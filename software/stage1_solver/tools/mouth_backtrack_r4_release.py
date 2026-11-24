"""Backtrack internal crystal displacements of the held R4 mouth after checkpoint 10."""
import argparse
from pathlib import Path

import numpy as np

from mouth_continuous_shape_relax import ContinuousRelease
from mouth_exact_melt_trial import ExactMeltTrial
from mouth_exact_r4_release import ExactRelease
from mouth_feasible_crystal_step import REPEATED, TINY_ALPHA, feasible_crystal_block
from mouth_joint_endpoint_audit import geometry
from mouth_ordered_j_bounds import ordered_j_readout, ordered_prism_point
from mouth_rim_pinned_comparison import digest, read
from mouth_whole_r4_relax import SHAPE_TOL, atomic_json, log, save_fields


class FeasibleOrderTrial(ExactMeltTrial):
    """At fixed crystal coordinates, the exact guard defines a face of the order box.

    A nonnegative spatial basis with positive weight anywhere on J <= 0 must
    have zero coefficient. Other free coefficients see the original energy.
    """
    def match(self, q, u, *args, **kwargs):
        previous = self.free
        self.matching_order = True
        try:
            return super().match(q, u, *args, **kwargs)
        finally:
            self.matching_order = False
            self.free = previous

    def crystal_load_for_match(self, q, u):
        load = super().crystal_load_for_match(q, u)
        if not getattr(self, 'matching_order', False):
            return load
        c = u[2, self.first]
        eligible = self.affected[(c[self.affected] == 0) & self.free[self.affected]]
        pending = [eligible] if len(eligible) else []
        blocked = []
        while pending:
            ids = pending.pop()
            # A certified union proves each member disjoint. A witnessed union
            # is subdivided down to singletons before constraining any node.
            indicator = np.isin(self.number, ids).astype(float)
            evidence = self.certificate(q, indicator, False)
            if evidence['status'] == 'certified':
                continue
            if evidence['status'] != 'active_fold':
                raise ArithmeticError(f"Undecided order-support mask: {evidence['status']}")
            if len(ids) == 1:
                blocked.append(int(ids[0]))
            else:
                cut = len(ids)//2
                pending.extend((ids[:cut], ids[cut:]))
        self.blocked = np.asarray(sorted(blocked), dtype=int)
        self.free = self.free.copy()
        self.free[self.blocked] = False
        return load

    def check_settled_fields(self, q, u, gradient):
        checks = super().check_settled_fields(q, u, gradient)
        if getattr(self, 'matching_order', False):
            derivatives = (self.pull.T @ gradient[2])[self.blocked]
            if self.blocked.size and np.any(u[2, self.first[self.blocked]] != 0.):
                raise ArithmeticError('Infeasible order on the fixed crystal')
            checks['excluded_support'] = dict(nodes=self.blocked.tolist(),
                lower_bound_dE_dchi=derivatives.tolist(), negative_count=int((derivatives < 0).sum()))
        return checks


class MaterialStall(ArithmeticError):
    """Both authorized update orders reached the recorded material stall."""


class BacktrackRelease(ExactRelease):
    trial_type = FeasibleOrderTrial

    def __init__(self, source, replay, controls, previous, folder):
        gate = read(previous/'coarse/full_trial_gate.json')
        failed = read(previous/'coarse/overlap.json')
        outcome = read(previous/'guard/outcome.json')
        if (not outcome['limitsVerified'] or outcome['exitCode'] != 1 or
                read(previous/'coarse/issue.json')['error'] != 'VerifiedOverlap' or
                digest(previous/'coarse/overlap_state.npz') != failed['saved_state_sha256'] or
                gate['replayed_sha'] != digest(replay/'trial.npz') or
                failed['checkpoint_sha'] != gate['checkpoint_sha']):
            raise ValueError('Original failed displacement or verified full shape step changed')
        super().__init__(source, replay, controls, folder)
        with np.load(previous/'coarse/overlap_state.npz', allow_pickle=False) as snap:
            if not np.array_equal(snap['x'], self.x) or not np.array_equal(snap['field'], self.u):
                raise ValueError('Prior displacement does not start from this rematched field')
            self.prior_overlap = dict(q=snap['q'].copy(), x=self.x.copy(), field=self.u,
                                      witness=failed['witness'])
        self.previous_pin = None
        self.previous_force = np.asarray(self.last['shape_force'])
        self.last_accepted = dict(iteration=10, sha256=self.last['sha256'], energy=self.last['energy'])
        self.block_reports = []
        self.block_stage = None

    def build(self, x, solve_fields=False):
        # Candidate certificates return their verdict to the line search. A
        # rematched positive-order fold still raises before any data are kept.
        return ContinuousRelease.build(self, x, solve_fields)

    def block(self, trial, q, u, steps, stage):
        self.block_stage = stage
        next_q, rep = feasible_crystal_block(trial, q, u, self.history, steps=steps,
            prior_overlap=self.prior_overlap,
            on_event=lambda row: log('crystal_search', checkpoint=self.iteration, stage=stage, **row))
        self.prior_overlap = None
        rep['stage'] = stage
        self.block_reports.append(rep)
        log('crystal_block', checkpoint=self.iteration, stage=stage, status=rep['status'],
            accepted=rep['steps'], min_alpha=rep['min_alpha'], energy_change=rep['delta_total_energy'])
        return next_q, rep

    def rematch(self, trial, q, u, stage, crystal_delta=0.):
        result = trial.match(q, u)
        field, e, _, load, meta = result
        # Compare with exactly the same assembled energy at this q. The large
        # common bath can obscure a few ulps in a field-only difference.
        trial.model.material_load = load
        initial = trial.transported(u)
        baseline, old_grad = trial.model.evaluate(initial)
        raw_change = float(e-baseline)
        roundoff = 32*np.spacing(max(1., abs(baseline)))
        change, error = raw_change, 0.
        if raw_change > roundoff:
            raise ArithmeticError(f'Field rematch increased total energy by {raw_change}')
        if crystal_delta+raw_change > 0:
            # Resolve work smaller than an ulp of the common bath. All fields
            # along this segment satisfy the box and exact support constraints.
            dx = field-initial
            slopes = []
            for fraction in (0., .25, .5, .75, 1.):
                grad = old_grad if fraction == 0. else trial.model.evaluate(initial+fraction*dx)[1]
                slopes.append(float(np.sum(grad.astype(np.longdouble)*dx, dtype=np.longdouble)))
            two = (slopes[0]+4*slopes[2]+slopes[4])/6
            four = (slopes[0]+4*slopes[1]+2*slopes[2]+4*slopes[3]+slopes[4])/12
            error = abs(four-two)+64*np.finfo(float).eps*sum(map(abs, slopes))
            change = four
            if abs(change-raw_change)>roundoff+error:
                raise ArithmeticError('Field work disagrees with the assembled energy')
            if crystal_delta+change+error > 0:
                # The incoming field is an admissible alternative only if it
                # already satisfies the field equations at the NEW crystal.
                mask = trial.free
                trial.free = mask.copy()
                trial.free[trial.blocked] = False
                trial.matching_order = True
                try:
                    chi_res = float(np.max(abs(trial.projected(initial[2,trial.first],trial.pull.T@old_grad[2]))))
                    carrier_res = max(float(np.max(abs(trial.model.solvers[j](
                        old_grad[j,trial.model.masks[j]])))) for j in (0,1))
                    if max(chi_res,carrier_res)>=2.e-7:
                        raise ArithmeticError('Combined descent unresolved and incoming fields are not stationary')
                    original_checks = trial.check_settled_fields(q, initial, old_grad)
                finally:
                    trial.free,trial.matching_order = mask,False
                solid, gq, _, _ = trial.b.evaluate(q, initial[2], load=False)
                exchange = float(solid-(load-trial.ext)@initial[2])
                if abs(exchange)>2.e-8*max(1.,abs(solid)):
                    raise ArithmeticError('Unbalanced work with the retained feasible fields')
                field, e = initial, float(baseline)
                meta = dict(field=dict(stationarity=max(chi_res,carrier_res),
                    carrier_metric_residual=carrier_res, order_projected_residual=chi_res,
                    retained_incoming_field=True),solid_energy=solid,exchange_error=exchange,
                    energy=e, **original_checks)
                result = field,e,gq,load,meta
                change,error = 0.,0.
        meta['field_update'] = dict(stage=stage, energy_before=float(baseline), delta_total_energy=change,
            raw_candidate_field_delta=raw_change, field_work_error_estimate=error, crystal_delta=crystal_delta,
            order51603_before=float(u[2, trial.first[51603]]),
            order51603_after=float(field[2, trial.first[51603]]))
        if not np.array_equal(field[2], u[2]):
            self.history.clear()
        return result

    def stalled(self, trial, q, u, reason, recent=()):
        folder = self.args.folder
        snapshot = folder/f'stalled_{self.iteration:05d}.npz'
        save_fields(snapshot, x=trial.x, q=q, solid=trial.b.unpack(q), field=u,
                    shape_params=trial.shape.params, shape_extra=trial.shape.extra,
                    state_sha=digest(self.args.state), radius=4, winding=1)
        entry = dict(reason=reason, attempted_checkpoint=self.iteration, last_accepted=self.last_accepted,
                     stall_state=snapshot.name, sha256=digest(snapshot), recent=list(recent),
                     material_blocks=self.block_reports, previous_pin=self.previous_pin)
        atomic_json(folder/'stall.json', entry)
        entry['ordered_J'] = ordered_j_readout(trial, q, u, folder)
        entry['historical_point'] = ordered_prism_point(trial, q, u, folder)
        last = next((b['last_overlap'] for b in reversed(self.block_reports) if b['last_overlap']), None)
        if last:
            entry['last_rejected_witness_on_feasible_state'] = ordered_prism_point(
                trial, q, u, folder, witness=last['witness'])
        atomic_json(folder/'stall.json', entry)
        raise MaterialStall(reason)

    def material_update(self, trial, q, u):
        self.block_reports = []
        self.interleaved = False
        # Keep the original update order unless the previous complete state
        # already demonstrated positive-order, machine-scale volume.
        if self.previous_pin is None:
            q, pure = self.block(trial, q, u, 36, 'fixed_order')
            fields = self.rematch(trial, q, u, 'after_fixed_order', pure['delta_total_energy'])
            pure['field_update'] = fields[4]['field_update']
            if pure['status'] not in ('tiny_stall','no_feasible_descent'):
                return q, pure['steps'], pure['status'], fields
            u = fields[0]
            reason = pure['status']
        else:
            fields = self.rematch(trial, q, u, 'order_first_after_small_j')
            u = fields[0]
            reason = 'positive-order machine-scale volume at previous checkpoint'
        self.interleaved = True
        self.history.clear()
        log('interleaving', checkpoint=self.iteration, reason=reason,
            order51603=float(u[2, trial.first[51603]]))
        recent = []
        for _ in range(36):
            q, row = self.block(trial, q, u, 1, 'interleaved')
            fields = self.rematch(trial, q, u, 'interleaved', row['delta_total_energy'])
            u = fields[0]
            change = fields[4]['field_update']
            row['field_update'] = change
            total_delta = row['delta_total_energy']+change['delta_total_energy']
            tiny = ((row['min_alpha'] is not None and row['min_alpha'] < TINY_ALPHA)
                    or row['status'] == 'no_feasible_descent')
            if row['status'] == 'material_tolerance':
                break
            if tiny and total_delta >= -row['negligible_drop']:
                recent.append(dict(alpha=row['min_alpha'], delta_total_energy=total_delta,
                                   last_attempt=row['last_attempt'], field_update=change))
                if len(recent) >= REPEATED:
                    self.stalled(trial, q, u, 'Repeated negligible interleaved energy descent', recent)
            else:
                recent.clear()
        return q, sum(b['steps'] for b in self.block_reports), 'interleaved', fields

    def begin_material_record(self, trial, q, u, energy, internal, info):
        # Preserve a recoverable state if a signed probe fails. It is never
        # published as a completed checkpoint until all six forces exist.
        path = self.args.folder/'pending.npz'
        save_fields(path, x=trial.x, solid=trial.b.unpack(q), q=q, field=u, energy=energy,
                    shape_params=trial.shape.params, shape_extra=trial.shape.extra,
                    state_sha=digest(self.args.state), radius=4, winding=1)

    def complete_material_record(self, trial, q, u, energy, internal, info, values):
        folder = self.args.folder
        current = values['current']
        checks = trial.check_material_domain(q, u)
        forces = np.asarray(values['shape_force'])
        if checks['status'] != 'certified' or forces.shape != (6,) or not np.isfinite(forces).all():
            raise ArithmeticError('Incomplete signed shape work or order certificate')
        alpha = [b['min_alpha'] for b in self.block_reports if b['min_alpha'] is not None]
        trend = (np.abs(forces)<np.abs(self.previous_force)).tolist()
        measure = ordered_j_readout(trial, q, u, folder)
        pin = measure['machine_zero_witness']
        name = f'state_{self.iteration:05d}.npz'
        record = dict(iteration=self.iteration, file=name, x=trial.x.tolist(), energy=energy,
            material=internal, orientation=checks, **geometry(trial.shape, trial.polar_shift),
            **values, ordered_J=measure, min_material_alpha=min(alpha, default=None),
            melted_front_depth=trial.shape.length*current['shaft_order']['terminal_from_fraction']['disk']['0.0']
                if current['shaft_order']['terminal_from_fraction']['disk']['0.0'] is not None else None,
            shape_magnitudes_decreased=trend, exit_and_middle_magnitudes_both_decreased=trend[1] and trend[2],
            material_blocks=self.block_reports, interleaved=self.interleaved)
        if energy >= self.last_accepted['energy'] and (values['material_steps'] or
                not np.array_equal(trial.x, self.last['x'])):
            raise ArithmeticError('Net checkpoint energy decrease is unresolved')
        old = self.last
        self.last = record
        try:
            self.check_progress()
        except Exception:
            self.last = old
            raise
        path = folder/name
        if path.exists():
            raise FileExistsError(path)
        (folder/'pending.npz').rename(path)
        record['sha256'] = digest(path)
        atomic_json(folder/path.with_suffix('.json').name, record)
        atomic_json(folder/'latest.json', record)
        self.last_accepted = dict(iteration=self.iteration, file=name, sha256=record['sha256'], energy=energy)
        self.previous_force = forces
        log('complete_checkpoint', iteration=self.iteration, energy=energy, shape_force=forces.tolist(),
            min_alpha=record['min_material_alpha'], current=current['feed_current'], J=measure['infimum_bounds'])
        if pin and self.previous_pin is not None and self.interleaved:
            self.stalled(trial, q, u, 'Positive-order machine-scale volume persisted after interleaving')
        self.previous_pin = dict(iteration=self.iteration, witness=measure['upper_witness']) if pin else None

    def shape_update(self, x, dx, q, u):
        before = self.x.copy()
        result = super().shape_update(x, dx, q, u)
        if not np.array_equal(before, result[0]):
            self.history.clear()
        return result

    def finish(self, trial, q, u, load, measured):
        if self.last['ordered_J']['machine_zero_witness'] or max(map(abs,self.last['half_shape_force']))>=SHAPE_TOL:
            raise ArithmeticError('Independent half-step force or ordered-volume check remains unresolved')
        if self.last['field'].get('excluded_support', {}).get('negative_count',0):
            raise ArithmeticError('One-sided order work on excluded support remains unresolved at endpoint')
        super().finish(trial, q, u, load, measured)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source','replay','controls','previous','output'):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    folder = args.output.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if any(folder.iterdir()):
        raise FileExistsError('Use an empty output folder and preserve previously computed states')
    try:
        run = BacktrackRelease(args.source.resolve(), args.replay.resolve(),
            args.controls.resolve(), args.previous.resolve(), folder)
        run.run(matched_energy=run.energy, last_x=run.parent_x, last_g=run.parent_g)
    except Exception as exc:
        record = dict(error=type(exc).__name__, detail=str(exc), last_accepted='parent checkpoint 10')
        if (folder/'latest.json').exists():
            record['last_accepted'] = read(folder/'latest.json')['file']
        atomic_json(folder/'issue.json', record)
        raise
