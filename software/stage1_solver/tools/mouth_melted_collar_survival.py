"""Continue the saved R4 coarse solution with solid orientation on ordered material.

The spatial order variation stays free. For an invalid placement confined to
zero order we test its most favorable possible reordering cost, W=0, only
when the former W is undefined. Any oriented crystal has W >= 0. A settled
zero with a nonnegative one-sided derivative for this lower bound is also
stationary for every nonnegative repair energy. Reactivation lacks that proof.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_active_prism_orientation import active_orientation
from mouth_continuous_depth_mesh import neutral_force
from mouth_continuous_order import norm
from mouth_continuous_shape_relax import ContinuousRelease, ShapeTrial
from mouth_exchange_carrier import BIAS
from mouth_free_shape_solid import InvalidPlacement
from mouth_joined_front_audit import Front
from mouth_last_rings_fluid import GRIDS
from mouth_r4_fold_order import describe
from mouth_rim_pinned_comparison import digest, read, validate_guard
from mouth_whole_r4_relax import FIELD_TOL, SOLID_TOL, STEP, atomic_json, log, material_norm, relax_material


START_SHA = '988d5c2e3d5124dd3caaf7b00f206031ee15e17bf81b2b2d64798d400e10ce12'


class ReorderingUnresolved(ArithmeticError):
    """The relaxed order requires an oriented reference or a sharper force bound."""


def brief(certificate):
    return {k: v for k, v in certificate.items() if k != 'potential_activation_nodes'}


class MeltTrial(ShapeTrial):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.b.melted_labels = True

    def certificate(self, q, order, collect):
        return active_orientation(self.b, self.b.unpack(q), order, collect_activation=collect)

    def check_material_domain(self, q, u):
        return brief(self.certificate(q, u[2], False))

    def fixed(self, q, u):
        if not self.check_material_domain(q, u)['all_active_cells_oriented']:
            raise InvalidPlacement('Infeasible ordered material orientation in virtual work')
        return super().fixed(q, u)

    def crystal_load_for_match(self, q, u):
        cert = self.certificate(q, u[2], True)
        if not cert['all_active_cells_oriented']:
            raise ReorderingUnresolved(f'Input order/placement: {brief(cert)}')
        _, _, load, details = self.b.lower_bound_load(q, u[2])
        ids = np.union1d(cert['potential_activation_nodes'], details['undefined_order_nodes']).astype(int)
        self.affected = np.unique(self.number[ids])
        self.invalid_points = details['nonpositive_J_points']
        return load

    def check_settled_fields(self, q, u, gradient):
        cert = self.certificate(q, u[2], True)
        if not cert['all_active_cells_oriented']:
            raise ReorderingUnresolved(f'Rematched order/placement: {brief(cert)}')
        affected = np.union1d(self.affected, self.number[np.asarray(cert['potential_activation_nodes'], dtype=int)])
        c = u[2, self.first]
        g = self.pull.T @ gradient[2]
        touched = affected[self.free[affected]]
        if np.any(c[affected] != 0.) or np.any(g[touched] < 0.):
            i = int(affected[np.argmax(c[affected])]) if np.any(c[affected] != 0.) else int(touched[np.argmin(g[touched])])
            raise ReorderingUnresolved(f'Cannot establish zero-order stability at joined node {i}: chi={c[i]!r}, lower-bound order force={g[i]!r}, nodes={len(affected)}, nonpositive quadrature sites={self.invalid_points}')
        free = self.free
        kkt = np.where(c == 0, np.minimum(g, 0.), np.where(c == 1., np.maximum(g, 0.), g))
        raw = {name: norm(gradient[i, self.model.masks[i]]) for i, name in enumerate(('amplitude', 'phase'))}
        raw.update(order_free=norm(g[free]), order_kkt=norm(kkt[free]),
                   order_interior=norm(g[free & (c>0) & (c<1)]))
        return dict(raw_residual=raw, zero_order_check=dict(nodes=len(affected), free_nodes=len(touched),
                    minimum_lower_bound_derivative=float(np.min(g[touched])) if touched.size else None,
                    invalid_quadrature_points=self.invalid_points), material_domain=brief(cert))

    def measure(self, u, load):
        result = super().measure(u, load)
        m = self.model
        region = u[2, m.ids[1, 0]]
        result['shaft_order'] = describe((m.rr[0], m.a, m.z[1]/self.shape.length), region)
        return result


class Survival(ContinuousRelease):
    trial_type = MeltTrial

    def __init__(self, args):
        self.args = args
        data = read(args.source/'coarse/latest.json')
        source = args.source/'coarse'/data['file']
        if data['iteration'] != 6 or digest(source) != START_SHA or data['sha256'] != START_SHA:
            raise ValueError('Expected the retained update 6')
        command = read(args.source/'coarse/invocation.json')['command']
        original, _, edges, folder, level = command[-5:]
        if Path(folder) != source.parent or level != 'coarse':
            raise ValueError('Source mesh does not match')
        prior = Path(original)/level
        proof = read(Path(edges))['results'][level]
        args.state = Path(read(prior/'invocation.json')['command'][5])
        if digest(args.state) != proof['preparation_sha256'] or digest(prior/'endpoint.json') != proof['endpoint_sha256']:
            raise ValueError('Previous radius-four preparation changed')
        self.discretization = proof['solid_discretization']
        if self.discretization != dict(n=8, k=24, angles=6, degree=3, subdivisions=2) or BIAS != .4510084282857874:
            raise ValueError('Incorrect mesh or reservoir phase drop')
        self.base = Front(args.state, GRIDS['double'], 1)
        self.reference_depth = read(prior/'endpoint.json')['depth']
        with np.load(source) as saved:
            if saved['winding'].item() != 1 or saved['radius'].item() != 4 or str(saved['state_sha']) != digest(args.state):
                raise ValueError('Incorrect winding, radius or preparation')
            self.q, self.x, self.u = (saved[key].copy() for key in ('q','x','field'))
            offsets = saved['solid']
        self.trial = self.build(self.x, True)
        if not np.array_equal(offsets, self.trial.b.unpack(self.q)):
            raise ValueError('Reference labels changed on loading checkpoint')
        self.iteration, self.last, self.history, self.shape_metric = 6, None, [], None
        self.previous = data

    def preflight(self):
        t, q, u = self.trial, self.q, self.u
        b, m = t.b, t.model
        results = []
        for active in (False, True):
            b.melted_labels = active
            e, force, load, _ = b.evaluate(q, u[2])
            if active:
                _, _, low, details = b.lower_bound_load(q, u[2])
                if details['undefined_order_nodes'] or not np.array_equal(load, low):
                    raise ArithmeticError('Update 6 order load differs from full positive-J deposition')
            m.material_load = load + t.ext
            total, gradients = m.evaluate(u)
            results.append((e, force, load, total, gradients, m.measure(u)))
        old, new = results
        differences = dict(solid_energy=abs(new[0]-old[0]), solid_force=norm(new[1]-old[1]),
                           order_load=norm(new[2]-old[2]), total_energy=abs(new[3]-old[3]),
                           full_field_gradient=norm(new[4]-old[4]),
                           port_currents=max(abs(new[5][k]-old[5][k]) for k in ('feed_current','bulk_current','weak_disk_current')))
        source_error = dict(total_energy=abs(new[3]-self.previous['energy']),
                            ring_residual=abs(material_norm(b.residual(q,new[1]))-material_norm(self.previous['material'])),
                            intake=abs(new[5]['feed_current']-self.previous['current']['feed_current']))
        validate_guard(self.args.reference.parent/'reference_guard_host')
        independent = read(self.args.reference.with_suffix('.json'))
        if independent['checkpoint_sha256'] != START_SHA:
            raise ValueError('Independent baseline uses another state')
        with np.load(self.args.reference) as ref:
            source_error.update(all_crystal_forces=norm(new[1]-ref['gq']),
                                all_crystal_to_order_loads=norm(new[2]-ref['inner']),
                                joined_order_first_variation=norm(t.pull.T@new[4][2]-ref['raw']),
                                all_field_gradients=norm(new[4]-ref['field_gradient']))
        cert = t.certificate(q,u[2],True)
        record = dict(source_sha256=START_SHA, equivalence=differences, archived_discrepancy=source_error,
                      total_energy=new[3], signed_ports=[new[5][k] for k in ('feed_current','bulk_current','weak_disk_current')],
                      source_material=b.residual(q,new[1]), domain=brief(cert),
                      raw_order_load_max=norm(new[2]), affected_nodes=len(cert['potential_activation_nodes']),
                      independent_reference_sha256=digest(self.args.reference))
        atomic_json(self.args.folder/'saved_state_comparison.json',record)
        if (max(differences.values())>1.e-10 or max(source_error.values())>1.e-8 or
                not cert['all_active_cells_oriented'] or new[5]['port_error']>3.e-5):
            raise ArithmeticError(f'Saved-state equivalence failed: {record}')
        log('saved_state_equivalence',**record)

    def initial_material(self, t, q, u):
        return dict(orientation=t.check_material_domain(q, u))

    def material_update(self, t, q, u):
        candidate, steps, why = relax_material(t, q, u, self.history, 36)
        delta = candidate-q
        scale = 1.
        while True:
            try:
                result = t.match(candidate, u)
                return candidate, steps, why, result
            except ReorderingUnresolved as exc:
                self.history.clear()
                log('rejected_material_rematch',scale=scale, detail=str(exc))
                scale *= .5
                candidate = q+scale*delta
                if np.array_equal(candidate,q):
                    return q, 0, str(exc), t.match(q,u)

    def shape_update(self, x, dx, q, u):
        scale = 1.
        while True:
            new_x = x+scale*dx
            t = self.build(new_x, True)
            try:
                return new_x, t, t.match(q,u)
            except ReorderingUnresolved as exc:
                log('rejected_shape_rematch',scale=scale,detail=str(exc))
                del t
                scale *= .5
                if np.array_equal(x+scale*dx,x):
                    raise ArithmeticError(f'No resolved ordered-domain shape step: {exc}') from exc

    def check_progress(self):
        row = self.last
        flow = row['current']
        if flow['port_error']>3.e-5 or abs(flow['disk_minus_feed_current'])>3.e-5:
            raise ArithmeticError('Port or weak disk flux mismatch on accepted state')
        sections = flow['shaft_order']['sections']
        if sections[0]['wall'][1] == 0.:
            raise ArithmeticError('Complete loss of order around the clamped entrance ring; inspect saved front trajectory')
        if row['minimum_radius']<=0. or flow['feed_current']<=0.:
            raise ArithmeticError('Positive draining branch lost; inspect geometry and current')

    def finish(self,t,q,u,load,measured):
        checks=t.check_material_domain(q,u)
        if (not checks['all_active_cells_oriented'] or measured['dense_total_subsonic_share']<=0. or
                material_norm(self.last['material'])>=SOLID_TOL or self.last['field']['field']['stationarity']>=FIELD_TOL):
            raise ArithmeticError('Coarse equilibrium evidence incomplete')
        record=dict(self.last,neutral_depth_null=neutral_force(t.shape.length,self.reference_depth,STEP/2,24))
        atomic_json(self.args.folder/'endpoint.json',record)
        log('coarse_endpoint',sha256=digest(self.args.folder/'endpoint.json'))

    def failure(self, exc):
        atomic_json(self.args.folder/'issue.json',dict(error_type=type(exc).__name__,error=str(exc),last=self.last))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('folder',type=Path)
    parser.add_argument('reference',type=Path)
    args=parser.parse_args()
    run=Survival(args)
    try:
        run.preflight()
        run.run()
    except (ArithmeticError,ValueError,RuntimeError) as exc:
        run.failure(exc)
        raise
