"""Feasible crystal descent at fixed shape, absolute spatial order and carrier fields."""
from collections import Counter
from math import fsum

import numpy as np

from mouth_whole_r4_relax import SOLID_TOL, inverse_step, material_norm


TINY_ALPHA = 1.e-4
NEGLIGIBLE_TOTAL_ULPS = 32
REPEATED = 2


def _same_bytes(a, b):
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def _background_energy(trial, u):
    """The part of ContinuousTrial.fixed independent of crystal coordinates."""
    model = trial.model
    previous = model.material_load
    try:
        model.material_load = trial.ext
        return float(model.evaluate(u)[0])
    finally:
        model.material_load = previous


def _trace(on_event, name, **row):
    if on_event is not None:
        on_event(dict(event=name, **row))


def feasible_crystal_block(trial, q, u, history, *, steps=36, prior_overlap=None, on_event=None):
    """Return the last feasible q and a JSON-serializable fixed-field search report.

    alpha is relative to the original trust-capped direction (max displacement
    0.025); solver_scale includes that cap. Energy differences use the varying
    crystal term alone: adding the same field energy before subtraction would
    erase resolvable descent near stationarity. The Armijo check is strict.

    A trusted archived rejection can be reused only for identical material
    coordinates, shape and all spatial fields. An undecided certificate or a
    disagreement between an exact certificate and native quadrature is an error.
    """
    if steps < 1:
        raise ValueError("A material block needs at least one search step")
    status = trial.certificate(q, u[2], False)
    if status["status"] != "certified":
        raise ArithmeticError(f"Uncertified material block input: {status['status']} {status['witness']}")
    crystal = trial.b
    field_energy = _background_energy(trial, u)
    energy, grad, _, _ = crystal.evaluate(q, u[2], load=False)
    initial_energy = float(energy)
    reference_total = fsum((field_energy, initial_energy))
    negligible = float(NEGLIGIBLE_TOTAL_ULPS*np.spacing(max(1., abs(reference_total))))
    archive = prior_overlap
    if archive is not None and not (_same_bytes(archive['x'], trial.x) and _same_bytes(archive['field'], u)):
        archive = None
    accepted = []
    reasons = Counter()
    first_overlap = None
    last_overlap = None
    last_attempt = None
    consecutive_tiny = 0
    status = 'advanced'

    for iteration in range(steps):
        if material_norm(crystal.residual(q, grad)) < SOLID_TOL:
            status = 'material_tolerance'
            break
        quasi = inverse_step(grad, history, crystal.metric)
        diagonal = grad/crystal.metric
        directions = [('lbfgs', quasi)]
        if not np.array_equal(quasi, diagonal):
            directions.append(('diagonal', diagonal))
        took_step = False

        for name, direction in directions:
            if float(grad@direction) <= 0:
                reasons['not_descent_direction'] += 1
                continue
            cap = min(1., .025/max(1.e-20, float(np.max(abs(direction)))))
            alpha = 1.
            while True:
                scale = cap*alpha
                step = -scale*direction
                candidate = q+step
                if np.array_equal(candidate, q):
                    reasons['unrepresentable_step'] += 1
                    break
                slope = float(grad@step)
                attempt = dict(iteration=iteration, direction=name, alpha=alpha,
                               solver_scale=scale, directional_work=slope)
                cached = archive is not None and _same_bytes(archive['q'], candidate)
                if cached:
                    verdict = dict(status='active_fold', witness=archive['witness'])
                    archive = None
                else:
                    verdict = trial.certificate(candidate, u[2], False)
                reason = verdict['status']
                if reason == 'certified':
                    new_energy, new_grad, _, _ = crystal.evaluate(candidate, u[2], load=False)
                    change = float(new_energy-energy)
                    attempt['delta_total_energy'] = change
                    if change < min(0., 1.e-4*slope):
                        took_step = True
                        break
                    reason = 'insufficient_strict_decrease'
                elif reason == 'active_fold':
                    witness = verdict['witness']
                    if not (witness['exact_order_positive'] and witness['exact_determinant_sign'] <= 0):
                        raise ArithmeticError(f'Unverified overlap: {witness}')
                    overlap = dict(**attempt, witness=witness, reused_archived_rejection=cached)
                    if first_overlap is None:
                        first_overlap = overlap
                    last_overlap = overlap
                    _trace(on_event, 'crystal_overlap_rejected', **overlap)
                elif reason != 'outside_chart':
                    raise ArithmeticError(f"Undecided material placement: {reason} {verdict['witness']}")
                reasons[reason] += 1
                last_attempt = dict(**attempt, reason=reason)
                # A tiny infeasible step can have a feasible half, even when
                # its predicted work is small. Search until the change itself
                # falls below machine resolution in the unit material chart.
                if np.max(abs(step)) <= np.finfo(float).eps and abs(slope) <= negligible:
                    reasons['below_chart_and_work_resolution'] += 1
                    break
                alpha *= .5
            if took_step:
                break

        if not took_step:
            status = 'no_feasible_descent'
            break
        delta = new_grad-grad
        product = float(delta@step)
        if product > 1.e-12*np.linalg.norm(step)*np.linalg.norm(delta):
            history.append((step, delta, 1./product))
            del history[:-16]
        consecutive_tiny = consecutive_tiny+1 if alpha < TINY_ALPHA and -change <= negligible else 0
        row = dict(**attempt, tiny_without_meaningful_drop=bool(consecutive_tiny),
                   crystal_energy=float(new_energy), energy=fsum((field_energy, new_energy)))
        accepted.append(row)
        _trace(on_event, 'crystal_step_accepted', **row)
        q, energy, grad = candidate, new_energy, new_grad
        if consecutive_tiny >= REPEATED:
            status = 'tiny_stall'
            break

    residual = crystal.residual(q, grad)
    if material_norm(residual) < SOLID_TOL:
        status = 'material_tolerance'
    return q, dict(status=status, accepted=accepted, steps=len(accepted),
                   energy_before=reference_total, energy_after=fsum((field_energy, energy)),
                   delta_total_energy=float(energy-initial_energy), negligible_drop=negligible,
                   tiny_alpha=TINY_ALPHA, tiny_streak=consecutive_tiny,
                   min_alpha=min((row['alpha'] for row in accepted), default=None),
                   last_attempt=last_attempt, rejection_counts=dict(reasons),
                   first_overlap=first_overlap, last_overlap=last_overlap,
                   material=residual)
