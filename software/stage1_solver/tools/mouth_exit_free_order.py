"""Fixed-shape relaxation of every unpinned order value on the occupied bonded branch."""
from dataclasses import dataclass

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_live_waist_reopen import reopen


TOLERANCE = 2e-5  # dE/d(sqrt(H_ii) * chi_i), as in the bounded-order benchmark.
FIELD_ERRORS = (ArithmeticError, ValueError, RuntimeError, np.linalg.LinAlgError)


@dataclass
class Opening:
    index: tuple
    ids: np.ndarray
    extension: np.ndarray
    slope: float
    lower: float | None


def opening_slopes(chi, free, gradient, wave):
    """True coordinate slopes and a sufficient bound for simultaneous openings."""
    omega, vector = wave.seed
    project = wave.project.tocsc()
    weights = project @ chi.T.ravel()
    touched = np.zeros(chi.shape, dtype=bool)
    if np.min(weights) > 0:
        return touched, [], []
    mass = wave.cut_m + wave.matrix(weights[:, None, None]*wave.mc)
    active = mass.diagonal() > 0
    if np.any(vector[~active] != 0):
        raise ArithmeticError("The pole has nonzero amplitude on silent wave nodes")
    cells = np.any(~np.r_[active, True][wave.ids], axis=1)
    touched = ((project[cells].sum(axis=0).A1 > 0).reshape(chi.T.shape).T & free)
    if not touched.any():
        return touched, [], []
    if np.any(chi[touched] != 0):
        raise ArithmeticError("Positive order touches a structurally silent wave node")
    edge = vector[wave.port.ids]
    derivative = (wave.port.block(omega+2e-6)-wave.port.block(omega-2e-6))/4e-6
    den = 2*omega*(vector @ (mass @ vector)) + edge @ (derivative @ edge)
    if not np.isfinite(den) or den == 0:
        raise ArithmeticError("Singular pole derivative at the reopening boundary")
    pencil = wave.kc - omega**2*wave.mc
    openings, failures = [], []
    for iz, ir in np.argwhere(touched):
        index = (int(iz), int(ir))
        column = ir*chi.shape[0] + iz
        section = slice(project.indptr[column], project.indptr[column+1])
        ids, coeff = project.indices[section], project.data[section]
        try:
            _, correction, new, ext, _ = reopen(ids, coeff, wave, pencil, vector, active)
            slope = float(gradient[index] + (correction/den).real)
            if not np.isfinite(slope):
                raise ArithmeticError("Nonfinite one-sided slope")
            # For positive combinations of openings the same silent wave nodes
            # must extend consistently. If H_j = Re A_j is positive on the new
            # nodes, |b^T A^-1 b| <= b* H^-1 b <= sum_j t_j b_j* H_j^-1 b_j.
            # Here A=sum t_j A_j, b=sum t_j b_j, and the second inequality is
            # the matrix fractional bound (Cauchy-Schwarz on the direct sum).
            # It controls mixed directional work without assuming additivity.
            block = np.zeros((len(new), len(new)), complex)
            for cell, factor in zip(ids, coeff):
                nodes = wave.ids[cell]
                keep = np.flatnonzero(np.isin(nodes, new))
                local = np.searchsorted(new, nodes[keep])
                block[np.ix_(local, local)] += factor*pencil[cell][np.ix_(keep, keep)]
            lower = float(gradient[index])
            if len(new):
                values = np.linalg.eigvalsh(block.real)
                if values[0] <= 1e-10*np.max(abs(values)):
                    lower = None
                else:
                    cross = -block @ ext
                    bound = np.vdot(cross, np.linalg.solve(block.real, cross)).real
                    lower -= float(bound/abs(den))
            openings.append(Opening(index, new, ext, slope, lower))
        except FIELD_ERRORS as exc:
            failures.append(dict(index=list(index), reason=str(exc)))
    return touched, openings, failures


def relax_order(chi, order, loads, elastic, record):
    """Return (accepted chi, info) and restore its pole and density on every exit.

    record(chi=array, seed=(omega, vector), **info) receives accepted states.
    inside_energy comprises interior order, static interior elastic work, carrier
    and occupied wave at the frozen shape. The receiving order is a constant.
    """
    chi = np.array(chi, dtype=float, copy=True)
    free = order.free
    if (chi.shape != free.shape or elastic.shape != chi.shape or not np.isfinite(chi).all()
            or not np.isfinite(elastic).all() or chi.min() < 0 or chi.max() > 1
            or np.max(abs(chi[~free]-order.seed[~free])) > 1e-10):
        raise ValueError("Inadmissible starting order or elastic field")
    scale = np.sqrt(2*order.stiff.diagonal()+1.5*order.lump).reshape(chi.shape)
    if np.min(scale[free]) <= 0:
        raise ValueError("Singular order scaling")
    wave = loads.wave

    def evaluate(c):
        eo, go = order.energy_gradient(c)
        ef, gf, fields = loads.evaluate(c, 1., BIAS)
        energy = float(eo + ef + np.sum(elastic*c))
        gradient = go + gf + elastic
        if (not np.isfinite(energy) or not np.isfinite(gradient).all()
                or not np.isfinite([fields[k] for k in ('mode_overlap', 'pole_residual',
                                         'current', 'terminal_current', 'wave_frequency')]).all()
                or fields['mode_overlap'] < .94 or fields['pole_residual'] > 1e-7
                or abs(fields['current']-fields['terminal_current']) >
                   1e-6*max(abs(fields['current']), 1e-12)):
            raise ArithmeticError("Unresolved energy, occupied pole or terminal current")
        return energy, gradient, fields

    def snapshot():
        return (wave.seed, wave.target, loads.carrier.f.copy())

    def restore(s, extended=None):
        wave.seed, wave.target, density = s
        loads.carrier.f = density.copy()
        if extended is not None:
            v = wave.seed[1].copy()
            v[extended.ids] = extended.extension
            wave.seed = (wave.seed[0], v)

    initial = snapshot()
    try:
        energy, grad, fields = evaluate(chi)
    except FIELD_ERRORS as exc:
        restore(initial)
        info = dict(stationary=False, reason=f"Initial state unresolved: {exc}",
                    accepted=False, iteration=0, inside_energy=None, order_residual=None)
        record(chi=chi, seed=wave.seed, **info)
        return chi, info

    origin = energy
    state = snapshot()
    history = []
    iteration = 0
    diagnostics = {}

    def report(reason, stationary=False):
        restore(state)
        info = dict(fields, **diagnostics, stationary=stationary, reason=reason,
                    accepted=True, iteration=iteration, inside_energy=energy,
                    energy_change=energy-origin)
        record(chi=chi, seed=wave.seed, **info)
        return chi, info

    while True:
        restore(state)
        try:
            singular, sites, failures = opening_slopes(chi, free, grad, wave)
        except FIELD_ERRORS as exc:
            diagnostics = dict(order_residual=None)
            return report(f"Order-opening audit failed: {exc}")
        smooth = free & ~singular
        projected = grad/scale
        projected[~smooth | ((chi == 0) & (grad >= 0)) | ((chi == 1) & (grad <= 0))] = 0.
        norm = float(np.max(abs(projected)))
        best = min(sites, key=lambda s: s.slope/scale[s.index], default=None)
        violation = 0. if best is None else max(0., -best.slope/scale[best.index])
        uncertified = [s for s in sites if s.lower is None or s.lower/scale[s.index] < -TOLERANCE]
        diagnostics = dict(order_residual=max(norm, violation) if not failures else None,
                           smooth_residual=norm, negative_opening_residual=float(violation),
                           lower_nodes=int(np.sum(free & (chi == 0))),
                           upper_nodes=int(np.sum(free & (chi == 1))),
                           silent_openings=int(singular.sum()), unresolved_openings=failures,
                           uncertified_joint_openings=len(uncertified),
                           minimum_opening_slope=None if best is None else best.slope,
                           minimum_scaled_opening_index=None if best is None else list(best.index))
        if iteration == 0 or iteration & (iteration-1) == 0 or iteration % 25 == 0:
            report("Accepted order field")
        opening = best if violation > max(norm, TOLERANCE) else None
        if max(norm, violation) <= TOLERANCE:
            # Check the energy itself on two scales before accepting first-order
            # stationarity, including a true resumed mode at the worst silent site.
            probes = []
            if smooth.any():
                idx = np.unravel_index(np.argmax(np.where(smooth, abs(projected), -1.)), chi.shape)
                sign = 1. if chi[idx] <= .5 else -1.
                probes.append((idx, sign, grad[idx]*sign, None))
            if best is not None:
                probes.append((best.index, 1., best.slope, best))
            checks = []
            for idx, sign, slope, extension in probes:
                secants = []
                try:
                    for h in (.0005, .00025):
                        c = chi.copy()
                        c[idx] += sign*h
                        restore(state, extension)
                        secants.append((evaluate(c)[0]-energy)/h)
                except FIELD_ERRORS as exc:
                    diagnostics['direct_checks'] = checks
                    return report(f"Unresolved signed order probe at {idx}: {exc}")
                error = abs(2*secants[1]-secants[0]-slope)
                checks.append(dict(index=list(map(int, idx)), sign=sign, slope=slope, secants=secants,
                                   extrapolation_error=error))
                if error > 2e-5*max(1., abs(slope)):
                    diagnostics['direct_checks'] = checks
                    return report("Direct order work disagrees with its one-sided slope")
            diagnostics['direct_checks'] = checks
            if failures or uncertified:
                return report("Individual order relaxation ended with uncertified simultaneous openings")
            if fields.get('wave_components', 1) > 1 and singular.any():
                return report("A separate shear component lacks a spectral-gap check on reopening")
            return report("Bounded first-order stationarity at fixed shape", stationary=True)

        if opening is not None:
            direction = np.zeros_like(chi)
            direction[opening.index] = min(.05, violation/scale[opening.index])
            search = [direction]
        else:
            # Limited-memory inverse curvature in the diagonally scaled field.
            q = projected[free].copy()
            coefficients = []
            for dx, dy, rho in reversed(history):
                a = rho*(dx @ q)
                coefficients.append(a)
                q -= a*dy
            for (dx, dy, rho), a in zip(history, reversed(coefficients)):
                q += dx*(a-rho*(dy @ q))
            direction = np.zeros_like(chi)
            direction[free] = -q/scale[free]
            direction[~smooth | ((chi == 0) & (grad >= 0))
                              | ((chi == 1) & (grad <= 0))] = 0.
            search = [direction, -projected/scale]

        accepted = False
        reason = "No resolved energy decrease"
        for direction in search:
            if not np.any(direction):
                continue
            step = min(1., .15/np.max(abs(direction)))
            while True:
                candidate = np.clip(chi + step*direction, 0., 1.)
                dx = candidate-chi
                predicted = (float(np.sum(grad*dx)) if opening is None
                             else float(opening.slope*dx[opening.index]))
                noise = 8*np.spacing(max(abs(energy), 1.))
                if predicted >= -noise or not np.any(dx):
                    break
                restore(state, opening)
                try:
                    new_energy, new_grad, new_fields = evaluate(candidate)
                    if new_energy < energy - max(noise, -1.e-4*predicted):
                        accepted = True
                        break
                except FIELD_ERRORS as exc:
                    reason = str(exc)
                step *= .5
            if accepted:
                break
        if not accepted:
            diagnostics['attempted_opening'] = None if opening is None else list(opening.index)
            return report(f"Order line search unresolved: {reason}")
        if (opening is not None or np.any((chi == 0) != (candidate == 0))
                or np.any((chi == 1) != (candidate == 1))):
            history.clear()
        else:
            dx = ((candidate-chi)*scale)[free]
            dy = ((new_grad-grad)/scale)[free]
            curvature = dx @ dy
            if curvature > 1e-12*np.linalg.norm(dx)*np.linalg.norm(dy):
                history.append((dx, dy, 1./curvature))
                history = history[-12:]
        state = snapshot()
        chi, energy, grad, fields = candidate, new_energy, new_grad, new_fields
        iteration += 1
