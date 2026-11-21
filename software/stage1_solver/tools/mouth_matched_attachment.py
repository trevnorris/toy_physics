"""Compare centered and offset attachments at the same local 0PN wall state."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from mouth_contact_crossing_force import REGISTRY_DRAG, REGISTRY_K, registry_force
from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import emit
from mouth_surviving_attachment import STEPS
from mouth_tagged_contact import SPAN


OFFSET = .125
ORDERS = (.75, .50, .30)


def force_pair(law, offset, sideways, chi):
    """Internal forces on the displaced feature. Work is per unit contact area."""
    r = np.hypot(offset, sideways)
    if np.any(r >= chi/2):
        raise ValueError('Comparison extends outside the shared intact region')
    if law == 'separable':
        return (registry_force(offset, chi, REGISTRY_K),
                registry_force(sideways, chi, REGISTRY_K))
    if law != 'central':
        raise ValueError('Unknown contact law')
    radial = registry_force(r, chi, REGISTRY_K)
    factor = np.divide(radial, r, out=np.zeros_like(r), where=r>0)
    return offset*factor, sideways*factor


def mobile_offsets(trace, start, gates):
    # The site sits upstream of the junction center. The transport clock belongs
    # to the junction, and the order is sampled at the upstream site.
    shift = -SPAN/2
    local_order = lambda w: trace.sample(trace.chi, w+shift)[0]
    def initial_balance(q):
        return OFFSET-q+registry_force(q, local_order(start), REGISTRY_K)
    initial = brentq(initial_balance, 0., OFFSET, xtol=1e-15)
    # Both start at the supported equilibrium, then the second loses its anchor.
    anchor_k = np.array([1., 0.])
    state = np.full(2, initial)
    edges = np.unique(np.r_[start, gates, trace.edges, trace.edges-shift])
    edges = edges[(edges>=start) & (edges<=gates[-1])]
    readings = {}
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        def velocity(w, q):
            clock_at = np.clip(w, np.nextafter(left, right), np.nextafter(right, left))
            f = anchor_k*(OFFSET-q)+registry_force(q, local_order(w), REGISTRY_K)
            return f/(REGISTRY_DRAG*trace.speed(clock_at))

        # q=0 and q=OFFSET bound both solutions by the sign of the force there.
        # A capture width >2*OFFSET over the segment therefore guarantees survival.
        mid = (left+right)/2
        orders = np.array([local_order(w) for w in (left, mid, right)])
        bernstein_middle = 2*orders[1]-(orders[0]+orders[2])/2
        if min(*orders, bernstein_middle)/2 <= np.hypot(OFFSET, max(STEPS)):
            raise ValueError('Cannot certify that both feature pairs stay attached')
        sol = solve_ivp(velocity, (left, right), state, method='DOP853', rtol=2e-12, atol=1e-14)
        if not sol.success:
            raise ArithmeticError(sol.message)
        state = sol.y[:, -1]
        if right in gates:
            readings[right] = state.copy()
    return initial, readings


def held_pair(law, chi, offset, step):
    r0, f0 = force_pair(law, offset, 0., chi)
    rm, fm = force_pair(law, offset, -step, chi)
    rp, fp = force_pair(law, offset, step, chi)
    center = force_pair(law, 0., step, chi)[1]
    if center >= 0 or fp >= 0:
        raise ValueError('Expected an attached restoring force in both preparations')
    return dict(step=step, centered_force=center, offset_force_plus=fp, offset_force_minus=fm,
                ratio=fp/center, center_slope=-center/step, offset_slope=-(fp-fm)/(2*step),
                baseline_sideways_force=f0, external_holding_force=-r0,
                holding_change_plus=r0-rp, holding_change_minus=r0-rm)


def mixed_work(chi, law):
    # Integrate the difference in sideways holding work along two paths in the
    # same live patch. The second path uses only the registry-force channel.
    nodes, weights = leggauss(4)
    step = STEPS[0]
    u = (1+nodes)/2
    df_side = force_pair(law, 0., step*u, chi)[1]-force_pair(law, OFFSET, step*u, chi)[1]
    df_row = force_pair(law, OFFSET*u, 0., chi)[0]-force_pair(law, OFFSET*u, step, chi)[0]
    return step*(weights @ df_side)/2, OFFSET*(weights @ df_row)/2


def main():
    for mesh in (34, 48):
        trace = WallTrace(mesh, .5)
        start = trace.crossing(.95)
        gates = np.array([trace.crossing(chi)+SPAN/2 for chi in ORDERS])
        initial, mobile = mobile_offsets(trace, start, gates)
        emit(test='preparation', mesh=mesh, radius=1., intake=trace.intake, head=trace.head,
             start_depth=start, initial_mobile_offset=initial, trial_anchor=OFFSET,
             anchor_stiffness=1., registry_drag=REGISTRY_DRAG, well_strength=REGISTRY_K,
             field_sha256=trace.hashes)
        for level, w in zip(ORDERS, gates, strict=True):
            chi = trace.sample(trace.chi, w-SPAN/2)[0]
            np.testing.assert_allclose(chi, level, rtol=0, atol=1e-13)
            elapsed = trace.transit(start, w, 8)
            for law in ('central','separable'):
                rows = [held_pair(law, chi, OFFSET, h) for h in STEPS]
                side_work, row_work = mixed_work(chi, law)
                emit(test='matched_steps', mesh=mesh, law=law, local_order=chi,
                     upstream_depth=w-SPAN/2, junction_depth=w, elapsed=elapsed,
                     registry=OFFSET, reads=rows, mixed_side_work=side_work,
                     mixed_registry_work=row_work)
                for reading in rows:
                    np.testing.assert_allclose(reading['offset_force_minus'], -reading['offset_force_plus'],
                                               rtol=0, atol=1e-20)
                    np.testing.assert_allclose(reading['holding_change_plus'], reading['holding_change_minus'],
                                               rtol=0, atol=1e-17)
                np.testing.assert_allclose(side_work, row_work, rtol=3e-6, atol=1e-19)
                if law == 'separable':
                    np.testing.assert_allclose([r['ratio'] for r in rows], 1., rtol=0, atol=2e-15)
                else:
                    predicted = -registry_force(OFFSET, chi, REGISTRY_K)/OFFSET
                    limit = (4*rows[1]['offset_slope']-rows[0]['offset_slope'])/3
                    np.testing.assert_allclose(limit, predicted, rtol=0, atol=1e-14)
                    np.testing.assert_allclose(rows[0]['holding_change_plus'], 4*rows[1]['holding_change_plus'],
                                               rtol=3e-6, atol=1e-17)
            for kind, q in zip(('elastic_anchor','released'), mobile[w], strict=True):
                reading = held_pair('central', chi, q, STEPS[0])
                emit(test='mobile_pair', mesh=mesh, support=kind, local_order=chi, elapsed=elapsed,
                     actual_offset=q, **reading)


if __name__ == '__main__':
    main()
