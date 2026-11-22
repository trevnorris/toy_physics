"""Match the flat-face bath densities and coherent phase to a spherical drain.

The planar path length adjusts. In a separate ensemble both individual bath
heads are fixed as well; their coherent quantum stresses are measured rather
than treating the finite planes as uniform reservoirs. The actuating body force
holds the planar front at the same location throughout a given comparison.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq

from mouth_free_front_wave_head import bath_head, emit
from mouth_held_order_current import measure


RECORD = (Path(__file__).resolve().parents[3] /
          'research/pde_ledger_v3/_measurements/mouth_coupled_conversion.json')


def curves():
    archive = json.loads(RECORD.read_text())
    for name in ('curved_pressure', 'curved_matched_spatial_phase'):
        record = archive['runs'][name]
        raw = record['literal_stdout']
        if (hashlib.sha256(raw.encode()).hexdigest() != record['stdout_sha256']
                or not record['limits_verified'] or record['guard_exit'] or record['worker_exit']):
            raise ValueError('Unverified curved reference')
        controls = [json.loads(line) for line in raw.splitlines()]
        row, = [item for item in controls if item['test']=='free_spherical_front'
                and item['trial_action_per_inner_area']==.001
                and (name!='curved_pressure' or item['brane_pressure_excess']==.0045)]
        yield name, row, record['stdout_sha256']


def matched(ref_name, curved, archived_sha, prescribe_heads):
    pressure = curved['brane_pressure_excess']
    n_bulk = (1-5*pressure)**.2
    j_inner = curved['current_per_inner_area']
    j_outer = j_inner*(curved['inner']/curved['outer'])**3
    curve_heads = np.array([bath_head(1.,j_inner), bath_head(n_bulk,j_outer)])
    phase = curved['phase_bulk_minus_brane']
    cached = {}
    def phase_error(length):
        row, sol, _ = measure(-pressure,.5,1.,length/2,.105,5e-10,return_fields=True,
                              endpoint_heads=curve_heads if prescribe_heads else None)
        x = np.linspace(-length/2, length/2, 8001)
        n = sol.sol(x)[0]**2
        actual_phase = sol.p[0]*simpson(1/n,x=x)
        cached[length] = row, sol, x, actual_phase
        return actual_phase-phase

    scale = phase/(3.*pressure)
    length = brentq(phase_error,.82*scale,1.17*scale,xtol=3.e-12)
    row, sol, x, actual_phase = cached[length]
    ends = x[[0,-1]]
    f, df, q, dq, heads, _ = sol.sol(ends)
    actual_ddf = sol.sol(ends,1)[1]
    j = sol.p[0]
    n = f*f
    u = j/n
    p_static = n**5/5
    convective = j*u
    quantum = .5*(df*df-f*actual_ddf)
    order = dq*dq - (q*(1.-q))**2
    normal = p_static+convective+quantum+order
    head_loss = heads[0]-heads[-1]
    order_on_grid = sol.sol(x)[2]
    emit(test='same_pressures_and_phase', reference=ref_name, reference_sha256=archived_sha,
         also_match_both_heads=prescribe_heads, pressure_brane_minus_bulk=pressure,
         bath_densities=n.tolist(), phase_zero_at_ordered_bath=0., target_phase=phase,
         flat_phase=actual_phase, curved_length=curved['outer']-curved['inner'], flat_length=length,
         curved_heads=curve_heads.tolist(), flat_heads=heads.tolist(),
         curved_port_current_per_area=[j_inner,j_outer], flat_current_per_area=j,
         curved_current_per_front_area=curved['current_at_front'],
         flat_over_curved_front=j/curved['current_at_front'],
         static_pressure_at_ports=p_static.tolist(), flat_velocity_at_ports=u.tolist(),
         advected_normal_flux=convective.tolist(), quantum_normal_stress=quantum.tolist(),
         order_normal_stress=order.tolist(), total_normal_flux=normal.tolist(),
         fixture_on_material=row['fixture_on_medium'], normal_flux_jump=normal[-1]-normal[0],
         head_per_constituent=head_loss, reservoir_power_per_front_area=j*head_loss,
         conversion_loss_per_front_area=row['conversion_power'], rms_residual=row['rms_residual_max'],
         max_order_rise=max(np.diff(order_on_grid)))
    if (max(abs(n-[1.,n_bulk]))>1e-12 or abs(actual_phase-phase)>2e-12
            or max(np.diff(order_on_grid))>2e-9):
        raise ArithmeticError('Phase, density, or order profile failed')
    if prescribe_heads and max(abs(heads-curve_heads))>1e-12:
        raise ArithmeticError('Imposed generalized heads failed')
    if abs(normal[-1]-normal[0]-row['fixture_on_medium'])>2e-8:
        raise ArithmeticError('External traction and port momentum discrepancy')
    if abs(row['conversion_power']-j*head_loss)>3e-11:
        raise ArithmeticError('External power and conversion discrepancy')


def main():
    for name, curved, sha in curves():
        for prescribe_heads in (False, True):
            matched(name,curved,sha,prescribe_heads)


if __name__=='__main__':
    main()
