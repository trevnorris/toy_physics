"""Measure the six global loads before and after fixed-wall angular equilibration.

Take signed variations on the same two prepared states, transporting the field
coefficients and physical material displacements on the established chart. Use
paired energy differences so the two large static receiver terms cancel.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_fixed_wall_angular import AngularBonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import emit, file_sha, check_guard


def run(state, source, equilibrium):
    prep, row, x, old_d, old_field = accepted(state,source)
    check_guard(equilibrium)
    rows=[json.loads(line) for line in (equilibrium.parent/'stdout').read_text().splitlines()]
    if rows[-1].get('sha256') != file_sha(equilibrium) or rows[-1].get('kind') != 'saved':
        raise ValueError('A completed guarded angular equilibrium is required')
    with np.load(equilibrium, allow_pickle=False) as z:
        if z['input_sha'].item() != file_sha(source) or z['state_sha'].item() != file_sha(state):
            raise ValueError('The paired preparations do not match')
        if not np.array_equal(z['x'], x):
            raise ValueError('Exterior shape changed between fixed-wall states')
        d,field,s = z['displacement'],z['field'],float(z['internal_amplitude'])
    if d.shape != old_d.shape or field.shape != old_field.shape:
        raise ValueError('Field or solid topology changed')
    base = Front(state, tuple(prep['field_grid']), 1)
    n,k = (d.shape[1]-1)//3, d.shape[0]-1
    emit(kind='input', state=file_sha(state), before=file_sha(source), after=file_sha(equilibrium),
         material_cells=[3*n,k], angular_coordinate=s)
    def comparison(y):
        t = Trial(base,state,y,n,k,False)
        m,b = t.model,t.body
        before_q = old_d[b.mask]*b.scale
        m.material_load = b.get_material_load(before_q,m)
        bs = b.energy(before_q)[0]
        bf = m.evaluate(old_field)[0]
        a = AngularBonded(state,n,k,shape=t.shape,full_depth=True,epsilon=y[5])
        a.internal_amplitude = s
        q = d[a.mask]*a.scale
        full_solid = a.energy(q)[0]
        m.material_load = a.get_material_load(q,m)
        full_field = m.evaluate(field)[0]
        ref_solid = b.energy(d[b.mask]*b.scale)[0]
        return np.array([t.static+bs+bf, (ref_solid-bs)+(full_field-bf), full_solid-ref_solid])
    center = comparison(x)
    emit(kind='center', original_energy=center[0], meridional_and_field_delta=center[1],
         internal_mode_delta=center[2], total_drop=-sum(center[1:]))
    if abs(center[0]-row['energy']) > 3e-8:
        raise ArithmeticError('Starting energy differs from the verified preparation')
    last = None
    for h in (1e-4,5e-5,2.5e-5):
        loads = np.array([(comparison(x-h*axis)-comparison(x+h*axis))/(2*h)
                          for axis in np.eye(6)])
        before,delta = loads[:,0],loads[:,1]+loads[:,2]
        after = before+delta
        peak = int(np.argmax(abs(after)))
        if h==1.e-4 and max(abs(before-row['global_forces'])) > 2e-5:
            raise ArithmeticError('Stored original six signed forces not reproduced')
        emit(kind='six_coordinate_work', h=h, before=before.tolist(), after=after.tolist(),
             meridional_and_field_change=loads[:,1].tolist(), internal_mode_change=loads[:,2].tolist(),
             paired_change=delta.tolist(), peak_coordinate=peak,
             peak_before=float(max(abs(before))),peak_after=float(max(abs(after))),
             peak_relieved=float(max(abs(before))-max(abs(after))),
             paired_step_difference=None if last is None else abs(last-delta).tolist())
        last=delta


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','source','equilibrium'):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    run(a.state,a.source,a.equilibrium)
