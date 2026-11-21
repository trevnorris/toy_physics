"""Reduce matched rim reactions and the uniform-coverage linear observable."""
import json
from pathlib import Path

import numpy as np

from mouth_order_wave_pulse import emit

BASE = Path('_scratch/s11c')


def checked(name):
    path = BASE/name
    outcome = json.loads((path/'outcome.json').read_text())
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or
        not json.loads((path/'limit-validation.json').read_text())['verified']):
        raise ValueError(f'Invalid guarded run {name}')
    return [json.loads(line) for line in (path/'stdout').read_text().splitlines()]


if __name__=='__main__':
    field_hash = checked('companion_patch_field_01')[-1]['sha256']
    projection = checked('companion_rim_projection_01')
    if projection[0]['sha256'] != field_hash:
        raise ValueError('Projected a different preparation')
    zero = checked('companion_rim_zero_grip_01')
    np.testing.assert_array_equal(zero[-1]['absolute_mean_force'],np.zeros(6))
    for law in ('central','layered'):
        readings = {}
        for variant, suffix in [('local','_01'),('coarse','_coarse_01'),('both','_both_01'),('angle','_angle01')]:
            rows = checked(f'companion_rim_{law}{suffix}')
            head,last = rows[0],rows[-1]
            expected = projection[-1]['sha256'] if variant=='angle' else field_hash
            if (head['flow_profile_sha256'] != expected or last['test']!='patch_result' or
                head['cases']!=['alone','toward8','mirror8','toward16','mirror16','half8']):
                raise ValueError('Input fields or channel ordering changed')
            arr = np.array(rows[-2]['rim_impulse'])/last['elapsed']
            delta = arr-arr[:,:1]
            np.testing.assert_allclose(delta,last['delta_force_by_seed'],rtol=1e-12)
            readings[variant] = delta.mean(axis=0)
        local,coarse,both,angle = [readings[k] for k in ('local','coarse','both','angle')]
        angular_odd = (angle[[1,3]]-angle[[2,4]])/2
        emit(test='force_comparison',law=law,
             local_toward=local[[1,3]].tolist(),ratio16_over_8=local[3]/local[1],
             both_planes_gain=(both[[1,3]]/local[[1,3]]).tolist(),
             coarse_over_fine=(coarse[[1,3]]/local[[1,3]]).tolist(),
             half_linearity=2*local[5]/local[1],zero_grip=zero[-1]['delta_mean_force'][1],
             projected_odd_force_per_uniform_patch=angular_odd.tolist(),
             projected_odd_ratio=angular_odd[1]/angular_odd[0])
