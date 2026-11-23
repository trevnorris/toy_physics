"""Constitutive force and signed joint work on the continuously matched held R4 crystal."""
import argparse
import json
from math import fsum
from pathlib import Path
import sys

import numpy as np

from mouth_r4_ring_work import coherent_energy, digest, joint_probe, measure_ring, read


def run(input_guard, output, build):
    from mouth_continuous_order import ContinuousTrial
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_rim_pinned_comparison import validate_guard
    from mouth_continuous_evidence import held_receipt

    receipt = held_receipt(input_guard)
    command = read(input_guard/'invocation.json')['command']
    runner, original, _, saved = map(Path,command[command.index('-u')+1:])
    record = read(saved/'record.json')
    parent = original/'fine'
    validate_guard(parent)
    state = Path(read(parent/'invocation.json')['command'][5])
    known = parent/read(parent/'endpoint.json')['file']
    controls = parent/'endpoint_controls.npz'
    model = runner.with_name('mouth_continuous_order.py')
    deps = original/'pin'/'software/stage1_solver/tools'
    for path,key in ((known,'fine'),(state,'state'),(controls,'controls'),(original/'sources.json','saved_original_pin'),
                     (runner,'code'),(model,'model')):
        if digest(path)!=record[key]:
            raise ValueError(f'Mismatched source or saved fine preparation: {path}')
    if set(record['solutions'])!={'rest','wound'}:
        raise ValueError('Both stationary preparations are required')
    with np.load(known) as f:
        x,q = f['x'],f['q']
    trial = ContinuousTrial(Front(state,GRIDS['double'],1),state,x,**record['solid'],precondition=False)
    common = Path(__file__).with_name('mouth_r4_ring_work.py')
    sources = (Path(__file__),common,common.with_suffix('.cpp'),deps/'mouth_free_shape_solid.cpp',
               original/'sources.json',known,controls,state,runner,model,saved/'record.json',receipt)
    result = dict(command=sys.argv, inputs={str(p):digest(p) for p in sources},
                  sign='negative derivative with respect to physical outward ring displacement',
                  functional='absolute merged order; reference-subtracted coherent vortex; frozen external crystal and material labels',
                  cases={})
    for name in ('rest','wound'):
        case = record['solutions'][name]
        file = saved/(name+'.npz')
        if digest(file)!=case['sha256']:
            raise ValueError(f'Unmatched {name} field')
        result['inputs'][str(file)] = digest(file)
        with np.load(file) as f:
            if not (np.array_equal(x,f['x']) and np.array_equal(q,f['q'])):
                raise ValueError('Preparation did not hold the same wall and lattice')
            u,load = f['field'],f['load']
        if np.max(abs(u[2]-u[2,trial.first][trial.number]))>2.e-13:
            raise ValueError('Order differs on coincident banks')
        trial.model.material_load = load
        ring = measure_ring(trial.b,q,u[2],deps,build)
        terms = coherent_energy(trial.model,u)
        terms.update(chi_W_inner=ring['full_domain_chi_W'],chi_W_exterior=float(trial.ext@u[2]))
        joint,field_forces,joint_force = joint_probe(trial.b,q,u,external_load=trial.ext)
        checks = dict(ring_matches_saved=abs(ring['full_domain_radial_force']-case['solid']['last_free_ring_angular_sum'][0])<2.e-10,
                      energy_matches_saved=abs(fsum(terms.values())-case['energy'])<2.e-8,
                      signed_joint_force=abs(joint_force-ring['full_domain_radial_force'])<2.e-5,
                      direct_field_forces=all(f==0. for f in field_forces.values()),
                      stationary=case['meta']['field']['stationarity']<2.e-7)
        result['cases'][name] = dict(ring=ring, energy_terms=terms, total_energy=fsum(terms.values()),
                    signed_joint=joint, full_joint_probe_force=joint_force, signed_direct_field_forces=field_forces,
                    checks={k: bool(v) for k,v in checks.items()})
        (build/(name+'.json')).write_text(json.dumps(result['cases'][name],allow_nan=False)+'\n')
        print(name, json.dumps(result['cases'][name],allow_nan=False),flush=True)
    result['accepted']=all(all(v['checks'].values()) and all(v['ring']['checks'].values())
                           for v in result['cases'].values())
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['accepted']:
        raise ArithmeticError('Continuous work readback checks failed. See the saved numbers.')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ('input_guard','output','build'):
        p.add_argument(arg,type=Path)
    args = p.parse_args()
    command = read(args.input_guard/'invocation.json')['command']
    runner, original = map(Path,command[command.index('-u')+1:command.index('-u')+3])
    for key,sha in read(original/'sources.json').items():
        if digest(original/'pin'/key)!=sha:
            raise ValueError(f'Changed original input: {key}')
    sys.path[:0] = [str(runner.parent),str(original/'pin'/'software/stage1_solver/tools')]
    run(args.input_guard,args.output,args.build)
