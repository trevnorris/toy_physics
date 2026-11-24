"""Collect the archived exact-melt and determinant-gauge checks on both R4 endpoints."""
import argparse
import json
from pathlib import Path

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_stage2_outcome import verify_guard


def assemble(scratch):
    checks = {
        'linear_reextension':('r4_slaved_melt01','static_guard','static.json'),
        'independent_harmonic_variation':('r4_slaved_harmonic_review01','guard','readout.json'),
        'independent_barrier_variation':('r4_slaved_harmonic_review01','barrier_cell_guard','barrier_cell.json'),
        'barrier':('r4_slaved_melt_barrier01','readout_guard','measured.json'),
        'order_source':('r4_slaved_melt01','trace_guard','trace.json'),
    }
    result = dict(radius=4,winding=1,starting_commit='12e02d90',experiments={})
    for name,(folder,guard,file) in checks.items():
        path=scratch/folder
        verify_guard(path/guard,0)
        doc=read(path/file)
        result['experiments'][name]=dict(command=read(path/guard/'invocation.json')['command'][2:],
                                         sha256=digest(path/file), data=doc)
    data=result['experiments']
    for case in ('primary','secondary'):
        state=data['linear_reextension']['data']['cases'][case]['state_sha256']
        barrier=data['barrier']['data']['cases'][case]
        review=data['independent_harmonic_variation']['data'][case]
        source=data['order_source']['data'][case]
        if (state!=barrier['sha256'] or state!=review['saved_state_sha256']
                or state!=source['source_sha256'] or barrier['fixed_trace_state_sha256']!=source['extension_sha256']):
            raise ValueError('Readbacks refer to different endpoints')
        if not barrier['inner']['stationary'] or not barrier['outer_probe']['inner']['stationary']:
            raise ValueError('The nonlinear inner solve did not converge')
        if barrier['outer_probe']['alpha'] < 1.e-4 or barrier['outer_probe']['reduced_energy_change'] >= 0:
            raise ValueError('The probed active step did not pass')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    record=assemble(args.scratch)
    with args.output.open('x') as fp:
        json.dump(record,fp,indent=2,allow_nan=False)
        fp.write('\n')
    for c,b in record['experiments']['barrier']['data']['cases'].items():
        print(c,b['fixed_trace_limiting_physical_Jinv'],b['outer_probe']['alpha'],
              b['outer_probe']['reduced_energy_change'],flush=True)
