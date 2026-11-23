"""Resample the saved end-point displacement on an independently integrated solid patch."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_internal_collar_check import read
from mouth_last_ring_solid import source
from mouth_last_rings_enriched import Patch, rim_patch
from mouth_last_ring_surface_review import FieldChart
from mouth_supercore_joined_front import check_guard, emit, file_sha


def run(state,wall,paired,endpoint,aligned):
    check_guard(paired)
    limit=json.loads((endpoint.parent/'limit-validation.json').read_text())
    guard=json.loads((endpoint.parent/'outcome.json').read_text())
    if not guard['limitsVerified'] or not limit['verified'] or int(limit['actual']['memory.swap.max']):
        raise ValueError('Unverified displacement containment')
    rows=[json.loads(line) for line in (endpoint.parent/'stdout').read_text().splitlines()]
    feedback=[row for row in rows if row.get('kind')=='feedback']
    if not feedback:
        raise ValueError('No completed paired feedback iterate to review')
    s,d,_,_=source(state,wall)
    _,_,inner,_=read(state)
    with np.load(endpoint) as end,np.load(paired) as f:
        if end['checkpoint_sha'].item()!=file_sha(wall) or end['field_sha'].item()!=file_sha(paired) or f['level'].item()!=end['grid'].item():
            raise ValueError('Input is not the reported, fixed geometry and matched field')
        model=FieldChart(s,None,1,f['grid'],precondition=False)
        ids=model.ids[1,0]
        order=np.einsum('raj,a->jr',end['field'][2,ids],angular_weights(f['angles']))
        physical_z=f['duct_z']
        grid=(physical_z/s.length,f['inside_x'],order)
        U=end['patch']
    emit(kind='snapshot',last_iterate=file_sha(endpoint),drive=file_sha(paired),preparation=file_sha(wall),
         primary_exit=guard['exitCode'],primary_sweep=feedback[-1]['sweep'],primary_min=feedback[-1]['min_stretch'],
         primary_force=feedback[-1]['max_force'],primary_field_residual=feedback[-1]['field_residual'])
    if aligned:
        from mouth_rim_compression_feedback import LocalMaterial, OFFSETS, collect
        for n in aligned:
            quad=LocalMaterial(state,wall,grid[1],physical_z,n)
            quad.mesh=8
            _,pt=collect(quad,build=True)
            pt.order=[inner,grid]
            p=Patch(pt,8,OFFSETS)
            values=p.sample(U)
            emit(kind='registered_resampling',subdivision=n,**values)
            del pt,p,quad
        return
    for n in (24,40):
        rr,offsets=rim_patch(state,s,d,[inner,grid],n,True)
        p=Patch(rr,8,offsets)
        if U.size!=p.zero.size:
            raise ValueError('Incompatible internal displacement coordinates')
        values=p.sample(U)
        emit(kind='uniform_resampling',subdivision=n,**values)
        del rr,p


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for nm in ('state','wall','fields','endpoint'):parser.add_argument(nm,type=Path)
    parser.add_argument('--registered',nargs='+',type=int)
    a=parser.parse_args()
    run(a.state,a.wall,a.fields,a.endpoint,a.registered)
