"""Rematch the driven and resting fixed R4 rim to one boundary-aligned crystal load.

Use the archived refined carrier axes and the same external momentum enclosures.
The direct disk displacement is a separate physical variation.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_closed_rim_endpoint import wall_reaction
from mouth_closed_rim_surface_review import Reader
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_rings_fluid import GRIDS, RimCarrier, SameWindows, read_source
from mouth_supercore_joined_front import check_guard, emit, file_sha


def checked_pair(path):
    check_guard(path)
    rows=[json.loads(line) for line in (path.parent/'stdout').read_text().splitlines()]
    records=[row for row in rows if row['kind']=='fields']
    if len(records)!=1 or records[0]['sha256']!=file_sha(path):
        raise ValueError('Field pair does not match the guarded measurement')
    return records[0]['sha256']


def material_load(model, archive, checkpoint):
    check_guard(archive)
    with np.load(archive,allow_pickle=False) as data:
        if (data['checkpoint_sha'].item()!=file_sha(checkpoint)
                or not np.array_equal(data['ax'],model.rr[0])
                or not np.array_equal(data['z'],model.z[1])):
            raise ValueError('Material load has different nodes or was generated on another wall')
        plane=data['plane']
    if plane.shape!=(len(model.rr[0]),len(model.z[1])) or not np.all(np.isfinite(plane)):
        raise ValueError('Invalid material load moments')
    weights=plane[:,None,:]*angular_weights(model.a)[None,:,None]
    return np.bincount(model.ids[1,0].ravel(),weights=weights.ravel(),minlength=model.size)


def from_archive(state, checkpoint, previous, precondition):
    x,d,_,s=read_source(state,checkpoint)
    source_hash=checked_pair(previous)
    base=Front(state,GRIDS['native'],1)
    tr=Trial(base,state,x,11,32,False)
    with np.load(previous,allow_pickle=False) as data:
        level=str(data['level'])
        if (level not in ('double','local') or data['checkpoint_sha'].item()!=file_sha(checkpoint)
                or not np.array_equal(data['grid'],GRIDS[level])):
            raise ValueError('Expected the completed double or rim-refined rejected-wall state')
        axes=tuple(data[key] for key in ('inside_x','outside_x','feed_z','duct_z','receiver_z'))
        old=[data[key] for key in ('wound','rest')]
        angles=data['angles']
    m=RimCarrier(tr.shape,base.c,1,GRIDS[level],polar_shift=x[5],axes=(axes[:2],axes[2:]),
                 terminal_profile=tr.model.terminal_profile,precondition=precondition)
    if not np.array_equal(m.a,angles) or any(v.shape!=(3,m.size) for v in old):
        raise ValueError('Field topology or angular coordinates changed')
    return base,tr,m,old,level,source_hash,d,s


def run(a):
    base,tr,m,old,level,sha,d,s=from_archive(a.state,a.checkpoint,a.previous,True)
    m.material_load=material_load(m,a.load,a.checkpoint)
    common=dict(checkpoint_sha=file_sha(a.checkpoint),prior_fields_sha=sha,load_sha=file_sha(a.load),
                grid=GRIDS[level],level=level,inside_x=m.rr[0],outside_x=m.rr[1],
                feed_z=m.z[0],duct_z=m.z[1],receiver_z=m.z[2],angles=m.a)
    emit(kind='input',checkpoint=common['checkpoint_sha'],source=file_sha(a.state),level=level,grid=GRIDS[level],
         prior_fields_sha=sha,load_sha=common['load_sha'],shape=[m.shape.length,m.shape.exit],
         internal_amplitude=s,nodes=m.size,free_order_nodes=int(m.masks[2].sum()))
    pair=[]
    for n,rest in enumerate((False,True)):
        init=old[n]
        initial_residual=np.max(abs(m.norm(m.vars(m.evaluate(init)[1]))))
        updated,info=m.solve(rest=rest,initial=init)
        measure=m.measure(updated)
        if not (np.isfinite(updated).all() and measure['port_error']<3e-5 and info['stationarity']<2e-7):
            raise ArithmeticError('Unresolved corrected field, phase, or current')
        diff=np.max(abs(updated-init),axis=1)
        pair.append(updated)
        emit(kind='matched',rest=rest,initial_corrected_residual=initial_residual,change_max=diff.tolist(),**info,**measure)
        if not rest:
            corr=updated[2,m.ids[1,0]]
            out=a.folder/'solid_load.npz'
            np.savez_compressed(out,duct_dc=corr,**common)
            emit(kind='solid_load',sha256=file_sha(out))
    target=a.folder/'fields.npz'
    np.savez_compressed(target,wound=pair[0],rest=pair[1],duct_dc=corr,**common)
    emit(kind='fields',sha256=file_sha(target))
    m.solvers=[]
    # No old-cubature crystal pin force enters this calculation. The material
    # stress difference vanishes at the taper because both order fields are fixed
    # there; the corrected loading was included in each field's rematch.
    study=SameWindows(m,None,d,pair,reference=tr.model,include_constraints=False)
    reader=Reader(m,base.c,*pair)
    for k in (2,3):
        taper,_,_,geo=study.lifts(k,4)
        low=m.shape.length-max(geo['lower_half_thickness'])
        banks=wall_reaction(reader,taper,low,10)
        check=wall_reaction(reader,taper,low,6)
        if np.max(abs(banks-check))>2e-6:
            raise ArithmeticError('Sealed-bank quadrature not resolved')
        for j in (4,5):
            row=study.audit(k,j)
            emit(**row,slit_banks_mass_axial=banks.tolist(),slit_force=banks[:,1].sum(),
                 closed_plus_slit=row['outward_momentum']+banks[:,1].sum(),
                 slit_quad_max=float(np.max(abs(banks-check))))
    # Reuse the matched reader, restricting to two bands wholly in the near
    # plateau. Compare the final part of the slit to the earlier archived study.
    span=m.h
    prev=np.zeros((2,2))
    for depth in (span/16,span/8,span/4):
        banks=wall_reaction(reader,lambda z:np.ones_like(z),m.shape.length-depth,10)
        emit(kind='strip',depth_from_exit=depth,bank_mass_axial=banks.tolist(),
             strip_axial=float(sum((banks-prev)[:,1])),cumulative_axial=float(sum(banks[:,1])))
        prev=banks


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','checkpoint','previous','load','folder'):
        p.add_argument(key,type=Path)
    run(p.parse_args())
