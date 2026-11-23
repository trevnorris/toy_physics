"""First-order axial force for each moving carrier extension on the held R4 disk.

Integrate the relative carrier/material work on its entire compact support,
including the central aperture and the full radial shoulder of the exterior
bank. Add the common physical disk work and the signed interelement work.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_closed_rim_surface_review import Reader, quadrature
from mouth_disk_pillbox import dot
from mouth_exchange_carrier import mapped_nodes
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_rim_pinned_comparison import digest, validate_guard
from mouth_shared_rim_extension_work import gradient
from mouth_shared_rim_shift import SharedRimChart
from mouth_supercore_joined_front import ORDER_INERTIA, state_from


def smooth_work(disk, chart, slab, side):
    model=disk.model
    view=SimpleNamespace(model=model,sh=disk.sh,reader=SimpleNamespace(vortex=disk.reader.vortex,
                         fields={(1,side):disk.reader.fields[slab,side]}))
    bank=Bank(view,side,None,0.)
    edges=model.rr[side]
    end=max(chart.outer_taper[-1],chart.edges[4][1],chart.edges[8][1]) if side else 1.
    nodes,weights=quadrature(np.unique(np.r_[edges[edges<end],end]),2)
    angles,aw=quadrature(model.a,2)
    angle_grid=np.meshgrid(nodes,angles,indexing='ij')
    x,angle=(p.ravel() for p in angle_grid)
    weights=(weights[:,None]*aw).ravel()
    lo,hi=chart.tents[[0,2] if slab==1 else [2,4]]
    zz=model.z[slab]
    depths,zw=quadrature(np.unique(np.r_[lo,zz[(zz>lo)&(zz<hi)],hi]),2)
    result=np.zeros((2,3))
    for z,wz in zip(depths,zw):
        b=disk.sh.radius(z)[0]
        width=28-b if side else b
        r=x*width+(b if side else 0.)
        volume=4*np.pi*weights*width*r*r*np.sin(angle)*wz
        physical=chart.physical(side,r,z)
        motions=[np.subtract(chart.carrier(side,r,z,n),physical) for n in (4,8)]
        for sign,before,src in zip((1.,-1.),bank.physical(x,angle,z),disk.reader.fields[slab,side]):
            (_,p,c,rho,_,_),amp,ga,F=before
            y=amp**2
            background=F**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20
            pc=dot(p,c)
            dn=rho+ORDER_INERTIA*dot(c,c)
            pj=-.5*ORDER_INERTIA*rho*pc**2/dn
            eps=1e-7
            gc_z=(gradient(disk,slab,side,(r,angle,z+eps),src.values[:,:,:,2])
                  -gradient(disk,slab,side,(r,angle,z-eps),src.values[:,:,:,2]))/(2*eps)
            for i,(eta,er,ez) in enumerate(motions):
                dp=np.column_stack((-er*p[:,2],np.zeros_like(r),-ez*p[:,2]))
                da=np.column_stack((-er*ga[:,2],np.zeros_like(r),-ez*ga[:,2]))
                dc=eta[:,None]*gc_z
                proj_derivative=(-ORDER_INERTIA*rho*pc/dn*(dot(dp,c)+dot(p,dc))
                                 -pj*2*ORDER_INERTIA*dot(c,dc)/dn)
                sectors=[ez*(.5*F**2*dot(ga,ga)+background)+F**2*dot(ga,da),
                         .5*rho*ez*dot(p,p)+rho*dot(p,dp),ez*pj+proj_derivative]
                result[i]-=sign*(np.array(sectors)@volume)
    return result


def run(a):
    for folder in (a.job,a.faces.parent,a.difference.parent):
        validate_guard(folder)
    fields_sha=digest(a.job/'zero.npz')
    previous=json.loads(a.manifest.read_text())
    c=next(c for c in previous['cases'] if c['inputs']['grid']==a.job.name)
    face=json.loads(a.faces.read_text())
    surface=next(s for s in json.loads(a.surface.read_text())['cases'] if s['grid']==a.job.name)
    diff=json.loads(a.difference.read_text())
    m=c['inputs']
    if (fields_sha!=c['zero_state_sha256'] or face['field']!=fields_sha or surface['zero_sha256']!=fields_sha
        or diff['zero_sha256']!=fields_sha or digest(a.axes)!=m['reference_axes_sha256']
        or digest(a.state)!=m['state_sha256'] or digest(a.wall)!=m['wall_sha256']
        or digest(a.original)!=m['field_sha256']):
        raise ValueError('The mesh, frozen state or physical geometry differs')
    disk=Disk(a.state,a.wall,a.original)
    with np.load(a.job/'zero.npz') as f:
        disk.reader=Reader(disk.model,state_from(a.state)[1],f['wound'],f['rest'])
    chart=object.__new__(SharedRimChart)
    chart.shape=disk.sh
    chart.outer_taper=mapped_nodes(13,disk.sh.length,28.,22.)[0][1][4:6]
    chart.tents=np.array(face['heights'])
    chart.window=np.array(disk.windows[0])
    chart.edges={n:(m['duct_bounds'][i],m['receiver_bounds'][i]) for i,n in enumerate((4,8))}
    fdata=next(r for r in face['readings'] if r['gauss']==2)['cases']
    sweeps=np.stack([np.array(fdata[key]['breakdown']).sum(axis=1) for key in ('rim_four','rim_eight')])
    regions=[]
    rel=np.zeros((2,2))
    for slab,side in np.ndindex(2,2):
        smooth=smooth_work(disk,chart,slab+1,side)
        strips=sweeps[:,:,slab,side]
        combined=smooth.sum(axis=1)[:,None]+strips
        rel+=combined
        regions.append(dict(slab=slab+1,outer=side,smooth_sectors=smooth.tolist(),
                            face_work=strips.tolist(),one_sided=combined.tolist()))
    resolved=next(r for r in diff['checks'] if r['gauss']==2)
    err=np.max(abs(np.subtract(rel[0],rel[1])-resolved['totals']))
    if err>1e-9:
        raise ArithmeticError(f'Full-support and difference-only extensions disagree by {err}')
    cap={r['window']:r for r in surface['results'] if r['degree']==2}
    common=cap['near']['common_reaction']
    return dict(grid=a.job.name,zero_sha256=fields_sha,face_sha=digest(a.faces),surface_sha=digest(a.surface),
                rule='Rows are 4 and 8 cells. One-sided columns are negative and positive displacement.',
                sectors=['amplitude and native density','phase','projection'],
                regions=regions,relative_one_sided=rel.tolist(),relative_centered=rel.mean(axis=1).tolist(),
                common_disk=common,near_outward=cap['near']['outward'][1],far_outward=cap['far']['outward'][1],
                near_minus_far_work=common-cap['far']['common_reaction'],predicted_one_sided=(common+rel).tolist(),
                predicted_centered=(common+rel.mean(axis=1)).tolist(),difference_check=err)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('state','wall','original','axes','job','manifest','faces','surface','difference'):
        p.add_argument(arg,type=Path)
    print(json.dumps(run(p.parse_args()),allow_nan=False))
