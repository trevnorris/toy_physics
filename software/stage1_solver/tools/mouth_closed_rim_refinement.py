"""Repeat the same four R4 closed work windows on a nested carrier and front grid."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_endpoint import wall_reaction
from mouth_closed_rim_surface_review import Reader
from mouth_closed_rim_work import Closed
from mouth_disk_mesh_rematch import remap
from mouth_fixed_wall_angular import AngularBonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_joined import PolarCoupled
from mouth_supercore_joined_front import check_guard, emit, file_sha


def run(a):
    for p in (a.state,a.endpoint,a.control):
        check_guard(p)
    rows=lambda p:[json.loads(l) for l in (p.parent/'stdout').read_text().splitlines()]
    facts=rows(a.control)
    if (facts[0]['controls_sha']!=file_sha(a.control) or facts[0]['state']!=file_sha(a.state)
            or facts[0]['endpoint']!=file_sha(a.endpoint)
            or rows(a.endpoint)[-1]['sha256']!=file_sha(a.endpoint)):
        raise ValueError('Unverified relaxed state or matched source control')
    with np.load(a.endpoint) as f:
        x,d,s=f['x'],f['displacement'],float(f['internal_amplitude'])
    with np.load(a.control) as f:
        if f['checkpoint_sha'].item()!=file_sha(a.endpoint):raise ValueError('Held shape mismatch')
        old=[f[k] for k in ('wound','rest')]
    base=Front(a.state,(13,10,14),1)
    n,k=(d.shape[1]-1)//3,d.shape[0]-1
    tr=Trial(base,a.state,x,n,k,False)
    b=AngularBonded(a.state,n,k,shape=tr.shape,epsilon=x[5],full_depth=True)
    b.internal_amplitude=s
    original=tr.model
    doubled=tuple(2*i for i in base.grid)
    m=PolarCoupled(tr.shape,base.c,1,doubled,polar_shift=x[5],terminal_profile=original.terminal_profile)
    if any(len(fin)!=2*len(old)-1 or max(abs(fin[::2]-old))>2e-13 for fin,old in
           zip((*m.z,*m.rr,m.a),(*original.z,*original.rr,original.a))):
        raise ValueError('Physical grids are not nested')
    m.material_load=b.get_material_load(d[b.mask]*b.scale,m)
    pair=[]
    emit(kind='source',state=file_sha(a.state),endpoint=file_sha(a.endpoint),
         control=file_sha(a.control),old_grid=base.grid,new_grid=doubled,
         held_internal_amplitude=s,solid_cells=[3*n,k])
    for is_rest,u in enumerate(old):
        init=remap(original,m,u)
        updated,result=m.solve(rest=bool(is_rest),initial=init)
        observed=m.measure(updated)
        emit(kind='field_equilibrium',rest=bool(is_rest),**result,**observed)
        if observed['port_error']>3e-5 or result['stationarity']>=2e-7:
            raise ArithmeticError('Nested mesh did not match both reservoirs')
        pair.append(updated)
    np.savez_compressed(a.save,wound=pair[0],rest=pair[1],checkpoint_sha=file_sha(a.endpoint),grid=doubled)
    emit(kind='matched_states',sha256=file_sha(a.save))
    study=Closed(m,b,d,pair,window_stride=2)
    reader=Reader(m,base.c,*pair)
    for ix in (2,3):
        ax,_,_,geo=study.lifts(ix,4)
        slit=wall_reaction(reader,ax,m.shape.length-geo['lower_half_thickness'][-1],10)
        for j in (4,5):
            row=study.audit(ix,j)
            old=next(f for f in facts if f.get('axial_band')==ix and f.get('lateral_band')==j)
            for key in ('lower_half_thickness','upper_half_thickness','lateral_radius_at_disk'):
                if not np.allclose(row[key],old[key],rtol=0.,atol=2e-13):
                    raise ArithmeticError('Physical test moved under refinement')
            emit(**row,held_slit_work_banks_mass_and_axial=slit.tolist(),
                 sharp_slit_total_axial=slit[:,1].sum(),outer_plus_slit_axial=row['outward_momentum']+slit[:,1].sum())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','endpoint','control','save'):
        p.add_argument(name,type=Path)
    run(p.parse_args())
