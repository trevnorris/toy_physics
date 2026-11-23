"""Read a hashed accepted R4 shape with its own bonded solid and both reservoirs.

Hold the wound displacement while solving rest at one winding and driven zero
winding. Side and disk tests stay in the unperturbed terminal-band coordinates.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import set_front
from mouth_disk_boundary_work import cap_node_work, crystal_cap, rim_order_work
from mouth_disk_pillbox import Pillbox
from mouth_joint_polar_shape import Trial
from mouth_joined_front_audit import Front
from mouth_joined_front_trace import face
from mouth_polar_joined import PolarCoupled
from mouth_supercore_joined_front import check_guard, emit, file_sha


def accepted(state,path):
    check_guard(state)
    # The solver publishes a separate hard link at each accepted shape. In an
    # ongoing solve, its signed record authenticates that immutable endpoint.
    # Decode exactly the bytes hashed, even if the caller chose the live name.
    data=path.read_bytes()
    sha=hashlib.sha256(data).hexdigest()
    limits=json.loads((path.parent/'limit-validation.json').read_bytes())
    if not limits['verified'] or int(limits['actual']['memory.swap.max'])!=0:
        raise ValueError('Unverified limit or swap on the producing solve')
    rows=[json.loads(s) for s in (path.parent/'stdout').read_text().splitlines()
          if s.startswith('{') and s.endswith('}')]
    prep=next(r for r in rows if r.get('kind')=='preparation')
    rec=next((r for r in rows if r.get('kind')=='accepted_shape' and r['sha256']==sha),None)
    if rec is None or prep['winding']!=1 or prep['state']!=file_sha(state):
        raise ValueError('Missing accepted wound checkpoint on supplied source')
    with np.load(io.BytesIO(data)) as f:
        x,d,u=f['x'],f['displacement'],f['field']
        if (int(f['iteration'])!=rec['iteration'] or f['state_sha'].item()!=prep['state'] or
                not np.array_equal(x,rec['coordinates'])):
            raise ValueError('Accepted checkpoint and record disagree')
    return prep,rec,x,d,u


def run(args):
    prep,rec,x,d,u=accepted(args.state,args.checkpoint)
    base=Front(args.state,tuple(prep['field_grid']),1)
    t=Trial(base,args.state,x,(d.shape[1]-1)//3,d.shape[0]-1,True)
    m,b=t.model,t.body
    q=d[b.mask]*b.scale
    m.material_load=b.get_material_load(q,m)
    set_front(b,m,u)
    b.energy(q)
    err=max(abs(m.norm(m.vars(m.evaluate(u)[1]))))
    material=max(abs(b.native[b.mask]))
    if err>=2e-7 or material>=3e-5:raise ArithmeticError(f'Checkpoint residuals: {err}, {material}')
    emit(kind='source',state_sha256=prep['state'],checkpoint_sha256=rec['sha256'],
         iteration=rec['iteration'],grid=base.grid,coordinates=x.tolist(),depth=m.shape.length,
         exit_radius=m.shape.exit,field_residual=err,material_residual=material)
    rest,rs=m.solve(rest=True,initial=u)
    zero=PolarCoupled(m.shape,base.c,0,base.grid,polar_shift=x[5],terminal_profile=m.terminal_profile)
    zero.material_load=b.get_material_load(q,zero)
    plain,zs=zero.solve(initial=u)
    drive_read=m.measure(u)
    rest_read=m.measure(rest)
    zero_read=zero.measure(plain)
    if any(v['port_error']>3e-5 for v in (drive_read,rest_read,zero_read)):
        raise ArithmeticError('Reservoir current did not close')
    np.savez_compressed(args.save,wound=u,rest=rest,unwound=plain,checkpoint_sha=rec['sha256'])
    emit(kind='matched_flow',sha256=file_sha(args.save),wound=drive_read,rest=rest_read,
         zero=zero_read,retained_fraction=drive_read['feed_current']/zero_read['feed_current'],
         rest_stationarity=rs['stationarity'],zero_stationarity=zs['stationarity'])
    pill=Pillbox(m,base.c,b,d,[u,rest])
    emit(kind='shared_cap_variation',banks=['duct','over_disk','outside_rim'],**cap_node_work(pill))
    for side in (-1,1):
        wt,n,ma,j,T=face(m,base.c,u,side,12)
        T-=face(m,base.c,rest,side,12)[-1]
        positive=np.maximum(j,0)
        dense=(n>.5)&(ma<1)
        solid=crystal_cap(base,pill,False,10) if side<0 else 0.
        emit(kind='disk',side=side,fluid_and_order=wt@T,crystal=solid,total=wt@T+solid,
             signed_current=wt@j,dense_subsonic_share=wt[dense]@positive[dense]/(wt@positive))
    for h in (.12*2/7,.12*4/7,.12):
        emit(kind='constrained_rim',h=h,gedot_chi_z_in=rim_order_work(pill,h,1),
             gedot_chi_z_out=rim_order_work(pill,h,2))
        for compact in (False,True):emit(**pill.audit(h,compact,8))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('state','checkpoint','save'):parser.add_argument(key,type=Path)
    run(parser.parse_args())
