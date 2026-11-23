"""Five whole-mouth virtual displacements at a reciprocally matched bonded R4 field."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import Coupled, get_material_load
from mouth_bonded_material_motion import Bonded
from mouth_internal_collar import CollarPair
from mouth_joined_front_audit import Front
from mouth_poly_shallow_release import CUT, Order, exterior
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit
from mouth_supercore_terminal_order import terminal


def run(a):
    for path in (a.state,a.joint):check_guard(path)
    rows=list(map(json.loads,(a.joint.parent/'stdout').read_text().splitlines()))
    digest=hashlib.sha256(a.state.read_bytes()).hexdigest()
    if rows[0]['source_sha256']!=digest or rows[-1]['sha256']!=hashlib.sha256(a.joint.read_bytes()).hexdigest():
        raise ValueError('Mismatched held equilibrium and base')
    base=Front(a.state,tuple(rows[0]['grid']),1)
    with np.load(a.joint) as z:
        u=z['field'];d=z['displacement']
    n=base.order_mesh;t=np.linspace(0,1,2*n+1);x=np.linspace(0,1,n+1)
    def energies(v):
        params=base.sh.params+np.array([1,.25,.25,.25])*v[:4]
        sh=Meridian(4.,params,extra=base.sh.extra+v[4]*.25)
        pair=CollarPair(sh,base.co,base.na,base.nt)
        out=pair.evaluate(base.labels,(t,x,base.c),((1-np.cos(np.pi*t))/2,.1*x+.9*x*x,base.outer))['outer_source']
        source=lambda wa,xo:out
        ext=exterior(sh,n,CUT,side_order_work=source)
        profile=terminal(sh,n,CUT,source)
        model=Coupled(sh,base.c,1,base.grid,terminal_profile=profile,precondition=False)
        body=Bonded(a.state,(d.shape[1]-1)//3,d.shape[0]-1,shape=sh,full_depth=True)
        q=d[body.mask]*body.scale
        model.material_load=get_material_load(body,q,model)
        # All inner crystal and bonded core strain work enters through this
        # body integral exactly once. Exterior.static_work is the other half.
        return np.array([Order(sh,n).energy_gradient(base.c)[0],ext['energy']+ext['static_work'],
                         body.energy(q)[0],model.evaluate(u)[0]])
    emit(kind='base',source_sha256=digest,energies=energies(np.zeros(5)).tolist(),
         coordinates=['log_depth','.25 log_exit','.25 middle','.25 upper_minus_lower','.25 waist'])
    for h in (.0003,.00015):
        drive=[]
        for i,e in enumerate(np.eye(5)):
            lo,hi=energies(-h*e),energies(h*e)
            f=(lo-hi)/(2*h);drive.append(f.sum())
            emit(kind='shape_work',coordinate=i,h=h,force=float(f.sum()),terms=f.tolist(),
                 values_minus=lo.tolist(),values_plus=hi.tolist())
        norm=np.linalg.norm(drive)
        if norm>0:
            dE=(energies(h*np.asarray(drive)/norm)-energies(np.zeros(5))).sum()
            emit(kind='joint_descent',h=h,direction=(np.array(drive)/norm).tolist(),energy_change=dE)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('joint',type=Path)
    run(p.parse_args())
