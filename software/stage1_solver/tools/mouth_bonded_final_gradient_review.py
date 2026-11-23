"""Reconstruct the two returned gradients of the work-conjugate held shape."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from mouth_bonded_front_adjoint import get_material_load, set_front
from mouth_bonded_front_adjoint_review import readout
from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_supercore_front import check_guard


def audit(a):
    for path in (a.state,a.flow,a.coupled):check_guard(path)
    base=Front(a.state,(13,10,14),1)
    model=base.joined(base.sh,True)
    # A common point-free source is used below for the constituent-order term.
    model.D[:,-1]=0.
    with np.load(a.coupled) as f:
        disp=f['displacement'];U=f['field'];nr=len(f['radial_fraction'])-1;nt=len(f['depth_fraction'])-1
    owner=Bonded(a.state,nr//3,nt,a.flow)
    q=disp[owner.mask]*owner.scale
    def output(**kw):print(json.dumps(kw),flush=True)
    def gradients(v,u):
        set_front(owner,model,u)
        _,gq=owner.energy(v)
        native_q=owner.native[owner.mask].copy()
        _,gu=model.evaluate(u)
        gu[2]+=get_material_load(owner,v,model)
        native_front=model.vars(gu)
        return gq,gu,native_q,max(abs(model.norm(native_front)))
    def energy(v,u):
        owner.front=None
        default=owner.energy(v)[0]
        motion=np.zeros_like(disp);motion[owner.mask]=v/owner.scale
        return default+model.evaluate(u)[0]+readout(model,owner,motion,u[2])
    gq,gu,native,gf=gradients(q,U)
    output(kind='reconstructed_final_residual',native_material=max(abs(native)),field_metric=gf)
    rng=np.random.default_rng(309)
    # Verify both gradients near the saved state and after a physical off-equilibrium
    # displacement, where their signed values cannot agree merely by being small.
    dq=np.zeros_like(disp);dq[owner.mask]=rng.normal(size=q.size)
    dq[~owner.mask]=0.;dq=dq[owner.mask]*owner.scale
    du=np.zeros_like(U);du[2,model.masks[2]]=rng.normal(size=model.masks[2].sum())
    du/=np.max(abs(du))
    for scale in (0.,4e-4):
        v=q+scale*dq
        u=U+scale*du
        gq,gu,*_=gradients(v,u)
        for axis in ('material','order'):
            dv=dq if axis=='material' else np.zeros_like(q)
            dU=du if axis=='order' else np.zeros_like(U)
            exact=gq@dv+np.sum(gu*dU)
            for eps in (2e-6,1e-6):
                fd=(energy(v+eps*dv,u+eps*dU)-energy(v-eps*dv,u-eps*dU))/(2*eps)
                output(kind='gradient_check',perturbation=scale,axis=axis,step=eps,
                       analytic=exact,direct=fd,error=fd-exact)
                if abs(fd-exact)>2e-5*max(1,abs(exact)):raise ArithmeticError('Signed total material/front work mismatch')
    output(kind='sources',**{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ('state','flow','coupled')})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('state','flow','coupled'):p.add_argument(k,type=Path)
    audit(p.parse_args())
