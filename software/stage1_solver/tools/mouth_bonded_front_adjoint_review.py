"""Check transposed material-to-order deposition against physical solid-angle quadrature."""
from pathlib import Path
import json,hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
from mouth_bonded_material_motion import Bonded
from mouth_bonded_front_transfer_review import Front,check_guard
def readout(model,owner,displacement,field):
    coeff=displacement[owner.mask]*owner.scale
    qdata=displacement.reshape(-1,2)[owner.id]
    motion=np.einsum('qi,qij->qj',owner.basis,qdata)
    grad=np.stack([np.einsum('qi,qij->qj',d,qdata) for d in (owner.dr,owner.dw)],axis=-1)
    r=owner.r+motion[:,0];w=owner.w+motion[:,1]
    B=np.eye(2)+grad
    hoop=r/owner.r*owner.h
    J=owner.j*np.linalg.det(B)*(r/owner.r)**2
    K=B@owner.K@B.transpose((0,2,1))
    W=.5*(np.trace(K,axis1=1,axis2=2)+2*hoop**2-4)-np.log(J)+.25*(J-1)**2
    # Integrate sin(theta) times 3-coordinate Eulerian hats at physical positions
    # by independent Gaussian quadrature in each angular subinterval.
    vals=field[model.ids[1,0]].transpose((0,2,1))
    interp=RegularGridInterpolator((model.rr[0],model.z[1]),vals)
    at=interp((r/model.shape.radius(w)[0],w))
    q,wt=leggauss(6);q=(q+1)/2;wt/=2
    angles=model.a[:-1];da=np.diff(model.a)
    integrals=np.stack([np.sin(angles+v*da)*a*da for v,a in zip(q,wt)])
    aw_right=q@integrals;aw_left=(1-q)@integrals
    chi_val=at[:,:-1]@aw_left+at[:,1:]@aw_right
    return (owner.measure*W)@chi_val


def main(args):
    from mouth_bonded_front_adjoint import get_material_load
    for path in (args.state,args.flow,args.motion):check_guard(path)
    base=Front(args.state,(13,10,14),1)
    model=base.joined(base.sh,False)
    with np.load(args.flow) as f:dc=f['fields'][2]
    with np.load(args.motion) as f:
        motion=f['displacement'];sx=f['radial_fraction'];st=f['depth_fraction']
    owner=Bonded(args.state,(len(sx)-1)//3,len(st)-1,args.flow)
    # Same scalar variation in theta and a separate angularly varying perturbation.
    # Shared cap node IDs and the free-order mask are respected explicitly.
    rng=np.random.default_rng(408)
    trial=np.zeros(model.size)
    trial[model.masks[2]]=rng.normal(scale=.07,size=np.sum(model.masks[2]))
    out=lambda **row:print(json.dumps(row),flush=True)
    for name,field in [('measured',dc),('signed_angular',trial)]:
        def energies(alpha):
            u=alpha*motion
            q=u[owner.mask]*owner.scale
            direct=readout(model,owner,u,field)
            loaded=get_material_load(owner,q,model)@field
            return np.array([direct,loaded])
        v=energies(1.)
        out(kind='value',field=name,independent=v[0],deposited=v[1],error=v[1]-v[0])
        if abs(v[1]-v[0])>2e-11:raise ArithmeticError('Absolute work mismatch')
        for h in (.0004,.0002):
            dif=(energies(1+h)-energies(1-h))/(2*h)
            out(kind='slope',field=name,step=h,independent=dif[0],deposited=dif[1],error=dif[1]-dif[0])
            if abs(dif[1]-dif[0])>4e-8:raise ArithmeticError('Directional transfer mismatch')
    out(kind='provenance',**{k:hashlib.sha256(getattr(args,k).read_bytes()).hexdigest() for k in ('state','flow','motion')})


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('state','flow','motion'):p.add_argument(k,type=Path)
    main(p.parse_args())
