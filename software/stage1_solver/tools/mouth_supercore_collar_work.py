"""Independent signed work for an interior redistribution of material in the radius-4 collar."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legder, leggauss, legval

from mouth_finite_bore_support import graph_values
from mouth_polycrystal_material import stable_energy
from mouth_sealed_shape_relax import Meridian


def bump(x, lo, hi):
    u = np.clip((x-lo)/(hi-lo),0.,1.)
    return 16*u**2*(1-u)**2, 32*u*(1-u)*(1-2*u)/(hi-lo)


def evaluate(shape, chi, co, degree):
    # Integrate in the unloaded material chart. Both the collar seam and its
    # core edge remain fixed along with the field at every spatial point.
    nodes, weights=leggauss(degree)
    nodes=(1+nodes)/2
    weights=weights/2
    nz,nx=np.subtract(chi.shape,1)
    tnodes=np.unique(np.r_[.10,np.arange(nz+1)/nz,.90])
    tnodes=tnodes[(tnodes>=.10)&(tnodes<=.90)]
    xnodes=np.arange(1,nx+1)/nx
    steps=np.array([-.0002,-.0001,0.,.0001,.0002])
    energies=np.zeros(len(steps))
    deriv=0.
    smallest=np.inf
    for t0,t1 in zip(tnodes,tnodes[1:]):
        for tn, wt in zip(nodes,weights):
            t=t0+(t1-t0)*tn
            kz=min(int(t*nz),nz-1)
            zeta=t*nz-kz
            row=chi[kz]*(1-zeta)+chi[kz+1]*zeta
            lp=legval(2*t-1,co)
            lp_t=2*legval(2*t-1,legder(co))
            # Keep the derivative quadrature inside individual order elements.
            cuts=-np.log(xnodes)/np.exp(lp)
            bounds=np.r_[.06,np.sort(cuts[(cuts>.06)&(cuts<.44)]),.44]
            da=np.diff(bounds)
            a=(bounds[:-1,None]+da[:,None]*nodes).ravel()
            weight=(da[:,None]*weights).ravel()*wt*(t1-t0)*4*np.pi*(shape.neck+a)**2
            r,rs,rz,hoop,jac,_=graph_values(a,t,shape,lp,lp_t)
            b=shape.profile(t)[0]
            B,Ba=bump(a,.06,.44)
            C,Ct=bump(t,.10,.90)
            h,ha,ht=B*C,Ba*C,B*Ct
            def interpolate(y):
                i=np.minimum((y*nx).astype(int),nx-1)
                u=y*nx-i
                return (1-u)*row[i]+u*row[i+1],nx*(row[i+1]-row[i])
            occupation,dxchi=interpolate(r/b)
            trace=rs**2+rz**2+shape.length**2+2*hoop**2
            work=stable_energy(trace,jac)
            drs=rs*h+r*ha
            drz=rz*h-r*ht
            dlogj=drs/rs+2*h
            dw=rs*drs+rz*drz+2*hoop**2*h+(-1+.5*jac*(jac-1))*dlogj
            deriv+=weight@(occupation*dw+dxchi*(r/b)*h*work)
            for k,e in enumerate(steps):
                scale=np.exp(e*h)
                new_rs,new_rz=scale*(rs+e*r*ha),scale*(rz-e*r*ht)
                new_h=scale*hoop
                new_j=-new_rs*shape.length*new_h**2
                smallest=min(smallest,float(new_j.min()))
                if new_j.min()<=0:
                    raise ArithmeticError('Virtual displacement changes orientation')
                c,_=interpolate(scale*r/b)
                energies[k]+=weight@(c*stable_energy(new_rs**2+new_rz**2+shape.length**2+2*new_h**2,new_j))
    fd=[(energies[3]-energies[1])/.0002,(energies[4]-energies[0])/.0004]
    return dict(degree=degree,analytic_slope=deriv,signed_force_secants=fd,
                changes_relative_to_zero=(energies-energies[2]).tolist(),log_displacement_steps=steps.tolist(),
                minimum_perturbed_J=smallest)


def main(field):
    parent=field.parent
    limits=json.loads((parent/'limit-validation.json').read_bytes())
    guard=json.loads((parent/'outcome.json').read_bytes())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError('Unverified radius-4 source')
    with np.load(field,allow_pickle=False) as f:
        if f['radius']!=4:
            raise ValueError('Expected radius four')
        shape=Meridian(4.,f['shape_params'],float(f.get('shape_extra',0.)))
        for d in (3,5):
            report=evaluate(shape,f['order'],f['registration'],d)
            print(json.dumps(dict(kind='material_variation',file_sha256=hashlib.sha256(field.read_bytes()).hexdigest(),
                                  entrance=4, depth=shape.length,exit=shape.exit,**report),allow_nan=False),flush=True)
            if abs(report['signed_force_secants'][0]-report['analytic_slope'])>.002*(1+abs(report['analytic_slope'])):
                raise ArithmeticError('Signed and analytic work disagree')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('field',type=Path)
    main(p.parse_args().field)
