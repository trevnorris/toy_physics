"""Direct finite angular-twist energy and cut-cell kinetic norm of a strained shallow pocket."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss, legder, legval
from scipy.optimize import minimize_scalar

from mouth_exchange_shear import DEGREE
from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import DELTA, Meridian
from mouth_shallow_puncture_fields import PORT, ROOT, read_case


def finite_twist():
    # Azimuthal rearrangement phi -> phi + eps*f(r,w)*P_ell'(cos(theta))
    # is exactly area-preserving on S² at every r,w.
    y,weights=leggauss(22)
    p=Legendre.basis(DEGREE)
    sn=np.sqrt(1-y*y)
    t0=p.deriv()(y)
    dt=-sn*p.deriv(2)(y)
    angular_norm=weights@((sn*t0)**2)
    probes=[]
    for mer,h in ((np.array([[-.8,1.2],[0.,-1.]]),.4),
                  (np.array([[-2.1,-.9],[.13,-.7]]),1.75)):
        r,u,fr,fw=.72,.85,.34,-.57
        F=np.zeros((4,4)); F[:2,:2]=mer; F[2,2]=F[3,3]=h
        B=mer@mer.T; J=np.linalg.det(F)
        def averaged_energy(eps):
            A=np.broadcast_to(np.eye(4),(len(y),4,4)).copy()
            A[:,3,0]+=eps*sn*t0*r*fr
            A[:,3,1]+=eps*sn*t0*r*fw
            A[:,3,2]+=eps*sn*dt*u
            deformed=A@F
            det=np.linalg.det(deformed)
            if np.max(abs(det/J-1))>1e-12:
                raise ArithmeticError("Torsional displacement changed volume")
            return weights@((np.sum(deformed**2,axis=(1,2))-4)/2 -
                            np.log(det)+(det-1)**2/4)
        step=.007
        measured=(averaged_energy(step)+averaged_energy(-step)-2*averaged_energy(0))/(step**2*angular_norm*J)
        grad=np.array([fr,fw])
        predicted=(r*r*(grad@B@grad)+(DEGREE*(DEGREE+1)-2)*h*h*u*u)/J
        probes.append(dict(J=J,direct=measured,predicted=predicted))
        if abs(measured/predicted-1)>2e-10:
            raise ArithmeticError("Finite twist fails large-strain angular Hessian")
    return probes


def mass_interval(op,shape,mesh,coeff,chi,v,zlo,zhi,ng):
    cell=op.project@chi.T.ravel()
    nt=len(cell)//mesh
    weights=cell.reshape(mesh,nt)
    nodes=np.r_[v,0j][op.ids].reshape(mesh,nt,4)
    x0=np.arange(mesh)[:,None]/mesh
    x1=x0+1/mesh
    z0=np.arange(nt)[None,:]*shape.length/nt
    zlow=np.maximum(zlo,z0)
    dz=np.maximum(0.,np.minimum(zhi,z0+shape.length/nt)-zlow)
    total=0.
    q,qw=leggauss(ng)
    for zb,wb in zip((q+1)/2,qw/2):
        z=zlow+dz*zb
        zb0=(z-z0)*nt/shape.length
        t=z/shape.length
        width=shape.radius(z)[0]
        p=np.exp(legval(2*t-1,coeff))
        mid=np.clip(np.exp(-p/2),x0,x1)
        for left,right,is_collar in ((x0,mid,False),(mid,x1,True)):
            for xa,wa in zip((q+1)/2,qw/2):
                x=left+(right-left)*xa
                a=(x-x0)*mesh
                r=x*width
                value=((1-a)*(1-zb0)*nodes[:,:,0]+a*(1-zb0)*nodes[:,:,1]
                       +(1-a)*zb0*nodes[:,:,2]+a*zb0*nodes[:,:,3])
                if is_collar:
                    # At degenerate quadrature regions x can coincide with zero.
                    s=shape.neck-np.log(np.maximum(x,1.e-60))/p
                    invj=s*s/(shape.length*p*np.maximum(r,1.e-60)**3)
                else:
                    invj=1.
                measure=wa*wb*(right-left)*width*dz*r**4*invj
                total+=np.sum(measure*weights*abs(value)**2)
    return total


def run(root,radius,n,path,check):
    # Verify the pole was emitted by the guarded calculation in this preparation.
    meta=json.loads((root/"outcome.json").read_text())
    if not meta["limitsVerified"] or meta["exitCode"] or meta["childOutcome"]["exitCode"]:
        raise ValueError("Missing validated pole")
    refrows=[json.loads(line) for line in (root/"stdout").read_text().splitlines()]
    prep,*_=read_case(radius,n)
    with np.load(path) as data:
        shape=Meridian(radius,data["coordinates"]);v=data["vector"]
        chi=data["order"];c=data["coeffs"];pole=complex(data["pole"].item())
    if check=="initial":
        row,=[r for r in refrows if r["kind"]=="strained"]
        reference=row["readout"]
        checksum=hashlib.sha256(path.read_bytes()).hexdigest()
    else:
        row,=[r for r in refrows if r["kind"]=="selected"]
        checksum=hashlib.sha256(path.read_bytes()).hexdigest()
        if row["sha256"]!=checksum: raise ValueError("Refitted collar checkpoint mismatch")
        reference=row["readout"]
    if abs(pole.real-reference["wave_frequency"])>1e-9:
        raise ValueError("Pole mismatch")
    m=prep["field_mesh"]
    op=OccupiedWave(shape,n,2*n,m,PORT)
    waist=minimize_scalar(lambda z:shape.radius(z)[0],bounds=(0,shape.length),method="bounded").x
    for ng in (3,5):
        edges=(0.,waist-2*DELTA,waist+2*DELTA,shape.length)
        areas=np.array([mass_interval(op,shape,m,c,chi,v,a,b,ng) for a,b in zip(edges[:-1],edges[1:])])
        total=mass_interval(op,shape,m,c,chi,v,0.,shape.length,ng)
        print(json.dumps(dict(radius=radius,mesh=n,quad=ng,source_sha256=checksum,
                              cavity=total,reference_cavity=reference["wave_pocket_fraction"],
                              fractions=(areas/total).tolist(),sum_error=sum(areas)/total-1,
                              waist=waist)),flush=True)
        if abs(total-reference["wave_pocket_fraction"])>2e-4 or abs(sum(areas)/total-1)>1e-4:
            raise ArithmeticError("Cavity kinetic partition or normalization mismatch")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--radius",type=int,required=True,choices=(1,2))
    ap.add_argument("--mesh",type=int,required=True,choices=(48,72))
    ap.add_argument("--check",choices=("initial","selected"),default="initial")
    args=ap.parse_args()
    print(json.dumps(dict(finite_angular_twist=finite_twist())),flush=True)
    run(args.root,args.radius,args.mesh,args.source,args.check)
