"""Incremental degree-four torsional wave of the stationary finite collar.

The meridional strain tensor is supplied on the current cavity. The central
pocket outside the registered material starts locally stress-free. The joined
torsional field has continuous displacement and natural traction.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder

from mouth_exchange_shear import DEGREE
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_fields import CoupledLoads
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT, ROOT, read_case
from mouth_unified_transverse import cell_difference


def core_factors(shape, coefficients, reference, t):
    """Stretch and meridional gradient for a bonded, initially unstressed core."""
    original, baseline = reference
    b, slope, _ = shape.profile(t)
    b0, slope0, _ = original.profile(t)
    p = np.exp(legval(2*t-1,coefficients))
    p0 = np.exp(legval(2*t-1,baseline))
    q = (b/b0)*np.exp((p0-p)/2)
    dlq = (shape.length*slope/b - original.length*slope0/b0
           - p*legval(2*t-1,legder(coefficients))
           + p0*legval(2*t-1,legder(baseline)))
    length_ratio = shape.length/original.length
    return q, dlq/original.length, length_ratio


def cell_tangents(wave, shape, coefficients, quadrature=3, degree=DEGREE, axial_release=0.,
                  core_reference=None):
    """Unit-reference-order K,M. Each cell later takes its own mean live order."""
    if core_reference is not None and (coefficients is None or axial_release):
        raise ValueError("The bonded-core chart uses the unshifted axial planes of both materials")
    mesh = wave.mesh
    nt = len(wave.ids)//mesh
    i,j = np.indices((mesh,nt)).reshape(2,-1)
    x0,x1 = i/mesh,(i+1)/mesh
    dz = shape.length/nt
    km,mm = np.zeros_like(wave.kc),np.zeros_like(wave.mc)
    q,qw = leggauss(quadrature)
    for zb,wb in zip((1+q)/2,qw/2):
        t=(j+zb)/nt
        width,slope,_=shape.profile(t)
        if coefficients is not None:
            p=np.exp(legval(2*t-1,coefficients))
            pt=2*legval(2*t-1,legder(coefficients))
            cut = np.exp(-.5*p)
            if axial_release:
                edge=t.copy()
                for _ in range(4):
                    edge-=(shape.length*(edge-t)+.5*axial_release*edge*(1-edge))/(
                        shape.length+.5*axial_release*(1-2*edge))
                edge_b=shape.profile(edge)[0]
                edge_p=np.exp(legval(2*edge-1,coefficients))
                cut=edge_b*np.exp(-.5*edge_p)/width
                if np.any((cut<=0)|(cut>=1)):
                    raise ValueError("Deformed collar crossed the wall or axis")
        else:
            cut=np.ones_like(width)
        if core_reference is not None:
            stretch, shear_per_radius, zstretch = core_factors(shape,coefficients,core_reference,t)
        # Split cells at the moving core-collar material interface.
        middle=np.clip(cut,x0,x1)
        for lo,hi,collar in ((x0,middle,False),(middle,x1,True)):
            if not np.any(hi>lo):
                continue
            for xa,wa in zip((1+q)/2,qw/2):
                x=lo+(hi-lo)*xa
                u=np.clip((x-x0)*mesh,0,1)
                r=(np.where(hi>lo,x,1.) if collar else x)*width
                n=np.stack(((1-u)*(1-zb),u*(1-zb),(1-u)*zb,u*zb),axis=1)
                nr=mesh*np.array([zb-1,1-zb,-zb,zb])[None,:]/width[:,None]
                nw=np.stack((u-1,-u,1-u,u),axis=1)/dz - (x*slope)[:,None]*nr
                if collar:
                    tag_t=t
                    a=np.log(width/r)/p
                    ws=0.
                    wz=-shape.length
                    p_here=p
                    pt_here=pt
                    bt_b=slope*shape.length/width
                    if axial_release:
                        # At fixed (r,w), invert the release on each cut cell.
                        tag_t=t.copy()
                        for k in range(5):
                            tag_b,tag_slope,_=shape.profile(tag_t)
                            p_here=np.exp(legval(2*tag_t-1,coefficients))
                            pt_here=2*legval(2*tag_t-1,legder(coefficients))
                            a=np.log(tag_b/r)/p_here
                            bt_b=shape.length*tag_slope/tag_b
                            err=shape.length*(tag_t-t)+axial_release*a*tag_t*(1-tag_t)
                            if np.max(abs(err))<1e-13:
                                break
                            da_dt=bt_b/p_here-a*pt_here
                            tag_t-=err/(shape.length+axial_release*(
                                da_dt*tag_t*(1-tag_t)+a*(1-2*tag_t)))
                        else:
                            raise ArithmeticError("Unresolved inverse of the released axial tag")
                        ws=axial_release*tag_t*(1-tag_t)
                        wz=-shape.length-axial_release*a*(1-2*tag_t)
                    h=r/(shape.neck+a)
                    rs=-p_here*r
                    rz=-r*(bt_b-a*p_here*pt_here)
                    jac=(rs*wz-rz*ws)*h*h
                    if np.min(jac)<=0:
                        raise ArithmeticError("Folded reference-volume map")
                    inv_j=1/jac
                    brr=rs*rs+rz*rz
                    brw=rs*ws+rz*wz
                    bww=ws*ws+wz*wz
                    h2=h*h
                elif core_reference is not None:
                    rz = r*shear_per_radius
                    inv_j = 1/(stretch**3*zstretch)
                    h2 = stretch**2
                    brr,brw,bww = h2+rz**2,rz*zstretch,zstretch**2
                else:
                    inv_j,brr,brw,bww,h2=1.,1.,0.,1.,1.
                measure=wa*wb*dz*width*(hi-lo)*r*r*inv_j
                def outer(a,b):
                    return a[:,:,None]*b[:,None,:]
                rgrad=(np.asarray(brr)[...,None]*nr + np.asarray(brw)[...,None]*nw)
                wgrad=(np.asarray(brw)[...,None]*nr + np.asarray(bww)[...,None]*nw)
                km+=(measure*r*r)[:,None,None]*(outer(nr,rgrad)+outer(nw,wgrad))
                km+=(measure*(degree*(degree+1)-2)*h2)[:,None,None]*outer(n,n)
                mm+=(measure*r*r)[:,None,None]*outer(n,n)
    return km,mm


def install(wave, shape, coefficients, *, quadrature=3, axial_release=0., core_reference=None):
    """Keep the feed and its exact outgoing boundary; replace the cavity pencil."""
    wave.kc,wave.mc=cell_tangents(wave,shape,coefficients,quadrature,axial_release=axial_release,
                                 core_reference=core_reference)
    return wave


def active_load(shape,n,mesh,order,pole,vector,coefficients,quadrature):
    loads=CoupledLoads(shape,n,2*n,mesh,PORT,(pole,vector),pole.real)
    original_k=loads.wave.kc-cell_difference(loads.wave,shape,mesh,DEGREE)
    original_m=loads.wave.mc.copy()
    identity_k,identity_m=cell_tangents(loads.wave,shape,None,3)
    error_k=np.max(abs(identity_k-original_k))/np.max(abs(original_k))
    error_m=np.max(abs(identity_m-original_m))/np.max(abs(original_m))
    if max(error_k,error_m)>1e-11:
        raise ArithmeticError("Finite-strain triangle disagrees with undeformed wave")
    ActiveWave.install(loads)
    wave=loads.wave
    mixed_k,mixed_m=cell_tangents(wave,shape,coefficients,quadrature)
    weights=wave.project @ order.T.ravel()
    # Continuation parameter interpolates pencils only to identify the old branch.
    progress=0.
    step=.08
    overlap_min=1.
    accepted=[]
    while progress<1.:
        end=min(1.,progress+step)
        wave.kc=original_k+(mixed_k-original_k)*end
        wave.mc=original_m+(mixed_m-original_m)*end
        old,target=wave.seed,wave.target
        try:
            out,_,_,frequency,overlap=wave.pole(weights)
            if overlap<.96 or abs(frequency.real/old[0].real-1)>.15:
                raise ArithmeticError("Unresolved mode correspondence")
        except (RuntimeError,ArithmeticError,ValueError) as exc:
            wave.seed,wave.target=old,target
            step/=2
            if step<.0001:
                raise ArithmeticError(f"Unresolved spectral branch at {progress}: {exc}") from exc
            continue
        accepted.append([end,frequency.real,overlap])
        progress=end
        step=min(.16,step*1.3)
        overlap_min=min(overlap_min,overlap)
    return loads, dict(undeformed_matrix_errors=[error_k,error_m],continuation=accepted,
                        minimum_overlap=overlap_min)


def validate(radius,n,quad):
    prep,_,_,_,_,_,orig=read_case(radius,n)
    mode="bounded" if radius==1 else "relaxed"
    folder=ROOT/f"_scratch/s11c/shallow_puncture_fields01/{mode}_r{radius}_m{n}"
    src=folder/"center.npz"
    guard=json.loads((folder/"outcome.json").read_text())
    if guard["exitCode"] or not guard["limitsVerified"] or hashlib.sha256(src.read_bytes()).hexdigest() not in (folder/"stdout").read_text():
        raise ValueError("Missing verified change-of-order seed")
    with np.load(src) as data:
        if orig!=str(data["original_snapshot_sha256"].item()): raise ValueError("Different order history")
        chi=data["order"]; x=data["coordinates"]; vector=data["vector"]; old=complex(data["pole"].item())
    shape=Meridian(radius,x)
    graph=ROOT/f"_scratch/s11c/finite_bore_graph01/fields/r{radius}_unit_g2_ordered.npz"
    archive=json.loads((ROOT/"research/pde_ledger_v3/_measurements/mouth_finite_bore_support.json").read_text())
    if hashlib.sha256(graph.read_bytes()).hexdigest()!=archive["central_fields"][graph.name]:
        raise ValueError("Invalid input collar")
    with np.load(graph) as data:
        if np.max(abs(data["coordinates"]-x))>1e-12: raise ValueError("Shape mismatch")
        c=data["coeffs"]
    loads,controls=active_load(shape,n,prep["field_mesh"],chi,old,vector,c,quad)
    wave=loads.wave
    freq,grad,metrics=wave.response(chi)
    print(json.dumps(dict(kind="strained",radius=radius,mesh=n,quadrature=quad,
                          former_frequency=old.real,readout=metrics,**controls),allow_nan=False),flush=True)
    base=wave.seed
    path=Path(f"_scratch/s11c/poly_wave01/fields/r{radius}n{n}q{quad}.npz")
    path.parent.mkdir(exist_ok=True,parents=True)
    np.savez_compressed(path,pole=base[0],vector=base[1],coordinates=x,order=chi,coeffs=c)
    # Interior order perturbation leaves zero and saturated sites unchanged.
    it,jt=np.indices(chi.shape)
    window=np.sin(np.pi*it/(2*n))**2*(1-(jt/n)**2)**2
    direction=window*chi*(1-chi)
    weights=wave.project@chi.T.ravel()
    for axis in ("order","radial","axial","spacing"):
        for h in (.001,.0005):
            poles=[]; systems=[]
            for s in (-1,1):
                trial_x=x.copy(); trial_c=c.copy(); trial_chi=chi
                if axis=="radial": trial_x[2]+=s*h/radius
                if axis=="axial": trial_x[0]+=np.log1p(s*h/shape.length)
                if axis=="order": trial_chi=chi+s*h*direction
                if axis=="spacing": trial_c[0]+=s*h
                trial=Meridian(radius,trial_x)
                obj=CoupledLoads(trial,n,2*n,prep["field_mesh"],PORT,base,base[0].real)
                install(obj.wave,trial,trial_c,quadrature=quad)
                ActiveWave.install(obj)
                v=obj.wave
                w=v.project@trial_chi.T.ravel()
                out,_,mass,omega,ov=v.pole(w)
                if ov<.96: raise ArithmeticError("Signed mode switch")
                poles.append(omega.real)
                stiff=v.cut_k+v.matrix(w[:,None,None]*v.kc)
                systems.append((stiff,mass))
            if axis=="order":
                predicted=float(np.sum(grad*direction))
            else:
                y=base[1]; om=base[0]
                edge=y[wave.port.ids]
                db=(wave.port.block(om+2.e-6)-wave.port.block(om-2.e-6))/4.e-6
                m=wave.cut_m+wave.matrix(weights[:,None,None]*wave.mc)
                den=2*om*(y@(m@y))+edge@(db@edge)
                ds=systems[1][0]-systems[0][0]-om**2*(systems[1][1]-systems[0][1])
                predicted=float(((y@(ds@y))/(2*h*den)).real)
            measured=(poles[1]-poles[0])/(2*h)
            print(json.dumps(dict(kind="slope",radius=radius,mesh=n,quadrature=quad,axis=axis,
                                  step=h,predicted=predicted,measured=measured)),flush=True)
            if abs(predicted-measured)>3e-4*max(1,abs(measured)):
                raise ArithmeticError("Signed pole work disagrees with reciprocal matrix variation")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--radius",type=int,required=True,choices=(1,2))
    ap.add_argument("--mesh",type=int,required=True,choices=(48,72))
    ap.add_argument("--quadrature",type=int,default=3,choices=(3,5))
    args=ap.parse_args()
    validate(args.radius,args.mesh,args.quadrature)
