"""Continue the strained mode to the stabilized static collar and release an internal label."""
import argparse
import hashlib
import json

import numpy as np

from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_fields import CoupledLoads
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT, ROOT, read_case
from mouth_strained_wave import install


def released_work(radius,n,quad,adopt):
    prep,*_=read_case(radius,n)
    old_dir=ROOT/f"_scratch/s11c/poly_wave01/r{radius}n{n}q{quad}"
    refit=ROOT/f"_scratch/s11c/polycrystal_material01/r{radius}g2"
    cases=[json.loads((p/"outcome.json").read_text()) for p in (old_dir,refit)]
    if any(r["exitCode"] or not r["limitsVerified"] for r in cases):
        raise ValueError("Both starting states must be validated")
    source=ROOT/f"_scratch/s11c/poly_wave01/fields/r{radius}n{n}q{quad}.npz"
    with np.load(source) as data:
        chi=data["order"]; x=data["coordinates"]; c0=data["coeffs"]
        pole=complex(data["pole"].item());v=data["vector"]
    file=refit/"coeffs.npz"
    sha=hashlib.sha256(file.read_bytes()).hexdigest()
    if sha not in (refit/"stdout").read_text():
        raise ValueError("Missing verified neutral-law material reference")
    with np.load(file) as data:
        c=data["ordered"]; ref_x=data["coordinates"]
    if np.max(abs(ref_x-x))>1e-12: raise ValueError("Static and wave shapes differ")
    shape=Meridian(radius,x)
    obj=CoupledLoads(shape,n,2*n,prep["field_mesh"],PORT,(pole,v),pole.real)
    ActiveWave.install(obj)
    wave=obj.wave
    weights=wave.project@chi.T.ravel()
    target=ROOT/f"_scratch/s11c/poly_wave01/fields/selected_r{radius}n{n}q{quad}.npz"
    if adopt:
        folder=ROOT/f"_scratch/s11c/poly_wave01/selected_r{radius}n{n}q{quad}"
        meta=json.loads((folder/"outcome.json").read_text())
        if (not meta["limitsVerified"] or meta["exitCode"]!=1
                or "operands could not be broadcast" not in (folder/"stderr").read_text()):
            raise ValueError("Wrong guarded failure to resume")
        row,=map(json.loads,(folder/"stdout").read_text().splitlines())
        if row["sha256"]!=hashlib.sha256(target.read_bytes()).hexdigest() or row["graph_sha256"]!=sha:
            raise ValueError("Stored center differs from successful stage")
        with np.load(target) as data:
            if not (np.array_equal(chi,data["order"]) and np.array_equal(c,data["coeffs"])
                    and np.array_equal(x,data["coordinates"])):
                raise ValueError("Mismatched saved material and pole")
            base=complex(data["pole"].item()),data["vector"]
        install(wave,shape,c,quadrature=quad)
        k=wave.cut_k+wave.matrix(weights[:,None,None]*wave.kc)
        m=wave.cut_m+wave.matrix(weights[:,None,None]*wave.mc)
        om,vec=base
        residual=(k-om**2*m)@vec
        residual[wave.port.ids]-=wave.port.block(om)@vec[wave.port.ids]
        error=np.linalg.norm(residual)/np.linalg.norm(k@vec)
        if error>1e-8 or abs(np.vdot(vec,m@vec).real-1)>1e-8:
            raise ArithmeticError("Archived vector fails updated pencil")
        print(json.dumps(dict(**row,recovered_relative_residual=error)),flush=True)
    else:
        for t in np.linspace(0.,1.,8)[1:]:
            install(wave,shape,(1-t)*c0+t*c,quadrature=quad)
            out,_,_,om,ov=wave.pole(weights)
            if ov<.95:
                raise ArithmeticError("Lost pole in mechanics update")
        freq,_,metrics=wave.response(chi)
        base=wave.seed
        np.savez_compressed(target,pole=base[0],vector=base[1],coordinates=x,order=chi,coeffs=c)
        print(json.dumps(dict(kind="selected",radius=radius,n=n,quad=quad,
                              source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                              graph_sha256=sha,sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                              readout=metrics)),flush=True)
    om,vec=base
    boundary=vec[wave.port.ids]
    db=(wave.port.block(om+2.e-6)-wave.port.block(om-2.e-6))/4.e-6
    mass=wave.cut_m+wave.matrix(weights[:,None,None]*wave.mc)
    den=2*om*(vec@(mass@vec))+boundary@(db@boundary)
    for h in (.002,.001):
        read=[]; pencils=[]
        for sign in (-1,1):
            install(wave,shape,c,quadrature=quad,axial_release=sign*h)
            wave.seed=base;wave.target=om.real
            out,_,mass,fp,ov=wave.pole(weights)
            if ov<.95: raise ValueError("Released pole mismatch")
            read.append(fp.real)
            k=wave.cut_k+wave.matrix(weights[:,None,None]*wave.kc)
            pencils.append(k-om**2*mass)
        raw=(read[1]-read[0])/(2*h)
        predicted=((vec@((pencils[1]-pencils[0])@vec))/(2*h*den)).real
        print(json.dumps(dict(kind="axial_label",step=h,measured=raw,reciprocal=predicted)),flush=True)
        if abs(raw-predicted)>1e-5*max(1.,abs(raw)):
            raise ArithmeticError("Wave load on the released axial displacement mismatched")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--radius",type=int,choices=(1,2),required=True)
    p.add_argument("--mesh",type=int,choices=(48,72),required=True)
    p.add_argument("--quadrature",type=int,choices=(3,5),default=3)
    p.add_argument("--adopt",action="store_true",help="Validate a saved refitted pole before signed work")
    a=p.parse_args()
    released_work(a.radius,a.mesh,a.quadrature,a.adopt)
