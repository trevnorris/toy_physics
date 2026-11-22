"""Verify that a formed central pocket follows the inner collar material surface."""
import argparse
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import legval

from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_fields import CoupledLoads
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT, ROOT, read_case
from mouth_strained_wave import core_factors, install


def cut(shape,c,t):
    return shape.profile(t)[0]*np.exp(-.5*np.exp(legval(2*t-1,c)))


def check_reference(original, old, current, coeff, t):
    stretch, shear, axial = core_factors(current,coeff,(original,old),t)
    initial = cut(original,old,t)
    final = cut(current,coeff,t)
    gap = np.max(abs(stretch*initial-final))
    d=1e-5
    def current_radius(t_shift):
        return .7*initial*cut(current,coeff,t_shift)/cut(original,old,t_shift)
    finite=(current_radius(t+d)-current_radius(t-d))/(2*d*original.length)
    shear_error = np.max(abs(finite-.7*final*shear))
    if gap>1e-12 or shear_error>1e-7 or np.min(stretch**3*axial)<=0:
        raise ArithmeticError("Central material and collar kinematics differ")
    return float(gap),float(shear_error)


def run(radius,n):
    prep,*_=read_case(radius,n)
    folder=ROOT/f'_scratch/s11c/poly_wave01/release_r{radius}n{n}q3'
    out=json.loads((folder/'outcome.json').read_text())
    row=json.loads((folder/'stdout').read_text().splitlines()[0])
    f=ROOT/f'_scratch/s11c/poly_wave01/fields/selected_r{radius}n{n}q3.npz'
    if out['exitCode'] or not out['limitsVerified'] or row['sha256']!=hashlib.sha256(f.read_bytes()).hexdigest():
        raise ValueError('Unverified stabilized reference')
    with np.load(f) as data:
        x,chi,c=data['coordinates'],data['order'],data['coeffs']
        base=(complex(data['pole'].item()),data['vector'])
    om,v=base
    shape0=Meridian(radius,x)
    mesh=prep['field_mesh']
    def system(shape, coeff, reference):
        op=CoupledLoads(shape,n,2*n,mesh,PORT,base,om.real)
        install(op.wave,shape,coeff,core_reference=reference)
        ActiveWave.install(op)
        wave=op.wave
        weights=wave.project@chi.T.ravel()
        K=wave.cut_k+wave.matrix(weights[:,None,None]*wave.kc)
        M=wave.cut_m+wave.matrix(weights[:,None,None]*wave.mc)
        return wave,weights,K,M
    reference=(shape0,c)
    src,_,K,M=system(shape0,c,reference)
    _,_,K0,M0=system(shape0,c,None)
    def rel(a,b):
        return np.linalg.norm((a-b).data)/np.linalg.norm(b.data)
    identity=[rel(K,K0),rel(M,M0)]
    print(json.dumps(dict(kind="bonded_center",radius=radius,mesh=n,equality=identity,frequency=om.real)),flush=True)
    if max(identity)>1e-11: raise ArithmeticError("Core state is prestrained at reference")
    ev=v[src.port.ids]
    db=(src.port.block(om+2e-6)-src.port.block(om-2e-6))/4e-6
    den=2*om*(v@(M@v))+ev@(db@ev)
    for axis in ('radial','axial','spacing'):
        for h in (.002,.001):
            poles=[];matrices=[];glued=[]
            for sign in (-1.,1.):
                xx=x.copy();co=c.copy()
                if axis=='radial': xx[2]+=sign*h/radius
                if axis=='axial': xx[0]+=np.log1p(sign*h/shape0.length)
                if axis=='spacing': co[0]+=sign*h
                shape=Meridian(radius,xx)
                check=check_reference(shape0,c,shape,co,np.linspace(.02,.98,51))
                wave,w,k,m=system(shape,co,reference)
                _,_,_,omega,overlap=wave.pole(w)
                if overlap<.98: raise ArithmeticError('Lost bonded mode on signed displacement')
                glued.append(check);poles.append(omega.real);matrices.append(k-om**2*m)
            direct=(poles[1]-poles[0])/(2*h)
            predicted=((v@((matrices[1]-matrices[0])@v))/(2*h*den)).real
            print(json.dumps(dict(kind='signed',axis=axis,radius=radius,mesh=n,step=h,
                                  direct=direct,reciprocal=predicted,cut_and_shear_errors=glued)),flush=True)
            if abs(direct-predicted)>3e-4*max(abs(direct),1.):
                raise ArithmeticError('Unresolved bonded angular work')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--mesh',type=int,choices=(48,72),required=True)
    a=p.parse_args();run(a.radius,a.mesh)
