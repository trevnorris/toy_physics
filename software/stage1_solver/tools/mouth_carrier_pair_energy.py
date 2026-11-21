"""Weak-phase, separated round mouths with shared feed and receiving reservoirs.

Diagonal conductance and added mass come from independent resolved throats.
Use the Neumann Green function of the laterally infinite 3D slab of thickness
t and that of the half-space in four spatial dimensions for mutual loading.
Only the monopole throughput varies; radial emission profiles are trial data.
"""
from pathlib import Path
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from mouth_flared_receiver import emit

BASE=Path('_scratch/s11c')
REFERENCE={.5:'h',1.:'one',2.:'two'}


def aperture(radius,count,profile):
    nodes,weights=leggauss(count)
    x=(nodes+1)/2
    law=np.ones(count) if profile=='uniform' else (1-x*x)**2
    weights=weights*x*x*law
    return radius*x,weights/weights.sum()


def offsets(radii,count,profile):
    a,wa=aperture(radii[0],count,profile)
    b,wb=aperture(radii[1],count,profile)
    mu,wm=leggauss(count)
    q=np.sqrt(a[:,None,None]**2+b[None,:,None]**2-2*a[:,None,None]*b[None,:,None]*mu)
    weight=wa[:,None,None]*wb[None,:,None]*wm/2
    return q,weight


def kernels(d,q,kind,thickness=1.):
    # The relative offset vector is isotropic. Average over its angle to
    # the line of centers by integrating s G(s) between d−q and d+q.
    lo,hi=d-q,d+q
    if lo.min()<=0:
        raise ValueError('Overlapping projected apertures')
    if kind=='bulk':
        integral=np.log1p(2*q/lo)/(2*np.pi**2)
        edges=(1/hi-1/lo)/(2*np.pi**2)
    else:
        k=np.pi/thickness
        integral=(2*q+2/k*(np.log1p(-np.exp(-k*hi))-np.log1p(-np.exp(-k*lo))))/(4*np.pi*thickness)
        edges=(1/np.tanh(k*hi/2)-1/np.tanh(k*lo/2))/(4*np.pi*thickness)
    mean=integral/(2*d*q)
    return mean, mean/d-edges/(2*d*q)


def pair(d,calibration,offsets_by_side,bias):
    mutual=[]
    slope=[]
    for kind,(q,weight) in zip(('feed','bulk'),offsets_by_side):
        g,s=kernels(d,q,kind)
        mutual.append(np.sum(weight*g))
        slope.append(np.sum(weight*s))
    matrix=np.diag([1/x['conductance_energy'] for x in calibration])
    matrix[0,1]=matrix[1,0]=sum(mutual)
    if np.linalg.eigvalsh(matrix)[0]<=0:
        raise ArithmeticError('Two-port energy ceased to be positive')
    rates=np.linalg.solve(matrix,bias*np.ones(2))
    energy=.5*bias*sum(rates)
    forces=np.prod(rates)*np.array(slope)
    return dict(energy=energy,rates=rates.tolist(),force_feed=forces[0],
                force_bulk=forces[1],force_toward=forces.sum(),mutual=mutual)


def read(radius):
    path=BASE/f'carrier_weak_wide_{REFERENCE[radius]}_01'
    if json.loads((path/'outcome.json').read_text())['exitCode'] or not json.loads((path/'limit-validation.json').read_text())['verified']:
        raise ValueError('Unverified isolated aperture')
    return json.loads((path/'stdout').read_text())


def verify_green():
    # Independently sum reflected 4D point sources separated along w by 2t.
    # Bound the positive omitted image tail by the decreasing integral.
    t=1.
    modes=np.arange(1,20001)
    for s in (2.5,8.,17.):
        partial=(1/s**2+2*np.sum(1/(s*s+4*t*t*modes**2)))/(2*np.pi**2)
        lower=partial+quad(lambda n:1/(s*s+4*t*t*n*n),modes[-1]+1,np.inf,epsabs=1e-14)[0]/np.pi**2
        upper=partial+quad(lambda n:1/(s*s+4*t*t*n*n),modes[-1],np.inf,epsabs=1e-14)[0]/np.pi**2
        exact=1/np.tanh(np.pi*s/(2*t))/(4*np.pi*t*s)
        emit(test='slab_images',r=s,lower=lower,formula=exact,upper=upper)
        if not lower<=exact<=upper:
            raise ArithmeticError('Slab image sum and kernel disagree')


if __name__=='__main__':
    verify_green()
    for radii in ((1.,1.),(.5,2.)):
        calibration=[read(r) for r in radii]
        for profile in ('uniform','center_weighted'):
            for count in (16,28):
                offsets_by_side=[offsets(np.array(radii)*k,count,profile) for k in (1.,1.5)]
                for d in (8.,16.,32.,64.):
                    bias=.05
                    result=pair(d,calibration,offsets_by_side,bias)
                    step=d*.0003
                    diff=(pair(d+step,calibration,offsets_by_side,bias)['energy']-
                          pair(d-step,calibration,offsets_by_side,bias)['energy'])/(2*step)
                    rates=result['rates']
                    normalizer=rates[0]*rates[1]/(4*np.pi*d*d)
                    accelerations=[result['force_toward']/r['mass_from_energy'] for r in calibration]
                    emit(test='pair',radii=radii,bias=bias,profile=profile,count=count,d=d,
                         energy_derivative=diff,normalized_feed_force=result['force_feed']/normalizer,
                         accelerations_diagonal_only=accelerations,**result)
                    np.testing.assert_allclose(diff,result['force_toward'],rtol=4e-7)
