"""Separate physical order-profile distortion from translation in the local loaded patch calculation."""
import argparse
import hashlib
import json

import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import null_space

from mouth_breathing_mobile_patch import WIDTH, evolution, reservoir, stationary_front
from mouth_poly_shallow_release import ROOT, emit


def measurement(j,omega,length,tol):
    f=stationary_front(j,length)
    g0=f(0.)[-1].item()
    domains=(-1,1)
    outgoing=[null_space(reservoir(omega,j,s)[0].conj().T).conj().T for s in domains]
    def equations(z, y):
        return np.vstack([s*np.einsum('nij,jn->in',evolution(s*z, f, omega,j),v)
                          for s,v in zip(domains,(y[:6],y[6:]))])
    def boundary(a,b):
        return np.r_[a[:5]-a[6:11],a[4]+g0,outgoing[0]@b[:6],outgoing[1]@b[6:]]
    zz=np.unique(np.r_[np.linspace(0,24*WIDTH,151),np.linspace(24*WIDTH,length,401)])
    y=np.zeros((12,len(zz)),dtype=complex)
    # The collocation equations are complex-linear. Supplying the exact Jacobians
    # avoids ill-scaled differencing across the short, order-bearing layer.
    def jac(z,y):
        d=np.zeros((12,12,len(z)),complex)
        for i,s in enumerate(domains):
            d[6*i:6*i+6,6*i:6*i+6]=s*evolution(s*z,f,omega,j).transpose(1,2,0)
        return d
    def bcjac(a,b):
        eye=np.eye(12)
        return (np.array([boundary(a+v,b)-boundary(a,b) for v in eye]).T,
                np.array([boundary(a,b+v)-boundary(a,b) for v in eye]).T)
    solution=solve_bvp(equations,boundary,zz,y,fun_jac=jac,bc_jac=bcjac,tol=tol,max_nodes=35000)
    if not solution.success:
        raise ArithmeticError(solution.message)
    z=np.unique(np.r_[np.linspace(0,24*WIDTH,4097),np.linspace(24*WIDTH,length,2049)])
    y=solution.sol(z)
    force=g0*(y[11,0]-y[5,0])
    result=[]
    for s,sl in zip(domains,(y[:6],y[6:])):
        g=f(s*z)[-1]
        co=sl[4]+g
        i=np.argmax(abs(co))
        result.append(dict(side=s,eulerian_max=abs(sl[4]).max(),moving_profile_max=abs(co[i]),
                           signed_location=s*z[i],unperturbed_order=f(s*z[i])[2].item()))
    return force,result,int(len(solution.x))


def main(radius):
    path=ROOT/f'_scratch/s11c/mouth_breathing_disk_load01/r{radius}'
    g=json.loads((path/'outcome.json').read_bytes())
    raw=(path/'stdout').read_bytes()
    if g['exitCode'] or g['childOutcome']['exitCode'] or not g['limitsVerified']:
        raise ValueError('Disk loading is not qualified')
    rows=list(map(json.loads,raw.splitlines()))
    center,=[v for v in rows if v['kind']=='local_loading' and v['offset']==0. and v['exit_step']==.006]
    ref,=[v for v in rows if v['kind']=='result' and v['step']==.006]
    omega=center['real_axis_root']
    j=center['mean_normal_current']
    checks=center['checks'][-1]
    # Total work and the weak-loss on-shell frequency encode Im(k_patch).
    cf_im=-2*checks['mean_power_per_amplitude_squared']/(omega*checks['disk_area_motion_weight'])
    depth=ref['displacement']['peak_depth_displacement']
    emit(kind='input',radius=radius,source_sha256=hashlib.sha256(raw).hexdigest(),omega=omega,current=j,
         conditional_peak_terminal_shift=depth)
    values=[]
    for L,tol in ((14.,2e-6),(21.,5e-7)):
        f,layers,nodes=measurement(j,omega,L,tol)
        cpeak=max(v['moving_profile_max'] for v in layers)
        values.append(cpeak)
        emit(kind='translated_order',length=L,tolerance=tol,nodes=nodes,force=[f.real,f.imag],
             work_relative_error=abs(f.imag/cf_im-1),sides=layers,
             eulerian_peak_at_load=max(v['eulerian_max'] for v in layers)*depth,
             peak_order_distortion_at_load=cpeak*depth)
        if abs(f.imag/cf_im-1)>.0004:
            raise ArithmeticError('Independent collocation disagrees with the dissipative traction')
    if abs(values[0]/values[1]-1)>.002:
        raise ArithmeticError('Translated-profile order peak changed under refinement')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('radius',type=int,choices=(1,2))
    main(cli.parse_args().radius)
