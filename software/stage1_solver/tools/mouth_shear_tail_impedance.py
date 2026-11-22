"""Scattering on one plane shear ray through the measured conversion profile.

Join a uniform ordered-side lead at chi=.99 to a measured axial line at r=.25R.
Continue the remote weak material with the local exponential decay length.
Compare residual and vanishing transverse inertia. This projects out the rim,
transverse gradients and other bulk channels.
"""
import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

from mouth_shear_ringdown_state import domain
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit, quadratic, tensor


def profile(radius,mesh):
    dom=domain(radius,mesh)
    if radius==1.:
        folder=ROOT/f'mouth_free_order_cubic_n{mesh}_d04'
        if mesh==48:
            # The verified deep state is the residual-checked continuation.
            folder=folder.with_name(folder.name+'_outer_residual')
        path=folder/'order_seed1.55.npz'
        expected=json.loads(ARCHIVE.read_text())['diffuse_orders_sha256'][str(path)]
        key='chi'
    else:
        tag={.5:'half',2.:'double'}[radius]
        folder=ROOT/f'shear_state_{tag}_{mesh}_01'
        if json.loads((folder/'outcome.json').read_text())['exitCode'] or not json.loads((folder/'limit-validation.json').read_text())['verified']:
            raise ValueError('Order preparation failed')
        path=folder/'state.npz'
        expected=json.loads((folder/'stdout').read_text().splitlines()[-1])['sha256']
        key='order'
    if digest(path)!=expected:
        raise ValueError('Invalid order profile')
    with np.load(path) as data:
        chi=data[key]
    def sample(z):
        z=np.atleast_1d(z)
        values=np.zeros_like(z,dtype=float)
        low=z<.5
        if np.any(low):
            values[low]=dom.read(chi,np.full(sum(low),.25),z[low])[0]
        if np.any(~low):
            w=z[~low]
            p,v,dofs=dom.blocks[2]
            i=np.clip(np.searchsorted(p,w)-1,0,len(p)-2)
            vq=1-.25/dom.radius(w)[0]
            j=np.clip(np.searchsorted(v,vq)-1,0,len(v)-2)
            a=quadratic(2*(w-p[i])/np.diff(p)[i]-1)[0]
            b=quadratic(2*(vq-v[j])/np.diff(v)[j]-1)[0]
            values[~low]=np.einsum('qi,qi->q',tensor(a,b),chi[dofs[i*(len(v)-1)+j]])
        return values
    z=np.linspace(.01,2.8,801)
    y=sample(z)
    def crossing(level):
        hits=np.flatnonzero((y[:-1]>level)&(y[1:]<=level))
        if len(hits)!=1:
            raise ValueError('Expected one falling order crossing')
        i=hits[0]
        return brentq(lambda w:sample(w)[0]-level,z[i],z[i+1])
    return sample,crossing,path


def tail_derivative(omega,ell,cut,floor):
    # u=exp(y/2ell)*H_order(2ell*omega*sqrt(floor/cut)*exp(y/2ell)).
    # For zero floor select the outgoing or finite-energy continuation directly.
    if floor==0.:
        q=np.sqrt(complex(omega**2-1/(4*ell**2)))
        return 1/(2*ell)+1j*q
    x=2*ell*omega*np.sqrt(floor/cut)
    nu=mp.sqrt(1-4*(ell*omega)**2*(1-floor))
    h=mp.hankel1(nu,x)
    return complex((1+x*(mp.hankel1(nu-1,x)-mp.hankel1(nu+1,x))/(2*h))/(2*ell))


def scatter(order,step,lead,cut,floor,omega,ell):
    rho=floor+(1-floor)*order
    wave_number=omega*np.sqrt(rho/order)
    impedance=order*wave_number
    phases=wave_number*step
    # Propagate (displacement, physical shear traction) backwards from the
    # asymptotic outgoing solution without placing a reflector in the bulk.
    u=1.+0j
    traction=cut*tail_derivative(omega,ell,cut,floor)
    outgoing=(u.conjugate()*traction).imag
    for cosine,sine,z in zip(np.cos(phases)[::-1],np.sin(phases)[::-1],impedance[::-1]):
        u,traction=cosine*u-sine*traction/z,z*sine*u+cosine*traction
    z0=omega*np.sqrt(lead*(floor+(1-floor)*lead))
    incident=(u+traction/(1j*z0))/2
    reflected=(u-traction/(1j*z0))/2
    reflect=abs(reflected/incident)**2
    transmit=outgoing/(z0*abs(incident)**2)
    return reflect,transmit,abs(reflect+transmit-1)


def main(a):
    for r in (.5,1.,2.):
        sample,crossing,path=profile(r,a.mesh)
        lead=.99
        start=crossing(lead)
        for cut in (.0001,.00001):
            end=crossing(cut)
            ell=(end-crossing(10*cut))/np.log(10.)
            edges=np.linspace(start,end,a.cells+1)
            order=sample((edges[:-1]+edges[1:])/2)
            if np.min(order)<=0:
                raise ValueError('Nonpositive shear order')
            emit(test='tail',radius=r,r_over_radius=.25,source=str(path),source_hash=digest(path),
                 mesh=a.mesh,ray_cells=a.cells,start_w=start,end_w=end,cut=cut,ell=ell,
                 zero_floor_cutoff=1/(2*ell))
            for f in (1.,.0625,0.):
                for om in (8.,12.2,16.):
                    refl,tr,err=scatter(order,edges[1]-edges[0],lead,cut,f,om,ell)
                    emit(test='ray',radius=r,cut=cut,inertia_floor=f,omega=om,
                         reflected_power=refl,outgoing_power=tr,flux_residual=err)
                    if not np.isfinite(tr) or err>1e-8 or tr < -1e-12:
                        raise ValueError('Lossless two-port work balance failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--cells',type=int,required=True)
    main(p.parse_args())
