"""Frictionless curved-guide benchmark, fixed-angle control and independently resolved carrier flux."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import prepare, trace as field_trace
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, emit


def guide(tr, cuts, straight):
    # x=C e + R n. To first order n=n0+q b P_n e. For no angular
    # traction d(R² db/dt)/dt=−R d²H/dt². K=R² db/dt+R dH/dt removes
    # derivatives of the archived speed even across finite-element interfaces.
    state=np.zeros(4)  # K, b, integrated H' x-speed, integrated H'^2 u
    edges=np.unique(np.r_[cuts,tr.edges[(tr.edges>cuts[0])&(tr.edges<cuts[-1])]])
    for a,end in zip(edges[:-1],edges[1:],strict=True):
        def values(w,state):
            w=np.array([np.clip(w,np.nextafter(a,end),np.nextafter(end,a))])
            radius,dr=tr.geometry.radius(w)
            if straight: radius,dr=tr.geometry.radius(cuts[:1])[0],np.zeros(1)
            h,hp=taper(w,cuts[0],cuts[2])
            R,Rp,H,Hw,u=(float(x[0]) for x in (radius,dr,1-h,-hp,tr.speed(w)))
            K,b=state[:2]
            vx=Hw*u/3+2*(Rp*u*b+K/R)/3
            return R,Rp,H,Hw,u,vx
        def ode(w,s):
            R,Rp,H,Hw,u,vx=values(w,s)
            return [Rp*Hw*u, s[0]/(u*R**2)-Hw/R,Hw*vx,Hw*Hw*u]
        sol=solve_ivp(ode,(a,end),state,method='DOP853',rtol=2e-11,atol=1e-13)
        if not sol.success: raise ArithmeticError(sol.message)
        state=sol.y[:,-1]
        if end in cuts:
            R,Rp,H,Hw,u,vx=values(end,state)
            conjugate=H*vx-state[2]
            if straight:
                np.testing.assert_allclose([vx,conjugate],[Hw*u/3,(H*Hw*u-state[3])/3],atol=3e-11)
            emit(test='hard_guide',straight=straight,depth=end,net_holder_slope=vx,
                 shape_holder_slope=conjugate,angular_dipole_slope=2*state[1],
                 fixed_angle_net_holder_slope=Hw*u,fixed_angle_shape_holder_slope=H*Hw*u-state[3])


def carrier_flux(dom,field,w,thickness):
    radius=dom.radius(w)[0]
    v_limit=min(float(thickness/radius),1.)
    radial=dom.blocks[2][1]
    cuts=np.r_[0.,radial[(radial>0)&(radial<v_limit)],v_limit]
    nodes,weights=leggauss(8)
    half=np.diff(cuts)[:,None]/2
    v=(cuts[:-1,None]+half*(1+nodes)).ravel()
    r,_,_,rv,_,_,f,_,fw=field_trace(dom,2,np.full_like(v,w),v,field[:,None])
    density=(f[:,0].conj()*fw[:,0]).imag
    return (half*weights).ravel()@(4*np.pi*r*r*abs(rv)*density),density.min()


def main(mesh):
    tr,rate=cubic_trace(mesh)
    lo,hi=1.175,tr.crossing(.5)
    cuts=np.array([lo,(lo+hi)/2,hi,tr.crossing(.25),tr.crossing(.1),tr.crossing(.01)])
    for straight in (False,True): guide(tr,cuts,straight)
    dom,field=prepare(SimpleNamespace(source=Path(f'_scratch/s11c/mouth_flared_cubic_n{mesh}')))
    # Only the inside of the guide belongs to this carrier cross-section.
    # The bonded-site calculation itself was launched at the potential minimum.
    half_well=PITCH*np.log1p(1/np.sqrt(2))/BETWEEN[1]
    for w in cuts:
        radius=dom.radius(w)[0]
        total,jmin=carrier_flux(dom,field,w,radius)
        check=dom.section(field,w)['integrated_jw']
        np.testing.assert_allclose(total,check,atol=1e-13)
        if jmin<=0: raise ArithmeticError('Axial flow reversed on a sampled section')
        for thickness in [half_well,PITCH,2*PITCH]:
            flux,_=carrier_flux(dom,field,w,thickness)
            emit(test='carrier_layer',depth=w,thickness=thickness,
                 first_radial_cell_width=radius*np.diff(dom.blocks[2][1])[0],
                 section_current=total,inlet_intake=rate,layer_current=flux,section_fraction=flux/total)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',required=True,type=int)
    main(parser.parse_args().mesh)
