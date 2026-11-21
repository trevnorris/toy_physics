"""Sideways GNLS momentum inside the neck and in bands next to its prescribed wall."""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_bend_band_momentum import ENTRANCE, MIDDLE
from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain


def measure(cells):
    point,_,_,flow,arch = settings()
    params = arch["parameters"]
    radius = arch["phase_A"]["q_star"][0]
    half = float(point["L"])/2
    domain = ReservoirDomain(radius,half,radius,cells,2,4,3,params)
    # The established 36-radius real domain and outgoing low-tone layer.
    def extend(x):
        return x+30*radius*(np.maximum(x/radius-2,0)/4)**3
    def stretch(x):
        return x+1j*np.maximum(x-9*radius,0)**3/(27*radius)**2
    domain.re = extend(domain.re)
    domain.ze = half+extend(domain.ze-half)
    metric = domain.metric(domain.re,domain.ze)
    domain.volume,domain.conductance,domain.inlet,domain.outlet,_,domain.operator = metric
    domain.set_current(float(flow["reference_inward_current"])/params["m"])
    vol,_,_,_,ang,lap = domain.metric(stretch(domain.re),half+stretch(domain.ze-half))
    p = domain.psi
    n = abs(p)**2
    weight = vol/3
    hc = 1.25*domain.eos_k
    h = lap/3+sparse.diags(domain.kinetic*ang+hc*(n**4-domain.bath_density**4)*weight)
    g = 8*hc*n**3*weight
    xy = g*p.real*p.imag
    pos = (domain.z_index>=domain.nfeed)&(domain.z_index<domain.neck_end)
    def neck_faces(v,dipole):
        vcell = np.empty((cells,cells),complex)
        vcell[domain.r_index[pos],domain.z_index[pos]-domain.nfeed] = v[pos]
        faces = np.empty((cells+1,cells),complex)
        faces[1:-1] = (vcell[:-1]+vcell[1:])/2
        faces[0] = 0 if dipole else (9*vcell[0]-vcell[1])/8
        faces[-1] = (3*vcell[-1]-vcell[-2])/2 if dipole else (9*vcell[-1]-vcell[-2])/8
        return faces
    base = neck_faces(p,False)
    dp = np.diff(base,axis=0)/domain.dr
    rface = domain.re[:cells+1,None]
    depth = (np.arange(cells)+.5)/cells
    for omega in (.01,.02):
        symp = 1j*omega*domain.hbar*weight
        mat = sparse.bmat([[h+sparse.diags(g*p.real**2),sparse.diags(xy-symp)],
                           [sparse.diags(xy+symp),h+sparse.diags(g*p.imag**2)]],format="csc")
        source = domain.spatial+1j*omega*domain.temporal
        response = spsolve(mat,source/2)
        power = omega/2*np.vdot(source,response).imag
        flux = domain.cuts(omega,response,round(1.25*cells),round(.75*cells))
        assert min(flux)>0 and abs(sum(flux)-power)<2e-7*power
        u,v = np.split(response,2)
        # Both real quadratures are phasors. Extract only the physical
        # coefficient in phase with bdot for exp(-i omega t).
        dip = neck_faces(-(u.imag+1j*v.imag)/omega,True)
        da = np.diff(dip,axis=0)/domain.dr
        # GNLS momentum density hbar Im(psi* grad psi). Angular projection
        # supplies 1/3 of the radial and 2/3 of the angular part.
        direct = np.zeros(dp.shape)
        volume_term = np.zeros_like(direct)
        for t in (.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)):
            r = rface[:-1]+t*domain.dr
            pval = base[:-1]+t*domain.dr*dp
            aval = dip[:-1]+t*domain.dr*da
            direct += domain.dr/2*(r**2*(aval.conj()*dp+pval.conj()*da)
                                    +2*r*pval.conj()*aval).imag
            volume_term += domain.dr*r**2*(aval.conj()*dp).imag
        boundary = np.diff(rface**2*(base.conj()*dip).imag,axis=0)
        assert np.allclose(direct,boundary+volume_term,atol=1e-13,rtol=2e-11)
        momentum = 4*np.pi*domain.hbar*domain.dw/3*direct
        for name,(lo,hi) in (("entrance",ENTRANCE),("middle",MIDDLE)):
            zz = (depth>=lo)&(depth<hi)
            for denom in (1,8,16,32):
                count = cells//denom
                assert count*denom==cells
                val = momentum[-count:,zz].sum()
                print(f"carrier cells={cells} omega={omega:g} band={name} shell_over_R={1/denom:g} "
                      f"Px_over_bdot={val:.12g}",flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,choices=(32,64),required=True)
    measure(parser.parse_args().cells)
