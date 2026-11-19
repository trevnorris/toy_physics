"""Project the actual bend-driven aperture trace onto straight-neck waves.

Each end uses a uniform local reference tube at the measured mean intake.
Evanescent interference can transmit work across a cut next to the junction.
Report this contribution separately from the traveling-wave return ratio.
"""
import argparse

import numpy as np

from mouth_aperture_return import IDENTITY,J,radial_basis,rotate
from mouth_open_reservoirs import ReservoirDomain, settings


class ApertureTrace(ReservoirDomain):
    def cuts(self,omega,u,cut_r,cut_w):
        self.trace_field=u
        return super().cuts(omega,u,cut_r,cut_w)


def project(d,omega,rate,side):
    upstream=side=="brane"
    inside=d.nfeed if upstream else d.neck_end-1
    other=inside-1 if upstream else inside+1
    ids=np.nonzero((d.z_index==inside)&(d.r_index<d.mouth_cells))[0]
    ext=np.nonzero((d.z_index==other)&(d.r_index<d.mouth_cells))[0]
    u=np.array(np.split(d.trace_field,2))
    wstep=(np.diff(d.ze)[other]+d.dw)/2
    # Extrapolated ghost preserves the exact finite-volume face derivative.
    ghost=u[:,ids]+(d.dw/wstep)*(u[:,ext]-u[:,ids])
    area,modes,basis=radial_basis(len(ids),d.dr)
    dc=-rate if upstream else rate
    delta=np.arcsin(dc*d.constituent_mass*d.dw/(3*sum(area)*d.bath_density*d.hbar))
    plane_phase=np.angle(np.dot(area,d.psi[ids]))
    first=rotate(-plane_phase)@(u[:,ids]*np.sqrt(area))@basis
    second=rotate(-plane_phase-delta)@(ghost*np.sqrt(area))@basis
    hop=d.kinetic/d.dw**2
    twist=rotate(delta)
    outward=inward=near=0.
    for j,lam in enumerate(modes):
        central=(2*hop*np.cos(delta)+d.kinetic*lam)*IDENTITY+np.diag([10*d.eos_k*d.bath_density**4,0])+1j*d.hbar*omega*J
        T=np.block([[np.zeros((2,2)),IDENTITY],[-twist.T@twist.T,twist.T@central/hop]])
        z,states=np.linalg.eig(T)
        v=states[:2]
        amplitudes=np.linalg.solve(np.vstack([v,v*z]),np.r_[first[:,j],second[:,j]])
        scale=omega*d.kinetic/d.dw
        evanescent=abs(abs(z)-1)>1e-7
        q=v[:,evanescent]@amplitudes[evanescent]
        qnext=v[:,evanescent]@(z[evanescent]*amplitudes[evanescent])
        near+=scale*np.vdot(q,twist@qnext).imag
        for p in np.nonzero(~evanescent)[0]:
            current=scale*abs(amplitudes[p])**2*(z[p]*np.vdot(v[:,p],twist@v[:,p])).imag
            outward+=max(0,current)
            inward-=min(0,current)
    return outward,inward,near


def calculate(args):
    point,_,bend,flow,archive=settings()
    params=archive["parameters"]
    rate=0. if args.rest else float(flow["reference_inward_current"])/params["m"]
    R=archive["phase_A"]["q_star"][0]
    d=ApertureTrace(R,float(point["L"])/2,R,args.cells,2.,2.,3.,params)
    d.set_current(rate)
    for tone in ("difference","sum"):
        omega=float(bend[tone])
        _,_,feed,bulk=d.response(omega)
        for side,flux in (("brane",feed),("bulk",bulk)):
            incident,returned,evanescent=project(d,omega,rate,side)
            residual=incident-returned+evanescent-flux
            print(f"aperture={side} tone={tone} cells={args.cells} incident_per_b2={incident:.12g} "
                  f"returned_per_b2={returned:.12g} near_field_per_b2={evanescent:.12g} "
                  f"net_per_b2={flux:.12g} return_over_incident={returned/incident:.12g} residual={residual:.3g}",flush=True)
            assert abs(residual)<1e-7*max(flux,1) and incident>0


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cells",type=int,required=True)
    p.add_argument("--rest",action="store_true")
    calculate(p.parse_args())
