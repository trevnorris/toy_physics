"""Scalar support response to a bend for the two retained wall sidebands.

Leave the selected trapped dipole explicit at the central frequency. The
lower and upper scalar satellites see all radial and axial dipole modes.
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, splu, spsolve
from scipy.special import expit

from mouth_carrier_bend_budget import settings
from mouth_smooth_dipole_relaxation import AxialModes
from mouth_wave_outgoing_sidebands import gate


def calculate(step,extent,strength):
    point,_,_,_,archive = settings()
    p = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    length = float(point["L"])
    mesh = AxialModes(radius,p,step,round(7.2/step)*step)
    pot,_,n = gate(mesh,length,step,p)
    h0 = mesh.kinetic+sparse.diags(pot.ravel())
    h1 = h0+sparse.diags(2*mesh.centrifugal)
    def lowest(h):
        lam,v = eigsh(h,k=1,sigma=0.,which="LM",v0=mesh.seed,tol=1e-11)
        return lam[0],v[:,0]
    ls,s = lowest(h0)
    _,d = lowest(h1)
    ws = np.sqrt(ls)
    scale = p["I_wave"]/(2*ws)
    f = np.where(mesh.w<length/2,np.sin(2*np.pi*mesh.w/length)**2,0.)
    axial = expit((length-2*np.sqrt(mesh.w**2+p["eps_w"]**2))/p["delta_parallel"])
    grad = mesh.depth*(2/p["delta_r"])*mesh.gr*(1-mesh.gr)
    source = -s*np.outer(grad,f*axial).ravel()/np.sqrt(3)
    overlap = d@source
    rest = source-d*overlap
    identity = sparse.eye(h1.shape[0],format="csc")
    stat = splu((h1-ls*identity).tocsc()).solve
    first = stat(source)
    mw = 2*scale*(first@first+4*ls*first@stat(first))
    print(f"adiabatic step={step:g} extent={extent:g} strength={strength:g} "
          f"overlap={abs(overlap):.12g} slow_wave_mass={mw:.12g}",flush=True)
    outer = AxialModes(radius,p,step,extent)
    pout,_,m = gate(outer,length,step,p)
    rout = step*np.arange(1,m)
    def absorber(x):
        return strength*(np.maximum(x-5,0)/(extent-5))**3
    cap = absorber(rout)[:,None]+absorber(outer.w)[None,:]
    ho = outer.kinetic+sparse.diags(pout.ravel()+2*outer.centrifugal-1j*cap.ravel())
    eye_out = sparse.eye(ho.shape[0],format="csc")
    src_out = np.zeros_like(pout)
    src_out[:n-1,:n] = source.reshape(n-1,n)
    # The continuum argument is absolute frequency anchored at the converged
    # reference. The explicit selected level is also taken from that reference.
    wd_ref, ws_ref = float(point["omega_d"]),float(point["omega_s"])
    for delta in (-.016,0.,.008,.016,.024,.032,.04,.05,.065,.08):
        z = wd_ref+delta
        # These meshes have eigenfrequencies distinct from the reference.
        # Impose d.v=0 by a scalar Schur complement.
        central = h1-z*z*identity
        fac = splu(central.tocsc())
        v = fac.solve(rest)
        correction = fac.solve(d)
        v -= correction*((d@v)/(d@correction))
        assert abs(d@v)<1e-9
        central_load = scale*(rest@v)
        low = scale*(source@spsolve(h1-(z-2*ws_ref)**2*identity,source))
        out = spsolve(ho-(z+2*ws_ref)**2*eye_out,src_out.ravel())
        high = scale*(src_out.ravel()@out)
        dissipated = scale*np.dot(cap.ravel(),abs(out)**2)
        assert high.imag>0 and abs(dissipated-high.imag)<1e-8
        print(f"wave offset={delta:g} central={central_load:.12g} lower={low:.12g} "
              f"upper_real={high.real:.12g} upper_imag={high.imag:.12g}",flush=True)


if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step",type=float,required=True)
    ap.add_argument("--extent",type=float,required=True)
    ap.add_argument("--strength",type=float,required=True)
    args = ap.parse_args()
    calculate(args.step,args.extent,args.strength)
