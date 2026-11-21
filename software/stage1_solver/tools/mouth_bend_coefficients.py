"""Restoring energy, projected inertia and carrier calibration for the sin² bend."""
import argparse
import math

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.sparse.linalg import eigsh, splu, spsolve
from scipy.special import expit

from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain
from mouth_smooth_dipole_relaxation import AxialModes
from mouth_wave_outgoing_sidebands import gate


def project(radius, half, params, step):
    # Both sides of the symmetric scalar well follow the same bend. Only one
    # side has a draining carrier. b is the peak displacement of each side.
    area = 4*np.pi*radius**2/3
    tension = 2*area*params["sigma"]*np.pi**2/(2*half)
    # A direct second difference of the exact translated-section area.
    def energy_increment(b):
        def angular(w):
            slope = np.pi/half*np.sin(2*np.pi*w/half)
            def increment(c):
                s = (b*slope*c)**2
                return s/(np.sqrt(1+s)+1)
            return quad(increment, -1, 1)[0]
        return 2*params["sigma"]*2*np.pi*radius**2*quad(angular, 0, half)[0]
    assert np.isclose(2*energy_increment(.0001)/.0001**2, tension, rtol=1e-6)
    # Conditional normal sheet inertia: all the archived breathing mass is
    # allocated to this sheet, with no cap contribution to radial breathing.
    breathing = params["rho0"]*4*np.pi*radius**3*(2*half)/3
    def moment(n):
        return math.comb(2*n,n)/4**n  # mean of sin^(2n)
    uniform = breathing*moment(2)/3
    # A smooth and strictly positive profile proportional to 0.02+f^4.
    central = breathing*(.02*moment(2)+moment(6))/(3*(.02+moment(4)))
    print(f"geometry radius={radius:.12g} half={half:.12g} tensile_one_arm={tension/2:.12g} "
          f"tensile_both={tension:.12g} breathing={breathing:.12g} "
          f"uniform_wall_mass={uniform:.12g} central_wall_mass={central:.12g}", flush=True)

    mesh = AxialModes(radius, params, step, round(7.2/step)*step)
    pot, _, n = gate(mesh, 2*half, step, params)
    radial = step*np.arange(1,n)
    axial = expit((2*half-2*np.sqrt(mesh.w**2+params["eps_w"]**2))/params["delta_parallel"])
    bend = np.where(mesh.w < half, np.sin(np.pi*mesh.w/half)**2, 0.)
    # Translate the smooth radial potential V(|x-b f e|, w).
    slope = mesh.depth*(2/params["delta_r"])*mesh.gr*(1-mesh.gr)
    curvature = (2/params["delta_r"])*(2*mesh.gr-1)*slope
    h0 = mesh.kinetic+sparse.diags(pot.ravel())
    h1 = h0+sparse.diags(2*mesh.centrifugal)
    def lowest(h):
        val, vec = eigsh(h, k=1, sigma=0., which="LM", v0=mesh.seed, tol=1e-11)
        return val[0], vec[:,0]
    lam, s = lowest(h0)
    lam_d, d = lowest(h1)
    omega_s = np.sqrt(lam)
    # Angular averages for Y00 are <cos²>=1/3 and <sin²>=2/3.
    direct_profile = np.outer((curvature+2*slope/radial)/3, axial*bend**2)
    direct = np.dot(s*s, direct_profile.ravel())
    source = -s*np.outer(slope, axial*bend).ravel()/np.sqrt(3)
    overlap = np.dot(d,source)
    h_shifted = h1-sparse.eye(h1.shape[0],format="csc")*lam
    response = spsolve(h_shifted,source)
    assert np.linalg.norm(h_shifted@response-source) < 1e-8
    relax = 2*np.dot(source,response)
    selected = 2*overlap**2/(lam_d-lam)
    scale = params["I_wave"]/(2*omega_s)
    print(f"support step={step:g} ws={omega_s:.12g} wd={np.sqrt(lam_d):.12g} "
          f"overlap={abs(overlap):.12g} direct_K={scale*direct:.12g} "
          f"relax_all={scale*relax:.12g} retained_relax={scale*selected:.12g} "
          f"static_wave_K={scale*(direct-relax):.12g} "
          f"background_wave_K={scale*(direct-relax+selected):.12g}", flush=True)


def carrier_mass(domain):
    # Finite real bath with far dipole perturbations fixed. Expand its force at
    # omega=0 before moving the far boundary out. No absorbing layer at DC.
    volume,_,_,_,angular,lap = domain.metric(domain.re,domain.ze)
    weight = volume/3
    p = domain.psi
    n = abs(p)**2
    hc = 1.25*domain.eos_k
    h = lap/3+sparse.diags(domain.kinetic*angular+hc*(n**4-domain.bath_density**4)*weight)
    c = 8*hc*n**3*weight
    xy = sparse.diags(c*p.real*p.imag)
    a = sparse.bmat([[h+sparse.diags(c*p.real**2),xy],
                     [xy,h+sparse.diags(c*p.imag**2)]],format="csc")
    solve = splu(a).solve
    s,t = domain.spatial,domain.temporal
    u = solve(s)
    x,y = np.split(u,2)
    corr = t-domain.hbar*np.r_[-weight*y,weight*x]
    added = .5*np.dot(corr,solve(corr))
    static = .5*np.dot(s,u)
    # Direct centered real-frequency susceptibility at fixed finite outer box.
    symp = sparse.diags(domain.hbar*weight)
    b = sparse.bmat([[None,-symp],[symp,None]],format="csc")
    probe = 1e-4
    src = s+1j*probe*t
    force = .5*np.vdot(src,spsolve(a+1j*probe*b,src))
    difference = (force.real-static)/probe**2
    assert abs(force.imag)<1e-12 and np.isclose(difference,added,rtol=2e-4)
    print(f"carrier cells={domain.mouth_cells} outer_R={domain.re[-1]/domain.radius:g} "
          f"static_induced={static:.12g} slow_added_mass={added:.12g} "
          f"finite_frequency_check={difference:.12g}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode",choices=("projection","carrier","extend"))
    ap.add_argument("--step",type=float,default=.03)
    ap.add_argument("--cells",type=int,default=32)
    ap.add_argument("--outer",type=float,default=4)
    args = ap.parse_args()
    point,_,bend,flow,archive = settings()
    p = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    half = float(point["L"])/2
    if args.mode=="projection":
        project(radius,half,p,args.step)
        return
    domain = ReservoirDomain(radius,half,radius,args.cells,args.outer-2,2,3,p)
    domain.set_current(float(flow["reference_inward_current"])/p["m"])
    if args.mode=="carrier":
        carrier_mass(domain)
        return
    # Positive offsets supplement the already measured negative offsets.
    for tone in ("difference","sum"):
        for delta in (.006,.009,.012,.015,.018,.021,.024,.028,.033,.04,.05,.065,.08):
            omega = float(bend[tone])+delta
            reaction,power,feed,bulk = domain.response(omega)
            assert min(power,feed,bulk)>0
            print(f"tone={tone} offset={delta:.8g} omega={omega:.12g} "
                  f"induced_real={reaction.real:.12g} power_per_b2={power:.12g}", flush=True)


if __name__=="__main__":
    main()
