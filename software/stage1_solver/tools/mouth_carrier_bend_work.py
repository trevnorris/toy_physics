"""Linear carrier work at a prescribed, translating reflecting sidewall.

One open sleeve and finite receiver have the same stationary fixed-current
solution as the bend control. The two real carrier perturbations have angular
factor cos(theta) and exp(-i omega t) phasors. A translating Neumann wall
contributes the boundary variation of the rotating-frame carrier action.
Uniform leads carry excess wave energy away from the finite domain.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from stage1_solver.config import TensorGridSpec
from stage1_solver.open_sleeve import OpenSleeve

IDENTITY = np.eye(2)
J = np.array([[0.,-1.],[1.,0.]])


def rotate(phi):
    return IDENTITY*np.cos(phi) + J*np.sin(phi)


def radial_basis(cells, dr):
    radius = np.arange(cells+1)*dr
    area = 4*np.pi*np.diff(radius**3)/9  # sphere integral of cos(theta)^2
    link = (4*np.pi/3)*radius[1:-1]**2/dr
    diagonal = (np.r_[0,link]+np.r_[link,0]+8*np.pi*dr/3)/area
    vals, vec = eigh_tridiagonal(diagonal,-link/np.sqrt(area[1:]*area[:-1]))
    return area, vals, vec


def surface_step(lam, omega, solver, drift):
    """Outgoing two-component ghost value for a unit boundary-cell value.

    The rotating phase of a uniform discrete current changes by drift per
    outward cell. Evanescent roots decay; propagating roots carry positive
    excess energy outwards. The material reservoir has a matched DC chemical
    offset for this Bloch background, with the same density and compressibility.
    """
    k = solver.kinetic
    hop = k/solver.dw**2
    forward = rotate(drift)
    compressibility = 5*solver.eos_k*solver.bath_density**4
    onsite = ((2*hop*np.cos(drift)+k*lam)*IDENTITY
              + np.diag([2*compressibility,0.]) + 1j*solver.hbar*omega*J)
    companion = np.block([[np.zeros((2,2)), IDENTITY],
                          [-forward.T@forward.T, forward.T@onsite/hop]])
    roots, cols = np.linalg.eig(companion)
    subspace = []
    for z, v in zip(roots, cols[:2].T):
        flux = np.imag(z*np.vdot(v,forward@v))
        if abs(z) < 1-1e-7 or (abs(abs(z)-1)<=1e-7 and flux>0):
            subspace.append((z,v))
    if len(subspace)!=2:
        raise ArithmeticError(f"Ambiguous outgoing subspace: lambda={lam}, roots={roots}")
    z, columns = zip(*subspace)
    modes = np.array(columns).T
    return (modes*z)@np.linalg.inv(modes)


def lead(solver, indices, omega, phase, current_out):
    area, lam, basis = radial_basis(len(indices),solver.dr)
    # Match number flux to the mean outward phase gradient of the lead.
    delta = np.arcsin(current_out*solver.constituent_mass*solver.dw /
                      (3*sum(area)*solver.bath_density*solver.hbar))
    angle = phase-delta/2
    dimension = 2*len(indices)
    matrix = np.zeros((dimension,dimension),complex)
    weights = np.sqrt(area)[:,None]*basis
    for mode, eigenvalue in enumerate(lam):
        step = surface_step(eigenvalue,omega,solver,delta)
        ghost = rotate(angle+delta)@step@rotate(-angle)
        matrix += np.kron(IDENTITY-ghost,np.outer(weights[:,mode],weights[:,mode]))
    return (solver.kinetic/solver.dw)*matrix


def wall_source(solver, psi, length, order):
    """Stationary wall action variation b(s.u+t.udot), per max displacement."""
    nc = len(psi)
    w = solver.sleeve_cells
    side = np.nonzero((solver.r_index==solver.mouth_cells-1)&(solver.w_index<w))[0]
    inner = np.nonzero((solver.r_index==solver.mouth_cells-2)&(solver.w_index<w))[0]
    assert np.array_equal(solver.w_index[side],np.arange(w))
    # Smooth stationary Neumann trace by even extrapolation in radius.
    weights = np.array([1.,0.]) if order==1 else np.array([9/8,-1/8])
    trace = weights[0]*psi[side]+weights[1]*psi[inner]
    f = np.sin(2*np.pi*(np.arange(w)+.5)*solver.dw/length)**2
    f_edge = np.sin(2*np.pi*np.arange(1,w)*solver.dw/length)**2
    area = 4*np.pi*(solver.mouth_cells*solver.dr)**2/3
    mu = 1.25*solver.eos_k*solver.bath_density**4
    enthalpy = 1.25*solver.eos_k*abs(trace)**8
    spatial = 2*area*solver.dw*f*(mu-enthalpy)*trace
    flux = 2*solver.kinetic*area*f_edge*np.diff(trace)/solver.dw
    spatial[:-1] += flux
    spatial[1:] -= flux
    temporal = -1j*solver.hbar*area*solver.dw*f*trace
    vectors=[]
    for value in (spatial, temporal):
        profile = np.zeros(nc,complex)
        profile[side] = weights[0]*value
        profile[inner] = weights[1]*value
        vectors.append(np.r_[profile.real,profile.imag])
    return vectors


def linear_operator(solver, psi, omega):
    vol = solver.volume/3
    n = abs(psi)**2
    hc = 1.25*solver.eos_k
    # Half the rotating-frame energy Hessian in real quadratures.
    h = (solver.operator/3 + sparse.diags((hc*(n**4-solver.bath_density**4))*vol
         + 8*np.pi*solver.kinetic*solver.dr*solver.dw/3*np.ones(len(psi))))
    stiff = 8*hc*n**3*vol
    cross = stiff*psi.real*psi.imag
    symp = 1j*omega*solver.hbar*vol
    return sparse.bmat([[h+sparse.diags(stiff*psi.real**2),sparse.diags(cross-symp)],
                        [sparse.diags(cross+symp),h+sparse.diags(stiff*psi.imag**2)]],format="csc")


def response(solver, psi, phi, rate, length, omega, ports, drift, trace_order):
    matrix = linear_operator(solver,psi,omega)
    nc = len(psi)
    at_port = []
    for name, z, coef, angle, flow in (("feed",0,solver.inlet,phi,-rate),
            ("receiver",solver.grid.nw-1,solver.outlet,0.,rate)):
        if name not in ports:
            continue
        ids = np.nonzero(solver.w_index==z)[0]
        idx = np.r_[ids,ids+nc]
        outgoing = lead(solver,ids,omega,angle,flow if drift else 0.)
        update = outgoing-np.diag(np.tile(solver.kinetic*coef[ids]/3,2))
        matrix += sparse.coo_matrix((update.ravel(),(np.repeat(idx,len(idx)),np.tile(idx,len(idx)))),
                                    shape=matrix.shape).tocsc()
        at_port.append((name,idx,outgoing))
    spatial, temporal = wall_source(solver,psi,length,trace_order)
    source = spatial+1j*omega*temporal
    # Unit displacement; multiply work and fluxes by physical amplitude squared.
    u = spsolve(matrix,source/2)
    traction = np.vdot(source,u)
    # Integrate the physical pressure variation independently on the wall.
    side=np.nonzero((solver.r_index==solver.mouth_cells-1)&(solver.w_index<solver.sleeve_cells))[0]
    inner=np.nonzero((solver.r_index==solver.mouth_cells-2)&(solver.w_index<solver.sleeve_cells))[0]
    coeff=0. if trace_order==1 else .125
    def trace(field):
        return (1+coeff)*field[side]-coeff*field[inner]
    base=trace(psi)
    a,d=trace(u[:nc]),trace(u[nc:])
    f=np.sin(2*np.pi*(solver.w_index[side]+.5)*solver.dw/length)**2
    fe=np.sin(2*np.pi*(solver.w_index[side][1:])*solver.dw/length)**2
    area=4*np.pi*(solver.mouth_cells*solver.dr)**2/3
    pressure=(1j*omega*solver.hbar*(base.real*d-base.imag*a) +
              2*1.25*solver.eos_k*(solver.bath_density**4-abs(base)**8)*(base.real*a+base.imag*d))
    direct=area*(solver.dw*sum(f*pressure)-2*solver.kinetic/solver.dw*
                 sum(fe*(np.diff(base.real)*np.diff(a)+np.diff(base.imag)*np.diff(d))))
    assert abs(direct-traction)<2e-9*max(1,abs(traction))
    work=.5*omega*direct.imag
    powers={"feed":0., "receiver":0.}
    for name, idx, terminator in at_port:
        powers[name]=-omega*np.vdot(u[idx],terminator@u[idx]).imag
    cut=(solver.grid.nw+solver.sleeve_cells)//2
    left=np.nonzero(solver.w_index==cut-1)[0]
    right=np.nonzero(solver.w_index==cut)[0]
    area=radial_basis(len(left),solver.dr)[0]
    cut_flux=omega*solver.kinetic/solver.dw*np.sum(area*(u[left].conj()*u[right]+u[left+nc].conj()*u[right+nc]).imag)
    scale=max(1,abs(work))
    assert abs(work-sum(powers.values()))<1e-8*scale
    assert abs(cut_flux-powers["receiver"])<1e-8*scale
    assert min(powers.values())>=-1e-9*scale
    return work, powers, cut_flux, direct


def main(args):
    archive=yaml.safe_load(args.archive.read_text())
    p=archive["parameters"]
    radius=archive["phase_A"]["q_star"][0]
    runs=json.loads(args.spectrum.read_text())["runs"]
    def fields(run, tag):
        line,=[line for line in runs[run]["stdout_lines"] if line.startswith(tag)]
        return dict(pair.split("=",1) for pair in line.split() if "=" in pair)
    length=float(fields("smooth_dipole_l015_01","point n=1 ")["L"])
    rate=float(fields("mouth_wave_depth_response01","archive ")["reference_inward_current"])/p["m"]
    bend=fields("mouth_wave_linewidth_budget01","bow h=0.015 ")
    tone=[float(bend[key]) for key in ("difference","sum")]
    cells=args.cells
    grid=TensorGridSpec(2*radius,2*cells,0.,length*(1+args.receiver_lengths)/2,(1+args.receiver_lengths)*cells)
    solver=OpenSleeve(grid,cells,cells,hbar=p["hbar"],constituent_mass=p["m"],eos_k=p["K_eos"],bath_density=p["rho0"]/p["m"])
    flows=[rate,0.] if args.controls else [rate]
    for flow in flows:
        psi,phase,res=solver.solve(flow)
        entry,thru,far,err=solver.budget(psi,phase)
        assert max(abs(np.array([entry,thru,far])-flow))<1e-10 and err<1e-8
        print(f"carrier cells={cells} receiver_lengths={args.receiver_lengths} current={flow:.12g} inlet={entry:.12g} exit={far:.12g} residual={res:.3g}",flush=True)
        cases=[("feed+receiver",True,2),("receiver",True,2)]
        if args.controls and flow:
            cases.extend([("none",True,2),("feed+receiver",False,2),("feed+receiver",True,1)])
        for ports,drift,order in cases:
            for omega in tone:
                work,powers,cut,traction=response(solver,psi,phase,flow,length,omega,ports,drift,order)
                print(f"work cells={cells} current={flow:.12g} ports={ports} drift={int(drift)} trace={order} omega={omega:.12g} "
                      f"work_per_b2={work:.12g} feed_per_b2={powers['feed']:.12g} receiver_per_b2={powers['receiver']:.12g} "
                      f"cut_per_b2={cut:.12g} effective_C={2*work/omega**2:.12g} reactive_F_per_b={traction.real:.12g}",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,required=True)
    parser.add_argument("--receiver-lengths",type=int,default=1)
    parser.add_argument("--controls",action="store_true")
    parser.add_argument("--archive",type=Path,default=Path("software/stage1_solver/reports/pathA_26_derrick_results.yaml"))
    parser.add_argument("--spectrum",type=Path,default=Path("research/pde_ledger_v3/_measurements/mouth_trapped_wave_spectrum.json"))
    main(parser.parse_args())
