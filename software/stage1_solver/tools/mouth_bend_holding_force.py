"""In-phase field loads for a held sin² bend, and a test of force identifiability.

The external force is conjugate to the common peak displacement b. Both well
halves are translated, one carrier arm drains, and all other imposed wall
harmonics vanish. No total holding-force sensor is assumed.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, spsolve
from scipy.special import expit

from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain
from mouth_smooth_dipole_relaxation import AxialModes
from mouth_wave_outgoing_sidebands import gate


TONES = (0., .01, .02)
AMPLITUDE = 1e-5
PREVIOUS = Path("research/pde_ledger_v3/_measurements/mouth_bend_coefficients.json")


def carrier_load(cells, outer):
    point,_,_,flow,archive = settings()
    p = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    half = float(point["L"])/2
    domain = ReservoirDomain(radius,half,radius,cells,2,4,3,p)
    # Preserve the neck, slab thickness and two real monitoring cuts. Grade
    # the far real cells to fit wavelengths much longer than the mouth radius.
    def extend(x):
        extra = np.maximum(x/radius-2,0)/4
        return x+(outer-6)*radius*extra**3
    domain.re = extend(domain.re)
    domain.ze = half+extend(domain.ze-half)
    metric = domain.metric(domain.re,domain.ze)
    domain.volume,domain.conductance,domain.inlet,domain.outlet,_,domain.operator = metric
    domain.set_current(float(flow["reference_inward_current"])/p["m"])
    def outgoing(x):
        start = outer*radius/4
        depth = outer*radius-start
        return x+1j*np.maximum(x-start,0)**3/depth**2
    ac = domain.metric(outgoing(domain.re),half+outgoing(domain.ze-half))
    for omega in TONES:
        domain.ac = metric if omega==0 else ac
        induced,power,feed,bulk = domain.response(omega)
        if omega:
            assert min(feed,bulk)>0 and abs(feed+bulk-power)<2e-7*power
        print(f"carrier cells={cells} outer_R={outer:g} omega={omega:g} "
              f"holding_per_b={-induced.real:.14g} holding_imag_per_b={-induced.imag:.14g} "
              f"power_per_b2={power:.12g} feed={feed:.12g} bulk={bulk:.12g}",flush=True)


def wave_load(step):
    point,_,_,_,archive = settings()
    p = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    length = float(point["L"])
    mesh = AxialModes(radius,p,step,round(7.2/step)*step)
    potential,_,n = gate(mesh,length,step,p)
    h0 = mesh.kinetic+sparse.diags(potential.ravel())
    h1 = h0+sparse.diags(2*mesh.centrifugal)
    eig,v = eigsh(h0,k=1,sigma=0.,which="LM",v0=mesh.seed,tol=1e-11)
    support = v[:,0]
    ws = np.sqrt(eig[0])
    scale = p["I_wave"]/(2*ws)
    f = np.where(mesh.w<length/2,np.sin(2*np.pi*mesh.w/length)**2,0.)
    axial = expit((length-2*np.sqrt(mesh.w**2+p["eps_w"]**2))/p["delta_parallel"])
    grad = mesh.depth*(2/p["delta_r"])*mesh.gr*(1-mesh.gr)
    curvature = (2/p["delta_r"])*(2*mesh.gr-1)*grad
    radial = step*np.arange(1,n)
    direct = scale*np.dot(support**2,np.outer((curvature+2*grad/radial)/3,axial*f**2).ravel())
    v1 = -support*np.outer(grad,axial*f).ravel()/np.sqrt(3)
    eye = sparse.eye(h1.shape[0],format="csc")
    def induced(nu):
        return scale*np.dot(v1,spsolve(h1-nu**2*eye,v1))
    for omega in TONES:
        # Only the pump-cycle-averaged, at-drive-frequency component enters
        # this clamped experiment. The pump also exerts a force at other tones.
        load = direct-induced(ws+omega)-induced(ws-omega)
        print(f"wave step={step:g} omega={omega:g} holding_per_b={load:.14g}",flush=True)


def fields(line):
    return dict(token.split("=",1) for token in line.split() if "=" in token)


def combine(carrier_path,wave_path):
    def read(path,prefix):
        rows = [fields(s) for s in path.read_text().splitlines() if s.startswith(prefix+" ")]
        assert [float(r["omega"]) for r in rows]==list(TONES)
        return np.array([float(r["holding_per_b"]) for r in rows])
    fluid = read(carrier_path,"carrier")
    wave = read(wave_path,"wave")
    known = fluid+wave
    prior = json.loads(PREVIOUS.read_text())["runs"]
    geo = fields(prior["bend_projection_wave015"]["stdout_lines"][0])
    threshold = fields(prior["bend_projection_fold64"]["stdout_lines"][0])
    k = float(geo["tensile_both"])
    design = np.column_stack((np.ones(3),-np.square(TONES)))
    assert np.linalg.matrix_rank(design)==2
    for i,omega in enumerate(TONES):
        print(f"correction b={AMPLITUDE:g} omega={omega:g} carrier_force={AMPLITUDE*fluid[i]:.12g} "
              f"wave_force={AMPLITUDE*wave[i]:.12g} combined_field_force={AMPLITUDE*known[i]:.12g}")
    # These are forward predictions using three explicitly supplied wall laws.
    # Fit only DC and the first tone, with the second held out. Independently
    # measured field data have zero sensitivity to either wall coefficient.
    for label,mass in (("uniform",float(geo["uniform_wall_mass"])),
                       ("spectral_fold",float(threshold["wall_mass_at_fold"])),
                       ("central",float(geo["central_wall_mass"]))):
        predicted = AMPLITUDE*(known+design@np.array([k,mass]))
        recovered = np.linalg.solve(design[:2],predicted[:2]/AMPLITUDE-known[:2])
        held_out = AMPLITUDE*(known[2]+design[2]@recovered)
        assert np.allclose(recovered,[k,mass],rtol=1e-10) and abs(held_out-predicted[2])<1e-15
        print(f"conditional sheet={label} F0={predicted[0]:.12g} F1={predicted[1]:.12g} "
              f"F2={predicted[2]:.12g} recovered_input_K={recovered[0]:.12g} "
              f"recovered_input_M={recovered[1]:.12g} F2_error={held_out-predicted[2]:.3g} "
              f"residual_drop1={AMPLITUDE*mass*TONES[1]**2:.12g} "
              f"residual_drop2={AMPLITUDE*mass*TONES[2]**2:.12g}")


if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode",required=True)
    c = sub.add_parser("carrier")
    c.add_argument("--cells",type=int,required=True)
    c.add_argument("--outer",type=float,required=True)
    w = sub.add_parser("wave")
    w.add_argument("--step",type=float,required=True)
    mix = sub.add_parser("combine")
    mix.add_argument("carrier",type=Path)
    mix.add_argument("wave",type=Path)
    args = ap.parse_args()
    if args.mode=="carrier":
        carrier_load(args.cells,args.outer)
    elif args.mode=="wave":
        wave_load(args.step)
    else:
        combine(args.carrier,args.wave)
