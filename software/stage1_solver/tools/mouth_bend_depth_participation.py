"""Resolve the fixed-intake carrier response in depth and bound a local sheet law."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh, splu
from scipy.special import expit

from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain
from mouth_projected_bend_spectrum import fields, spectrum
from mouth_smooth_dipole_relaxation import AxialModes
from mouth_wave_outgoing_sidebands import gate

FOLDER = Path("research/pde_ledger_v3/_measurements")


def carrier(cells, outer):
    point,_,_,flow,archive = settings()
    params = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    half = float(point["L"])/2
    domain = ReservoirDomain(radius,half,radius,cells,outer-2,2,3,params)
    domain.set_current(float(flow["reference_inward_current"])/params["m"])
    # Same finite, real bath as the zero-frequency inertia calibration.
    volume,links,inlet,outlet,angular,lap = domain.metric(domain.re,domain.ze)
    weight = volume/3
    p = domain.psi
    density = abs(p)**2
    hc = 1.25*domain.eos_k
    potential = hc*(density**4-domain.bath_density**4)*weight
    h = lap/3+sparse.diags(domain.kinetic*angular+potential)
    slope = 8*hc*density**3*weight
    xy = sparse.diags(slope*p.real*p.imag)
    a = sparse.bmat([[h+sparse.diags(slope*p.real**2),xy],
                     [xy,h+sparse.diags(slope*p.imag**2)]],format="csc")
    solve = splu(a).solve
    static = solve(domain.spatial)
    x,y = np.split(static,2)
    rhs = domain.temporal-domain.hbar*np.r_[-weight*y, weight*x]
    response = solve(rhs)
    x,y = np.split(response,2)
    direct_mass = .5*(rhs@response)

    # Local allocation of 1/2 response^T A response. Split a face-gradient
    # term equally between the adjacent cells. Sum includes both reservoirs.
    amplitude2 = x*x+y*y
    local = .5*((domain.kinetic*(angular+(inlet+outlet)/3)+potential)*amplitude2
                +slope*(p.real*x+p.imag*y)**2)
    l,r = domain.left,domain.right
    face_term = domain.kinetic*links*((x[l]-x[r])**2+(y[l]-y[r])**2)/6
    np.add.at(local,l,face_term/2)
    np.add.at(local,r,face_term/2)
    assert np.isclose(local.sum(),direct_mass,rtol=1e-10)
    z = domain.z_index
    feed_mass = local[z<domain.nfeed].sum()
    bulk_mass = local[z>=domain.neck_end].sum()
    neck = (z>=domain.nfeed)&(z<domain.neck_end)
    by_depth = np.bincount(z[neck]-domain.nfeed,weights=local[neck],minlength=cells)

    # Eulerian velocity component in phase with b_dot after removing the
    # static shape variation. For exp(-i omega t), d phase / d b_dot has
    # a minus sign and the source in the quadrature equation is divided by 2.
    phi = -domain.hbar*(p.real*y-p.imag*x)/(2*domain.constituent_mass*density)
    vpot = np.empty((cells,cells))
    vpot[domain.r_index[neck],z[neck]-domain.nfeed] = phi[neck]
    # At the moving reflecting wall v_r also contains advection of the bent
    # surface by the *steady* axial stream. That term is in phase with b and
    # drops out of the b_dot derivative. Extrapolate the potential linearly.
    wall_phi = 1.5*vpot[-1]-.5*vpot[-2]
    # For Phi=phi(r,w)cos(theta), integral_B d_x Phi / volume = phi(R,w)/R.
    # This is the volume-mean transverse carrier velocity, per peak wall speed.
    vmean = wall_phi/radius
    normalized_depth = (np.arange(cells)+.5)/cells
    for depth,speed,mass in zip(normalized_depth,vmean,by_depth):
        print(f"slice w_over_H={depth:.10g} mean_vx_over_bdot={speed:.12g} "
              f"carrier_added_mass_in_slice={mass:.12g}")
    for i,ids in enumerate(np.array_split(np.arange(cells),4)):
        print(f"quarter index={i+1} mean_vx_over_bdot={vmean[ids].mean():.12g} "
              f"carrier_added_mass={by_depth[ids].sum():.12g}")
    print(f"carrier cells={cells} outer_R={outer:g} slab={feed_mass:.12g} "
          f"neck={sum(by_depth):.12g} bulk={bulk_mass:.12g} total={direct_mass:.12g}")


def sheet_bound():
    old = json.loads((FOLDER/"mouth_bend_coefficients.json").read_text())["runs"]
    geo = fields(old["bend_projection_wave015"]["stdout_lines"][0])
    breathing = float(geo["breathing"])
    threshold = float(fields(old["bend_projection_fold64"]["stdout_lines"][0])["wall_mass_at_fold"])
    # Every nonnegative nonincreasing line density is a positive mixture of
    # top hats [0,t]. Its normalized f² moment is a weighted average of theirs.
    def integral(t):
        return 3*t/8-np.sin(2*np.pi*t)/(4*np.pi)+np.sin(4*np.pi*t)/(32*np.pi)
    def f2(t):
        return np.sin(np.pi*t)**4
    depth = brentq(lambda t:t*f2(t)-integral(t),.5,1.)
    moment_max = integral(depth)/depth
    mass_max = breathing*moment_max/3
    required = 3*threshold/breathing
    # At most 1/4 of the maximal weight lies at any point outside the middle
    # half of the arm. This yields a necessary bound on mass inside that half.
    min_middle_share = (required-.25)/.75
    assert np.isclose(integral(1),.375) and moment_max<required<1
    print(f"sheet peak_threshold_mass={threshold:.12g} required_f2_moment={required:.12g} "
          f"monotone_optimal_depth={depth:.12g} monotone_max_moment={moment_max:.12g} "
          f"monotone_max_mass={mass_max:.12g} necessary_middle_half_share={min_middle_share:.12g}")
    spectral_readout(mass_max,"maximal_monotone")


def spectral_readout(mass,label):
    old = json.loads((FOLDER/"mouth_bend_coefficients.json").read_text())["runs"]
    geo,wave = map(fields,old["bend_projection_wave015"]["stdout_lines"])
    point,level,bend,_,arch = settings()
    ws,wd = float(point["omega_s"]),float(point["omega_d"])
    g = float(bend["overlap"])*np.sqrt(2*arch["parameters"]["I_wave"]/ws)
    gap = abs(float(level["U_over_action"]))
    prior = json.loads((FOLDER/"mouth_open_reservoirs.json").read_text())["runs"]
    for cells,key,ext,wkey in [(48,"open_reservoir_lines48","rad48","wresp05"),
                               (64,"final_open_lines64","rad64","wresp03")]:
        rows = prior[key]["stdout_lines"]+old["bend_projection_"+ext]["stdout_lines"]
        wr = [fields(s) for s in old["bend_projection_"+wkey]["stdout_lines"] if s.startswith("wave ")]
        offsets = [float(r["offset"]) for r in wr]
        kernels = [PchipInterpolator(offsets,[float(r[k]) for r in wr],extrapolate=False)
                   for k in ("central","lower","upper_real","upper_imag")]
        direct = float(wave["direct_K"])
        shift,width,asym,count = spectrum(rows,ws,wd,g,gap,mass,float(geo["tensile_both"])+direct,
                                          (direct,*kernels,1))
        print(f"{label}_spectrum cells={cells} shift={shift:.12g} width={width:.12g} "
              f"width_over_gap={width/gap:.12g} asymmetry={asym:.12g} peaks={count}")


def pressure_attached(step):
    # Hypothesis: the inertia of additional, normally entrained material is
    # proportional to the prepared cycle-averaged normal pressure. Only its
    # radial normalization is imported; the bend profile is a prediction.
    point,_,_,_,arch = settings()
    p = arch["parameters"]
    radius = arch["phase_A"]["q_star"][0]
    half = float(point["L"])/2
    mesh = AxialModes(radius,p,step,round(7.2/step)*step)
    pot,_,n = gate(mesh,2*half,step,p)
    h = mesh.kinetic+sparse.diags(pot.ravel())
    def ground(operator):
        lam,vec = eigsh(operator,k=1,sigma=0.,which="LM",v0=mesh.seed,tol=1e-11)
        return np.sqrt(lam[0]),vec[:,0]
    ws,s = ground(h)
    axial = expit((2*half-2*np.sqrt(mesh.w**2+p["eps_w"]**2))/p["delta_parallel"])
    grad = (2*mesh.depth/p["delta_r"])*mesh.gr*(1-mesh.gr)
    cell_load = p["I_wave"]/(2*ws)*(grad@(s.reshape(n-1,n)**2))*axial
    # The eigenvector on the even half domain carries unit Euclidean norm;
    # these cell-integrated reactions are for both mirrored halves combined.
    left = np.maximum(mesh.w-step/2,0)
    right = mesh.w+step/2
    dz = np.maximum(0,np.minimum(right,half)-left)
    f2 = np.where(mesh.w<half,np.sin(np.pi*mesh.w/half)**4,0.)
    # Check the *local* force by probing a nonuniform radial dilation f².
    radial = step*np.arange(1,n)
    def energy(shift):
        gate_r = expit(2*(radius+shift*f2[None,:]-radial[:,None])/p["delta_r"])
        gate_v = p["mu_out"]**2-mesh.depth*gate_r*axial[None,:]
        return p["I_wave"]*ground(mesh.kinetic+sparse.diags(gate_v.ravel()))[0]
    epsilon = 2e-4
    oracle = (energy(-epsilon)-energy(epsilon))/(2*epsilon)
    assert np.isclose(oracle,cell_load@f2,rtol=3e-6)
    total = cell_load@(dz/step)
    normalized = cell_load/total
    moment = normalized@(dz/step*f2)
    radial_mass = p["rho0"]*(4*np.pi/3)*radius**3*2*half
    mass = radial_mass*moment/3
    shares = []
    for i in range(4):
        widths = np.maximum(0,np.minimum(right,(i+1)*half/4)-np.maximum(left,i*half/4))
        shares.append(normalized@(widths/step))
    assert np.isclose(sum(shares),1) and np.all(cell_load>=0)
    print(f"pressure step={step:g} radial_force={total:.12g} f2_weighted_force={oracle:.12g} "
          f"moment={moment:.12g} wall_mass={mass:.12g} "
          f"quarter_shares={','.join(f'{v:.9g}' for v in shares)}")
    for w,frac,width in zip(mesh.w,normalized,dz):
        if width>0:
            print(f"attached w_over_H={(w-step/2+width/2)/half:.10g} "
                  f"radial_mass_share={frac*width/step:.12g}")
    spectral_readout(mass,"pressure_attached")


def draw(carrier_path,pressure_path,destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def rows(path,tag):
        return [fields(s) for s in path.read_text().splitlines() if s.startswith(tag+" ")]
    slices = rows(carrier_path,"slice")
    points = rows(pressure_path,"attached")
    x = [float(r["w_over_H"]) for r in slices]
    mean_velocity = [float(r["mean_vx_over_bdot"]) for r in slices]
    # A cell's reported x is the midpoint of its overlap with the physical arm.
    edges = [0.]
    for row in points:
        edges.append(2*float(row["w_over_H"])-edges[-1])
    assert abs(edges[-1]-1)<1e-7
    line_density = np.array([float(r["radial_mass_share"]) for r in points])/np.diff(edges)

    fig,axes = plt.subplots(1,2,figsize=(9.0,3.5),layout="constrained")
    grid = np.linspace(0,1,301)
    axes[0].plot(grid,np.sin(np.pi*grid)**2,color="0.45",ls="--",label="Prescribed wall motion")
    axes[0].plot(x,mean_velocity,color="#126878",lw=2.5,label="Carrier, section mean")
    axes[0].set_title("Sideways velocity at fixed intake")
    axes[0].set_ylabel("Fraction of peak wall speed")
    axes[0].legend(frameon=False,fontsize=9)
    axes[1].stairs(line_density,edges,color="#934714",lw=2,baseline=None)
    axes[1].axhline(1,color="0.55",lw=1,ls="--",label="Uniform attachment")
    axes[1].set_title("If attachment follows wave pressure")
    axes[1].set_ylabel("Attached mass per depth / arm mean")
    axes[1].legend(frameon=False,fontsize=9)
    for ax in axes:
        ax.set(xlim=(0,1),ylim=(0,None),xlabel="Depth along the draining arm, w/H",
               xticks=[0,.25,.5,.75,1])
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.18)
    destination.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(destination,dpi=160)
    plt.close(fig)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="mode",required=True)
    c = subs.add_parser("carrier")
    c.add_argument("--cells",type=int,required=True)
    c.add_argument("--outer",type=int,required=True)
    subs.add_parser("sheet")
    p = subs.add_parser("pressure")
    p.add_argument("--step",type=float,required=True)
    plot = subs.add_parser("plot")
    for arg in ("carrier","pressure","destination"):
        plot.add_argument(arg,type=Path)
    args = parser.parse_args()
    if args.mode=="carrier":
        carrier(args.cells,args.outer)
    elif args.mode=="pressure":
        pressure_attached(args.step)
    elif args.mode=="plot":
        draw(args.carrier,args.pressure,args.destination)
    else:
        sheet_bound()
