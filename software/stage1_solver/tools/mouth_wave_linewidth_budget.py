"""Project an anchored bend onto the smooth modes and budget measured wall work.

The optional single-coordinate fast response is an extrapolation of the
archived mechanical parameters. All other rate coefficients are kept in
terms of the independently measurable displacement response of this shape.
"""
import argparse
from pathlib import Path
import numpy as np
import yaml
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from scipy.special import expit

from mouth_smooth_dipole_relaxation import ARCHIVE, AxialModes
from mouth_wave_outgoing_sidebands import gate


def row(path, prefix):
    line=next(line for line in path.read_text().splitlines() if line.startswith(prefix))
    return dict(item.split("=",1) for item in line.split()[1:] if "=" in item)


def odd_overlap(p,radius,length,step):
    mesh = AxialModes(radius,p,step,round(7.2/step)*step)
    pot, _, n = gate(mesh,length,step,p)
    gw = expit((length-2*np.sqrt(mesh.w**2+p["eps_w"]**2))/p["delta_parallel"])
    axial_bow = np.where(mesh.w<length/2, np.sin(2*np.pi*mesh.w/length)**2, 0.)
    # Use a side displacement beta * cos(theta) * axial_bow on both ends.
    derivative = -mesh.depth*np.outer((2/p["delta_r"])*mesh.gr*(1-mesh.gr),gw*axial_bow)
    vectors=[]
    frequencies=[]
    for ell in (0,1):
        h = mesh.kinetic + diags(pot.ravel()+ell*(ell+1)*mesh.centrifugal)
        val,vec=eigsh(h,k=1,sigma=0,which="LM",v0=mesh.seed,tol=1e-11)
        vectors.append(vec[:,0])
        frequencies.append(np.sqrt(val[0]))
    overlap=abs(np.vdot(vectors[0],derivative.ravel()*vectors[1]))/np.sqrt(3)
    return overlap, frequencies


def main(ladder, sideband):
    cfg=yaml.safe_load(ARCHIVE.read_text())
    p=cfg["parameters"]
    point=row(ladder,"point n=1 ")
    levels=row(ladder,"levels ")
    high=row(sideband,"sideband kind=high ")
    radius=cfg["phase_A"]["q_star"][0]
    length=float(point["L"])
    width=abs(float(levels["U_over_action"]))
    action=p["I_wave"]
    q=p["hbar"]
    occupation=action/q
    for step in (.06,.03,.015):
        k,(ws,wd)=odd_overlap(p,radius,length,step)
        f=q*k/(2*np.sqrt(ws*wd))
        coeff=2*f*f/q
        print(f"bow h={step:.4g} overlap={k:.12g} force_matrix={f:.12g} "
              f"difference={wd-ws:.12g} sum={wd+ws:.12g} "
              f"loss_per_ImChi_difference={coeff*(occupation+1):.12g} "
              f"loss_per_ImChi_sum={coeff*occupation:.12g} "
              f"equal_ImChi_one_spacing={width/(coeff*(2*occupation+1)):.12g}",flush=True)
    ws,wd=map(float,(point["omega_s"],point["omega_d"]))
    omega=2*ws
    area=4*np.pi*radius**3/3
    mass=p["rho0"]*area*length/4
    k_static=float(levels["support_stiffness_at_one"])
    force=action*abs(float(point["domega_s_dL"]))
    power_per_unit=float(high["power_per_action_eps2"])
    # This coefficient is evaluated on the finest outgoing grid, with its
    # own sideband frequency; retain the same normalization for the rate.
    nu=float(high["out_frequency"])
    for damper in (0.,p["c_damp_sample"]):
        amplitude=force/abs(k_static-mass*omega**2-1j*damper*omega)
        gamma=power_per_unit*amplitude**2/nu
        print(f"coherent_sample damping={damper:g} mass={mass:.12g} stiffness={k_static:.12g} "
              f"force_amplitude={force:.12g} drive_frequency={omega:.12g} "
              f"depth_amplitude={amplitude:.12g} power_per_dipole_action={power_per_unit*amplitude**2:.12g} "
              f"upper_action_loss_rate={gamma:.12g} ratio_to_spacing={gamma/width:.12g}")
        assert amplitude < .01*length and gamma > 0
    length_slope=abs(float(point["domega_d_dL"]))
    print(f"target linewidth={width:.12g} inverse_width={1/width:.12g} "
          f"rms_depth_for_one_spacing={width/length_slope:.12g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("ladder",type=Path)
    parser.add_argument("upper_sideband_results",type=Path)
    args=parser.parse_args()
    main(args.ladder,args.upper_sideband_results)
