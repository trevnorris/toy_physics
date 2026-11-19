"""Bulk scalar power radiated per imposed sinusoidal depth displacement squared.

Use the frozen symmetric smooth well at the relaxed one-increment geometry.
The reduced radial field is r A. Outer absorbers approximate outgoing waves
in both open coordinates. Independently compare interior flux, absorbed power
and work on the source; quote refinement sensitivity for the open response.
"""
import argparse
from pathlib import Path
import re

import numpy as np
import yaml
from scipy import sparse
from scipy.sparse.linalg import eigsh, spsolve
from scipy.special import expit

from mouth_smooth_dipole_relaxation import AxialModes
from mouth_trapped_wave_ports import ARCHIVE


def gate(mesh, length, step, p):
    n = len(mesh.w)
    gw = expit((length - 2*np.sqrt(mesh.w**2 + p["eps_w"]**2))/p["delta_parallel"])
    profile = np.outer(mesh.gr, gw)
    potential = p["mu_out"]**2 - mesh.depth*profile
    derivative = -mesh.depth * np.outer(mesh.gr,gw*(1-gw)/p["delta_parallel"])
    return potential, derivative, n


def prepared(radius, length, p, step):
    inner = AxialModes(radius,p,step,round(6/step)*step)
    pot, slope, n = gate(inner,length,step,p)
    modes = []
    for ell in (0,1):
        h = inner.kinetic + sparse.diags(pot.ravel()+ell*(ell+1)*inner.centrifugal)
        val, vec = eigsh(h,k=1,sigma=0.,which="LM",v0=inner.seed,tol=1e-11)
        modes.append((np.sqrt(val[0]),slope*vec[:,0].reshape(n-1,n)))
    return modes


def calculate(step, extent, absorber_strength, length, closed):
    config = yaml.safe_load(ARCHIVE.read_text())
    radius = config["phase_A"]["q_star"][0]
    p = config["parameters"]
    modes = prepared(radius,length,p,step)
    ws,wd = [row[0] for row in modes]
    outer = AxialModes(radius,p,step,extent)
    pot, _, n = gate(outer,length,step,p)
    r = step*np.arange(1,n)
    w = outer.w
    start = 6.
    def profile(x):
        return absorber_strength*(np.maximum(0., x-start)/(extent-start))**3
    absorber = profile(r)[:,None]+profile(w)[None,:]
    tests = [(1,wd+2*ws,"high"),(0,2*wd-ws,"low")]
    if closed:
        tests.append((0,1.,"closed"))
    results = []
    for ell,nu,label in tests:
        source = np.zeros_like(pot)
        win,patch = modes[ell]
        source[:patch.shape[0],:patch.shape[1]] = patch
        # Five-point nearest-neighbor stencil, with Neumann even symmetry at w=0.
        diagonal = pot.ravel() + ell*(ell+1)*outer.centrifugal - nu**2 - 1j*absorber.ravel()
        operator = outer.kinetic + sparse.diags(diagonal,format="csc")
        out = spsolve(operator,source.ravel()).reshape(pot.shape)
        work = np.vdot(source,out).imag
        sink = np.sum(abs(out)**2*absorber)
        # The cut excludes the absorber and contains essentially the entire source.
        radial_cells = np.count_nonzero(r<5.)
        axial_cells = np.count_nonzero(w<5.)
        radial_edge = np.vdot(out[radial_cells-1,:axial_cells],out[radial_cells,:axial_cells]).imag
        axial_edge = np.vdot(out[:radial_cells,axial_cells-1],out[:radial_cells,axial_cells]).imag
        flux = (radial_edge+axial_edge)*p["c_w"]**2/step**2
        print(f"sideband kind={label} ell={ell} h={step:g} extent={extent:g} absorber={absorber_strength:g} "
              f"depth={length:.12g} input_frequency={win:.12g} out_frequency={nu:.12g} "
              f"Im_response={work:.12g} shell_flux={flux:.12g} absorbed={sink:.12g} "
              f"power_per_action_eps2={nu/(4*win)*work:.12g}",flush=True)
        assert work>=-1e-13 and abs(work-sink)<2e-10 and abs(work-flux)<2e-10
        results.append(nu/(4*win)*work)
    return np.array(results)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("shape_result",type=Path)
    parser.add_argument("--step",type=float,required=True)
    parser.add_argument("--extent",type=float,required=True)
    parser.add_argument("--absorber",type=float,required=True)
    parser.add_argument("--closed-control",action="store_true")
    args=parser.parse_args()
    line=next(s for s in args.shape_result.read_text().splitlines() if "n=1 " in s)
    length=float(re.search(r"\bL=([0-9.eE+-]+)",line).group(1))
    calculate(args.step,args.extent,args.absorber,length,args.closed_control)
