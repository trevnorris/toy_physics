"""First and second action increments of a dipole over an occupied monopole.

The scalar uses the symmetric finite potential gate from the freeze sheet.
The sharp analytic archived volume and tension terms are retained, as is
the radius. The length can move on a local equilibrium branch. We exclude
an implicit change of zero point: n=0,1,2 add n times the trial action step.
"""
import argparse
from functools import lru_cache

import numpy as np
import yaml
from scipy import sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh
from scipy.special import expit

from mouth_trapped_wave_ports import ARCHIVE


class AxialModes:
    def __init__(self, radius, params, step, extent):
        self.params = params
        cells = round(extent / step)
        assert np.isclose(cells * step, extent)
        r = step * (1 + np.arange(cells - 1))
        self.w = step * (.5 + np.arange(cells))
        radial = sparse.diags((-np.ones(cells-2),2*np.ones(cells-1),-np.ones(cells-2)),(-1,0,1))
        axial_diagonal = 2*np.ones(cells)
        axial_diagonal[0], axial_diagonal[-1] = 1, 3
        axial = sparse.diags((-np.ones(cells-1),axial_diagonal,-np.ones(cells-1)),(-1,0,1))
        speed2 = params["c_w"]**2
        self.kinetic = (speed2/step**2)*sparse.kronsum(axial, radial, format="csc")
        self.centrifugal = np.broadcast_to((speed2/r**2)[:,None], (len(r),cells)).ravel()
        self.gr = expit(2*(radius-r)/params["delta_r"])
        self.depth = params["mu_out"]**2-params["mu_in"]**2
        self.seed = (r[:,None]*np.exp(-r[:,None]**2-self.w[None,:]**2)).ravel()

    @lru_cache(maxsize=128)
    def eigendata(self, length, ell):
        p = self.params
        gw = expit((length-2*np.sqrt(self.w**2+p["eps_w"]**2))/p["delta_parallel"])
        gate = np.outer(self.gr, gw)
        potential = p["mu_out"]**2 - self.depth*gate.ravel() + ell*(ell+1)*self.centrifugal
        operator = self.kinetic + sparse.diags(potential, format="csc")
        # The shift is below all eigenvalues by positivity of the discrete K and V.
        values, vectors = eigsh(operator, k=1, sigma=0., which="LM", v0=self.seed, ncv=12, tol=1.e-11)
        omega = np.sqrt(values[0])
        vector = vectors[:,0]
        dgate = np.outer(self.gr, gw*(1-gw)/p["delta_parallel"])
        slope = -self.depth*np.dot(dgate.ravel(), vector**2)/(2*omega)
        residual = np.linalg.norm(operator@vector-values[0]*vector)
        return omega, slope, residual


def calculate(step, extent):
    archive = yaml.safe_load(ARCHIVE.read_text())
    p = archive["parameters"]
    radius, _ = archive["phase_A"]["q_star"]
    modes = AxialModes(radius, p, step, extent)
    cross_section = 4*np.pi*radius**3/3
    area_coefficient = 4*np.pi*radius**2
    cap_cost = 2*p["sigma"]*cross_section
    f_wall = p["alpha_volume"]*cross_section + p["sigma"]*area_coefficient
    action, increment = p["I_wave"],p["hbar"]

    def gradient(length, added):
        mono = modes.eigendata(length,0)[1]
        dip = modes.eigendata(length,1)[1] if added else 0.
        return f_wall + action*mono + added*dip

    def point(count):
        added = count*increment
        bracket = (1.8,2.5)
        g_left,g_right = (gradient(length,added) for length in bracket)
        if not g_left < 0 < g_right:
            bracket = (1.6,3.2)
            g_left,g_right = (gradient(length,added) for length in bracket)
        length = brentq(gradient,*bracket,args=(added,),xtol=1.e-10)
        mono,dm,err_m = modes.eigendata(length,0)
        dip,dd,err_d = modes.eigendata(length,1)
        energy = cap_cost + f_wall*length + action*mono + added*dip
        print(f"point n={count} step={step:g} extent={extent:g} L={length:.12g} "
              f"E={energy:.12g} omega_s={mono:.12g} omega_d={dip:.12g} "
              f"domega_s_dL={dm:.12g} domega_d_dL={dd:.12g} "
              f"cutoff_margin={p['mu_out']-dip:.10g} grad={gradient(length,added):.3g} "
              f"bracket_gradients={g_left:.6g},{g_right:.6g} eig_residual={max(err_m,err_d):.3g}",flush=True)
        assert max(err_m,err_d)<1.e-6
        assert mono < dip < p["mu_out"] and g_left<0<g_right
        return length,energy,dip,dd

    data = [point(n) for n in range(3)]
    lengths,energies,frequencies,slopes = np.array(data).T
    first,second = np.diff(energies)
    splitting = second-first
    # For a smooth local branch E_min'' = -omega_d'^2 / d_L^2 E.
    # This differentiation check is independent of differencing total energies.
    mid_length = lengths[1]
    delta_l = .0005
    stiff = (gradient(mid_length+delta_l,increment)-gradient(mid_length-delta_l,increment))/(2*delta_l)
    midpoint_prediction = -increment**2 * slopes[1]**2 / stiff
    print(f"levels step={step:g} extent={extent:g} rigid_spacing_at_empty={increment*frequencies[0]:.12g} "
          f"relaxed_first={first:.12g} relaxed_second={second:.12g} U={splitting:.12g} "
          f"U_over_action={splitting/increment:.12g} support_stiffness_at_one={stiff:.12g} "
          f"midpoint_U_prediction={midpoint_prediction:.12g} "
          f"prediction_residual={abs(splitting-midpoint_prediction):.5g}",flush=True)
    assert stiff>0 and np.all(np.diff(lengths)>0) and first>second>0
    for n in range(2):
        slope_secant = (energies[n+1]-energies[n])/increment
        assert frequencies[n+1]-2.e-8 <= slope_secant <= frequencies[n]+2.e-8


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=float, required=True)
    parser.add_argument("--extent", type=float, required=True)
    args = parser.parse_args()
    calculate(args.step,args.extent)
