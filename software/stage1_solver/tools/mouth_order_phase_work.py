"""Test specific order force and dissipation on the stored radius-one 0PN mouth.

The carrier is a measured input. A stationary passive layer requires the work
one-form (mu/n) dchi to be exact if all remaining stationary specific forces
are gradients. Test that requirement on two nearby regular chi-levels.
"""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline

from mouth_curved_phase import Domain


SIGMA = 1 / 3
NU = 0.4
ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_curved_phase.json")


def rectangle_work(domain, field, surface, delta):
    # An explicit smooth normal-width logistic profile near the measured graph.
    # At fixed chi its height is h + delta*g*log((1-chi)/chi), g=sqrt(1+h'^2).
    # Choose kappa/M by keeping sigma=kappa/(6*delta) and nu=M*kappa fixed.
    # The two constant-chi legs have exactly zero integral of (mu/n) dchi.
    nodes, weights = leggauss(64)
    chi = 0.5 + 0.05 * nodes
    log_odds = np.log((1 - chi) / chi)
    work = []
    flat_controls = []
    for rad in (0.25, 0.75):
        r = np.full_like(chi, rad)
        hp = surface(rad, 1)
        g = np.hypot(1, hp)
        gp = hp * surface(rad, 2) / g
        z = surface(rad) + delta * g * log_odds
        assert np.all(np.hypot(domain.rc-r, domain.rounding-z) > domain.rounding)
        f, fr, fw = domain.read(field, r, z)
        density = abs(f)**2
        jr, jw = (f.conj()*fr).imag, (f.conj()*fw).imag
        level_slope = hp + delta * gp * log_odds
        b = (6*SIGMA/NU) * chi*(1-chi) * (jw - level_slope*jr) / (g*density**2)
        work.append(float(0.05 * weights @ b))
        # Same orientation, uniform flat flow and constant density.
        flat_controls.append(float(0.05 * weights @ (chi*(1-chi) * (6*SIGMA/NU) * 0.1)))
    return dict(delta=delta, chi_bounds=[0.45, 0.55], radii=[0.25, 0.75],
                leg_work=work, signed_closed_work=work[0]-work[1],
                flat_control=flat_controls[0]-flat_controls[1])


def measure(folder, hashes):
    input_path = folder / "field.npz"
    with input_path.open("rb") as stream:
        assert file_digest(stream, "sha256").hexdigest() == hashes[str(input_path)]
    state = json.loads((folder/"stdout").read_text().splitlines()[0])
    domain = Domain(state["n"], state["far"])
    with np.load(input_path) as saved:
        np.testing.assert_array_equal(domain.points, saved["points"])
        radius, height, field = saved["r"], saved["h_0.4"], saved["psi"]
        q = saved["q_0.4"]
    surface = CubicSpline(radius, height)
    val, _, _ = domain.read(field, radius, height)
    density = abs(val)**2
    g = np.hypot(1, surface(radius, 1))
    un = q / (g*density)
    # For a thin, locally equilibrated logistic profile the midpoint coefficient
    # |grad chi|/M is 6*sigma/(4*nu). The exact edge test only needs regularity.
    mid_force = (6*SIGMA/(4*NU)) * un / density
    width_weight = 4*np.pi*radius**2
    heat = (SIGMA/NU) * simpson(un**2 * g * width_weight, x=radius)
    slope = 0.25/(4*np.pi*NU)
    axial_tension = 4*np.pi*SIGMA*slope/np.hypot(1,slope)
    rings = [0, .25, .5, .75, .99, 1]
    print(json.dumps(dict(mesh=state["n"], far=state["far"], radii=rings,
              q=[float(np.interp(r, radius, q)) for r in rings],
              midpoint_mu_per_n=[float(np.interp(r, radius, mid_force)) for r in rings],
              centered_density=density[0], sampled_min_normal_speed=un[:-1].min(),
              predicted_thin_layer_heat=heat, linearized_required_axial_load=SIGMA*state["rate"]/NU,
              finite_slope_tension_on_same_graph=axial_tension,
              loop_tests=[rectangle_work(domain, field, surface, delta) for delta in (.02,.04)])))
    assert abs(simpson(q*width_weight,x=radius) - state["rate"]) < 0.0001
    assert np.all(density>0) and heat>=0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="+", type=Path)
    args = parser.parse_args()
    fields = json.loads(ARCHIVE.read_text())["raw_fields_sha256"]
    for folder in args.folders:
        measure(folder, fields)
