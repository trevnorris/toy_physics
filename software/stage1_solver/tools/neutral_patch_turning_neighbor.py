"""Torque on a turning two-layer patch from an overdamped neighboring annulus.

Angles are in radians. Positive torque is the torque the actuator supplies
to the exterior. The annulus exerts its negative on the patch. Material
units, the isolated same-plane Morse law and unit site drag are those of the
fixed-order patch. The homogeneous, bonded exterior and its modulus law
are hypotheses. A planar stress-free nearest-neighbor shear coupon fixes
the coefficient at the entry gate without using patch turns or forces.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh_tridiagonal, solve
from scipy.special import j1, kve, y1

from neutral_patch_resting_order_hold import LEVELS
from neutral_patch_sites import COUNT, WITHIN, morse


LAYERS = 2
ENTRY = .95
AREA_PER_SITE = np.sqrt(3.) * WITHIN[2]**2 / 2
RADIUS = np.sqrt(COUNT * AREA_PER_SITE / np.pi)
DRAG_PER_AREA = LAYERS / AREA_PER_SITE
RATIOS = (2, 4, 8, 16)
LAWS = ('curl', 'symmetric')
# Fixed-order stationary patch trial, Q=0.24999999999821823.
TRANSIT = 5.296673641978323


def sheet_modulus(order):
    """Entry shear curvature per sheet area for both planes, times chi/entry."""
    depth, inverse_range, pitch = WITHIN
    bond_curvature = 2 * depth * inverse_range**2
    # Of the six neighbors per site, count three distinct bonds. The sum
    # of (n_x n_y)^2 over their three directions is 3/8.
    return LAYERS * bond_curvature * pitch**2 * (3 / 8) / AREA_PER_SITE * (order / ENTRY)


def shear_coupon():
    """Affine simple shear of a periodic triangular Bravais primitive cell."""
    directions = np.arange(3) * (np.pi / 3)
    bonds = WITHIN[2] * np.column_stack((np.cos(directions), np.sin(directions)))

    def density(shear):
        shifted = bonds.copy()
        shifted[:, 0] += shear * shifted[:, 1]
        return LAYERS * morse(np.linalg.norm(shifted, axis=1), *WITHIN)[0].sum() / AREA_PER_SITE

    exact = sheet_modulus(ENTRY)
    readings = []
    for h in (2e-3, 1e-3, 5e-4, 2.5e-4):
        curvature = (density(h) - 2*density(0.) + density(-h)) / h**2
        readings.append(dict(step=h, curvature=curvature, relative_error=abs(curvature/exact - 1)))
    if readings[-1]['relative_error'] > 3e-7:
        raise ArithmeticError('Finite shear-energy curvature disagrees with the lattice sum')
    return dict(morse_depth_range_rest=list(WITHIN), layers=LAYERS, sites_per_layer=COUNT,
                area_per_site=AREA_PER_SITE, inner_radius=RADIUS, drag_per_area=DRAG_PER_AREA,
                entry_order=ENTRY, two_plane_entry_modulus=exact, curvature_readings=readings)


def annulus_matrices(order, ratio, law='curl', elements=128):
    """Return radial coordinates, K and D. The fixed outer angle is eliminated.

    Row zero is the patch angle. Rows 1...elements-1 are free material
    angles. The physical tangential displacement is r*theta. D contains
    only the exterior drag, including the inner node's share of the first
    exterior element. Couple by pulling back K,D with the same map from
    site displacements to the inner angle. dE/dtheta is a torque.
    """
    if law not in LAWS or not 0 < order <= 1 or ratio <= 1 or elements < 2:
        raise ValueError('Invalid annulus law, order, radius ratio or mesh size')
    radii = np.geomspace(RADIUS, RADIUS*ratio, elements+1)
    inner2, outer2 = radii[:-1]**2, radii[1:]**2
    width2 = outer2 - inner2
    scale = 4 * np.pi * sheet_modulus(order) / width2
    # On each cell the angle is affine in 1/r². u_phi = A*r+B/r is the
    # exact static radial interpolant. The elastic energy is the integral
    # of (mu/2)*(r theta' [+ 2 theta])² against 2*pi*r dr.
    off = -scale * inner2 * outer2
    if law == 'curl':
        left = scale * inner2**2
        right = scale * outer2**2
    else:
        left = right = -off
    diagonal = left.copy()
    diagonal[1:] += right[:-1]
    stiffness = np.diag(diagonal) + np.diag(off[:-1], 1) + np.diag(off[:-1], -1)
    # Full element row sums of 2*pi*zeta*integral r³ N_i N_j dr. Compute
    # these before fixing the outer ring to keep the physical area weight.
    lump = np.pi * DRAG_PER_AREA * width2 / 2
    drag = lump * inner2
    drag[1:] += lump[:-1] * outer2[:-1]
    return radii[:-1], stiffness, np.diag(drag)


def static_torque(order, ratio, law):
    scale = 4 * np.pi * sheet_modulus(order) * RADIUS**2
    return scale * (1 if law == 'curl' else ratio**2) / (ratio**2-1)


def field_energy(radii, theta, modulus, law):
    """Independent quadrature of strain energy for virtual-work testing."""
    left, right = radii[:-1], radii[1:]
    z, weights = np.polynomial.legendre.leggauss(12)
    r = (left+right)[:, None]/2 + (right-left)[:, None]*z/2
    fluctuation = (theta[:-1]-theta[1:]) / (left**-2-right**-2)
    constant = theta[:-1] - fluctuation/left**2
    strain = 2*constant[:, None] if law == 'curl' else -2*fluctuation[:, None]/r**2
    return np.sum(np.pi * modulus * (right-left)[:, None]/2 * weights*r*strain**2)


def measure(order, ratio, law, elements, transit):
    nodes, k, d = annulus_matrices(order, ratio, law, elements)
    static = np.r_[1., solve(k[1:, 1:], -k[1:, 0], assume_a='pos')]
    exact = static_torque(order, ratio, law)
    measured = k[0] @ static
    b = RADIUS*ratio
    profile = (nodes**-2 - b**-2)/(RADIUS**-2 - b**-2)
    radii = np.r_[nodes, b]
    perturb = np.sin(np.linspace(.4, 3.1, elements))
    # Off equilibrium, probe virtual work along a direction that also
    # changes the angle at the interface. The far ring remains fixed.
    trial = static + .13 * np.sin(3.5*np.log(nodes/RADIUS))
    h = 1e-4
    high, low = [field_energy(radii, np.r_[trial + sign*h*perturb, 0.],
                              sheet_modulus(order), law) for sign in (1., -1.)]
    work = perturb @ k @ trial
    err_work = abs((high-low)/(2*h) - work) / max(abs(work), exact)
    err_static = abs(measured/exact - 1)
    err_shape = max(abs(static-profile))
    # The interior relaxation is symmetric in drag-weighted coordinates.
    # A held unit angle from rest has reaction K_static + sum w exp(-rate*t)
    # at t>0. The inner node's lumped drag also gives an impulse at t=0.
    root_drag = np.sqrt(np.diag(d)[1:])
    rates, modes = eigh_tridiagonal(np.diag(k)[1:]/root_drag**2,
                                   np.diag(k, 1)[1:] / (root_drag[:-1]*root_drag[1:]))
    drive = k[1:, 0] / root_drag
    weights = (modes.T @ drive)**2 / rates
    step = measured + weights @ np.exp(-rates*transit)
    # A second driver turns through the same angle in T/4 with a half
    # cosine and holds thereafter. At T its angle is one and velocity zero.
    duration = transit/4
    omega = np.pi/duration
    ramp_memory = .5*omega**2 * (1+np.exp(-rates*duration))/(rates**2+omega**2)
    smooth = measured + weights @ (ramp_memory * np.exp(-rates*(transit-duration)))
    if max(err_shape, err_static, err_work) > 3e-8 or min(rates) <= 0:
        raise ArithmeticError('Annulus did not reproduce static shape, work, or positive rates')
    return dict(order=order, outer_over_inner_radius=ratio, law=law, elements=elements,
                static_exact=exact, static_discrete=measured, step_reaction_at_T=step,
                smooth_reaction_at_T=smooth, static_relative_error=err_static,
                static_field_max_abs_error=err_shape, virtual_work_relative_error=err_work)


def run():
    calibration = shear_coupon()
    rows = []
    for order in LEVELS:
        for ratio in RATIOS:
            history = []
            offset = static_torque(order, ratio, 'symmetric') - static_torque(order, ratio, 'curl')
            for elements in (64, 128, 256):
                curl, symmetric = [measure(order, ratio, law, elements, TRANSIT) for law in LAWS]
                for key in ('static_discrete', 'step_reaction_at_T', 'smooth_reaction_at_T'):
                    if abs((symmetric[key]-curl[key])/offset - 1) > 1e-8:
                        raise ArithmeticError('Two exterior laws violate their boundary identity')
                history.append((curl, symmetric))
            for coarse, medium, fine in zip(*history, strict=True):
                fine['mesh_readings'] = [
                    {key: row[key] for key in ('elements', 'step_reaction_at_T', 'smooth_reaction_at_T')}
                    for row in (coarse, medium, fine)]
                for key in ('step_reaction_at_T', 'smooth_reaction_at_T'):
                    if abs(fine[key]/medium[key]-1) > .002:
                        raise ArithmeticError('Unresolved exterior torque at transit time')
                rows.append(fine)
    return dict(transit=TRANSIT, intake=.24999999999821823, calibration=calibration,
                torque_convention='Actuator torque delivered to the annulus per radian of turning',
                smooth_turn='(1-cos(pi*t/(T/4)))/2 to T/4, then held through T',
                measurements=rows)


def unbounded_exterior(finite_path):
    """Continuum step torque and a Laplace-transform check, with no outer ring.

    Set p=s*zeta*a²/mu. The decaying solution has u_phi proportional to
    K1(sqrt(p)*r/a). Its curl traction gives 2*pi*mu*a²*sqrt(p)*K0/K1
    per unit boundary angle in Laplace space. Dividing by p for a step,
    the jump of K0/(sqrt(p)*K1) on the negative real axis gives the
    relaxation-rate density below (rate x² in sheet diffusion units).
    """
    finite_rows = json.loads(finite_path.read_text())['measurements']

    def density(x):
        return (4 / np.pi**2) / (x * (j1(x)**2 + y1(x)**2))

    transform_checks = []
    for p in (.1, 1., 10.):
        z = np.sqrt(p)
        exact = kve(0, z)/(z*kve(1, z))
        quadrature = quad(lambda x: density(x)/(p+x*x), 0, np.inf,
                          epsabs=1e-13, epsrel=1e-12)[0]
        error = abs(quadrature/exact-1)
        if error > 2e-10:
            raise ArithmeticError('Continuum spectral density failed its Bessel transform')
        transform_checks.append(dict(p=p, exact=exact, quadrature=quadrature, relative_error=error))

    readings = []
    for order in LEVELS:
        mu = sheet_modulus(order)
        diffusion_time = mu*TRANSIT/(DRAG_PER_AREA*RADIUS**2)
        dimensional = 2*np.pi*mu*RADIUS**2
        # Scale y=x*sqrt(diffusion_time) to resolve the Gaussian for every gate.
        factor = np.sqrt(diffusion_time)
        integral = quad(lambda y: np.exp(-y*y)*density(y/factor)/factor, 0, np.inf,
                        epsabs=1e-13, epsrel=1e-12)[0]
        curl = dimensional*integral
        for law in LAWS:
            step = curl if law == 'curl' else curl + 2*dimensional
            measured = [r for r in finite_rows if r['order'] == order and r['law'] == law]
            if len(measured) != len(RATIOS):
                raise ArithmeticError('Finite-radius comparison missing')
            comparisons = []
            for row in measured:
                coarse, fine = row['mesh_readings'][-2:]
                # For a positive held turn the heat-equation comparison
                # principle bounds the reaction below by the unbounded
                # reaction. Allow measured spatial discretization error.
                mesh_change = abs(fine['step_reaction_at_T'] - coarse['step_reaction_at_T'])
                gap = row['step_reaction_at_T'] - step
                if gap < -mesh_change-1e-9*step:
                    raise ArithmeticError('Finite-ring torque below continuum lower bound')
                comparisons.append(dict(outer_over_inner_radius=row['outer_over_inner_radius'],
                                        finite_minus_unbounded=gap, last_mesh_change=mesh_change))
            readings.append(dict(order=order, law=law, sheet_diffusion_time=diffusion_time,
                                 unbounded_step_reaction_at_T=step, comparisons=comparisons))
    return dict(transit=TRANSIT, transform_checks=transform_checks, measurements=readings)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--unbounded-check', type=Path, metavar='FINITE_JSON',
                        help='compute the infinite-sheet comparison with an existing finite-ring run')
    args = parser.parse_args()
    result = unbounded_exterior(args.unbounded_check) if args.unbounded_check else run()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
