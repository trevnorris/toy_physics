"""Contact-free 0PN controls for a registry holding-force measurement."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_contact_tracking import WallTrace
from mouth_curved_phase_stability import weighted_gram
from mouth_free_order_diffuse import NU, read
from mouth_native_order_contrast import rate_field
from mouth_order_wave_pulse import ROOT, emit


def shear_profile(w_edges, probe):
    # Compatible displacements u_x=g*b(w), u_y=h*b(w). With piecewise linear
    # b' on native element edges, its square belongs to the Q2 order space.
    i = np.searchsorted(w_edges, probe)
    knots = w_edges[i+np.array([-8, -3, 3, 8, 13, 18])]
    pos_area = (knots[1]-knots[0])/2 + knots[2]-knots[1] + (knots[3]-knots[2])/2
    height = -2*pos_area/(knots[-1]-knots[3])
    values = [0., 1., 1., 0., height, 0.]
    if not knots[1] < probe < knots[2]:
        raise ValueError('Continuum probe outside the uniform shear plateau')
    # The signed areas cancel, so b returns to zero outside the drive region.
    return lambda w: np.interp(w, knots, values, left=0., right=0.)


def continuum(mesh):
    wall = WallTrace(mesh, .5)
    folder = ROOT/f'mouth_flared_delayed_n{mesh}'
    order_folder = ROOT/f'contact_supply_n{mesh}_01'
    dom, psi = prepare(SimpleNamespace(source=folder))
    record = [json.loads(s) for s in (order_folder/'stdout').read_text().splitlines()]
    calibration = next(s for s in record if s['test']=='calibration')
    if wall.intake != calibration['intake']:
        raise ValueError('Calibration and traced carrier differ')
    with np.load(folder/'field.npz') as carrier:
        carrier_keys = list(carrier.files)
    with np.load(order_folder/'order_0.50000.npz') as saved:
        chi = saved['chi']
        order_keys = list(saved.files)
    emit(test='state', mesh=mesh, intake=wall.intake, head=wall.head, fields=wall.hashes,
         carrier_arrays=carrier_keys, order_arrays=order_keys, mobility=calibration['mobility'])

    start = wall.crossing(.95)
    entry_radius = wall.geometry.radius(start)[0]
    for level in (.75, .50, .30):
        w = wall.crossing(level)
        radius = wall.geometry.radius(w)[0]
        emit(test='carrier_markers', mesh=mesh, wall_order=level, w=w,
             elapsed=wall.transit(start, w, 8), arc_spacing_over_entry=radius/entry_radius)

    # Interior volume stress, measured with no contact model or wall spring.
    radius_at = lambda w: .8*dom.radius(w)[0]
    wc = brentq(lambda w: read(dom, chi, radius_at(w), w)[0]-.5, wall.crossing(.95)-.12,
                wall.crossing(.01)+.12)
    rc = radius_at(wc)
    local_order = float(read(dom, chi, rc, wc)[0])
    profile = shear_profile(wall.edges, wc)
    nq = dom.N @ np.ones_like(chi)
    curl = profile((dom.N @ dom.points[:, 1])/nq)
    weight = dom.weight*np.exp(-np.angle(dom.N @ psi)/NU)
    free_basis = dom.N[:, dom.free].multiply((1/nq)[:, None]).tocsr()
    mass = splu(weighted_gram(free_basis, weight).tocsc())
    mu = 1.
    mobility = calibration['mobility']

    def probe(g, h):
        # Subtract the common no-strain residual BEFORE assembly. All copies
        # start with the same carrier, order gradients, free energy and supply.
        source = free_basis.T @ (weight*mobility*mu*(g*g+h*h)*curl**2/2)
        initial_order_rate = rate_field(dom, mass, source)
        order_rate = float(read(dom, initial_order_rate, rc, wc)[0])
        return mu*local_order*g, mu*g*order_rate

    for g in (.001, -.001, 0.):
        reference, base_rate = probe(g, 0.)
        for h in (1e-5, 5e-6):
            positive, plus_rate = probe(g, h)
            negative, minus_rate = probe(g, -h)
            even_rate = (plus_rate+minus_rate)/2-base_rate
            predicted = -mobility*mu**2*g*h**2/2
            emit(test='contact_free_bulk', mesh=mesh, position=[rc, wc], local_order=local_order,
                 preload_strain=g, transverse_strain=h, frozen_stress=[negative, reference, positive],
                 order_induced_stress_rates_relative_to_unstrained=[minus_rate, base_rate, plus_rate],
                 even_frozen=(positive+negative)/2-reference, odd_frozen=(positive-negative)/2,
                 even_initial_rate=even_rate, odd_initial_rate=(plus_rate-minus_rate)/2,
                 local_continuum_even_rate=predicted)
            np.testing.assert_allclose(even_rate, predicted, atol=1e-22, rtol=2e-8)
            if positive != negative or positive != reference or plus_rate != minus_rate:
                raise ArithmeticError('Broken control symmetry')


def identifiability():
    import sympy as s
    q, y, a, h, lam, u = s.symbols('q y a h lambda u', real=True, nonzero=True)
    # Smooth interior perturbation in an occupied patch; it is not used as a
    # material law. The envelope and its first derivative vanish at the edges.
    phi = u**2*(1-u**2)**2*(u**2-s.Rational(1, 3))
    energy_difference = lam*y**2*phi.subs(u, q/a)/2
    holding = s.diff(energy_difference, q)
    gradients = [s.diff(energy_difference, v) for v in (q, y)]
    axes = [s.simplify(f.subs(z, 0)) for f in gradients for z in (q, y)]
    even = s.simplify((holding.subs(y,h)+holding.subs(y,-h))/2-holding.subs(y,0))
    odd = s.simplify((holding.subs(y,h)-holding.subs(y,-h))/2)
    mean = s.integrate(phi, (u, -1, 1))/2
    half = s.simplify(even.subs({q:s.Rational(1,8), a:s.Rational(1,4)}))
    emit(test='undetermined_contact_readout', same_uniaxial_force_differences=list(map(str, axes)),
         change_in_uniform_mean_stiffness=str(lam*mean), change_in_odd_read=str(odd),
         change_in_even_read_at_half_order_eighth_pitch=str(half))
    if axes != [0]*4 or mean != 0 or odd != 0 or s.simplify(half-lam*h*h/4) != 0:
        raise ArithmeticError('Mixed response witness failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--mesh', choices=(34,48), type=int)
    choice.add_argument('--identifiability', action='store_true')
    args = parser.parse_args()
    if args.identifiability:
        identifiability()
    else:
        continuum(args.mesh)
