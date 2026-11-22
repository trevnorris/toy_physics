"""Tangential force from a site's order disturbance in a mean-held local 4D sample.

The same isotropic periodized Gaussian samples order for both surface laws and
depths. Each force is a response coefficient per (site mass * specific slope)^2.
The unconstrained converting front and the total mouth force are outside this test.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_free_order_diffuse import NU, SIGMA
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import emit


DATA = Path('research/pde_ledger_v3/_measurements')


def archived_rows(path, name):
    run = json.loads(path.read_text())['guarded_runs'][name]
    assert run['guard_exit'] == run['child_exit'] == 0 and run['limits_verified']
    raw = run['literal_stdout']
    assert hashlib.sha256(raw.encode()).hexdigest() == run['stdout_sha256']
    return [json.loads(row) for row in raw.splitlines()]


def input_states(mesh):
    energy_file = DATA/'mouth_bound_site_transfer.json'
    case = 'primary' if mesh == 34 else 'mesh26'
    energy = archived_rows(energy_file, case)
    sample, = [r for r in energy if r['test'] == 'sample']
    transfers = [r for r in energy if r['test'] == 'transfer']
    flow = archived_rows(DATA/'mouth_neutral_face_and_slip.json', f'carrier{mesh}')
    states = [r for r in flow if r['test'] == 'local_state']
    assert len(transfers) == len(states) == 2 and 0 < sample['scale']
    for state, key in zip(states, ('entry', 'destination'), strict=True):
        for row in transfers:
            np.testing.assert_allclose(row[key]['chi'], state['order'], rtol=0, atol=1e-10)
        state['slopes'] = {r['law']: r[key]['order_derivative']*sample['scale'] for r in transfers}
    path = Path(f'_scratch/s11c/mouth_free_order_cubic_n{mesh}_d04/stdout')
    domain = json.loads(path.read_text().splitlines()[0])
    assert domain['profile'] == 'cubic' and domain['n'] == mesh
    np.testing.assert_allclose([domain['source_flux'],domain['nu'],domain['sigma']],
                               [.25, NU, SIGMA], atol=1.e-8, rtol=0.)
    np.testing.assert_allclose(domain['mobility'], NU/(6*SIGMA*domain['delta']))
    return states, domain


def spectrum(length, modes):
    # Retain every signed wave vector in [-modes,modes]^4 with its multiplicity.
    # The last three components enter only through their squared lengths.
    modes1 = np.arange(-modes, modes+1)
    transverse = sum(np.ix_(*([modes1**2]*3)))
    shell, count = np.unique(transverse, return_counts=True)
    axial, shell = np.broadcast_arrays(modes1[:,None], shell)
    nonzero = axial**2+shell > 0
    factor = 2*np.pi/length
    return (factor*axial[nonzero], factor**2*(axial**2+shell)[nonzero],
            np.broadcast_to(count, axial.shape)[nonzero])


def stable_mode_bound(speed, b):
    # Supremum of lambda*k_x^2/(lambda^2+k_x^2*speed^2) over all stable modes.
    # b=(1-6 chi (1-chi))/delta^2; k_x^2 <= k^2. For b<0 use y=1+b/k^2.
    if b >= 0:
        return 1/NU
    q = np.sqrt(-b)
    if speed == 0:
        return None
    if speed >= NU*q:
        return 1/NU
    return q/(speed*(2-speed/(NU*q)))


def responses(state, domain, args):
    u, order = state['tangent_speed'], state['order']
    length, radius = args.length, args.width
    kx, k2, copies = spectrum(length, args.modes)
    volume, mobility = length**4, domain['mobility']
    b = (1-6*order*(1-order))/domain['delta']**2
    rate = NU*(k2+b)
    if min(rate) <= 0:
        raise ValueError('Uniform background unstable at a retained nonzero wave number')
    shape = np.exp(-.5*radius**2*k2)
    source = shape/volume
    static = -mobility*source/rate
    weighted = abs(copies @ (shape*static))
    gamma0 = (copies @ (kx**2*shape**2/rate**2))*mobility/volume
    map_gamma0 = volume/mobility*(copies @ (kx**2*static**2))
    np.testing.assert_allclose(gamma0, map_gamma0, rtol=2e-13)
    wx = (2*np.pi/length)**2
    box_bound = max(1., wx/(wx+b))/NU if b < 0 else 1/NU
    emit(test='static', cut=state['cut'], chi=order, speed=u, b=b, minimum_nonzero_rate=min(rate),
         site_specific_order_slopes=state['slopes'], static_overlap_per_abs_source=weighted,
         zero_slip_coefficient_per_source_squared=gamma0,
         static_map_coefficient_per_source_squared=map_gamma0, zero_slip_ratio=gamma0/weighted,
         all_slip_box_shape_bound=box_bound, iso_gaussian_dc_box_bound=box_bound/4,
         unheld_uniform_decay_rate=NU*b)
    # D = velocity of transporting material minus velocity of the site.
    # Zero, counter-drives, half drives, and a common boost test the two protocols.
    for carrier, host, site in [(0.,0.,0.), (u,u,u), (2*u,u,u), (u,2*u,u),
                                *((a*u,0.,0.) for a in (1.,-1.,.5,-.5)),
                                *((0.,a*u,0.) for a in (1.,-1.,.5,-.5))]:
        for frame, d in [('carrier', carrier-site), ('sheet', host-site)]:
            mode = -mobility*source/(rate+1j*kx*d)
            # mu is projected off the held uniform mode. Intrinsic stress integrates
            # to zero on the periodic sample. The transporting component feels
            # integral mu grad(eta). Site force comes from translating its footprint.
            mu = (rate/mobility)*mode + source
            grad = 1j*kx*mode
            site_force = -(copies @ (shape*grad.real))
            medium_force = volume*(copies @ (mu.conj()*grad).real)
            heat = mobility*volume*(copies @ abs(mu)**2)
            np.testing.assert_allclose(medium_force, -site_force, atol=2e-10, rtol=2e-9)
            np.testing.assert_allclose(heat, d*site_force, atol=1e-14, rtol=2e-9)
            coefficient = gamma0 if d == 0 else site_force/d
            ratio = coefficient/weighted
            all_boxes = stable_mode_bound(abs(d), b)
            if not (0 < coefficient <= gamma0*(1+1e-12) and ratio <= box_bound/4*(1+1e-12)):
                raise ArithmeticError('Invalid response relative to static overlap')
            ceilings = {}
            for law, slope in state['slopes'].items():
                # Each archived slope is negative, so a positive footprint samples a
                # positive static perturbation. The site-weighted order cannot exceed 1.
                if slope >= 0:
                    raise ArithmeticError('Order headroom must be evaluated for the actual slope sign')
                headroom = 1-order
                ceilings[law] = dict(rate_per_unit_static_overlap=abs(slope)*ratio,
                    max_rate_in_gaussian_trial=abs(slope)*ratio*headroom,
                    max_rate_any_stable_footprint=(None if all_boxes is None else abs(slope)*headroom*all_boxes))
            emit(test='drive', cut=state['cut'], advected_by=frame, carrier=carrier, sheet=host,
                 site=site, transport_minus_site=d, site_force_per_source_squared=site_force,
                 medium_force_per_source_squared=medium_force, heat_per_source_squared=heat,
                 ratio_to_static_overlap=ratio, all_stable_mode_shape_bound=all_boxes, rates=ceilings)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26,34), default=34)
    parser.add_argument('--length', type=float, default=.12)
    parser.add_argument('--width', type=float, default=PITCH)
    parser.add_argument('--modes', type=int, default=14)
    args=parser.parse_args()
    if min(args.length,args.width,args.modes)<=0:
        parser.error('Use a positive length, width and cutoff')
    states,domain=input_states(args.mesh)
    emit(test='setup', mesh=args.mesh, length=args.length, width=args.width,
         max_mode_per_axis=args.modes, delta=domain['delta'], nu=NU, mobility=domain['mobility'],
         normalization='site force and dissipation per (site mass times specific order derivative) squared')
    for state in states:
        responses(state, domain, args)


if __name__=='__main__':
    main()
