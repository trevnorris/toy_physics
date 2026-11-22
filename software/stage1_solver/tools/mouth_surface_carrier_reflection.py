"""Register a fixed neutral sheet against a normally incident carrier density wave.

The extra hypothesis assigns each carrier constituent the analytic plane mean
of the underlying Morse interaction with its far tail left untapered. Only the fixed-sheet
carrier problem is solved; two identical connected carrier reservoirs are assumed.
Length, energy and wavenumber below use the carrier's mouth units.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_bvp

from mouth_attachment_incidence import HalfOrderBand
from mouth_curved_site_front import sha
from mouth_segment_contact_hold import PITCH
from mouth_surface_site_transfer import sheet_mean
from neutral_patch_sites import emit


def stationary(mesh, potential, nb, tol):
    far = np.sqrt(nb)

    def ode(x, y):
        f = y[0]
        return np.array([y[1], 2 * (potential(x) + (f**8-nb**4)/4) * f])

    def jac(x, y):
        j = np.zeros((2, 2, len(x)))
        j[0, 1] = 1
        j[1, 0] = 2 * potential(x) + (9*y[0]**8-nb**4)/2
        return j

    def bc(left, right):
        return np.array([left[1], right[0]-far])

    guess = np.tanh(mesh+.002)
    result = solve_bvp(ode, bc, mesh, far*np.array([guess, 1-guess**2]),
                       fun_jac=jac, tol=tol, max_nodes=50000)
    if not result.success or np.min(result.y[0])<=0:
        raise ArithmeticError(f'Positive stationary state unresolved: {result.message}')
    return result


def parity_channel(base, nb, k, odd, tol):
    # Factoring the positive base out of each perturbation enforces its
    # exact global phase zero mode and avoids subtracting the large barrier.
    p = k*k/2
    omega = np.sqrt(p*(p+2*nb**4))
    decay = np.sqrt(k*k+4*nb**4)
    sound_ratio = p/omega            # (amplitude, phase) = (p/omega, 1)
    evan_ratio = -2*omega/decay**2   # decaying mode = (1, evan_ratio)

    def jac(x, _):
        f, g = base.sol(x)
        matrix = np.zeros((4, 4, len(x)))
        matrix[0, 2] = matrix[1, 3] = 1
        matrix[2, 0] = 4*f**8
        matrix[2, 1] = matrix[3, 0] = -2*omega
        matrix[2, 2] = matrix[3, 3] = -2*g/f
        return matrix

    def ode(x, y):
        return np.einsum('ijk,jk->ik', jac(x, y), y)

    def far_modes(y):
        a, b, da, db = y
        den = 1-evan_ratio*sound_ratio
        return ((b-evan_ratio*a)/den, (db-evan_ratio*da)/den,
                (a-sound_ratio*b)/den, (da-sound_ratio*db)/den)

    def bc(left, right):
        prop, dprop, evan, devan = far_modes(right)
        center = left[:2] if odd else left[2:]
        return np.r_[center, prop-1, devan+decay*evan]

    wave = solve_bvp(ode, bc, base.x, np.zeros((4,len(base.x))),
                     fun_jac=jac, tol=tol, max_nodes=50000)
    if not wave.success:
        raise ArithmeticError(f'Wave channel unresolved: {wave.message}')
    prop, slope, _, _ = far_modes(wave.y[:,-1])
    # A real standing wave is the sum of incoming exp(-ikz) and outgoing exp(ikz).
    return np.exp(-2j*k*wave.x[-1])*(1j*k*prop+slope)/(1j*k*prop-slope)


def run(args):
    previous = json.loads(args.record.read_bytes())
    assert previous['input_and_source_sha256'][str(args.section)] == sha(args.section)
    band = HalfOrderBand.load(args.section)
    trace = band.launch(0.,0.,0.)
    nb = float(trace['j_w']/trace['v_w'])
    mesh = np.unique(np.r_[np.linspace(0, .06, 121), np.linspace(.06,1,95),
                           np.linspace(1,args.far,1+int(16*(args.far-1)))])
    emit(test='setup', section=str(args.section), section_sha256=sha(args.section),
         record=str(args.record), record_sha256=sha(args.record), pitch=band.pitch,
         microscopic_length=PITCH, energy_scale=band.scale, scalar_order=band.scalar_order,
         background_density=nb, local_mouth_depth=band.depth, far=args.far,
         tolerance=args.tol, density_backing='identical carriers on both sides; fixed neutral sheet')
    for strength in args.strength:
        def potential(x):
            # The retained-mean continuum candidate keeps the analytic far tail.
            # The discrete part of the tag potential is tapered at 4.5–5 native lengths.
            return strength*band.scale*sheet_mean(x/PITCH, band.pitch)[0]
        base = stationary(mesh,potential,nb,args.tol)
        def excess(x):
            return nb/base.sol(x)[0]**2-1
        extra = 2*sum(quad(excess, a, b, epsabs=3e-10)[0]
                      for a,b in zip((0.,.01,.1,1.),(.01,.1,1.,args.far),strict=True))
        emit(test='base', strength=strength, center_density=base.y[0,0]**2,
             center_fraction=base.y[0,0]**2/nb, excess_phase_length=extra,
             far_field_gradient=base.y[1,-1], stationary_nodes=len(base.x))
        for k in args.k:
            even, odd = [parity_channel(base,nb,k,o,args.tol) for o in (False,True)]
            r, t = (even+odd)/2,(even-odd)/2
            reflect, transmit = abs(r)**2, abs(t)**2
            emit(test='scattering', strength=strength, k=k, omega=k*np.sqrt(nb**4+k*k/4),
                 r=[r.real,r.imag], t=[t.real,t.imag], reflected_energy=reflect,
                 transmitted_energy=transmit, identical_baths_flux_sum=reflect+transmit)
            assert abs(reflect+transmit-1)<1e-10
            if strength==0:
                assert abs(r)<2e-6 and abs(t-1)<2e-6


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section', type=Path, required=True)
    parser.add_argument('--record', type=Path, required=True)
    parser.add_argument('--far', type=float, default=18.)
    parser.add_argument('--tol', type=float, default=3e-7)
    parser.add_argument('--k', type=float, nargs='+', default=[.2,.1,.05,.025,.0125,.00625,.003125,.001])
    parser.add_argument('--strength', type=float, nargs='+', default=[0.,1.])
    args = parser.parse_args()
    if args.far<10 or args.tol<=0 or min(args.k)<=0 or min(args.strength)<0:
        parser.error('Use positive wavenumbers and tolerance, nonnegative strength, far >= 10')
    run(args)
