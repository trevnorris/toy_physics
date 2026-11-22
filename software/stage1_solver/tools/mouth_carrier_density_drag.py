"""Cold linear density-wave drag on a weak prescribed probe in a uniform 4D carrier.

Only the local density and tangential speed come from the archived mouth.
The outgoing-wave calculation assumes infinite homogeneous material in its rest
frame. The Fourier transform of the test potential is g exp(-k^2/2).
"""
import argparse

import numpy as np
from scipy.integrate import quad

from mouth_middle_tag_grip import cubic_trace
from neutral_patch_sites import emit


def radiation_force(density, relative_speed):
    """Two quadratures of force on the probe, divided by g^2."""
    speed = abs(relative_speed)
    sound = density**2
    if speed <= sound:
        return 0., 0.
    # Retarded density susceptibility has residue n k^2 / (2 Omega).
    # Integrate the two frequency poles using first angle, then wave number.
    end = 2*np.sqrt(speed**2-sound**2)
    factor = density * np.pi * (4*np.pi) / (2*np.pi)**4

    def radial(k):
        cone = np.sqrt(max(0., (end-k)*(end+k))) / (2*speed)
        return k**4*np.exp(-k*k)*cone

    def angular(mu):
        k = 2*np.sqrt(max(0., (speed*mu-sound)*(speed*mu+sound)))
        return 4*mu*k**3*np.exp(-k*k)*np.sqrt(max(0., 1-mu**2))

    fk = factor/speed**2*quad(radial, 0., end, epsabs=1.e-14, epsrel=1.e-10)[0]
    fa = factor*quad(angular, sound/speed, 1., epsabs=1.e-14, epsrel=1.e-10)[0]
    np.testing.assert_allclose(fk, fa, rtol=1.e-9, atol=1.e-15)
    if min(fk, fa) <= 0:
        raise ArithmeticError('Above-threshold control failed to radiate')
    return -np.sign(relative_speed)*fk, -np.sign(relative_speed)*fa


def check_spectrum(density):
    # Differentiate the native enthalpy h(n)=(n^4-1)/4 and diagonalize the
    # complex field linearization instead of using the hydrodynamic dispersion.
    delta = 1.e-5*density
    derivative = ((density+delta)**4-(density-delta)**4)/(8*delta)
    compressibility_term = density*derivative
    errors = []
    for k in np.geomspace(.02, 5., 21):
        diag = k*k/2 + compressibility_term
        matrix = [[diag, compressibility_term], [-compressibility_term, -diag]]
        eigen = np.sort(np.linalg.eigvals(matrix))
        predicted = k*np.sqrt(density**4+k*k/4)*np.array([-1., 1.])
        np.testing.assert_allclose(eigen, predicted, atol=1.e-11, rtol=2.e-10)
        errors.append(max(abs(eigen-predicted)))
    return max(errors)


def run(mesh):
    trace, intake = cubic_trace(mesh)
    np.testing.assert_allclose(intake, .25, atol=1.e-8, rtol=0.)
    cuts = [1.175, trace.crossing(.5)]
    # Include both one-sided values of the elementwise phase gradient.
    interior = trace.edges[(cuts[0] < trace.edges) & (trace.edges < cuts[1])]
    knots = np.r_[cuts[0], interior, cuts[1]]
    path = np.unique(np.r_[np.concatenate([np.linspace(a, b, 65) for a, b
                                          in zip(knots[:-1], knots[1:])]),
                           np.nextafter(interior, -np.inf)])
    density = abs(trace.sample(trace.psi, path)[0])**2
    velocity = trace.speed(path)*np.hypot(1., trace.geometry.radius(path)[1])
    emit(test='path_sample', mesh=mesh, intake=intake, depths=cuts, points=len(path),
         density_range=[min(density), max(density)], max_sampled_mach=max(velocity/density**2))
    for label, i in [('upstream', 0), ('half_order', -1)]:
        n, u = density[i], velocity[i]
        c = n**2
        emit(test='local_state', mesh=mesh, cut=label, depth=path[i],
             order=trace.sample(trace.chi, path[i])[0], density=n, tangent_speed=u,
             sound=c, mach=u/c, spectrum_max_error=check_spectrum(n))
        for magnitude in [0., u/2, u, c, 1.25*c]:
            for rel in sorted({magnitude, -magnitude}):
                force_k, force_mu = radiation_force(n, rel)
                emit(test='steady_probe', mesh=mesh, cut=label, probe_minus_carrier=rel,
                     mach=rel/c, force_per_g2_radial=force_k,
                     force_per_g2_angular=force_mu, radiated_power_per_g2=-rel*force_k)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26, 34), required=True)
    run(parser.parse_args().mesh)
