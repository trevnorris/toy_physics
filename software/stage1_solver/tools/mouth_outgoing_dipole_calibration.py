"""Homogeneous acoustic dipole standards for a candidate outgoing reservoir.

Use exp(+i omega t), fluid velocity grad(phi), pressure -i rho omega phi,
and unit rho, sound speed and sphere radius. This is a full sphere at rest
in d dimensions. A mirror half-space halves the integrals for a tangential
dipole. The two drive normalizations have different low-frequency limits.
"""

import json

import mpmath as mp


def pair(z):
    return [float(mp.re(z)), float(mp.im(z))]


def calibrate(dim, omega):
    alpha = mp.mpf(dim - 2) / 2
    order = alpha + 1
    area = 2 * mp.pi**(mp.mpf(dim) / 2) / mp.gamma(mp.mpf(dim) / 2)
    angular = area / dim  # integral of n_x^2
    k = omega

    def radial(kind, radius):
        return radius**(-alpha) * kind(order, k * radius)

    def radiation_derivative(radius):
        x = k * radius
        # dH_n/dx = H_(n-1) - n H_n/x.
        return (x * mp.hankel2(order - 1, x) / mp.hankel2(order, x)
                - order - alpha) / radius

    def measured_flux(coefficient, radius):
        potential = lambda r: coefficient * radial(mp.hankel2, r)
        pressure = -1j * omega * potential(radius)
        velocity = mp.diff(potential, radius)
        return angular * radius**(dim - 1) * mp.re(pressure * mp.conj(velocity))

    a = mp.mpf(1)
    R1 = radial(mp.hankel2, a)
    derivative = mp.diff(lambda r: radial(mp.hankel2, r), a)
    rigid_coefficient = 1 / derivative
    specific = -1j * omega * R1 / derivative  # p/n_x for V=1
    force = angular * specific
    B = derivative  # a=1
    specific_wronskian = 2 * omega / (mp.pi * abs(B)**2)
    rigid_flux = [measured_flux(rigid_coefficient, r)
                  for r in (a, mp.mpf(18), 8 / k)]

    # A steady point source Q at X has static potential perturbation
    # -X Q n_x/(area*r^(d-1)). Set Q=1 and unit X, so |V|=omega.
    # Fix the singular dipole coefficient at the origin, never impose unit
    # normal velocity on a finite sphere for this second source.
    singular_coefficient = -1 / area
    translated_coefficient = (singular_coefficient * mp.pi * k**order
                              / (1j * 2**order * mp.gamma(order)))
    translated_flux = [measured_flux(translated_coefficient, r) / omega**2
                       for r in (a, mp.mpf(18), 8 / k)]
    translated_wronskian = 2 * angular * abs(translated_coefficient)**2 / (mp.pi * omega)

    # Independently impose the outgoing logarithmic derivative at each
    # matching radius on an annulus with real J,Y fundamental solutions.
    # Normalize velocity on the sphere, then recover its force.
    matching = []
    for radius in map(mp.mpf, ('9', '13.5', '18')):
        at_wall = [mp.diff(lambda r: radial(kind, r), a) for kind in (mp.besselj, mp.bessely)]
        dtn = radiation_derivative(radius)
        at_match = [mp.diff(lambda r: radial(kind, r), radius) - dtn * radial(kind, radius)
                    for kind in (mp.besselj, mp.bessely)]
        cj, cy = mp.lu_solve(mp.matrix([at_wall, at_match]), mp.matrix([1, 0]))
        potential = cj * radial(mp.besselj, a) + cy * radial(mp.bessely, a)
        recovered = angular * mp.re(-1j * omega * potential)
        matching.append(dict(radius=float(radius), radius_dtn=pair(radius * dtn),
                             resistance=float(recovered)))
        if not mp.almosteq(recovered, mp.re(force), rel_eps=mp.mpf('1e-35'), abs_eps=0):
            raise ArithmeticError(f'Failed annulus transfer at {dim=}, {omega=}, {radius=}')

    power = mp.re(force)
    results = dict(test='ideal_dipole', dimension=dim, omega=float(omega),
                   specific_impedance=pair(specific), integrated_impedance=pair(force),
                   specific_wronskian=float(specific_wronskian), rigid_power=map_float(rigid_flux),
                   translated_Q1_power_per_speed_squared=map_float(translated_flux),
                   translated_wronskian=float(translated_wronskian), match=matching)
    print(json.dumps(results, allow_nan=False))
    for actual, expected in [(specific_wronskian, mp.re(specific)),
                             *((f, power) for f in rigid_flux),
                             *((f, translated_wronskian) for f in translated_flux)]:
        if not mp.almosteq(actual, expected, rel_eps=mp.mpf('1e-35'), abs_eps=0):
            raise ArithmeticError(f'Failed independent power or Wronskian at {dim=}, {omega=}')
    if power <= 0 or translated_wronskian <= 0:
        raise ArithmeticError('Outgoing power must be positive')
    return power, translated_wronskian, specific


def map_float(values):
    return [float(v) for v in values]


def main():
    mp.mp.dps = 70
    print(json.dumps(dict(test='preparation', rho=1, sound_speed=1, sphere_radius=1,
                          phasor='exp(+i omega t)', power='2*mean outward power/|V|^2',
                          rigid_velocity=1, translated_static_flux=1, translated_displacement=1,
                          probe_radii=['1', '18', '8/k'], digits=mp.mp.dps)))
    frequencies = list(map(mp.mpf, ('.000225', '.000637', '.002', '.005', '.01', '.02')))
    for dim in [3, 4]:
        rows = [calibrate(dim, w) for w in frequencies]
        slopes = [[mp.log(rows[i+1][j] / rows[i][j]) / mp.log(frequencies[i+1] / frequencies[i])
                   for i in range(len(frequencies) - 1)] for j in (0, 1)]
        order = mp.mpf(dim) / 2
        leading = 2 * mp.pi / ((dim - 1)**2 * 2**(2 * order) * mp.gamma(order)**2)
        observed = mp.re(rows[0][2]) / frequencies[0]**(dim+1)
        measured_mass = mp.im(rows[0][2]) / frequencies[0]
        print(json.dumps(dict(test='low_frequency', dimension=dim, rigid_slopes=map_float(slopes[0]),
                              translated_slopes=map_float(slopes[1]), leading_specific=float(leading),
                              observed_specific=float(observed), reactive_mass_per_area=float(measured_mass))))
        if (max(abs(s-(dim+1)) for s in slopes[0]) > mp.mpf('.002')
                or max(abs(s-(dim-1)) for s in slopes[1]) > mp.mpf('1e-35')
                or abs(observed/leading-1) > mp.mpf('1e-6')
                or abs((dim-1)*measured_mass-1) > mp.mpf('1e-6')):
            raise ArithmeticError(f'Unexpected low-frequency dipole calibration in {dim}D')


if __name__ == '__main__':
    main()
