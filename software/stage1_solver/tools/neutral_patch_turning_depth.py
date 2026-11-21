"""Static torque on a prescribed rotating inclusion in a homogeneous shell.

Lengths are in inner radii and the shear modulus is unity. The planar annulus
is per unit invariant transverse length, the spherical shell is in three
dimensions. Restricting the four-dimensional brane to w-independent fields
also leaves a factor of participating w depth. All torques are actuator torques.
"""
import json

import sympy as sp


def require_zero(expression):
    if sp.simplify(expression) != 0:
        raise ArithmeticError(str(expression))


def test_shell(dim):
    r, b = sp.symbols('r b', positive=True)
    t, to, z = sp.symbols('theta theta_outer cos_polar')
    f = sp.Function('f')(r)
    if dim == 2:
        normal = sp.Matrix([1, 0])
        generator = sp.Matrix([[0, -1], [1, 0]])
    else:
        normal = sp.Matrix([sp.sqrt(1-z**2), 0, z])
        generator = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    # At one azimuth, grad u comes directly from u = f(r) Jx.
    grad = f*generator + f.diff(r)*r*generator*normal*normal.T
    require_zero(grad.trace())
    curl = grad-grad.T
    strain = (grad+grad.T)/2
    densities = {'curl': (curl.T*curl).trace()/4, 'strain': (strain.T*strain).trace()}
    results = {}
    radial_profile = to + (t-to)*(r**(-dim)-b**(-dim))/(1-b**(-dim))
    for law, density in densities.items():
        angular = sp.expand(density)
        if dim == 3:
            angular = sp.integrate(angular, (z, -1, 1))
        lagrangian = 2*sp.pi*r**(dim-1)*angular
        radial_momentum = lagrangian.diff(f.diff(r))
        euler = sp.simplify(radial_momentum.diff(r)-lagrangian.diff(f))
        require_zero(euler.subs(f, radial_profile).doit())
        require_zero(radial_profile.subs(r, 1)-t)
        require_zero(radial_profile.subs(r, b)-to)
        # Check the unrestricted bulk Euler equation as well. In either
        # material it reduces to the vector Laplacian since div u vanishes.
        require_zero(radial_profile.diff(r, 2)+(dim+1)*radial_profile.diff(r)/r)
        momentum = radial_momentum.subs(f, radial_profile).doit()
        inner = -momentum.subs(r, 1)
        outer = momentum.subs(r, b)
        radial_energy = sp.expand(lagrangian.subs(f, radial_profile).doit())
        energy = sp.factor(sp.integrate(radial_energy, (r, 1, b)))
        require_zero(energy.diff(t)-inner)
        require_zero(energy.diff(to)-outer)
        core_density = lagrangian.subs(f, t).doit()
        core_energy = sp.integrate(core_density, (r, 0, 1))
        results[law] = (energy, core_energy, sp.factor(inner.diff(t)))

    curl_energy, curl_core, kc = results['curl']
    strain_energy, strain_core, ks = results['strain']
    # Equal moduli across a filled core conserve the fixed-rim energy difference.
    require_zero((curl_energy+curl_core-strain_energy-strain_core).subs(to, 0))
    return {
        'dimensions': dim,
        'stiffness_normalization': 'mu*a^dimension',
        'b_symbol': 'outer radius / inner radius',
        'inner_curl_K': str(kc),
        'inner_strain_K': str(ks),
        'filled_core_extra_curl_K': str(curl_core.diff(t, 2)),
        'outer_restraining_torque_per_inner_angle': {
            law: str(sp.factor(energy.diff(to, t))) for law, (energy, _, _) in results.items()
        },
        'samples_including_unbounded_fixed_far_frame': [
            {'b_over_a': str(size),
             'shell_curl': str(sp.factor(sp.limit(kc, b, size))),
             'shell_symmetric': str(sp.factor(sp.limit(ks, b, size)))}
            for size in (2, 4, 8, sp.oo)
        ],
    }


if __name__ == '__main__':
    for dimension in (2, 3):
        print(json.dumps(test_shell(dimension)), flush=True)
