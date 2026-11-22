"""Independent zero-curvature sensitivity equations for a symmetric two-face strip."""
import json

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.special import expit


def susceptibility(width, extent):
    half = width / 2.
    susceptibility_half = 0.
    axial_weight_half = 0.
    for inner, (lo, hi) in enumerate(((0., half), (half, half + extent))):
        x = np.linspace(lo, hi, int(20 * (hi - lo)) + 1)
        kink = expit(half - x)
        guess = np.array([kink, -kink*(1.-kink), np.zeros(x.size), np.zeros(x.size)])

        def eq(x, y):
            c, cp, h, hp = y
            return np.stack((cp, c*(1.-c)*(1.-2.*c), hp, (1.-6.*c+6.*c*c)*h - cp))

        def ends(left, right):
            if inner == 0:
                return (left[1], right[0] - .5, left[2], right[2])
            return (left[0] - .5, right[0], left[2], right[2])

        result = solve_bvp(eq, ends, x, guess, tol=1.e-10, max_nodes=30000)
        if not result.success:
            raise ArithmeticError(result.message)
        def integrand(s):
            c, cp, h, hp = result.sol(s)
            return cp*h
        susceptibility_half += quad(integrand, lo, hi, epsabs=2.e-13)[0]
        axial_weight_half += quad(lambda s: s*s*result.sol(s)[0]**2, lo, hi, epsabs=2.e-12)[0]
    return -4.*susceptibility_half, 2.*axial_weight_half


if __name__ == '__main__':
    B, I = susceptibility(25., 20.)
    print(json.dumps(dict(H_over_delta=25, B_over_a_delta3=B, B_over_one_sigma_delta2=3.*B,
                          area_weighted_second_moment_delta3=I, weighted_moment=I*.04**3,
                          B_neutral_native=75.*.04**3*B, threshold_Y=-75.*B/I,
                          homogeneous_threshold_Y_over_a=-12.*B/25**3)))
