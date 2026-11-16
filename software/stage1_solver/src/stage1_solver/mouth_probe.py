"""Apply virtual work and read the center of a star-shaped 3-D aperture.

The rim is sampled on a sphere at fixed w. The angular samples represent an
axisymmetric signal about the x force axis; weights integrate d(cos theta).
"""

import numpy as np


def aperture_center_x(radius, mu, weights):
    """Volume centroid of the aperture, integrating radial rays to the rim."""
    radius, mu, weights = (np.asarray(v, dtype=float) for v in (radius, mu, weights))
    if (radius.ndim != 1 or mu.shape != radius.shape or weights.shape != radius.shape
            or not np.all(np.isfinite([radius, mu, weights])) or np.any(radius <= 0)
            or np.any(abs(mu) > 1) or np.any(weights <= 0)):
        raise ValueError("Invalid aperture radii or angular quadrature")
    # The radial integrals for volume and its first moment are R^3/3 and R^4/4.
    return 0.75 * np.dot(weights * mu, radius**4) / np.dot(weights, radius**3)


def spherical_rim_probe(radius, mu, force):
    """Signed normal load per rim area at the initial spherical mouth.

    External work is force times variation of aperture_center_x. At a sphere
    the variation of the centroid is 3/(4*pi) integral mu*delta_R dOmega.
    """
    mu = np.asarray(mu, dtype=float)
    if (not np.isfinite(radius) or radius <= 0 or not np.isfinite(force)
            or not np.all(np.isfinite(mu)) or np.any(abs(mu) > 1)):
        raise ValueError("Expected a positive radius and finite force and direction cosines")
    return 3 * force * mu / (4 * np.pi * radius**2)
