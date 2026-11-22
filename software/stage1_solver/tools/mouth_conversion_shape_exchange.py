"""Free shape perturbations at fixed axial length and exactly matched volume.

A wave with any smooth energy E(V,L) can exert lateral pressure and can
restore a uniform change of volume. Test shape redistribution invisible to
that energy, using the 4D isotropic membrane, on normalized turning profiles.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.optimize import brentq

from mouth_conversion_shape_mode import profile


def test_shape(radius, ratio, functions, order):
    p_unit, _, unit_length, _, solution = profile(ratio)
    length = radius*unit_length
    x, q = leggauss(order)
    w = length*(1+x)/2
    dw = q*length/2
    b, slope = solution(w/radius)
    b = radius*b
    arc = np.hypot(1., slope)
    darea = 4*np.pi*b*b*arc*dw
    k_angular = 1/(b*arc)
    k_meridian = p_unit/radius-2*k_angular
    k_squared = k_meridian**2 + 2*k_angular**2
    modes = np.arange(1, functions+1)[:, None]
    k = modes*np.pi/length
    mean = (np.cos(k*w)@darea)/sum(darea)
    basis = np.cos(k*w)-mean[:,None]
    derivative = -k*np.sin(k*w)
    gram = (basis*darea)@basis.T
    stiffness = ((derivative*(darea/arc**2))@derivative.T
                 -(basis*(darea*k_squared))@basis.T)
    values, vectors = eigh(stiffness, gram, subset_by_index=[0,0])
    vector = vectors[:,0]
    # Fix normalization to unit outward amplitude at the waist.
    vector /= vector@(1-mean)
    f = vector@basis
    df = vector@derivative
    stiffness_of_f = vector@stiffness@vector
    np.testing.assert_allclose(f@darea, 0., atol=2e-10)

    def surface(step, offset):
        h = step*f+offset
        dh = step*df
        new_b = b+h/arc
        axial_jacobian = 1+h*k_meridian-dh*slope/arc
        metric = np.hypot(arc*(1+h*k_meridian), dh)
        if min(new_b)<=0 or min(axial_jacobian)<=0:
            raise ArithmeticError('Crossed axis or lost regular meridian')
        return (dw@(4*np.pi*new_b**2*metric),
                dw@((4*np.pi/3)*new_b**3*axial_jacobian))

    area, volume = surface(0., 0.)
    controls = []
    for fraction in (.01, .005):
        energies = []
        errors = []
        offsets = []
        step = fraction*radius
        for sign in (-1,1):
            h = step*sign
            correction = brentq(lambda c: surface(h,c)[1]-volume, -abs(h), abs(h), xtol=1e-14)
            moved_area, moved_volume = surface(h,correction)
            energies.append(moved_area-area)
            errors.append((moved_volume-volume)/volume)
            offsets.append(correction)
        controls.append(dict(fraction=fraction, volume_errors=errors, offsets=offsets,
                             area_changes=energies, curvature=sum(energies)/step**2))
    if max(abs(e) for c in controls for e in c['volume_errors'])>1e-12:
        raise ArithmeticError('Signed readout left the constant-volume manifold')
    np.testing.assert_allclose(controls[-1]['curvature'],stiffness_of_f,rtol=.01,atol=1e-5)
    return dict(waist=radius,exit_to_waist=ratio,free_length=length,area=area,
                volume=volume,cosines=functions,quadrature=order,
                lowest_Ritz_quotient=values[0],area_stiffness=stiffness_of_f,
                waist_displacement=1.,exit_displacement=float(vector@((-1.)**modes[:,0]-mean)),
                coefficients=vector.tolist(),controls=controls)


if __name__=='__main__':
    rows=[test_shape(1.,ratio,6,96) for ratio in (1.25,1.5,2.)]
    rows.extend(test_shape(r,1.5,10,144) for r in (1.,2.))
    print(json.dumps(rows,indent=2,allow_nan=False))
