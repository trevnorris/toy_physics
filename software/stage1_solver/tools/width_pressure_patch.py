"""Static calibration for equal-and-opposite pressure on two flat faces."""

import sympy as sp


def main():
    radius, length, stiffness, gradient, tension = sp.symbols("r L K D sigma", positive=True)
    width_slope, center_slope, eps = sp.symbols("b_x Z_x eps", real=True)
    # A slope along x suffices for isotropic area tension. Faces sit at Z +/- W/2.
    slopes = (center_slope + width_slope/2, center_slope - width_slope/2)
    area_energy = sum(tension * (sp.sqrt(1 + eps**2*slope**2) - 1) for slope in slopes)
    quadratic = sp.expand(sp.diff(area_energy, eps, 2).subs(eps, 0)/2)
    d_width = sp.diff(quadratic, width_slope, 2)
    print("area energy to quadratic order =", quadratic)
    print("width gradient coefficient =", d_width)
    assert d_width == tension/2

    # Use the radial Green function, normalized by the inward flux of its
    # gradient at a small surrounding sphere, in three brane dimensions.
    decay = sp.sqrt(stiffness/gradient)
    green = sp.exp(-decay*radius)/(4*sp.pi*gradient*radius)
    normal_flux = -4*sp.pi*radius**2*gradient*sp.diff(green, radius)
    assert sp.simplify(sp.limit(normal_flux, radius, 0, dir="+")-1) == 0
    response = sp.simplify(sp.integrate(4*sp.pi*radius**2*green, (radius, 0, length)))
    expected = (1-(1+length*decay)*sp.exp(-length*decay))/stiffness
    print("center width per unit pressure =", response)
    soft = sp.simplify(sp.limit(response, stiffness, 0, dir="+"))
    print("zero local stiffness =", soft)
    print("large patch at positive local stiffness =", sp.limit(response, length, sp.oo))
    assert sp.simplify(response-expected) == 0 and soft == length**2/(2*gradient)
    assert sp.limit(response, length, sp.oo) == 1/stiffness
    doubled = sp.simplify(soft.subs(length, 2*length)/soft)
    print("width ratio at radii 2L and L for tension alone =", doubled)
    assert doubled == 4


if __name__ == "__main__":
    main()
