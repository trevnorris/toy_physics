"""First variation of sideways momentum accounting on a fixed brane cut."""

import sympy as sp


def main():
    eps, t, x, y, z = sp.symbols("epsilon t x y z", real=True)
    coords = (x, y, z)
    def pair(name):
        base = sp.Function(name + "0")(*coords)
        change = sp.Function("d" + name)(t, *coords)
        return base + eps*change, base, change

    g, _, dg = pair("g")
    j, j0, dj = pair("j")
    v, v0, dv = zip(*(pair(name) for name in ("Vx", "Vy", "Vz")))
    source, _, dsource = pair("S")
    flux, _, dflux = zip(*(pair("C"+k) for k in ("x","y","z")))
    stress, _, dstress = zip(*(pair("sigma"+k) for k in ("x","y","z")))

    def div(vector):
        return sum(sp.diff(f, coord) for f, coord in zip(vector, coords))

    def variation(expr):
        return sp.diff(expr, eps).subs(eps, 0)

    # C is spatial transport excluding tensile stress; S excludes the net
    # rim traction (applied on the boundary) and the outward carried loss.
    residual = sp.diff(g,t) + div([f-p for f,p in zip(flux,stress)]) + j*v[0] - source
    expected = (sp.diff(dg,t) + div([f-p for f,p in zip(dflux,dstress)])
                + j0*dv[0] + v0[0]*dj - dsource)
    general = sp.expand(variation(residual) - expected)
    print("Fixed-cut first variation minus storage + transport - stress + carried loss - other supply:", general)

    # Conditional test only: ordinary advected mass momentum on a flat fixed
    # measure, with carried mass leaving at its local material velocity.
    rho, rho0, drho = pair("rho")
    mass = sp.diff(rho,t) + div([rho*w for w in v]) + j
    momentum = sp.diff(rho*v[0],t) + div([rho*v[0]*w for w in v]) + j*v[0]
    acceleration0 = sum(w*sp.diff(v0[0],c) for w,c in zip(v0,coords))
    dacceleration = (sp.diff(dv[0],t)
                     + sum(w*sp.diff(dv[0],c) for w,c in zip(v0,coords))
                     + sum(w*sp.diff(v0[0],c) for w,c in zip(dv,coords)))
    conditional = sp.expand(variation(momentum-v[0]*mass) - drho*acceleration0 - rho0*dacceleration)
    print("Conditional advected-momentum tangent minus density, material and background-gradient terms:", conditional)

    # Two example incremental relaxation laws, not parameter choices for the drain.
    s, tau, modulus = sp.symbols("s tau G_inf", positive=True)
    def shear_tangent(relaxation):
        return modulus*s*relaxation/(1+s*relaxation)
    slopes = [shear_tangent(relaxation) for relaxation in (tau, 2*tau)]
    low = [sp.limit(value,s,0) for value in slopes]
    high = [sp.limit(value,s,sp.oo) for value in slopes]
    finite = [sp.factor(value.subs(s,1/tau)) for value in slopes]
    print("Shear-stress increments per strain increment: relaxed, instantaneous, at s=1/tau:", low, high, finite)
    if general != 0 or conditional != 0 or low != [0,0] or high != [modulus]*2 or finite[0] == finite[1]:
        raise ArithmeticError("Tangent identities or example time-domain distinction failed")


if __name__ == "__main__":
    main()
