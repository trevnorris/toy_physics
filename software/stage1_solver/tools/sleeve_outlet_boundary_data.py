"""Local bulk carrier on the discharge side of a flat section, at 0PN.

Find decaying, steady perturbations with a common downstream state and a
prescribed normal number current at the section. The background drains out
along +w. These dimensionless bulk solutions do not choose a throat wall.
"""
import sympy as sp


x, w = sp.symbols("x w", real=True)
eps, U = sp.symbols("epsilon U", real=True)
beta = sp.Symbol("beta", positive=True)  # hbar**2/(4*m**2)
k, lam, v, c2 = sp.symbols("k lambda V c_s_squared", positive=True)
f, g = sp.symbols("f g")


def gradient(field):
    return sp.Matrix([sp.diff(field, coord) for coord in (x, w)])


def laplacian(field):
    return sum(sp.diff(field, coord, 2) for coord in (x, w))


# n0=m=1; h(n)/m=c_s^2*n^4/4. The other two brane coordinates are spectators.
def linear_equations(fractional_density, potential_change):
    density = 1 + eps*fractional_density
    potential = v*w + eps*potential_change
    velocity = gradient(potential)
    current = density*velocity
    continuity = sp.diff(current[0], x) + sp.diff(current[1], w)
    bernoulli = (velocity.dot(velocity)/2 + c2*density**4/4
                 - 2*beta*laplacian(sp.sqrt(density))/sp.sqrt(density))
    return sp.Matrix([sp.diff(eq, eps).subs(eps, 0) for eq in (continuity, bernoulli)])


def main():
    params = {v: 1, c2: sp.Rational(35, 24), beta: sp.Rational(1, 24)}
    local_params = {**params, k: 1}
    mode = sp.sin(k*x)*sp.exp(-lam*w)
    equations = linear_equations(f*mode, g*mode)
    matrix = equations.applyfunc(lambda term: sp.cancel(term/mode)).jacobian([f, g])
    print("normal mode matrix [fractional density, potential] =", matrix)
    print("sample characteristic determinant =", sp.factor(matrix.det().subs(local_params)))

    def decaying_modes(wave_number):
        chosen = matrix.subs({**params, k: wave_number})
        roots = sp.solve(chosen.det(), lam)
        modes = []
        for root in roots:
            if root.is_positive:
                vector, = chosen.subs(lam, root).nullspace()
                modes.append((root, vector/vector[0]))
        return modes

    # k=0 uses a brane-uniform perturbation with the same matrix at k=0.
    (uniform_rate, uniform_vector), = decaying_modes(0)
    uniform = uniform_vector*sp.exp(-uniform_rate*w)
    iso_residual = linear_equations(*uniform).subs(params).applyfunc(sp.simplify)
    print("isolated unit density mode [density, potential] =", list(uniform))
    print("isolated [continuity, Bernoulli, normal current change] =",
          list(iso_residual) + [sp.simplify(uniform[0]+sp.diff(uniform[1], w))])

    modes = decaying_modes(1)
    # Choose a density amplitude convention for one mode. Solve the normal
    # flux constraint for the other amplitude.
    coefficients = sp.symbols(f"a0:{len(modes)}")
    density, potential = sum(
        (coeff*vec*sp.exp(-rate*w)*sp.sin(x)
         for coeff, (rate, vec) in zip(coefficients, modes)), sp.zeros(2, 1))
    flux_at_cut = (density + sp.diff(potential, w)).subs(w, 0)
    amplitude, = sp.linsolve([flux_at_cut/sp.sin(x), coefficients[-1]-8], coefficients)
    normalization = dict(zip(coefficients, amplitude))
    density = density.subs(normalization)
    potential = potential.subs(normalization)
    # Add the incident constant sideways velocity U. The unknown response b
    # vanishes with the companion forcing when multiplied by epsilon.
    b = sp.Symbol("b", real=True)
    perturbation = U*x + b*potential
    residual = linear_equations(b*density, perturbation).subs(params).applyfunc(sp.simplify)
    velocity = gradient(perturbation)
    print("paired unit decaying mode [density, potential] =", [density, potential])
    print("paired [continuity, Bernoulli] =", list(residual))
    fields = [b*density, *velocity, b*density + velocity[1]]
    print("at cut [density, sideways velocity, normal velocity, normal current] change =",
          [sp.simplify(field.subs(w, 0)) for field in fields])
    print("downstream changes =", [sp.limit(field, w, sp.oo) for field in fields])
    if any(iso_residual) or any(residual):
        raise ValueError("A linearized field equation has a nonzero residual.")


if __name__ == "__main__":
    main()
