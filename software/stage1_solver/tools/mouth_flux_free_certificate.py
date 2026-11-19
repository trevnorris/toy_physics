"""Independent free-angular-flux comparison for a three-dimensional collar.

All energies are in 4*pi*K**2/(u*R). The radial variable is t=r/R and
the density unit is K/(u*R**2). Angular brackets denote unit-sphere means.
The field coefficient gamma*u/K**2 multiplies the full exterior Maxwell
energy, including the Coulomb contribution (gamma*u/K**2)*n**2/4.
"""

import sympy as sp


def main():
    Q = sp.Rational
    aspect, strong, weak, gamma = Q(5, 3), Q(9, 5), Q(1, 10), Q(1, 2)
    t, z = sp.symbols("t z", positive=True)
    length, volume, inverse = [sp.integrate(w, (t, 1, aspect))
                               for w in (1, t**2, t**-2)]
    c = length/volume

    # At critical Maxwell stiffness the completed square subtracts the
    # variance of the signed coherence density from the local quartic. A
    # helicity/sign reversal can always be absorbed by changing the sign of b.
    b, p, mean_p = sp.symbols("b p mean_p", real=True)
    rs, rw = sp.symbols("rho_s rho_w", nonnegative=True)
    square_residual = sp.expand(
        gamma*b**2/t**2 + b*(p-mean_p) + t**2*(p-mean_p)**2/2
        - (b+t**2*(p-mean_p))**2/(2*t**2))
    density_surplus = sp.expand((rs+rw)**2-(rs-rw)**2)
    print("angular square residual, two-component quartic surplus:",
          square_residual, density_surplus)

    a = sp.symbols("a", positive=True)
    # Require (t^2 v)'=a*(c*t^2-1) and both free-edge calibration values zero.
    tau = sp.symbols("tau", positive=True)
    radial_shift = sp.factor(a*sp.integrate(c*tau**2-1, (tau, 1, t))/t**2)
    # 4*aspect*(t-1)*(aspect-t)/t^2 <= (aspect-1)^2; multiply by the
    # maximum of (t+aspect+1)/(1+aspect+aspect^2) for the shift bound.
    interval_square = sp.factor((aspect-1)**2*t**2 -
                                4*aspect*(t-1)*(aspect-t))
    radial_max = (a*(aspect-1)**2/(4*aspect) *
                  (2*aspect+1)/(1+aspect+aspect**2))
    print("radial moments (length, volume, inverse), c:", length, volume, inverse, c)
    print("calibration p_a, boundary values, (t^2*p_a)'/a:",
          radial_shift, radial_shift.subs(t, 1), radial_shift.subs(t, aspect),
          sp.factor(sp.diff(t**2*radial_shift, t)/a))
    print("interval square, upper bound on |p_a|:", interval_square, radial_max)

    # Opposite shifts on strong and weak radial kinetic terms give the
    # signed mean U=x-y an unconstrained square and give y a positive cost.
    x, y = sp.symbols("x y", nonnegative=True)
    d, v = sp.symbols("d v", real=True)
    U = x-y
    dp = a*(c*t**2-1)
    material_bound = t**2*U**2/2-a*U+d*y -dp*U-t**2*v**2*(x+y)
    completed = (t**2*(U-a*c-v**2)**2/2-t**2*(a*c+v**2)**2/2 +
                 (d-2*t**2*v**2)*y)
    square_check = sp.expand(material_bound-completed)
    print("calibrated signed-mean square residual:", square_check)

    # Choose orientation so the strong native degree is +2. The reflected
    # problem maps n to -n. Only a trial for the selected sign is needed.
    floors = {}
    weak_costs = []
    no_radial_floors = {}
    for n in (0, -1):
        deg_s, deg_w = n-2, n+2
        gain_s = strong-Q(abs(deg_s), 2)
        gain_w = weak-Q(abs(deg_w), 2)
        leftover = -gain_s-gain_w
        p_max = radial_max.subs(a, gain_s)
        weak_cost = leftover-2*aspect**2*p_max**2
        floors[n] = gamma*n**2/4-volume*(gain_s*c+p_max**2)**2/2
        no_radial_floors[n] = gamma*n**2/4-inverse*gain_s**2/2
        print(f"n={n} strong/weak degrees, gains, d, weak cost lower, E floor:",
              deg_s, deg_w, gain_s, gain_w, leftover, weak_cost, floors[n])
        weak_costs.append(weak_cost)

    # For |n|>=2 discard the positive angular gradient. Keep both material
    # gradients with the same shift and use the full unscreened quartic with
    # Jensen on the total density. Uniform Coulomb Maxwell cost is minimized
    # by |n|=2 in this class.
    p_tail = radial_max.subs(a, strong)
    floors[2] = gamma-volume*(strong*c+p_tail**2)**2/2
    no_radial_floors[2] = gamma-inverse*strong**2/2
    print("|n|>=2 E floor, all floors without radial calibration:",
          floors[2], *(no_radial_floors[n] for n in (0, -1, 2)))

    # On the favored unit sector a degree-one section in the uniform radial
    # monopole connection has chart amplitude sqrt(z)*(1+(1-z)/3). The sign
    # of its effective degree is handled by conjugation. Keep it constant in
    # radius, minimize just its overall scale, and leave the weak field zero.
    f = sp.sqrt(z)*(1+(1-z)/3)
    angular_term = z*(1-z)*sp.diff(f, z)**2+(1-z)*f**2/(4*z)
    T, M, P4 = [sp.integrate(sp.cancel(h), (z, 0, 1))
                for h in (angular_term, f**2, f**4)]
    amplitude = (strong*M-T)*length/(volume*P4)
    trial = gamma/4+length*(T-strong*M)*amplitude+volume*P4*amplitude**2/2
    print("favored n=1 trial (T, M, P4, amplitude squared, E):",
          T, M, P4, amplitude, trial)
    print("certified gaps (neutral, wrong unit, |n|>=2):",
          *(floors[n]-trial for n in (0, -1, 2)))

    assert square_check == square_residual == 0 and density_surplus == 4*rs*rw
    assert interval_square == sp.factor(((aspect+1)*t-2*aspect)**2)
    assert radial_shift.subs(t, 1) == radial_shift.subs(t, aspect) == 0
    assert sp.factor(sp.diff(t**2*radial_shift, t)-dp) == 0
    assert all(cost > 0 for cost in weak_costs)
    assert amplitude > 0 and trial < min(floors.values())


if __name__ == "__main__":
    main()
