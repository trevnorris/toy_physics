"""
em_breathing_gauss_check.py

Purpose:
    Numerically verify that the static, spherically symmetric breathing-mode
    potential outside a compact core behaves like a Coulomb field:

        phi(r) ~ Q / (4 * pi * eps0 * r),

    and that the electric flux

        Phi_E(R) = 4 * pi * R^2 * E_r(R),  E_r = -dphi/dr,

    is independent of R and equal to Q / eps0 (Gauss's law).

Model:
    Outside the defect core, the static breathing-mode potential Phi(r)
    solves the l=0 Laplace equation:

        0 = ∇^2 Phi = Phi''(r) + (2/r) Phi'(r),

    whose nontrivial solution is A + B/r. Dropping the constant A, this
    gives the usual 1/r profile. We integrate this ODE numerically with
    initial conditions chosen to match the analytic Coulomb solution at
    a small radius r0.

Dependencies:
    mpmath  (usually comes with SymPy)
"""

import mpmath as mp


def make_radial_solution(Q=1.0, eps0=1.0, r0=0.1):
    """
    Construct a numerical solution Phi(r) of

        Phi''(r) + 2/r * Phi'(r) = 0

    for r >= r0, with initial data chosen so that

        Phi(r) = Q / (4*pi*eps0*r)

    at r = r0.

    Returns:
        sol(r): a callable that returns (phi(r), dphi/dr).
    """

    def ode(r, y):
        phi, dphi = y
        # Laplace equation in spherical symmetry: phi'' + 2/r phi' = 0
        return [dphi, -2 * dphi / r]

    # Analytic Coulomb-like initial conditions at r0:
    phi0 = Q / (4 * mp.pi * eps0 * r0)
    dphi0 = -Q / (4 * mp.pi * eps0 * r0**2)

    sol = mp.odefun(ode, r0, (phi0, dphi0))
    return sol


def sample_and_print(Q=1.0, eps0=1.0, r0=0.1):
    """
    Sample the numerical solution at a few radii and print:

        r, phi(r), E_r(r) = -dphi/dr, and
        Phi_E(r) = 4*pi*r^2*E_r(r).

    For a perfect Coulomb field we should have:

        phi(r) = Q / (4*pi*eps0*r),
        Phi_E(r) = Q / eps0   (independent of r).
    """
    sol = make_radial_solution(Q=Q, eps0=eps0, r0=r0)

    radii = [0.2, 0.5, 1.0, 2.0, 5.0]
    print("Breathing-mode / Coulomb check")
    print(f"Q = {Q}, eps0 = {eps0}, r0 = {r0}\n")
    print("{:>6}  {:>14}  {:>14}  {:>14}".format("r", "phi(r)", "E_r(r)", "Phi_E(r)"))
    print("-" * 60)

    for r in radii:
        phi, dphi = sol(r)
        E_r = -dphi
        flux = 4 * mp.pi * r**2 * E_r
        # Cast to float for printing
        print("{:6.3f}  {:14.8e}  {:14.8e}  {:14.8e}".format(
            float(r),
            float(phi),
            float(E_r),
            float(flux),
        ))

    print("\nExpected Gauss-law flux: Q/eps0 = {:.8e}".format(float(Q / eps0)))
    print("If Phi_E(r) is approximately constant and equal to Q/eps0, "
          "the numeric solution behaves like a Coulomb field.\n")


if __name__ == "__main__":
    # You can tweak Q, eps0, r0 if you like; keep r0 small but not too
    # close to 0 to avoid numerical issues from the 2/r term.
    sample_and_print(Q=1.0, eps0=2.0, r0=0.1)


"""
Output:

Breathing-mode / Coulomb check
Q = 1.0, eps0 = 2.0, r0 = 0.1

     r          phi(r)          E_r(r)        Phi_E(r)
------------------------------------------------------------
 0.200  1.98943679e-01  9.94718394e-01  5.00000000e-01
 0.500  7.95774715e-02  1.59154943e-01  5.00000000e-01
 1.000  3.97887358e-02  3.97887358e-02  5.00000000e-01
 2.000  1.98943679e-02  9.94718394e-03  5.00000000e-01
 5.000  7.95774715e-03  1.59154943e-03  5.00000000e-01

Expected Gauss-law flux: Q/eps0 = 5.00000000e-01
If Phi_E(r) is approximately constant and equal to Q/eps0, the numeric solution behaves like a Coulomb field.
"""
