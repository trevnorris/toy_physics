"""Angular momentum and charge bookkeeping for mirror-pair formation.

Take both spin vectors in common brane laboratory coordinates. The orthogonal
matrix M sends the label of a patch on mouth 1 to a patch on mouth 2. It is
not a component transformation of the spin being summed. The constants l_i
are spin per unit solid angle or include the matched patch measure. The spin
test concerns this texture only. The charge law n_i=-s_i*chi_i is supplied.
"""

from itertools import product

import sympy as sp


def main():
    x, phi, theta = sp.symbols("x phi theta", real=True)
    rhat = sp.Matrix([sp.sqrt(1-x**2)*sp.cos(phi),
                      sp.sqrt(1-x**2)*sp.sin(phi), x])

    def average(expr):
        azimuth = sp.integrate(expr, (phi, 0, 2*sp.pi))
        return sp.integrate(azimuth, (x, -1, 1))/(4*sp.pi)

    first = rhat.applyfunc(average)
    second = (rhat*rhat.T).applyfunc(average)
    print("unit-sphere direction mean, second moment:", first.T, second)
    l1, l2 = sp.symbols("l1 l2", nonnegative=True)
    # In axes along its rotation axis, any real orthogonal 3x3 matrix is a
    # proper rotation times its determinant. The moments are isotropic, so a
    # simultaneous change of lab basis leaves both the average and trace fixed.
    rotation = sp.Matrix([[sp.cos(theta), -sp.sin(theta), 0],
                          [sp.sin(theta), sp.cos(theta), 0], [0, 0, 1]])
    checks = [first == sp.zeros(3, 1), second == sp.eye(3)/3]
    for det_sign, handed_product in product((1, -1), repeat=2):
        M = det_sign*rotation
        # Factor out the unit overall orientation and handed sign on mouth 1.
        # The minus sign comes from the reverse orientation of mouth 2.
        difference = l1*sp.eye(3)-handed_product*l2*M
        mean = difference*first
        mismatch = sp.trigsimp(sp.trace(difference*second*difference.T))
        trace_form = l1**2+l2**2-2*l1*l2*handed_product*sp.trace(M)/3
        mismatch_matrix = sp.eye(3)-handed_product*M
        nonnegative_form = ((l1-l2)**2 +
                            l1*l2*sp.trace(mismatch_matrix.T*mismatch_matrix)/3)
        print(f"det M={det_sign} chi1*chi2={handed_product}, trace M, "
              "global residual, patch mean square:", sp.trace(M), mean.T,
              sp.factor(mismatch))
        checks.extend((sp.simplify(M.T*M) == sp.eye(3), mean == sp.zeros(3, 1),
                       sp.trigsimp(mismatch-trace_form) == 0,
                       sp.trigsimp(mismatch-nonnegative_form) == 0))

    # Enumerate the exact pair charges of the patch-cancelling branch with
    # l1=l2>0, choosing s1=+1, s2=-1; an overall flip reverses both charges.
    for label, M in (("identical", sp.eye(3)), ("antipodal", -sp.eye(3))):
        for chi2 in (1, -1):
            chi1, s1, s2 = 1, 1, -1
            spin_map = s1*chi1*l1*sp.eye(3)+s2*chi2*l2*M
            mismatch = sp.factor(sp.trace(spin_map*second*spin_map.T))
            n1, n2 = -s1*chi1, -s2*chi2
            print(f"{label} patches chi2={chi2}, mean square, n1, n2, total n:",
                  mismatch, n1, n2, n1+n2)
            same_relation = M == chi1*chi2*sp.eye(3)
            checks.append((mismatch.subs(l2, l1) == 0) == same_relation)

    # Charge isolation of each initially neutral birth, with no other charge
    # deposits, imposes a separate conservation law on its two daughters.
    def units(orientation, handedness):
        return -orientation*handedness

    neutral_pairs = [(units(1, left), units(-1, right))
                     for left, right in product((1, -1), repeat=2)
                     if units(1, left)+units(-1, right) == 0]
    relative_signs = set()
    for birth1, birth2 in product(neutral_pairs, repeat=2):
        same_side_product = birth1[0]*birth2[0]
        relative_signs.add(same_side_product)
        print("two neutral births (A+, A-), (B+, B-), same-side charge product:",
              birth1, birth2, same_side_product)
    checks.append(relative_signs == {1, -1})
    assert all(checks)


if __name__ == "__main__":
    main()
