"""Linear response of a mirror-symmetric slab in a fixed tangential lab frame.

The physical reflection w -> -w interchanges the two faces, without an extra
internal conjugation. Apply the face transformation to transverse connections
(or their curls); pure gauge shifts of a connection drop out of a flux probe.
All examples refer to nonzero in-brane wave number and a *local* curvature
component, measurable by a signed oriented open patch. A source-free wave has
zero total flux through a closed surface.
"""

import sympy as sp


def parity_blocks():
    # Work at an arbitrary fixed (k, omega), for one tangential polarization.
    # Rows are linear readouts and columns drives. No reciprocity, Hermiticity
    # or locality constraint is imposed on these complex response functions.
    swap = sp.Matrix([[0, 1], [1, 0]])
    mirror = sp.diag(swap, swap)              # (u_top, u_bottom, A_top, A_bottom)
    raw = sp.Matrix(4, 4, sp.symbols("g0:16"))
    response = (raw + mirror * raw * mirror.T) / 2
    parity = sp.Matrix([[1, 0, 1, 0], [1, 0, -1, 0],
                        [0, 1, 0, 1], [0, 1, 0, -1]]) / sp.sqrt(2)
    projected = sp.simplify(parity.T * response * parity)
    print("mirror eigenvalues in (u+,A+,u-,A-) order:", list((parity.T * mirror * parity).diagonal()))
    print("arbitrary response even-to-odd, odd-to-even:", projected[2:, :2], projected[:2, 2:])
    print("arbitrary response surviving coefficients per 2x2 diagonal block:",
          sum(c != 0 for c in projected[:2, :2]), sum(c != 0 for c in projected[2:, 2:]))
    assert (parity.T * mirror * parity) == sp.diag(1, 1, -1, -1)
    assert projected[2:, :2] == projected[:2, 2:] == sp.zeros(2)


def degenerate_residue():
    # At a coincident pole a rotated eigenbasis can give both drives nonzero
    # overlap with an individual basis vector. The projector onto the full
    # pole is invariant, hence the sum of cross residues is the observable.
    k = sp.symbols("k", positive=True)
    angle = sp.symbols("angle", real=True)
    rotated = sp.rot_axis3(angle)[:2, :2]
    drive = sp.Matrix([1, 0])
    flux = sp.Matrix([0, k])  # one curl helicity, eigenvalue +k
    residues = [(flux.dot(v) * v.dot(drive)).expand(trig=True) for v in rotated.columnspace()]
    print("coincident pole cross weights at pi/4, sum at general angle:",
          [r.subs(angle, sp.pi/4) for r in residues], sp.trigsimp(sum(residues)))
    assert sp.trigsimp(sum(residues)) == 0


def slab_profile():
    H, w, mu_w, rho = sp.symbols("H w mu_w rho", positive=True)
    even = sp.Integer(1)
    odd = sp.sin(sp.pi * w / H)
    interval = (w, -H/2, H/2)
    cross = sp.integrate(even * odd, interval)
    gap_squared = sp.simplify(mu_w * sp.integrate(sp.diff(odd,w)**2, interval)
                              / (rho * sp.integrate(odd**2, interval)))
    print("free homogeneous slab even-odd overlap, odd face slopes, odd squared gap:",
          cross, [sp.diff(odd, w).subs(w, face) for face in interval[1:]], gap_squared)
    assert cross == 0 and gap_squared == sp.pi**2 * mu_w / (H**2 * rho)


def mechanical_control():
    # Real amplitudes of a transverse curl helicity with eigenvalue q = +/-|k|.
    # u is even tangential translation, v an odd internal shear amplitude,
    # A the odd connection. The illustrative curl coupling works when v and
    # curl A transform alike under the other symmetries. Its positive energy
    # gives the strain-gradient energy to A and a thickness gap to v.
    u, v, A = sp.symbols("u v A", real=True)
    k, mass, ce, cv, ca = sp.symbols("k m c_e c_v c_A", positive=True)
    coupling, z = sp.symbols("lambda z", real=True)
    bias = sp.symbols("b", real=True)
    for q in (k, -k):
        energy = (ce**2 * k**2 * u**2 + cv**2 * k**2 * v**2
                  + ca**2 * k**2 * A**2 + mass**2 * (v - coupling * q * A)**2) / 2
        stiffness = sp.hessian(energy, (u, v, A))
        dynamic = stiffness - z * sp.eye(3)    # positive unit inertias, z=omega^2
        source = sp.Matrix([1, bias, 0])
        probe = sp.Matrix([0, 0, q])
        odd_det = dynamic[1:, 1:].det()
        transfer = sp.factor((probe.T * dynamic.inv() * source)[0])
        numerator = sp.factor(transfer * odd_det)
        sample = {k: sp.Rational(1,2), mass: 2, ce: 1, cv: 1, ca: 1,
                  coupling: sp.Rational(1,2), z: 0}
        sample_readings = [transfer.subs({**sample, bias:b}) for b in (0, 1, -1)]
        print("helicity, whole-slab transfer, biased transfer times odd determinant, "
              "sample (whole,upper,lower):", q/k, transfer.subs(bias,0), numerator, sample_readings)
        assert transfer.subs(bias, 0) == 0
        assert numerator == bias * coupling * k**2 * mass**2
        assert sample_readings == [0, sp.Rational(4,9), -sp.Rational(4,9)]

    strain, field, g = sp.symbols("S B g", real=True)
    energy = g * strain * field**2
    mixed = sp.diff(energy, strain, field)
    print("even-strain odd-flux-squared energy, mixed derivative at zero and finite flux:",
          energy, mixed.subs(field,0), mixed)
    assert energy == energy.subs(field, -field) and mixed.subs(field, 0) == 0


if __name__ == "__main__":
    parity_blocks()
    degenerate_residue()
    slab_profile()
    mechanical_control()
