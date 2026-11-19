"""Local pairing on the two faces, its possible exterior response, and finite pair transport."""

import numpy as np
import sympy as sp


def paired_faces():
    angle = sp.symbols("angle", real=True)
    top = sp.Matrix([1, 0])
    bottom = sp.Matrix([sp.cos(angle/2), sp.sin(angle/2)])
    epsilon = sp.Matrix([[0, 1], [-1, 0]])
    pair = (top.T*epsilon*bottom)[0]
    overlap = top.dot(bottom)
    assert sp.trigsimp(pair**2-(1-sp.cos(angle))/2) == 0
    assert sp.trigsimp(overlap**2-(1+sp.cos(angle))/2) == 0
    print("squared pair and overlap at parallel, perpendicular, antiparallel:",
          [[sp.simplify(f.subs(angle,a)**2) for f in (pair,overlap)]
           for a in (0,sp.pi/2,sp.pi)])

    # Positive, identical transverse kinetic coefficient b and curl coefficient k on each face.
    # Sigma (cross-face pair) and delta (ordinary face coherence) are distinct stiffnesses.
    b, k, q = sp.symbols("b k q", positive=True)
    sigma, delta = sp.symbols("sigma delta", nonnegative=True)
    same, opposite = sp.Matrix([1,1]), sp.Matrix([1,-1])
    kernel = (k*q**2*sp.eye(2) + sigma*same*same.T + delta*opposite*opposite.T)/b
    waves = [sp.expand((vector.T*kernel*vector)[0]/2) for vector in (same,opposite)]
    print("squared frequencies, common face phase then opposite:", waves)
    assert sp.simplify(b*waves[0]-k*q**2-2*sigma) == 0
    assert sp.simplify(b*waves[1]-k*q**2-2*delta) == 0
    for s in (1,-1):
        flux = sp.Matrix([s,-s])
        print("face flux, pair obstruction, ordinary obstruction, top/bottom/equal probes:",
              list(flux), same.dot(flux), (-opposite).dot(flux),
              [flux[0], flux[1], same.dot(flux)])

    # Normalize a single-face monopole to total flux 2 pi s. Use
    # integral grad(1/r1).grad(1/r2) = 4 pi/d in the point-source exterior.
    gamma, d = sp.symbols("gamma d", positive=True)
    s1, s2 = sp.symbols("s1 s2")
    single = 2*gamma*(s1/2)*(s2/2)*4*sp.pi/d
    interaction = sp.simplify(single + single.subs({s1:-s1, s2:-s2}, simultaneous=True))
    print("two-face exterior interaction, outward force:", interaction, -sp.diff(interaction,d))


def pair_transfer():
    # L=patch 0, R=patch 1. A forward elementary hop carries exp(i phase).
    a, b, at, ab = sp.symbols("t_top t_bottom a_top a_bottom", real=True)
    gap, energy = sp.symbols("gap E", positive=True)
    for character in (1,-1):
        def hop(t,phase):
            return -t*sp.Matrix([[0,sp.exp(-sp.I*phase)], [sp.exp(sp.I*phase),0]])
        h = (sp.kronecker_product(hop(a,at),sp.eye(2))
             + sp.kronecker_product(sp.eye(2),hop(b,character*ab))
             + sp.diag(0,gap,gap,0))
        low, high = [0,3], [1,2]
        # Eliminate LR and RL exactly at E; in the strong-binding expansion E is small.
        eff = h.extract(low,high)*h.extract(high,low)/(energy-gap)
        forward = sp.simplify(eff[1,0])
        target = 2*a*b*sp.exp(sp.I*(at+character*ab))/(energy-gap)
        print("bottom character, virtual forward pair hop, diagonal shift:", character,
              forward, sp.simplify(eff[0,0]))
        assert sp.simplify(forward-target) == 0


def ring_energy(length, binding, top_twist, bottom_twist):
    # One component on each face. Separate patches pay the same binding gap.
    def motion(total_twist):
        shift = np.roll(np.eye(length), 1, axis=0)*np.exp(1j*total_twist/length)
        return -(shift+shift.conj().T)
    h = np.kron(motion(top_twist),np.eye(length)) + np.kron(np.eye(length),motion(bottom_twist))
    h += np.diag((binding*(1-np.eye(length))).ravel())
    return np.linalg.eigvalsh(h)[0]


def twist_experiment():
    twist = np.pi/3
    for length, binding in ((6,0),(6,4),(4,16),(6,16),(8,16),(10,16)):
        ground = ring_energy(length,binding,0,0)
        common = ring_energy(length,binding,twist,twist) - ground
        opposite = ring_energy(length,binding,twist,-twist) - ground
        print(f"periodic ring N={length} gap={binding} common={common:.12g}"
              f" opposite={opposite:.12g} ratio={opposite/common:.9g}")
        assert common > 0 and opposite > 0
        assert np.isclose(ring_energy(length,binding,twist,-twist),
                          ring_energy(length,binding,-twist,twist))
        if binding == 0:
            assert np.isclose(common,opposite)


if __name__ == "__main__":
    paired_faces()
    pair_transfer()
    twist_experiment()
