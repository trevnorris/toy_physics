"""Subthreshold polarization and circulating current of one bound pair.

The loop response is computed as a ground energy curvature, a signed
probability current, and a cancellation of instantaneous/virtual responses.
"""
import numpy as np


def ground_ring(n, binding, flux, uniform_strain=0.):
    eye = np.eye(n)
    # A common real hopping perturbation on a particular bond, on both faces.
    t = np.ones(n)
    t[0] += uniform_strain
    shift = np.roll(eye,1,axis=0) @ np.diag(t)
    def face(sign):
        hop = -np.exp(1j*sign*flux/n)*shift
        h = hop+hop.conj().T
        dh = 1j*sign*(hop-hop.conj().T)/n
        return h, dh, -h/n**2
    top, bot = face(1), face(-1)
    h, j, h2 = [np.kron(ti,eye)+np.kron(eye,bi) for ti,bi in zip(top,bot)]
    h += np.diag((binding*(1-eye)).ravel())
    energy, states = np.linalg.eigh(h)
    ground = states[:,0]
    current = -(ground.conj()@j@ground).real
    weights = abs(states[:,1:].conj().T @ (j@ground))**2
    virtual = -2*sum(weights/(energy[1:]-energy[0]))
    instant = (ground.conj()@h2@ground).real
    return energy[0], current, instant, virtual


def main():
    for binding in (4.,16.):
        splitting = np.hypot(binding,4.)
        omega = (binding+splitting)/2
        # Matrix element of Q=n_top,R-n_bottom,R from even ground into odd state.
        q2 = (1-binding/splitting)/2
        susceptibility = -2*q2/omega
        # Direct occupation response at small opposite site-frequency offsets.
        t = np.array([[0.,-1.],[-1.,0.]])
        q = np.diag([0.,-1.,1.,0.])
        base = np.kron(t,np.eye(2))+np.kron(np.eye(2),t)+np.diag([0.,binding,binding,0.])
        def occupation(v):
            _, eig = np.linalg.eigh(base+v*q)
            return eig[:,0]@q@eig[:,0]
        eps = .0002
        direct = (occupation(eps)-occupation(-eps))/(2*eps)
        print(f"dimer binding={binding:g} dQ/d(v+hbar*adot)={susceptibility:.12f} "
              f"direct={direct:.12f} following_clock_drive=0")
        assert abs(direct-susceptibility)<1e-9
    for n, binding in ((4,4.),(4,8.),(4,16.),(6,16.)):
        e0, current0, instant, virtual = ground_ring(n,binding,0.)
        step = .01
        eplus, iplus, _, _ = ground_ring(n,binding,step)
        eminus, iminus, _, _ = ground_ring(n,binding,-step)
        di = (iplus-iminus)/(2*step)
        curvature = instant+virtual
        dd = (eplus+eminus-2*e0)/step**2
        # All hoppings remain real under common strain at zero flux.
        _, strain_current, _, _ = ground_ring(n,binding,0.,uniform_strain=.01)
        print(f"loop N={n} binding={binding:g} instantaneous={instant:.12f} "
              f"virtual={virtual:.12f} curvature={curvature:.12f} "
              f"energy_fd={dd:.12f} d_current_d_flux={di:.12f} "
              f"real_strain_current={strain_current:.3g}")
        assert abs(current0)<1e-12 and abs(strain_current)<1e-12
        assert curvature>0 and abs(di+curvature)<2e-5*curvature
        assert abs(dd-curvature)<3e-5*curvature+1e-9


if __name__ == "__main__":
    main()
