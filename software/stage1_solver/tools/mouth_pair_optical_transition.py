"""Time-dependent material spectroscopy on the existing two-site bound pair.

hbar=1. A rightward hop has phase +a on the top and -a on the bottom.
The site-frequency difference is a separate observable. Also test the
specific local material clock induced by projection onto a moving spinor.
"""
import numpy as np
import sympy as sp


def pair_test():
    t, gap, a = sp.symbols("t Delta a", positive=True, real=True)
    h0 = sp.Matrix([[0,-t],[-t,0]])
    ident = sp.eye(2)
    base = (sp.kronecker_product(h0,ident)+sp.kronecker_product(ident,h0)
            +sp.diag(0,gap,gap,0))
    q = sp.diag(0,-1,1,0)
    def hop(angle):
        return sp.Matrix([[0,-t*sp.exp(-sp.I*angle)],[-t*sp.exp(sp.I*angle),0]])
    physical = (sp.kronecker_product(hop(a),ident)+sp.kronecker_product(ident,hop(-a))
                +sp.diag(0,gap,gap,0))
    u = sp.diag(*(sp.exp(sp.I*a*v) for v in q.diagonal()))
    equivalence = sp.simplify(physical-u*base*u.conjugate().T)
    angular = physical.diff(a).subs(a,0)
    # Ground/even relative mixture: order LL, LR, RL, RR.
    columns = sp.Matrix([[1,0,1,0],[0,1,0,1],[0,1,0,-1],[1,0,-1,0]])/sp.sqrt(2)
    block = columns.T*base*columns
    print("basis=bound+,separated+,bound-,separated- H0=",block.tolist())
    print("Q in same basis=",(columns.T*q*columns).tolist())
    print("physical link/current commutator and frame identities=",
          equivalence==sp.zeros(4), sp.simplify(angular-sp.I*(q*base-base*q))==sp.zeros(4))
    assert equivalence==sp.zeros(4) and angular==sp.I*(q*base-base*q)
    for binding in (4,16):
        subs = {t:1,gap:binding}
        energies, vectors = np.linalg.eigh(np.array(base.subs(subs),float))
        ground = vectors[:,0]
        phase = np.array(angular.subs(subs),complex)
        number = np.array(q,float)
        # Fractional real hopping modulation from strain; its calibration
        # to a physical bond length is d ln(t)/d strain, not supplied here.
        strength = np.array((t*base.diff(t)).subs(subs),float)
        for j in range(1,4):
            frequency = energies[j]-energies[0]
            state = vectors[:,j]
            # e^{-i omega t}, scalar voltage = -eta * dot(a).
            transition = [abs(state.conj() @ (phase+1j*eta*frequency*number) @ ground)**2
                          for eta in (0,.5,1)]
            length_signal = abs(state.conj()@strength@ground)**2
            print(f"Delta/t={binding} final={j} gap/t={frequency:.9f} "
                  f"real_hop_strength={length_signal:.9f} phase_fixed_sites={transition[0]:.9f} "
                  f"phase_half_follow={transition[1]:.9f} phase_full_follow={transition[2]:.3g}")
            assert transition[2]<1e-24
            if j==1:
                assert max(length_signal,transition[0])<1e-24
            if j==2:
                assert transition[0] > .1 and length_signal < 1e-24


def projected_clock():
    # Tilt the left normal from +z toward x. The fixed right normal tilts
    # from +z toward y. Each local time clock is -i z^dagger dot(z).
    for beta in (0., .4):
        eps = 1e-5
        left = lambda theta: np.array([np.cos(theta/2),np.sin(theta/2)],complex)
        right = np.array([np.cos(beta/2),1j*np.sin(beta/2)])
        for mirror in (1,-1):
            def registered(theta, rephase):
                zi, zj = left(theta), right.copy()
                if mirror == -1:
                    zi, zj = np.array([zi[1].conj(),-zi[0].conj()]), np.array([zj[1].conj(),-zj[0].conj()])
                if rephase:
                    zi, zj = zi*np.exp(1j*.37*theta), zj*np.exp(-1j*.61*theta)
                return zi, zj
            values = []
            for rephase in (False,True):
                i_minus,j_minus = registered(-eps,rephase)
                i_plus,j_plus = registered(eps,rephase)
                w_minus,w_plus = np.vdot(j_minus,i_minus),np.vdot(j_plus,i_plus)
                phase_rate = np.angle(w_minus.conjugate()*w_plus)/(2*eps)
                clock_i = np.angle(np.vdot(i_minus,i_plus))/(2*eps)
                clock_j = np.angle(np.vdot(j_minus,j_plus))/(2*eps)
                values.append(phase_rate+clock_j-clock_i)
            predicted = -mirror*np.tan(beta/2)/2
            print(f"neighbor_tilt={beta:g} mirror={mirror:+d} bond_drive_per_tilt_rate={values[0]:.12f} "
                  f"rephased={values[1]:.12f} prediction={predicted:.12f}")
            assert max(abs(np.array(values)-predicted))<5e-11


if __name__ == "__main__":
    pair_test()
    projected_clock()
