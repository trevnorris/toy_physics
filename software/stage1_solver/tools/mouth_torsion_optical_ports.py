"""Reciprocal mechanical probes embedded in two infinite transverse leads.

Two uniform instrumented planes make this a normal-incidence problem in x.
q_top and q_bottom are individual contact piston readouts. Operator and
loads are constructed in the physical face basis before any projection.
"""
from math import sin, sqrt

import numpy as np

K, EDGE, CORNER = 8., 2., 1.
MU, RHO = 2., 4.
ARM, DIAGONAL, H, BRIDGE, MASS = .2, .75, 1.2, .4, 12.
DRAG = .6
EVEN = np.array([1., 1.]) / sqrt(2)
ODD = np.array([1., -1.]) / sqrt(2)


def lead(omega, angular):
    s = omega**2
    if angular:
        r2 = -CORNER*s/(2*K-CORNER*s)
        r4 = -CORNER*s/(4*K-CORNER*s)
        c = K*(1-2*r2+r4)/2
        z = 8*K*r2-2*EDGE*s
    else:
        c, z = MU, -RHO*s
    # The two eliminated semi-infinite leads carry outgoing positive flux.
    cos_k = 1 + z / (2*c)
    if not -1 < cos_k < 1 or c <= 0:
        raise ValueError("Frequency outside the tested acoustic passband")
    return c, z, np.arccos(cos_k)


def segment(omega, separation, drive, read, source_at_right=False):
    sites = separation + 7
    matrix = np.zeros((2*sites+4, 2*sites+4), dtype=complex)
    cuts = (3, 3+separation)
    channels = []
    for number in range(2):
        c, z, k = lead(omega, angular=bool(number))
        sl = slice(number*sites, (number+1)*sites)
        idx = np.arange(sites-1)
        local = np.eye(sites, dtype=complex) * z
        local[idx, idx] += c
        local[idx+1, idx+1] += c
        local[idx, idx+1] -= c
        local[idx+1, idx] -= c
        local[[0,-1],[0,-1]] += c*(1-np.exp(1j*k))
        matrix[sl, sl] = local
        channels.append((sl, c, k))
    probes = []
    for number, cut in enumerate(cuts):
        pair = np.arange(2*sites+2*number, 2*sites+2*number+2)
        probes.append(pair)
        for piston, sign in zip(pair, (1, -1)):
            matrix[piston,piston] += H+4*K*ARM**2 - MASS*omega**2
            # The ordinary diagonal sensor and the four signed spring torques.
            for offset, coupling in ((0,-H*DIAGONAL), (sites,sign*K*ARM)):
                nodes = offset + np.array([cut,cut+1])
                gradient = np.array([-1.,1.])
                matrix[piston,nodes] += coupling * gradient
                matrix[nodes,piston] += coupling * gradient
                if offset == 0:
                    matrix[np.ix_(nodes,nodes)] += H*DIAGONAL**2*np.outer(gradient,gradient)
        matrix[np.ix_(pair,pair)] += BRIDGE*np.array([[1.,-1.],[-1.,1.]])
    src, det = probes[::-1] if source_at_right else probes
    # Apply the same passive two-channel mechanical termination to both
    # probes. 'read' picks which detector dial's dissipated power is quoted.
    for pair in probes:
        matrix[pair,pair] -= 1j*omega*DRAG
    rhs = np.zeros(len(matrix), complex)
    rhs[src] = drive
    q = np.linalg.solve(matrix,rhs)
    incoming = omega*(rhs.conj()@q).imag/2
    absorption = [DRAG*omega**2*np.vdot(q[p],q[p]).real/2 for p in probes]
    radiation = [omega*c*sin(k)*sum(abs(q[sl][[0,-1]])**2)/2 for sl,c,k in channels]
    measured = DRAG*omega**2*abs(read @ q[det])**2/2
    error = abs(incoming-sum(absorption)-sum(radiation))/incoming
    return read @ q[det], measured / incoming, np.array(radiation)/incoming, error


def static_loading():
    # Independent energy Hessians of the square after relaxing its corners.
    # Coordinates are frame strain e, primitive circulation F and both pistons.
    hess = np.diag([MU, 0., 0., 0.])
    for piston, sign in ((2,1),(3,-1)):
        spring = np.array([0., sign, 0., 0.])
        spring[piston] = 4*ARM
        strap = np.array([-DIAGONAL, 0., 0., 0.])
        strap[piston] = 1
        hess += K/4*np.outer(spring,spring) + H*np.outer(strap,strap)
    relative = np.array([0.,0.,1.,-1.])
    hess += BRIDGE*np.outer(relative,relative)
    rest = hess[:2,:2] - hess[:2,2:]@np.linalg.solve(hess[2:,2:],hess[2:,:2])
    expected = [MU+8*H*K*ARM**2*DIAGONAL**2/(H+4*K*ARM**2),
                (K/2)*(H+2*BRIDGE)/(H+2*BRIDGE+4*K*ARM**2)]
    print(f"decorated strain/curl moduli={rest.diagonal().tolist()} mixed={rest[0,1]:.3g} "
          f"bare={[MU,K/2]} missing_bare_shear_residual={rest[0,0]-MU:.12f}")
    assert np.allclose(rest,np.diag(expected))


def main():
    static_loading()
    for frequency in (.23, .38, .54):
        for source, drive in (("common",EVEN),("opposite",ODD),("one_face",np.array([1.,0.]))):
            result = [segment(frequency,12,drive,port) for port in (EVEN,ODD)]
            detected = [r[1] for r in result]
            radiated = result[0][2]
            balance = max(r[3] for r in result)
            swaps = [segment(frequency,12,port,drive,True)[0] for port in (EVEN,ODD)]
            reciprocity = max(abs(r[0]-sw) for r,sw in zip(result,swaps))
            print(f"omega={frequency:g} source={source} detector_common={detected[0]:.12g} "
                  f"detector_opposite={detected[1]:.12g} lead_ordinary={radiated[0]:.9g} "
                  f"lead_torsion={radiated[1]:.9g} balance={balance:.3g} swap_error={reciprocity:.3g}")
            assert balance < 1e-11 and reciprocity < 1e-11
            if source != "one_face":
                bright = int(source == "opposite")
                assert detected[bright] > 1e-8 and detected[1-bright] < 1e-24
                assert radiated[bright] > 0 and radiated[1-bright] < 1e-24
            else:
                assert min(detected) > 1e-8

if __name__ == "__main__":
    main()
