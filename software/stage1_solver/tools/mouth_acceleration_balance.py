"""Read the momentum transfer, force and full volume balance of the crossing pair.

Positive z momentum flux flows from the feed into the mouth. Q denotes the
field per unit sine velocity amplitude, and -Q/omega its response per
acceleration at the equal-speed crossings. All momentum uses the lab frame.
"""

import argparse
from pathlib import Path

import numpy as np

from mouth_acceleration_response import load_preparation
from mouth_carrier_momentum import integrate
from mouth_feed_weak_check import dx4


def mouth_readings(patch, base, j, h):
    """Flux derivative on the flat aperture, with the density-slope rim term zero."""
    face, change = base[-1], h[-1]
    n = abs(face)**2
    dn = 2*(face.conj()*change).real
    dx, dz = patch.dx, patch.dz
    area = dx**3
    dj = j*(dn/n - np.sum(j*dn/n)/j.sum())
    vx = (face.conj()*dx4(face,dx)).imag/n
    dvx = ((face.conj()*dx4(change,dx) + change.conj()*dx4(face,dx)).imag - dn*vx)/n
    slope = -2*patch.restoring*(1-1/abs(face))
    slope_change = -patch.restoring*dn/abs(face)**3
    carried = area*np.sum(dj*vx+j*dvx)
    stress = area*np.sum(slope_change*dx4(n,dx) + slope*dx4(dn,dx))/4

    # Reconstruct both normal derivatives ignoring the imposed mouth law.
    def from_below(f):
        return (11*f[-1]-18*f[-2]+9*f[-3]-2*f[-4])/(6*dz)
    fz = from_below(base)
    j_fd = (face.conj()*fz).imag
    dj_fd = (change.conj()*fz + face.conj()*from_below(h)).imag
    slab_n = abs(base[-4:])**2
    slab_dn = 2*(base[-4:].conj()*h[-4:]).real
    slope_fd = from_below(np.log(slab_n))
    slope_change_fd = from_below(slab_dn/slab_n)
    mask = patch.mouth
    carried_fd = area*np.sum((dj_fd*vx+j_fd*dvx)[mask])
    stress_fd = -area*np.sum((dn*dx4(slope_fd,dx)+n*dx4(slope_change_fd,dx))[mask])/4
    return carried, stress, carried_fd, stress_fd, vx, dvx, area*np.sum(dj)


def check_equations(patch, f, current, h, other, v0, omega, cosine):
    """Real-space ghost, PDE and conserved input checks on the archived harmonics."""
    n = abs(f)**2
    dn = 2*(f.conj()*h).real
    jf = current*(dn[-1]/n[-1] - np.sum(current*dn[-1]/n[-1])/current.sum())
    velocity_z = current/n[-1]
    speed_change = (jf - velocity_z*dn[-1])/n[-1]
    amp = abs(f[-1])
    normal = -patch.restoring*(1-1/amp) + 1j*velocity_z
    dnormal = -patch.restoring*dn[-1]/(2*amp**3) + 1j*speed_change
    ghost = h[-2]+2*patch.dz*(normal*h[-1]+dnormal*f[-1])
    lap = np.diff(np.concatenate([h, ghost[None]]), n=2, axis=0)/patch.dz**2
    for axis in (1,2,3):
        lap += (np.roll(h[1:],1,axis)+np.roll(h[1:],-1,axis)-2*h[1:])/patch.dx**2
    residual = (-lap/2 + (n[1:]**4-1)*h[1:]/4 + n[1:]**3*dn[1:]*f[1:]
                + 1j*v0*patch.x_derivative(h[1:]))
    if cosine:
        residual += 1j*patch.x_derivative(f[1:]) - 1j*omega*other[1:]
    else:
        residual += 1j*omega*other[1:]
    return float(np.max(abs(residual)))


def momentum_budget(patch, f, p, q, omega):
    # At sine crossing the momentum time derivative per acceleration is
    # the cosine stored-momentum amplitude. Independently derive the
    # semidiscrete conservation identity with the same antisymmetric D_x
    # as the PDE. The lattice potential chain-rule error is explicit.
    grad0 = patch.x_derivative(f)
    gradp = patch.x_derivative(p)
    stored = integrate((p.conj()*grad0+f.conj()*gradp).imag, patch.dx, patch.dz)
    stored4 = integrate((p.conj()*dx4(f.swapaxes(0,1),patch.dx).swapaxes(0,1)
                         + f.conj()*dx4(p.swapaxes(0,1),patch.dx).swapaxes(0,1)).imag,
                        patch.dx, patch.dz)
    n = abs(f)**2
    dn = 2*(f.conj()*q).real
    local = (n**4-1)/4
    base_potential = local*f
    dpotential = local*q + n**3*dn*f
    lattice_sine = 2*integrate((dpotential.conj()*grad0 +
                                base_potential.conj()*patch.x_derivative(q)).real, patch.dx, patch.dz)
    return stored, stored4, lattice_sine


def main(directory, output):
    for file in sorted(directory.glob("omega*.npz")):
        with np.load(file) as data:
            level, law = int(data["level"]), str(data["law"])
            omega, v0, amplitude = (float(data[k]) for k in ("omega","v0","amplitude"))
            p, q = data["p"], data["q"]
        patch, f, j = load_preparation(level,law)
        stored, stored4, lattice_sine = momentum_budget(patch, f, p, q, omega)
        carry, stress, carry_fd, stress_fd, _, _, net = mouth_readings(patch,f,j,q)
        # The integrated mixed flux equals Re integral f_z^* D_x f:
        # the derivative-of-density stress is integrated by parts on the
        # periodic plane, with its given slope zero outside and at the rim.
        face, change = f[-1], q[-1]
        n = abs(face)**2
        dn = 2*(face.conj()*change).real
        normal = -patch.restoring*(1-1/abs(face)) + 1j*j/n
        normal_change = (-patch.restoring*dn/(2*abs(face)**3)
                         - 1j*(j/n)*np.sum(j*dn/n)/j.sum())
        fz, qz = normal*face, normal*change+normal_change*face
        mixed2 = patch.dx**3*np.sum((qz.conj()*patch.x_derivative(face)
                                     +fz.conj()*patch.x_derivative(change)).real)
        # At the uniform feed the tangential derivatives vanish; the net
        # feed stress is the periodic integral of a tangential derivative.
        # Read it independently from the one-sided normal density variations.
        feed_dn = 2*(f[:4].conj()*q[:4]).real/abs(f[:4])**2
        feed_ds = (-11*feed_dn[0]+18*feed_dn[1]-9*feed_dn[2]+2*feed_dn[3])/(6*patch.dz)
        feed = -patch.dx**3 * dx4(feed_ds,patch.dx).sum()/4
        residuals = [check_equations(patch,f,j,h,other,v0,omega,is_p)
                     for h,other,is_p in ((p,q,True),(q,p,False))]
        # Coefficients of the brane load relative to A, so a negative
        # coefficient is a reaction opposing the prescribed acceleration.
        mouth = -(carry+stress)/omega
        mouth2 = -mixed2/omega
        correction = -lattice_sine/omega
        balance = stored + mouth2 - correction  # feed is zero in the discrete identity
        scale = max(abs(stored),abs(mouth2),1e-12)
        acceleration = amplitude*omega
        if np.max(abs(p[0]))+np.max(abs(q[0])) != 0:
            raise ValueError("Lab feed phase is perturbed")
        print(f"n={level} depth={f.shape[0]-1} law={law} omega={omega:g} "
              f"V_cross={v0:g} A_size={acceleration:.8g} "
              f"carried_per_A={-carry/omega:.12g} stress_per_A={-stress/omega:.12g} "
              f"mouth_per_A={mouth:.12g} mouth_direct_fd_per_A={-(carry_fd+stress_fd)/omega:.12g} "
              f"stored_rate_per_A={stored4:.12g} feed_per_A={-feed/omega:.3g} "
              f"brane_drive_half_difference={-acceleration*mouth:.12g} "
              f"mixed2_per_A={mouth2:.12g} stored2_per_A={stored:.12g} "
              f"lattice_body_per_A={correction:.12g} discrete_balance_residual={balance:.3g} "
              f"relative_balance_residual={balance/scale:.3g} "
              f"max_harmonic_residual={max(residuals):.3g} net_rate_variation={net:.3g}", flush=True)
        if (not np.isfinite([mouth,stored4,balance,*residuals]).all() or abs(net)>1e-11
                or max(residuals)>2e-8 or abs(feed)>1e-12 or abs(balance)>3e-7*scale):
            raise ArithmeticError("Harmonic equations, feed or integrated momentum balance failed")
        output.append((level,law,omega))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories",type=Path,nargs="+")
    args=parser.parse_args()
    cases=[]
    for directory in args.directories:
        main(directory,cases)
    if len(cases)!=len(set(cases)) or len(cases)!=2*len(args.directories):
        raise ValueError("Incomplete or duplicate frequency pair")
