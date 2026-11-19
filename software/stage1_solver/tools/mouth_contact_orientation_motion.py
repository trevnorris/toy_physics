"""Dipolar tilt of two orientable faces crossing a radial ordering texture.

The tested site dynamics add equal first-order Berry actions to the prior
static pin and overlap alignment. Each face can cant and roll independently.
An optional antialignment energy lets one vary the even cant stiffness.
All surface densities below refer to one paired site per unit area. Lengths
are physical (the angular connection is phase per length). Time uses exp(-iwt).
"""

import numpy as np


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def ray(direction):
    """Use regular north/south sections of the normalized two-component ray."""
    x, y, z = np.moveaxis(direction, -1, 0)
    north = np.stack((1 + z, x + 1j*y), axis=-1)
    south = np.stack((x - 1j*y, 1 - z), axis=-1)
    return unit(np.where((z >= 0)[..., None], north, south))


def overlap(a, b):
    return np.sum(a.conj()*b, axis=-1)


def loop(states):
    phase = overlap(states, np.roll(states, -1, axis=0))
    return np.prod(phase / abs(phase))


def moving_latitude(radius, theta, eta, carry, count, dt=0.0001):
    """Gauge-invariant time/azimuth quadrilaterals on tagged, translating sites."""
    phi = 2*np.pi*np.arange(count)/count
    sphere = radius * np.column_stack((np.sin(theta)*np.cos(phi),
                                       np.sin(theta)*np.sin(phi), np.full(count, np.cos(theta))))
    old, new = [ray(eta*unit(sphere + np.array([0, 0, (carry-1)*time])))
                * np.exp(1j*(0.7*np.sin(3*phi) + time*np.cos(2*phi)))[:, None]
                for time in (-dt/2, dt/2)]
    mixed = (overlap(old, new)*overlap(new, np.roll(new, -1, axis=0))
             * overlap(np.roll(old, -1, axis=0), np.roll(new, -1, axis=0)).conj()
             * overlap(old, np.roll(old, -1, axis=0)).conj())
    length = 2*np.pi*radius*np.sin(theta)
    measured = np.angle(mixed).sum()/(dt*length)
    fidelity_rate = np.mean((1-abs(overlap(old, new))**2)/dt**2)
    return measured, fidelity_rate


def site_form(eta, berry, step):
    """Obtain the four real spin-coordinate symplectic coefficients from loops."""
    def pair(q):
        # At longitude zero the equator has e_r=x, e_theta=-z, e_phi=y.
        return ray(unit(np.array([[eta, q[1], -q[0]], [-eta, q[3], -q[2]]])))

    form = np.zeros((4, 4))
    for i, ei in enumerate(np.eye(4)*step):
        for j, ej in enumerate(np.eye(4)*step):
            if j > i:
                states = np.stack([pair(q/2) for q in (-ei-ej, ei-ej, ei+ej, -ei+ej)])
                form[i, j] = berry*np.angle(loop(states[:, 0])*loop(states[:, 1]))/step**2
    return form-form.T


class Ring:
    """Evaluate energy per paired site on an equatorial circle of an l=1 tilt.

    Meridional tilts give the azimuthal coexact link. An azimuthal tilt's
    first-order meridional link is an exact gradient for this axisymmetric
    dipole and both faces can roll to relax it. Integrated azimuthal loops
    discard exact gradients, including the choice of local section.
    """

    def __init__(self, radius, eta, carry, count, pin, exchange, bond):
        self.radius, self.eta, self.carry = radius, eta, carry
        self.pin, self.exchange, self.bond = pin, exchange, bond
        phi = 2*np.pi*np.arange(count)/count
        self.radial = np.column_stack((np.cos(phi), np.sin(phi), np.zeros(count)))
        self.azimuth = np.column_stack((-np.sin(phi), np.cos(phi), np.zeros(count)))
        self.meridian = np.array([0., 0., -1.])

    def state(self, q):
        tp, up, tm, um, translation, _ = q
        a = unit(self.eta*self.radial + tp*self.meridian + up*self.azimuth)
        b = unit(-self.eta*self.radial + tm*self.meridian + um*self.azimuth)
        center_offset = (1-self.carry)*translation/self.radius
        native = self.eta*unit(self.radial + center_offset*self.meridian)
        phases = [loop(ray(p))*loop(ray(n)).conj() for p, n in ((a, native), (b, -native))]
        shifts = np.angle(phases)/(2*np.pi*self.radius)
        return a[0], b[0], native[0], shifts

    def energy(self, q):
        a, b, native, shifts = self.state(q)
        pinning = self.pin*(sum((a-native)**2) + sum((b+native)**2))/2
        canting = self.exchange*sum((a+b)**2)/2
        slip = np.array([q[-1], -q[-1]]) - shifts
        return pinning + canting + self.bond*(slip@slip)/4


def response(radius, eta, carry, count, pin, exchange, bond, berry, omega):
    ring = Ring(radius, eta, carry, count, pin, exchange, bond)
    step = 0.0002
    basis = step*np.eye(6)
    # Hessian of the actual normalized-ray spring energies, including both drives.
    hessian = np.array([[(ring.energy(i+j) - ring.energy(i-j)
                          - ring.energy(-i+j) + ring.energy(-i-j))/(4*step**2)
                         for j in basis] for i in basis])
    form = site_form(eta, berry, step)
    q = -np.linalg.solve(hessian[:4, :4]+1j*omega*form, hessian[:4, 4:])
    # Columns: imposed translation X at y=0; imposed relative angular slip y at X=0.
    full_q = np.vstack((q, np.eye(2)))
    shift_gradient = np.stack([(ring.state(e)[3]-ring.state(-e)[3])/(2*step)
                               for e in basis], axis=1)
    measured = (np.array([1., -1.]) @ shift_gradient @ full_q)/2
    cant = (q[1]+q[3])/2
    # Independently eliminate the common cant in the continuum dipole action.
    # Paired-site inertia for the actual odd direction is berry^2/(2(pin+2exchange)).
    inertia = 2*radius**2*berry**2/(pin+2*exchange)
    stiffness = 8*radius**2*pin
    geometric = eta*(1-carry)/(2*radius**2)
    predicted = np.array([inertia*omega**2*geometric, bond])/(stiffness+bond-inertia*omega**2)
    error = np.max(abs(measured-predicted))
    projection = np.vstack((np.eye(2), -np.eye(2)))
    conjugate_form = np.max(abs(projection.T@form@projection))
    print(f"R={radius:g} eta={eta:+d} site_carry={carry:g} J={exchange:g} omega={omega:g} "
          f"numerical_s_per_X={measured[0].real:.9f} reduced={predicted[0]:.9f} "
          f"numerical_s_per_y={measured[1].real:.9f} reduced={predicted[1]:.9f} "
          f"even_azimuthal_cant_per_X_imag={cant[0].imag:.9f} "
          f"max_error={error:.3g} locked_Berry_form={conjugate_form:.3g}")
    return error, conjugate_form


def main():
    radius, theta = 1.7, 0.91
    for eta, carry in ((1, 0.), (-1, 0.), (1, 0.4), (1, 1.)):
        geometry = eta*(1-carry)*np.sin(theta)/(2*radius**2)
        metric = ((1-carry)*np.sin(theta)/(2*radius))**2
        coarse = moving_latitude(radius, theta, eta, carry, 128)
        fine = moving_latitude(radius, theta, eta, carry, 512)
        print(f"eta={eta:+d} site_carry={carry:g} "
              f"mixed_time_phi={fine[0]:.9f} curvature_transport={geometry:.9f} "
              f"fidelity_metric={fine[1]:.9f} predicted={metric:.9f} "
              f"loop_error_N128={abs(coarse[0]-geometry):.3g} "
              f"loop_error_N512={abs(fine[0]-geometry):.3g}")
        assert abs(fine[0]-geometry) < 1e-5 and abs(fine[1]-metric) < 1e-6

    for exchange, omega in ((0., 0.7), (0.8, 0.7), (0.8, 1.9)):
        for eta, carry in ((1, 0.), (-1, 0.), (1, 0.4), (1, 1.)):
            err, locked = response(radius, eta, carry, 512, 1.2, exchange, 5., 1.6, omega)
            assert err < 1e-4 and locked < 1e-8


if __name__ == "__main__":
    main()
