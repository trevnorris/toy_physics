"""Radiation from two prescribed cuts with a staggered two-state pair on each patch.

Spin operators are dimensionless. One pair per unit physical area, no bath,
exp(-i omega t), independent free site rolls, and bonds following the dominant
local face ray. All sources translate along the line of centers. Surface
lengths in the exterior solver are measured in units of each mouth's radius.
"""
import numpy as np
import sympy as sp

from mouth_pair_translating_cuts import AMPLITUDE, MODULUS, SPEED, Pair, hankel, polynomial
from mouth_translation_load import components


PIN, EXCHANGE, HBAR, BOND = 0.10, 0.10, 10., 8.
MOTION = np.array([AMPLITUDE, -AMPLITUDE])
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
TOP = np.array([np.kron(matrix, np.eye(2))/2 for matrix in PAULI])
BOTTOM = np.array([np.kron(np.eye(2), matrix)/2 for matrix in PAULI])
STAGGERED, TOTAL = TOP - BOTTOM, TOP + BOTTOM


def symbolic_check():
    # In the singlet and central triplet block the normalized ground-state
    # components are proportional to (E,h), E the gap to either outer triplet.
    E, pin, energy = sp.symbols("E h epsilon", positive=True)
    exchange = (E**2-pin**2)/E
    block = sp.Matrix([[-exchange, -pin], [-pin, 0]])  # measured from outer triplet
    state = sp.Matrix([E, pin])
    pol = 2*E*pin/(E**2+pin**2)
    chi = E**2/(E**2+pin**2)*(1/(E-energy)+1/(E+energy))
    residual = list(sp.simplify((block+E*sp.eye(2))*state))
    residual.append(sp.factor(chi/(pol/pin)-E**2/(E**2-energy**2)))
    print("symbolic ground-state and transverse-susceptibility residuals", residual)
    assert residual == [0, 0, 0]


def analytic_site(exchange):
    half_width = np.hypot(exchange/2, PIN)
    return PIN/half_width, exchange/2 + half_width


def response(H, drive, energy):
    """Solve the four-state linearized von Neumann equation in the product basis."""
    _, states = np.linalg.eigh(H)
    rho = np.outer(states[:, 0], states[:, 0].conj())
    identity = np.eye(4)
    commutator = np.kron(identity, H) - np.kron(H.T, identity)
    source = rho@drive - drive@rho  # perturbation of the Hamiltonian is -f drive
    drho = np.linalg.solve(energy*np.eye(16)-commutator, source.ravel(order="F"))
    return rho, drho.reshape((4, 4), order="F")


def site(exchange, eta, energy):
    H = sum(exchange*a@b for a, b in zip(TOP, BOTTOM))-eta*PIN*STAGGERED[2]
    rho, change = response(H, STAGGERED[0], energy)
    polarization = eta*np.trace(rho@STAGGERED[2]).real
    transverse = np.trace(change@STAGGERED[0])
    common = np.trace(change@TOTAL[1])
    return polarization, transverse, common


def compare_local():
    for exchange in (0., EXCHANGE, 3*EXCHANGE):
        pol, gap = analytic_site(exchange)
        for eta in (1, -1):
            energy = .62*gap
            obs, susceptibility, common = site(exchange, eta, energy)
            expected = pol/PIN/(1-(energy/gap)**2)
            # Heisenberg equation for common angular momentum about y.
            residual = -1j*energy*common - eta*(pol - PIN*susceptibility)
            print(f"J={exchange:g} face={eta:+d} polarization={obs:.9g} "
                  f"transverse_chi={susceptibility.real:.9g} common_y_imag={common.imag:.9g} "
                  f"gap={gap:.9g} torque_residual={abs(residual):.3g}")
            assert abs(obs-pol)<1e-13 and abs(susceptibility-expected)<1e-11
            assert abs(residual)<1e-12


def geometric_drive(lmax, flux, carry):
    pref = np.zeros(lmax)
    pref[0] = (1-carry)*flux/(4*np.pi)
    return pref


def coupled(pair, k, radius, carry, energy, surface):
    """Surface field plus total spin tilt U solved jointly from spin matrix response.

    U is the total angular connection shift on the tagged material sites.
    Relative angular slip of the reference is U-bX. The source driving the
    spin is the torque from the bond plus that from turning its local pin.
    """
    count = pair.lmax*len(surface)
    matrix = np.zeros((2*count, 2*count), complex)
    rhs = np.zeros(2*count, complex)
    bond = radius*BOND/MODULUS
    for i, ((theta, y, t, v), x, q) in enumerate(zip(surface, MOTION, pair.charges)):
        rows = np.arange(i*pair.lmax, (i+1)*pair.lmax)
        pol, chi, _ = site(EXCHANGE, int(np.sign(q)), energy)
        susceptibility = MODULUS/radius*chi/(2*radius*pol)**2
        b = geometric_drive(pair.lmax, q, carry)
        # External torque in units of C/R. The pin force times susceptibility
        # equals PIN*chi/pol times the physical tilt of the material normal.
        matrix[rows, :count] = theta + bond*y
        matrix[rows, rows+count] = -bond
        rhs[rows] = -x*(t+bond*(v+b))
        matrix[rows+count, :count] = -susceptibility*bond*y
        matrix[rows+count, rows+count] = 1+susceptibility*bond
        rhs[rows+count] = x*(susceptibility*bond*(v+b) + PIN*chi/pol*b)
    return np.linalg.solve(matrix, rhs)[:count]


def effective(pair, radius, carry, ratio, surface):
    """Use the analytic one-pair spectrum and spring constants, with no moving fit."""
    pol, _ = analytic_site(EXCHANGE)
    h = 4*radius**2*pol*PIN
    inverse = h*(1-ratio**2)
    g = radius/MODULUS*BOND*inverse/(BOND+inverse)
    motion_torque = -radius/MODULUS*BOND*h*ratio**2/(BOND+inverse)
    blocks, drives = [], []
    for (theta, y, t, v), x, q in zip(surface, MOTION, pair.charges):
        blocks.append(theta+g*y)
        drives.append(-x*(t+g*v+motion_torque*geometric_drive(pair.lmax, q, carry)))
    return np.linalg.solve(np.vstack(blocks), np.concatenate(drives))


def equator(pair, out, k):
    # r a_phi exp(-ikr) referenced to positive-mouth peak displacement at t=0.
    value = 0j
    for block in out.reshape(-1, pair.lmax):
        for degree, coef in enumerate(block, 1):
            value += coef*(-1j)**(degree+1)*polynomial(degree, 0.)[1]/(k*hankel(degree,k))
    return value


def pair_readings(radius, ratio, degree):
    pol, gap = analytic_site(EXCHANGE)
    k = ratio*gap/HBAR*radius/SPEED
    h = 4*radius**2*pol*PIN
    g_static = radius/MODULUS*BOND*h/(BOND+h)
    pair = Pair(np.array([-2., 2.]), degree, g_static)
    static_tilt = max(np.max(abs(pair.background(z, 0.)[1])) for z in pair.centers)*2*MODULUS/(radius*h)
    print(f"R={radius:g} L={degree} sampled_static_tilt_rad={static_tilt:.9g}")
    _, surface = components(pair, k)
    readings = []
    for carry in (0., .4, 1.):
        direct = coupled(pair, k, radius, carry, ratio*gap, surface)
        predicted = effective(pair, radius, carry, ratio, surface)
        far = equator(pair, direct, k)
        phase_check = abs(equator(pair, predicted, k)/far-1)
        coeff_error = np.linalg.norm(direct-predicted)/np.linalg.norm(direct)
        power = pair.outgoing(direct, k)
        work = sum(pair.flux(direct, z, 1.3, k) for z in pair.centers)
        work_error = abs(work/power-1)
        print(f"R={radius:g} kR={k:.8g} frequency_over_gap={ratio:g} L={degree} carry={carry:g} "
              f"far_angle_abs={abs(far):.10g} phase_rad={np.angle(far):.10g} "
              f"power_in_radius_units={power:.10g} complex_error={phase_check:.3g} "
              f"coeff_error={coeff_error:.3g} work_error={work_error:.3g}")
        assert max(phase_check, coeff_error)<1e-10 and work_error<1e-5 and power>0
        readings.append((far, power))
    interpolation = (readings[0][0]*.6+readings[-1][0]*.4)
    assert abs(interpolation/readings[1][0]-1)<1e-11
    relative = readings[0][0]/readings[-1][0]
    print(f"R={radius:g} frequency_over_gap={ratio:g} L={degree} ambient_over_carried "
          f"amplitude={abs(relative):.9g} phase_rad={np.angle(relative):.9g} "
          f"power={readings[0][1]/readings[-1][1]:.9g}")
    return np.array([reading[0] for reading in readings])


def main():
    symbolic_check()
    compare_local()
    for radius in (5., 10.):
        for ratio in (.30, .65, 1.15):
            coarse = pair_readings(radius, ratio, 6)
            fine = pair_readings(radius, ratio, 10)
            difference = max(abs(coarse/fine-1))
            print(f"R={radius:g} frequency_over_gap={ratio:g} multipole_refinement={difference:.3g}")
            assert difference<2e-4


if __name__ == "__main__":
    main()
