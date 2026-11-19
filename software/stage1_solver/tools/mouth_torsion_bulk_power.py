"""Power through a proposed normal torsion contact, per face and polarization.

The earlier cubic springs act along the brane. Here a separate routing
hypothesis couples adjacent layers along w. In unit pitch the twist g = du/dw
loads four identical series springs. In a parallel fraction chi of complete
units the spring references are held. Each reference in the remaining units
has inertia and can slide, optionally against a viscous housing. The normal
angular inertia density is independently specified. All housings are held.
The linear calculation is for a welded sharp interface at normal incidence.
"""

import cmath
import math

import numpy as np


# Dimensionless one-face parameters. C, H and D below refer to total slip,
# equal to four times the slip at a single reference in symmetric loading.
K = 8.0
C = K / 4
RHO = 4.0
H = 8192.0
W0 = math.sqrt(C / H)
D_LOSSY = 0.25 * math.sqrt(C * H)
Z_INC = math.sqrt(C * RHO)


def free_contacts(twist, omega, drag):
    """Solve for four openings, four reference slips and the common torque."""
    matrix = np.zeros((9, 9), dtype=complex)
    forcing = np.zeros(9, dtype=complex)
    for j in range(4):
        # Identical series torque K(opening - slip) for each spring.
        matrix[j, j] = K
        matrix[j, j + 4] = -K
        matrix[j, 8] = -1
        # Inertial/viscous torque on one reference equals the spring torque.
        matrix[4 + j, 4 + j] = -4 * (H * omega**2 + 1j * drag * omega)
        matrix[4 + j, 8] = -1
    matrix[8, :4] = 1
    forcing[8] = twist
    if omega == 0:
        # Several static partitions of the opening leave all springs unloaded.
        # The minimum-norm representative has the permutation symmetry.
        result = np.linalg.lstsq(matrix, forcing, rcond=None)[0]
    else:
        result = np.linalg.solve(matrix, forcing)
    return result[:4], result[4:8], result[8]


def material_modulus(chi, omega, drag):
    dynamic = -H * omega**2 - 1j * drag * omega
    return chi * C + (1 - chi) * C * dynamic / (C + dynamic)


def wave_number(modulus, omega):
    k = omega * cmath.sqrt(RHO / modulus)
    if k.imag < 0 or (k.imag == 0 and k.real < 0):
        k = -k
    return k


def power(torque, velocity):
    """Positive direction is from the ordered incident lead into the receiver."""
    return -(torque * velocity.conjugate()).real / 2


def scatter(label, chi, ratio, drag, *, connected=True):
    omega = ratio * W0
    modulus = material_modulus(chi, omega, drag)
    k = wave_number(modulus, omega)
    impedance = modulus * k / omega
    if connected:
        # Unit incident angular displacement and no incoming far-bulk wave.
        reflection = (Z_INC - impedance) / (Z_INC + impedance)
        transmission = 1 + reflection
    else:
        # Open/clutch-disengaged end of lead; undriven receiver.
        reflection, transmission = 1.0, 0.0

    twist = 1j * k * transmission
    velocity = -1j * omega * transmission
    opening, slip, free_torque = free_contacts(twist, omega, drag)
    torque = chi * C * twist + (1 - chi) * free_torque
    incident = omega**2 * Z_INC / 2
    incoming = power(torque, velocity) / incident
    reflected = abs(reflection)**2

    # Positive instantaneous energy is kinetic + spring energy. Mechanical
    # loss at each reference is (4D) slip_velocity**2. For exp(i k w) the
    # cycle average decays as exp(-2 Im(k) w); integrate over the half-space.
    free_energy = (K * sum(abs(opening - slip)**2)
                   + 4 * H * omega**2 * sum(abs(slip)**2)) / 4
    stored = (RHO * abs(velocity)**2 + chi * C * abs(twist)**2) / 4
    stored += (1 - chi) * free_energy
    heat_density = (1 - chi) * (4 * drag) * omega**2 * sum(abs(slip)**2) / 2
    absorbed = heat_density / (2 * k.imag * incident) if k.imag > 0 else 0.0
    far = 0.0 if k.imag > 0 else incoming

    # Check tractions separately from power, including the open control.
    lead_torque = 1j * omega * Z_INC * (1 - reflection)
    traction_error = abs(lead_torque - torque)
    print(f"{label} chi={chi:g} omega/omega0={ratio:g} D={drag:g} "
          f"k=({k.real:.9f},{k.imag:.9f}) abs_u0={abs(transmission):.9f} "
          f"R={reflected:.9f} entry={incoming:.9f} far={far:.9f} "
          f"heat={absorbed:.9f} E0_over_Pinc={stored/incident:.9f} "
          f"traction_res={traction_error:.2g} "
          f"entry_balance={1-reflected-incoming:.2g} "
          f"bulk_balance={incoming-far-absorbed:.2g}")
    assert traction_error < 1e-12
    assert abs(torque - modulus * twist) < 1e-12
    assert abs(1 - reflected - incoming) < 1e-12
    assert abs(incoming - far - absorbed) < 1e-12
    assert incoming >= -1e-12 and stored >= 0 and absorbed >= 0

    if connected and drag == 0 and k.imag == 0:
        step = omega * 1e-5
        k_plus, k_minus = (wave_number(material_modulus(chi, w, drag), w)
                           for w in (omega + step, omega - step))
        group_velocity = 2 * step / (k_plus - k_minus).real
        energy_velocity = incoming * incident / stored
        print(f"{label} work_per_energy={energy_velocity:.12f} "
              f"dispersion_slope={group_velocity:.12f} "
              f"difference={energy_velocity-group_velocity:.2g}")
        assert abs(group_velocity - energy_velocity) < 1e-8


def main():
    print(f"one_face per_contact_K={K:g} C={C:g} rho_w={RHO:g} "
          f"per_reference_inertia={4*H:g} omega0={W0:g} "
          f"per_reference_drag_lossy={4*D_LOSSY:g} pitch=1")
    for chi in (0.0, 0.25, 0.5, 1.0):
        opening, slip, torque = free_contacts(1.0, 0.0, D_LOSSY)
        energy = chi * C / 2 + (1 - chi) * K * sum(abs(opening - slip)**2) / 2
        stress = chi * C + (1 - chi) * torque
        print(f"static chi={chi:g} unit_twist_torque={stress.real:.9f} "
              f"twice_contact_energy={2*energy:.9f} chi_C={chi*C:.9f}")
        assert abs(2 * energy - stress) < 1e-12 and abs(stress - chi * C) < 1e-12

    scatter("held", 1.0, 0.5, 0.0)
    for ratio in (0.5, 2.0):
        scatter("free_lossless", 0.0, ratio, 0.0)
        scatter("free_damped", 0.0, ratio, D_LOSSY)
    scatter("open", 0.0, 2.0, 0.0, connected=False)


if __name__ == "__main__":
    main()
