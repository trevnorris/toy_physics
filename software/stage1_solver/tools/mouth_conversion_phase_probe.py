"""Independent patch observable for a common shift of both transit levels.

The observable compares phase-referenced two-channel transmission matrices
on the same patch, before and after preparation. No circuit twist enters.
Unknown spin-dependent rotation drops out of their determinant ratio.
All inputs and readings below are synthetic, in the conversion test units.
"""
import numpy as np
from scipy.linalg import expm

from mouth_collar_order_preparation import PAULI, calibration, predicted


def exposure(axis, eta, retention, encounters):
    return (1 - retention**encounters) * (1 - eta * axis[2]) / 2


def transmission(axis, eta, retention, encounters, dwell, coupling):
    orientation = predicted(axis[None, :], eta, retention, encounters)[0]
    deposited = exposure(axis, eta, retention, encounters)
    energy = np.diag([0.0, 1.0]) - 2 * np.einsum("k,kij->ij", orientation, PAULI)
    return expm(-1j * dwell * (energy + coupling * deposited * np.eye(2)))


def phase_shift(axis, eta, retention, encounters, dwell, coupling):
    reference = transmission(axis, eta, retention, 0, dwell, coupling)
    prepared = transmission(axis, eta, retention, encounters, dwell, coupling)
    relative = reference.conj().T @ prepared
    determinant = np.linalg.det(relative)
    # Follow the zero-dwell branch. Chosen inputs keep the double phase in (-pi, pi).
    assert 2 * abs(coupling * exposure(axis, eta, retention, encounters) * dwell) < np.pi
    return np.angle(determinant) / 2


def main():
    retention = calibration()
    axis = np.array([0.8, 0.0, -0.6])
    for coupling in (0.0, 1.0, -1.0):
        calibrated_phase = phase_shift(axis, 1, retention, 6, 1.0, coupling)
        estimate = -calibrated_phase / exposure(axis, 1, retention, 6)
        print(f"training input={coupling:+g} phase={calibrated_phase:.12g} inferred={estimate:.12g}")
        for encounters, dwell, eta in ((3, 1.0, 1), (6, 0.5, 1), (6, 1.0, -1)):
            direction = eta * axis
            measured = phase_shift(direction, eta, retention, encounters, dwell, coupling)
            forecast = -estimate * exposure(direction, eta, retention, encounters) * dwell
            print(f"held_out encounters={encounters} dwell={dwell} eta={eta:+} "
                  f"observed_phase={measured:.12g} predicted_phase={forecast:.12g} "
                  f"absolute_error={abs(measured-forecast):.3g}")
            assert abs(measured - forecast) < 1e-12


if __name__ == "__main__":
    main()
