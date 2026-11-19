"""Interferometric phase of the existing resonant order-preparation contact.

One path interacts with an isolated site, the reference shares its free
transit and bypasses the interaction. An independent incoming parcel can
interrogate a receiver kept after an earlier encounter. Previous parcels
are traced out. The retained receiver is a diagnostic change of preparation.
"""
import numpy as np

from mouth_collar_order_preparation import COUPLING, DURATION, INPUT, VACUUM, event, full_block


def interference_matrix(propagator, receiver):
    # K_ab = sum_{ij} <a,i|U|b,j> rho_ji. Tracing an unobserved
    # receiver on the interacting arm gives an amplitude, not a spin channel.
    return np.einsum("aibj,ji->ab", propagator.reshape(2, 4, 2, 4), receiver)


def run(axis, eta):
    local_axis = eta * axis
    preparation = full_block(local_axis, DURATION)
    spin = INPUT.copy()
    for _ in range(6):
        spin, last_receiver = event(spin, preparation)
    coherent = np.array([1, 1j, -0.3, 0.5 + 0.2j])
    coherent /= np.linalg.norm(coherent)
    states = {"fresh": VACUUM, "sixth_out": last_receiver, "coherent": np.outer(coherent, coherent.conj())}
    for dwell in (DURATION / 2, DURATION, 2 * DURATION):
        assert np.sqrt(2) * COUPLING * dwell < np.pi / 2
        propagator = full_block(local_axis, dwell)
        reference = interference_matrix(propagator, VACUUM)
        for name, receiver in states.items():
            deposited = np.trace(receiver[1:, 1:]).real
            for shift in (0.0, 1.0):
                # The added scalar shift exists only on the interacting arm.
                amplitude = interference_matrix(np.exp(-1j * shift * deposited * dwell) * propagator, receiver)
                relative_det = np.linalg.det(reference.conj().T @ amplitude)
                phase = np.angle(relative_det) / 2
                predicted = -shift * deposited * dwell
                contrast = np.linalg.svd(amplitude, compute_uv=False)[-1]
                print(f"native eta={eta:+} receiver={name} dwell={dwell:g} shift={shift:g} "
                      f"excitation={deposited:.12g} min_amplitude={contrast:.12g} "
                      f"det_modulus={abs(relative_det):.12g} common_phase={phase:.12g} "
                      f"prediction_error={abs(phase-predicted):.3g}")
                assert abs(predicted) < np.pi / 2
                assert abs(phase - predicted) < 1e-12 and contrast > 0


if __name__ == "__main__":
    axis = np.array([0.8, 0.0, -0.6])
    for orientation in (1, -1):
        run(orientation * axis, orientation)
