"""Apply a weak harmonic rim load in a conditional brane-only force balance.

Linearize about the steadily translating, inward-draining mouth. The steady
drive is held fixed. The test sets the *incremental* total external force
plus the integrated carrier load to zero. A physical free-mouth response
would require an independent check of this closure.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_acceleration_balance import mouth_readings
from mouth_acceleration_response import load_preparation
from mouth_release_input_check import manifest_fields
from stage1_solver.mouth_probe import spherical_rim_probe


def run(level, folder, hashes):
    law = "restoring_slope"
    omega = 0.1
    radius = 1.0
    patch, base, current = load_preparation(level, law)
    rate = patch.dx**3 * current.sum()
    if not np.isclose(rate, 0.25, atol=1e-12, rtol=0):
        raise ValueError("Mouth preparation differs from the reference intake")

    path = folder / "fields/omega0.1.npz"
    if sha256(path.read_bytes()).hexdigest() != hashes[path]:
        raise ValueError(f"Unverified harmonic response: {path}")
    with np.load(path, allow_pickle=False) as packet:
        if (int(packet["level"]), str(packet["law"]), float(packet["omega"]),
                float(packet["v0"])) != (level, law, omega, 0.01):
            raise ValueError("Harmonic field and base preparation do not match")
        reference_amplitude = float(packet["amplitude"])
        load = []
        direct = []
        for key in ("p", "q"):
            # Each field is per unit velocity amplitude. The second readout
            # estimates the normal traces from the interior field planes.
            carried, stress, carried_fd, stress_fd, _, _, net = mouth_readings(
                patch, base, current, packet[key])
            if not np.isfinite([carried, stress, carried_fd, stress_fd, net]).all() or abs(net) > 1e-11:
                raise ArithmeticError("Invalid harmonic load or varying total intake")
            load.append(carried + stress)
            direct.append(carried_fd + stress_fd)

    # Re[V exp(i omega t)] gives cosine coefficient Re V and sine
    # coefficient -Im V. Z converts velocity to compensating force.
    impedance = -complex(load[0], -load[1])
    direct_impedance = -complex(direct[0], -direct[1])
    if abs(impedance) == 0:
        raise ArithmeticError("The brane-only force equation cannot be inverted")
    print(f"grid={level} radius={radius:g} inward_rate={rate:.12g} V0=0.01 omega={omega:g} "
          f"Z_real={impedance.real:.12g} Z_imag={impedance.imag:.12g}")

    mu, weights = leggauss(32)
    for force in (0.0, 2.5e-8, 5e-8):
        traction = spherical_rim_probe(radius, mu, force)
        applied = 2*np.pi*radius**2 * np.dot(weights, mu*traction)
        if not np.isclose(applied, force, atol=0, rtol=2e-14):
            raise ArithmeticError("Rim traction integral does not reproduce the applied force")
        velocity = applied / impedance
        acceleration = 1j*omega*velocity
        if abs(velocity) > reference_amplitude:
            raise ArithmeticError("Response exceeds the archived small-motion amplitude")
        direct_force = direct_impedance * velocity
        discrepancy = abs(direct_force - applied)
        print(f"  F_amp={applied:.12g} V_amp={abs(velocity):.12g} A_amp={abs(acceleration):.12g} "
              f"A_cos={acceleration.real:.12g} A_sin={-acceleration.imag:.12g} "
              f"direct_force_cos={direct_force.real:.12g} direct_force_sin={-direct_force.imag:.12g} "
              f"max_direct_force_mismatch={discrepancy:.12g}")


if __name__ == "__main__":
    hashes = manifest_fields(Path("docs/_measurements/mouth_acceleration_response.json"))
    for level, resolution in ((24, "coarse"), (32, "fine")):
        run(level, Path(f"_scratch/s11c/mouth_accel_{resolution}_restoring01"), hashes)
