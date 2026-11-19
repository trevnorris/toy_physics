"""Does a retained excitation give a common conservative energy during transit?

The rim plus semi-infinite bulk supports one discrete level. Two off-band
parcel choices let both the initially occupied and empty parcel states return
after a smooth passage. Track the discrete energies and test the prediction
by evolving the original open path at finite passage duration.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from mouth_rim_delayed_transit import BAND_CENTER, Patch


RIM = 5.5
BOUND_ENERGY = 37 / 6
Z = 5 / 9


def projected_level(parcel_energy, coupling):
    def secular(energy):
        offset = energy - BAND_CENTER
        surface = 2 / (offset + np.sqrt(offset**2 - 4))
        return energy - RIM - surface - coupling**2 / (energy - parcel_energy)

    right = np.nextafter(parcel_energy, -np.inf) if parcel_energy > BOUND_ENERGY else 9.
    return brentq(secular, 6., right, xtol=1.e-14)


def follow(parcel_energy, strength, duration, sites):
    patch = Patch(sites, RIM, parcel_energy)
    states = np.column_stack((patch.parcel, patch.bound))
    centered = patch.diagonal - BAND_CENTER
    def derivative(t, flat):
        v = flat.reshape(states.shape)
        h_v = centered[:, None] * v
        g = strength * np.sin(np.pi * t / duration)**2
        h_v[2:] += v[1:-1]
        h_v[1:-1] += v[2:]
        h_v[0] += g * v[1]
        h_v[1] += g * v[0]
        return (-1j * h_v).ravel()

    result = solve_ivp(derivative, (0., duration), states.astype(complex).ravel(),
                       method="DOP853", rtol=4.e-11, atol=2.e-12)
    assert result.success
    final = result.y[:, -1].reshape(states.shape)
    free_phase = np.exp(1j * (np.array([parcel_energy, BOUND_ENERGY])-BAND_CENTER)*duration)
    overlap = free_phase[:, None] * (states.T @ final)
    a = overlap[0, 0]
    ground = overlap[1, 1]
    excited = np.linalg.det(overlap) / a
    average_level_shift = quad(lambda s: projected_level(parcel_energy, strength*np.sin(np.pi*s)**2)
                               - BOUND_ENERGY, 0, 1, epsabs=2.e-13)[0]
    prediction = -duration * average_level_shift
    # Compare phases modulo a whole turn. The chosen sweeps stay on the initial branch.
    errors = np.angle(np.array([ground, excited]) * np.exp(-1j*prediction))
    bulk_leakage = np.maximum(0., 1 - np.sum(abs(overlap)**2, axis=0))
    energy_in = np.array([parcel_energy, BOUND_ENERGY])
    energy_out = np.array([patch.energy(v) for v in final.T])
    print(f"pass parcel={parcel_energy:g} dwell={duration:g} sites={sites} "
          f"ground_phase={np.angle(ground):.12g} excited_phase={np.angle(excited):.12g} "
          f"adiabatic_prediction={prediction:.12g} max_phase_error={max(abs(errors)):.8g} "
          f"bound_leakage={bulk_leakage[1]:.8g} parcel_leakage={bulk_leakage[0]:.8g} "
          f"bound_energy_change={energy_out[1]-energy_in[1]:.8g} "
          f"parcel_energy_change={energy_out[0]-energy_in[0]:.8g}")
    return np.r_[np.angle(ground),np.angle(excited),errors,bulk_leakage,energy_out-energy_in]


def main():
    strength = .3
    for e in (1.5, 7.5):
        for g in (.3, .15, .075, .0375):
            displacement = projected_level(e, g) - BOUND_ENERGY
            print(f"energy parcel={e:g} G={g:g} shift={displacement:.12g} "
                  f"shift_over_G2={displacement/g**2:.12g} weak_limit={Z/(BOUND_ENERGY-e):.12g}")
        series = [follow(e, strength, duration, 256) for duration in (16, 48, 128)]
        repeat = follow(e, strength, 128, 352)
        size_error = np.max(abs(repeat-series[-1]))
        print(f"check parcel={e:g} max_size_error={size_error:.3g} "
              f"level_approx_error_per_time={max(abs(series[-1][2:4]))/128:.4g}")
        assert size_error < 1e-7
        assert max(series[-1][4:6]) < max(series[0][4:6])
        assert max(abs(series[-1][2:4]))/128 < max(abs(series[0][2:4]))/16


if __name__ == "__main__":
    main()
