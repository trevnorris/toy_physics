"""Normal pressure pulse in the neutral two-sided coherent carrier, per unit area.

The zero-strain, unoccupied neutral order energy is additive in the native
carrier Hamiltonian.  The equilibrium half-order surface is held at s=0 for
this fixed-plane diagnostic.  No order transport kinetics are invented here.
The even-image control is a separate, reflecting half-space with s<0 at rest.
"""

import json

import numpy as np
from scipy.integrate import simpson


def sample_peak(values):
    return float(values[np.argmax(np.abs(values))])


def run(cells, dt, amplitude, reflecting):
    length = 384.0
    dx = length / cells
    s = np.arange(-cells // 2, cells // 2) * dx
    k = 2 * np.pi * np.fft.fftfreq(cells, d=dx)
    propagation = np.exp(-0.5j * dt * k**2)
    bulk_bump = np.exp(-0.5 * ((s - 48) / 4)**2)
    if reflecting:
        bulk_bump += np.exp(-0.5 * ((s + 48) / 4)**2)
    # Inverse of p(n)-p(1)=(n^5-1)/5. The initially resting bulk
    # splits this perturbation into inward and outward travelling pulses.
    psi = (1 + 5 * amplitude * bulk_bump).astype(complex) ** 0.1
    face = 6.0
    control = (s >= (0.0 if reflecting else -face)) & (s <= face)
    probes = [np.flatnonzero(s == position)[0] for position in (-face, 0, face)]
    cadence = round(0.2 / dt)
    samples = []
    energy = {}
    for step in range(round(96 / dt) + 1):
        if step % cadence == 0:
            time = step * dt
            spectrum = np.fft.fft(psi)
            psi_s = np.fft.ifft(1j * k * spectrum)
            current = (psi.conjugate() * psi_s).imag
            density = np.abs(psi)**2
            delta = density - 1
            density_ss = np.fft.ifft(-k**2 * np.fft.fft(delta)).real
            # Convective + thermodynamic + quantum normal momentum flux.
            pressure = np.expm1(5 * np.log1p(delta)) / 5
            stress = np.abs(psi_s)**2 + pressure - density_ss / 4
            left, middle, right = probes
            momentum = simpson(current[control], x=s[control])
            samples.append((time, 0 if reflecting else stress[left], stress[right],
                            momentum, current[middle], stress[middle]))
            if step in (0, round(24 / dt), round(96 / dt)):
                local_energy = (0.5 * np.abs(psi_s)**2
                                + delta**2 * (10 + delta * (10 + delta * (5 + delta))) / 20)
                energy[f"total_{time:g}"] = dx * np.sum(local_energy)
                if step == 0:
                    windows = {}
                elif time < 48:
                    windows = {"incident": (0, 48)}
                else:
                    windows = {"transmitted": (-88, -10), "reflected": (10, 88)}
                for label, (lo, hi) in windows.items():
                    window = (lo <= s) & (s <= hi)
                    energy[label] = (0.0 if reflecting and label == "transmitted"
                                     else simpson(local_energy[window], x=s[window]))
        if step == round(96 / dt):
            break
        # Strang evolution of i psi_t = -psi_ss/2 + ((|psi|^8 - 1)/4) psi.
        psi *= np.exp(-0.125j * dt * (np.abs(psi)**8 - 1))
        psi = np.fft.ifft(propagation * np.fft.fft(psi))
        psi *= np.exp(-0.125j * dt * (np.abs(psi)**8 - 1))

    t, left, right, momentum, current, center_stress = np.array(samples).T
    interval = slice(2, -2)
    momentum_rate = (momentum[:-4] - 8 * momentum[1:-3] + 8 * momentum[3:-1]
                     - momentum[4:]) / (12 * cadence * dt)
    transferred = (left - right)[interval] - momentum_rate
    independent_wall = -center_stress[interval] if reflecting else np.zeros_like(transferred)
    error_scale = max(np.max(np.abs(right)), np.finfo(float).tiny)
    return {
        "cells": cells, "dt": dt, "type": "reflecting_control" if reflecting else "continuous",
        "peak_bulk_pressure_preparation": amplitude,
        "energy_window_ratio_R": float(energy["reflected"] / energy["incident"]) if amplitude else None,
        "energy_window_ratio_T": float(energy["transmitted"] / energy["incident"]) if amplitude else None,
        "incident_window_energy": float(energy["incident"]),
        "max_relative_energy_drift": float(max(abs(energy[f"total_{t}"] / energy["total_0"] - 1)
                                                for t in (24, 96))) if amplitude else None,
        "ordered_cut_peak_stress": sample_peak(left), "bulk_cut_peak_stress": sample_peak(right),
        "peak_storage_rate": sample_peak(momentum_rate),
        "fixed_plane_peak_current": sample_peak(current),
        "fixed_plane_number_passing": float(simpson(current, x=t)),
        "peak_residual_force": sample_peak(transferred),
        "residual_force_impulse": float(simpson(left - right, x=t) - (momentum[-1] - momentum[0])),
        "independent_peak_force": sample_peak(independent_wall),
        "momentum_audit_error_over_peak_bulk_cut_stress":
            float(np.max(np.abs(transferred - independent_wall)) / error_scale) if amplitude else None,
    }


def main():
    print(json.dumps({"model": "native carrier, neutral flat held order profile at s=0",
                      "coordinates": "s increases from ordered brane to disordered bulk",
                      "prepared": "bulk-side pressure bump centered at 48 with sigma 4, initially zero velocity",
                      "box": [-192, 192], "cut": [-6, 6], "first_passage_end": 96,
                      "sample_spacing": 0.2, "pressure_scale": "n0 = c = 1"}), flush=True)
    for cells, dt, amplitudes in [(1536, 0.02, (0.001,)), (3072, 0.005, (0, 0.001, -0.001, 0.0005, -0.0005))]:
        for amplitude in amplitudes:
            for reflecting in (False, True):
                if reflecting and amplitude not in (0.001, -0.001):
                    continue
                result = run(cells, dt, amplitude, reflecting)
                print(json.dumps(result, allow_nan=False), flush=True)
                if amplitude:
                    assert result["max_relative_energy_drift"] < 1e-4
                    assert result["momentum_audit_error_over_peak_bulk_cut_stress"] < 1e-3
                    assert abs(result["energy_window_ratio_T"] + result["energy_window_ratio_R"] - 1) < 1e-3
                    if reflecting:
                        assert abs(result["fixed_plane_peak_current"]) < abs(amplitude) * 1e-6
                else:
                    assert abs(result["fixed_plane_peak_current"]) < 1e-10
                    assert abs(result["peak_residual_force"]) < 1e-10


if __name__ == "__main__":
    main()
