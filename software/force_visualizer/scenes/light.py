"""Animated two-polarization transverse wave and weak-field lensing ray."""

from __future__ import annotations

from math import pi
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from ..params import DEFAULT_PARAMS, ModelParameters, Provenance
from ..physics import light as core
from .data import build_light_data
from .shared import add_provenance_box, breathing_size, save_gif


def build_animation(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    frames: int = 90,
):
    """Build the light scene with a toggleable stray longitudinal overlay."""

    params.validate()
    data = build_light_data(params)
    wavenumber = data.wavenumber
    first_evolution = data.first_evolution
    second_evolution = data.second_evolution
    coordinate = first_evolution.coordinates
    ray = data.ray
    departure = data.departure
    longitudinal_speed = data.longitudinal_speed

    figure, (wave_axis, ray_axis) = plt.subplots(1, 2, figsize=(11.2, 5.1), constrained_layout=True)
    figure.suptitle("Light: MacCullagh transverse waves and calibrated 1PN refractive lens")
    wave_axis.set(xlim=(coordinate[0], coordinate[-1]), ylim=(-1.45, 1.45), xlabel="brane coordinate", ylabel="displacement")
    wave_axis.set_title("[DERIVED-FORM] exactly 2 transverse polarizations")
    wave_axis.text(
        0.02, 0.12, "transverse wave = brane shear field",
        transform=wave_axis.transAxes, fontsize=8, color="#0077b6"
    )
    first_line, = wave_axis.plot([], [], color="#0077b6", lw=2, label="T₁")
    second_line, = wave_axis.plot([], [], color="#90e0ef", lw=2, label="T₂ (phase shifted)")
    longitudinal_line, = wave_axis.plot([], [], "--", color="#d00000", lw=1.4, label="stray longitudinal")
    wave_axis.legend(loc="lower left", fontsize=8)
    add_provenance_box(wave_axis, params, ("c_gamma", "lambda_gamma", "rho_br", "mu_R"))

    ray_axis.set_aspect("equal")
    ray_axis.set(xlim=(-0.5, 4.0), ylim=(-20.5, 20.5), xlabel="impact coordinate x", ylabel="propagation z")
    ray_axis.set_title("[DERIVED-FORM] n(r)=1+2GM/(cγ²r); bends toward mass")
    lens_body, = ray_axis.plot(0.0, 0.0, "o", color="#f4a261", ms=13)
    ray_axis.plot([2.5, 2.5], [-20.0, 20.0], ":", color="0.65", lw=1, label="unlensed")
    ray_line, = ray_axis.plot([], [], color="#6a4c93", lw=2, label="Hamiltonian ray")
    ray_axis.legend(loc="lower right", fontsize=8)
    add_provenance_box(ray_axis, params, ("G", "c_gamma"))
    ray_axis.text(
        0.02,
        0.02,
        f"[{Provenance.DERIVED_FORM.value}] signed bend={ray.signed_deflection:.4g} rad",
        transform=ray_axis.transAxes,
        fontsize=8,
    )
    if show_departure:
        wave_axis.text(
            0.98,
            0.98,
            f"DEPARTURE: {departure.code}\n"
            f"[{Provenance.DERIVED_FORM.value}] 1 longitudinal DOF; cL={longitudinal_speed:.4g}\n"
            f"[{Provenance.CALIBRATED_MAGNITUDE.value}] display amplitude="
            f"{params.longitudinal_display_fraction:.4g}\nMaxwell locus only BY TUNING: B_eff→0",
            transform=wave_axis.transAxes,
            va="top",
            ha="right",
            fontsize=7.2,
            color="#9d0208",
            bbox={"facecolor": "#fff1e6", "alpha": 0.9, "edgecolor": "#d00000"},
        )

    wave_samples = np.linspace(0, len(first_evolution.times) - 1, frames, dtype=int)
    ray_samples = np.linspace(0, len(ray.positions) - 1, frames, dtype=int)

    def update(frame_index: int):
        wave_index = wave_samples[frame_index]
        time = first_evolution.times[wave_index]
        first_line.set_data(coordinate, first_evolution.displacements[wave_index])
        second_line.set_data(coordinate, second_evolution.displacements[wave_index])
        if show_departure:
            longitudinal_line.set_data(
                coordinate,
                params.longitudinal_display_fraction
                * np.sin(wavenumber * coordinate - longitudinal_speed * wavenumber * time),
            )
        else:
            longitudinal_line.set_data([], [])
        index = ray_samples[frame_index]
        ray_line.set_data(ray.positions[: index + 1, 0], ray.positions[: index + 1, 1])
        lens_body.set_markersize(breathing_size(13.0, 2.0 * pi * frame_index / 24.0))
        return first_line, second_line, longitudinal_line, ray_line, lens_body

    animation = FuncAnimation(figure, update, frames=frames, interval=45, blit=False)
    return figure, animation


def save(
    output: str | Path,
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    frames: int = 90,
    fps: int = 20,
) -> Path:
    """Build and save the headless light GIF."""

    figure, animation = build_animation(params, show_departure=show_departure, frames=frames)
    return save_gif(animation, figure, output, fps)
