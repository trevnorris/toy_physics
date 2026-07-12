"""Animated localized Coulomb interaction with toggleable Yukawa correction."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from ..params import DEFAULT_PARAMS, ModelParameters, Provenance
from ..physics import charge as core
from .data import build_charge_data
from .shared import add_provenance_box, breathing_size, normalized_field, save_gif, sign_color


def build_animation(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    frames: int = 90,
):
    """Build the default opposite-charge attraction and force-distance law."""

    params.validate()
    orientations = (1.0, -1.0)
    data = build_charge_data(params, orientations=orientations, show_departure=show_departure)
    trajectories = data.trajectories
    departure = data.departure
    sample = np.linspace(0, len(trajectories.times) - 1, frames, dtype=int)
    distances = data.distances
    coulomb_force = data.coulomb_force
    gap = float(departure.diagnostics["mass_gap"])
    corrected_force = data.corrected_force

    figure, (motion_axis, law_axis) = plt.subplots(1, 2, figsize=(10.8, 4.9), constrained_layout=True)
    figure.suptitle("Electric charge: odd throat orientation, localized Coulomb form")
    motion_axis.set_aspect("equal")
    limit = max(2.6, float(np.max(np.abs(trajectories.positions)))) * 1.05
    motion_axis.set(xlim=(-limit, limit), ylim=(-1.5, 1.5), xlabel="x", ylabel="y")
    motion_axis.set_title("[DERIVED-FORM] like repel / unlike attract")
    grid_x, grid_y = np.meshgrid(
        np.linspace(-limit, limit, 17), np.linspace(-1.35, 1.35, 11)
    )
    grid_points = np.stack((grid_x, grid_y), axis=-1)
    initial_field = normalized_field(
        core.electric_field(
            grid_points,
            trajectories.positions[0],
            orientations,
            params.Q_E,
            params.ell,
            params.mouth_half_width_b,
        )
    )
    field_arrows = motion_axis.quiver(
        grid_x, grid_y, initial_field[..., 0], initial_field[..., 1],
        color="#6a4c93", alpha=0.42, pivot="mid", scale=23, width=0.004, zorder=0
    )
    left_track, = motion_axis.plot([], [], color=sign_color(orientations[0]), lw=1.2)
    right_track, = motion_axis.plot([], [], color=sign_color(orientations[1]), lw=1.2)
    left_body, = motion_axis.plot([], [], "o", color=sign_color(orientations[0]), ms=12, label="positive (+w)")
    right_body, = motion_axis.plot([], [], "o", color=sign_color(orientations[1]), ms=12, label="negative (-w)")
    motion_axis.legend(loc="lower right", fontsize=8)
    add_provenance_box(motion_axis, params, ("Q_E", "ell", "mouth_half_width_b", "N0"))
    motion_axis.text(
        0.02,
        0.04,
        "electric field (throat-body interaction)\n— a field, not a medium flow",
        transform=motion_axis.transAxes,
        fontsize=8,
        color="#5a189a",
    )

    law_axis.loglog(distances, coulomb_force, color="#2a9d8f", lw=2, label="[DERIVED-FORM] Coulomb p=2")
    if show_departure:
        law_axis.loglog(distances, corrected_force, "--", color="#e76f51", lw=2, label="Coulomb + Yukawa")
    law_axis.set(xlabel="separation R", ylabel="force magnitude", title="Force falloff and short-range partner")
    law_axis.legend(loc="lower left", fontsize=8)
    add_provenance_box(law_axis, params, ("yukawa_mass", "yukawa_fraction"))
    if show_departure:
        law_axis.text(
            0.98,
            0.98,
            f"DEPARTURE: {departure.code}\n"
            f"[{Provenance.DERIVED_FORM.value}] exp(-√3 R/ℓ)/R; gap={gap:.4g}\n"
            f"[{Provenance.CALIBRATED_MAGNITUDE.value}] relative residue="
            f"{params.yukawa_fraction:.4g}",
            transform=law_axis.transAxes,
            va="top",
            ha="right",
            fontsize=7.3,
            color="#9d0208",
            bbox={"facecolor": "#fff1e6", "alpha": 0.9, "edgecolor": "#e76f51"},
        )

    def update(frame_index: int):
        index = sample[frame_index]
        positions = trajectories.positions
        field = normalized_field(
            core.electric_field(
                grid_points,
                positions[index],
                orientations,
                params.Q_E,
                params.ell,
                params.mouth_half_width_b,
            )
        )
        field_arrows.set_UVC(field[..., 0], field[..., 1])
        left_track.set_data(positions[: index + 1, 0, 0], positions[: index + 1, 0, 1])
        right_track.set_data(positions[: index + 1, 1, 0], positions[: index + 1, 1, 1])
        left_body.set_data([positions[index, 0, 0]], [positions[index, 0, 1]])
        right_body.set_data([positions[index, 1, 0]], [positions[index, 1, 1]])
        phase = 2.0 * np.pi * frame_index / 24.0
        left_body.set_markersize(breathing_size(12.0, phase))
        right_body.set_markersize(breathing_size(12.0, phase + np.pi))
        return left_track, right_track, left_body, right_body, field_arrows

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
    """Build and save the headless charge GIF."""

    figure, animation = build_animation(params, show_departure=show_departure, frames=frames)
    return save_gif(animation, figure, output, fps)
