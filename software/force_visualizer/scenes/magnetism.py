"""Animated parallel-current attraction and scalar-current admixture."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from ..params import DEFAULT_PARAMS, ModelParameters, Provenance
from ..physics import magnetism as core
from .data import build_magnetism_data, following_transverse_grid
from .shared import (
    add_provenance_box,
    breathing_size,
    current_color,
    current_marker,
    normalized_field,
    save_gif,
)


def build_animation(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    current_directions: tuple[float, float] = (1.0, 1.0),
    frames: int = 90,
):
    """Build the end-on-current scene with core-driven transverse motion."""

    params.validate()
    data = build_magnetism_data(
        params,
        current_directions=current_directions,
        show_departure=show_departure,
    )
    trajectories = data.trajectories
    sample = np.linspace(0, len(trajectories.times) - 1, frames, dtype=int)
    departure = data.departure
    ratio = data.scalar_ratio
    directions = data.current_directions
    positions = trajectories.positions
    relation = "parallel currents attract" if directions[0] * directions[1] > 0.0 else "antiparallel currents repel"

    figure, (motion_axis, channel_axis) = plt.subplots(1, 2, figsize=(10.8, 4.9), constrained_layout=True)
    figure.suptitle("Magnetism: end-on currents and transverse core-force response")
    motion_axis.set_aspect("equal")
    x_padding = 0.65
    x_limits = (
        float(positions[:, :, 0].min() - x_padding),
        float(positions[:, :, 0].max() + x_padding),
    )
    y_center = float(positions[:, :, 1].mean())
    motion_axis.set(
        xlim=x_limits,
        ylim=(y_center - 1.0, y_center + 1.0),
        xlabel="transverse coordinate x₁",
        ylabel="transverse coordinate x₂",
    )
    motion_axis.set_title(f"[DERIVED-FORM] {relation}; point force p=2")
    grid_points = following_transverse_grid(positions[0])
    initial_field = normalized_field(
        core.magnetic_field_on_brane(
            grid_points, positions[0], directions, params.N_u,
            params.aT, params.aL, params.mu_R, params.B_eff, show_departure
        )
    )
    field_arrows = motion_axis.quiver(
        grid_points[..., 0], grid_points[..., 1], initial_field[..., 0], initial_field[..., 1],
        color="#2a9d8f", alpha=0.42, pivot="mid", scale=23, width=0.004, zorder=0
    )
    first_track, = motion_axis.plot([], [], color=current_color(directions[0]), lw=1.4)
    second_track, = motion_axis.plot([], [], color=current_color(directions[1]), lw=1.4)
    first_body, = motion_axis.plot(
        [], [], marker=current_marker(directions[0]), linestyle="None",
        color=current_color(directions[0]), ms=13
    )
    second_body, = motion_axis.plot(
        [], [], marker=current_marker(directions[1]), linestyle="None",
        color=current_color(directions[1]), ms=13
    )
    motion_axis.text(
        0.98,
        0.96,
        "⊙ out of page   ⊗ into page\ncurrent is normal to this transverse plane",
        transform=motion_axis.transAxes,
        ha="right",
        va="top",
        fontsize=7.4,
    )
    motion_axis.text(
        0.02,
        0.04,
        "magnetic field on the brane (moving throat's swirl, felt via localization)\n"
        "literal swirl lives in the throat 4D body",
        transform=motion_axis.transAxes,
        fontsize=7.4,
        color="#006d77",
    )
    add_provenance_box(motion_axis, params, ("aT", "aL", "mu_R", "B_eff"))

    channel_axis.set(xlim=(0.0, 1.0), ylim=(0.0, 1.0), xticks=[], yticks=[])
    channel_axis.set_title("Derived field content")
    channel_axis.text(
        0.08,
        0.76,
        "[DERIVED-FORM] transverse vector channel\n"
        "U_T ∝ -(D+A)/(μ_R R)\n"
        "2 transverse physical DOF",
        fontsize=10,
        color="#264653",
    )
    channel_axis.text(
        0.08,
        0.50,
        "[DERIVED-FORM] longitudinal scalar-current channel\n"
        "U_L ∝ -(D-A)/(B_eff R)\n"
        "same channel sign; cannot cancel\n"
        "(attractive for parallel currents)",
        fontsize=10,
        color="#9d0208" if show_departure else "0.6",
        alpha=1.0 if show_departure else 0.45,
    )
    add_provenance_box(channel_axis, params, ("c_E", "c_gamma", "C_hu"), location=(0.08, 0.30))
    if show_departure:
        channel_axis.text(
            0.08,
            0.06,
            f"DEPARTURE: {departure.code}\n"
            f"[{Provenance.DERIVED_FORM.value}] scalar/vector magnitude ratio form\n"
            f"[{Provenance.CALIBRATED_MAGNITUDE.value}] transverse=1, scalar={ratio:.4g}\n"
            "transverse channel leads; scalar remains visible\n"
            f"[{Provenance.FREE_UNREDUCED.value}] C_hu={params.C_hu:.4g}; total physical DOF=4",
            fontsize=8,
            color="#9d0208",
            bbox={"facecolor": "#fff1e6", "alpha": 0.9, "edgecolor": "#d00000"},
        )

    def update(frame_index: int):
        index = sample[frame_index]
        grid_points = following_transverse_grid(positions[index])
        field = normalized_field(
            core.magnetic_field_on_brane(
                grid_points, positions[index], directions, params.N_u,
                params.aT, params.aL, params.mu_R, params.B_eff, show_departure
            )
        )
        field_arrows.set_offsets(grid_points.reshape(-1, 2))
        field_arrows.set_UVC(field[..., 0], field[..., 1])
        first_track.set_data(positions[: index + 1, 0, 0], positions[: index + 1, 0, 1])
        second_track.set_data(positions[: index + 1, 1, 0], positions[: index + 1, 1, 1])
        first_body.set_data([positions[index, 0, 0]], [positions[index, 0, 1]])
        second_body.set_data([positions[index, 1, 0]], [positions[index, 1, 1]])
        phase = 2.0 * np.pi * frame_index / 24.0
        first_body.set_markersize(breathing_size(13.0, phase))
        second_body.set_markersize(breathing_size(13.0, phase + np.pi))
        return first_track, second_track, first_body, second_body, field_arrows

    animation = FuncAnimation(figure, update, frames=frames, interval=45, blit=False)
    return figure, animation


def save(
    output: str | Path,
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    current_directions: tuple[float, float] = (1.0, 1.0),
    frames: int = 90,
    fps: int = 20,
) -> Path:
    """Build and save the headless magnetism GIF."""

    figure, animation = build_animation(
        params,
        show_departure=show_departure,
        current_directions=current_directions,
        frames=frames,
    )
    return save_gif(animation, figure, output, fps)
