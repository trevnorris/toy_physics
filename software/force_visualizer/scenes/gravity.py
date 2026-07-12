"""Animated Newtonian-versus-1PN orbit and gravity return residual."""

from __future__ import annotations

from math import pi
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

from ..params import DEFAULT_PARAMS, ModelParameters, Provenance
from ..physics import gravity as core
from .data import build_gravity_data
from .shared import add_provenance_box, breathing_size, normalized_field, save_gif


def build_animation(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    include_25pn_benchmark: bool = False,
    frames: int = 90,
):
    """Build the gravity scene; departures are visible by default and toggleable."""

    params.validate()
    data = build_gravity_data(
        params, include_25pn_benchmark=include_25pn_benchmark
    )
    newtonian, relativistic = data.newtonian, data.relativistic
    body1, body2 = data.body1, data.body2
    semi_major_axis, eccentricity = data.semi_major_axis, data.eccentricity
    steps = len(relativistic.times) - 1
    sample = np.linspace(0, steps, frames, dtype=int)
    expected_shift = data.expected_shift
    departure = data.departure

    figure, axis = plt.subplots(figsize=(7.4, 6.2), constrained_layout=True)
    axis.set_aspect("equal")
    limit = semi_major_axis * (1.0 + eccentricity) * 1.05
    axis.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="x", ylabel="y")
    axis.set_title("Gravity: calibrated effective orbit (not a parent-PDE simulation)")
    grid_x, grid_y = np.meshgrid(
        np.linspace(-limit, limit, 17), np.linspace(-limit, limit, 17)
    )
    grid_points = np.stack((grid_x, grid_y), axis=-1)
    initial_sources = np.stack((body1[0], body2[0]))
    inflow = normalized_field(
        core.drain_inflow_field(grid_points, initial_sources, [0.03, 1.0], params.G)
    )
    inflow_arrows = axis.quiver(
        grid_x,
        grid_y,
        inflow[..., 0],
        inflow[..., 1],
        color="#2a9d8f",
        alpha=0.30,
        pivot="mid",
        scale=25,
        width=0.003,
        zorder=0,
    )
    tracer_angles = np.linspace(0.0, 2.0 * pi, 28, endpoint=False)
    tracers = np.column_stack((np.cos(tracer_angles), np.sin(tracer_angles))) * (0.94 * limit)
    tracer_artist, = axis.plot(
        tracers[:, 0], tracers[:, 1], ".", color="#006d77", ms=3.0,
        alpha=0.72, label="inward medium tracers", zorder=1
    )
    axis.plot(
        newtonian.positions[:, 0],
        newtonian.positions[:, 1],
        linestyle="--",
        color="0.60",
        linewidth=1.0,
        label="[DERIVED-FORM] 0PN closed reference",
    )
    pn_line, = axis.plot([], [], color="#7b2cbf", lw=1.8, label="[DERIVED-FORM] EIH 1PN orbit")
    orbiting_body, = axis.plot([], [], "o", color="#264653", ms=7)
    central_body, = axis.plot([], [], "o", color="#f4a261", ms=12)
    pulse = Circle((0.0, 0.0), 0.1, fill=False, color="#e63946", alpha=0.0, lw=1.3)
    axis.add_patch(pulse)
    axis.legend(loc="lower right", fontsize=8)
    axis.text(
        0.02,
        0.90,
        "medium inflow (gravity = the drain)\none-way sinks; tracers terminate at masses",
        transform=axis.transAxes,
        fontsize=8,
        color="#006d77",
    )
    add_provenance_box(axis, params, ("G", "c_gamma", "epsilon0"))

    order_label = "2.5PN benchmark" if include_25pn_benchmark else "1PN"
    shift_text = axis.text(
        0.02,
        0.02,
        f"[DERIVED-FORM] {order_label}; Δϖ₁PN={expected_shift:.4g} rad/orbit\n"
        f"[{Provenance.DERIVED_FORM.value}] p=2 attractive localization",
        transform=axis.transAxes,
        va="bottom",
        fontsize=8,
    )
    del shift_text
    if show_departure:
        axis.text(
            0.98,
            0.98,
            f"DEPARTURE: {departure.code}\n"
            f"[{Provenance.DERIVED_FORM.value}] ℓ=0/1 powers 1/3; bounded drain factor\n"
            f"[{Provenance.CALIBRATED_MAGNITUDE.value}] ε₀/(1+ε₀)="
            f"{departure.diagnostics['epsilon0/(1+epsilon0)']:.4g}",
            transform=axis.transAxes,
            va="top",
            ha="right",
            fontsize=7.5,
            color="#9d0208",
            bbox={"facecolor": "#fff1e6", "alpha": 0.9, "edgecolor": "#e63946"},
        )

    def update(frame_index: int):
        index = sample[frame_index]
        sources = np.stack((body1[index], body2[index]))
        field = normalized_field(
            core.drain_inflow_field(grid_points, sources, [0.03, 1.0], params.G)
        )
        inflow_arrows.set_UVC(field[..., 0], field[..., 1])
        tracer_field = normalized_field(
            core.drain_inflow_field(tracers, sources, [0.03, 1.0], params.G)
        )
        tracers[:] += 0.10 * tracer_field
        sink_distance = np.min(
            np.linalg.norm(tracers[:, None, :] - sources[None, :, :], axis=2), axis=1
        )
        consumed = sink_distance < 0.16
        if np.any(consumed):
            ids = np.flatnonzero(consumed)
            angles = 2.0 * pi * (ids + frame_index / max(frames, 1)) / len(tracers)
            tracers[ids] = 0.94 * limit * np.column_stack((np.cos(angles), np.sin(angles)))
        tracer_artist.set_data(tracers[:, 0], tracers[:, 1])
        pn_line.set_data(relativistic.positions[: index + 1, 0], relativistic.positions[: index + 1, 1])
        orbiting_body.set_data([body1[index, 0]], [body1[index, 1]])
        central_body.set_data([body2[index, 0]], [body2[index, 1]])
        breath_phase = 2.0 * pi * frame_index / 24.0
        orbiting_body.set_markersize(breathing_size(7.0, breath_phase))
        central_body.set_markersize(breathing_size(12.0, breath_phase + pi))
        if show_departure:
            phase = 2.0 * pi * frame_index / max(frames - 1, 1)
            pulse.set_radius(0.25 + (0.5 + 0.5 * np.sin(phase)) * 1.1)
            pulse.set_alpha(0.15 + 0.30 * (0.5 + 0.5 * np.cos(phase)))
        return pn_line, orbiting_body, central_body, pulse, inflow_arrows, tracer_artist

    animation = FuncAnimation(figure, update, frames=frames, interval=45, blit=False)
    return figure, animation


def save(
    output: str | Path,
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    show_departure: bool = True,
    include_25pn_benchmark: bool = False,
    frames: int = 90,
    fps: int = 20,
) -> Path:
    """Build and save the headless gravity GIF."""

    figure, animation = build_animation(
        params,
        show_departure=show_departure,
        include_25pn_benchmark=include_25pn_benchmark,
        frames=frames,
    )
    return save_gif(animation, figure, output, fps)
