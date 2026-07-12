"""Live interactive Matplotlib front end for the four verified force sectors.

Run locally with ``python -m force_visualizer.app`` from the repository's
``software/`` directory. Scenario setup is shared with the GIF scenes and
every trajectory/value is produced by the existing render-agnostic physics core.
"""

from __future__ import annotations

import os
import sys
from itertools import count
from math import pi
from typing import Any

import numpy as np

from .params import DEFAULT_PARAMS, ModelParameters


def _has_display() -> bool:
    """Return whether this process has a plausible local graphical display."""

    if sys.platform.startswith(("win", "darwin")):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def _select_interactive_backend() -> str:
    """Select an interactive backend before pyplot is imported."""

    if not _has_display():
        raise RuntimeError(
            "No graphical display was detected (DISPLAY/WAYLAND_DISPLAY is unset). "
            "Run the live app in a local desktop session; use render_all for headless GIF export."
        )
    import matplotlib

    interactive = {name.lower() for name in matplotlib.rcsetup.interactive_bk}
    if matplotlib.get_backend().lower() not in interactive:
        candidates = ["MacOSX"] if sys.platform == "darwin" else ["TkAgg", "QtAgg", "GTK3Agg"]
        errors: list[str] = []
        for candidate in candidates:
            try:
                matplotlib.use(candidate, force=True)
                return candidate
            except (ImportError, RuntimeError, ValueError) as error:
                errors.append(f"{candidate}: {error}")
        raise RuntimeError(
            "A display is present, but no interactive Matplotlib backend could be loaded ("
            + "; ".join(errors)
            + ")."
        )
    return str(matplotlib.get_backend())


class _SectorView:
    """Small common controller for one live sector and its Matplotlib widgets."""

    name = "sector"

    def __init__(self, app: "ForceVisualizerApp") -> None:
        self.app = app
        self.figure = app.figure
        self.params = app.params
        self.axes: list[Any] = []
        self.control_axes: list[Any] = []
        self.widgets: dict[str, Any] = {}
        self.departure_enabled = True
        self.revision = 0
        self.frame = 0

    def _axis(self, rectangle: list[float], **kwargs: Any):
        axis = self.figure.add_axes(rectangle, **kwargs)
        self.axes.append(axis)
        return axis

    def _slider(self, name: str, label: str, limits: tuple[float, float], initial: float, **kwargs: Any):
        from matplotlib.widgets import Slider

        row = len(self.widgets)
        axis = self.figure.add_axes([0.12, 0.205 - row * 0.052, 0.62, 0.028])
        self.control_axes.append(axis)
        widget = Slider(axis, label, limits[0], limits[1], valinit=initial, **kwargs)
        widget.on_changed(self._controls_changed)
        self.widgets[name] = widget
        return widget

    def _departure_widget(self):
        from matplotlib.widgets import CheckButtons

        axis = self.figure.add_axes([0.82, 0.54, 0.16, 0.09])
        self.control_axes.append(axis)
        widget = CheckButtons(axis, ["departures"], [self.departure_enabled])
        widget.on_clicked(self._departure_changed)
        self.departure_toggle = widget
        return widget

    def _controls_changed(self, _value: float) -> None:
        self.rebuild()

    def _departure_changed(self, _label: str) -> None:
        self.departure_enabled = bool(self.departure_toggle.get_status()[0])
        self.rebuild()

    def rebuild(self) -> None:
        raise NotImplementedError

    def update(self, frame: int) -> tuple[Any, ...]:
        raise NotImplementedError

    def snapshot(self) -> np.ndarray:
        """Return representative core-produced animated values for smoke tests."""

        raise NotImplementedError

    def set_visible(self, visible: bool) -> None:
        for axis in self.axes + self.control_axes:
            axis.set_visible(visible)

    def _finish_rebuild(self) -> None:
        self.frame = 0
        self.revision += 1
        self.update(0)
        self.figure.canvas.draw_idle()


class _GravityView(_SectorView):
    name = "Gravity"

    def __init__(self, app: "ForceVisualizerApp") -> None:
        super().__init__(app)
        from matplotlib.patches import Circle
        from .scenes.shared import normalized_field

        self.orbit_axis = self._axis([0.06, 0.30, 0.43, 0.63])
        self.radius_axis = self._axis([0.54, 0.30, 0.24, 0.63])
        self.orbit_axis.set(aspect="equal", xlabel="x", ylabel="y")
        self.orbit_axis.set_title("Gravity: 0PN reference and EIH 1PN")
        self.radius_axis.set(xlabel="time", ylabel="separation")
        self.radius_axis.set_title("Live orbital radius")
        self.reference_line, = self.orbit_axis.plot([], [], "--", color="0.65", lw=1)
        self.track, = self.orbit_axis.plot([], [], color="#7b2cbf", lw=1.8)
        self.body1, = self.orbit_axis.plot([], [], "o", color="#264653", ms=7)
        self.body2, = self.orbit_axis.plot([], [], "o", color="#f4a261", ms=12)
        grid_x, grid_y = np.meshgrid(np.linspace(-5.0, 5.0, 15), np.linspace(-5.0, 5.0, 15))
        self.field_grid = np.stack((grid_x, grid_y), axis=-1)
        self.inflow_arrows = self.orbit_axis.quiver(
            grid_x, grid_y, np.zeros_like(grid_x), np.zeros_like(grid_y),
            color="#2a9d8f", alpha=0.30, pivot="mid", scale=23, width=0.003, zorder=0
        )
        tracer_angles = np.linspace(0.0, 2.0 * pi, 24, endpoint=False)
        self.tracers = 4.5 * np.column_stack((np.cos(tracer_angles), np.sin(tracer_angles)))
        self.tracer_artist, = self.orbit_axis.plot(
            self.tracers[:, 0], self.tracers[:, 1], ".", color="#006d77", ms=3, alpha=0.72
        )
        self.inflow_label = self.orbit_axis.text(
            0.02, 0.88, "medium inflow (gravity = the drain)\none-way sink tracers",
            transform=self.orbit_axis.transAxes, fontsize=8, color="#006d77"
        )
        self.pulse = Circle((0.0, 0.0), 0.1, fill=False, color="#e63946", lw=1.4)
        self.orbit_axis.add_patch(self.pulse)
        self.radius_line, = self.radius_axis.plot([], [], color="#4361ee")
        self.status = self.radius_axis.text(0.04, 0.96, "", transform=self.radius_axis.transAxes, va="top", fontsize=8)
        self._slider("mass1", "m₁ [calibrated]", (0.01, 0.30), 0.03)
        self._slider("mass2", "m₂ [calibrated]", (0.4, 2.0), 1.0)
        self._slider("eccentricity", "e", (0.05, 0.70), 0.35)
        self._departure_widget()
        self.rebuild()

    def rebuild(self) -> None:
        from .scenes.data import build_gravity_data

        self.data = build_gravity_data(
            self.params,
            mass1=self.widgets["mass1"].val,
            mass2=self.widgets["mass2"].val,
            eccentricity=self.widgets["eccentricity"].val,
            periods=1.5,
            steps_per_period=260,
        )
        limit = self.data.semi_major_axis * (1.0 + self.data.eccentricity) * 1.08
        self.orbit_axis.set(xlim=(-limit, limit), ylim=(-limit, limit))
        self.reference_line.set_data(
            self.data.newtonian.positions[:, 0], self.data.newtonian.positions[:, 1]
        )
        radii = np.linalg.norm(self.data.relativistic.positions, axis=1)
        self.radius_axis.set_xlim(0.0, self.data.relativistic.times[-1])
        self.radius_axis.set_ylim(max(0.0, radii.min() * 0.90), radii.max() * 1.10)
        residual = self.data.departure.diagnostics["epsilon0/(1+epsilon0)"]
        self.status.set_text(
            f"[DERIVED-FORM] Δϖ={self.data.expected_shift:.3g} rad/orbit\n"
            f"departure ℓ=0/1 residual={residual:.3g}"
        )
        angles = np.linspace(0.0, 2.0 * pi, len(self.tracers), endpoint=False)
        self.tracers[:] = 0.94 * limit * np.column_stack((np.cos(angles), np.sin(angles)))
        self._finish_rebuild()

    def update(self, frame: int) -> tuple[Any, ...]:
        from .physics import gravity as core
        from .scenes.shared import breathing_size, normalized_field

        orbit = self.data.relativistic
        index = frame % len(orbit.times)
        sources = np.stack((self.data.body1[index], self.data.body2[index]))
        masses = [self.widgets["mass1"].val, self.widgets["mass2"].val]
        field = normalized_field(core.drain_inflow_field(self.field_grid, sources, masses, self.params.G))
        self.inflow_arrows.set_UVC(field[..., 0], field[..., 1])
        tracer_field = normalized_field(core.drain_inflow_field(self.tracers, sources, masses, self.params.G))
        self.tracers[:] += 0.09 * tracer_field
        distances = np.min(
            np.linalg.norm(self.tracers[:, None, :] - sources[None, :, :], axis=2), axis=1
        )
        consumed = distances < 0.16
        if np.any(consumed):
            ids = np.flatnonzero(consumed)
            limit = self.data.semi_major_axis * (1.0 + self.data.eccentricity) * 1.02
            angles = 2.0 * pi * (ids + frame / max(len(orbit.times), 1)) / len(self.tracers)
            self.tracers[ids] = limit * np.column_stack((np.cos(angles), np.sin(angles)))
        self.tracer_artist.set_data(self.tracers[:, 0], self.tracers[:, 1])
        self.track.set_data(orbit.positions[: index + 1, 0], orbit.positions[: index + 1, 1])
        self.body1.set_data([self.data.body1[index, 0]], [self.data.body1[index, 1]])
        self.body2.set_data([self.data.body2[index, 0]], [self.data.body2[index, 1]])
        breath_phase = 2.0 * pi * frame / 24.0
        self.body1.set_markersize(breathing_size(7.0, breath_phase))
        self.body2.set_markersize(breathing_size(12.0, breath_phase + pi))
        radii = np.linalg.norm(orbit.positions[: index + 1], axis=1)
        self.radius_line.set_data(orbit.times[: index + 1], radii)
        if self.departure_enabled:
            phase = 2.0 * pi * index / max(len(orbit.times) - 1, 1)
            self.pulse.set_radius(0.25 + 0.8 * (0.5 + 0.5 * np.sin(phase)))
            self.pulse.set_alpha(0.18 + 0.30 * (0.5 + 0.5 * np.cos(phase)))
        else:
            self.pulse.set_alpha(0.0)
        self.frame = frame
        return (
            self.track, self.body1, self.body2, self.radius_line, self.pulse,
            self.inflow_arrows, self.tracer_artist,
        )

    def snapshot(self) -> np.ndarray:
        residual = float(self.data.departure.diagnostics["epsilon0/(1+epsilon0)"])
        return np.r_[self.data.relativistic.positions[-1], residual if self.departure_enabled else 0.0]


class _LightView(_SectorView):
    name = "Light"

    def __init__(self, app: "ForceVisualizerApp") -> None:
        super().__init__(app)
        self.wave_axis = self._axis([0.06, 0.30, 0.45, 0.63])
        self.ray_axis = self._axis([0.56, 0.30, 0.22, 0.63])
        self.wave_axis.set(xlim=(0.0, 2.0 * pi), ylim=(-1.45, 1.45), xlabel="brane coordinate", ylabel="displacement")
        self.wave_axis.set_title("Exactly two transverse polarizations")
        self.ray_axis.set(aspect="equal", xlim=(-0.5, 4.0), ylim=(-20.5, 20.5), xlabel="x", ylabel="z")
        self.ray_axis.set_title("Refractive lens")
        self.first_line, = self.wave_axis.plot([], [], color="#0077b6", lw=2, label="T₁")
        self.second_line, = self.wave_axis.plot([], [], color="#90e0ef", lw=2, label="T₂")
        self.longitudinal_line, = self.wave_axis.plot([], [], "--", color="#d00000", lw=1.4, label="stray L")
        self.wave_axis.legend(loc="lower left", fontsize=8)
        self.lens_body, = self.ray_axis.plot(0.0, 0.0, "o", color="#f4a261", ms=12)
        self.ray_axis.plot([2.5, 2.5], [-20.0, 20.0], ":", color="0.65")
        self.ray_line, = self.ray_axis.plot([], [], color="#6a4c93", lw=2)
        self.status = self.ray_axis.text(0.04, 0.03, "", transform=self.ray_axis.transAxes, fontsize=8)
        self.wave_axis.text(
            0.02, 0.12, "transverse wave = brane shear field",
            transform=self.wave_axis.transAxes, fontsize=8, color="#0077b6"
        )
        self._slider("k", "k [derived dispersion]", (1.0, 5.0), 2.0, valstep=1.0)
        self._slider("c_gamma", "cγ [calibrated]", (5.0, 15.0), self.params.c_gamma)
        self._departure_widget()
        self.rebuild()

    def rebuild(self) -> None:
        from .scenes.data import build_light_data

        self.data = build_light_data(
            self.params,
            wavenumber=self.widgets["k"].val,
            c_gamma=self.widgets["c_gamma"].val,
            spatial_points=180,
            periods=1.0,
            ray_steps=600,
        )
        self.status.set_text(
            f"ω=cγk={self.data.c_gamma * self.data.wavenumber:.3g}\n"
            f"bend={self.data.ray.signed_deflection:.3g} rad"
        )
        self._finish_rebuild()

    def update(self, frame: int) -> tuple[Any, ...]:
        from .scenes.shared import breathing_size

        first = self.data.first_evolution
        index = frame % len(first.times)
        coordinate = first.coordinates
        self.first_line.set_data(coordinate, first.displacements[index])
        self.second_line.set_data(coordinate, self.data.second_evolution.displacements[index])
        if self.departure_enabled:
            self.longitudinal_line.set_data(
                coordinate,
                self.params.longitudinal_display_fraction
                * np.sin(
                    self.data.wavenumber * coordinate
                    - self.data.longitudinal_speed * self.data.wavenumber * first.times[index]
                ),
            )
        else:
            self.longitudinal_line.set_data([], [])
        ray_index = min(len(self.data.ray.positions) - 1, int(index * (len(self.data.ray.positions) - 1) / max(len(first.times) - 1, 1)))
        self.ray_line.set_data(
            self.data.ray.positions[: ray_index + 1, 0], self.data.ray.positions[: ray_index + 1, 1]
        )
        self.lens_body.set_markersize(breathing_size(12.0, 2.0 * pi * frame / 24.0))
        self.frame = frame
        return self.first_line, self.second_line, self.longitudinal_line, self.ray_line, self.lens_body

    def snapshot(self) -> np.ndarray:
        amplitude = self.params.longitudinal_display_fraction if self.departure_enabled else 0.0
        return np.array([
            self.data.c_gamma * self.data.wavenumber,
            self.data.ray.signed_deflection,
            amplitude,
        ])


class _ChargeView(_SectorView):
    name = "Charge"

    def __init__(self, app: "ForceVisualizerApp") -> None:
        super().__init__(app)
        self.motion_axis = self._axis([0.06, 0.30, 0.43, 0.63])
        self.law_axis = self._axis([0.54, 0.30, 0.24, 0.63])
        self.motion_axis.set(aspect="equal", xlabel="x", ylabel="y")
        self.motion_axis.set_title("Odd throat orientations")
        self.law_axis.set(xscale="log", yscale="log", xlabel="R", ylabel="|F|")
        self.law_axis.set_title("Coulomb + partner")
        self.left_track, = self.motion_axis.plot([], [], color="#0077b6", lw=1.2)
        self.right_track, = self.motion_axis.plot([], [], color="#d62828", lw=1.2)
        self.left_body, = self.motion_axis.plot([], [], "o", color="#0077b6", ms=12)
        self.right_body, = self.motion_axis.plot([], [], "o", color="#d62828", ms=12)
        grid_x, grid_y = np.meshgrid(np.linspace(-3.0, 3.0, 17), np.linspace(-0.75, 0.75, 9))
        self.field_grid = np.stack((grid_x, grid_y), axis=-1)
        self.field_arrows = self.motion_axis.quiver(
            grid_x, grid_y, np.zeros_like(grid_x), np.zeros_like(grid_y),
            color="#6a4c93", alpha=0.42, pivot="mid", scale=23, width=0.004, zorder=0
        )
        self.field_label = self.motion_axis.text(
            0.02, 0.04,
            "electric field (throat-body interaction)\n— a field, not a medium flow",
            transform=self.motion_axis.transAxes, fontsize=8, color="#5a189a"
        )
        self.coulomb_line, = self.law_axis.plot([], [], color="#2a9d8f", lw=2, label="p=2")
        self.yukawa_line, = self.law_axis.plot([], [], "--", color="#e76f51", lw=2, label="Yukawa")
        self.law_axis.legend(loc="lower left", fontsize=8)
        self.status = self.motion_axis.text(0.04, 0.96, "", transform=self.motion_axis.transAxes, va="top", fontsize=8)
        self._slider("sign1", "sign q₁", (-1.0, 1.0), 1.0, valstep=[-1.0, 1.0])
        self._slider("sign2", "sign q₂", (-1.0, 1.0), -1.0, valstep=[-1.0, 1.0])
        self._slider("Q_E", "Q_E [calibrated]", (0.5, 3.5), self.params.Q_E)
        self._departure_widget()
        self.rebuild()

    def rebuild(self) -> None:
        from .scenes.data import build_charge_data
        from .scenes.shared import sign_color

        orientations = (self.widgets["sign1"].val, self.widgets["sign2"].val)
        self.data = build_charge_data(
            self.params,
            orientations=orientations,
            Q_E=self.widgets["Q_E"].val,
            show_departure=self.departure_enabled,
            dt=0.003,
            steps=360,
        )
        positions = self.data.trajectories.positions
        limit = max(1.2, float(np.max(np.abs(positions[:, :, 0]))) * 1.12)
        self.motion_axis.set(xlim=(-limit, limit), ylim=(-0.8, 0.8))
        self.coulomb_line.set_data(self.data.distances, self.data.coulomb_force)
        if self.departure_enabled:
            self.yukawa_line.set_data(self.data.distances, self.data.corrected_force)
        else:
            self.yukawa_line.set_data([], [])
        self.law_axis.relim()
        self.law_axis.autoscale_view()
        for artist, sign in (
            (self.left_track, orientations[0]),
            (self.left_body, orientations[0]),
            (self.right_track, orientations[1]),
            (self.right_body, orientations[1]),
        ):
            artist.set_color(sign_color(sign))
        behavior = "repel" if orientations[0] * orientations[1] > 0 else "attract"
        self.status.set_text(f"[DERIVED-FORM] {behavior}; p=2\nQ_E={self.widgets['Q_E'].val:.3g}")
        self._finish_rebuild()

    def update(self, frame: int) -> tuple[Any, ...]:
        from .physics import charge as core
        from .scenes.shared import breathing_size, normalized_field

        positions = self.data.trajectories.positions
        index = frame % len(positions)
        orientations = (self.widgets["sign1"].val, self.widgets["sign2"].val)
        field = normalized_field(
            core.electric_field(
                self.field_grid, positions[index], orientations, self.widgets["Q_E"].val,
                self.params.ell, self.params.mouth_half_width_b
            )
        )
        self.field_arrows.set_UVC(field[..., 0], field[..., 1])
        self.left_track.set_data(positions[: index + 1, 0, 0], positions[: index + 1, 0, 1])
        self.right_track.set_data(positions[: index + 1, 1, 0], positions[: index + 1, 1, 1])
        self.left_body.set_data([positions[index, 0, 0]], [positions[index, 0, 1]])
        self.right_body.set_data([positions[index, 1, 0]], [positions[index, 1, 1]])
        phase = 2.0 * pi * frame / 24.0
        self.left_body.set_markersize(breathing_size(12.0, phase))
        self.right_body.set_markersize(breathing_size(12.0, phase + pi))
        self.frame = frame
        return self.left_track, self.right_track, self.left_body, self.right_body, self.field_arrows

    def snapshot(self) -> np.ndarray:
        return self.data.trajectories.positions[-1].ravel().copy()


class _MagnetismView(_SectorView):
    name = "Magnetism"

    def __init__(self, app: "ForceVisualizerApp") -> None:
        super().__init__(app)
        from .scenes.data import following_transverse_grid
        from .scenes.shared import current_marker

        self.motion_axis = self._axis([0.06, 0.30, 0.43, 0.63])
        self.channel_axis = self._axis([0.54, 0.30, 0.24, 0.63])
        self.motion_axis.set(
            aspect="equal",
            xlabel="transverse coordinate x₁",
            ylabel="transverse coordinate x₂",
        )
        self.motion_axis.set_title("End-on currents: core-driven transverse force")
        self.first_track, = self.motion_axis.plot([], [], color="#4361ee", lw=1.4)
        self.second_track, = self.motion_axis.plot([], [], color="#f72585", lw=1.4)
        self.first_body, = self.motion_axis.plot(
            [], [], marker=current_marker(1.0), linestyle="None", color="#4361ee", ms=13
        )
        self.second_body, = self.motion_axis.plot(
            [], [], marker=current_marker(1.0), linestyle="None", color="#f72585", ms=13
        )
        self.field_grid = following_transverse_grid(np.array([[-1.3, 0.0], [1.3, 0.0]]))
        self.field_arrows = self.motion_axis.quiver(
            self.field_grid[..., 0], self.field_grid[..., 1],
            np.zeros(self.field_grid.shape[:2]), np.zeros(self.field_grid.shape[:2]),
            color="#2a9d8f", alpha=0.42, pivot="mid", scale=23, width=0.004, zorder=0
        )
        self.current_key = self.motion_axis.text(
            0.98,
            0.96,
            "⊙ out of page   ⊗ into page\ncurrent is normal to the transverse plane",
            transform=self.motion_axis.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
        )
        self.field_label = self.motion_axis.text(
            0.02, 0.04,
            "magnetic field on the brane (moving throat's swirl, felt via localization)\n"
            "literal swirl lives in the throat 4D body",
            transform=self.motion_axis.transAxes, fontsize=7.5, color="#006d77"
        )
        self.channel_axis.set(xlim=(0.0, 1.0), ylim=(0.0, 1.0), xticks=[], yticks=[])
        self.channel_axis.set_title("Field content")
        self.channel_text = self.channel_axis.text(0.07, 0.92, "", va="top", fontsize=9)
        self._slider("aL", "aL [calibrated]", (0.0, 0.05), self.params.aL)
        self._slider("direction1", "I₁ (+out/−in)", (-1.0, 1.0), 1.0, valstep=[-1.0, 1.0])
        self._slider("direction2", "I₂ (+out/−in)", (-1.0, 1.0), 1.0, valstep=[-1.0, 1.0])
        self._departure_widget()
        self.rebuild()

    def rebuild(self) -> None:
        from .scenes.data import build_magnetism_data
        from .scenes.shared import current_color, current_marker

        directions = (self.widgets["direction1"].val, self.widgets["direction2"].val)
        self.data = build_magnetism_data(
            self.params,
            current_directions=directions,
            aL=self.widgets["aL"].val,
            show_departure=self.departure_enabled,
            dt=0.25,
            steps=360,
        )
        positions = self.data.trajectories.positions
        x_center = float(positions[:, :, 0].mean())
        x_half_extent = max(
            1.5,
            float(np.max(np.abs(positions[:, :, 0] - x_center))) + 0.55,
        )
        y_center = float(positions[:, :, 1].mean())
        self.motion_axis.set(
            xlim=(x_center - x_half_extent, x_center + x_half_extent),
            ylim=(y_center - 0.95, y_center + 0.95),
        )
        for artist, direction in (
            (self.first_track, directions[0]),
            (self.first_body, directions[0]),
            (self.second_track, directions[1]),
            (self.second_body, directions[1]),
        ):
            artist.set_color(current_color(direction))
        self.first_body.set_marker(current_marker(directions[0]))
        self.second_body.set_marker(current_marker(directions[1]))
        relation = "parallel: attract" if directions[0] == directions[1] else "antiparallel: repel"
        self.motion_axis.set_title(f"End-on currents — {relation}; transverse core force")
        scalar = self.data.scalar_ratio if self.departure_enabled else 0.0
        self.channel_text.set_text(
            "[DERIVED-FORM]\ntransverse vector channel\n2 transverse DOF\n\n"
            f"{relation}\npoint force p=2\n\n"
            f"scalar-current ratio={scalar:.3g}\n"
            + (
                "uncancelable scalar admixture\n(attractive for parallel currents)"
                if self.departure_enabled
                else "counterfactual scalar hidden"
            )
        )
        self._finish_rebuild()

    def update(self, frame: int) -> tuple[Any, ...]:
        from .physics import magnetism as core
        from .scenes.data import following_transverse_grid
        from .scenes.shared import breathing_size, normalized_field

        positions = self.data.trajectories.positions
        index = frame % len(positions)
        directions = (self.widgets["direction1"].val, self.widgets["direction2"].val)
        self.field_grid = following_transverse_grid(positions[index])
        field = normalized_field(
            core.magnetic_field_on_brane(
                self.field_grid, positions[index], directions, self.params.N_u,
                self.params.aT, self.widgets["aL"].val, self.params.mu_R,
                self.params.B_eff, self.departure_enabled
            )
        )
        self.field_arrows.set_offsets(self.field_grid.reshape(-1, 2))
        self.field_arrows.set_UVC(field[..., 0], field[..., 1])
        self.first_track.set_data(positions[: index + 1, 0, 0], positions[: index + 1, 0, 1])
        self.second_track.set_data(positions[: index + 1, 1, 0], positions[: index + 1, 1, 1])
        self.first_body.set_data([positions[index, 0, 0]], [positions[index, 0, 1]])
        self.second_body.set_data([positions[index, 1, 0]], [positions[index, 1, 1]])
        phase = 2.0 * pi * frame / 24.0
        self.first_body.set_markersize(breathing_size(13.0, phase))
        self.second_body.set_markersize(breathing_size(13.0, phase + pi))
        self.frame = frame
        return self.first_track, self.second_track, self.first_body, self.second_body, self.field_arrows

    def snapshot(self) -> np.ndarray:
        scalar = self.data.scalar_ratio if self.departure_enabled else 0.0
        return np.r_[self.data.trajectories.positions[-1].ravel(), scalar]


class ForceVisualizerApp:
    """Constructed live application; pass ``animate=False`` for Agg smoke tests."""

    def __init__(self, params: ModelParameters = DEFAULT_PARAMS, *, animate: bool = True) -> None:
        from matplotlib import pyplot as plt
        from matplotlib.animation import FuncAnimation
        from matplotlib.widgets import RadioButtons

        params.validate()
        self.params = params
        self.figure = plt.figure(figsize=(12.8, 7.8))
        self.figure.suptitle("Analog Four-Force Phenomenology Visualizer — calibrated effective forms", fontsize=14)
        self.sectors: dict[str, _SectorView] = {
            "gravity": _GravityView(self),
            "light": _LightView(self),
            "charge": _ChargeView(self),
            "magnetism": _MagnetismView(self),
        }
        self.selector_axis = self.figure.add_axes([0.82, 0.70, 0.16, 0.20])
        self.selector_axis.set_title("Sector")
        self.selector = RadioButtons(self.selector_axis, tuple(self.sectors))
        self.selector.on_clicked(self.select_sector)
        self.active_name = "gravity"
        self.select_sector(self.active_name)
        self.animation = None
        if animate:
            self.animation = FuncAnimation(
                self.figure,
                self._tick,
                frames=count(),
                interval=40,
                blit=False,
                cache_frame_data=False,
            )

    @property
    def active_sector(self) -> _SectorView:
        return self.sectors[self.active_name]

    def select_sector(self, name: str) -> None:
        key = name.lower()
        if key not in self.sectors:
            raise ValueError(f"unknown sector {name!r}")
        self.active_name = key
        for sector_name, sector in self.sectors.items():
            sector.set_visible(sector_name == key)
        self.selector_axis.set_visible(True)
        self.figure.canvas.draw_idle()

    def _tick(self, frame: int) -> tuple[Any, ...]:
        return self.active_sector.update(frame)

    def close(self) -> None:
        from matplotlib import pyplot as plt

        if self.animation is not None:
            self.animation.event_source.stop()
        plt.close(self.figure)


def create_app(
    params: ModelParameters = DEFAULT_PARAMS, *, animate: bool = True
) -> ForceVisualizerApp:
    """Construct the app using the caller-selected backend (Agg is test-safe)."""

    return ForceVisualizerApp(params, animate=animate)


def main() -> int:
    """Select a desktop backend and run the local interactive window."""

    try:
        _select_interactive_backend()
        from matplotlib import pyplot as plt

        app = create_app(animate=True)
        plt.show()
    except (ImportError, RuntimeError) as error:
        print(f"force_visualizer live app unavailable: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
