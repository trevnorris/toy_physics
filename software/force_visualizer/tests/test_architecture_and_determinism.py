"""Hard purity, provenance, deterministic, and scene-build acceptance checks."""

from __future__ import annotations

import ast
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib import pyplot as plt

from software.force_visualizer.params import (
    DEFAULT_PARAMS as P,
    PARAMETER_INFO,
    Provenance,
    parameter_names,
)
from software.force_visualizer.physics import charge, gravity
from software.force_visualizer.scenes import charge as charge_scene
from software.force_visualizer.scenes import gravity as gravity_scene
from software.force_visualizer.scenes import light as light_scene
from software.force_visualizer.scenes import magnetism as magnetism_scene


ROOT = Path(__file__).resolve().parents[1]


def test_physics_package_imports_no_rendering_module() -> None:
    forbidden_roots = {"matplotlib", "seaborn", "plotly", "bokeh", "scenes", "rendering"}
    for source in (ROOT / "physics").glob("*.py"):
        tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported = {alias.name.split(".")[0] for alias in node.names}
                assert imported.isdisjoint(forbidden_roots), (source, imported)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_root = node.module.lstrip(".").split(".")[0]
                assert imported_root not in forbidden_roots, (source, imported_root)


def test_every_shared_parameter_has_an_explicit_status_and_consistent_identity() -> None:
    for name in parameter_names():
        assert name in PARAMETER_INFO
        assert PARAMETER_INFO[name].status in set(Provenance)
    for derived_name in ("c_gamma", "B_eff", "N0", "yukawa_mass"):
        assert PARAMETER_INFO[derived_name].status is Provenance.DERIVED_FORM
    np.testing.assert_allclose(P.c_gamma**2, P.mu_R / P.rho_br, atol=1e-14)
    np.testing.assert_allclose(P.B_eff, P.rho_B0**2 / P.chi_c, atol=1e-14)
    P.validate()


def test_same_inputs_produce_bitwise_identical_trajectories() -> None:
    position, velocity = gravity.kepler_periapsis_state(2.5, 0.2, 1.2, P.G)
    first = gravity.integrate_relative_orbit(
        position, velocity, 1.0, 0.2, P.G, P.c_gamma, 0.002, 300, pn_order=1.0
    )
    second = gravity.integrate_relative_orbit(
        position, velocity, 1.0, 0.2, P.G, P.c_gamma, 0.002, 300, pn_order=1.0
    )
    assert np.array_equal(first.positions, second.positions)
    assert np.array_equal(first.velocities, second.velocities)

    kwargs = dict(
        positions=[[-1.0, 0.0], [1.0, 0.0]],
        velocities=[[0.0, 0.0], [0.0, 0.0]],
        masses=[1.0, 1.0],
        orientations=[1.0, -1.0],
        Q_E=P.Q_E,
        ell=P.ell,
        mouth_half_width=P.mouth_half_width_b,
        dt=0.002,
        steps=200,
    )
    assert np.array_equal(charge.integrate_charges(**kwargs).positions, charge.integrate_charges(**kwargs).positions)


def test_all_scenes_build_headlessly_with_departure_toggle() -> None:
    assert matplotlib.get_backend().lower() == "agg"
    for module in (gravity_scene, light_scene, charge_scene, magnetism_scene):
        figure, animation = module.build_animation(P, show_departure=True, frames=2)
        assert isinstance(animation, FuncAnimation)
        animation._draw_next_frame(0, blit=False)  # exercise the update and mark it rendered
        plt.close(figure)
