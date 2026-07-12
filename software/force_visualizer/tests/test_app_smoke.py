"""Agg smoke coverage for the live app's widgets and core-backed callbacks."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np
from matplotlib.animation import FuncAnimation

from software.force_visualizer.app import create_app
from software.force_visualizer.scenes.shared import (
    FORWARD_CURRENT_COLOR,
    NEGATIVE_COLOR,
    POSITIVE_COLOR,
    REVERSE_CURRENT_COLOR,
    current_marker,
)


def test_live_sector_widgets_update_core_quantities_headlessly() -> None:
    cases = (
        ("gravity", "mass2", 1.35),
        ("light", "k", 3.0),
        ("charge", "Q_E", 2.75),
        ("magnetism", "aL", 0.035),
    )
    app = create_app(animate=True)
    try:
        assert isinstance(app.animation, FuncAnimation)
        charge_sector = app.sectors["charge"]
        assert charge_sector.widgets["sign1"].val == 1.0
        assert charge_sector.widgets["sign2"].val == -1.0
        assert charge_sector.left_body.get_color() == POSITIVE_COLOR
        assert charge_sector.right_body.get_color() == NEGATIVE_COLOR
        charge_sector.widgets["sign2"].set_val(1.0)
        assert charge_sector.right_body.get_color() == POSITIVE_COLOR
        charge_sector.widgets["sign2"].set_val(-1.0)
        assert charge_sector.right_body.get_color() == NEGATIVE_COLOR

        magnetic_sector = app.sectors["magnetism"]
        assert magnetic_sector.first_body.get_color() == FORWARD_CURRENT_COLOR
        assert magnetic_sector.first_body.get_marker() == current_marker(1.0)
        magnetic_sector.widgets["direction1"].set_val(-1.0)
        assert magnetic_sector.first_body.get_color() == REVERSE_CURRENT_COLOR
        assert magnetic_sector.first_body.get_marker() == current_marker(-1.0)
        magnetic_sector.update(len(magnetic_sector.data.trajectories.positions) - 1)
        positions = magnetic_sector.data.trajectories.positions[-1]
        offsets = magnetic_sector.field_arrows.get_offsets()
        assert offsets[:, 0].min() < positions[:, 0].min()
        assert offsets[:, 0].max() > positions[:, 0].max()
        x_limits = magnetic_sector.motion_axis.get_xlim()
        assert x_limits[0] < positions[:, 0].min() < positions[:, 0].max() < x_limits[1]
        magnetic_sector.widgets["direction1"].set_val(1.0)

        for sector_name, slider_name, new_value in cases:
            app.select_sector(sector_name)
            sector = app.sectors[sector_name]
            initial_revision = sector.revision
            initial = sector.snapshot()

            sector.widgets[slider_name].set_val(new_value)
            after_slider = sector.snapshot()
            assert sector.revision == initial_revision + 1
            assert not np.allclose(after_slider, initial)

            sector.departure_toggle.set_active(0)
            after_toggle = sector.snapshot()
            assert sector.revision == initial_revision + 2
            assert not np.allclose(after_toggle, after_slider)

            # Advance the live FuncAnimation itself, then inspect the artists
            # returned by its update callback.
            assert app.animation._step()
            artists = app.animation._drawn_artists
            assert artists
    finally:
        app.close()
