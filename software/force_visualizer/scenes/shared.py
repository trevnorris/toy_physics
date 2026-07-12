"""Shared rendering helpers.  This module is intentionally outside physics/."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/force_visualizer-mpl")

from matplotlib import pyplot as plt
from matplotlib.animation import Animation, PillowWriter

from ..params import ModelParameters, labeled_value


POSITIVE_COLOR = "#d62828"
NEGATIVE_COLOR = "#277da1"
FORWARD_CURRENT_COLOR = "#e76f51"
REVERSE_CURRENT_COLOR = "#4361ee"


def sign_color(value: float) -> str:
    """Return the required red-positive / blue-negative charge color."""

    return POSITIVE_COLOR if value > 0.0 else NEGATIVE_COLOR


def current_color(value: float) -> str:
    """Use distinct colors for the two current directions."""

    return FORWARD_CURRENT_COLOR if value > 0.0 else REVERSE_CURRENT_COLOR


def current_marker(value: float) -> str:
    """Return the conventional end-on marker: circled dot out, cross into page."""

    return r"$\odot$" if value > 0.0 else r"$\otimes$"


def normalized_field(field):
    """Normalize a core-sampled vector field for truthful direction display."""

    vectors = np.asarray(field, dtype=float)
    magnitude = np.linalg.norm(vectors, axis=-1, keepdims=True)
    return np.divide(vectors, magnitude, out=np.zeros_like(vectors), where=magnitude > 0.0)


def breathing_size(base_size: float, phase: float) -> float:
    """Subtle non-quantitative breathing cue (plus/minus six percent)."""

    return float(base_size * (1.0 + 0.06 * np.sin(phase)))


def provenance_text(params: ModelParameters, names: Iterable[str]) -> str:
    """Format the mandatory status tag for every numerical scene value."""

    return "\n".join(labeled_value(params, name) for name in names)


def add_provenance_box(axis, params: ModelParameters, names: Iterable[str], *, location: Tuple[float, float] = (0.02, 0.98)):
    """Attach a compact provenance legend in axes coordinates."""

    return axis.text(
        location[0],
        location[1],
        provenance_text(params, names),
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=7,
        family="monospace",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )


def save_gif(animation: Animation, figure, output: str | Path, fps: int = 20) -> Path:
    """Save an animation headlessly with Pillow and close its figure."""

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    animation.save(path, writer=PillowWriter(fps=fps))
    plt.close(figure)
    return path
