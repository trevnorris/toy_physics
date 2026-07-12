"""Render-agnostic records for the four characterized model departures.

Honest-scope contract: ``software/force_visualizer/notes/build_spec.md`` §5.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class CharacterizedDeparture:
    """A named departure with explicit form/magnitude epistemic labels.

    Contract source: ``software/force_visualizer/notes/build_spec.md`` §5.
    """

    code: str
    description: str
    derived_form: str
    calibrated_magnitude: str
    diagnostics: Mapping[str, float | str]
