"""Pure deterministic physics core; rendering dependencies are forbidden.

Architecture source: ``software/force_visualizer/notes/build_spec.md`` §2.
Sector law sources are cited in every submodule and public function.
"""

from . import charge, gravity, integrators, light, magnetism

__all__ = ["charge", "gravity", "integrators", "light", "magnetism"]
