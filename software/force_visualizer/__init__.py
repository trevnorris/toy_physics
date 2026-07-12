"""Reduced, calibrated four-sector phenomenology visualizer.

This package visualizes effective force laws.  It does not solve the parent
four-dimensional PDE or claim that the displayed forces emerge dynamically.
"""

from .params import DEFAULT_PARAMS, ModelParameters, Provenance

__all__ = ["DEFAULT_PARAMS", "ModelParameters", "Provenance"]
