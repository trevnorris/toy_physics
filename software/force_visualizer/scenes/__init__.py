"""Matplotlib scene layer, loaded lazily so the live app can choose a backend."""

from __future__ import annotations

from importlib import import_module

__all__ = ["charge", "gravity", "light", "magnetism"]


def __getattr__(name: str):
    if name in __all__:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(name)
