"""Headless command-line renderer for all four force-sector animations."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/force_visualizer-mpl")

import matplotlib

matplotlib.use("Agg", force=True)

from .params import DEFAULT_PARAMS
from .scenes import charge, gravity, light, magnetism


def render_all(
    output_directory: str | Path,
    *,
    frames: int = 90,
    fps: int = 20,
    show_departures: bool = True,
    include_25pn_benchmark: bool = False,
) -> list[Path]:
    """Render all scenes with the Agg backend and return their GIF paths."""

    if frames < 2 or fps < 1:
        raise ValueError("frames must be >=2 and fps must be >=1")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    paths = [
        gravity.save(
            output / "gravity.gif",
            DEFAULT_PARAMS,
            show_departure=show_departures,
            include_25pn_benchmark=include_25pn_benchmark,
            frames=frames,
            fps=fps,
        ),
        light.save(
            output / "light.gif",
            DEFAULT_PARAMS,
            show_departure=show_departures,
            frames=frames,
            fps=fps,
        ),
        charge.save(
            output / "charge.gif",
            DEFAULT_PARAMS,
            show_departure=show_departures,
            frames=frames,
            fps=fps,
        ),
        magnetism.save(
            output / "magnetism.gif",
            DEFAULT_PARAMS,
            show_departure=show_departures,
            frames=frames,
            fps=fps,
        ),
    ]
    return paths


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
        help="directory for the four GIF files",
    )
    parser.add_argument("--frames", type=int, default=90)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument(
        "--hide-departures",
        action="store_true",
        help="counterfactual display toggle; default scenes expose all characterized departures",
    )
    parser.add_argument(
        "--include-25pn-benchmark",
        action="store_true",
        help="show the calibrated Burke--Thorne benchmark (native toy normalization is blocked)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""

    arguments = _parser().parse_args(argv)
    paths = render_all(
        arguments.output_dir,
        frames=arguments.frames,
        fps=arguments.fps,
        show_departures=not arguments.hide_departures,
        include_25pn_benchmark=arguments.include_25pn_benchmark,
    )
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
