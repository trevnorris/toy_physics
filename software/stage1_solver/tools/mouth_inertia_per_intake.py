"""Read slow brane reaction per isolated intake from the fixed-inward-speed radius study."""

import hashlib
import json
from decimal import Decimal
from pathlib import Path


ARCHIVE = Path("docs/_measurements/lateral_mouth_radius.json")


def read_fields(lines, prefix):
    line, = (line for line in lines if line.startswith(prefix))
    return dict(word.split("=") for word in line.split())


def main():
    runs = json.loads(ARCHIVE.read_text())["runs"]
    by_grid = {}
    for run in runs:
        lines = run["stdout"]
        digest = hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()
        if digest != run["stdout_sha256"] or not run["limits_verified"] or run["guard_exit"] or run["worker_exit"]:
            raise ValueError(f"Unverified archived run at {run['grid']}, radius {run['radius']}")
        isolated = read_fields(lines, "state_d=0 ")
        slow = read_fields(lines, "omega=0 ")
        states = by_grid.setdefault(run["grid"], {})
        radius = Decimal(str(run["radius"]))
        if radius in states:
            raise ValueError(f"Duplicate radius in {run['grid']}")
        states[radius] = tuple(Decimal(value) for value in (
            isolated["mean_speed"], isolated["target_rate"], slow["reaction"], slow["direct"]
        ))

    if set(by_grid) != {"e24_c4", "e24_c6", "e48_c4"}:
        raise ValueError("Expected the three archived grids")
    for grid, states in sorted(by_grid.items()):
        if sorted(states) != [Decimal("0.75"), Decimal("1"), Decimal("1.5")]:
            raise ValueError(f"Incomplete set of radii in {grid}")
        if len({speed for speed, _, _, _ in states.values()}) != 1:
            raise ValueError(f"Different mean inward speeds in {grid}")
        normalized = {}
        for radius, (_, rate, reaction, direct) in sorted(states.items()):
            normalized[radius] = reaction / rate, direct / rate
            print(f"{grid} radius={radius} D={rate} B={reaction} B/D={reaction/rate:.10g} "
                  f"direct_B/D={direct/rate:.10g}", flush=True)
        small = normalized[min(states)]
        large = normalized[max(states)]
        print(f"{grid} ratio_large/small={large[0]/small[0]:.10g} "
              f"required_delta_other_per_D={small[0]-large[0]:.10g} "
              f"direct_ratio={large[1]/small[1]:.10g} "
              f"direct_required_delta={small[1]-large[1]:.10g}", flush=True)


if __name__ == "__main__":
    main()
