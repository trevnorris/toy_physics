"""Fixed-phase carrier response when only the terminal conversion disk conducts."""
import argparse
import json
from pathlib import Path

from mouth_exchange_carrier import BIAS
from mouth_exchange_carrier_density import solve
from mouth_exchange_geometry import ExchangeShape


RECORD = (Path(__file__).resolve().parents[3]
          / "research/pde_ledger_v3/_measurements/mouth_exchange_carrier.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, required=True, choices=(40, 64))
    mesh = parser.parse_args().mesh
    previous = json.loads(RECORD.read_text())["finite"][f"mesh{mesh}_drive1"]
    old_center, = (row for row in previous["states"] if row["step"] == 0)
    steps = (-.005, 0., .005) if mesh == 40 else (-.010, -.005, 0., .005, .010)
    states = {}
    for step in steps:
        row = solve(step, mesh, BIAS, ExchangeShape, sealed_side=True)
        cuts = row["current_by_cut"]
        if set(cuts) != {"cap_layer"} or abs(cuts["cap_layer"] - row["current"]) > 1e-7:
            raise ArithmeticError("Terminal flux fails the sealed duct conservation check")
        states[step] = row
        print(json.dumps({"state": row}, allow_nan=False), flush=True)

    for h in sorted(set(s for s in steps if s > 0)):
        low, center, high = states[-h], states[0], states[h]
        stiffness = (high["energy"] + low["energy"] - 2 * center["energy"]) / h**2
        work_stiffness = (high["work"] - low["work"]) / (2 * h)
        slope = (high["energy"] - low["energy"]) / (2 * h)
        old, = (s for s in previous["signed"] if s["half_width"] == h)
        print(json.dumps(dict(
            mesh=mesh, half_step=h, phase_bias=BIAS, open_side_current=old_center["current"],
            sealed_side_current=center["current"], sealed_side_terminal_current=center["current_by_cut"]["cap_layer"],
            current_ratio=center["current"]/old_center["current"],
            open_side_slope=old["virtual_slope"], sealed_side_slope=center["work"],
            sealed_energy_slope=slope, open_side_stiffness=old["energy_stiffness"],
            sealed_side_stiffness=stiffness, sealed_stiffness_from_work=work_stiffness,
        ), allow_nan=False), flush=True)
        if abs(slope-center["work"]) > .0005 or abs(stiffness-work_stiffness) > .0005:
            raise ArithmeticError("Two-sided work and stationary energy disagree")


if __name__ == "__main__":
    main()
