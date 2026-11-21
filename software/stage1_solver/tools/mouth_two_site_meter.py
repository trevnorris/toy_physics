"""Read a matched two-pad experiment using only positions and gauge calibration."""
import argparse
import json
from decimal import Decimal as D, getcontext
from pathlib import Path


getcontext().prec = 55
SITE_IDS = ("wall", "ordered")
POSITION_TOL = D("1e-33")


def read_site(row, name):
    pad = row["pads"][name]
    position = [D(x) for x in pad["xy"]]
    stage = [D(x) for x in pad["stage_xy"]]
    stiffness = D(pad["gauge_k"])
    if stiffness <= 0:
        raise ValueError("Invalid gauge calibration")
    return position, [stiffness*(c-r) for c, r in zip(stage, position, strict=True)]


def paired(on, off, header):
    corrected = []
    forces = []
    mismatch = D(0)
    for row in (on, off):
        if set(row["pads"]) != set(SITE_IDS):
            raise ValueError("Missing or ambiguous material site")
        a, fa = read_site(row, "wall")
        b, fb = read_site(row, "ordered")
        gap = D(row["pads"]["ordered"]["z"]) - D(row["pads"]["wall"]["z"])
        targets = [D(row["q_target"]), D(row["y_target"])]
        for target, left, right in zip(targets, a, b, strict=True):
            mismatch = max(mismatch, abs(right-left-target), abs((left+right)/2))
        if gap != D(header["normal_gap"]) or mismatch > POSITION_TOL:
            raise ValueError("Measured midpoint, offset, sideways step or gap missed")
        forces.append((fa, fb))
    for name in SITE_IDS:
        p, c = on["pads"][name], off["pads"][name]
        if (p["gauge_k"], p["z"]) != (c["gauge_k"], c["z"]):
            raise ValueError("Unmatched force sensors or normal locations")
        mismatch = max(mismatch, *(abs(D(x)-D(y)) for x, y in zip(p["xy"], c["xy"], strict=True)))
    if mismatch > POSITION_TOL:
        raise ValueError("Control and connected pad coordinates differ")
    if len(on["strand_lengths"]) != 2 or off["strand_lengths"]:
        raise ValueError("Invalid attachment/control pairing")
    if min(map(D, on["strand_lengths"])) <= D(header["rest_length"]):
        raise ValueError("Trial attachment lost tensile preload")
    for f_on, f_off in zip(*forces, strict=True):
        corrected.append([u-v for u, v in zip(f_on, f_off, strict=True)])
    imbalance = max(abs(u+v) for u, v in zip(*corrected, strict=True))
    if imbalance > D("1e-31"):
        raise ValueError("Residual external force after subtracting controls")
    h_on, h_off = [(fb[0]-fa[0])/2 for fa, fb in forces]
    return h_on, h_off, h_on-h_off, mismatch, imbalance


def main(path):
    header, *rows = [json.loads(line) for line in path.read_text().splitlines()]
    q = D(header["registry"])
    full = D("0.00001")
    targets = {(q, y) for y in (D(0), -full, full, -full/2, full/2)}
    targets.update({(q-full, D(0)), (q+full, D(0))})
    expected = {(width, x, y, state) for width in (D(".25"), D(1))
                for x, y in targets for state in ("connected", "brane_only")}
    indexed = {(D(r["span"]), D(r["q_target"]), D(r["y_target"]), r["fixture"]): r
               for r in rows}
    if expected != indexed.keys() or len(rows) != len(expected):
        raise ValueError("Incomplete or duplicated experiment")
    for width in (D(".25"), D(1)):
        data = {}
        for x, y in sorted(targets):
            pair = [indexed[width, x, y, state] for state in ("connected", "brane_only")]
            data[x, y] = paired(*pair, header)
        slope = (data[q+full, 0][2]-data[q-full, 0][2])/(2*full)
        for h in (full, full/2):
            trio = [data[q, y] for y in (-h, 0, h)]
            on, off, net = [[r[i] for r in trio] for i in (0, 1, 2)]
            even = (net[0]+net[2])/2 - net[1]
            print(json.dumps({"span": width, "step": h, "step_order": ["-", "0", "+"],
                              "H_connected": on, "H_brane": off, "H_attachment": net,
                              "delta_even": even, "delta_odd": (net[2]-net[0])/2,
                              "brane_even": (off[0]+off[2])/2-off[1],
                              "registry_slope_centered": slope,
                              "equivalent_registry_drift": abs(even/slope) if slope else None,
                              "max_position_discrepancy": max(r[3] for r in data.values()),
                              "max_two_pad_force_imbalance": max(r[4] for r in data.values())},
                             default=str))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("readings", type=Path)
    main(parser.parse_args().readings)
