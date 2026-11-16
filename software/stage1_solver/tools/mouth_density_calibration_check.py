"""Inventory the existing radius-one preparations and check a conditional rest-density calibration."""

from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np
import sympy as sp

from held_crossflow_budget import archived_packet


def feed_packets():
    for speed, manifest in [(0, "mouth_feed_weak_run.json"),
                            (0.01, "mouth_intake_redistribution.json"),
                            (0.02, "mouth_intake_redistribution.json")]:
        path = archived_packet(manifest, f"n32_restoring_slope_u{speed:g}.npz")
        with np.load(path, allow_pickle=False) as packet:
            field = packet["field"]
            feed = abs(field[0])**2
            end = (len(field)-1)*float(packet["dz"])
            rate = float(packet["dx"])**3 * packet["current"].sum()
            print(f"feed patch: U={float(packet['speed']):g} feed density=[{feed.min():.12g},{feed.max():.12g}] "
                  f"last z={end:g} planes={len(field)} inward rate={rate:.12g}")
            print("  keys:", ", ".join(sorted(packet.files)))
            assert np.all(feed == 1) and end == 1 and float(packet["speed"]) == speed
            np.testing.assert_allclose(rate, .25, rtol=0, atol=1e-12)


def receiver_packets():
    record = Path("docs/_measurements/open_sleeve_current.txt").read_text()
    files = re.findall(r"^([a-f0-9]{64})  (\S*/level_32\.npz)$", record, re.M)
    expected_rates = {"open_sleeve_current01": 0.125, "open_sleeve_zero01": 0.0}
    preparations = [Path(name).parent.parent.name for _, name in files]
    if sorted(preparations) != sorted(expected_rates):
        raise ValueError("Expected the draining and zero-rate receiver preparations")
    for digest, name in files:
        path = Path(name)
        run = json.loads((path.parent.parent/"outcome.json").read_text())
        assert run["limitsVerified"] and run["exitCode"] == run["childOutcome"]["exitCode"] == 0
        assert sha256(path.read_bytes()).hexdigest() == digest
        with np.load(path, allow_pickle=False) as packet:
            width = float(packet["dw"])
            bath_density = float(packet["bath_density"])
            rate = float(packet["rate"])
            first_center = (packet["w_index"].min() + 0.5) * width
            receiver_join = int(packet["sleeve_cells"]) * width
            print(f"open sleeve: shared bath density={bath_density:g} "
                  f"prescribed rate={rate:g} first center w={first_center:g} "
                  f"receiver join w={receiver_join:g}")
            # This archived mesh spans 3.7 in 64 cells, with the sleeve ending halfway.
            expected = (1.0, expected_rates[path.parent.parent.name], 3.7 / 128, 1.85)
            np.testing.assert_allclose((bath_density, rate, first_center, receiver_join), expected,
                                       rtol=0, atol=1e-14, err_msg="Receiver preparation or geometry differs")


def finite_density_difference():
    # q is the local density difference from its reference at a fixed position.
    # The two smooth jumps coincide at the center and both sample points.
    q, h, slope, gamma, offset = sp.symbols("q h slope gamma offset", real=True)
    first = offset + slope*q
    second = first + gamma*q*(q-h)*(q+h)
    differences = [sp.factor((second-first).subs(q, x)) for x in (-h, 0, h)]
    slope_change = sp.factor(sp.diff(second-first, q).subs(q, 0))
    secant = sp.simplify((second.subs(q,h)-second.subs(q,-h))/(2*h))
    print("jump changes at local density offsets -h, 0, h:", differences)
    print("center derivative change:", slope_change, "centered secant:", secant)
    assert differences == [0, 0, 0] and slope_change == -gamma*h**2 and secant == slope


if __name__ == "__main__":
    feed_packets()
    receiver_packets()
    finite_density_difference()
