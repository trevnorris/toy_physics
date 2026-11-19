"""Inspect the paired drain records and test the impulse readout on the separate collar."""

import hashlib
import json
from pathlib import Path

import numpy as np


RUNS = Path("_scratch/s11c")
MEASUREMENTS = Path("docs/_measurements")
CASES = {"positive": .02, "negative": -.02, "positive_half": .01, "negative_half": -.01}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def carrier_pair():
    pins = set((MEASUREMENTS/"mouth_prescribed_acceleration.txt").read_text().splitlines())
    rows = {}
    for name, acceleration in CASES.items():
        folder = RUNS/f"mouth_prescribed_accel_{name}02"
        for path in (folder/"outcome.json", folder/"stdout"):
            if f"{digest(path)}  {path}" not in pins:
                raise ValueError(f"Unverified archived input {path}")
        outcome = json.loads((folder/"outcome.json").read_text())
        if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
            raise ValueError(f"Incomplete drive {name}")
        packet = json.loads((folder/"stdout").read_text())
        rows[name] = packet["32"]
        assert packet["acceleration"] == acceleration
        for row in rows[name]:
            np.testing.assert_allclose(row["inward_rate"], .25, atol=1e-13, rtol=0)
    print("carrier sample columns:", sorted(set.union(*(set(row) for r in rows.values() for row in r))))
    for suffix in ("", "_half"):
        pos, neg = rows["positive"+suffix], rows["negative"+suffix]
        assert [r["t"] for r in pos] == [r["t"] for r in neg]
        plus, minus = pos[-1], neg[-1]
        odd = {key: (plus[key]-minus[key])/2 for key in ("carried","stress","feed_load")}
        print(f"feeding face {suffix or 'full'}: t={plus['t']} intake={plus['inward_rate']} "
              f"carried={odd['carried']:.12g} stress={odd['stress']:.12g} "
              f"total_rate={odd['feed_load']:.12g}")
        np.testing.assert_allclose(odd["carried"]+odd["stress"], odd["feed_load"], rtol=1e-12, atol=0)


def benchmark_pair(cells, suffix, pins):
    samples = []
    base = f"benchmark_continued_{cells}"
    for sign in ("positive", "negative"):
        path = RUNS/"signed_mouth_collar01"/"fields"/f"{base}_{sign}{suffix}.npz"
        if digest(path) != pins[str(path)]:
            raise ValueError(f"Unverified benchmark {path}")
        with np.load(path, allow_pickle=False) as data:
            assert str(data["scope"]) == base and float(data["acceleration"]) == CASES[sign+suffix]
            assert data["times"][-1] == 1.
            z, cut = data["z"], float(data["cut_depth"])
            assert z[0] == 0 and np.count_nonzero(z == cut) == 1
            # Index 1 is the wave envelope. It receives wall-to-wave exchange,
            # has zero rim traction and exports momentum at the measuring cut.
            samples.append((z[z <= cut], data["velocity"][1, z <= cut, -1],
                            data["source_impulse"][1,-1], data["cut_impulse"][1,-1],
                            data["stored"][1,-1]))
    pos, neg = samples
    np.testing.assert_array_equal(pos[0], neg[0])
    v, source, exported, stored = [(a-b)/2 for a,b in zip(pos[1:], neg[1:], strict=True)]
    net = source - exported
    print(f"collar n={cells} {suffix or 'full'}: odd_source={source:.12g} "
          f"odd_export={exported:.12g} impulse_difference={net:.12g} odd_stored={stored:.12g} "
          f"velocity_range=[{v.min():.12g},{v.max():.12g}]")
    np.testing.assert_allclose(net, stored, rtol=2e-10, atol=0)
    assert net > 0 and np.min(v) > 0
    lower, upper = net/v.max(), net/v.min()
    # External check: integrate the archive generator's specified continuum
    # line mass cos^2(pi*z/(2*1.85)). It is never used to form the bounds.
    known = cut/2 + 1.85/(2*np.pi)*np.sin(np.pi*cut/1.85)
    print(f"collar n={cells} {suffix or 'full'}: integrated_mass_interval=[{lower:.12g},{upper:.12g}] "
          f"supplied_continuum_mass={known:.12g}")
    assert lower < known < upper
    return v, net, lower, upper


if __name__ == "__main__":
    carrier_pair()
    pins = json.loads((MEASUREMENTS/"signed_mouth_collar.json").read_text())["output_sha256"]
    for cells in (64,128):
        full = benchmark_pair(cells, "", pins)
        half = benchmark_pair(cells, "_half", pins)
        for i, (a, b) in enumerate(zip(full, half, strict=True)):
            np.testing.assert_allclose(b, a/2 if i < 2 else a, rtol=2e-7, atol=1e-13)
