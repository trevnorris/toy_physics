"""Measure finite-box decay of the signed energy currents in the frozen push data."""

import hashlib
import json
import math
from pathlib import Path


MEASUREMENTS = Path(__file__).resolve().parents[3] / "research/pde_ledger_v3/_measurements"


def interval_exponent(a, b, key):
    return math.log(abs(a[key] / b[key])) / math.log(b["far"] / a["far"])


def main():
    archive = json.loads((MEASUREMENTS / "mouth_rigid_drag_table.json").read_text())
    rows = [r for r in archive["table"] if r["mesh"] == 34]
    for w in sorted({r["omega"] for r in rows}):
        a, b, c = sorted((r for r in rows if r["omega"] == w), key=lambda r: r["far"])
        assert [r["far"] for r in (a, b, c)] == [9, 13.5, 18]
        assert all(r["feed_resistance"] > 0 > r["bulk_resistance"] for r in (a, b, c))
        result = {"test": "finite_box_scaling", "omega": w}
        for key, dim in [("feed_resistance", 3), ("bulk_resistance", 4)]:
            # The smallest box fixes the coefficient. The other boxes provide checks.
            result[key] = dict(
                measured=[r[key] for r in (a, b, c)],
                adjacent_exponents=[interval_exponent(a, b, key), interval_exponent(b, c, key)],
                endpoint_exponent=interval_exponent(a, c, key),
                normalized_r_to_d=[r[key] * r["far"]**dim for r in (a, b, c)],
                predicted_from_smallest=[a[key] * (a["far"]/r["far"])**dim for r in (b, c)],
                observed_over_prediction=[r[key]/(a[key]*(a["far"]/r["far"])**dim) for r in (b,c)])
        result.update(
            root_under_assumed_powers=[-r["far"]*r["bulk_resistance"]/r["feed_resistance"]
                                       for r in (a,b,c)],
            direct_over_net_largest=-c["omega_im_inverse_chi"]/c["net_port_resistance"])
        print(json.dumps(result, allow_nan=False))

    checks = json.loads((MEASUREMENTS / "mouth_rigid_drag_checks.json").read_text())["runs"]
    steady = checks["steady"]
    assert hashlib.sha256(steady["literal_stdout"].encode()).hexdigest() == steady["stdout_sha256"]
    dc = json.loads(steady["literal_stdout"].splitlines()[-1])
    low = min((r for r in rows if r["far"] == 18), key=lambda r: r["omega"])
    static = dc["FoverV_wall_via_domain"][0]
    print(json.dumps(dict(test="steady_comparison", omega=low["omega"],
                          steady_reservoir=static, harmonic_reservoir=low["net_port_resistance"],
                          relative_difference=low["net_port_resistance"]/static-1,
                          direct_wall_steady=dc["FoverV_direct_wall"][0],
                          direct_wall_harmonic=-low["omega_im_inverse_chi"])))


if __name__ == "__main__":
    main()
