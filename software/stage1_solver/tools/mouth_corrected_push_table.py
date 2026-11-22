"""Apply the corrected momentum force to saved unit-speed fields in the original boxes."""

from itertools import groupby
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_carrier_free_center import momentum_force_rows, motion_source, wall_row
from mouth_carrier_lateral_response import prepare
from mouth_occupied_push_baths import digest, load_bath


ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_rigid_drag_table.json")


def main():
    archive = json.loads(ARCHIVE.read_text())
    old = {(p["mesh"], p["far"], p["omega"]): p for p in archive["table"]}
    assert len(old) == len(archive["runs"]) == 22
    seen = set()
    for prep, cases in groupby(sorted(archive["runs"].values(), key=lambda r: r["preparation"]),
                               key=lambda r: r["preparation"]):
        cases = list(cases)
        command = cases[0]["command"]
        flag = "--bath" if "--bath" in command else "--source"
        folder = Path(command[command.index(flag) + 1])
        dom, base = load_bath(folder) if flag == "--bath" else prepare(SimpleNamespace(source=folder))
        base_hash = digest(folder / "field.npz")
        static, temporal = momentum_force_rows(dom, base, motion_source(dom, base))
        direct, rate, boost = wall_row(dom, base)
        for run in cases:
            key = run["mesh"], run["far"], run["omega"]
            w = run["omega"]
            saved = Path(run["command"][-1])
            assert run["guard_exit"] == run["worker_exit"] == run["swap_cap_bytes"] == 0
            assert run["limits_verified"] and run["source_field_sha256"] == base_hash
            assert digest(saved) == run["response_field_sha256"] and key not in seen
            seen.add(key)
            with np.load(saved, allow_pickle=False) as field:
                np.testing.assert_array_equal(field["points"], dom.points)
                assert field["background_sha256"].item() == base_hash and field["omega"].item() == w
                h = np.r_[field["harmonic_real"], field["harmonic_imag"]]
                expected_direct = field["impedance"].item()
            z = (static + 1j*w*temporal) @ h
            old_z = -((direct + 1j*w*rate) @ h + boost)
            a_per_f = 1j*w/z
            port = old[key]["net_port_resistance"]
            windows = [old[key]["window1"], old[key]["window2"]]
            print(json.dumps(dict(test="corrected_push", prep=prep, mesh=key[0], far=key[1], omega=w,
                                  response_sha256=run["response_field_sha256"],
                                  trace_impedance=[old_z.real, old_z.imag],
                                  impedance=[z.real, z.imag], a_per_F=[a_per_f.real, a_per_f.imag],
                                  signed_push=w*(1/a_per_f).imag, twice_outward_power_per_v2=port,
                                  independent_windows=windows, minus_port=z.real-port,
                                  minus_windows=[z.real-x for x in windows]), allow_nan=False))
            np.testing.assert_allclose(old_z, expected_direct, rtol=1e-8, atol=2e-12)
            np.testing.assert_allclose(z.real, port, rtol=1e-8, atol=2e-12)
    assert seen == set(old)


if __name__ == "__main__":
    main()
