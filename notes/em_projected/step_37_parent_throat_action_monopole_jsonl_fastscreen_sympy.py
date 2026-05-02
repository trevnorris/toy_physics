#!/usr/bin/env python3
"""Self-test the single_throat_monopole JSONL fast-screen classifier."""
from __future__ import annotations

import hashlib
import json
import pathlib
import tempfile

from single_throat_monopole_jsonl_fastscreen import classify_events, load_events


def write_jsonl(path: pathlib.Path, events: list[dict]) -> None:
    path.write_text("".join(json.dumps(event) + "\n" for event in events), encoding="utf-8")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_monopole_jsonl_fastscreen",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "partial fast-screen for the existing single_throat_monopole JSON diagnostics",
    }

    init = {"event": "init", "EOS": {"n": 5}}
    good_diag = {
        "event": "diag",
        "mach_max": 0.18,
        "fits": {"dP_slope": -1.08, "dP_npts": 20, "geff_slope": -1.97, "geff_npts": 20},
    }
    bad_diag = {
        "event": "diag",
        "mach_max": 0.84,
        "fits": {"dP_slope": -0.32, "dP_npts": 20, "geff_slope": -1.11, "geff_npts": 20},
    }
    weak_diag = {
        "event": "diag",
        "mach_max": 0.11,
        "fits": {"dP_slope": -1.02, "dP_npts": 4, "geff_slope": -1.99, "geff_npts": 5},
    }
    boundary_diag = {
        "event": "diag",
        "mach_max": 0.6,
        "fits": {"dP_slope": -0.65, "dP_npts": 8, "geff_slope": -1.65, "geff_npts": 8},
    }
    just_outside_diag = {
        "event": "diag",
        "mach_max": 0.6001,
        "fits": {"dP_slope": -0.649, "dP_npts": 8, "geff_slope": -1.649, "geff_npts": 8},
    }
    done = {"event": "done", "t_final": 3.5}

    good_verdict = classify_events([init, good_diag, done])
    bad_verdict = classify_events([init, bad_diag, done])
    weak_verdict = classify_events([init, weak_diag, done])
    boundary_verdict = classify_events([init, boundary_diag, done])
    outside_verdict = classify_events([init, just_outside_diag, done])

    if good_verdict["status"] != "PASS":
        raise AssertionError(f"good monopole log should PASS, got {good_verdict}")
    if bad_verdict["status"] != "FAIL":
        raise AssertionError(f"bad monopole log should FAIL, got {bad_verdict}")
    if weak_verdict["status"] != "INCOMPLETE":
        raise AssertionError(f"weak monopole log should be INCOMPLETE, got {weak_verdict}")
    if boundary_verdict["status"] != "PASS":
        raise AssertionError(f"threshold-boundary monopole log should PASS, got {boundary_verdict}")
    if outside_verdict["status"] != "FAIL":
        raise AssertionError(f"just-outside-threshold monopole log should FAIL, got {outside_verdict}")

    with tempfile.TemporaryDirectory(prefix="em_projected_monopole_jsonl_") as tmp:
        log_path = pathlib.Path(tmp) / "monopole.log"
        write_jsonl(log_path, [init, good_diag, done])
        roundtrip = load_events(log_path)
        if len(roundtrip) != 3 or roundtrip[1]["fits"]["dP_slope"] != good_diag["fits"]["dP_slope"]:
            raise AssertionError("JSONL loader roundtrip failed")
        malformed_path = pathlib.Path(tmp) / "malformed.log"
        malformed_path.write_text(json.dumps(init) + "\nnot json\n" + json.dumps(good_diag) + "\n", encoding="utf-8")
        malformed_verdict = classify_events(load_events(malformed_path))
        if malformed_verdict["status"] != "INCOMPLETE":
            raise AssertionError(f"malformed JSONL should be INCOMPLETE, got {malformed_verdict}")
        if not any("malformed JSONL line" in warning for warning in malformed_verdict["warnings"]):
            raise AssertionError(f"malformed JSONL warning missing: {malformed_verdict}")

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "threshold_cases": {
            "good": good_verdict["status"],
            "bad": bad_verdict["status"],
            "weak": weak_verdict["status"],
            "boundary": boundary_verdict["status"],
            "outside": outside_verdict["status"],
            "malformed": malformed_verdict["status"],
        },
        "thresholds": good_verdict["thresholds"],
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    print("STEP 37 MONOPOLE JSONL FASTSCREEN SELF-TEST")
    print("Validated the partial fast-screen for single_throat_monopole.py JSON diagnostics.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Fast-screen thresholds:")
    print("  dP_slope target = -1 +/- 0.35")
    print("  geff_slope target = -2 +/- 0.35")
    print("  mach_max <= 0.6")
    print("  min_fit_points = 8")
    print("Classifier verdicts:")
    print("  good log =", good_verdict["status"])
    print("  bad log =", bad_verdict["status"])
    print("  weak log =", weak_verdict["status"])
    print("  threshold-boundary log =", boundary_verdict["status"])
    print("  just-outside-threshold log =", outside_verdict["status"])
    print("  malformed JSONL log =", malformed_verdict["status"])
    print("Representative failure reasons:")
    print("  bad log:", "; ".join(bad_verdict["failures"]))
    print("  weak log warnings:", "; ".join(weak_verdict["warnings"]))
    print("  malformed log warnings:", "; ".join(malformed_verdict["warnings"]))
    print("CLI:")
    print("  python single_throat_monopole_jsonl_fastscreen.py monopole.log --output-json monopole_verdict.json")
    print("Interpretation:")
    print("  This is a partial early-kill screen for the existing monopole solver output. It does not replace the full snapshot-based runtime monitor, but it lets us reject clearly wrong 1/r and 1/r^2 behavior immediately.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
