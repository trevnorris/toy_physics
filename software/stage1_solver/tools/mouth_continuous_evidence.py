"""Validate the terminal record for the saved radius-four held order fields."""
from pathlib import Path

from mouth_rim_pinned_comparison import digest, read, validate_guard


def held_receipt(folder):
    folder = Path(folder)
    native = folder/"outcome.json"
    if native.is_file():
        validate_guard(folder)
        return native
    path = folder/"terminal_worker_verification.json"
    adoption = read(path)
    inv = read(folder/"invocation.json")
    child = read(folder/"child-outcome.json")
    limits = read(folder/"limit-validation.json")
    duration = read(folder/"duration-validation.json")
    hashes = adoption["evidence_sha256"]
    required = {"invocation.json", "child-outcome.json", "limit-validation.json",
                "duration-validation.json", "stdout", "../saved/record.json",
                "../saved/rest.npz", "../saved/wound.npz"}
    if not (folder.resolve() == Path(inv["logDirectory"]).resolve()
            and inv["unit"] == "s11c-guard-3b3b69784cb8"
            and adoption["unit"] == inv["unit"]+".service"
            and adoption["adopted_authorized_worker"] is True
            and adoption["guard_exit"] is None
            and adoption["limits_verified"] is True
            and child == adoption["worker_outcome"] and child["exitCode"] == 0
            and child["guardReason"] is None and limits["verified"] and duration["verified"]
            and limits["actual"]["memory.swap.max"] == "0"
            and duration["wallDeadlineSeconds"] is None and duration["inactivityDeadlineSeconds"] is None
            and required.issubset(hashes)
            and all(digest(folder/name) == sha for name, sha in hashes.items())):
        raise ValueError("Unverified terminal worker or mismatched radius-four held fields")
    return path
