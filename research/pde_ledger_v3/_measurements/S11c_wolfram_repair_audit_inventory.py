"""Inventory completed native audit evidence and optionally publish new outputs.

The Wolfram programs compute the mathematical operands. This file validates
their source pins, stream structure, grade bounds and machine-readable census.
It does not evaluate a second CAS engine or repair a physics producer.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]


def pin(path):
    return {"bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def inventory(run, name, publish):
    manifest = json.loads((run / "manifest.json").read_text())
    summary = json.loads((run / "summary.json").read_text())
    if manifest["status"] != "completed" or manifest["exitCode"] != 0:
        raise ValueError(f"incomplete run: {run}")
    if pin(run / "full.out") != manifest["output"] or (run / "stderr.txt").stat().st_size:
        raise ValueError(f"output integrity: {run}")
    for relative, expected in manifest["sources"].items():
        if pin(ROOT / relative) != expected or pin(run / "sources" / relative) != expected:
            raise ValueError(f"source pin changed: {relative}")
    for relative, expected in manifest["priorOutputs"].items():
        if pin(ROOT / relative) != expected:
            raise ValueError(f"historical output changed: {relative}")
    lines = (run / "full.out").read_text().splitlines()
    keys = []
    physical_records = 0
    for line in lines:
        match = re.match(r"([a-z][A-Za-z0-9]*) = (.*)$", line)
        if not match:
            raise ValueError(f"unframed output: {line[:100]}")
        key, payload = match.groups()
        if key in keys or '"dimensionLTM" ->' not in payload or '"epsilonOrder" ->' not in payload:
            raise ValueError(f"key or metadata: {key}")
        keys.append(key)
        if '"retained" ->' in payload:
            physical_records += 1
            retained = payload.split('"retained" ->', 1)[1].split(', "discarded" ->', 1)[0]
            if re.search(r"(?:etaBg|sigmaW)\^(?:[2-9]|[1-9][0-9])", retained):
                raise ValueError(f"out-of-rectangle retained power: {key}")
            if re.search(r"\b(?:Indeterminate|ComplexInfinity|Series|SeriesCoefficient)\b", retained):
                raise ValueError(f"uncomputed/nonfinite retained operand: {key}")
    output = ROOT / "mathematica/out" / f"S11c_wolfram_{name}_audit.out"
    if publish:
        if output.exists() or output.is_symlink():
            if pin(output) != manifest["output"]:
                raise ValueError(f"refusing to overwrite existing evidence: {output}")
        else:
            shutil.copyfile(run / "full.out", output)
        if pin(output) != manifest["output"]:
            raise ValueError("publication hash mismatch")
    return {"runDirectory": str(run), "manifest": manifest, "manifestPin": pin(run / "manifest.json"),
            "summary": summary, "summaryPin": pin(run / "summary.json"),
            "tagCount": len(keys), "physicalRecordCount": physical_records,
            "keyDigest": hashlib.sha256("\n".join(keys).encode()).hexdigest(),
            "publication": {"path": str(output.relative_to(ROOT)), **manifest["output"]}}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    mechanical = inventory(args.run_root / "mechanical-final", "mechanical_repair", args.publish)
    trace = inventory(args.run_root / "trace-final", "pressure_trace", args.publish)
    ms, ts = mechanical["summary"], trace["summary"]
    for prefix in ("wlRepairInertia", "wlRepairStoredProbe", "wlRepairFaceWork", "wlRepairFaceRow", "wlRepairFacePower"):
        if ms[prefix]["nonzeroResiduals"] or ms[prefix]["records"] != ms[prefix]["zeroResiduals"]:
            raise ValueError(f"mechanical residual census needs inspection: {prefix}")
    for key in ("c1JoinResiduals", "nativeResponseJoinResiduals", "independentAffineResiduals"):
        if ts[key] != ["0"] * 8:
            raise ValueError(f"trace control census needs inspection: {key}")
    if any(value == "0" for value in ts["witnessResiduals"]):
        raise ValueError("witness census changed; inspect before classifying")
    checkpoint = {
        "recordedUtc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "AUDITED_DEFECT_REQUIRES_USER_APPROVAL", "published": args.publish,
        "leftPrerequisiteCommit": "b271c71b4601adf3aaa1457a4ac04a44f2f57995",
        "mechanical": mechanical, "pressureTrace": trace,
        "inventorySource": pin(Path(__file__)),
        "issues": [
            {"id": "conservative_inertia_sign", "status": "ABSENT_ON_AUDITED_DOMAIN",
             "scope": "b native kinetic terms and assembly coefficients in raw/frozen/live constructors; both density representatives; both constrained stored-variation routes on the stated energy probe"},
            {"id": "mechanical_load_normalization", "status": "ABSENT_ON_AUDITED_DOMAIN",
             "scope": "b prescribed face work, all five generalized components and power, both faces/routes/anchorings/densities with symbolic stationary width; not full closed c2 power"},
            {"id": "reference_physical_pressure_trace", "status": "PRESENT",
             "scope": "c2 N6 native image assignments and affine trace on the exact constant-height diagonal locus; c1 pressure independently joined to the boundary solution",
             "earliestAffectedProducer": "mathematica/S11c_c2_N6_mathematica_audit.wl",
             "unexecutedCoverage": "nonconstant-profile ordered two/three-leg eta/sigma inverse; no absence or global spectral claim"},
            {"id": "thickness_kinetic_coordinate", "status": "ABSENT_ON_AUDITED_DOMAIN",
             "scope": "b native kinetic terms against the physical deltaW=W0*eW action with symbolic nonuniform width, both densities, all three constructors"}
        ],
        "nextAction": "Obtain user approval for S11c_wolfram_pressure_trace_repair_plan.md before editing the native c2 producer. REFERENCE/scattering remains paused at this decision.",
        "unchanged": ["native b/c1/c2 Wolfram producers and prior transcripts", "all SymPy producers/exports/accepted normalization outputs", "physics authorities", "S10/Lean", "retained solver/export contract"],
        "developmentRuns": {
            "mechanical-01": "Completed exploratory mechanical run; final adds computed summary and action-level negative-kinetic mutation.",
            "trace-01": "Completed counterexample; final uses explicit SeriesCoefficient selection to remove beyond-order rational remainders from retained metadata. Raw evidence preserved."
        }
    }
    path = ROOT / "_measurements/S11c_wolfram_repair_audit_checkpoint.json"
    path.write_text(json.dumps(checkpoint, indent=2) + "\n")
    print(json.dumps({"checkpoint": str(path), "mechanicalTags": mechanical["tagCount"],
                      "traceTags": trace["tagCount"], "published": args.publish}))


if __name__ == "__main__":
    main()
