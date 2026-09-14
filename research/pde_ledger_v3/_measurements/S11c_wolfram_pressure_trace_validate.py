"""Validate source-pinned focused/native Wolfram pressure-repair evidence.

Only inert text/JSON/sparse-array data are read here. The four balanced text
parser functions are reused from the existing native harness without importing
or executing its run driver, SymPy engine, review legs or comparator.
"""
from __future__ import annotations

import argparse
import ast
from collections import OrderedDict
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]


def pin(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            result.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": result.hexdigest()}


def parser_functions():
    source = ROOT / "scripts/S11c_c2_N6_ablation_harness_wl.py"
    names = {"split_top", "association", "literal", "sparse_row_counts"}
    module = ast.parse(source.read_text())
    selected = [node for node in module.body if isinstance(node, ast.FunctionDef) and node.name in names]
    if {node.name for node in selected} != names:
        raise ValueError("native text parser census")
    namespace = {"ast": ast, "OrderedDict": OrderedDict}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(source), "exec"), namespace)
    return namespace, {"path": str(source.relative_to(ROOT)), **pin(source)}


def check_manifest(run):
    manifest = json.loads((run / "manifest.json").read_text())
    if manifest["status"] != "completed" or manifest["exitCode"] != 0:
        raise ValueError("run incomplete or failed")
    if (run / "stderr.txt").stat().st_size or pin(run / "full.out") != manifest["output"]:
        raise ValueError("stderr/output integrity")
    if not all(manifest["sourceStability"].values()) or not all(manifest["priorOutputStability"].values()):
        raise ValueError("run changed a frozen source or prior output")
    for relative, expected in manifest["sources"].items():
        if pin(ROOT / relative) != expected or pin(run / "sources" / relative) != expected:
            raise ValueError(f"source integrity: {relative}")
    if manifest.get("sourceTransform"):
        if pin(Path(manifest["command"][-1])) != manifest["sourceTransform"]["executedPin"]:
            raise ValueError("executed case restriction changed")
    return manifest


def focused(run):
    summary = json.loads((run / "summary.json").read_text())
    values = [v for row in summary["residuals"].values() for v in row]
    result = {"summary": summary, "summaryPin": pin(run / "summary.json"),
              "residualScalarCount": len(values), "nonzeroResidualScalars": sum(v != "0" for v in values)}
    if result["nonzeroResidualScalars"]:
        raise ValueError("focused residual requires inspection")
    if not summary["mutationNonzeroComponents"] or any(n == 0 for n in summary["mutationNonzeroComponents"].values()):
        raise ValueError("missing mutation sensitivity")
    if any(n == 0 for n in summary["orderedMomentumNonzeroComponents"]):
        raise ValueError("missing intermediate-momentum sensitivity")
    return result


def native(run):
    parser, source = parser_functions()
    association, literal, sparse_counts = (parser[k] for k in ("association", "literal", "sparse_row_counts"))
    # The adopted N6 disposition is complete-operator covariance (Reading B),
    # not equality of coefficients in different field-variable frames. Keep
    # the raw R_N6 diagnostic and check its independently built channel split.
    residual_tags = {
        "WL_S11CC2_N6RC_SPLIT_CHECK", "WL_S11CC2_N6COV_R_COV", "WL_S11CC2_N6COV_R_COV_INCREMENT",
        "WL_S11CC2_N6RC_CARRIER_BRIDGE_RESIDUAL", "WL_S11CC2_N6RC_CARRIER_CHANNEL", "WL_S11CC2_N6RC_CROSS_CHANNEL",
        "WL_S11CC2_N6_SLOT_GUARD_RESIDUAL", "WL_S11CC2_N6_CLOSURE_GUARD_RESIDUAL",
        "wlS11cc2PressureTraceReconstructionResidual",
    }
    diagnostic_tags = {"WL_S11CC2_N6RC_R_N6", "WL_S11CC2_N6RC_SOURCE_CHANNEL"}
    split_top = parser["split_top"]
    cases, tags, bounds, checks, hashes, trace_dimensions, numeric_checks = {}, [], {}, {}, {}, {}, {}
    trace_tags = {"wlS11cc2PressureTrace" + name for name in
                  ("PhysicalResponse", "ReferenceResponse", "ReconstructedPressure", "ReconstructionResidual")}
    def denominator_shape(text):
        if not text.startswith("SparseArray[") or not text.endswith("]"):
            raise ValueError("denominator sparse head")
        args = list(split_top(text[len("SparseArray["):-1]))
        if len(args) != 4 or args[0].strip() != "Automatic":
            raise ValueError("denominator sparse storage")
        shape, default, storage = map(literal, args[1:])
        if len(shape) != 2 or default != 1 or storage[0] != 1:
            raise ValueError("denominator dimensions/default")
        pointers, columns = storage[1]
        values = storage[2]
        if len(pointers) != shape[0] + 1 or pointers[-1] != len(values) or len(columns) != len(values) or any(v == 0 for v in values):
            raise ValueError("denominator sparse index or zero at a retained sample")
        return shape
    required = residual_tags | diagnostic_tags | {"WL_S11CC2_N6_LOCAL_PROBE", "WL_S11CC2_N6_LOCAL_INVENTORY",
        "WL_S11CC2_N6RC_DIMENSIONS", "wlS11cc2PressureTracePhysicalResponse",
        "wlS11cc2PressureTraceReferenceResponse", "wlS11cc2PressureTraceReconstructedPressure"}
    # Native emissions are one complete assignment per line. Individual sparse
    # numerator arrays are decoded; arithmetic circuits remain inert strings.
    with (run / "full.out").open() as stream:
        for line in stream:
            tag, separator, payload = line.rstrip("\n").partition(" = ")
            if not separator or tag in tags:
                raise ValueError(f"native stream framing/duplicate: {tag[:80]}")
            tags.append(tag)
            hashes[tag] = hashlib.sha256(line.encode()).hexdigest()
            if tag in residual_tags | trace_tags | diagnostic_tags:
                case_records = {}
                for case, entries in association(payload).items():
                    records = association(entries)
                    scalars = nonzero = 0
                    sample_rows = set()
                    trace_keys, dimension_components = set(), 0
                    for key, component in records.items():
                        fields = association(component)
                        count, shape = sparse_counts(fields['"PROBE_NUMERATORS"'])
                        if denominator_shape(fields['"PROBE_DENOMINATORS"']) != shape:
                            raise ValueError("numerator/denominator sample shape join")
                        scalars += shape[0] * shape[1]
                        nonzero += sum(count)
                        sample_rows.add(shape[0])
                        if tag in trace_tags:
                            family, face, grade = literal(key)
                            trace_keys.add((family, face, tuple(grade)))
                            dimension_text = fields['"DIMENSIONS"'].strip()
                            dimensions = list(split_top(dimension_text[1:-1]))
                            if len(dimensions) != shape[1]:
                                raise ValueError("trace dimension/component census")
                            for dimension in dimensions:
                                d = association(dimension)
                                coefficient = literal(d['"COEFFICIENT_SUPPORT"'])
                                measure = literal(d['"MEASURE_SUPPORT"'])
                                pairing = literal(d['"TRIAL_TEST_SUPPORT"'])
                                restored = literal(d['"PAIRED_SUPPORT"'])
                                if len(coefficient) > 1 or len(measure) != 1 or any(len(u) != 3 for u in coefficient + measure + [pairing]):
                                    raise ValueError("trace unit support is unresolved/inconsistent")
                                if restored != [[u[i] + measure[0][i] + pairing[i] for i in range(3)] for u in coefficient]:
                                    raise ValueError("trace restored dimension join")
                                dimension_components += 1
                    if tag in trace_tags:
                        expected = {(family, face, (0, eta, sigma))
                                    for family in ("FOURIER_KOUT_Y", "FOURIER_KOUT_KIN_Y", "FOURIER_KOUT_KIN_Y_MIDDLE")
                                    for face in (-1, 1, "SUM") for eta in (0, 1) for sigma in (0, 1)}
                        if trace_keys != expected:
                            raise ValueError("trace family/face/grade census")
                        trace_dimensions.setdefault(tag, {})[case] = {"keys": len(trace_keys), "components": dimension_components,
                            "restoredDimensionJoinNonzero": 0, "grades": [[0, e, s] for e in (0, 1) for s in (0, 1)]}
                    if len(sample_rows) != 1:
                        raise ValueError(f"inconsistent joint sample count: {tag} {case}")
                    case_records[case] = {"entries": len(records), "sampleRows": next(iter(sample_rows)),
                                          "numeratorScalars": scalars, "nonzeroNumerators": nonzero}
                if tag in residual_tags:
                    checks[tag] = case_records
                numeric_checks[tag] = case_records
            elif tag == "WL_S11CC2_N6_LOCAL_PROBE":
                for case, value in association(payload).items():
                    fields = association(value)
                    index = literal(fields['"SAMPLE_INDEX"'])
                    cases[case] = {"sampleRows": len(index), "cells": len({row[0] for row in index}),
                                   "primes": sorted({row[1] for row in index}),
                                   "familyCardinality": int(fields['"FAMILY_CARDINALITY"'])}
                    bounds[case] = fields['"CONDITIONAL_UNION_BOUND"']
    if required - set(tags):
        raise ValueError(f"missing native tags: {required - set(tags)}")
    nonzero = 0
    for tag, entries in numeric_checks.items():
        if set(entries) != set(cases):
            raise ValueError(f"case census: {tag}")
        for case, record in entries.items():
            if record["sampleRows"] != cases[case]["sampleRows"]:
                raise ValueError(f"joint sample join: {tag} {case}")
            if tag in residual_tags:
                nonzero += record["nonzeroNumerators"]
    result = {"parserSource": source, "tagCount": len(tags), "tags": tags, "cases": cases,
              "residuals": checks, "nonzeroCertifiedNumerators": nonzero,
              "rawRepresentationDiagnostics": {tag: numeric_checks[tag] for tag in sorted(diagnostic_tags)},
              "interpretationSources": {name: pin(ROOT / "_measurements" / name) for name in
                  ("S11c_c2_N6_RESOLVED.md", "S11c_c2_N6_reconcile_disposition_PATH_B.md")},
              "conditionalUnionBounds": bounds,
              "traceDimensions": trace_dimensions, "emissionSha256": hashes,
              "sampleDenominatorZeroScalars": 0,
              "scope": "native recorded real-axis cells and joint samples, excluding stated singular/diagonal loci; not global spectral coverage"}
    if nonzero:
        (run / "native-unresolved-inventory.json").write_text(json.dumps(result, indent=2) + "\n")
        raise ValueError("native residual requires inspection")
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--run-directory", type=Path, required=True)
    parser.add_argument("--kind", choices=("focused", "native"), required=True)
    parser.add_argument("--expected-cases", type=int, required=True)
    parser.add_argument("--publish", type=Path)
    args = parser.parse_args()
    run = args.run_directory.resolve()
    manifest = check_manifest(run)
    result = focused(run) if args.kind == "focused" else native(run)
    actual_cases = result["summary"]["nativeRowCases"] // 4 if args.kind == "focused" else len(result["cases"])
    if actual_cases != args.expected_cases:
        raise ValueError("requested case census")
    publication = None
    if args.publish:
        destination = (ROOT / args.publish).absolute()
        destination.relative_to(ROOT / "mathematica/out")
        if destination.suffix != ".out":
            raise ValueError("publication must be a Mathematica .out")
        # Replacing a symlink leaves its historical annex object untouched.
        staging = run / "publication.out"
        shutil.copyfile(run / "full.out", staging)
        if pin(staging) != manifest["output"]:
            raise ValueError("staged publication hash")
        staging.replace(destination)
        publication = {"path": str(destination.relative_to(ROOT)), **pin(destination)}
    checkpoint = {"kind": args.kind, "runDirectory": str(run), "manifest": manifest,
                  "manifestPin": pin(run / "manifest.json"), "validationSource": pin(Path(__file__)),
                  "result": result, "publication": publication}
    name = f"S11c_wolfram_pressure_trace_{args.kind}_checkpoint.json"
    if args.expected_cases != 4:
        name = name.replace("_checkpoint", "_development_checkpoint")
    (ROOT / "_measurements" / name).write_text(json.dumps(checkpoint, indent=2) + "\n")
    (run / "validation.json").write_text(json.dumps(checkpoint, indent=2) + "\n")
    print(json.dumps({"checkpoint": name, "kind": args.kind, "cases": actual_cases,
                      "publication": publication}))


if __name__ == "__main__":
    main()
