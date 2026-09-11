#!/usr/bin/env python3
"""Compare the focused anisotropic stratum reruns, allowing different witnesses.

The coverage contract (parallel and perpendicular, D=3,4) comes from the Lean
direction classification. Engine root ordering and stratum numbering are not
assumed to agree. Compare root/|k|² and integer counts, not matrices at different
wavevectors. Check each engine's matrices and bases against its own witness.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import sympy as sp

from S10_cross_engine_comparator import (
    load_pair, normalize, parse_sympy_payload, parse_wolfram_payload,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def expression_equal(left, right) -> bool:
    return sp.cancel(left - right) == 0


def matrix_zero(matrix: sp.Matrix) -> bool:
    return all(expression_equal(x, 0) for x in matrix)


def wolfram_list_length(raw: str) -> int:
    """Count outer list elements without interpreting Reduce's Boolean syntax."""
    require(raw.startswith("{") and raw.endswith("}"), "expected a Wolfram list")
    body = raw[1:-1].strip()
    if not body:
        return 0
    # Association delimiters participate in nesting, just like brackets.
    tokens = re.findall(r'"(?:\\.|[^"\\])*"|<\||\|>|[{}\[\](),]|[^{}\[\](),"<>]+|.', body)
    stack, count = [], 1
    pairs = {"}": "{", "]": "[", ")": "(", "|>": "<|"}
    for token in tokens:
        if token in pairs.values():
            stack.append(token)
        elif token in pairs:
            require(bool(stack) and stack.pop() == pairs[token], "unbalanced Wolfram list")
        elif token == "," and not stack:
            count += 1
    require(not stack, "unbalanced Wolfram list")
    return count


def read_engine(transcript) -> dict:
    parse = parse_sympy_payload if transcript.engine == "PY" else parse_wolfram_payload
    require(not transcript.duplicates, f"{transcript.engine}: duplicate tags")
    allowed_message = 'no tagged-row grammar: Solve::svars: Equations may not give solutions for all "solve" variables.'
    require(not transcript.format_issues or (transcript.engine == "WL" and
            all(issue.endswith(allowed_message) for issue in transcript.format_issues)),
            f"{transcript.engine}: unexpected untagged output")

    def get(name):
        require(name in transcript.records, f"{transcript.engine}: missing {name}")
        return normalize(parse(transcript.records[name].raw))

    run_pairs = get("S10_RUN_PAIRS")
    require(len(run_pairs) == 2, f"{transcript.engine}: expected both anisotropic dimensions")
    require(not get("S10_SKIPPED_PAIRS"), f"{transcript.engine}: skipped run pairs")
    result = {}
    for dimension in (3, 4):
        base = f"S10_XFORM_ANISO_D{dimension}_"
        allowed_name = base + "Q8_ALLOWED_STRATA"
        require(allowed_name in transcript.records, f"{base}: missing stratum inventory")
        allowed_count = (len(get(allowed_name)) if transcript.engine == "PY" else
                         wolfram_list_length(transcript.records[allowed_name].raw))
        pattern = re.compile(re.escape(base) + (
            r"Q8_STRATUM([0-9]+)_POINT" if transcript.engine == "PY"
            else r"STRATUM([0-9]+)_Q8_POINT"))
        point_names = [(name, pattern.fullmatch(name)) for name in transcript.records]
        point_names = [(name, match) for name, match in point_names if match]
        require(len(point_names) == allowed_count, f"{base}: a discovered stratum lacks a rerun")
        for name, match in point_names:
            index = int(match.group(1))
            prefix = base + ("Q8_" if transcript.engine == "PY" else "") + f"STRATUM{index}_"
            substitutions = {str(rule.args[0]): rule.args[1] for rule in get(name)}
            k = sp.Matrix([substitutions[f"k{i}"] for i in range(1, dimension + 1)])
            require(all(x.is_Rational for x in k), f"{prefix}: witness is not exact rational")
            norm = (k.T * k)[0]
            require(norm > 0, f"{prefix}: zero wavevector witness")
            parallel = all(x == 0 for x in k[1:])
            perpendicular = k[0] == 0
            require(parallel != perpendicular, f"{prefix}: witness is outside the two proved strata")
            direction = "parallel" if parallel else "perpendicular"
            key = f"D{dimension}_{direction}"
            require(key not in result, f"{prefix}: two reruns cover the same direction")
            roots = get(prefix + "ROOT_ORDERING")
            require(get(prefix + "Q3_ROOT_COUNT") == len(roots), f"{prefix}: root-count mismatch")
            determinant = get(prefix + "Q3_DETERMINANT")
            spectral = sp.Symbol("omegaSquared")
            solved = sp.solve(determinant, spectral)
            require(len(solved) == len(roots) and
                    all(any(expression_equal(r, s) for s in roots) for r in solved),
                    f"{prefix}: emitted roots do not exhaust the emitted determinant")
            modes = []
            for root_index, root in enumerate(roots, 1):
                rp = prefix + f"ROOT{root_index}_"
                matrix = sp.Matrix(get(rp + "N1_MATRIX"))
                stacked = sp.Matrix(get(rp + "N3_STACKED_MATRIX"))
                require(matrix.shape == (dimension, dimension), f"{rp}: matrix shape")
                require(matrix_zero(stacked - matrix.col_join(k.T)), f"{rp}: incorrect N3 stack")
                n2 = int(get(rp + "N2_NULLITY"))
                n3 = int(get(rp + "N3_TRANSVERSE_NULLITY"))
                rank = int(get(rp + "N2_RANK"))
                stacked_rank = int(get(rp + "N3_STACKED_RANK"))
                require(rank == matrix.rank() and n2 == dimension - rank, f"{rp}: N2 mismatch")
                require(stacked_rank == stacked.rank() and n3 == dimension - stacked_rank,
                        f"{rp}: N3 mismatch")
                require(get(rp + "N4_NULLITY_DIFFERENCE") == n2 - n3, f"{rp}: N4 mismatch")
                vectors = [sp.Matrix(v) for v in get(rp + "N6_NULLSPACE_BASIS")]
                basis = sp.Matrix.hstack(*vectors) if vectors else sp.zeros(dimension, 0)
                require(matrix_zero(matrix * basis) and basis.rank() == n2, f"{rp}: basis mismatch")
                require(get(rp + "N7_BASIS_COUNT") == len(vectors) == n2, f"{rp}: N7 mismatch")
                modes.append({"root": str(root), "root_over_norm": str(sp.cancel(root / norm)),
                              "N2": n2, "N3": n3, "N4": n2 - n3,
                              "N2_rank": rank, "N3_rank": stacked_rank, "N7": len(vectors)})
            result[key] = {"point": [str(x) for x in k], "norm_squared": str(norm),
                           "prefix": prefix, "modes": modes}
    require(set(result) == {f"D{d}_{direction}" for d in (3, 4)
                            for direction in ("parallel", "perpendicular")},
            f"{transcript.engine}: incomplete Lean direction coverage")
    return result


def compare(paths: list[Path]) -> dict:
    py, wl = load_pair(paths)
    left, right = read_engine(py), read_engine(wl)
    joined = []
    for key in sorted(left):
        lm, rm = left[key]["modes"], right[key]["modes"]
        require(len(lm) == len(rm), f"{key}: distinct-root counts differ")
        for mode in lm:
            matches = [other for other in rm if expression_equal(
                sp.sympify(mode["root_over_norm"]), sp.sympify(other["root_over_norm"]))]
            require(len(matches) == 1, f"{key}: no unique normalized-root match")
            other = matches[0]
            require(all(mode[field] == other[field] for field in
                        ("N2", "N3", "N4", "N2_rank", "N3_rank", "N7")),
                    f"{key}: cross-engine count mismatch at {mode['root_over_norm']}")
            joined.append({"case": key, "root_over_norm": mode["root_over_norm"],
                           "N2": mode["N2"], "N3": mode["N3"]})
    return {"status": "PASS", "scope": "sampled strata: normalized roots and N2/N3/N4/N7 counts",
            "inputs": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            "diagnostics": {"WL_Solve_svars_messages": len(wl.format_issues)},
            "PY": left, "WL": right, "joined": joined}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_files", nargs=2, type=Path)
    args = parser.parse_args()
    try:
        result = compare(args.out_files)
    except (ValueError, KeyError, TypeError, AttributeError) as error:
        print(json.dumps({"status": "FAIL", "reason": str(error)}, indent=2))
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
