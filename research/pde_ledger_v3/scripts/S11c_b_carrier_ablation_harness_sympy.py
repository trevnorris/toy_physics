#!/usr/bin/env python3
"""S11c-b SymPy pressure-carrier ablation harness.

Manifest (the orchestrator-owned directive fixes these sites and forms):
K_A: build_operator, per-face closure_residuals before closure_residual_sum:
     expand each closure_shape_deriv and delete Lambda_A_0-bearing addends.
K_T: build_operator, final U_BODY_BALANCE and E_W_BALANCE construction:
     remove all four face_u / face_e additions.
K_W: substrate_substitutions, returned substitution dictionary:
     delta_p_minus -> delta_p_plus and d_w_delta_p_minus -> d_w_delta_p_plus.

1. The harness may **PRINT** computed objects (carriers, diffs). It may ⛔ NOT state conclusions — no `PASS`, no
   verdict, no "bites"/"fails". Interpretation is the orchestrator's, from the printed triples.
2. **Print operand and residual, then guard.** Emit `baseline`, `corrupted`, and their `diff`; a residual asserted
   zero/nonzero carries no information.
3. Interpretation belongs to the review / step record, ⛔ not the script.
"""

from __future__ import annotations

import argparse
import ast
from contextlib import redirect_stdout
from dataclasses import dataclass
import difflib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import pickle
import shlex
import subprocess
import sys
import tempfile
import traceback

import sympy as sp


SCRIPT = Path(__file__).resolve()
ENGINE = SCRIPT.with_name("S11c_b_brane_operator_sympy_audit.py")
DIRECTIVE = SCRIPT.parent.parent / "directives/S11c_b_carrier_ablation_harness_directive.md"
CASE = ("MATERIAL_ADVECTED", "RHO4_CONSTANT", "EULERIAN")
PRESSURE_NAMES = (
    "delta_p_plus", "d_w_delta_p_plus", "delta_p_minus", "d_w_delta_p_minus",
)
ROW_NAMES = (
    "U_BODY_BALANCE[1]", "U_BODY_BALANCE[2]", "U_BODY_BALANCE[3]",
    "E_W_BALANCE", "THETA_BALANCE",
)
SCHEMA = tuple((row, slot) for row in ROW_NAMES for slot in PRESSURE_NAMES)


@dataclass(frozen=True)
class Mutation:
    name: str
    function: str
    site: str
    form: str
    replacements: tuple[tuple[str, str, int], ...]


CLOSURE_BLOCK = '''    closure_residuals = tuple(
        sp.sympify(
            bind_mu_theta_operand(
                selected_substrate_axes(
                    faces,
                    "closure_shape_deriv",
                    (branch, face, "DELTA_W", representative),
                ),
                branch,
                mu_theta_amplitude,
            )
        )
        for face in FACES
    )
    closure_residual_sum = sp.Add(*closure_residuals)
'''
CLOSURE_FILTER = '''    Lambda_A_0 = next(atom for atom in DECLARED_SYMBOLS if str(atom) == "Lambda_A_0")
    closure_residuals = tuple(
        sp.Add(*(t for t in sp.Add.make_args(sp.expand(closure_shape_deriv)) if not t.has(Lambda_A_0)))
        for closure_shape_deriv in closure_residuals
    )
'''
CLOSURE_RESCALE = '''    Lambda_A_0 = next(atom for atom in DECLARED_SYMBOLS if str(atom) == "Lambda_A_0")
    closure_residuals = tuple(
        sp.Add(*(2*t if t.has(Lambda_A_0) else t for t in sp.Add.make_args(sp.expand(closure_shape_deriv))))
        for closure_shape_deriv in closure_residuals
    )
'''
COLLAPSE = '''    pressure_atoms = {str(atom): atom for atom in DECLARED_SYMBOLS}
    substitutions.update({
        pressure_atoms["delta_p_minus"]: pressure_atoms["delta_p_plus"],
        pressure_atoms["d_w_delta_p_minus"]: pressure_atoms["d_w_delta_p_plus"],
    })
    return substitutions
'''
MINUS_RESCALE = '''    pressure_atoms = {str(atom): atom for atom in DECLARED_SYMBOLS}
    substitutions[pressure_atoms["delta_p_minus"]] = 2 * pressure_atoms["delta_p_minus"]
    return substitutions
'''


def closure_replacement(insertion: str) -> str:
    return CLOSURE_BLOCK.replace(
        "    closure_residual_sum = sp.Add(*closure_residuals)\n",
        insertion + "    closure_residual_sum = sp.Add(*closure_residuals)\n",
    )


KNIVES = (
    Mutation(
        "K_A", "build_operator", "per-face closure_residuals before their response fold",
        "Expand every closure_shape_deriv; delete Lambda_A_0-bearing additive terms.",
        ((CLOSURE_BLOCK, closure_replacement(CLOSURE_FILTER), 1),),
    ),
    Mutation(
        "K_T", "build_operator", "four face_u / face_e additions in final balance rows",
        "Remove both face_u additions and both face_e additions.",
        (("                        + face_u[a]\n", "", 2),
         (" + face_e)", ")", 2)),
    ),
    Mutation(
        "K_W", "substrate_substitutions", "returned substitution dictionary",
        "Add delta_p_minus -> delta_p_plus and d_w_delta_p_minus -> d_w_delta_p_plus.",
        (("    return substitutions\n", COLLAPSE, 1),),
    ),
)
SELF_TESTS = (
    Mutation(
        "SELF_DEAD_PATH", "build_operator", "e_W kinetic addend at the reduced EXPANDED row",
        "Delete the -e_kinetic additive term.",
        (('            sp.expand(reduced_e["EXPANDED"] - e_kinetic),\n',
          '            sp.expand(reduced_e["EXPANDED"]),\n', 1),),
    ),
    Mutation(
        "SELF_RESCALE_K_A", "build_operator", "same per-face closure_residuals as K_A",
        "Expand every closure_shape_deriv; multiply Lambda_A_0-bearing addends by 2.",
        ((CLOSURE_BLOCK, closure_replacement(CLOSURE_RESCALE), 1),),
    ),
    Mutation(
        "SELF_RESCALE_K_T", "build_operator", "same four face_u / face_e additions as K_T",
        "Multiply each of the four face_u / face_e additions by 2.",
        (("                        + face_u[a]\n", "                        + 2 * face_u[a]\n", 2),
         (" + face_e)", " + 2 * face_e)", 2)),
    ),
    Mutation(
        "SELF_RESCALE_K_W", "substrate_substitutions", "same substitution dictionary as K_W",
        "Add delta_p_minus -> 2 * delta_p_minus, without face collapse.",
        (("    return substitutions\n", MINUS_RESCALE, 1),),
    ),
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def patch_construction(source: str, mutation: Mutation) -> tuple[str, dict]:
    """Replace exact fragments inside exactly the directive's named function."""
    nodes = [node for node in ast.parse(source).body
             if isinstance(node, ast.FunctionDef) and node.name == mutation.function]
    if len(nodes) != 1:
        raise RuntimeError(f"construction function mismatch: {mutation.function}")
    node = nodes[0]
    lines = source.splitlines(keepends=True)
    before = "".join(lines[node.lineno - 1:node.end_lineno])
    after = before
    for old, new, count in mutation.replacements:
        if after.count(old) != count:
            raise RuntimeError(f"construction site mismatch: {mutation.name}, {old!r}")
        after = after.replace(old, new)
    patched = "".join(lines[:node.lineno - 1]) + after + "".join(lines[node.end_lineno:])
    ast.parse(patched)
    return patched, {
        "name": mutation.name,
        "function": mutation.function,
        "canonical_function_lines": [node.lineno, node.end_lineno],
        "site": mutation.site,
        "form": mutation.form,
        "source_sha256": sha256(patched.encode()),
    }


def import_engine(path: Path):
    # Temporary engine copies resolve the canonical engine's unchanged incoming ledger.
    sys.path.insert(0, str(ENGINE.parent))
    spec = importlib.util.spec_from_file_location("s11cb_carrier_live_engine", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import engine: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def pressure_atoms(engine) -> tuple:
    slots = []
    for name in PRESSURE_NAMES:
        found = [atom for atom in engine.DECLARED_SYMBOLS if str(atom) == name]
        if len(found) != 1:
            raise RuntimeError(f"native atom table mismatch: {name}")
        slots.append(found[0])
    return tuple(slots)


def observed_rows(engine, operator) -> tuple:
    u_rows = engine.named_tuple_row(
        engine.named_tuple_row(operator, "U_BODY_BALANCE"), "EXPANDED"
    )
    if not isinstance(u_rows, sp.Tuple) or len(u_rows) != 3:
        raise RuntimeError("U_BODY_BALANCE EXPANDED row shape mismatch")
    return (
        *u_rows,
        engine.named_tuple_row(engine.named_tuple_row(operator, "E_W_BALANCE"), "EXPANDED"),
        engine.named_tuple_row(engine.named_tuple_row(operator, "THETA_BALANCE"), "EXPANDED"),
    )


def carrier(rows: tuple, slots: tuple, *, zero_first: bool = False) -> tuple:
    zero = dict.fromkeys(slots, sp.S.Zero)
    return tuple(
        (row_name, slot_name,
         sp.diff(expression.xreplace(zero), atom) if zero_first
         else sp.diff(expression, atom).xreplace(zero))
        for row_name, expression in zip(ROW_NAMES, rows, strict=True)
        for slot_name, atom in zip(PRESSURE_NAMES, slots, strict=True)
    )


def worker(engine_path: Path, output_path: Path) -> None:
    os.environ["S11CB_PROJECTION_WORKERS"] = "1"
    engine = import_engine(engine_path)
    operator, _, _ = engine.build_operator(*CASE)
    retained = engine.retained_grade(operator)
    rows = observed_rows(engine, retained)
    slots = pressure_atoms(engine)
    result = {
        "operator_sha256": sha256(engine.render(retained).encode()),
        "carrier": carrier(rows, slots),
        "zero_first_carrier": carrier(rows, slots, zero_first=True),
    }
    with output_path.open("wb") as handle:
        pickle.dump(result, handle, protocol=pickle.HIGHEST_PROTOCOL)


def run_engine(path: Path, output_path: Path, label: str) -> dict:
    print(f"BEGIN_ENGINE {label}", flush=True)
    environment = dict(os.environ, S11CB_PROJECTION_WORKERS="1", PYTHONHASHSEED="0")
    result = subprocess.run(
        [sys.executable, "-u", str(SCRIPT), "--worker-engine", str(path),
         "--worker-output", str(output_path)],
        env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n", flush=True)
    if result.returncode:
        raise RuntimeError(f"engine subprocess exited {result.returncode}: {label}")
    with output_path.open("rb") as handle:
        data = pickle.load(handle)
    print(f"END_ENGINE {label}", flush=True)
    return data


def emit_carrier(label: str, values: tuple) -> None:
    print(label)
    for row, slot, expression in values:
        print(json.dumps({"row": row, "slot": slot, "value": sp.srepr(expression)}, sort_keys=True))


def guard_schema(values: tuple) -> None:
    if tuple((row, slot) for row, slot, _ in values) != SCHEMA:
        raise RuntimeError("carrier row/slot schema mismatch")


def emit_triple(label: str, baseline: tuple, corrupted: tuple) -> None:
    print(f"BEGIN_TRIPLE {label}")
    emit_carrier("baseline", baseline)
    emit_carrier("corrupted", corrupted)
    diff = tuple(
        (row, slot, sp.expand(changed - original))
        for (row, slot, original), (_, _, changed) in zip(baseline, corrupted, strict=True)
    )
    emit_carrier("diff", diff)
    print(f"END_TRIPLE {label}", flush=True)
    guard_schema(baseline)
    guard_schema(corrupted)
    guard_schema(diff)


def run() -> None:
    source = ENGINE.read_text()
    source_digest = sha256(source.encode())
    print("S11CB_SYMPY_PRESSURE_CARRIER_HARNESS")
    print(json.dumps({
        "branch": CASE[0], "density": CASE[1], "route": CASE[2],
        "pressure_slot_order": PRESSURE_NAMES, "row_order": ROW_NAMES,
        "python": sys.version, "sympy": sp.__version__, "projection_workers": 1,
        "engine_source_sha256": source_digest,
        "harness_source_sha256": sha256(SCRIPT.read_bytes()),
        "directive_source_sha256": sha256(DIRECTIVE.read_bytes()),
        "incoming_ledger_sha256": sha256(ENGINE.with_name("S11c_a_exports.py").read_bytes()),
    }, sort_keys=True))
    patched_sources = []
    for mutation in (*KNIVES, *SELF_TESTS):
        patched, manifest = patch_construction(source, mutation)
        print("MANIFEST " + json.dumps(manifest, sort_keys=True))
        print(f"BEGIN_SOURCE_PATCH {mutation.name}")
        print("".join(difflib.unified_diff(
            source.splitlines(keepends=True), patched.splitlines(keepends=True),
            fromfile=ENGINE.name, tofile=f"{mutation.name}/{ENGINE.name}",
        )), end="")
        print(f"END_SOURCE_PATCH {mutation.name}")
        patched_sources.append((mutation, patched))
    print("EXTRACTOR_SELF_TEST " + json.dumps({
        "name": "SELF_EXTRACTOR_ORDER",
        "form": "Apply all four native pressure atoms -> 0 before differentiation.",
    }, sort_keys=True), flush=True)

    # subprocess.run is blocking: at most one production CAS construction is active.
    with tempfile.TemporaryDirectory(prefix="s11cb_sympy_carrier_") as directory:
        scratch = Path(directory)
        baseline = run_engine(ENGINE, scratch / "baseline.pkl", "CANONICAL")
        canonical_carrier = baseline["carrier"]
        unchanged = scratch / "unablated" / ENGINE.name
        unchanged.parent.mkdir()
        unchanged.write_text(source)
        copied = run_engine(unchanged, scratch / "unablated.pkl", "UNABLATED_COPY")
        emit_triple("DRIFT_UNABLATED_COPY", canonical_carrier, copied["carrier"])
        digest_comparison = {
            "baseline": baseline["operator_sha256"],
            "corrupted": copied["operator_sha256"],
            "equal": baseline["operator_sha256"] == copied["operator_sha256"],
        }
        print("RETAINED_OPERATOR_DIGEST_COMPARISON " + json.dumps(digest_comparison, sort_keys=True))
        print("CARRIER_COPY_COMPARISON " + json.dumps({
            "equal": canonical_carrier == copied["carrier"],
        }, sort_keys=True), flush=True)
        if not digest_comparison["equal"] or canonical_carrier != copied["carrier"]:
            raise RuntimeError("canonical/unablated-copy object drift")

        for mutation, patched in patched_sources[:len(KNIVES)]:
            path = scratch / mutation.name / ENGINE.name
            path.parent.mkdir()
            path.write_text(patched)
            result = run_engine(path, scratch / f"{mutation.name}.pkl", mutation.name)
            emit_triple(mutation.name, canonical_carrier, result["carrier"])

        emit_triple("SELF_EXTRACTOR_ORDER", canonical_carrier, baseline["zero_first_carrier"])

        for mutation, patched in patched_sources[len(KNIVES):]:
            path = scratch / mutation.name / ENGINE.name
            path.parent.mkdir()
            path.write_text(patched)
            result = run_engine(path, scratch / f"{mutation.name}.pkl", mutation.name)
            emit_triple(mutation.name, canonical_carrier, result["carrier"])

    final_digest = sha256(ENGINE.read_bytes())
    print("CANONICAL_SOURCE_COMPARISON " + json.dumps({
        "baseline": source_digest, "after": final_digest, "equal": source_digest == final_digest,
    }, sort_keys=True), flush=True)
    if source_digest != final_digest:
        raise RuntimeError("canonical engine source changed during the harness run")


class Tee:
    def __init__(self, stream, capture):
        self.stream = stream
        self.capture = capture

    def write(self, text):
        self.capture.write(text)
        return self.stream.write(text)

    def flush(self):
        self.stream.flush()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="write invocation, digests, and literal stdout as Markdown")
    parser.add_argument("--worker-engine", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_engine is not None:
        if args.worker_output is None or args.record is not None:
            parser.error("worker requires an output path and cannot record")
        worker(args.worker_engine, args.worker_output)
        return 0
    if args.worker_output is not None:
        parser.error("worker output requires an engine")
    capture = io.StringIO()
    status = 0
    with redirect_stdout(Tee(sys.stdout, capture)):
        try:
            run()
        except Exception:
            traceback.print_exc(file=sys.stdout)
            status = 1
    if args.record is not None:
        transcript = capture.getvalue()
        invocation = shlex.join([sys.executable, "-u", *sys.argv])
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(
            "# S11c-b SymPy carrier ablation harness transcript\n\n"
            f"Working directory: `{Path.cwd()}`\n\n"
            f"Invocation:\n\n```sh\n{invocation}\n```\n\n"
            f"Exit code: `{status}`\n\n"
            f"Harness source SHA-256: `{sha256(SCRIPT.read_bytes())}`\n\n"
            f"Engine source SHA-256: `{sha256(ENGINE.read_bytes())}`\n\n"
            f"Directive source SHA-256: `{sha256(DIRECTIVE.read_bytes())}`\n\n"
            f"Literal stdout SHA-256 (UTF-8, including final newline): `{sha256(transcript.encode())}`\n\n"
            f"Literal stdout:\n\n```text\n{transcript}```\n"
        )
    return status


if __name__ == "__main__":
    raise SystemExit(main())
