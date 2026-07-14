#!/usr/bin/env python3
"""Out-of-process mutation and metamorphic harness for Phase B1.

Every case edits a temporary copy and invokes the unchanged production
executable.  A case is accepted only when its guarded artifact object changes,
the process exits 1, and the expected non-sentinel assertion fires.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_FIXTURES = HERE / "u1_body_mechanics_fixtures.yaml"
DEFAULT_INPUT = HERE / "u1_body_mechanics_inputs.yaml"
DEFAULT_ARTIFACTS = HERE / "reports/u1_body_dynamics_artifacts"
RUNTIME = HERE / "_scratch/u1_b1_mutations"
ASSERTION = re.compile(r"^ASSERT_FAIL:([^:\s]+):", re.MULTILINE)
# Wolfram launches are serialized: concurrent activation handshakes can race and
# return 255 before either .wl file executes.  One seat is below the two-seat
# ceiling and prevents a license race from masquerading as a failed assertion.
WOLFRAM_SEATS = threading.Semaphore(1)
WOLFRAM_COLD_START = threading.Lock()
WOLFRAM_READY = threading.Event()


def load(path: Path) -> dict[str, Any]:
    text = path.read_text()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return yaml.safe_load(text)


def sha(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def at(value: Any, path: list[Any]) -> Any:
    cur = value
    for key in path:
        cur = cur[int(key)] if isinstance(cur, list) else cur[key]
    return cur


def dotted_at(value: Any, path: str) -> Any:
    return at(value, path.split("."))


def leaf_path_value(value: Any, path: str) -> Any:
    return dotted_at(value, path)


def parent(value: Any, path: list[Any]) -> tuple[Any, Any]:
    return at(value, path[:-1]), path[-1]


def mutate(value: dict[str, Any], operation: dict[str, Any]) -> None:
    op = operation["op"]
    if op == "set":
        host, key = parent(value, operation["path"])
        if isinstance(host, list): host[int(key)] = operation["value"]
        else: host[key] = operation["value"]
    elif op == "add":
        host, key = parent(value, operation["path"])
        if isinstance(host, list): host[int(key)] += operation["value"]
        else: host[key] += operation["value"]
    elif op == "delete_value":
        at(value, operation["path"]).remove(operation["value"])
    elif op == "delete_by_id":
        rows = at(value, operation["path"])
        rows[:] = [row for row in rows if row.get("id") != operation["id"]]
    elif op == "set_by_id":
        rows = at(value, operation["path"])
        row = next(row for row in rows if row.get("id") == operation["id"])
        row[operation["field"]] = operation["value"]
    elif op == "copy":
        host, key = parent(value, operation["path"])
        copied = copy.deepcopy(at(value, operation["from"]))
        if isinstance(host, list): host[int(key)] = copied
        else: host[key] = copied
    elif op == "wrap_expression":
        host, key = parent(value, operation["path"])
        old = host[int(key)] if isinstance(host, list) else host[key]
        wrapped = f"{operation['factor']}*({old})"
        if isinstance(host, list): host[int(key)] = wrapped
        else: host[key] = wrapped
    else:
        raise ValueError(f"unsupported fixture operation {op}")


def dump_yaml(path: Path, value: dict[str, Any]) -> None:
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=180))


def dump_artifact(path: Path, value: dict[str, Any]) -> None:
    if path.suffix == ".json": path.write_text(json.dumps(value, indent=2))
    else: dump_yaml(path, value)


def run(command: list[str], log_path: Path) -> tuple[int, str]:
    try:
        proc = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
        text = proc.stdout
        code = proc.returncode
    except subprocess.TimeoutExpired as exc:
        text = (exc.stdout or "") + "\nASSERT_FAIL:MUTATION_TIMEOUT:600s"
        code = 124
    log_path.write_text(text)
    return code, text


def run_engine(engine: str, command: list[str], log_path: Path) -> tuple[int, str]:
    if engine == "Mathematica":
        if not WOLFRAM_READY.is_set():
            with WOLFRAM_COLD_START:
                if not WOLFRAM_READY.is_set():
                    result = run(command, log_path)
                    WOLFRAM_READY.set()
                    return result
        with WOLFRAM_SEATS:
            return run(command, log_path)
    return run(command, log_path)


def asserted(text: str) -> str | None:
    match = ASSERTION.search(text)
    return match.group(1) if match else None


def engine_command(engine: str, input_path: Path, output_path: Path, phase_path: Path, source_path: Path | None = None) -> list[str]:
    if engine == "SymPy":
        return [sys.executable, str(source_path or HERE / "u1_body_mechanics_sympy.py"), "--input", str(input_path), "--output", str(output_path), "--phase-artifact", str(phase_path)]
    return ["wolframscript", "-file", str(source_path or HERE / "u1_body_mechanics_dual.wl"), "--input", str(input_path), "--output", str(output_path), "--phase-artifact", str(phase_path)]


def write_case_inputs(case_dir: Path, base_input: dict[str, Any], base_phase_input: dict[str, Any], fixture: dict[str, Any]) -> Path:
    data, phase_input = copy.deepcopy(base_input), copy.deepcopy(base_phase_input)
    if fixture.get("target") == "phase_input": mutate(phase_input, fixture["operation"])
    elif fixture.get("target") == "mechanics_input": mutate(data, fixture["operation"])
    phase_input_path = case_dir / "phase_a_input.yaml"
    mechanics_path = case_dir / "mechanics_input.yaml"
    dump_yaml(phase_input_path, phase_input)
    data["substrate"]["phase_a_input"] = str(phase_input_path.relative_to(ROOT))
    dump_yaml(mechanics_path, data)
    return mechanics_path


def expression(text: Any) -> sp.Expr:
    normalized = str(text).replace("b1", "").replace("ZZ", "_").replace("Pi", "pi").replace("^", "**")
    return sp.sympify(normalized, locals={"pi": sp.pi})


def canonical_expr(rows: list[dict[str, Any]]) -> sp.Expr:
    total: sp.Expr = sp.Integer(0)
    for row in rows:
        term: sp.Expr = sp.Float(row["coefficient"], 17)
        for name, power in row["powers"].items(): term *= sp.Symbol(name) ** int(power)
        total += term
    return sp.expand(total)


def has_parity(new: Any, old: Any, sign: int) -> bool:
    delta = sp.simplify(expression(new) - sign * expression(old))
    return bool(delta == 0)


def matrix_entry_parity(new: dict[str, Any], old: dict[str, Any], sign: int) -> bool:
    delta = canonical_expr(new["g4_control"]["Omega_physical_per_sheet_area"][0][1]) - sign * canonical_expr(old["g4_control"]["Omega_physical_per_sheet_area"][0][1])
    return abs(float(sp.N(delta.subs({symbol: 1 for symbol in delta.free_symbols}), 14))) < 1e-10


def force_parity(new: dict[str, Any], old: dict[str, Any], sign: int) -> bool:
    for name in old["g4_control"]["force_cases"]:
        before = old["g4_control"]["force_cases"][name]["F_extract"]
        after = new["g4_control"]["force_cases"][name]["F_extract"]
        if not all(has_parity(a, b, sign) for a, b in zip(after, before)): return False
    return True


def parity_check(new: dict[str, Any], old: dict[str, Any], parity: dict[str, int]) -> dict[str, bool]:
    ng, og = new["g4_control"], old["g4_control"]
    checks = {
        "k": has_parity(ng["computed_winding"], og["computed_winding"], parity["k"]),
        "Gamma": has_parity(ng["Gamma_computed"], og["Gamma_computed"], parity["Gamma"]),
        "epsilon": has_parity(ng["epsilon_matrix"][0][1], og["epsilon_matrix"][0][1], parity["epsilon"]),
        "Omega": matrix_entry_parity(new, old, parity["Omega"]),
        "F": force_parity(new, old, parity["F"]),
        "sigma": has_parity(ng["computed_sigma"], og["computed_sigma"], parity["sigma"]),
    }
    return checks


def route_residual_nonzero(artifact: dict[str, Any]) -> bool:
    return any(row for matrix_row in artifact["berry"]["equivalence_residual"] for row in matrix_row)


def run_engine_mutations(fixtures: dict[str, Any], base_input: dict[str, Any], base_phase_input: dict[str, Any],
                         baselines: dict[str, dict[str, Any]], phase_paths: dict[str, Path]) -> list[dict[str, Any]]:
    def one_case(fixture: dict[str, Any]) -> dict[str, Any]:
        case_dir = RUNTIME / fixture["id"]; case_dir.mkdir(parents=True)
        mechanics_path = write_case_inputs(case_dir, base_input, base_phase_input, fixture)
        def execute(engine: str) -> dict[str, Any]:
            output = case_dir / ("sympy.json" if engine == "SymPy" else "mathematica.json")
            phase_path = phase_paths[engine]
            if fixture.get("target") == "phase_artifact":
                phase_copy = copy.deepcopy(load(phase_path)); mutate(phase_copy, fixture["operation"])
                phase_path = case_dir / f"{engine.lower()}_phase_a.json"; dump_artifact(phase_path, phase_copy)
            mutated_source: Path | None = None
            if fixture.get("target") == "engine_source":
                original = HERE / ("u1_body_mechanics_sympy.py" if engine == "SymPy" else "u1_body_mechanics_dual.wl")
                replacement = fixture["source_replacements"][engine]
                source_text = original.read_text()
                if source_text.count(replacement["old"]) != 1:
                    raise RuntimeError(f"source mutation cardinality {fixture['id']}:{engine}")
                mutated_source = HERE / f".mutation_{fixture['id']}_{original.name}"
                mutated_source.write_text(source_text.replace(replacement["old"], replacement["new"]))
            try:
                code, log = run_engine(engine, engine_command(engine, mechanics_path, output, phase_path, mutated_source), case_dir / f"{engine.lower()}.log")
            finally:
                if mutated_source and mutated_source.exists(): mutated_source.unlink()
            tooth = asserted(log); artifact = load(output) if output.exists() else {}
            before = sha(at(baselines[engine], fixture["guarded_path"]))
            after = sha(at(artifact, fixture["guarded_path"])) if artifact else before
            observation = True; observation_record: dict[str, Any] = {"kind": "guard_digest_only"}
            if fixture.get("required_observation") == "sigma_changes":
                baseline_sigma = baselines[engine]["g4_control"]["computed_sigma"]
                mutated_sigma = artifact.get("g4_control", {}).get("computed_sigma")
                observation = mutated_sigma != baseline_sigma
                observation_record = {"kind": "sigma_changes", "baseline_sigma": baseline_sigma, "mutated_sigma": mutated_sigma}
            elif fixture.get("required_observation") == "route_equivalence_nonzero":
                observation = bool(artifact) and route_residual_nonzero(artifact)
                observation_record = {"kind": "route_equivalence_nonzero", "residual_digest": sha(artifact["berry"]["equivalence_residual"]) if artifact else None}
            passed = code == 1 and tooth == fixture["expected_assert"] and tooth != "MUTATION_NOOP" and before != after and observation
            return {"engine": engine, "exit_code": code, "asserted": tooth, "guarded_digest_changed": before != after,
                    "required_observation_passed": observation, "observation": observation_record, "status": "PASS" if passed else "FAIL"}
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            engine_rows = list(pool.map(execute, ("SymPy", "Mathematica")))
        passed = all(row["status"] == "PASS" for row in engine_rows)
        return {"id": fixture["id"], "kind": "engine_input_or_source", "expected_assert": fixture["expected_assert"],
                "guarded_path": fixture["guarded_path"], "guarded_digest_changed": all(r["guarded_digest_changed"] for r in engine_rows),
                "assertion_fired": all(r["asserted"] == fixture["expected_assert"] for r in engine_rows),
                "engine_results": engine_rows, "status": "PASS" if passed else "FAIL"}
    # The semaphore, rather than the number of case workers, enforces the two
    # Wolfram-seat ceiling.  Additional workers let independent SymPy jobs run
    # while Wolfram cases wait for a seat, keeping the whole script below its
    # 600-second cap.
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(8, len(fixtures["engine_mutations"]))) as pool:
        rows = list(pool.map(one_case, fixtures["engine_mutations"]))
    for row in rows:
        print(f"B1_EXTERNAL_MUTATION:{row['id']} {row['status']}")
        if row["status"] != "PASS": raise RuntimeError(f"external mutation failed: {row['id']}: {row['engine_results']}")
    return rows


def run_metamorphic(fixtures: dict[str, Any], base_input: dict[str, Any], base_phase_input: dict[str, Any],
                    baselines: dict[str, dict[str, Any]], phase_paths: dict[str, Path]) -> list[dict[str, Any]]:
    def one_control(fixture: dict[str, Any]) -> dict[str, Any]:
        case_dir = RUNTIME / fixture["id"]; case_dir.mkdir(parents=True)
        data = copy.deepcopy(base_input)
        for operation in fixture["operations"]: mutate(data, operation)
        mechanics_path = case_dir / "mechanics_input.yaml"; dump_yaml(mechanics_path, data)
        def execute(engine: str) -> dict[str, Any]:
            output = case_dir / ("sympy.json" if engine == "SymPy" else "mathematica.json")
            code, log = run_engine(engine, engine_command(engine, mechanics_path, output, phase_paths[engine]), case_dir / f"{engine.lower()}.log")
            artifact = load(output) if output.exists() else {}
            parity = parity_check(artifact, baselines[engine], fixture["parity"]) if artifact else {}
            passed = code == 0 and asserted(log) is None and parity and all(parity.values())
            return {"engine": engine, "exit_code": code, "parity_checks": parity, "status": "PASS" if passed else "FAIL"}
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            engine_rows = list(pool.map(execute, ("SymPy", "Mathematica")))
        passed = all(row["status"] == "PASS" for row in engine_rows)
        return {"id": fixture["id"], "expected_parity": fixture["parity"], "engine_results": engine_rows, "status": "PASS" if passed else "FAIL"}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        controls = list(pool.map(one_control, fixtures["metamorphic_controls"]))
    for row in controls:
        print(f"B1_METAMORPHIC:{row['id']} {row['status']}")
        if row["status"] != "PASS": raise RuntimeError(f"metamorphic control failed: {row['id']}: {row['engine_results']}")
    return controls


def run_comparator_mutations(fixtures: dict[str, Any], base_input_path: Path, baselines: dict[str, dict[str, Any]],
                             phase_paths: dict[str, Path]) -> list[dict[str, Any]]:
    def one(fixture: dict[str, Any]) -> dict[str, Any]:
        case_dir = RUNTIME / fixture["id"]; case_dir.mkdir(parents=True)
        targets = {"sympy_artifact": copy.deepcopy(baselines["SymPy"]), "math_artifact": copy.deepcopy(baselines["Mathematica"]),
                   "phase_artifact": load(phase_paths["SymPy"])}
        target = targets[fixture["target"]]; before = sha(at(target, fixture["guarded_path"])); mutate(target, fixture["operation"]); after = sha(at(target, fixture["guarded_path"]))
        sym_path, math_path, phase_path = case_dir / "sympy.yaml", case_dir / "math.yaml", case_dir / "phase.json"
        dump_artifact(sym_path, targets["sympy_artifact"]); dump_artifact(math_path, targets["math_artifact"]); dump_artifact(phase_path, targets["phase_artifact"])
        command = [sys.executable, str(HERE / "u1_body_mechanics_compare.py"), "--input", str(base_input_path), "--sympy-artifact", str(sym_path),
                   "--math-artifact", str(math_path), "--phase-a-artifact", str(phase_path), "--verify-only"]
        code, log = run(command, case_dir / "comparator.log"); tooth = asserted(log)
        passed = code == 1 and tooth == fixture["expected_assert"] and tooth != "MUTATION_NOOP" and before != after
        return {"id": fixture["id"], "kind": "comparator_artifact", "expected_assert": fixture["expected_assert"],
                "guarded_path": fixture["guarded_path"], "guarded_digest_changed": before != after, "assertion_fired": tooth == fixture["expected_assert"],
                "exit_code": code, "asserted": tooth, "status": "PASS" if passed else "FAIL"}

    with concurrent.futures.ThreadPoolExecutor(max_workers=min(6, len(fixtures["comparator_mutations"]))) as pool:
        rows = list(pool.map(one, fixtures["comparator_mutations"]))
    for row in rows:
        print(f"B1_EXTERNAL_MUTATION:{row['id']} {row['status']}")
        if row["status"] != "PASS": raise RuntimeError(f"comparator mutation failed: {row['id']}: {row}")
    return rows


def mutate_leaf(data: dict[str, Any], path: str) -> tuple[str, str, bool]:
    parts: list[Any] = [int(part) if part.isdigit() else part for part in path.split(".")]
    host, key = parent(data, parts)
    old = host[int(key)] if isinstance(host, list) else host[key]
    if isinstance(old, bool): new: Any = not old
    elif isinstance(old, int): new = old + 1
    elif isinstance(old, float): new = old + 0.5
    elif old is None: new = "__LIVENESS_NON_NULL__"
    elif isinstance(old, str): new = old + "__LIVENESS_MUTATION__"
    elif isinstance(old, list) and not old: new = ["__LIVENESS_MEMBER__"]
    else: raise TypeError(f"cannot mutate leaf {path}: {type(old).__name__}")
    if isinstance(host, list): host[int(key)] = new
    else: host[key] = new
    discrete_tokens = ("schema_version", "directive_version", "spec_version", "spatial_dimension", "radial_dimension",
                       "harmonic_degree", "dimensions.", "coordinate_dimensions", "derivative_order_bound")
    semantically_valid = isinstance(old, (int, float)) and not isinstance(old, bool) and not any(token in path for token in discrete_tokens)
    return sha(old), sha(new), semantically_valid


def run_input_liveness(base_input: dict[str, Any], base_phase_input: dict[str, Any],
                       baselines: dict[str, dict[str, Any]], phase_paths: dict[str, Path],
                       artifacts: Path, batch_index: int = 0, batch_count: int = 1) -> list[dict[str, Any]]:
    baseline = baselines["SymPy"]
    all_rows = baseline["input_liveness"]["rows"]
    if batch_count < 1 or not 0 <= batch_index < batch_count: raise ValueError("invalid liveness batch")
    rows = [(index, row) for index, row in enumerate(all_rows) if index % batch_count == batch_index]
    baseline_by_path = {row["path"]: row for row in all_rows}
    math_by_path = {row["path"]: row for row in baselines["Mathematica"]["input_liveness"]["rows"]}
    live_root = RUNTIME / "input_liveness"; live_root.mkdir(parents=True, exist_ok=True)
    dual_sample_indices: set[int] = set()
    for index, row in rows:
        probe = copy.deepcopy(base_input)
        _, _, semantically_valid = mutate_leaf(probe, row["path"])
        if semantically_valid and not row["path"].startswith("phase_a_protection."):
            dual_sample_indices.add(index)
            if len(dual_sample_indices) == 2:
                break

    def one(index_and_row: tuple[int, dict[str, Any]]) -> dict[str, Any]:
        index, row = index_and_row; path = row["path"]
        case_dir = live_root / f"{index:04d}"; case_dir.mkdir(parents=True, exist_ok=True)
        data = copy.deepcopy(base_input); before_value, after_value, semantically_valid = mutate_leaf(data, path)
        mechanics_path = case_dir / "mechanics_input.yaml"; dump_yaml(mechanics_path, data)
        output = case_dir / "sympy.json"
        if path.startswith("phase_a_protection."):
            command = [sys.executable, str(HERE / "u1_body_mechanics_compare.py"), "--input", str(mechanics_path),
                       "--sympy-artifact", str(artifacts / "sympy_phase_b1.yaml"),
                       "--math-artifact", str(artifacts / "mathematica_phase_b1.yaml"),
                       "--phase-a-artifact", str(phase_paths["SymPy"]), "--verify-only"]
            code, log = run(command, case_dir / "comparator.log"); artifact = None; executor = "Comparator"
        else:
            code, log = run(engine_command("SymPy", mechanics_path, output, phase_paths["SymPy"]), case_dir / "sympy.log")
            artifact = load(output) if output.exists() else None; executor = "SymPy"
        tooth = asserted(log)
        baseline_sink = sha(dotted_at(baseline, row["semantic_sink"]))
        changed_sink = False
        if artifact is not None:
            try: changed_sink = sha(dotted_at(artifact, row["semantic_sink"])) != baseline_sink
            except (KeyError, IndexError, TypeError, ValueError): changed_sink = True
        guarded_assertion = code == 1 and tooth is not None and tooth.startswith("B1_") and tooth != "MUTATION_NOOP"
        engine_results = [{"engine": executor, "exit_code": code, "asserted": tooth,
                           "sink_digest_changed": changed_sink, "guarded_assertion_fired": guarded_assertion}]
        dual_passed = True
        if index in dual_sample_indices:
            math_output = case_dir / "mathematica.json"
            math_code, math_log = run_engine("Mathematica", engine_command("Mathematica", mechanics_path, math_output,
                                                  phase_paths["Mathematica"]), case_dir / "mathematica.log")
            math_artifact = load(math_output) if math_output.exists() else None
            math_tooth = asserted(math_log)
            math_sink = math_by_path[path]["semantic_sink"]
            math_baseline_sink = sha(dotted_at(baselines["Mathematica"], math_sink))
            math_changed = False
            if math_artifact is not None:
                try: math_changed = sha(dotted_at(math_artifact, math_sink)) != math_baseline_sink
                except (KeyError, IndexError, TypeError, ValueError): math_changed = True
            math_assertion = math_code == 1 and math_tooth is not None and math_tooth.startswith("B1_") and math_tooth != "MUTATION_NOOP"
            dual_passed = math_code != 124 and (math_changed if semantically_valid else (math_changed or math_assertion))
            engine_results.append({"engine": "Mathematica", "exit_code": math_code, "asserted": math_tooth,
                                   "sink_digest_changed": math_changed, "guarded_assertion_fired": math_assertion})
        passed = (code != 124 and before_value != after_value and
                  (changed_sink if semantically_valid else (changed_sink or guarded_assertion)) and dual_passed)
        return {"path": path, "liveness_batch_index": batch_index, "liveness_batch_count": batch_count,
                "engine": "SymPy+Mathematica" if index in dual_sample_indices else executor,
                "dual_engine_sampled": index in dual_sample_indices, "engine_results": engine_results,
                "typed_role": baseline_by_path[path]["typed_role"],
                "semantic_sink": row["semantic_sink"], "input_digest_changed": before_value != after_value,
                "sink_digest_changed": changed_sink, "guarded_assertion_fired": guarded_assertion,
                "mutation_class": "semantically_valid" if semantically_valid else "declared_malformed",
                "asserted": tooth, "exit_code": code, "out_of_process": True, "status": "PASS" if passed else "FAIL"}

    workers = min(16, max(2, os.cpu_count() or 2))
    completed = 0; output_rows: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(one, item) for item in rows]
        for future in concurrent.futures.as_completed(futures):
            result = future.result(); output_rows.append(result); completed += 1
            if completed % 25 == 0 or completed == len(rows): print(f"B1_INPUT_LIVENESS: {completed}/{len(rows)}", flush=True)
    output_rows.sort(key=lambda row: row["path"])
    failed = [row for row in output_rows if row["status"] != "PASS"]
    if failed: raise RuntimeError(f"input liveness failures ({len(failed)}): {failed[:4]}")
    return output_rows


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES); parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--artifacts", type=Path, default=DEFAULT_ARTIFACTS); parser.add_argument("--output", type=Path)
    parser.add_argument("--core-only", action="store_true"); parser.add_argument("--core-results", type=Path)
    parser.add_argument("--core-batch"); parser.add_argument("--core-batch-results", type=Path, nargs="+")
    parser.add_argument("--liveness-batch"); parser.add_argument("--liveness-results", type=Path, nargs="+")
    args = parser.parse_args(); output = args.output or args.artifacts / "b1_mutation_results.yaml"
    try:
        if args.core_only and args.core_results is not None: raise ValueError("--core-only and --core-results are mutually exclusive")
        if args.core_batch and args.core_batch_results: raise ValueError("--core-batch and --core-batch-results are mutually exclusive")
        if (args.core_batch or args.core_batch_results) and not args.core_only: raise ValueError("core batching requires --core-only")
        if (args.liveness_batch or args.liveness_results) and args.core_results is None: raise ValueError("liveness batching requires --core-results")
        if args.liveness_batch and args.liveness_results: raise ValueError("--liveness-batch and --liveness-results are mutually exclusive")
        fixtures = yaml.safe_load(args.fixtures.read_text()); base_input = yaml.safe_load(args.input.read_text())
        base_phase_input = yaml.safe_load((ROOT / base_input["substrate"]["phase_a_input"]).read_text())
        phase_paths = {"SymPy": args.artifacts / "sympy_phase_a.json", "Mathematica": args.artifacts / "mathematica_phase_a.json"}
        baseline_paths = {"SymPy": args.artifacts / "sympy_phase_b1.yaml", "Mathematica": args.artifacts / "mathematica_phase_b1.yaml"}
        baselines = {name: load(path) for name, path in baseline_paths.items()}
        baseline_digests = {name: sha(value) for name, value in baselines.items()}
        if fixtures.get("schema_version") != "U1_PHASE_B1_EXTERNAL_FIXTURES_V3": raise ValueError("fixtures are not V3")
        if args.core_results is None and RUNTIME.exists(): shutil.rmtree(RUNTIME)
        RUNTIME.mkdir(parents=True, exist_ok=True)
        required = {f"B1_R{i}" for i in range(1, 10)} | {"B1_C_ENGINE_MATH", "B1_C_RECORD_MAP", "B1_C_PHASE_A", "B1_C_REPORT_SCHEMA", "B1_C_INPUT_LIVENESS"}

        if args.core_batch_results:
            batches = [load(path) for path in args.core_batch_results]
            counts = {batch.get("batch_count") for batch in batches}
            indices = {batch.get("batch_index") for batch in batches}
            if (len(counts) != 1 or counts != {len(batches)} or indices != set(range(len(batches))) or
                    any(batch.get("schema_version") != "U1_PHASE_B1_MUTATION_CORE_BATCH_V3" or batch.get("status") != "PASS" or
                        batch.get("baseline_digests") != baseline_digests for batch in batches)):
                raise ValueError("mutation core batch artifacts are incomplete, invalid, or stale")
            ordered = sorted(batches, key=lambda item: item["batch_index"])
            cases = [row for batch in ordered for row in batch["cases"]]
            controls = [row for batch in ordered for row in batch["metamorphic_controls"]]
            covered = {row["expected_assert"] for row in cases}
            core = {"schema_version": "U1_PHASE_B1_MUTATION_CORE_V3", "status": "PASS" if required <= covered else "FAIL",
                    "production_executables_mutation_unaware": True, "sentinel_rejected": all(row["status"] == "PASS" for row in cases),
                    "covered_assertions": sorted(covered), "baseline_digests": baseline_digests,
                    "cases": cases, "metamorphic_controls": controls,
                    "source_batches": len(batches)}
            output.parent.mkdir(parents=True, exist_ok=True); dump_yaml(output, core)
            if core["status"] != "PASS": raise RuntimeError(f"missing mutation coverage {required-covered}")
            print(f"B1_MUTATION_CORE: PASS batches={len(batches)} cases={len(cases)} controls={len(controls)}"); return 0

        if args.core_batch:
            parts = args.core_batch.split("/", 1)
            if len(parts) != 2: raise ValueError("core batch must be INDEX/COUNT")
            batch_index, batch_count = map(int, parts)
            if batch_count < 1 or not 0 <= batch_index < batch_count: raise ValueError("invalid core batch")
            local = copy.deepcopy(fixtures)
            local["engine_mutations"] = [row for index, row in enumerate(fixtures["engine_mutations"])
                                         if index % batch_count == batch_index]
            local["comparator_mutations"] = [row for index, row in enumerate(fixtures["comparator_mutations"])
                                              if index % batch_count == batch_index]
            local["metamorphic_controls"] = [row for index, row in enumerate(fixtures["metamorphic_controls"])
                                              if index % batch_count == batch_index]
            cases = run_engine_mutations(local, base_input, base_phase_input, baselines, phase_paths)
            controls = run_metamorphic(local, base_input, base_phase_input, baselines, phase_paths)
            cases += run_comparator_mutations(local, args.input, baselines, phase_paths)
            batch = {"schema_version": "U1_PHASE_B1_MUTATION_CORE_BATCH_V3",
                     "status": "PASS" if all(row["status"] == "PASS" for row in cases + controls) else "FAIL",
                     "batch_index": batch_index, "batch_count": batch_count, "baseline_digests": baseline_digests,
                     "cases": cases, "metamorphic_controls": controls}
            output.parent.mkdir(parents=True, exist_ok=True); dump_yaml(output, batch)
            if batch["status"] != "PASS": raise RuntimeError(f"mutation core batch failed {batch_index}/{batch_count}")
            print(f"B1_MUTATION_CORE_BATCH: PASS batch={batch_index}/{batch_count} cases={len(cases)} controls={len(controls)}"); return 0

        if args.core_results is not None:
            core = load(args.core_results)
            if (core.get("schema_version") != "U1_PHASE_B1_MUTATION_CORE_V3" or core.get("status") != "PASS"
                    or core.get("baseline_digests") != baseline_digests):
                raise ValueError("mutation core artifact is invalid or stale")
            cases = core["cases"]; controls = core["metamorphic_controls"]
        else:
            cases = run_engine_mutations(fixtures, base_input, base_phase_input, baselines, phase_paths)
            controls = run_metamorphic(fixtures, base_input, base_phase_input, baselines, phase_paths)
            cases += run_comparator_mutations(fixtures, args.input, baselines, phase_paths)
            covered = {row["expected_assert"] for row in cases}
            core = {"schema_version": "U1_PHASE_B1_MUTATION_CORE_V3", "status": "PASS" if required <= covered else "FAIL",
                    "production_executables_mutation_unaware": True, "sentinel_rejected": all(row["status"] == "PASS" for row in cases),
                    "covered_assertions": sorted(covered), "baseline_digests": baseline_digests,
                    "cases": cases, "metamorphic_controls": controls}
            if args.core_only:
                output.parent.mkdir(parents=True, exist_ok=True); dump_yaml(output, core)
                if core["status"] != "PASS": raise RuntimeError(f"missing mutation coverage {required-covered}")
                print(f"B1_MUTATION_CORE: PASS cases={len(cases)} controls={len(controls)}"); return 0

        core_digest = sha(core)
        if args.liveness_results:
            batches = [load(path) for path in args.liveness_results]
            counts = {batch.get("batch_count") for batch in batches}
            indices = {batch.get("batch_index") for batch in batches}
            if (len(counts) != 1 or counts != {len(batches)} or indices != set(range(len(batches)))
                    or any(batch.get("schema_version") != "U1_PHASE_B1_INPUT_LIVENESS_BATCH_V3" or batch.get("status") != "PASS"
                           or batch.get("baseline_digests") != baseline_digests or batch.get("core_results_sha256") != core_digest
                           for batch in batches)):
                raise ValueError("liveness batch artifacts are incomplete, invalid, or stale")
            liveness_cases = [row for batch in sorted(batches, key=lambda item: item["batch_index"])
                              for row in batch["input_liveness_cases"]]
        else:
            batch_index, batch_count = 0, 1
            if args.liveness_batch:
                parts = args.liveness_batch.split("/", 1)
                if len(parts) != 2: raise ValueError("liveness batch must be INDEX/COUNT")
                batch_index, batch_count = map(int, parts)
            liveness_cases = run_input_liveness(base_input, base_phase_input, baselines, phase_paths, args.artifacts,
                                                batch_index, batch_count)
            if args.liveness_batch:
                batch = {"schema_version": "U1_PHASE_B1_INPUT_LIVENESS_BATCH_V3", "status": "PASS",
                         "batch_index": batch_index, "batch_count": batch_count, "baseline_digests": baseline_digests,
                         "core_results_sha256": core_digest, "input_liveness_cases": liveness_cases}
                output.parent.mkdir(parents=True, exist_ok=True); dump_yaml(output, batch)
                print(f"B1_INPUT_LIVENESS_BATCH: PASS batch={batch_index}/{batch_count} cases={len(liveness_cases)}"); return 0

        covered = {row["expected_assert"] for row in cases}
        declared_paths = set(baselines["SymPy"]["input_liveness"]["declared_leaf_paths"])
        passed_paths = [row["path"] for row in liveness_cases if row["status"] == "PASS"]
        live_ok = set(passed_paths) == declared_paths and len(passed_paths) == len(declared_paths) == len(liveness_cases)
        result = {"schema_version": "U1_PHASE_B1_OUT_OF_PROCESS_MUTATIONS_V3", "status": "PASS" if required <= covered and live_ok else "FAIL",
                  "production_executables_mutation_unaware": True, "sentinel_rejected": all(row["status"] == "PASS" for row in cases),
                  "covered_assertions": sorted(covered), "cases": cases, "metamorphic_controls": controls,
                  "input_liveness_cases": liveness_cases, "input_liveness_declared_count": len(declared_paths),
                  "dual_engine_liveness_sample_count": sum(bool(row.get("dual_engine_sampled")) for row in liveness_cases),
                  "core_results_sha256": core_digest}
        output.parent.mkdir(parents=True, exist_ok=True); dump_yaml(output, result)
        if result["status"] != "PASS": raise RuntimeError(f"missing mutation coverage {required-covered}")
        print(f"B1_OUT_OF_PROCESS_HARNESS: PASS cases={len(cases)} controls={len(controls)}"); return 0
    except Exception as exc:
        print(f"ASSERT_FAIL:MUTATION_HARNESS:{type(exc).__name__}:{exc}", file=sys.stderr); return 1


if __name__ == "__main__":
    raise SystemExit(main())
