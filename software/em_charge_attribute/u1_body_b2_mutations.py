#!/usr/bin/env python3
"""Out-of-process production comparator ablations; engines stay mutation-unaware."""

from __future__ import annotations

import argparse
import copy
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable

from u1_body_b2_common import digest, dump_yaml, load_yaml, require, sha256_file


Mutation = Callable[[dict[str, Any], dict[str, Any]], None]


def cases() -> list[tuple[str, str, Mutation]]:
    def one(path: list[Any], value: Any) -> Mutation:
        def mutate(sym: dict[str, Any], _: dict[str, Any]) -> None:
            obj: Any = sym
            for key in path[:-1]: obj = obj[key]
            obj[path[-1]] = value
        return mutate
    return [
        ("engine_identity", "B2_C_ENGINE", one(["engine"], "Mathematica")),
        ("route_independence", "B2_C_INDEPENDENCE", lambda s, m: m.__setitem__("independent_route", s["independent_route"])),
        ("grid_omission", "B2_C_GRID", lambda s, m: s["grid"]["cells"].pop()),
        ("C_mdot_ablation", "B2_C_CELL", one(["grid", "cells", 0, "C_mdot", "status"], "ZERO(padded)")),
        ("velocity_alias_ablation", "B2_C_CELL", one(["grid", "cells", 0, "velocity_block", "alias"], "independent_D")),
        ("pdot_omission", "B2_C_CELL", one(["grid", "cells", 0, "pdot_generalized_velocity_remainder", "status"], "OMITTED")),
        ("double_count", "B2_C_CELL", one(["grid", "cells", 0, "velocity_block", "count_in_P_local"], 2)),
        ("operator_coefficient", "B2_C_OPERATORS", one(["operator_inventory", 0, "O"], "1")),
        ("partition_owner", "B2_C_PARTITION", one(["partition", "terminal_owner_enum"], ["M", "pending"])),
        ("radiation_total", "B2_C_RADIATION", one(["radiation", "totals", "F_rad"], "ZERO(padded)")),
        ("phase_C_execution", "B2_C_PHASE_C", one(["phase_C", "G8"], "PASS")),
        ("G9_router", "B2_C_G9", lambda s, m: [row["G9"]["energy"].__setitem__("causes", ["return_energy_closure"]) for artifact in (s, m) for row in artifact["grid"]["cells"]]),
        ("resolvent_identity", "B2_C_RESOLVENT", lambda s, m: [artifact["radiation"]["per_cell_resolvent_identity"][0]["retarded_green_operator"].__setitem__("contact_residual", "1") for artifact in (s, m)]),
        ("known_nonzero_control", "B2_C_CONTROL", lambda s, m: [artifact["radiation"]["per_cell_known_nonzero_control"][0].__setitem__("status", "ZERO(padded)") for artifact in (s, m)]),
        ("obligation_missing", "B2_C_OBLIGATIONS_EXPECTED_MINUS_REACHABLE", lambda s, m: [(artifact["typed_dag"]["nodes"].pop("mode_coverage_residual"), artifact["typed_dag"]["nodes"]["report_headline"]["depends_on"].remove("mode_coverage_residual")) for artifact in (s, m)]),
    ]


def main() -> int:
    p = argparse.ArgumentParser(); p.add_argument("--sympy", type=Path, required=True); p.add_argument("--mathematica", type=Path, required=True); p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True); p.add_argument("--work", type=Path, required=True); p.add_argument("--output", type=Path, required=True); a = p.parse_args()
    try:
        base_s, base_m = load_yaml(a.sympy), load_yaml(a.mathematica); a.work.mkdir(parents=True, exist_ok=True); rows = []
        comparator = Path(__file__).parent / "u1_body_b2_compare.py"
        for name, tooth, mutate in cases():
            s, m = copy.deepcopy(base_s), copy.deepcopy(base_m); before = digest({"s": s, "m": m}); mutate(s, m); require(before != digest({"s": s, "m": m}), "MUTATION_NOOP", name)
            work = a.work / name; work.mkdir(parents=True, exist_ok=True); spath, mpath = work / "sympy.yaml", work / "mathematica.yaml"; dump_yaml(spath, s); dump_yaml(mpath, m)
            seal = work / "producer.yaml"
            dump_yaml(seal, {"schema_version": "U1_PHASE_B2_PRODUCER_SEAL_V1", "status": "PASS", "producer": f"production_mutation:{name}", "outputs": [{"path": str(spath.resolve()), "sha256": sha256_file(spath), "producer": f"production_mutation:{name}"}, {"path": str(mpath.resolve()), "sha256": sha256_file(mpath), "producer": f"production_mutation:{name}"}]})
            proc = subprocess.run([sys.executable, str(comparator), "--sympy", str(spath), "--mathematica", str(mpath), "--stage0", str(a.stage0.resolve()), "--stage0-contract-digest", a.stage0_contract_digest, "--producer-record", str(seal), "--output", str(work / "agreement.yaml")], cwd=comparator.parent, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
            marker = f"ASSERT_FAIL:{tooth}:"; require(proc.returncode == 1 and marker in proc.stdout, "MUTATION_NOOP", f"{name}:{proc.returncode}:{proc.stdout[-300:]}")
            rows.append({"id": name, "expected_assert": tooth, "killed_at_own_assert": True, "sink_digest_changed": True})
        dump_yaml(a.output, {"schema_version": "U1_PHASE_B2_MUTATIONS_V1", "status": "PASS", "out_of_process": True, "engines_mutation_unaware": True, "MUTATION_NOOP_sentinel": "armed", "case_count": len(rows), "cases": rows})
        print(f"B2_MUTATIONS: PASS cases={len(rows)}"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__": raise SystemExit(main())
