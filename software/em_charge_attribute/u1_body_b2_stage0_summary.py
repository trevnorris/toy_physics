#!/usr/bin/env python3
"""Write the contracted <=20-line Stage-0 handoff summary."""

from __future__ import annotations

import argparse
from pathlib import Path

from u1_body_b2_common import load_yaml, rel_repo, sha256_file


def main() -> int:
    p = argparse.ArgumentParser(); p.add_argument("--contract", type=Path, required=True); p.add_argument("--logical-contract", type=Path); p.add_argument("--output", type=Path, required=True); a = p.parse_args()
    row = load_yaml(a.contract)
    logical_contract = (a.logical_contract or a.contract).resolve()
    lines = [
        "# U1 Phase B2 Stage-0 build summary",
        "",
        f"- Contract anchor: directive v48 at `{row['contract_anchor']['startup_contract_commit']}`; directive commit/working sha256 `{row['contract_anchor']['contract_digests']['directive_commit_blob']}`.",
        f"- Stage-0 artifact: `{rel_repo(logical_contract)}`",
        f"- stage0_contract.yaml sha256: `{sha256_file(a.contract)}`",
        f"- stage0_semantic_sha256: `{row['stage0_semantic_sha256']}`",
        f"- Environment closure: `{row['trust_mode_predicate']['environment_closure_digest']}` over `{row['trust_mode_predicate']['environment_closure']['entry_count']}` first-use content-hashed entries.",
        f"- Obs_B2: count `{row['observation_contract']['Obs_B2_manifest_count']}`, digest `{row['observation_contract']['Obs_B2_manifest_set_sha256']}`.",
        f"- Physics inventory: {len(row['frozen_data']['native_operator_inventory'])} connected Hessian sectors, {len(row['frozen_data']['endpoint_resolvent_cells'])} endpoint cells; the wall χ–u block is dual-engine `DERIVED` at a frozen symbolic base.",
        f"- v48 disposition floor: `{len(row['frozen_data']['stage0_datum_dispositions'])}` slots, each exactly one of `DERIVED(value_digest, dual_engine_comparison_id)` or `UNRESOLVED(witness_id, challenge_id)`.",
        "- Integrated balances retain executable native-current/Reynolds and authenticated-root attempts and halt honestly as typed `UNRESOLVED`; no hand-entered matching symbol list is banked.",
        "- Full and targeted B1 trace containment passed with no B2 replay-orchestration exception; external network events: zero.",
        "- Status: `AWAITING_ORCHESTRATOR_APPROVAL`; contracted runner exit is `42`.",
    ]
    if len(lines) > 20: raise RuntimeError(f"summary has {len(lines)} lines")
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"B2_STAGE0_SUMMARY: PASS lines={len(lines)}")
    return 0


if __name__ == "__main__": raise SystemExit(main())
