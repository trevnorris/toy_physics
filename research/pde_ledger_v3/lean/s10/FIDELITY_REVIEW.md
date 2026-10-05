# S10 Lean contract completion and fidelity review

Completed 2026-09-11 under [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).
**The scoped S10 Lean task is complete.** The deliverable is the six-family
mathematical classification, compact action/operator identification, exhaustive
coverage contract and meaningful mutation controls. This does not clear the
separate S10 CAS production, comparator/export or ledger/paper work.

## Fixed revision and independent reviews

Codex authored the contract. Two fresh non-author review sessions received the
same 157-file packet, with no access to one another's findings. The
original prompt (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_review_prompt.txt`) asked
for statement fidelity, explicit domains, action/parameter normalization,
exceptional cases, full subspaces and meaningful controls. It excluded
systematic per-output bridge expansion.

- Revision: `S10 compact fidelity contract v1`.
- Packet SHA-256: `96d2c9d90a32b26d1c031e43ccc8580a8f8a102e18425491b231d4ac0cbc3b35`.
- Policy commit: `c2bdb663`; canonical proof checkpoint: `f21e6459`.
- Manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_review_packet.json`): each
  file hash; aggregate is SHA-256 of Python `json.dumps(files, sort_keys=True)`
  encoded as UTF-8. The author recomputed the aggregate and all 157 file hashes.
- Exact reviewed contract (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_contract_v1.md`):
  SHA-256 `013cf4df7b9cfcb55881550dca3bc9d9effef417343d488c0cda917123c1c5db`.
  The live [COVERAGE.md](COVERAGE.md) includes the editorial and status changes
  described below. Canonical proofs and the control instrument are unchanged.

At closure, 153 packet files still match the live tree. The reviewed contract
is preserved in the snapshot above; the three edited navigation documents
(`lean/README.md`, `lean/s10/README.md`, `lean/s10/CAS_CHECKPOINT.md`) have their
reviewed contents in policy commit `c2bdb663`. All 157 reviewed file hashes were
also checked using those durable originals after the status updates.

| Reviewer | Session | Verdict and evidence |
|---|---|---|
| Claude Opus 5 (`claude-opus-5[1m]` in CLI metadata) | `cad2e599-63c7-4417-878d-5bd5c3f7aed7` | **CLEAR**; [report](../../_measurements/S10_lean_fidelity_claude_v1.md), raw response (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_claude_v1.json`), run (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_claude_v1_run.json`) |
| Grok 4.6 (`grok-4.6-build`) | `d419314e-5241-48c0-9a9e-9fc04c8a891b` | **CLEAR**; [report](../../_measurements/S10_lean_fidelity_grok_v1.md), raw final response (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_grok_v1_completion.json`), completion run (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_grok_v1_completion_run.json`) |

Grok's first run (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_grok_v1_run.json`)
returned `stopReason: cancelled` after requesting a shell hash check. Its
partial response (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_grok_v1.json`) is not
review clearance. The same independent session resumed with read/search tools
only and this prompt (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_grok_resume_prompt.txt`),
without Claude's findings. The completed response returned `end_turn` and the
full report above. Both completed reports were collected before adjudication.

## Findings and disposition

Neither reviewer identified a substantive blocker. All comments were editorial
or optional evidence improvements; none requires changing a mathematical claim
or extending the bridge.

| Comment | Disposition |
|---|---|
| Claude E1: distinguish the first amplitude component from a coordinate label | Clarified displacement/velocity component, Lean index `0`. |
| Claude E2: state ANISO `sigma` is dimensionless | Added the declaration already present in both CAS constructors and the dimensional proofs. |
| Claude E3 / Grok: `c != 1` concerns distinctness from MAIN | Retained the nontrivial-control domain and stated that root/count formulas also hold at `c=1`; both reviews confirm that existing coverage. |
| Claude E4 / Grok: make implicit static N3 anchors explicit | Added scalar `zero_counts` at `c=1`, `longitudinal_inf_transverse`, and `div_zero_space` with `T ∩ T = T`. |
| Grok: scalar stationarity inherits the MAIN theorems | Named `actionStationary_eq` and the inheritance explicitly in the evidence table. |
| Claude E5: some mutants also emit linter errors; rejection regex scans whole output | No instrument change needed for this run: both reviewers inspected the mathematical failures inside the required declarations. ANISO fails `kinetic_eq`; ignored scaling fails `lagrangian_eq`. The linter errors are additional diagnostics, not the evidence of rejection. Linters remain enabled. A future instrument revision can narrow diagnostic matching if needed. |
| Claude E6: older CAS mutation JSON stores hashes rather than full mutant text | Retained the committed record. The existing instrument (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_check.py`) contains the mutation constructions; its diagnostics and positive controls were reviewed. No unrelated evidence-format migration is required for closure. |

These clarifications implement the reviewers' source-backed observations. They
do not alter definitions, hypotheses of formal theorems, the classification, or
control execution. Status changes and navigation links record completion rather
than constituting a new physics claim.

## Completion evidence and stopping point

| Item | Closed by |
|---|---|
| C1: theorem and exhaustive coverage contract | Existing general six-family proofs and the contract's theorem map, including parallel/perpendicular/oblique ANISO cases, root coincidence, and full/transverse counts. Both reviews clear the mapping and domain. |
| C2: compact fidelity link | Actual action selectors and parameter/index maps in both CAS constructors; Lean action/operator identities; retained, explicitly limited ANISO D3 imported-matrix identity. Both reviews clear the stated connection. |
| C3: meaningful controls | [Closure controls](../../_measurements/S10_lean_contract_checks.json): seven mathematical rejections and six passing controls, plus the mapped committed spectrum/stratum, matrix and Q6/Q7 records. Both reviews found the inspected failures substantive and hypotheses nonvacuous. |
| C4: independent review and verification | Both CLEAR reports, the dispositions above, and [CONTRACT_BUILD_VERIFICATION.txt](CONTRACT_BUILD_VERIFICATION.txt): all five libraries build; 1962 selected axiom audits; no admissions/custom axioms; all 135 original artifact hashes match. All 115 canonical source hashes still match the tested revision. |

The reviewers inspected source and recorded diagnostics; they did not rerun
Lean, the CAS engines or mutation instruments, and did not independently
recompute the packet aggregate. Claude directly spot-checked selected imported
D3 matrix entries against the transcript; neither review certifies the entire
translation generator. The generator was not in the packet. Its existing
provenance/check evidence remains separate from the kernel proof.

No further Lean implementation is required for this contract. Production Q7,
CAS stratum handling, broad comparator/export integration and ledger/paper
reconciliation remain separately open in [COVERAGE.md](COVERAGE.md). S9's
original pilot remains complete within its declared scope. Future work starts
from a new named contract obligation, not the superseded bridge expansion plan.
