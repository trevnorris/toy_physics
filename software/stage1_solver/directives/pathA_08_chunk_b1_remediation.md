# Directive pathA_08 — Chunk B1 remediation (strengthen the validation gates; resolve a schema divergence)

**Owner:** Codex (codes). Claude reviewed B1 (`pathA_07`): the extraction ALGEBRA is term-by-term FAITHFUL (two
independent fidelity reads: eigenproblem vs analytic `κ̂`; coupling/observable vs decision-11 §5 + v2_21/v2_22a
source). The adversarial pass found the **validation gates are weaker than claimed** + one undisclosed structural
divergence. Fix those; do NOT change the (correct) physics algebra. This is the B1 analog of the 1c remediation
(`pathA_05`).

## Background (what's already good — do not regress)
- The 3 genuine physics gates are real and tight; leave them intact: (a) eigenproblem `K/τ→κ̂=7+(10π/37)²` at O(h²)
  over 4 levels with a residual bound; (b) `K`-independence-of-`R0` via a real perturbation; (c) the manufactured
  field→coefficient MMS with hand-derived exacts at O(h²). The v2_21/v2_22a fixture reproductions are genuine
  separate-implementation matches (0.0 diff is legitimate).

## Fixes required

### 1. [MUST] Make the test diff helper unable to silently skip anything
`tests/test_patha_extraction.py` `_numeric_diff` currently recurses over `set(a) & set(b)` (intersection) and
`zip`s lists — so a key present on only one side, or a list-length mismatch, passes vacuously, and string fields
no-op. Fix so the comparison is total:
- At each dict level, assert the key SETS are EQUAL (fail on any symmetric-difference key), then recurse on all keys.
- Before comparing two lists, assert EQUAL LENGTHS, then compare elementwise.
- Add a string arm: assert string equality (not a no-op).
- Acceptance for this item: a deliberately mutated module output — an extra key, a missing key, a changed string,
  or a different list length — MUST make the test fail. State in the report which mutation each new assertion catches.

### 2. [MUST] Resolve the `input_hash` divergence
The module's `extract_branch` drops `input_hash`, which BOTH oracles (`stage_v2_21`, `stage_v2_22a`) return. With
fix #1 enforcing key-equality this becomes a hard failure, which is correct — it surfaces a real schema divergence.
Resolve it the faithful way: have the module emit `input_hash` computed the SAME way the oracles compute it (same
object, same hashing) so the fixture reproductions stay 0.0-diff. If there is a principled reason to omit it
instead, document that reason in code AND exclude it explicitly + symmetrically in the test (not via a silent skip).
Prefer matching the oracles.

### 3. [SHOULD] Actually compare the gate sub-structures
With #1 in place, ensure `open_gate` (the per-component flags `R_exit_positive`, `boundary_class_open_impedance`,
`hard_cap_forbidden`, …) and `name` are inside the compared subtree for BOTH fixture tests, so a per-component gate
bug is caught (today only the `all()` aggregate is checked). The two oracles differ in a couple of keys
(`axisymmetric_P0_pattern_pass` exists only in v2_22a; `constant_prefactor_N2/N4_residual` handling differs) — handle
those explicitly/symmetrically so each is compared against the oracle that has it, rather than relying on the
coincidence of differing key sets.

### 4. [SHOULD] Upgrade the coupling/observable validation from "transcription-faithful" to "independently validated"
The directive `pathA_07` claimed an INDEPENDENT oracle for these layers, but the fixture tests only prove the new
module transcribes the same algebra as the v2_21/v2_22a scripts. The genuinely independent check is those scripts'
own `run_symbolic_audit` (the SymPy series-expansion derivation of `P0/P2/P4` and the `u2/u4` identities). Add a B1
test that exercises an independent symbolic derivation of the observable algebra and asserts the module's numeric
`{D0,D2,D4,P0,P2,P4,R_norm,R_pole}` agree with it on the fixture inputs — either by invoking the oracle scripts'
`run_symbolic_audit` and asserting it passes, or by replicating the SymPy derivation directly against the module's
formulas. This makes a shared *conceptual* error (a wrong formula copied identically into both files) catchable.

## Acceptance criteria (Codex iterates until ALL pass, exit 0; any script `timeout 600`)
1. All pre-existing B1 gates still pass unchanged (eigen/κ̂, K-indep-R0, manufactured MMS, fixture reproductions).
2. `_numeric_diff` now fails on any key-set difference, list-length difference, or string mismatch (demonstrate the
   3 mutation catches in the report, then revert the mutations).
3. `input_hash` divergence resolved (module emits it matching the oracles, OR documented + symmetric explicit
   exclusion); fixture reproductions remain 0.0-diff.
4. `open_gate` per-component flags + `name` are compared in both fixture tests.
5. The new independent symbolic-audit gate passes and would catch a wrong P2/P4/R_norm/R_pole formula.
6. Full `pytest` for the patha_* suite green. Firewall (`research/pde_audit/simulation/`,
   `physical_export_permitted`) untouched. No physical frozen solve / `R_norm(τ)=0` root-find (still B2). No
   git add/commit (orchestrator commits).

## Report back
What changed, the 3 mutation-catch demonstrations (item 2), how `input_hash` was resolved, the new symbolic-audit
gate + what wrong formula it would catch, and the full test results.
