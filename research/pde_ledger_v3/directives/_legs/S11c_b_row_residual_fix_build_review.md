# Independent build review — S11c-b row-residual instrument, coupling weak-verdict FIX

You are reviewing a **Codex fix** to a physics CAS instrument. A prior build was leg-verified sound EXCEPT for
the coupling verdict path; this fix reworks only the coupling classification and the certificate-witness
determinism. Verify the fix is correct AND that nothing else changed. A prose re-derivation is worth nothing —
**run ablations and report literal stdout**. Work on `/tmp` copies; never modify the working tree. A leg that
finds nothing is weak evidence — probe hard.

## Artifact
`research/pde_ledger_v3/scripts/S11c_b_row_residual.py` (run:
`python3 research/pde_ledger_v3/scripts/S11c_b_row_residual.py`; emits ~16MB, redirect to a file). Imports the
committed layer `scripts/S11c_b_adjudicated_comparison.py` (do NOT modify it) and reads the committed
transcripts under `scripts/out/` and `mathematica/out/`.

## What the fix must have done (the reviewed construction — verify each step is implemented, on the PRE-bridge residual)
For each coupling case (both cross-sector blocks + both relabelled adjoints), with `A =
S11c_b_adjudicated_comparison`, `T = requested_truncation`, `R_pre = A._arithmetic_residual(left_pre_bridge_d,
right_pre_bridge_d)`:
- `FULL_PREBRIDGE_ROUTE`, `EULER_SIGNATURE` from `A.classify_total_divergence(R_pre, A.PRODUCTION_FIELD_REGISTRY,
  apply_bridge_d=True)` — the full-order classification runs on the **pre-bridge** residual (independent jets).
- `ROW_RESIDUAL = T(A._bridge_d(R_pre))` — the requested-truncated bridged cross-engine **A−B density** (A, B,
  A−B all emitted).
- `IN_SCOPE_WEAK_REMAINDER = T(A._bridge_d(A._normalise_exact(R_pre − A.formal_divergence(A._homotopy_vector(
  R_pre, A.PRODUCTION_FIELD_REGISTRY), A.PRODUCTION_FIELD_REGISTRY))))` — the homotopy remainder computed **in
  the instrument** so it exists even when the layer short-circuits (`RESIDUAL_BULK`); on `_homotopy_vector`
  throw, a `NO_CLEAN_QUOTIENT` flag + fall back to `ROW_RESIDUAL`.
- The Euler / homotopy operators must run on the **pre-bridge** residual; the truncation on the **bridged**
  result; the witness canonicalized deterministically (presentation only, never applied to
  `ROW_RESIDUAL`/remainder/V).

## Independent derivation FIRST (save script + literal stdout to named /tmp paths)
Confirm, with your own small SymPy check, that `_bridge_d` does NOT commute with `formal_divergence` (so running
Euler on the bridged residual is wrong), using the layer's fixture or a `div(W_bg u_T)`-type example: show that
`div(W_bg u_T)` classifies REPRESENTATIONAL on the pre-bridge residual but produces a spurious nonzero Euler
signature on the bridged residual. This is the defect the fix must remove.

## MANDATORY form ablations (report literal before/after stdout; a FORM change tests physics)
1. **Pre-bridge vs post-bridge Euler.** On a /tmp copy, change the coupling classification input from the
   pre-bridge `R_pre` to the post-bridge (bridged) residual (and/or flip `apply_bridge_d`). The
   `FULL_PREBRIDGE_ROUTE` / `EULER_SIGNATURE` for at least one coupling case MUST move (the wrong convention
   manufactures false bulk). If nothing moves, the fix did not actually switch to the pre-bridge convention.
2. **Homotopy divergence-peeling.** On a /tmp copy, skip the `R_pre − div(V)` step (set `IN_SCOPE_WEAK_REMAINDER
   = T(bridge(R_pre))` — i.e. do not peel the divergence). Construct a synthetic coupling residual that is a
   pure in-plane divergence with genuine bulk only at `η²` (e.g. `div(W_bg u_T) + (W_bg−W_0)² u_T`): with
   peeling, its `IN_SCOPE_WEAK_REMAINDER` MUST be 0; without peeling it is nonzero. Report both. (If you cannot
   inject a synthetic case cleanly, at minimum confirm the two emitted objects `ROW_RESIDUAL` and
   `IN_SCOPE_WEAK_REMAINDER` differ for some real coupling case, proving the peeling is wired.)
3. **Determinism.** Run the instrument twice (or the coupling path twice) and diff the emitted `EULER_SIGNATURE`
   witnesses byte-for-byte — they MUST be identical. Then confirm the canonicalization is NOT applied to
   `ROW_RESIDUAL`/`IN_SCOPE_WEAK_REMAINDER` (those must be the raw computed expressions).

## Regression: non-coupling output unchanged
Diff the full stdout of the fixed build against the pre-fix build (or against a re-run) restricted to the
strong slab rows, mass row, MU_THETA, and admissibility emissions — they MUST be **byte-identical**. Any change
outside the coupling family is a regression. Report the literal diff summary.

## Probes (report the line number that computes each, or report uncomputed)
- No `assert` before any residual/route emit; emission conditional only on row/case/engine/block, never on a
  value. No expected route/residual hardcoded.
- The homotopy remainder is computed on the **pre-bridge** residual (not the bridged one); the truncation is on
  the bridged output.
- The `NO_CLEAN_QUOTIENT` throw fallback exists and is exercised or unreachable (state which).
- ⛔ The committed layer `S11c_b_adjudicated_comparison.py` is not modified (diff it against HEAD).

## Physics filter / sandbox
Report a finding only if it catches a way the coupling verdict is computed wrong, a regression outside coupling,
or a leaked/hardcoded answer. Copy to /tmp and ablate the copy; save each ablation script + stdout to named
/tmp paths and report them.
