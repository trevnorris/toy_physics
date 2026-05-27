---
unit_id: 091
batch: IV.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 091

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
The orchestrator (post-user-resolution per `batch_IV1_paper_alignment.md` Cluster A direction (a)) annotated the SymPy script's module docstring at
`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py:1-16`
with a "Carry-forward annotations" block naming stage 094 (Isotropic Geometry-Decoupling Theorem) as the upstream verifier of the `l=0` / `l=2` orthogonality Check. No script assertion was added, no paper edit, no notes edit.

**Assessment:**
Correct. The directive routed F1 to the user and the user picked direction (a) (docstring carry-forward annotation, no paper-side edit, no new script check). The applied edit is exactly that: a multi-line docstring annotation that explicitly names the upstream stage (094) authorising the scalar-`Kcons(omega)` ansatz this stage uses. No collateral edit. No new identity introduced. Not tautological because it makes no math claim — it is a documentation pointer. The Mathematica `.wl` file has no analogous block, but the cluster directive specifies a single-script docstring annotation, so this is consistent rather than partial.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl:70-78` inserts the prescribed independent partial-fraction recombination block:

```
kConsBranchDirect = 3*kPole + kPole/(1 - omega^2/omegaQ^2);
k0BranchDirect = 4*kPole;
yHatRecomb = Together[kConsBranchDirect/k0BranchDirect];
yHatTargetRecomb = Together[3/4 + (1/4)/(1 - omega^2/omegaQ^2)];
expectZero["Yhat partial-fraction recombination", yHatRecomb - yHatTargetRecomb];
```

The trailing `Print[""]; Print["Stage 091 Mathematica audit passed."]; Exit[0];` (lines 80-83) is preserved unchanged. The banner sweep `STAGE 074 → STAGE 091` mentioned in the `## Applied: F2` block is reflected at line 26 (`banner["STAGE 091 — GROUPED-P2 + STATIC-GEOMETRY DERIVATION"]`).

**Assessment:**
Match against the directive verbatim (modulo a wrapped comment header that does not change semantics). The new check is non-tautological: it builds `kConsBranchDirect` directly from the branch-identity result `K_geom = 3*K_pole` without going through `Series` + `Coefficient` + `Solve`, then verifies `Together[kConsBranchDirect/(4*kPole)] == Together[3/4 + 1/4/(1-omega^2/omegaQ^2)]`. If the bottom-line identity were wrong, `Together` of the difference would not collapse to 0. Independent path achieved.

## Exec log assessment

**SymPy:** exit=0 (run-time inferred from log: no Python exceptions, script completes through `FINAL LEDGER` banner at the tail; `expect_zero` raises `AssertionError` on non-zero, none seen). Notable lines:

- `K_geom - 3*K_pole = 0`
- `Yhat - [3/4 + 1/4/(1-omega^2/OmegaQ^2)] = 0`
- `rho_alpha - 4/3 = 0`
- `zeta_req - 1/3 = 0`

**Mathematica:** exit=0 (`Exit[0]` reached; would `Exit[1]` on any `fail[...]` call). 9 `PASS:` lines (expected: 8 original + 1 new from F2). Notable lines:

- `PASS: K_geom - 3*K_pole`
- `PASS: K0 - 4*K_pole on branch`
- `PASS: Yhat - [3/4 + 1/4/(1-omega^2/OmegaQ^2)]`
- `PASS: rho_alpha - 4/3`
- `PASS: zeta_req - 1/3`
- `PASS: Yhat partial-fraction recombination`   <- new check from F2 present
- `Stage 091 Mathematica audit passed.`

**Output freshness:** SymPy script mtime `May 27 11:12`, SymPy output mtime `May 27 14:28` (output newer). Mathematica script mtime `May 27 11:13`, Mathematica output mtime `May 27 14:29` (output newer). Both outputs were re-generated post-fix. MANIFEST.yaml confirms `status: verified, iteration_count: 1`.

## Material-change assessment

`material_change`: false.

The SymPy edit is a docstring-only annotation — no symbolic content, no assertion change. The Mathematica edit adds a redundant independent check verifying the same identity already established by the existing assertions; it adds zero new derived results. Downstream units cannot depend on intermediates that did not change. No `upstream_stale` propagation warranted on substantive grounds.

## Side observations (non-blocking)

- The F1 docstring annotation appears in the SymPy `.py` only; the Mathematica `.wl` carries no parallel annotation. This is consistent with the orchestrator's stated direction (a) action and is not a missed-mirror issue, but a future "mirror-policy completeness" sweep might consider adding the same docstring to the `.wl`. Non-blocking.
- The Mathematica comment block at lines 70-73 is wrapped differently than the directive's literal text (semantics unchanged). Non-blocking.

## Verdict justification

Both findings are resolved. F1 reflects the user's chosen direction (Cluster A (a)) via a script docstring carry-forward annotation naming stage 094 as the upstream verifier; F2's independent partial-fraction recombination block is inserted exactly where the directive specified (between the `zeta_req` assertion and the closing `Print[""]`), produces a non-tautological check, and shows up in the saved Mathematica transcript with the expected `PASS: Yhat partial-fraction recombination` line. Both engines' saved outputs are fresher than the corresponding scripts; SymPy and Mathematica both reach normal termination with all assertions passing. Verdict: `verified`, `material_change: false`.
