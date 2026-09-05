# S11c-c1 WL repair-2 directive — decision-leg review (rule-7 gate, BEFORE the build)

This is the gate the prior repair directive skipped. Here it is run **before** the repair-2 build, as rule 7
requires. Artifact reviewed: `directives/S11c_c1_wl_repair2_directive.md` (draft, pre-fold). Prompt (committed):
`directives/_legs/S11c_c1_wl_repair2_directive_review.md`. Both legs: orchestrator-written directive → **Codex +
Grok**.

## Commands (literal)

```bash
# Codex leg
codex exec -c model_reasoning_effort=xhigh "$(</var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_c1_wl_repair2_directive_review.md)" \
  > .../scratchpad/codex_repair2_dir.log 2>&1        # exit 0; 112,260 tokens

# Grok leg
grok --prompt-file /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_c1_wl_repair2_directive_review.md \
  --cwd /var/projects/toy_physics --model grok-4.6 --effort high --permission-mode bypassPermissions \
  --output-format plain > .../scratchpad/grok_repair2_dir.log 2>&1        # exit 0
```

## Verdicts (literal key lines)

Both legs: **do not build against this directive** (convergent). Neither leg was given the SymPy engine, step
records, or prior-art; each formed its view from the engine code + spec.

- **Codex** (`codex_repair2_dir.log`): "Decision: **do not build against this directive.** I found three MUST-FIX
  issues." — (1) "R2 is not a real defect; the prescribed rewrite would corrupt a correct port-basis
  transformation" (the congruence `Aᵀ·diag(P₊,P₋)·A`, `A=faceToParityMatrix=[[1/2,1],[1/2,-1]]`, given
  `V_s=s∂_tζ_s`); (2) "completeness misses the actual additional parity-emission defect" —
  `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (`.wl:1713-1736`) emits `portCase["HERMITIAN_FORM"]` (per-face) under both
  `PARITY_*` keys, parity iterator unused; (3) "R1 contains a forbidden lazy escape through a 'full-operator
  diagonal symbol'" — `k=k′` selects a diagonal slice, a fresh dummy renames the freeze; remove it + the +1/+2
  tag allowance.
- **Grok** (`grok_repair2_dir.log`): "**Not sound. Do not build against this directive.** R1 is real and the
  prescribed two-leg object is right. R2 is not real *as described*, and the fix it names would attach `DELTA_W`
  to the wrong object." Decisive disproof: "if `M₊ = M₋ = M` (flat order, equal faces), a thickness drive still
  sees port `M`. The directive's `DELTA_W = (plus−minus)` is then **0** — the thickness port vanishes." Both
  `DELTA_W`/`ZETA_C` blocks are even combinations (congruence); the odd combination is the coupling — "not a
  swap." Grok NIT on R1: the re-freeze control should be constructor-level (`rightLegMomentum → momentumOutput`,
  as `.wl:1329-1332`), not a post-hoc substitution (the `Count` probe is syntactic). Grok also lists
  `UNIFORM_LIMIT_*` (`.wl:2017-2039`) as the same dead-parity-axis class on a flat object.

## Reconciliation (verified against code myself — rule 13)

The two legs **converge**, and they **correct the retroactive finding's location**:

1. **R1 REAL** (both legs). `fredholmFunctionSpaceOperator` `.wl:580-597`: both `nZero`/`gZero` on `momentumOutput`,
   `momentumInput` absent, profile source dropped; consumed at `.wl:1537-1549`. Fix = the two-leg `Z`; only
   legitimate diagonal symbol is the already-emitted flat one; re-freeze is constructor-level. Verified.
2. **R2 (as I wrote it) WRONG** (both legs). `PERMEABLE_PORT_HERMITIAN` is **correct**: `deltaWHermitian =
   (plus+minus)/4`, `zetaCHermitian = plus+minus`, coupling `(plus−minus)/2` are exactly the congruence
   `Aᵀ·diag·A` blocks under `A=[[1/2,1],[1/2,-1]]`, correct given `V_s=s∂_tζ_s` (`S11c_a:58`). The naive
   `DELTA_W=(ζ₊−ζ₋)` is a category error (matrices ≠ face values) that vanishes the thickness port at equal faces.
   **Building my original R2 would have corrupted a correct object.** ⇒ the retroactive plan's confirmed-defect #3
   ("PORT_HERMITIAN swapped") is **refuted and relocated**.
3. **The real propagated parity defect** = the **dead parity axis** in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (curved
   per-face memory Hermitian form under both parity keys; should be the parity-combination — spec §3b `:308-320`).
   `UNIFORM_LIMIT_*` is the same pattern on a flat (parity-independent) object — a redundant axis, NIT-class.

## Fold (one pass, then go — rule 7)

- **R1**: keep; require the two-leg `Z` (reuse `operatorCompositionFromDerivation`); **remove** the diagonal-symbol
  escape + tag allowance (tag count stays 51); re-freeze control **constructor-level**; make `.wl:577-579` comment
  true or remove.
- **R2**: **replaced** — drop the PORT_HERMITIAN "swap" fix (return it to byte-identical), fix the dead parity axis
  in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (emit the parity-combination via the existing change-of-basis; distinct
  per parity; one-sided per-face corruption moves them differently); `UNIFORM_LIMIT_*` as a lighter NIT preserving
  `UNIFORM_LIMIT_RESIDUAL`.
- **Deferred (not folded):** PORT_HERMITIAN's blocks are hand-written to match the congruence rather than computed
  from the stored `portParityTransformation` (Grok control-strength note), and the similarity-vs-congruence
  consistency with `FACE_RESPONSE` parity — to the step record / T7 cross-check.

Fold applied to `directives/S11c_c1_wl_repair2_directive.md`; re-leak-gated clean; **not re-legged** (rule 7 one
pass). NEXT: detached Mathematica build → 2 re-review legs.
