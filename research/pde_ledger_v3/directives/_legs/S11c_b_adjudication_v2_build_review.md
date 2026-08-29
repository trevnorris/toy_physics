# Independent build review — S11c-b ADJUDICATION layer v2 (Bridge D + IBP classification)

## Artifact
`research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py` (Codex-extended, UNCOMMITTED; the v1 layer at
`2e5c6755` was both-legs-SOUND — this reviews the v2 DIFF: `git diff 2e5c6755 -- <that file>`). v2 adds Bridge D
(the engine's `PROFILE_GRADE_SUBS` background expansion) + an IBP/total-in-plane-divergence classifier for the
COUPLING_KERNEL density, strong-operator EXACT comparison, protected-atom gating, and an explicit route order.
Directive: `directives/S11c_b_adjudication_v2_build_directive.md`. Its run reported: `MATCH=38, FLAG=12,
RESIDUAL_BULK=8, REPRESENTATIONAL_DIVERGENCE=0, DIVERGENCE_INCOMPLETE=0, PROTECTED_UNREDUCED=32,
STRUCTURE_INCOMPLETE=57, COVERAGE=84` (231 total); JET 231/0. This layer DECIDES whether the S11c-b operator/
coupling cross-engine differences are representational or genuine findings, so the failures that matter are (i) a
FALSE MATCH / FALSE REPRESENTATIONAL_DIVERGENCE (over-reduction → manufactured agreement) and (ii) a FALSE
RESIDUAL_BULK / FALSE FLAG from a BROKEN classifier or a mis-applied Bridge D (under-reduction → a manufactured
"finding"). Prove which of the two the `REPRESENTATIONAL_DIVERGENCE=0` is.

## Method + constraints
- Derive/ablate independently; SAVE every script + literal stdout to named absolute paths and report them. Prose
  "I checked" is discarded.
- ⛔ Copy anything you ablate to /tmp; NEVER modify the working tree / commit / `git stash`. Pure Python (reads
  committed `.out`). ⛔ No Mathematica.
- ⭐ Run in the FOREGROUND; do NOT background or spawn a monitor loop. The full layer run is ~11-18 min — the
  harness's native background+notify is acceptable ONLY if you cannot fit a foreground call, but never a
  setsid/&/poll loop. Do not stop until your full verdict is written.

## Required checks — a finding only if it catches a real way the comparison could be wrong

1. **The divergence classifier — is `REPRESENTATIONAL_DIVERGENCE=0` REAL or a broken-classifier false negative?**
   This is the crux. (a) Confirm the 4 fixtures pass (product-rule bulk→RESIDUAL_BULK; `a_d1·φ+a·φ_d1`→
   REPRESENTATIONAL_DIVERGENCE with verified V; strong-family divergence→ineligible; the anchored non-commutation
   certificate). (b) Then take a REAL coupling residual the layer classified `RESIDUAL_BULK` and INDEPENDENTLY
   attempt to reconstruct a flux `V` with `R − Σ∂_iV_i == 0` (the coupling residual has the curl-like structure
   `A_T·coef·(∂_2 v·∂_3 w − ∂_3 v·∂_2 w)`; try `V` built from the operands). If a real coupling residual IS a
   total divergence but the layer said RESIDUAL_BULK, that is a broken classifier (FALSE finding) — BLOCKING. If
   you also cannot find a V, corroborate RESIDUAL_BULK. (c) Ablate: construct a KNOWN real-scale total divergence
   in the coupling family and confirm the classifier finds its V (proves it is capable at scale, not just on
   toy fixtures).

2. **No FALSE MATCH (v1 property preserved).** `--drop-bridge-a` (38 unchanged? they carry no bRho),
   `--drop-rename` surgical, `--collapse-jet` → JET_LOST. Confirm the 38 MATCH are still the rename-level
   admissibility/count agreements and Bridge D created none by over-reduction.

3. **Bridge D correctness + non-commutation.** Confirm Bridge D references the imported `PROFILE_GRADE_SUBS`
   (all 12 entries; `sigma_W ≢ W_0·eta_bg`), retains all jets (JET 231/0). Confirm HeldDiv on strong rows is
   EXPANDED (`Σ formal_∂_i(V_i)`) NOT dropped (`--drop-bridge-d` resurfaces the residual; check no
   `_drop_held_divergences`/`reduce_divergence=True` on strong families). Verify the V-certificate uses
   `BridgeD(R_preD − Σ∂_iV_i)`, NEVER `∂_i(BridgeD(V_i))` — and that raw `case.value` is consumed, no
   `DIVERGENCE_REDUCED` context reaches the classifier.

4. **The 12 FLAG are genuine strong-operator differences.** Pick 2-3 and independently reduce operand_A − B
   under {renames, Bridge A, Bridge D, HeldDiv-expand, integral linearity} only; confirm a genuine nonzero bulk
   remains (not an artifact of a mis-applied Bridge D or an un-expanded HeldDiv).

5. **Protected + accounting.** The 32 `PROTECTED_UNREDUCED` each carry a protected atom (`07/10`, gamma-DivGrad);
   none were Bridge-D-folded or divergence-reduced; ENERGY_BASIS never entered Bridge D/divergence;
   `ENERGY_BASIS_COUNT` compared exactly. Case-ID multiset == emitted `join+py_only+wl_only` (231); duplicate-key
   raises; ablations exit 2 on bad args. Hygiene: no assert/PASS/FAIL/VERDICT/target on measured payloads.

## Report
Numbered findings (BLOCKING vs non-blocking), each with file+line, the ablation/reduction command + its literal
stdout path, and a concrete correction. Overall verdict (SOUND / NOT SOUND). CRITICALLY: state your independent
determination of whether `REPRESENTATIONAL_DIVERGENCE=0` is a real result (coupling differences are genuinely
bulk) or a broken-classifier false negative, and whether the 8 RESIDUAL_BULK + 12 FLAG are genuine.
