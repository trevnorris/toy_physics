# P2b — §3a coefficient bridge (expression-valued scale + bijection certificate) — decision list v2

⚠ **v2 — folded BOTH decision legs (Codex + Grok, computation-backed, convergent), rule 7 one pass.** Record
`directives/_measurements/S11c_b_p2b_gamma_bridge_directive.md`. Targets: `scripts/S11c_b_handcoded_comparison.py` +
`scripts/S11c_b_adjudicated_comparison.py` (+ `scripts/S11c_b_cross_engine_comparator.py` `_extra_basic` for B4).
Physics-bearing — the manufactured-agreement danger zone. Both engines now enumerate the same 15/source
O(3)-Kronecker §3a family (#89/#89a) but with different NAMES and different NORMALIZATION. ⚠ Sequencing: this lands
BEFORE P2a's final `row_residual` validation (P2a's parse passes through B4's collapse). ⛔ never blanket-collapse.

## The measured facts this must honor (both legs, orchestrator rule-13 spot-checked)
- `I_PY = W_0·I_WL` (W-family) and `I_PY = μ_R·I_WL` (μ-family): `EXACT_UNIQUE 0 / SCALED_UNIQUE 30`. WL spurion
  `∇W/W_0`, `∇μ/μ_R` (`…audit.wl:721-722`); PY raw jets `grad_W`/`grad_mu` (`sympy_audit.py:182`). ⇒ energy terms
  `γ·I` equal ⟺ `γ_WL = W_0·γ_PY` (resp. `μ_R·γ_PY`).
- The committed bridge apply path can ONLY emit a `Symbol` (`_reconcile_basic` → `sp.Symbol(active_map[name])`
  `handcoded:279/310`; `_extra_basic` → `sp.Symbol(target)` `comparator:396`). Writing `"W_0*gamma_…"` into the
  string table makes a NEW ATOM, not a product. The mechanisms that DO express a scale as `Mul(scale, symbol)` are
  Bridge A (`adjudicated_comparison.py:85-94`) and `s11ca.field_symbol` (`S11c_a_cross_engine_comparator.py:818-823`).
- Positional emission order CURRENTLY coincides with the invariant pairing (`POSITIONAL_EQ_SCALED True`), and WL does
  NOT emit in quotient order (`enumerateO3KroneckerBlueprints` + `DeleteDuplicatesBy`, `…audit.wl:540-616`) — so a
  positional table passes today but is not the pairing.
- `ENERGY_BASIS_NEW_INVARIANTS` is routed to `STRUCTURE_INCOMPLETE` (`adjudicated_comparison.py:67-72,927`); the
  §3a coefficients live inside the operator/kernel rows.
- The applied→bare collapse `widthBackground[x] − widthBackground[chi] → 0` happens at the FULL residual parse
  (`s11ca.canonical_value` keeps args since `widthBackground ∉ s11ca FIELD`, then `_extra_basic`/`BARE_APPLIED` drop
  them); same for `modulusBackground`, `eWBackground`. `thetaWave` is stripped by s11ca canon (a DIFFERENT,
  legitimate S11c-a fold — must not be retargeted). `INERT_APPLIED` today = only `{mu_theta_L, mu_theta_M}`.
  Jet decode is gated on `BARE_APPLIED` membership (`handcoded:274`).

## What must be TRUE
B1. **Expression-valued scale bridge.** Each WL coefficient's applied IMAGE is the EXPRESSION `W_0·gamma_s11cb_w_bg_NN`
    (W-family) / `μ_R·gamma_s11cb_mu_r_bg_NN` (μ-family) — a `Mul(scale, PY_symbol)`, using an expression-valued
    substitution (the Bridge-A / `field_symbol` mechanism), ⛔ NOT the string rename table (which can only make a bare
    atom; a string rename leaves a stray `1/W_0` = manufactured disagreement). ⛔ The scale is EMITTED as the family
    factor; ⛔ NEVER folded out of `I_WL` (that fold is the blanket collapse — and it must be made UNREPRESENTABLE,
    not merely forbidden in prose). ⛔ Not `COEFFICIENT_STANDARD_NAME` (a disjoint no-op name).
B2. **Pairing by the invariant OBJECT, as a METHOD (not a zip).** Each WL coefficient pairs to the PY coefficient
    whose UNFOLDED invariant satisfies `I_PY = factor·I_WL` (factor the B1 family symbol), matched under the
    already-complete profile-jet/DOF renames. ⛔ NEVER positional — the directive must NOT claim WL emits in
    "quotient order" (no such object); a PERMUTATION of either engine's emitted record order must leave the computed
    pairing + certificate UNCHANGED (a permutation-invariance control), and duplicate emitted coefficient keys are
    rejected.
B3. **Bijection CERTIFICATE bound to the RESIDUAL apply path.** Walking the `ENERGY_BASIS_NEW_INVARIANTS` records
    (un-gate them from `STRUCTURE_INCOMPLETE` as needed), the certificate PRINTS per coefficient, per branch: the two
    coefficient names; the two UNFOLDED invariants `I_WL`, `I_PY`; the family factor LOCKED ∈ `{W_0, μ_R}` by source
    (⛔ computed AS that symbol — not `I_PY/I_WL`, not `1`); the invariant residual `I_PY − factor·I_WL`; AND the
    ENERGY-TERM residual `γ_PY·I_PY − subst(γ_WL)·I_WL` under the SAME expression substitution the operator/kernel
    residual applies (⛔ not a sidecar). Bijection is checked BOTH directions: every WL degree 1, every PY reverse
    degree 1, 15/15 per source, 30/30 per branch, 60/60 overall. ⛔ Decisive failures (verified by both legs): a
    SWAP → invariant residual ≠ 0; a FOLD of `W_0` / a missing coefficient scale → term residual ≠ 0 OR the factor
    lock fails; a positional mismatch under permutation → the pairing changes.
B4. **Argument-sensitive applied→bare, on the FULL parse, pinned to the S11c-b tables.** A fixture through the FULL
    residual parse (`canonical_value` → `_extra_basic` → `_reconcile_basic`) asserts base-vs-shifted stays DISTINCT
    for `widthBackground` AND `modulusBackground` AND `eWBackground` (all live-background anchors). Fix the S11c-b
    tables (`_extra_basic`/`EXTRA_HEAD`, handcoded `BARE_APPLIED`) — ⛔ NOT `s11ca.canonical_value`/`FIELD` (that
    strips `thetaWave` legitimately and changing it regresses S11c-a). Classify EVERY `BARE_APPLIED` head: inert
    evaluation (removable), live-background anchor (args KEPT + argument-sensitive), trial/test/support with a bare
    PY image (args removable). ⛔ Keep the jet decode working (it is gated on `BARE_APPLIED` membership at two sites).
    A positive fixture proves an explicitly inert head still reduces.
B5. **Retire the stale protection.** `PROTECTED_ATOM_NAMES` drops the dead `gammaWidth*`/`gammaModulus*`; re-adjudicate
    PY `07/10` (both engines carry all 15/source now, so a bridged coefficient must not stay protected — protection
    diverts a nonzero residual to `PROTECTED_UNREDUCED` `:953` and excludes it from ablation-touch `:1478`, hiding
    disagreement). B3 exercises `07/10` with a planted residual → it must NOT route to `PROTECTED_UNREDUCED`.
B6. Output states operands + residual; no verdict token, no assert on the residual (rule 5). S11c-a NOT regressed
    (B4 pinned to the S11c-b tables; the s11ca canon/`FIELD` folds are untouched).

## Legs
Codex build ⇒ 2 build legs (fresh Claude agent + Grok). Mandatory FORM ablations that must BITE (as OBJECTS, per B3/
B4): a wrong pairing → B3 invariant residual ≠ 0 with the factor locked; fold `W_0` out (factor→1, string rename) →
B3 term residual ≠ 0 or factor-lock fails; a base-vs-shifted pair through the FULL parse → nonzero for width AND
modulus AND eW; leave `07/10` protected → a planted residual routes to `PROTECTED_UNREDUCED`; permute an engine's
emission order → the computed pairing is unchanged.
