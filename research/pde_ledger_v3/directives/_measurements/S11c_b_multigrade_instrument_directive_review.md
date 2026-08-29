# S11c-b multigrade-instrument directive — DECISION-LIST review, two legs, consolidated (folded once)

**Artifact reviewed:** `directives/S11c_b_multigrade_instrument_build_directive.md` (orchestrator-written ⇒ two
legs = Codex + Grok). **Leg prompt:** `directives/_legs/S11c_b_multigrade_instrument_directive_review.md`.
**Raw leg outputs (outside repo, tree hygiene):** Codex
`~/.s11_build/S11c_b_multigrade_dirreview_codex.txt` (116,331 tok), Grok
`~/.s11_build/S11c_b_multigrade_dirreview_grok.txt`.
**Ground truth handed:** the layer `scripts/S11c_b_adjudicated_comparison.py`, comparator
`scripts/S11c_b_cross_engine_comparator.py`, engine `scripts/S11c_b_brane_operator_sympy_audit.py`, spec
§2a/§3a/§3d, and the fresh adjudication run `~/.s11_build/S11c_b_adjud_fullres_run.out` (601.88s clean run;
`ROUTE_ACCOUNTING MATCH=38 FLAG=12 RESIDUAL_BULK=8 PROTECTED_UNREDUCED=32 STRUCTURE_INCOMPLETE=57 COVERAGE=84`,
`EMITTED=CLASSIFIED=231`, `CASE_ID_MULTISET_EQUAL=true`, `JET 231/0` — reproduces the v2 build review exactly).

Both legs read the real `srepr` operands (not prose) and verdicted **NOT SOUND to build as written**. Six
distinct defects, all verified by me against the code (rule 13); the two legs independently converged on the
rename-map and kinetic-alignment defects. Folded ONCE into the directive (rule 7: one leg pass → fold → go).

## Findings (verified) and the fold

1. **Wrong alignment path (Codex #1, Grok #2).** Directive said `active_names = dict(C.…WL_TO_PY_RENAME)` and
   `transform(bridge_d=True)` then `_bridge_d`. VERIFIED: `WL_TO_PY_RENAME` is defined in
   `S11c_b_handcoded_comparison.py:72`, re-exported as `S11c_b_adjudicated_comparison.py:35`
   (`= H.WL_TO_PY_RENAME`); comparator `C` has **0** occurrences. `main` (`:1356,:1403`) uses
   `active_names = dict(WL_TO_PY_RENAME)` and `transform(..., bridge_d=False)` then `_bridge_d` **once**.
   FOLD: prescribe `A.WL_TO_PY_RENAME`, `transform(bridge_a=True, bridge_d=False, collapse=None)`,
   `A._bridge_d(...)` once.
2. **Kinetic residual + leaf alignment under-specified (Codex #2, Grok #1, Grok #3).** VERIFIED on the real
   KINETIC operands: SymPy `left` is a nested tuple (paths `0/0..0/2`, `1`), Wolfram `right` is a labelled
   `Association` (`U_MOMENTUM_ROWS/…`, `THICKNESS_ROW`); they share no container path. `_kinetic_residual`
   consumes `_kinetic_pairs(family,key,left,right)` output, not raw operands. Grading A/B/residual by native
   container path would emit three incomparable path schemes and void `GRADE_DIFFERENCE` on the 4 kinetic
   cases. FOLD: use `A._kinetic_pairs(...)` (fail if `None`) → `A._kinetic_residual(pairs)`, and key all three
   kinetic multigrades by the SAME `_kinetic_pairs` semantic labels; acceptance now enforces this.
3. **Grading symbols not pinned (Codex #3).** VERIFIED: `P.eta_bg`, `P.sigma_W` are `Symbol(real=True)`;
   `sp.Symbol("eta_bg") == P.eta_bg` is `False`. A builder using plain symbols mis-grades every bookkeeper
   occurrence into the `(0,0)` cell while the guards still pass. FOLD: bind `eta_bg=P.eta_bg`,
   `sigma_W=P.sigma_W`; assert identity at import; require every coefficient free of those two atoms.
4. **Finite-`sp.series` incomplete + nondeterministic (Codex #4; Grok confirms rational leaves).** VERIFIED:
   `COUPLING_KERNEL` operands carry `1/(1 + eta_bg*w1_profile)` (Poly fails; formal Taylor order unbounded);
   `sp.series` truncation-order semantics are ambiguous. FOLD: replace with EXACT derivative-coefficient
   `c[a,b] = ∂ᵃ_η∂ᵇ_σ leaf|₀/(a!b!)` on a fixed inclusive window `[0,N]²` (`N=4`, above the observed ≤2
   degree), plus an EXACT remainder `R = cancel(together(leaf − Σ_window))` and R's computed leading grade —
   nothing dropped, fully deterministic.
5. **Acceptance can pass a wrong extraction (Codex #5, Grok #3).** VERIFIED as design logic: guards required
   *printing*, not normalising to zero; a builder could set `REMAINDER=leaf` with no coefficients, or grade in
   wrong variables, and satisfy every clause. FOLD: strengthen to mechanical exact-zero on `RECONSTRUCTION`,
   `WINDOW_CLEAN` (R has no window content — anti-gaming), `GRADE_DIFFERENCE`, `REMAINDER_DIFFERENCE`; plus
   coefficient-symbol-freeness and kinetic leaf-path-set match. Pure algebraic self-consistency ⇒ no rule-5 leak.
6. **Output key/nesting format under-pinned (Grok #4).** VERIFIED as a determinism defect: nested Association
   vs flat composite key, and `(0,1)` vs `"0,1"` vs `"eta_bg**0*sigma_W**1"` spellings — two faithful builders
   emit incomparable stdout. FOLD: one nesting convention (leaf-path → grade-pair) and one grade-key spelling
   (`"a,b"`), leaf path = `"ROOT"` (scalar) or the `_kinetic_pairs` label (kinetic).

## Checked clean by both legs (no finding)
Twenty-case set exact and `DIVERGENCE_ROUTE_FIXTURE` correctly excluded; the reuse path is the routed one and
none of the twenty families are energy-gated; the object (`(a,b)↦c[a,b]` in `eta_bg,sigma_W`) is the right
measurement for the downstream background-order adjudication; no expected grade population / engine verdict /
coefficient target leaked; `RECONSTRUCTION`/`GRADE_DIFFERENCE` are decomposition self-checks, not physics
acceptances.

## Disposition
Directive folded once (all six). Leak-gate re-run clean (only the rule-5 prohibition text and guard mechanics
match probes). NEXT: Codex builds `scripts/S11c_b_background_multigrade.py` (Codex-written ⇒ build legs =
fresh Claude agent + Grok). No re-leg of the directive (rule 7).
