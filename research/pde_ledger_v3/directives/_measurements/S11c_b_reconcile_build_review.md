# S11c-b hand-coded reconcile layer — build review, two legs, consolidated (INSTRUMENT SOUND, MAP COMPLETE)

**Artifact:** `scripts/S11c_b_handcoded_comparison.py` (Codex-built; imports the committed comparator; 80-entry
verified WL→PY rename map + re-check-zero; scoped to the sympy-parsed CORE families; controls →
`NAMESPACE_INCOMPLETE`). **Directive:** `directives/S11c_b_handcoded_reconcile_directive.md`. **Legs:** fresh
Claude agent + Grok (Codex-written ⇒ fresh-Claude + Grok), prompt `_legs/S11c_b_reconcile_build_review.md`.
**Raw:** Grok `~/.s11_build/S11c_b_reconcile_review_grok.txt` + `/tmp/s11cb_handcoded_review/`
(ablation `.stdout` + `REVIEW_REPORT.md`); Claude-agent
`/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/`
(`harness.py`, `probe_match.py`, `completeness_probe.py`, `drop_*.out`, `base_*.out`). Default run: **2 MATCH,
14 FLAG, 0 COVERAGE**; rename-map size 80; run 798s (Codex). Verdicts are the orchestrator's after
adjudication (rule 13); both legs' own verdict was NO FINDING under the false-MATCH filter.

## The ONE thing that mattered — no FALSE MATCH is possible (both legs, ablation-verified)
- **`--drop-rename` is surgical, not blanket** — dropping `forceHoldEw` flips exactly the 4 `DOF=E_W`
  ADMISSIBILITY_SUPPORT cases MATCH→FLAG (THETA/U/traction stay); `sigmaW`/`etaBg`/`WZero`/`muR` each move
  ONLY the residuals of the families where that variable occurs (`muR` on SLAB: `mu_R` 96→0 in
  U_MOMENTUM_ROWS 4/64; a no-op where the variable is absent, legitimately). Each rename is load-bearing =
  the same variable, and the drop punches through BOTH the comparator prepass and the reconcile map.
- **Map INJECTIVITY (the highest false-MATCH risk) verified faithful:** two WL names map to one PY target in
  two spots — `{transverseTrial,transverseTest}→A_T_s11cb_i` and `{eWave,eWField}→e_W`. Both correct: WL's own
  `reverseRelabelRules` substitutes `transverseTest→transverseTrial` before forming its residual (WL
  L1102-1113) and PY uses one symbol `A_T_s11cb_i` (PY L1906); `eWField = eWave[xOne,xTwo,xThree,time]` is
  literally the same field (WL L240). The map deliberately UNDER-folds (keeps theta/e_W/longitudinal
  trial≠test distinct) — the safe direction. No two DIFFERENT quantities fold to one name.
- **The 2 MATCH are genuine:** `ENERGY_BASIS_COUNT` = `Integer(26)==Integer(26)` (rename-independent, survives
  the empty map; a real but low-information ⭐bare-integer agreement); `ADMISSIBILITY_SUPPORT_OPERAND` 20/20 =
  supplied-premise consistency (both engines re-emit the declared `f_hold_*`/`t_hold_*` support bundle after
  spelling reconciliation; breaks surgically when the hold renames are dropped). Neither is an over-broad
  collapse.

## Correctly NOT folded (must stay FLAGGED) — both legs source-verified
- **`bRho` vs PY `B_rho_3`** — a genuine `W_0` SCALE difference (`B_rho_3 = B_rho·W_0`, "uniform integrated
  compression modulus", WL L472 vs PY L133/L1135), NOT a spelling. Folding it would have been a scale-collapse
  false MATCH. Left unmapped ⇒ FLAG.
- **`gamma{Width,Modulus}DivGrad{Theta,Ew}` (4)** — the `DIVU_GRADTHETA`/`DIVU_GRADEW` invariants = PY's
  OMITTED candidates (08/11), not the selected representatives. Correctly unmapped. (The SELECTED gamma
  representatives — `U_THETA→_04`, shear→`_06`, etc. — ARE in the map and source-checked against both engines'
  invariant expressions.)
- **Energy-basis quotient / coupling adjointness** — no representative-fold and no extra adjointness collapse
  in the map; ENERGY_BASIS_* stay FLAG; adjointness uses the comparator's already-IBP-reduced
  `*_DIVERGENCE_REDUCED` operands with no further fold.
- No `assert`, no `PASS/FAIL/VERDICT/target`; no blanket `camelCase→snake_case`. Decides no physics.

## COMPLETENESS — MAP COMPLETE; the 14 FLAGs are genuine, NOT naming noise (both legs)
After renames, sampled FLAG residuals (MU_THETA, ADMISSIBILITY_*, SLAB, COUPLING) contain NO leftover *mapped*
WL keys. Every remaining camelCase atom is one of three LEGITIMATE classes, none a missing spelling pair:
1. intentional physics non-folds (`bRho`, `gamma*DivGrad*`);
2. WL-coordinate-space vs PY-jet-space residue with NO PY sibling (`xOne/xTwo/xThree` — WL carries fields as
   `f(coords,time)`, PY is jet-space with ZERO coordinate symbols; held-operator heads `HeldDiv`,
   `CompactInPlaneSupport`);
3. genuinely different one-engine-only terms already in snake_case — e.g. `ADMISSIBILITY_OPERATOR_OPERAND`'s
   4 THETA cases have WL-only=[] and PY-only = all snake_case (`L_W,W_0,eta_bg,kappa_theta_W,sigma_W,
   w1_profile*`) — PY's first-background-jet admissibility contribution with no WL counterpart.
⇒ The FLAGs are READY for orchestrator physics adjudication; the map needs no second pass.

## DoD note (compute-time, not a soundness hole)
`COUPLING_KERNEL`'s full FLAG residual could not be materialized in a leg's foreground window (its adjointness
`sp.simplify` alone >470s). The false-MATCH question is closed for it STRUCTURALLY (it is a FLAG = the safe
direction; its two risk folds — transverse trial/test, IBP-reduced adjointness — are verified faithful; raw
scan shows no missing spelling pairs). The full committed-layer run (798s, completes) produces the residual
for adjudication.

## Verdict
Reconcile layer PASSES both build legs; instrument SOUND, rename map COMPLETE (verified, source-cited,
ablation-load-bearing, no false MATCH, under-folds in the safe direction, decides no physics). NEXT: commit;
run the committed layer → the 14 FLAG residuals; ADJUDICATE each (rule 13) under the further VERIFIED bridges
(bRho `W_0`-scale, gamma-index, inert-coordinate) — what reduces to 0 = engines AGREE on that object (the
coupling-kernel number); what survives = a genuine cross-engine finding (rule 1/6). ENERGY_BASIS = non-unique
quotient, never fold a representative. Then step record + S11c roll-up card.
