# Independent review — S11c-b ADJUDICATION build directive v3 (decision list, pre-build; round 3)

## Artifact
`research/pde_ledger_v3/directives/S11c_b_adjudication_build_directive.md` — orchestrator-written build brief.
DECISION-LIST review BEFORE the builder launches (no code yet). Two prior rounds removed an ill-posed
computed-collapse gate, a wrong Bridge-A form, a subset-matching container hole, a Boolean-as-algebraic hole,
and repeated rule-5 leaks; v3 restructures around operand-sort invariants + accounting rather than a per-family
schema list. Do NOT assume it is clean — re-derive. The worst failure to catch is any route or fold that
MANUFACTURES false cross-engine agreement.

## Context you are handed (read all of it)
- The directive above (v3).
- v1 reconcile `scripts/S11c_b_handcoded_comparison.py` (committed, legs SOUND) + comparator
  `scripts/S11c_b_cross_engine_comparator.py` (jet decoding via `import S11c_a_cross_engine_comparator as s11ca`;
  `canon_wl_basic`, `FIELD`, `EXTRA_HEAD`, `canon_jet_name`).
- Engine sources `scripts/S11c_b_brane_operator_sympy_audit.py` (PY) and
  `mathematica/S11c_b_brane_operator_mathematica_audit.wl` (WL); spec `directives/S11c_b_SHARED_PHYSICS.md`.
- The committed run `~/.s11_build/S11c_b_reconcile_run.out`.

## Required checks — report a finding only if it changes what gets built or what may be claimed

1. **Bridge A form + reach.** v3 states the SYMBOL substitution `bRho ↦ B_rho_3/W_0` (refs WL:472/1621-1630,
   PY:1130-1140, spec:102). Confirm from source the identity is exactly `B_rho_3 = bRho·W_0` (align the θ²
   coefficients; show them). Then check the substitution is neither too narrow (does it rewrite bare `bRho` and
   `bRho·W_bg`?) nor too broad (could `bRho ↦ B_rho_3/W_0` wrongly rewrite an occurrence that is NOT the
   compression modulus?). ⛔ Reading + algebra; do not spawn Mathematica.

2. **Routing invariant — can a non-arithmetic operand leak into the ALGEBRAIC route, or a CONTAINER MATCH arise
   from a subset?** v3 permits ALGEBRAIC only for arithmetic `sp.Expr`/same-shape arithmetic matrices/tuples,
   and a CONTAINER zero-verdict ONLY via a source-cited TOTAL BIJECTION over ALL leaves. Check: (a) is the
   arithmetic-sort test tight enough that a Boolean (`Equivalent`/`Not`/relational) or `TextAtom` cannot enter
   ALGEBRAIC (see the admissibility operands, PY:868-890, run:160-175)? (b) Can any container path yield a zero
   verdict without a total bijection — subset match, dropped/duplicate/unlabeled leaf, or differently-labelled
   pairing? (c) Could a total-bijection adapter still identify an ENERGY_BASIS quotient representative or a
   PROTECTED 07/10 term?

3. **Deferred-vs-adaptable families.** v3 hands the builder the KINETIC total-bijection adapter as the pattern
   (PY:1573 / WL:851) and requires it to enumerate every source-backed total correspondence and defer the rest
   as STRUCTURE_INCOMPLETE. Check the run's FLAG families: is any family that is genuinely algebraically
   comparable (under a source-fixed total correspondence) at risk of being wrongly deferred (hiding a real
   residual), or any family with NO faithful total correspondence at risk of being force-adapted (manufacturing
   agreement)? Name any family where the boundary is wrong.

4. **Jet conservation.** Semantic jet ID `(canonical_base_after_rename, multiindex)`; the substitution map must
   be jet-order-preserving; `JET_CONSERVED/JET_LOST` computed from before/after multisets. Check: (a) does it
   correctly NOT false-alarm on the legitimate rename `theta_d1→grad_theta_1` (order preserved)? (b) Does it
   still fire on a genuine order-reducing collapse? (c) Does `--collapse-jet` decisively let a leg show a
   background-gradient jet is load-bearing?

5. **One-engine protection.** From PY (~1273) and WL:417-435 confirm selected `{1,4,5,6,7,9,10,13}`, WL explicit
   → PY `{1,4,5,6,8,9,11,13}`, so WL-only DivGrad = 08/11 and PY-only selected = 07/10. Verify nothing in the
   routing or adapters can absorb 07/10 or the DivGrad pair into a match. Flag any unprotected one-engine term.

6. **Leaks (rule 5) + ablation decisiveness + accounting DoD.** Flag any leaked value, measured classification,
   or predicted ablation outcome. Confirm each ablation hook is decisive (unknown/non-occurring argument →
   operational error, nonempty touched set, syntactically changed operand, before/after residuals). Confirm the
   accounting DoD (per-route counts sum to the case count; no silent drop) is value-free and able to FAIL.

## Method
Read the comparator, engine sources, spec, and run FIRST; form your own view; then judge the directive. For
Bridge A show the coefficients. A prose "looks fine" is worth nothing — for each check show what you read/computed
or say you could not verify it. Numbered findings (blocking vs non-blocking), file+line, concrete correction. If
sound, say so and name what you verified.
