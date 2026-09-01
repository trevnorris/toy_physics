# S11c-b #89 PY §3a-repair engine — build-review clearance (2 legs, both CLEAR)

Artifact reviewed: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` at the
pre-migration checkpoint commit (`f655ea65`, rewritten from `fce14c1a`) + its annexed output
`research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out` (basis 40).

## Legs (rule 7: engine is Codex-written + orchestrator lever-C ⇒ legs = fresh Claude agent + Grok)
- **Leg A — fresh Claude general-purpose agent.** VERDICT CLEAR. Evidence dir
  `~/.s11_build/S11c_b_89_buildleg_claude/`: saved `.py`+`.out` for A1 (A1_independent_basis_count +
  A1_form_ablation_driver + engine_ablated.py), A2/A3/**A4** (A2_A3_paths — carries the σ_W 2nd-jet), and
  B5/B6/**B7** (B5_bounded [carries the `srepr-formdiffs=17` reduced-form result], B5_B6_leverC,
  B6_fallback_ordering, B6_ordering_small). ⚠ Leg A saved NO separate C8 script — it reported C8 soundness
  in its VERDICT only; C8 is fully evidenced by Leg B (`C8_soundness.py/.out`) with Leg A concurring.
- **Leg B — Grok (grok-4.5, high).** VERDICT CLEAR. Evidence dir `~/.s11_build/S11c_b_89_buildleg_grok2/`
  (A1_independent_basis_count, A1_form_ablation + engine_ablated.py, A2_A3_A4_paths, A3_material_pullback,
  A3_operator_from_density, B5_leverC_numeric, B6_fallback_ordering, B7_reduced_form_comparator, C8_soundness).
  ⚠ Grok's FIRST attempt was rejected (bare CLEAR, no saved script+stdout, part B abandoned — a bare Grok
  clear is weak evidence); a second attempt died spuriously at startup (not OOM); this third run, under a
  strict completion+evidence contract (`_legs/..._grok2.md`), is the leg of record.
- Prompts: `directives/_legs/S11c_b_89_sympy_3a_repair_build_review.md` (both legs) +
  `directives/_legs/S11c_b_89_sympy_3a_repair_build_review_grok2.md` (Grok re-run contract).

## Physics claims and the computations behind them (rule 2 — command + literal output)
Both legs derived every count/identity from first principles (own scripts), then matched the engine; the two
engines are independent (fresh Claude + Grok) and agree.

- **§3a basis is 40 = 10 + 15(∂W_bg) + 15(∂μ_R,bg), reduces to frozen 26.**
  Independent from-scratch enumeration: LIVE spurion carried WITH its own first jet (Hessian) in the
  EL/divergence quotient → rank/new-invariants **15 per source, nullity 0**; FROZEN (∂g≡0) → **8 per source,
  nullity 7**. Totals 10+15+15 = **40** live, 10+8+8 = **26** frozen. (Leg A `A1_independent_basis_count.out`;
  Leg B `A1_independent_basis_count.out`.) Engine emits `S11CB_ENERGY_BASIS_COUNT … Integer(40)` for both
  LAB_HELD and MATERIAL_ADVECTED, `NEW_INVARIANT_COUNT = W_BG:15, MU_R_BG:15` (committed `.out`).
- **The repair is load-bearing (MANDATORY form ablation).** In a `/tmp` copy with the Hessian /
  second-background-jet map zeroed (`engine_ablated.py`), the LIVE basis drops **40 → 26** and the ablated
  strong rows become structurally identical to the frozen rows (residual 0 for U/THETA/E_W). The engine's own
  frozen-depth switch also yields 26. (Both legs, `A1_form_ablation*.out`.)
- **Spurion enters the divergence map, not the variation-field set.** `basis_euler_signatures` fields =
  {bu_1,bu_2,bu_3,btheta,be_local}; `bg` ∉ fields (`any bg in fields → False`, L1294); it enters only via
  `derivative_maps[i][bg[j]] = Hessian`. No `δ/δ(spurion)` EL variation is taken. (Both legs, A2.)
- **All four frozen paths repaired.** `total_derivative(grad_W[0],·,depth=2) = σ_W·w1_profile_d1d1/L_W`
  (depth-1 → 0); coupling cascade depth-3 generates a genuine THIRD jet `σ_W·w1_profile_d1d1d1/L_W²`;
  `material_pullback` depth-2 minus depth-1 carries 2nd-jet atoms (6); `operator_from_density` E_W/THETA
  Hessian-atom counts 12, U 8. (Both legs, A2_A3*, A3_material_pullback, A3_operator_from_density.)
- **Generated higher jets are σ_W¹ and retained** (not truncated by the η≤1 ∧ σ_W≤1 filter); σ_W² and η²
  drop to 0. The emitted SLAB_OPERATOR retains 2nd-jet atoms. (Both legs, A4.)

## Lever-C fast projector (performance change) is value-exact, reduced-form only
- Numeric physics-identity: `reference(x) − fast(x) == 0` for EVERY sampled operator+kernel scalar —
  Leg A 262/262 (max|Δ|=0), Leg B 190/190 (nonzero_bugs=0). (`B5_bounded.out` / `B5_leverC_numeric.out`.)
- Form-only difference: **17** scalars differ in srepr (reduced `X/D` vs reference `2X/(2D)`) but are
  value-identical; the within-engine comparator is value-based (`cancel(together(expand(·)))` / `.equals`),
  so the reduced form cannot shift a comparison — `exact_rational_residual(2X/(2D), X/D)=0`. (Both legs, B7.)
- Integral/Derivative-bearing scalars route to the exact reference (srepr byte-identical); serial vs parallel
  `retained_grade` reassembly residual 0. (Both legs, B6.)

## Soundness (C8) — Leg B saved (`C8_soundness.py/.out`); Leg A concurred in its verdict (no saved C8 script)
Primary `count = len(basis_rows)+len(new_rows)` is computed, not asserted; no `Integer(40)` hardcode. The only
hardcoded count is `Integer(26)` (L3540) = the HESSIAN_FREEZE **control reference** (its residual emitted as an
operand), which does NOT gate the live primary. Working tree == committed.

## ⚠ Provenance caveat — controls are verified OUT-OF-BAND, not in the committed `.out`
The committed `.out` was produced with `S11CB_PRIMARIES_ONLY=1` (deferred build-cost: the full-controls run
OOMs — the MATERIAL coupling-kernel BUILDS dominate, see [[feedback_measure_the_slow_phase_before_optimizing]]).
So `RUN_TASKS` = the 4 primaries only; `SKIPPED_TASKS` includes `HESSIAN_FREEZE`, `PROJECTION_EQUIVALENCE`,
and `FORM`. Confirmed independently: the emitted top-level tags contain no HESSIAN_FREEZE / PROJECTION_EQUIVALENCE
/ FORM result. The frozen-limit regression (40→26) and the lever-C value-identity are therefore proven ONLY by
the two legs' out-of-band ablations/numerics above (plus `~/.s11_build/S11c_b_89_leverC_numeric2.out` = 0 real
diffs). A reader of the `.out` alone has no in-band evidence of these; running the controls in-band remains a
deferred item pending the kernel-build optimization. This is a provenance note, not an engine defect (both
legs CLEAR).

## Disposition
Two independent legs CLEAR with agreeing from-scratch derivations. #89 PY §3a repair (+ lever-C) is cleared.
