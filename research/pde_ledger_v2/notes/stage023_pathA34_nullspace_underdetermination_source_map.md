# II-G5b (ledger_stage023) source map — pathA_34 native nullspace underdetermination departure (FAIL_UNDERDETERMINED_NOT_PREDICTIVE 2/2, the COMPLETING leg that DELIVERS the FAIL)

> Running-start prep captured 2026-07-10 (post stage022 = pathA_34 II-G5a cross-ℓ fingerprints DONE, commit `4e22959b`),
> before authoring the II-G5b reshape directive, so the directive can be written without re-discovery. **All line refs below
> VERIFIED against the current sources 2026-07-10** by a full orchestrator read of every 023-owned slice in BOTH engines
> (`pathA_34_cross_l_unification_sympy.py` = **1693 lines**; `pathA_34_cross_l_unification.wl` = **388 lines**;
> `reports/pathA_34_cross_l_unification.md` = 48 lines) PLUS an independent fresh-agent structural distillation that
> corroborated every `.py` mechanic — **no line drift** from the ranges the ▶ NEXT entry cites.
> Companion: `part2_gravity_atomic_split.md` (rows **022/023** L42–43 = the pathA_34 2-way split + the pathA_34 trip-ups
> **L90–92** + the cross-stage flows L104–109 + the ▶ NEXT stage023 entry L496–506) and the **stage022 pair + source map +
> reshape directive** (the SIBLING leg — 023 is the OTHER half of the SAME pathA_34 source; 022 is the EARNED-first slice,
> 023 the FAIL-delivering completing slice) + the **stage021 pair + directive** (the freshest full μ̂₀-free dimensional-gate
> exemplar — 023 reuses the exact `dim_of` (L,M,T)-triple machinery + the free-carrier control; 021's "which corruptions MAKE
> IT FIRE, not a single negative control" lesson) + **stage009/010** (the pathA_29 return residual `A_res ∝ ε_ℓ/(1+ε_ℓ)`,
> `Z=−M0·ε0/(1+ε0)` — 023's `A0/A1` are the cross-ℓ/port continuation of exactly these). Build-order id **023**, Part II.
> Source top-line (verbatim, report `reports/pathA_34_cross_l_unification.md:3`): **`FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** —
> **023 DELIVERS this FAIL** and completes the joint (022 landed the 1/2 EARNED-first PARTIAL). Proposed target stem:
> `ledger_stage023_nullspace_underdetermination` (confirm slug at directive authoring).

## ⭐ The FIVE headline points (READ FIRST)

1. **⚠ 023 is the COMPLETING, FAIL-DELIVERING leg of the pathA_34 2-way split (022/023).** 022 built the cross-ℓ
   `−(ℓ+1)/Λ_ℓ` fingerprints + the Gate-4 non-regression (the EARNED-first 1/2, script exit 0, `LOCAL_AUDIT_VERDICT =
   CROSS_L_FINGERPRINT_OK`). **023 builds the native nullspace underdetermination that makes the JOINT gate FAIL** — and
   UNLIKE 022 (whose script emits a local PASS token and only PRINTS the joint FAIL string), **023's script's own verdict IS
   the joint verdict `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** (the earned FAIL headline, computed from a genuine rank audit).
   023 exits 0 with the FAIL as its earned result (the pathA_36 / stage003 pattern: a script that PASSES by correctly
   COMPUTING a characterized FAIL departure — the FAIL is the earned physics, not a script error).
2. **⭐⭐ THE EARNED CONTENT = a GENUINELY-COMPUTED native nullspace, NOT a rigged verdict.** `build_rank_audit` (`.py`
   L615–699) computes, over **11 genuine generator dofs** `[OmegaU, OmegaW, Rmix, gU, gW, D0, K0c, K_eta, T_Omega, Z0_ret,
   Z1_ret]`, the constraint-Jacobian **rank 3** → **native nullspace dimension 8** → augmenting with the return-transfer
   gradients `∂T0, ∂T1` gives **return-augmented rank 5** → **return-moving nullity 2** (`return_aug_rank − rank0`) → the two
   return admittances `Z0_ret, Z1_ret` are genuinely UNPINNED by every collected Gate-5 named constraint. `native_moves =
   (return_moving_nullity > 0)` ⟹ `underdetermined` ⟹ `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`. **This is a real
   `sp.Matrix(rows).rank()` on symbolic partial-derivative rows** (`.wl` = native `MatrixRank`), with explicit
   nullspace-direction WITNESSES (`example_return_moving_directions` L639–657 builds the `Z0_ret`/`Z1_ret` unit vectors and
   confirms they preserve every constraint yet move `T0`/`T1`). **The v1 REJECTION locus is 023's and the current source is
   CLEAN of it** (see headline 5).
3. **⭐ THE ABLE-TO-FAIL CONTROL = the derived selector equation → `CROSS_L_RESIDUAL_PREDICTION`.** The gate is genuinely
   able-to-fail because a **real derived selector** (`Z0_ret = K0c`, `Z1_ret = K_eta + 2·T_Omega`, `selector_equations` `.py`
   L306–315) collapses the return-moving nullity to **0** → `native_moves = False` → the verdict FLIPS to
   `CROSS_L_RESIDUAL_PREDICTION` (the predictive/PASS-like token). So the DEFAULT verdict IS the predictive token
   (`base_verdict` L1008–1023 returns `CROSS_L_RESIDUAL_PREDICTION` unless a gate trips) — **`FAIL_UNDERDETERMINED` is EARNED**
   from the computed nullity, not baked. **Honest landing: the linear theory cannot pin the ℓ=0/1 return `{Z0_ret, Z1_ret}`;
   Gate-6's nonlinear closure (sim-deferred) must supply the selector — the first concrete, proven Gate-6 input.** (This is
   exactly what `docs/conceptual_foundation.md` §5 records: "Gate 5 (pathA_34) = FAIL_UNDERDETERMINED — linear theory can't pin
   the ℓ=0/1↔ℓ=2 return admittances, so Gate 6's nonlinear closure must supply that selector".)
4. **⭐ THE RESIDUALS `A0/A1` CONSUME 022's `{1, 1/2}` — a CHECKABLE, FORWARD-built (never back-solved) derive-vs-typed.**
   `build_residuals` (`.py` L464–508) reads 022's ℓ=0/1 radiative coefficients `v0=1, v1=1/2` and builds the residual
   amplitudes FORWARD: `A0 = i·v0·(aω/c_s)·M0·(1−T0)`, `A1 = i·v1·(aω/c_s)³·D1·(1−T1)` (the `(1−T_ℓ)` return-transmission
   multiplier = `ε_ℓ/(1+ε_ℓ)`), then CHECKS them against the **independently-formed pathA_29 residual form**
   `expected_A0 = i·aω·M0·ε0/(c_s(1+ε0))`, `expected_A1 = i·a³ω³·D1·ε1/(2c_s³(1+ε1))`. ⭐ **The `v1=1/2` is exactly what makes
   `A1 = expected_A1` (the `2·c_s³` in `expected_A1`'s denominator absorbs the `1/2`) — so a corrupted ℓ=1 coefficient fires
   the `A1_form` residual.** This is the checkable 022→023 consumption (the stage022-of-018 pattern). `ε_eff = Z_ret/K` is
   built FORWARD from the transfer definition (`z/k`), **NEVER back-solved from a residual** (the `FAIL_TAUTOLOGICAL` firewall).
5. **⭐⭐ THE pathA_34 v1 REJECTION locus IS 023's (the sharpest in Part II) — the current source is CLEAN; the reshape must
   PRESERVE that (part2 L90–92).** v1 was tri-review-REJECTED for **(a) a rigged-to-UNDERDETERMINED zero-padded constraint**
   (nullity forced by a hand-assembled zero-padded matrix), **(b) flag-driven probes** (verdict tokens stamped, not computed),
   and **(c) a headline-only `.wl`**. The CURRENT source fixes all three: (a) nullity is a genuine `sp.Matrix(rows).rank()` /
   native `MatrixRank` on real symbolic Jacobian rows (no zero-padding); (b) every probe re-runs `run_gate`/`gateVerdictFor`
   and reads a COMPUTED verdict (no stamped token); (c) the `.wl` computes `MatrixRank`/`Series`/`FullSimplify` natively and
   hard-`Exit[1]`s on a native `headlineOk`. ⚠ **Back-solving `ε_eff`/`Z` from residuals is FORBIDDEN** — `build_provenance`
   (`.py` L899–970) is the firewall: `epsilon_eff_nonzero_value` is CLASSIFIED `deferred_branch_data` while underdetermined
   (its magnitude "not computed because native nullspace leaves Z0_ret/Z1_ret free"), and the `class_matches_computed` check
   fires `FAIL_TAUTOLOGICAL` if a derived thing is emitted as asserted (or an ε magnitude is emitted as derived).

## §1 The 023 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED 2026-07-10)

The whole computation is `build_*` helpers assembled by `run_gate` (L1027–1079) + `main` (L1670–1689). **023 owns the
transfers/residuals/rank-audit/dimensional-checker/provenance/decoupling/verdict ladder + probes R1/3a/3b/3c/3d/3f/3g/3h;
022 (DONE `4e22959b`) owns the fingerprint core + the Gate-4 non-regression + probe 3e.** The 023-owned cuts:

- **`selector_equations`/`selector_subs`/`selector_provenance` L306–358** — the nullspace-collapsing control.
  `selector_equations(mutation)` L306–315 returns `[]` for `"none"`, else `[Eq(Z0ret, K0c), Eq(Z1ret, Keta+2*TOmega)]`
  (both `"derived_pde_admissibility"` and `"asserted_unproven"`). `selector_subs` L318–325 → the `{Z0ret: K0c, Z1ret:
  Keta+2*TOmega}` dict (asserts each lhs is a bare Symbol). `selector_provenance` L332–358 carries the tags the firewall
  reads: `"none"`→`absent/derived=False/tautological=False`; `"derived_pde_admissibility"`→`derived=True/tautological=False`
  (source: a counterfactual Gate-6 PDE-admissibility equation); `"asserted_unproven"`→`derived=False/**tautological=True**`
  (an asserted selector with no named PDE provenance — the `3g` `FAIL_TAUTOLOGICAL` hook).
- **`positive_bounded_transfer` L361–376 + `build_transfers(mutation)` L379–461** — `T_dc`, `ε_eff`, boundedness.
  `k0=K0c`, `k1=Keta+2*TOmega`; `z0/z1 = selector_subs.get(Z*ret, Z*ret)`; `t0=k0/(k0+z0)`, `t1=k1/(k1+z1)`,
  `eps0=z0/k0`, `eps1=z1/k1` (⭐ **ε_eff FORWARD from `z/k`, NOT back-solved**). Mutation hooks on `z0/z1`:
  `perfect_return`→`0/0`, `wrong_sign_return`→`−Z0ret/−Z1ret`, `no_consistent_return`→`−2K0c/−2(Keta+2TOmega)`,
  `inject_null`→`z+=eta_null·k`, `decouple_knobs`→`t*=gain`. `eps_relation` L408–421 verifies `T_ℓ=1/(1+ε_ℓ)`,
  `1−T_ℓ=ε_ℓ/(1+ε_ℓ)`. Derived: `overcancel = (eps0==0 ∧ eps1==0)` L424, `no_consistent = (eps0==−2 ∧ eps1==−2 ∧ ¬admissible)`
  L425–429, `admissible_branch_exists` from `positive_bounded_transfer` (ε>0, 0<T<1). ⚠ All mutation-driven — no stamped verdict.
- **`build_residuals(fingerprints, transfers, mutation)` L464–508** — the `A0/A1` residuals vs pathA_29 (headline 4).
  `v0 = fingerprints["outgoing"][0]["radiative_coeff_z"]` (=1, **022**), `v1 = …[1]…` (=1/2, **022**);
  `a0 = i·v0·(a·ω/c_s)·M0·(1−t0)` L473, `a1 = i·v1·(a·ω/c_s)³·D1·(1−t1)` L474; `expected_a0 = i·a·ω·M0·eps0/(c_s(1+eps0))`
  L475, `expected_a1 = i·a³ω³·D1·eps1/(2c_s³(1+eps1))` L476; `A0/A1_residual_to_pathA29_form = a0−expected_a0` L482–483;
  `match` L485–493 = `{A0_form, A1_form, A0_order(ω¹), A1_order(ω³), positive_bounded, nonzero_epsilon,
  admissible_branch_exists}`; `ok` L494–502. ⭐ **The consumption is checkable derive-vs-typed** (corrupt `v1` → `A1_form` fires).
- **`GENERATOR_DOFS` L511–523 + the rank machinery L526–612 + `build_rank_audit` L615–699** — THE CENTRAL ARTIFACT (headline 2).
  - `GENERATOR_DOFS = [OmegaU, OmegaW, Rmix, gU, gW, D0, K0c, Keta, TOmega, Z0ret, Z1ret]` (11).
  - `rank_row(expr, dofs) = [compact(sp.diff(expr, dof)) for dof in dofs]` L526–527 (a genuine Jacobian row);
    `matrix_rank(rows) = int(sp.Matrix(rows).rank())` L530–533 (**native SymPy rank — the anti-v1 core**).
  - `named_constraints` L536–587 documents the 3 named Gate-5 constraints (`ell0_collective_gate2_stiffness`→`K0c` [handoff
    §8.2 / Gate-2 collective, = stage013]; `ell1_section9_4_harmonic_stiffness`→`Keta+2*TOmega` [handoff §9.4, = stage017's
    `K_η`+`T_Ω`]; `ell2_section10_port_kernel`→`port["P0_raw"]` [handoff §10.2–10.3, = stage017's grouped-P2 port kernel])
    with `touches_generator_dofs` + `linearized_row`.
  - `pathA29_premise_citation` L590–597: `Z_is_premise:True, boundary_dof:none` (the provenance for why `Z` survives —
    pathA_29 supplies no boundary dof or PDE equation that fixes the return admittance). `selector_constraint_rows` L600–612.
  - **`build_rank_audit` L615–699 (the nullity=8 / return-moving-nullity=2 computation):** `base_constraints =
    [port["P0_raw"], K0c, Keta+2*TOmega]` L621; `constraint_exprs = base_constraints + selector_rows_exprs` (0 or 2 selector
    rows); `rows = [rank_row(e, dofs)]`; `rank0 = matrix_rank(rows)` (=3 baseline); `nullity = len(dofs)−rank0` (=8);
    `grad_t0/grad_t1 = rank_row(t0/t1, dofs)`; `return_aug_rank = matrix_rank(rows + [grad_t0, grad_t1])` (=5);
    `return_moving_nullity = return_aug_rank − rank0` (=2); `native_moves = return_moving_nullity > 0`. `untouched_return_dofs`
    L634–638 (per `zret`, ANY constraint has nonzero derivative? → the survivors). `example_return_moving_directions`
    L639–657 (the explicit unit-vector witnesses). Returns `native_nullspace_dimension`, `return_moving_nullity`,
    `native_null_moves_return`, `underdetermined_not_predictive = native_moves`, `branch_selector = selector_provenance`,
    `injected_null_probe`, `why_Z_survives`. ⭐ **With the selector present**: the 2 extra rows `Z0ret−K0c`,
    `Z1ret−(Keta+2TOmega)` raise `rank0` to 5 AND make `t0=t1=1/2` (const, zero gradient) → `return_moving_nullity=0` →
    `determined` → `CROSS_L_RESIDUAL_PREDICTION`. Genuine flip.
- **The dimensional checker L702–896** — the residual dim gate (headline: 021's `dim_of` machinery, reused). `DimError`
  L702–703; `dim_of` L718–744 (recursive `Mul/Pow/Add` LMT-triple; raises on missing symbol / non-numeric exponent / **sum
  mismatch** L741–742 — the able-to-fail core, IDENTICAL discipline to stage021). `SOURCED_DIMS` L775–801
  (`M0=(0,1,−1)`, `D1=(1,1,−1)`, `Z0ret=Z1ret=K0c=Keta=TOmega=(0,1,−2)`, `q_free=ZERO_DIM`, gains/eta_null `ZERO_DIM`).
  `EXPECTED_DIMS` L803–817 (`A0=(0,1,−1)`, `A1=(1,1,−1)`, `T0=T1=ε0=ε1=P0_physical=ZERO_DIM`). `run_dimension_check`
  L833–896: checks `[A0],[A1]` (from `residuals["leading"]`), `[T0],[T1],[ε0],[ε1],[N0],[D0]`,
  `[P0_physical]=(c_s/a)²·P0_raw` + `K0/K1`; two controls: `corrupt_sourced` (`dims[M0]+=(1,0,0)` → `[A0]` mismatch →
  `FAIL_DIMENSIONAL`) L842–843, `corrupt_free_carrier` (`dims[q_free]=(7,0,0)` — but `q_free` is in NO checked expression →
  still passes → `NO_FAIL`, the **free-carrier-independence control**, L844–845 + `free_carrier_independence…` L893–895);
  `verdict = "NO_FAIL"/"FAIL_DIMENSIONAL"` L890.
- **`build_provenance(mutation, rank)` L899–970** — the `FAIL_TAUTOLOGICAL` firewall (headline 5). Each item carries
  `computed_class` vs `emitted_class`; `class_matches_computed = (emitted == computed)` L943; `ok = all(...)` L951. ⭐ The
  anti-back-solve: `epsilon_eff_nonzero_value.computed_class = "deferred_branch_data" if underdetermined else "derived_in_gate"`
  L917–925, with `magnitude_note = "not computed because native nullspace leaves Z0_ret/Z1_ret free"`. Two ways to trip:
  `assert_not_derive` L939–940 (forces `Y_l_out_fingerprints.emitted_class="asserted_literal"` vs computed `"derived_in_gate"`
  → mismatch — the `3g` hook); selector `tautological_assertion` L944–950 (injects a `branch_selector` item with
  emitted `"asserted_literal"` vs computed `"derived_in_gate"`).
- **`detect_decoupling(port, transfers, mutation)` L973–1005** — the `FAIL_DECOUPLED` probe (the anti-rig that catches
  introducing a free knob to dial ε after P0 is fixed). Introduced set `[gain0, gain1]`; `independently_moves_return` L989–994
  (gain0 moves T0 not T1, gain1 moves T1 not T0); `p0_unaffected` (neither gain in `port["P0_raw"].free_symbols`); `decoupled`
  L1002–1004. All from `sp.diff` — genuine.
- **`base_verdict(conditions)` L1008–1023** — the ordered FAIL ladder, DEFAULT `CROSS_L_RESIDUAL_PREDICTION`:
  `decoupled→FAIL_DECOUPLED, tautological→FAIL_TAUTOLOGICAL, quad_regression→FAIL_QUAD_REGRESSION, dimensional→FAIL_DIMENSIONAL,
  no_consistent_return→FAIL_NO_CONSISTENT_RETURN, overcancel→FAIL_OVERCANCEL, epsilon_mismatch→FAIL_EPSILON_MISMATCH,
  underdetermined→FAIL_UNDERDETERMINED_NOT_PREDICTIVE, able_to_fail_bad→FAIL_TAUTOLOGICAL` → else `CROSS_L_RESIDUAL_PREDICTION`.
  ⚠ **`quad_regression` (022's Gate-4 non-regression) sits in the joint ladder** — see §3 for the 023 cut decision.
- **`run_gate(mutation)` L1027–1079** — the assembly; `conditions` L1045–1059 wires every gate from computed objects
  (`tautological = ¬provenance["ok"] ∨ selector tautological`; `quad_regression = ¬gate4["ok"]` [consumes 022];
  `underdetermined = rank["underdetermined_not_predictive"]`; `epsilon_mismatch` from `residuals["pathA_29_comparison"]`;
  `able_to_fail_bad = False` here — the real able-to-fail aggregate is applied in `build_final_payload` L1501–1502).
- **`ablation(baseline, mutated, expected_fail)` L1082–1099** — the DYNAMIC two-verdict re-run (`fail_suppressed =
  with_mutation==expected_fail ∧ without_mutation≠expected_fail` L1096; carries both `conditions` dicts). Re-runs `run_gate`
  on both — a genuine two-verdict re-run (the pathA_34 v1 constant-`self_ablation` trip-up avoided).
- **`build_counterfactuals(actual_baseline)` L1102–1244** — the probe suite. `clean` baseline = selector present L1103–1107.
  023-owned probes:

  | Probe | `.py` lines | Mutates | expected_fail |
  |---|---|---|---|
  | **R1** port-kernel dependency | L1111–1123 | `corrupt_port_kernel` (`OmegaU→2·OmegaU`) | asserts `P0_changes` AND `ell2_determinacy_row_changes` (not a FAIL token) |
  | **3a** decouple knobs | L1124–1131 | `decouple_knobs` on `clean` | `FAIL_DECOUPLED` |
  | **3b** null direction / selector flip | L1132–1162 | selector control + `inject_null` | flip: without-selector→`FAIL_UNDERDETERMINED`, with-selector→`CROSS_L_RESIDUAL_PREDICTION`; + injected-null detector |
  | **3c** wrong sign antilocalizing | L1163–1170 | `wrong_sign_return` | `FAIL_EPSILON_MISMATCH` |
  | **3d** perfect return | L1171–1178 | `perfect_return` | `FAIL_OVERCANCEL` |
  | **3f** corrupt dimension | L1187–1203 | `corrupt_dimension` + `q_free` free-carrier | `FAIL_DIMENSIONAL` (+ free-carrier `NO_FAIL`) |
  | **3g** assert-not-derive | L1204–1215 | `selector_equation_set="asserted_unproven"` | `FAIL_TAUTOLOGICAL` |
  | **3h** no consistent return | L1216–1223 | `no_consistent_return` | `FAIL_NO_CONSISTENT_RETURN` |

  `flags` L1225–1239 gate on `self_ablation["fail_suppressed"]` (+ 3b's `verdict_flips_to_prediction`, 3f's free-carrier
  `NO_FAIL`, R1's both-changes); `able_to_fail_ok = all(flags.values())` L1243. **⚠ 3e (`break_gate4` → `FAIL_QUAD_REGRESSION`,
  L1179–1186) is 022's** — NOT rebuilt as an 023 tooth.

- **⭐ Shared helpers 023 uses (NOT cut boundaries):** `compact` L84, `bool_zero` L108–113, `hstr`/`equation_text` (rendering).
  Symbol decls L44–63 — 023 uses `M0, D1, R0, R1, K0c, Keta, TOmega, Z0ret, Z1ret, q_free, gain0, gain1, eta_null, OmegaU,
  OmegaW, Rmix, gU, gW, D0` (+ `a, c_s, omega` shared from 022). ⚠ **`Delta`/`Sport` L57–58 are declared-but-NEVER-referenced
  (vestigial) — DROP them in the reshape.**

- **⭐ CLEAN CUT — 023 owns L306–358 (selector) + L361–508 (transfers+residuals) + L511–699 (rank audit) + L702–896 (dim
  checker) + L899–1005 (provenance firewall + decoupling) + L1008–1079 (verdict ladder + assembly) + L1082–1099 (ablation) +
  the probes R1/3a/3b/3c/3d/3f/3g/3h L1111–1223. It touches NONE of 022's fingerprint core (L127–222) NOR 019's prefactor:**
  - **022 (DONE `4e22959b`) — CONSUME as PROVENANCE/typed-input, do NOT rebuild:** the fingerprint core `spherical_j`/`spherical_y`/
    `dtn_branch`/`build_fingerprints` L127–222 (023 needs ONLY the ℓ=0/1 radiative coefficients `{1, 1/2}` → cite typed, §1c);
    the Gate-4 non-regression + probe `3e` (§3 — the `quad_regression` condition is CONSUMED-False from 022, not re-derived).
  - **019 (DONE `f1c426f9`) — STRIP entirely from 023:** the squared-denominator prefactor `P(ω)=D₀N/D^cons²` (`.py` inside
    the old `build_gate4_non_regression`; `.wl` L88–113) — 023 has NO prefactor content. (The port kernel `build_port_kernel_for`
    L281–303 IS used by 023 — it supplies `P0_raw`, the ℓ=2 named-constraint row + the R1 probe target — but the PREFACTOR
    algebra on top of it is 019's, stripped.)

## §1b The `.wl` 023 slice (VERIFIED) — ⚠ THE KEY RESHAPE DECISION: genuine native engine, but is it a TRANSLITERATION?

⭐ **The pathA_34 `.wl` 023 region (L115–298) is a GENUINE native engine** (native `MatrixRank`, `D`, `Series`,
`FullSimplify`, recursive `dimOf`; hard-`Exit[1]` on a native `headlineOk`; does NOT `Import`/`Get` the `.py`). The fresh-agent
distillation recommended **keep-native**. **BUT** — ⚠ apply the stage022 lesson: **the 023 `.wl` region is STRUCTURALLY PARALLEL
to the `.py`** (same `rankDofs` list + order L128, same `base_constraints = {P0port, K0c, Keta+2 TOmega}` L132, same
`augRank − rank0` return-moving-nullity subtraction L136, same `gateVerdictFor` `Which`-ladder L168–177 mirroring
`base_verdict`, same `dimOf` recursion mirroring `dim_of`). By the `MATHEMATICA_MIRROR_POLICY.md` transliteration screen —
which is NECESSARY-NOT-SUFFICIENT to lack `Import`/`Get` — **this is the SAME concern that forced stage022's `branchData`
re-authoring.** ⭐ **The directive must RESOLVE this (route it through the Codex→Grok→Codex bookend), not rubber-stamp
"keep-native":**
- **Recommendation (leaning RE-AUTHOR the nullspace core, per reshape-spec §5 point 3 + the stage022 precedent):** a
  materially-different Wolfram route for the underdetermination DOES exist — instead of the `.py`'s `augRank − rank0`
  rank-subtraction, exhibit the nullspace CONSTRUCTIVELY: `NullSpace[constraintJacobian]` (or `RowReduce`/`Reduce[ForAll]`)
  to get the null directions directly, then show the `{Z0_ret, Z1_ret}` return-admittance subspace is 2-dimensional within it
  and moves `T0/T1` (a constructive witness, not a rank difference). That is genuinely independent corroboration (a different
  algorithm answering the same question). The verdict ladder can stay native `Which`, but the nullity/return-moving-nullity
  numbers should come from a different decomposition than the `.py`'s.
- **Alternative (if Codex/Grok judge native `MatrixRank` materially independent enough):** the rank computation is canonical
  (there is no "hand-built vs built-in" fork like 022's spherical Bessel), so keep-native `MatrixRank` MAY be acceptable — but
  ONLY with the transliteration screen EXPLICITLY applied (dof-name/decomposition parallelism argued non-load-bearing) and the
  `gateVerdictFor`/`dimOf` re-structured enough to not read as a line-by-line mirror. ⚠ Default to RE-AUTHOR unless the review
  affirmatively clears keep-native.
- The 023 `.wl` slice, whichever route: transfers/residuals L115–126 (consume 022's `out0/out1["radiativeCoeff"]` = `{1,1/2}`,
  or cite typed §1c), rank audit L128–149, verdict machinery L151–177, dim checker L179–265, probes + `headlineOk` L267–298.
- **⚠ `.wl` STRIP (022/019 DONE):** the fingerprint block `j0..y2`/`branchData`/`out*`/`in*`/`fingerprintOk` L35–86 (022's — 023
  cites `{1,1/2}` typed OR consumes them; it does NOT re-run the full fingerprint battery) + the prefactor `Pport`/`N0port`/
  `Nomega`/`prefObj`/`p*`/`resP*`/`gate4Ok` L88–113 (019/022's — but `Pport`/`N0port`/`P0port` L88–91 ARE needed as `P0_raw`
  for the rank-audit ℓ=2 constraint row + R1 probe; keep the P0 construction, drop the prefactor series L93–113).
- **⚠ `.wl` SEVER (the bridge):** `scratchDir`/`yamlOut` setup L15–22 + the whole `lines={…}` YAML assembly L300–383 +
  `Export[yamlOut,…]` L385. **Zero file I/O**, print-only + `fail[msg_]` (L5, already present) on failure, `Exit[0]` on the
  native `headlineOk`. **Dual-engine agreement = transcript-level** (both engines print the SAME `nullity=8`,
  `return_moving_nullity=2`, `baseline_verdict=FAIL_UNDERDETERMINED_NOT_PREDICTIVE`, `selector_verdict=CROSS_L_RESIDUAL_PREDICTION`,
  the residual `A0/A1` residual-zero, the dim `dimensional_ok`/corrupt-sourced-FAIL/corrupt-free-NO_FAIL, and the 3c/3d/3f/3h
  probe tokens). Arity discipline (def/call scan + unevaluated-leakage transcript scan; the `.wl` has `Module`s in
  `rankAuditFor`/`gateVerdictFor`/`dimOf`).

## §1c The consumption resolution (READ — the checkable derive-vs-typed vs the provenance-only cites)

⭐ **UNLIKE 022 (whose ONLY checkable consumption was the stage018 non-regression), 023 has TWO checkable cross-stage
relations (the 022 fingerprint coefficients feeding `A0/A1`, and the pathA_29 residual form) PLUS several provenance cites.**

- **⭐ 022's ℓ=0/1 radiative coefficients `{1, 1/2}` — a CHECKABLE derive-vs-typed feeding `A0/A1`.** 023 reads `v0=1, v1=1/2`
  (stage022's earned exports) and builds `A0/A1` forward, checking them against the pathA_29 form. **The genuine tooth: corrupt
  the cited `v1` (say `1/2→1/3`) → `A1 ≠ expected_A1` → the `A1_form` residual fires** (the `2·c_s³` in `expected_A1` is exactly
  the `1/2`). ⚠ **Decision for the directive:** does 023 (i) CITE `{1, 1/2}` TYPED (the stage022-of-018 pattern — 022 OWNS the
  fingerprints, so 023 cites its earned exports and does NOT re-run the fingerprint battery), or (ii) re-derive them
  self-contained from its own Hankel? **Recommend (i) CITE-typed** (clean split: 022 owns the fingerprints; the checkable
  consumption is the `A1_form` residual, not a re-derivation). ⚠ Ensure `expected_A0/A1` are built from the INDEPENDENT pathA_29
  residual form (stage009's `A_res ∝ ε_ℓ/(1+ε_ℓ)` + the ℓ-order `(aω/c_s)^(2ℓ+1)` + the dipole `1/2`), NOT typed to match
  `A0/A1` — else it degrades to X≡X. The two sides must have independent provenance (022's `{1,1/2}`; pathA_29's residual form).
- **⭐ The pathA_29 residual form (stage009/010) — the SECOND checkable relation.** `expected_A0 = i·aω·M0·ε0/(c_s(1+ε0))`,
  `expected_A1 = i·a³ω³·D1·ε1/(2c_s³(1+ε1))` are the cross-ℓ/port continuation of stage009's `RETURN_RESIDUAL_PREDICTION`
  (`A_res ∝ ε_ℓ/(1+ε_ℓ)`, `Z = −M0·ε0/(1+ε0)`) and stage010's return channels. Cite `stage009`/`stage010` — the `A0/A1`
  residual-form match IS the checkable relation (the return residual survives the cross-ℓ/port language). ⚠ The `ε_ℓ` here
  is `Z_ret/K` (FORWARD from the transfer), the SAME `ε_ℓ` family stage009 parameterized — but now UNPINNED (the FAIL).
- **008's `R0=−M0`/`R1=−D1` cancellation targets + `M0`/`D1` moments (PROVENANCE + the residual amplitudes).** `M0`/`D1` enter
  `A0/A1` as the source amplitudes; `R0`/`R1` appear only in the dim checker (`[R0]=[M0]`, `[R1]=[D1]`) and as typed provenance
  tags in `build_provenance` (`external_R0_equals_minus_M0`, `external_R1_equals_minus_D1`, class `external_bridge_input`).
  Cite `stage008` (the ℓ=0/1 constraint targets + the `M0`/`D1` moments). No dual-site on `R0`/`R1` (they are external-bridge tags).
- **pathA_29's `Z_is_premise=True, boundary_dof=none` (the KEYSTONE PROVENANCE).** `pathA29_premise_citation` L590–597 is WHY
  `Z0_ret`/`Z1_ret` survive: pathA_29 supplies no boundary dof / PDE equation fixing the return admittance. Cite `pathA_29`
  (via the pathA_29 fold = stages 009/010) as the premise that makes the underdetermination EARNED (not an omission). ⚠ This is
  the load-bearing "the linear theory genuinely cannot pin it" justification — the FAIL is honest BECAUSE `Z` is a premise.
- **017's grouped-P2 port kernel (`P0_raw`) — a CONSUMED structural input** (the ℓ=2 named-constraint row + the R1 probe target).
  `port["P0_raw"] = (OmegaU²·gW + Rmix·gU)²/(OmegaU²·OmegaW² − Rmix²)²/D0`. Cite `stage017` (the exported ℓ=2 port kernel; the
  R1 probe genuinely perturbs it → the ℓ=2 determinacy row changes). The port scalars `{B̃,Z̃}` were tracked-downstream-pinned
  at 017 (register row 180) — 023 is one of the downstream pins.
- **013's ℓ=0 collective stiffness (`K0c`) + 017's ℓ=1 harmonic stiffness (`K_eta+2·T_Omega`) — CONSUMED as the ℓ=0/1 named
  constraints.** `K0c` ← stage013's Gate-2 collective (a,L) reduction (§8.2); `K_eta` = R29-derived `T_wβ²` (013), `T_Omega` =
  017's counted `T_Ω`, so `K1 = K_eta + 2·T_Omega` (the ℓ=1 angular stiffness, `λ_m=ℓ(ℓ+1)=2` for ℓ=1 — 016's covariance).
  Cite `stage013`/`stage017`. ⭐ **These are the register crux (§6): likely DERIVED manifestations, NOT new counted knobs.**

## §2 The 023 claim-set (derive + assert; report/directive quotes)

- **(a) ⭐ The native nullspace underdetermination (EARNED — the FAIL-delivering headline; report `:5`, `:12`, `:47`).** Over 11
  genuine generator dofs, the 3 collected Gate-5 named constraints `{P0_raw, K0c, K_eta+2·T_Omega}` have constraint-Jacobian
  **rank 3** → **native nullspace dimension 8**; augmenting with `∂T0, ∂T1` gives **return-augmented rank 5** → **return-moving
  nullity 2** → the return admittances `{Z0_ret, Z1_ret}` are UNPINNED (`untouched_return_dofs = ['Z0_ret', 'Z1_ret']`), with
  explicit unit-vector witnesses that preserve every constraint yet move `T0/T1`. `native_null_moves_return = True` ⟹
  `underdetermined_not_predictive = True` ⟹ **`FAIL_UNDERDETERMINED_NOT_PREDICTIVE`**. ⭐ COMPUTED via `sp.Matrix(rows).rank()`
  / native `MatrixRank` (report `:12`: "Native nullspace dimension: 8; return-moving nullity: 2; moves return: True; selector
  present: False"). The pathA_29 premise (`Z_is_premise=True, boundary_dof=none`) is the provenance for why `Z` survives.
- **(b) ⭐ The selector-equation control → `CROSS_L_RESIDUAL_PREDICTION` (the ABLE-TO-FAIL; report `:15`, `:37`).** The derived
  selector `{Z0_ret = K0c, Z1_ret = K_eta + 2·T_Omega}` collapses the return-moving nullity to 0 → `determined` →
  `CROSS_L_RESIDUAL_PREDICTION`. So the gate genuinely DISTINGUISHES: default (no selector) → FAIL; selector supplied → the
  predictive token. **Honest landing: Gate-6's nonlinear closure (sim-deferred) must supply the `{Z0_ret, Z1_ret}` selector —
  the first concrete, proven Gate-6 input.**
- **(c) ⭐ The `A0/A1` scalar/dipole residuals vs pathA_29 (EARNED — form/sign/order; report `:19–21`, `:47`).** `A0 =
  i·M0·Z0_ret·aω/(c_s(K0c+Z0_ret))`, `A1 = i·D1·Z1_ret·a³ω³/(2c_s³(K_eta+2T_Omega+Z1_ret))` (report `:19–20`), `ε0_eff =
  Z0_ret/K0c`, `ε1_eff = Z1_ret/(K_eta+2T_Omega)` (report `:21`) — the cross-ℓ/port continuation of stage009's
  `RETURN_RESIDUAL_PREDICTION`, consuming 022's `{1, 1/2}`, form/order matched to the pathA_29 residual (`A0_order=ω¹`,
  `A1_order=ω³`), conditional on a positive bounded branch. ⭐ FORWARD-built, NEVER back-solved.
- **(d) ⭐ The μ̂₀-free-style dimensional gate (EARNED — report `:24–27`).** `[A0]=(0,1,−1)`, `[A1]=(1,1,−1)`,
  `[T]=[ε]=[P0_physical]=0` from the sourced dims via the recursive `dim_of` (021's machinery); `dimensional_ok=True`;
  the sourced corruption (`[M0]` wrong) → `FAIL_DIMENSIONAL`; the unconstrained free-carrier (`q_free`) corruption → `NO_FAIL`
  (the free-carrier-independence control, report `:27`).
- **(e) The provenance firewall (EARNED — the `FAIL_TAUTOLOGICAL` anti-back-solve).** `class_matches_computed` for every
  provenance item; `ε_eff` magnitude CLASSED `deferred_branch_data` while underdetermined (not emitted as derived); the
  `assert_not_derive` + tautological-selector mutations trip `FAIL_TAUTOLOGICAL`.
- **(f) The full probe suite (EARNED able-to-fail).** R1 (port-kernel dependency: `P0_changes` ∧ `ell2_row_changes`), 3a
  (`FAIL_DECOUPLED`), 3b (inject_null detector + the selector flip), 3c (`FAIL_EPSILON_MISMATCH`), 3d (`FAIL_OVERCANCEL`), 3f
  (`FAIL_DIMENSIONAL` + free-carrier `NO_FAIL`), 3g (`FAIL_TAUTOLOGICAL`), 3h (`FAIL_NO_CONSISTENT_RETURN`) — each a DYNAMIC
  two-verdict `ablation` re-run (report `:33–43`).
- **(g) The 023 landing (the COMPLETING slice — DELIVERS the joint FAIL).** 023's script verdict IS
  `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (the earned FAIL headline, computed). The report's `Deferred:` list `:47` is 023's:
  "the scalar/dipole return magnitude and nonzero prediction, because the native nullspace leaves `ε_eff` free at the linear
  Gate-5 level." ⭐ The SCRIPT exits 0 (correctly computing the earned FAIL + every tooth firing — the stage003/pathA_36 pattern).

## §3 Reshape cost (the bridge to sever) + the shared-machinery + the 023 verdict (⭐ the KEY 022/023 seam decision)

Same family as pathA_30–34 / 018–022 (the cross-script scratch-YAML reshape, driven by the presence of the MMA scratch file,
`.py` L1670–1689). **Reshape = sever ALL file I/O both directions:** `.py` — drop `main`'s `yaml_write(SYM_YAML)`/
`yaml_read(MMA_YAML)`/`yaml_write(RESULTS_YAML)`/`REPORT_MD.write_text` (L1670–1689) + `yaml_write`/`yaml_read` (L93–107) +
`engine_summary`/`build_engine_probe_summary`/`sympify_engine_number`/`compare_engines`/`build_final_payload`/`build_report`
(L1259–1667) + the path constants (L36–39); `.wl` — drop the `scratchDir`/`yamlOut` setup L15–22 + the `lines={…}` YAML
assembly L300–383 + `Export` L385. Each engine → standalone: print-only, `expect_zero`/`expect_bool`/`expect_fail`-style
asserts (`.py` local ledger idioms — `banner`/`subbanner`/`_record_pass`/`_record_fail`, a `Verdict labels:` block, tallies,
`OVERALL PASS`/nonzero exit — the stage022 template verbatim), `fail[]`/`Exit[1]` on failure (`.wl`). **Zero file I/O.** Arity
discipline. **⚠ The `.wl` route (keep-native vs RE-AUTHOR the nullspace) is §1b — the KEY genuineness decision.**

**⭐⭐ THE 022/023 VERDICT SEAM (the KEY reshape decision for 023).** UNLIKE 022 (which built a 022-LOCAL verdict and only
PRINTED the joint FAIL string), **023's OWN verdict IS the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** — 023 builds the joint
`base_verdict` ladder (it is the FAIL-delivering leg). ⚠ **But the joint ladder's `quad_regression` condition is 022's**
(`conditions["quad_regression"] = ¬gate4["ok"]`, from 022's Gate-4 non-regression). Resolve cleanly:
- **⭐ 023 CONSUMES 022's Gate-4 result as a CITED-True input (`quad_regression = False`, provenance from stage022), it does
  NOT re-run the Gate-4 non-regression** (that is 022's owned EARNED content, DONE). So 023's `base_verdict` ladder either (i)
  DROPS the `quad_regression` rung (cite 022 established it False) so 023's ladder is `{decoupled, tautological, dimensional,
  no_consistent_return, overcancel, epsilon_mismatch, underdetermined, able_to_fail_bad}` → default `CROSS_L_RESIDUAL_PREDICTION`;
  OR (ii) keeps `quad_regression` as a consumed-False cited condition (022 provenance) NOT recomputed here. **Recommend (i)**
  (cleanest split — 022 owns Gate-4, 023 owns the return/nullspace ladder; the joint token is unaffected since Gate-4 passed).
  ⚠ Whichever, 023 must NOT rebuild the Gate-4 fingerprint non-regression or probe `3e` (022's DONE tooth).
- **023's verdict is the EARNED FAIL:** the baseline would return `CROSS_L_RESIDUAL_PREDICTION` but for `underdetermined=True`
  from the computed rank audit. The selector control (a separate `run_gate`/`gateVerdictFor` call with the selector) returns
  `CROSS_L_RESIDUAL_PREDICTION` — the able-to-fail proof. Print BOTH (`baseline_verdict` = FAIL, `selector_control_verdict` =
  prediction), the stage022-of-021 discipline but now the baseline IS the earned FAIL.
- **⚠ 023 does NOT print a 022-style "1/2 PARTIAL" line** — 023 COMPLETES the joint. The landing: `VERDICT:
  FAIL_UNDERDETERMINED_NOT_PREDICTIVE (2/2, COMPLETING — the native nullspace underdetermination departure that DELIVERS the
  FAIL; 022 landed the EARNED-first cross-ℓ fingerprints + Gate-4 non-regression 1/2)`.

**Acceptance (dual-engine, both exit 0, CWD-independent):**
- Run each engine from the **repo root** AND from a **foreign CWD** (e.g. `/tmp`), both print-only, both exit 0, no files
  written (verify with `find` for new files; a `.wl` `Export` slip is the classic leak).
- Both engines emit the same transcript: `native_nullspace_dimension=8`, `return_moving_nullity=2`,
  `native_null_moves_return=True`, `baseline_verdict=FAIL_UNDERDETERMINED_NOT_PREDICTIVE`,
  `selector_control_verdict=CROSS_L_RESIDUAL_PREDICTION`; the residual `A0/A1` residual-to-pathA29-form = 0 (with `A0_order=ω¹`,
  `A1_order=ω³`); `dimensional_ok=True`, corrupt-sourced → `FAIL_DIMENSIONAL`, corrupt-free-carrier → `NO_FAIL`; the probe
  tokens 3a `FAIL_DECOUPLED` / 3c `FAIL_EPSILON_MISMATCH` / 3d `FAIL_OVERCANCEL` / 3f `FAIL_DIMENSIONAL` / 3g
  `FAIL_TAUTOLOGICAL` / 3h `FAIL_NO_CONSISTENT_RETURN` / 3b selector-flip; R1 `P0_changes ∧ ell2_row_changes`.
- **All able-to-fail teeth fire at their own assert** (per-tooth ablation, §5) — the rank-audit nullity, the return-moving
  nullity, the selector flip, the `A0/A1` form/order (incl. the corrupted-`v1` consumption tooth), the dim gate (sourced-FAIL,
  free-carrier-NO_FAIL), the provenance firewall (`assert_not_derive`, tautological-selector), the decoupling probe, each
  DYNAMIC `ablation` two-verdict re-run; the `.wl` independent-route + arity.

## §4 Consumed / exported

- **Consumes:**
  - **⭐ stage022's ℓ=0/1 radiative coefficients `{1, 1/2}` — the CHECKABLE derive-vs-typed feeding `A0/A1`** (§1c): corrupt
    `v1` → `A1_form` fires. Cite `stage022`. ⚠ Cite-typed (022 owns the fingerprints), NOT re-derived.
  - **⭐ the pathA_29 residual form (stage009/010) — the SECOND checkable relation** (`A_res ∝ ε_ℓ/(1+ε_ℓ)`, cross-ℓ/port
    continuation). Cite `stage009`/`stage010`. The `A0/A1` residual-form match IS the check.
  - **pathA_29's `Z_is_premise=True, boundary_dof=none`** — the keystone premise for why the underdetermination is EARNED.
    Cite the pathA_29 fold.
  - **008's `R0=−M0`/`R1=−D1` + `M0`/`D1` moments** — the source amplitudes (`M0/D1`) + external-bridge tags (`R0/R1`). Cite
    `stage008` PROVENANCE.
  - **017's grouped-P2 port kernel `P0_raw`** — the ℓ=2 named-constraint row + the R1 probe target. Cite `stage017`.
  - **013's ℓ=0 collective stiffness `K0c` + 017's ℓ=1 harmonic stiffness `K_eta+2·T_Omega`** — the ℓ=0/1 named constraints.
    Cite `stage013`/`stage017`. (Register crux §6.)
  - **`c_s`/`a`** — R1 units carrier (stage005) + `CONV` pin (`z=aω/c_s` realizes `z^(2ℓ+1)→ω^(2ℓ+1)` in `A0/A1`). Distinct
    from the frozen-wall `c_S` (011–017).
- **Exports (→ Cluster C 024–027 + Part VII):** the native nullspace departure (dim-8/return-nullity-2) + the **Gate-6
  `{Z0_ret, Z1_ret}` selector need** (the first concrete, proven Gate-6 input — → the Part-VII open-items register + the
  Gate-6 caveat on 024/027/028) + the `A0/A1` residual class + the completed pathA_34 joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`
  (022∧023). Per the cross-stage flows (`part2` L104–108). ⭐ **This COMPLETES the pathA_34 fold + closes Cluster A/the Gate-1–5
  gravity ladder** — next is Cluster C (024–029), then the MIDWAY KNOB AUDIT.

## §5 Teeth candidates (023-specific, per-tooth ablation MANDATORY — mutate the named object, confirm exit-1 AT its own assert)

1. **⭐⭐ The native nullspace teeth (`native_nullspace_dimension=8`, `return_moving_nullity=2`, `native_null_moves_return`).**
   COMPUTED via `sp.Matrix(rows).rank()` / native `MatrixRank` on the real constraint Jacobian. Per-tooth: mutate a constraint
   (e.g. inject a constraint that touches `Z0ret` → the return-moving nullity drops → `native_moves` flips → the underdetermined
   assert changes) or the `inject_null` probe (adds `eta_null` → a spurious null direction the detector catches). ⚠ **The
   firewall: the nullity MUST be a genuine rank of symbolic Jacobian rows — forbidden to hardcode `8`/`2` or a zero-padded
   constraint matrix (the v1 REJECTION locus).** Mutating `Z0ret`'s constraint coverage must change the computed nullity.
2. **⭐ The selector-flip tooth (`selector → CROSS_L_RESIDUAL_PREDICTION`; the able-to-fail control).** The derived selector
   `{Z0_ret=K0c, Z1_ret=K_eta+2·T_Omega}` collapses the nullity → 0 → `CROSS_L_RESIDUAL_PREDICTION`. Per-tooth: with the
   selector, `return_moving_nullity=0` and the verdict is the prediction; WITHOUT it, `=2` and the verdict is the FAIL — the
   two-verdict `ablation` (probe 3b). Neuter the selector (make it not touch `Z*ret`) → no flip → the control flags not-able-to-fail.
3. **⭐ The `A0/A1` residual form/order + the `{1,1/2}` consumption tooth.** `A0/A1` FORWARD-built from 022's `{1,1/2}`,
   checked vs the pathA_29 form. Per-tooth: (a) corrupt the cited `v1` (`1/2→1/3`) → `A1_form` residual ≠ 0 → fires (the
   checkable consumption); (b) corrupt the order-realization (`(aω/c_s)³→²`) → `A1_order` fires. ⚠ Derive-vs-typed vs the
   INDEPENDENT pathA_29 form, NEVER a back-solve of `ε_eff` (the `FAIL_TAUTOLOGICAL` firewall).
4. **⭐ The dimensional gate (021's machinery): sourced-corruption FIRES, free-carrier does NOT.** `[A0]/[A1]` from the sourced
   dims. Per-tooth: corrupt `[M0]` (`+(1,0,0)`) → `[A0]` mismatch → `FAIL_DIMENSIONAL` (probe 3f); corrupt the free carrier
   `[q_free]` (`(7,0,0)`) → `NO_FAIL` (the free-carrier-independence control — `q_free` is in no checked expression). ⚠ The
   021 lesson: which corruptions MAKE IT FIRE (the sourced ones) proves able-to-fail, not a single negative control.
5. **⭐ The provenance firewall (`FAIL_TAUTOLOGICAL`) — the anti-back-solve.** `class_matches_computed`. Per-tooth: (a)
   `assert_not_derive` → `Y_l_out.emitted_class="asserted_literal"` ≠ computed `"derived_in_gate"` → `FAIL_TAUTOLOGICAL` (probe
   3g); (b) the tautological selector (`asserted_unproven`) → the `branch_selector` item mismatches → `FAIL_TAUTOLOGICAL`. ⚠
   Confirm `epsilon_eff_nonzero_value.computed_class` is `deferred_branch_data` while underdetermined (the ε magnitude is NOT
   emitted derived).
6. **⭐ The transfer/return-family probes (3c/3d/3h) + the decoupling probe (3a).** `wrong_sign_return` → `FAIL_EPSILON_MISMATCH`
   (3c), `perfect_return` (`Z=0`) → `FAIL_OVERCANCEL` (3d), `no_consistent_return` (`Z=−2K`) → `FAIL_NO_CONSISTENT_RETURN`
   (3h), `decouple_knobs` (gain0/gain1 dial ε after P0 fixed) → `FAIL_DECOUPLED` (3a). Per-tooth: each mutation → its own token
   via a DYNAMIC `ablation` two-verdict re-run; neuter → no fire.
7. **The R1 port-kernel-dependency tooth.** `corrupt_port_kernel` (`OmegaU→2·OmegaU`) → `P0_raw` changes AND the ℓ=2
   determinacy (rank) row changes. Per-tooth: the perturbation genuinely propagates to `P0_raw` + the rank audit.
8. **The `.wl` INDEPENDENT-ROUTE + arity integrity (§1b).** Confirm the `.wl` nullspace is a MATERIALLY different route
   (constructive `NullSpace`/`Reduce`, OR native `MatrixRank` with the transliteration screen cleared) — NOT a line-by-line
   mirror of the `.py` `rankAuditFor`/`gateVerdictFor`/`dimOf`; def/call arity scan + unevaluated-leakage transcript scan (the
   `Module`s in `rankAuditFor`/`gateVerdictFor`/`dimOf`).

⚠ **NOT 023 (do not rebuild — 022 owns; 019 the prefactor):** the fingerprint core `spherical_j`/`spherical_y`/`dtn_branch`/
`build_fingerprints` (022), the Gate-4 non-regression + probe `3e`/`FAIL_QUAD_REGRESSION` (022 — consumed-True as provenance,
NOT re-run), the squared-denominator prefactor `P(ω)=D₀N/D^cons²` (019 DONE, STRIPPED). ⚠ **Vestigial `Delta`/`Sport` symbols —
DROP.**

## §6 Register expectation — ⭐ THE KEY 023 QUESTION (the ℓ=0/1 stiffnesses + the return admittances — CONFIRM + Codex-verify)

This is the sharpest register step in Part II. 023 introduces four new symbols into the counted-knob question — `K0c`, `K1 =
K_eta + 2·T_Omega`, `Z0_ret`, `Z1_ret`. The honest pre-read (⚠ CONFIRM at the register step + Codex-verify against the scripts):

- **⭐ `K0c` (ℓ=0 collective stiffness) — likely `DERIVED` (a manifestation of stage013's breathing packet), NOT a new counted
  knob.** Sourced to "handoff §8.2 / Gate-2 collective (delta_a,delta_L) reduction" = stage013's calibrated `(a,L)` collective
  closure; `K0c` is the ℓ=0 (monopole/breathing) scalar reduction of stage013's `K_AB`, so it is a DERIVED manifestation of the
  already-counted `{μ_η, T_w, β}` (edge R30-family). ⚠ CONFIRM it is not a genuinely-new independent stiffness.
- **⭐ `K1 = K_eta + 2·T_Omega` (ℓ=1 harmonic stiffness) — `DERIVED`, NOT new.** `K_eta = K_η` is R29-derived (`T_wβ²`, from
  {T_w,β}, 013); `T_Omega = T_Ω` is counted CALIB at 017; the `2` is `λ_m=ℓ(ℓ+1)=2` for ℓ=1 (016's covariance). So `K1` is a
  DERIVED combination of already-counted knobs (the ℓ=1 analogue of 017's `K₂=K̃+6·T̃_Ω`). A new DERIVED-manifestation edge.
- **⭐⭐ `Z0_ret`, `Z1_ret` (the ℓ=0/1 return admittances) — the CRUX: `FREE-UNREDUCED` reduction-debt WITH a named Gate-6
  route, NOT calibrated (they are UNDETERMINED — the whole point of the FAIL).** These are the 2 return admittances the native
  nullspace leaves free (`Z_is_premise=True`, pathA_29). They CANNOT be calibrated-to-a-benchmark (they are unpinned by the
  linear theory) — so they are NOT a clean CALIB headline knob. They are **reduction debt with a named-PENDING route: the Gate-6
  nonlinear closure supplying the selector `{Z0_ret=K0c, Z1_ret=K_eta+2·T_Omega}`** (the first concrete, proven Gate-6 input).
  ⚠ **This is exactly the "debt-WITH-a-route" the register's health metric tracks** (parameter_register §"two kinds of free" —
  vs the route-less `{ρ_B0,χ_c,C_hu}`). Register them as `FREE-UNREDUCED` (2 dofs) flagged with the Gate-6 selector route
  (a new edge, call it **R42** — confirm the next free number; R41 was 022's).
- **⭐ So 023 likely adds ZERO new *counted CALIB* knobs (K0c/K1 = DERIVED manifestations) but DOES add the `{Z0_ret, Z1_ret}`
  return-admittance FREE-UNREDUCED reduction-debt (2 dofs, Gate-6 route) + the DERIVED-manifestation rows (K0c, K1) + a new
  structural/obligation edge R42** (the cross-ℓ nullspace underdetermination departure + the Gate-6 `{Z0_ret, Z1_ret}` selector
  need + the `CROSS_L_RESIDUAL_PREDICTION` selector control; the honest FAIL landing). Part-II CALIB *headline* set likely stays
  = 6; the sim-deferred selector debt is the honest new liability-with-a-route. ⚠ **Codex-verify the classification** (the
  register verify is the gate that catches an over-count that falsely inflates, or a mislabel — e.g. `Z_ret` dressed as
  `DERIVED` when it is genuinely undetermined debt, or `K0c` counted afresh when it is a 013 manifestation — that would falsely
  shrink or inflate the irreducible codimension count).
- **Cited provenance (NOT re-counted):** stage022's `{1,1/2}`; the pathA_29 fold (009/010) residual form + the `Z_is_premise`;
  008's `M0/D1/R0/R1`; 017's port kernel + `T_Ω`; 013's `{μ_η,T_w,β}`; `c_s` (R1); `a` (`CONV`).
- **Control/tracked-not-counted:** `q_free` (free-carrier dim probe symbol — like 021's back-solve `μ̂₀`), `eta_null`
  (inject_null probe dof), `gain0/gain1` (decouple-knob control symbols) — all control-construction, tracked-not-counted (like
  `k_warp`/`α`).

⚠ **Do NOT let 023 silently (a) count `K0c`/`K1` as new CALIB (they are 013/017 DERIVED manifestations), nor (b) dress
`Z0_ret/Z1_ret` as `DERIVED`/calibrated (they are genuinely-undetermined reduction-debt — the FAIL's whole content), nor (c)
count the control symbols.** Resolve each and Codex-verify.

## Verdict tokens + honest scope

023 DELIVERS the joint **`FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** (2/2, COMPLETING; 022 landed the EARNED-first 1/2): the native
nullspace departure (dim-8 / return-nullity-2 over 11 genuine generator dofs, COMPUTED via `sp.Matrix(rows).rank()` / native
`MatrixRank` — no zero-padding), leaving the ℓ=0/1 return admittances `{Z0_ret, Z1_ret}` UNPINNED (`Z_is_premise`, pathA_29),
so `ε_eff` is free at the linear Gate-5 level → the return magnitude + nonzero prediction are DEFERRED (report `:47` `Deferred:`).
EARNED = the nullspace departure + the `A0/A1` residual form/sign/order (conditional on a positive bounded branch) + the
dimensional gate + the provenance firewall + the full able-to-fail probe suite; the selector-equation control
(`{Z0_ret=K0c, Z1_ret=K_eta+2·T_Omega}` → `CROSS_L_RESIDUAL_PREDICTION`) is what makes the gate able-to-fail and names the
first concrete Gate-6 input. ⭐ **023's SCRIPT passes (exit 0) by correctly COMPUTING the earned FAIL departure + every tooth
firing** (the stage003/pathA_36 earned-content-with-a-characterized-FAIL pattern). Consumes stage022's `{1,1/2}` (checkable
derive-vs-typed), the pathA_29 fold's residual form + `Z_is_premise` (keystone premise), 008's `M0/D1/R0/R1`, 017's port kernel,
013/017's ℓ=0/1 stiffnesses, `c_s`/`a` (units). Caveats: the return magnitude is genuinely undetermined (the FAIL); Gate-6's
nonlinear closure (sim-deferred) must supply the `{Z0_ret, Z1_ret}` selector; the `{Z0_ret, Z1_ret}` freedom is reduction-debt
with a named Gate-6 route (register §6). ⭐ **The pathA_34 v1 REJECTION locus is 023's (the sharpest) — the current source is
CLEAN (genuine rank, computed probes, native `.wl`); the reshape MUST preserve that: no zero-padded constraint, no hardcoded
`8`/`2` nullity, no stamped verdict token, no ε_eff/Z back-solve (the `FAIL_TAUTOLOGICAL` firewall), and the `.wl` must not
degrade from a native-rank engine into a token-asserting mirror.**

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G5b reshape directive (§1 the clean 023 slice / 2-way cut + the CONSUME-022-fingerprints-typed + the
STRIP-019-prefactor + §1b the `.wl` KEEP-native-vs-RE-AUTHOR decision [the KEY genuineness call — route through the bookend] +
§1c the checkable `{1,1/2}`→`A0/A1` consumption + the pathA_29-form + the `Z_is_premise` keystone + §2 faithful cover + §3 the
bridge-strip incl. sever-YAML + the 022/023 verdict-seam [023 owns the joint FAIL ladder, `quad_regression` consumed-from-022] +
transcript-level agreement + §5 the native-nullspace / selector-flip / `A0/A1` / dim-gate / firewall / probe teeth with per-tooth
ablation + §6 the register crux [K0c/K1 DERIVED, `{Z0_ret,Z1_ret}` FREE-UNREDUCED Gate-6-route-debt, R42 edge]) → **Codex xhigh
design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ FINAL Grok-4.5 headless compute-verify pass** (Grok SymPy-verifies the rank
audit `nullity=8`/`return-moving-nullity=2`, the selector flip to `CROSS_L_RESIDUAL_PREDICTION`, the `A0/A1` form + the `{1,1/2}`
consumption, the dim gate + free-carrier control, the firewall; it caught the 016 volume-vs-line + the 020 rule-inversion + the
022 pole-order-vs-inert mutant — so watch the nullspace-genuineness [genuine rank, not zero-padded], the selector-flip logic, the
back-solve firewall, and the `.wl` transliteration-screen decision) → assess + independently verify each catch → fold → **Codex
confirm-pass on the folds** → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`, background,
xhigh) → dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run → full tri-review (fidelity + adversarial-with-**per-
tooth ablation**; ⭐ hunt the nullspace genuine-rank-vs-rigged, the selector-flip able-to-fail, the `A0/A1` derive-vs-typed [not
a back-solve], the firewall genuine, a mirror-`.wl`, any vacuous able-to-fail) → remediate → fresh-agent re-verify → bump counts
22→23 → parameter register (⭐ K0c/K1 DERIVED, `{Z0_ret,Z1_ret}` FREE-UNREDUCED Gate-6-route-debt, R42 edge, control symbols
tracked-not-counted) + Codex-verify → note/card/`\input{stages/stage_023}` + registration → rebuild PDF → commit + docs/memory
sync (keep STATUS ▶ RESUME HERE THIN; append per-stage detail to part2 Progress). ⭐ **This COMPLETES the pathA_34 fold (022∧023)
+ closes the Gate-1–5 gravity ladder** → next = Cluster C (024–029), then the scheduled MIDWAY KNOB AUDIT. Orchestrator authors
notes/cards/LaTeX/registration/register; Codex codes. Target stem: `ledger_stage023_nullspace_underdetermination` (confirm slug
at directive authoring).
