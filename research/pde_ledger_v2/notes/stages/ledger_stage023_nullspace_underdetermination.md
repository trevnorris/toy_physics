# ledger_stage023 — the native nullspace underdetermination departure (Check II-G5b)

**Part / anchor.** Part II — Gravity (the frozen-throat cross-ℓ return sector). The COMPLETING, FAIL-delivering leg of a 2-way
split of `pathA_34`: stage 022 earned the cross-ℓ `−(ℓ+1)/Λ_ℓ` fingerprints + the Gate-4 non-regression (the EARNED-FIRST 1/2);
**this stage builds the native nullspace underdetermination that DELIVERS the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** — a
genuinely-computed nullspace over 11 generator dofs that leaves the ℓ=0/1 return admittances `{Z0_ret, Z1_ret}` unpinned, so the
linear theory cannot fix the ℓ=0/1↔ℓ=2 return. ⭐ **This COMPLETES the pathA_34 fold (022∧023) and closes the Gate-1–5 gravity
ladder** — the only remaining gravity item, Gate 6 (the nonlinear return-selector), is sim-deferred (a Part-VII open-item, R42).

**Verdict.** `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (JOINT, 2-stage) — **DELIVERED here** (2/2, COMPLETING; 022 landed the
EARNED-first 1/2). ⭐ **FAIL-headline stage (cf. stage 003 = earned transverse photons + the `FAIL_CAUCHY` departure): the
top-line verdict is a `FAIL_*`, but the earned content — a genuine native nullspace + the able-to-fail selector control — is the
headline, and the FAIL is the first-class characterized landing (no softening).** The script prints two distinct labels:
```
AUDIT_STATUS=PASS                                              (the script ran, all teeth fired, exit 0)
PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE (2/2, COMPLETING)   (the earned characterized-FAIL result)
```
The script PASSES (exit 0) by correctly COMPUTING the earned FAIL — the FAIL is the earned physics, not an execution failure.

**Status.** Symbolic / exact / float-free: the constraint Jacobian, its rank/nullspace, the return-transfer gradients, the
`A0/A1` residuals, and the `(L,M,T)`-triple dimensional algebra are exact `sympy` (`sp.Matrix(rows).rank()`) / native Wolfram
(`NullSpace`, `MatrixRank`) computations; `expect_zero`/`expect_bool`/`expect_fail` residual asserts, no `scipy`/`numpy`/floats/
tolerances. Dual-engine SymPy **111 PASS** / Mathematica **117 PASS**, exit 0, CWD-independent (post-remediation; the build's
116/123 became 111/117 after an honest de-count of 4 vacuous/bookkeeping teeth — §5).

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_34_cross_l_unification_{sympy.py,.wl}` +
> `reports/pathA_34_cross_l_unification.md` (the 023 slice = report :3, :5, :12, :15, :19–21, :24–27, :47) + the original
> directive `directives/pathA_34_cross_l_unification.md`. The report/directive are cited for provenance only; the derivation
> below is self-contained.

---

## 1. What this stage earns

The ℓ=0/1 brane↔bulk return (pathA_28/29 — `R0=−M0`, `R1=−D1`) and the ℓ=2 quadrupole (the grouped-P2 port kernel `P0_raw`) are
sectors of one linear return law. The earned result is that this law is **genuinely underdetermined**: a real nullspace computation
shows the ℓ=0/1 return admittances survive every collected Gate-5 constraint. A derived counterfactual selector collapses that
freedom, so the gate is able-to-fail.

### 1.1 The native nullspace (a genuine rank, not a rigged construction)
The linearized Gate-5 sector has **11 genuine generator dofs**
`[Ω_U, Ω_W, R_mix, g_U, g_W, D₀, K0c, K_eta, T_Omega, Z0_ret, Z1_ret]`. The three collected named constraints are the ℓ=2 port
kernel `P0_raw = (Ω_U²g_W + R_mix g_U)² / (Ω_U²Ω_W² − R_mix²)² / D₀` (handoff §10.2–10.3), the ℓ=0 collective stiffness `K0c`
(§8.2 / Gate-2), and the ℓ=1 harmonic stiffness `K1 = K_eta + 2·T_Omega` (§9.4). Their Jacobian over the 11 dofs (the row for
each constraint is `[∂c/∂dof]`) has
```
rank0 = rank( J ) = 3          (a genuine sp.Matrix(rows).rank() on symbolic ∂-rows — NO zero-padding, NO hardcoded 8/2)
native_nullspace_dimension = 11 − 3 = 8
```
Augmenting `J` with the two return-transfer gradients `∂T0, ∂T1` — where `T_ℓ = K_ℓ / (K_ℓ + Z_ℓ,ret)` is the DC return
transmission — gives
```
return_augmented_rank = rank( J ∪ {∂T0, ∂T1} ) = 5
return_moving_nullity = 5 − 3 = 2      ⟹  native_null_moves_return = True  ⟹  FAIL_UNDERDETERMINED_NOT_PREDICTIVE
```
So the two return admittances `{Z0_ret, Z1_ret}` are **untouched** by every collected named constraint (each constraint's
`∂/∂Z_ℓ,ret = 0`), yet they move `T0/T1` — the return is a genuine 2-parameter freedom. This is witnessed CONSTRUCTIVELY: the
unit vectors `e_{Z0_ret}` and `e_{Z1_ret}` each lie in the constraint nullspace (they preserve every constraint row) yet give
`ΔT0 = −K0c/(K0c+Z0_ret)² ≠ 0` and `ΔT1 = −(K_eta+2T_Omega)/(…)² ≠ 0` respectively. **Why the return survives (the keystone
premise):** pathA_29 records `Z_is_premise = True`, `boundary_dof = none` — the brane↔bulk return admittance is a *premise* of
the flat-slab family, not something the linear theory supplies an equation to fix.

### 1.2 The counterfactual selector control → `CROSS_L_RESIDUAL_PREDICTION` (the able-to-fail)
The gate is genuinely able-to-fail because a derived selector equation collapses the freedom. Adding the two selector rows
`{Z0_ret = K0c, Z1_ret = K_eta + 2·T_Omega}` as constraints raises the constraint rank 3→5 (native nullity 8→**6**) and makes
`T0 = T1 = 1/2` (constants with zero gradient), so
```
return_moving_nullity: 2 → 0      ⟹  native_moves = False  ⟹  CROSS_L_RESIDUAL_PREDICTION
```
Because `base_verdict` DEFAULTS to the predictive token `CROSS_L_RESIDUAL_PREDICTION` and only returns
`FAIL_UNDERDETERMINED_NOT_PREDICTIVE` when the computed `return_moving_nullity > 0`, **the FAIL is EARNED from the computed
nullity, not baked in.** ⚠ The selector `{Z0_ret=K0c, Z1_ret=K_eta+2·T_Omega}` is a **COUNTERFACTUAL RANK-COLLAPSE WITNESS** — it
is merely typed when the control mutation is chosen (`derived_from_named_pde = False`, `control_only = True`, `tautological =
False`); it is NOT a proven Gate-6 selector. The earned export is the *need*: Gate 6 must supply two independent equations
fixing the two return directions (or an equivalent return law).

### 1.3 The `A0/A1` scalar/dipole residuals vs pathA_29 (forward-built; consuming 022's `{1, 1/2}`)
The ℓ=0/1 return residual amplitudes are the cross-ℓ / port continuation of stage 009's `RETURN_RESIDUAL_PREDICTION`
(`A_res ∝ ε_ℓ/(1+ε_ℓ)`). They are built FORWARD from stage 022's earned radiative coefficients `v₀ = 1`, `v₁ = 1/2` and the
return transmission `(1 − T_ℓ)`:
```
A0 = i · v₀ · (aω/c_s)   · M0 · (1 − T0),      A1 = i · v₁ · (aω/c_s)³ · D1 · (1 − T1),
```
with `ε_ℓ = Z_ℓ,ret / K_ℓ` built FORWARD from the transfer definition (`1 − T_ℓ = ε_ℓ/(1+ε_ℓ)`) — **never** back-solved from a
residual. These are checked against the INDEPENDENT pathA_29 residual form
```
expected_A0 = i·aω·M0·ε0/(c_s(1+ε0)),      expected_A1 = i·a³ω³·D1·ε1/(2c_s³(1+ε1)),
```
and `A0 − expected_A0 = 0`, `A1 − expected_A1 = 0`, with `A0` at order `ω¹` and `A1` at order `ω³`. ⭐ **The consumption is
CHECKABLE:** the factor `1/2` in `v₁` is exactly what the `2·c_s³` in `expected_A1` encodes — algebraically
`A1 − expected_A1 = i·(a³ω³/c_s³)·D1·(ε1/(1+ε1))·(v₁ − 1/2)`, which is zero iff `v₁ = 1/2`. So corrupting the cited coefficient
(`v₁: 1/2 → 1/3`) fires the `A1_form` residual. The two sides have independent provenance — 022's fingerprint on one side, the
pathA_29 dipole form on the other — so the match is not an `X≡X`. The residual class is EARNED (form / sign / order, conditional
on a positive bounded branch); the *magnitude* is DEFERRED because `ε_eff` is left free by the nullspace.

### 1.4 The dimensional gate (the stage-021 `(L,M,T)`-triple machinery)
A recursive `dim_of` over `Mul/Pow/Add` (raising on a sum-mismatch) certifies, from the sourced dims (`[M0]=(0,1,−1)`,
`[D1]=(1,1,−1)`, `[Z_ℓ,ret]=[K_ℓ]=(0,1,−2)`, `[a]=(1,0,0)`, `[c_s]=(1,0,−1)`, `[ω]=(0,0,−1)`):
```
[A0] = (0,1,−1),   [A1] = (1,1,−1),   [T0]=[T1]=[ε0]=[ε1]=[P0_physical] = (0,0,0),   where [P0_physical] = (c_s/a)²·[P0_raw].
```
Corrupting a **sourced** dimension (`[M0] += (1,0,0)`) propagates to `[A0] → (1,1,−1)` → `FAIL_DIMENSIONAL`; corrupting an
**unconstrained free carrier** (`[q_free] = (7,0,0)`) gives `NO_FAIL` because `q_free` appears in no checked expression — the
free-carrier-independence control (which corruptions MAKE IT FIRE proves able-to-fail, per the stage-021 lesson).

### 1.5 The `FAIL_TAUTOLOGICAL` firewall (the anti-back-solve)
Each provenance item carries a `computed_class` and an `emitted_class`; the firewall is `class_matches_computed = (emitted ==
computed)`. While the sector is underdetermined, the `ε_eff` magnitude is CLASSED `deferred_branch_data` (not computed) — so
emitting it as derived is forbidden. Two mutations fire `FAIL_TAUTOLOGICAL`: `assert_not_derive` (emitting a genuinely-023-derived
object — the forward `T0/T1` transfer map — as an asserted literal) and a dedicated `emit_epsilon_magnitude_as_derived` mutation
(leaving the computed class `deferred_branch_data` but emitting a concrete `ε` magnitude). The baseline structurally asserts no
concrete `ε` magnitude is emitted. This is the firewall against the v1 rejection's `ε_eff/Z` back-solve.

---

## 2. The able-to-fail battery (023-owned)

023's own PHYSICS verdict IS the joint `base_verdict` ladder (the `quad_regression` rung is 022's — consumed `False`, not
re-run; the inert `able_to_fail_bad` rung was removed so a failed probe exits nonzero). Every tooth is a computed residual or a
DYNAMIC two-verdict `ablation` re-run:

| tooth | mutation → outcome | notes |
|---|---|---|
| native nullspace (nullity 8 / return-moving 2) | hardcode the rank, or zero-pad the Jacobian, or add a constraint touching `Z_ℓ,ret` → the computed nullity/return-moving changes → fires | genuine `sp.Matrix(rows).rank()` / `NullSpace`; forbidden to hardcode `8`/`2` |
| isolated rank teeth | `inject_null` (add `eta_null`) → native nullity 8→9; ONE selector row → return-moving 2→1 (native_moves STILL True); TWO independent rows → native_moves True→False | a single return row does NOT flip the verdict |
| selector-flip control | selector present → `return_moving_nullity=0` → `CROSS_L_RESIDUAL_PREDICTION`; neuter it (not touching `Z_ℓ,ret`) → no flip → fires | the counterfactual rank-collapse witness |
| `A0/A1` form + `{1,1/2}` consumption | corrupt cited `v₁` (1/2→1/3) → `A1_form` residual ≠ 0 → fires; corrupt the ω-order → `A1_order` fires | derive-vs-typed vs the INDEPENDENT pathA_29 form; never a back-solve |
| dimensional gate | corrupt `[M0]` → `FAIL_DIMENSIONAL`; corrupt `[q_free]` → `NO_FAIL` | sourced-corruption fires, free-carrier does not |
| firewall (`FAIL_TAUTOLOGICAL`) | `assert_not_derive` on the 023-derived T0/T1 map, asserted-selector, and `emit_epsilon_magnitude_as_derived` — 3 SEPARATE named assertions → each fires | the anti-back-solve; `ε` magnitude classed `deferred_branch_data` |
| transfer/return probes | `wrong_sign_return`→`FAIL_EPSILON_MISMATCH` (3c), `perfect_return`→`FAIL_OVERCANCEL` (3d), `no_consistent_return`→`FAIL_NO_CONSISTENT_RETURN` (3h), `decouple_knobs`→`FAIL_DECOUPLED` (3a) | each a DYNAMIC two-verdict re-run |
| R1 port-kernel dependency | `corrupt_port_kernel` (Ω_U→2·Ω_U) → `P0_raw` changes AND the ℓ=2 rank row changes → fires | proves the port kernel genuinely enters the audit |
| verdict read-set | wire a forbidden object into the verdict → the read-set assert fires | the ladder reads exactly the stage023 physics conditions |

De-counted (labeled prints, NOT verdict teeth): the provenance/premise documentation (`Z_is_premise`, the `cited_earned_input`
class of 022's `{1,1/2}`, the `ell2_P0_map` tags) and the `forward_relations_ok` T/ε self-consistency identity (§5).

---

## 3. Honest scope

- **The FAIL is the earned physics.** The genuinely-computed native nullspace (dim 8 / return-moving nullity 2, via a real rank)
  is the result: the LINEAR theory cannot pin the ℓ=0/1 return `{Z0_ret, Z1_ret}`. This is a first-class characterized departure
  (the stage-003 pattern), not a softened negative.
- **Gate 6 (sim-deferred) must supply the selector.** The honest export is the *need*: Gate 6's nonlinear closure must supply two
  independent equations fixing the two return directions (or an equivalent return law). The selector `{Z0_ret=K0c,
  Z1_ret=K_eta+2·T_Omega}` is a COUNTERFACTUAL rank-collapse WITNESS that the gate is able-to-fail — NOT a proven Gate-6 selector
  (deriving those exact equalities is separate Gate-6 work).
- **The raw nullity 8 is partly bookkeeping; the verdict rides return-moving 2.** Two of the three collected constraints (`K0c`,
  `K1`) are self-constraint rows (each self-pins its own stiffness), so the raw native nullity 8 carries some bookkeeping — using
  `P0_raw` alone gives rank 1 / nullity 10, but the **return-moving nullity is still 2**. The verdict-bearing quantity is
  `return_moving_nullity`; the raw nullity is a reported diagnostic.
- **The return magnitude is DEFERRED.** The `A0/A1` residual *class* (form / sign / order, conditional on a positive bounded
  branch) is earned; the *magnitude* and the nonzero prediction are deferred because the native nullspace leaves `ε_eff` free at
  the linear Gate-5 level. `ε_eff` is FORWARD-built (`ε = Z/K`), never back-solved (the `FAIL_TAUTOLOGICAL` firewall).

---

## 4. Consumed / exported

- **Consumed — two CHECKABLE relations + PROVENANCE cites.**
  - **stage 022's ℓ=0/1 radiative coefficients `{1, 1/2}` — the CHECKABLE derive-vs-typed feeding `A0/A1`** (§1.3): corrupting the
    cited `v₁` fires the `A1_form` residual (the `2c_s³` in `expected_A1` encodes the `1/2`). Cited typed — 022 owns the
    fingerprints; 023 does not re-run the fingerprint battery.
  - **the pathA_29 residual form (stage 009/010, `A_res ∝ ε_ℓ/(1+ε_ℓ)`, `Z=−M0·ε0/(1+ε0)`)** — the second checkable relation;
    `expected_A0/A1` are built from this INDEPENDENT form, not typed to match `A0/A1`.
  - **pathA_29's `Z_is_premise = True`, `boundary_dof = none`** — the keystone premise for why the underdetermination is EARNED
    (the linear theory supplies no fixing equation).
  - **008's `R0=−M0` / `R1=−D1` targets + the `M0`/`D1` moments** — `M0/D1` are the source amplitudes in `A0/A1`; `R0/R1` are
    external-bridge dim tags. Cited PROVENANCE.
  - **017's grouped-P2 port kernel `P0_raw`** — the ℓ=2 named-constraint row + the R1 probe target.
  - **013's ℓ=0 collective + 017's ℓ=1 harmonic sectors** — cited as provenance CONTEXT for the effective stiffnesses `K0c`/`K1`
    (their reduction is PENDING — see register).
  - **`c_s`** (stage 005 R1 `c_s²=5Kρ⁴/m`) — the units carrier in `z = aω/c_s`; **`a`** — the `CONV` pin. (Distinct from the
    frozen-wall Helmholtz speed `c_S`, 011–017.)
- **Register.** **ZERO new counted CALIB knobs (set stays 6); `Z0_ret/Z1_ret` add ZERO new free dofs (aliases); but `K0c/K1`
  add COUNTED `FREE-UNREDUCED` PENDING reduction-debt.** The return admittances `{Z0_ret, Z1_ret}` are COORDINATE
  ALIASES of the existing `ε0/ε1` FREE-UNREDUCED debt (register row for stage 009; `ε_ℓ = Z_ℓ,ret/K_ℓ` invertible once `K_ℓ`
  fixed) — not two new freedoms, no double-count. ⚠ The effective ℓ=0/1 stiffnesses `K0c` and the ℓ=1 sector `{K_eta, T_Omega}`
  (appearing as `K1 = K_eta + 2·T_Omega`) are pathA_34-convention scalars with dims `M T⁻²` (the scripts' `(L,M,T)`-tuple
  `(0,1,−2)`); they are classed **`FREE-UNREDUCED`, `PENDING` scalar-reduction, and COUNTED as reduction debt** (per the
  register's rule pending debt stays counted until DERIVED) — NOT `DERIVED`, NOT new `CALIB`: their dims
  do NOT match registered stage 013 `K_η=T_wβ²` (`M L⁻¹T⁻²`) or stage 017 `T_Ω` (`M L⁻³T⁻²`), and stage 017 records
  `K_η=T_wβ²` as non-transferable across the volume-vs-line convention (the stage-016 lesson) — so they are NOT identified with
  the raw 013/017 densities (an explicit profile+measure scalar-reduction to the wall packet would be needed to earn `DERIVED`
  and discharge the debt). The control symbols `q_free`/`eta_null`/`gain0`/`gain1` are tracked-not-counted. New obligation edge
  **R42** (the cross-ℓ nullspace underdetermination departure + the sharpened Gate-6 return-selector obligation; it SHARPENS the
  existing `ε0/ε1` R24-family debt — "the linear return law leaves exactly a 2-dim return-admittance nullspace" — and adds no
  free dofs). Part-II counted CALIB set unchanged at `{μ_η, T_w, β}`(013) + `{Vp0/ℓ_c}`(015) + `{T_Ω, β₂}`(017) = 6. ⚠ The Codex
  register-verify (a 4th check, post-build) caught that "zero new free dofs" over-claimed for `K0c/K1` (pending reduction-debt
  stays counted) + a dim-tuple-convention nit — both folded → `REGISTER_CLEAN`.
- **Exported.** The native nullspace departure (dim-8 / return-moving-nullity-2) + the **Gate-6 need for two independent return
  equations** fixing `{Z0_ret, Z1_ret}` (the sharpened R42 obligation → the Part-VII open-items register + the Gate-6 caveat on
  Cluster-C 024/027/028) + the `A0/A1` residual class + the completed pathA_34 joint (022∧023). ⭐ **This COMPLETES the pathA_34
  fold and closes the Gate-1–5 gravity ladder.**

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. The `.wl` is a **genuinely independent route,
RE-AUTHORED** (not a kept transliteration): the source `.wl`'s `rankAuditFor` block was structurally parallel to the `.py`'s
`augRank−rank0` (same dof order, same constraint decomposition), so it was discarded per the mirror-policy transliteration
screen. The re-authored `.wl` proves the nullspace CONSTRUCTIVELY with native Wolfram primitives — `Length[NullSpace[Jbase]] = 8`,
`MatrixRank[basis·Transpose[Greturn]] = 2` (the return-moving dimension read directly off the null basis), the explicit
`Z0_ret`/`Z1_ret` unit directions lie in `NullSpace[Jbase]`, and (with the selector) `Length[NullSpace[Jselector]] = 6` with the
original unsubstituted return gradients giving `Greturn · Transpose[NullSpace[Jselector]] = 0` — a materially different algorithm
from the `.py`'s rank-subtraction. Agreement is transcript-level (both engines emit `native_nullspace_dimension=8`,
`return_moving_nullity=2`, `baseline_verdict=FAIL_UNDERDETERMINED_NOT_PREDICTIVE`, `selector_control_verdict=CROSS_L_RESIDUAL_PREDICTION`,
the `A0/A1` residual-zero, the dim gate's sourced-FAIL / free-carrier-NO_FAIL, and the 3a/3c/3d/3f/3g/3h + R1 probe tokens). The
stage-007 unevaluated-leakage failure mode is guarded (arity self-check + transcript scan).

**Directive review** used the Codex→Grok→Codex bookend — and this stage is the pathA_34 v1 REJECTION locus (the sharpest in
Part II: v1 was rejected for a rigged zero-padded nullity, flag-driven probes, and a headline-only `.wl`). Codex's design-review
returned **7 BLOCKING**, all genuine and folded: the `.wl` must be decisively RE-AUTHORED via constructive `NullSpace`; the
selector collapses the RETURN-MOVING nullity to 0, not "the nullity" (the full native nullity goes 8→6, rank 3→5) — a directive
math fix, with corrected isolated ablation teeth; the selector is a counterfactual rank-collapse WITNESS, not a proven Gate-6
selector; the provenance cut must class 022's `{1,1/2}` as `cited_earned_input` and rewire `assert_not_derive` to a
genuinely-023-derived object (dropping the `gate4_prefactor` tag); the firewall must add the `emit_epsilon_magnitude_as_derived`
tooth, de-count the literal `rerun_gate_logic`, and fix the inert `able_to_fail_bad`; the register must class `K0c/K1` as
`PENDING` scalar-reduction (their dims `(0,1,−2)` do not match 013/017 — the stage-016 convention trap), NOT `DERIVED`; and it
must register `Z0_ret/Z1_ret` as coordinate aliases of the existing `ε0/ε1` debt, not two new free dofs. Codex confirm-passes
folded the remaining consistency-sweep gaps. A **Grok-4.5 compute-verify** of the folded directive returned `DIRECTIVE_CLEAN`,
independently confirming (with its own SymPy) the rank 3 / nullity 8 / return-moving 2, the selector flip to
`CROSS_L_RESIDUAL_PREDICTION` with native nullity 8→6, `A1 − expected_A1 = 0` iff `v₁ = 1/2`, the dims, and the `K0c/K1`
dim-conflict + `Z_ret` alias conventions; it validated the `.wl` constructive route (including `Greturn·NullSpace[Jselector] = 0`
as a genuine identity) and added one non-blocking honest-scope note (the raw nullity 8 includes `K0c/K1` self-constraint
bookkeeping; the verdict rides return-moving 2). A closing Codex confirm closed the bookend `DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents (arbiter re-run reproduced the build SymPy 116 / Mathematica 123, exit 0, CWD-independent):
`FIDELITY_CLEAN` (an independent read hand-re-derived the rank audit, the selector flip, the `A0/A1` forward-build + the `v₁=1/2`
consumption, and the dims with its own SymPy — all faithful; the rank is a genuine `sp.Matrix(rows).rank()` not zero-padded,
`ε_eff` is forward-built with no back-solve, and the `.wl` is a materially-different constructive engine) + `ADVERSARIAL_CLEAN`
(per-tooth ablation, 15 mutations across both engines: hardcoding the rank, zero-padding the Jacobian, and faking the `.wl`
`NullSpace` basis all make the audit FAIL; the `emit_epsilon` firewall rides a real class-mismatch; the `{1,1/2}` consumption and
the independent `expected_A1` form both fire; the `.wl` genuinely computes `NullSpace`, not a token-asserting mirror — the four
KEY anti-rig properties confirmed; **4 non-blocking de-count nits**). Codex remediated: **2 make-genuine** (the witness-preservation
assert now recomputes each Jacobian-row dot product from the stored witness vector; the neutralized-mutation meta-test now uses a
cache-distinct inert context with an independence check, defeating the `compare=False` name-collapse) + **2 de-count** (the
provenance/premise documentation asserts → labeled `PROVENANCE` prints; the `forward_relations_ok` T/ε identity → a labeled
`SELF-CONSISTENCY` check). Fresh-agent `REVERIFY_CLEAN` by the coupling meta-test (each made-genuine tooth fires under a mutation
of its object and goes vacuous when neutered, both engines; the de-counts keep coverage; no KEY earned tooth regressed). Tallies
116/123 → **111/117** (net −5/−6 per engine from the honest de-counts). Symbolic per-tooth ablation, mutations on copies. ⭐ With
this stage the pathA_34 fold is COMPLETE (022∧023) and the Gate-1–5 gravity ladder is closed.
