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

### 1.6 Dimension-object enumeration (step (a), frozen before .wl emission)

**Membership rule and working.** One object is counted for each stable source binding or addressable association member whose
value is an `(L,M,T)` exponent vector. Distinct symbol keys remain distinct declarations even when their values alias `zeroDim`.
A container is not counted again on top of its clean per-key members. A local alias is represented by the stable returned/top-level
path that survives it, not counted twice. A mutation map or mutation result is counted once as a whole state-specific object,
including its replacement binding and its unchanged copied members, rather than re-counting per-key members. The memoized
`caseFor[mode]` result is one conditional object family with all reached modes stated, because it is written at one source locus
and its provenance changes by `mode`. A physical expression and a transient numeric atom visited while walking it are not stable
dimension objects; the vector returned by the walk is. This cut is deliberately reviewable: a reviewer can reject, for example,
the choice to treat the 14 memoized mode instances as one conditional family rather than 14 rows.

**Search predicates fixed before the results.** The `.wl` was searched at its complete 769-line pre-emission
revision; it now stands at **816 lines**, the 47 added being the step-(d) dimension-record emitter alone —
`dimensionAxisSlots`/`dimensionAxesLabel`/`dimensionComponents`/`printDimRecord` (`.wl:212-220`),
`emitDimensionRecords[]` (`.wl:298-332`) and its `runAll` call (`.wl:789`) — which adds no `baseDims` or
`expectedDims` entry and no new walked result, only re-reading the 29 records proposed below. The four finite
predicates were:
`P1`, every exact three-component integer/rational list literal; `P2`, every identifier containing case-insensitive `dim`;
`P3`, every exact occurrence of the seven `dimensionAudit` record keys; and `P4`, every identifier assignment to an exact
integer/rational scalar, including candidate coefficients, counters, exponents and print-text false positives. The re-runnable
search and every match are in
`research/pde_ledger_v2/_scratch/stage023/enum/candidate_ledger.md`. It produced 179 match-level hits: `P1=23`, `P2=113`,
`P3=25`, `P4=18`. Exactly 52 hits are disposed as `PROMOTE` to one of the 42 row keys below and 127 as `RULE_OUT`; redundant
predicates explain
why promoted-hit count exceeds row count. Every one of those 179 hits has exactly one disposition, and every one of the 42 rows
below has at least one promoting hit; the scope of both universals is the `P1`–`P4` hit set recorded in that ledger and the row
table in this section, not the `.wl`'s dimension content at large.

This ledger makes omission reviewable: a reviewer can re-run the predicates, see every matched locus and challenge any
disposition. It does **not** prove the predicate set exhaustive. It would systematically miss a dimension vector assembled
without an exact three-literal and stored under a neutral identifier containing no `dim`, an opaque helper return assigned to
such a name, or an unlabelled vector imported from elsewhere.

**Classes defined from this file.** `neutral helper vector` is the walker's stable identity/alias source; `sourced declaration`
is a clean symbol→vector entry; `expected target` is an independently typed comparison/fallback oracle; `clean walked result` is
a vector computed from a live baseline expression; `mutation-only declaration map` and `mutation-only walked result set` are
deliberate probe states; and `memoized conditional walked result family` is the returned per-mode set whose provenance depends on
control flow. `Transient numeric fall-through` is a rule-out mechanism, not a row class: `.wl:224` returns `zeroDim` for the
evaluated numeric atom the walker visits, not for a surviving source coefficient binding. These cuts distinguish how a stable
vector came into existence (declaration, target, or walk) and whether its lifetime is clean, helper-only, mutation-only, or
conditional; collapsing those kinds would hide self-comparison and corruption paths.

<!-- B4_DETERMINATION_BEGIN -->
**B4 determination — `transient numeric fall-through`, complete membership and disposition.** The row membership of this
rule-out class is empty. `v0`, `v1`, and `coeff1` are all `RULE_OUT`, consistently: `v0=1` at `.wl:197` is absorbed before the
walk (the committed expression at
`research/pde_ledger_v2/mathematica/out/ledger_stage023_nullspace_underdetermination_mathematica_audit.out:86` contains no `1`
factor); `v1=1/2` at `.wl:198` is fused with `I` into the numeric atom `I/2` shown at the same transcript `:87`; and local
`coeff1=1/2` at `.wl:432` feeds the structurally identical `I coeff1` product at `.wl:453` before `dimensionAudit` walks it at
`.wl:461-463`. Thus the mechanism that actually occurs is that `dimOf` visits an evaluated/fused numeric atom and assigns that
atom `zeroDim` at `.wl:224`; none of the three scalar bindings survives as a stable exponent-vector object. The search still
records all three as `P4-H003`–`P4-H005`, and the source-context disposition derivation applies the same rule to each.
<!-- B4_DETERMINATION_END -->

**Quoted source reconstruction (the prior survey is not needed for any of the 42 rows tabled below).**

```text
.wl:209       zeroDim = {0, 0, 0};
.wl:222-225   dimOf[expr_, dims_Association] := Module[... Which[
                  TrueQ[expr == 0] || NumericQ[expr], zeroDim,
                  AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr], ...
.wl:247-255   baseDims = <|
                  a -> {1,0,0}, cs -> {1,0,-1}, omega -> {0,0,-1},
                  M0 -> {0,1,-1}, D1 -> {1,1,-1}, R0 -> {0,1,-1}, R1 -> {1,1,-1},
                  D0 -> {-1,1,-2}, K0c -> {0,1,-2}, Keta -> {0,1,-2},
                  TOmega -> {0,1,-2}, Z0ret -> {0,1,-2}, Z1ret -> {0,1,-2},
                  OmegaU -> {0,0,-1}, OmegaW -> {0,0,-1}, Rmix -> {0,0,-2},
                  gU -> {-1/2,1/2,-2}, gW -> {-1/2,1/2,-2},
                  etaNull -> zeroDim, gain0 -> zeroDim, gain1 -> zeroDim, qfree -> zeroDim |>;
.wl:257-262   expectedDims = <|"A0"->{0,1,-1}, "A1"->{1,1,-1},
                  "T0"->zeroDim, "T1"->zeroDim, "epsilon0"->zeroDim,
                  "epsilon1"->zeroDim, "P0Physical"->zeroDim|>;
.wl:270-280   computed = AssociationMap[
                  If[TrueQ[expressions[name] == 0], expectedDims[name],
                     dimOf[expressions[name], dims]] &, Keys[expressions]];
                ... "Computed" -> computed, "Expected" -> expectedDims, ...
.wl:286-288   baselineDimAudit = dimensionAudit[
                  baseDims, A0lead, A1lead, T0dc, T1dc, eps0, eps1,
                  baselinePort["P0Physical"]];
.wl:289-296   corruptSourcedDims = Join[KeyDrop[baseDims,{M0}], <|M0->{1,1,-1}|>];
                corruptSourcedDimAudit = dimensionAudit[corruptSourcedDims, ...];
                corruptFreeDims = Join[KeyDrop[baseDims,{qfree}], <|qfree->{7,0,0}|>];
                corruptFreeDimAudit = dimensionAudit[corruptFreeDims, ...];
.wl:424,461-463,490-492
                caseFor[mode_String] := caseFor[mode] = Module[...;
                  dimAudit = dimensionAudit[dims,a0,a1,t0,t1,e0,e1,
                    baselinePort["P0Physical"]]; ...;
                  "Dimension" -> dimAudit, ...];
```

**B1 determination — group-hazard completeness and honesty.** Every one of the **29 `PROPOSE EMIT` rows** carries an explicit
`B1 flag`; the 13 `PROPOSE NOT EMIT` rows carry none, because nothing is proposed for emission there to group. The universal is
scoped to the emission-proposed subset of the 42-row table, not to all 42 rows. The flagged
set is `sourced_dims.a`, `sourced_dims.c_s`, `sourced_dims.M0`, `sourced_dims.D1`, `sourced_dims.R0`, `sourced_dims.R1`,
`sourced_dims.D0`, `sourced_dims.K_eta`, `sourced_dims.T_Omega`, `computed_dims.epsilon0`,
`computed_dims.epsilon1`, and `computed_dims.P0_physical`; the live hazards are `a`, `K_eta`, and `T_Omega`, and the other nine
are forward hazards. The remaining 17 names are explicitly `NOT FLAGGED`.
Each flagged row names the existing or future member, its evidence, and the group it will create. This section records only the
proposal, the resulting group, and this stage's own identity text; it makes no tenability ruling. The physics leg in §1.7 owns
that determination. Canonical evidence is a FROZEN PRE-EMISSION SNAPSHOT of `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md`
**at revision `8c4b25b0`** — that version's census (30-stage corpus, **6 of 30** converted, **93** quantity rows, 5 candidate
groups, 2 `NEEDS_ADJUDICATION`) plus its full quantity-row and candidate-group tables; it predates this stage's own 29 rows and
therefore omits them. ⚠ The file **has since been regenerated** (122 rows, 7 of 30), so no live line range reaches that
snapshot — cite the revision. Stage evidence is the pre-existing opening/§4 and the new §1.7.

**Per-object record.**

<!-- ENUM_ROW_TABLE_BEGIN -->
| row key | object | definition locus | class | consumer locus | artifact status | reason |
|---|---|---|---|---|---|---|
| `zeroDim` | `zeroDim` | `.wl:209` | neutral helper vector | `.wl:224,237,254,259-261` | PROPOSE NOT EMIT | Private neutral element and alias source, not a declaration for a distinct quantity or a clean walked result. |
| `sourced.a` | `baseDims[a]` | `.wl:247-248` | sourced declaration | `.wl:116,199-202,225,286-295,462` | PROPOSE EMIT as `sourced_dims.a` | Record denotes the clean declared dimension of `a`; both engines spell the identifier `a`. B1 flag: GROUP HAZARD; it will join stage016 `dim_rules.a` (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:127`) and stage018 `a` (`:133`), already grouped at `:187` (**the `a` candidate-key group**); this stage calls `a` the `CONV` pin (§4). The physics leg in §1.7 governs the identity decision. |
| `sourced.cs` | `baseDims[cs]` | `.wl:247-248` | sourced declaration | `.wl:116,199-202,225,286-295,462` | PROPOSE EMIT as `sourced_dims.c_s` | Record denotes the clean declared dimension of Wolfram `cs`; the SymPy symbol name is `c_s`, proposed as join spelling. B1 flag: FORWARD GROUP HAZARD; future stages 011–017 `c_S` normalizes to the same `cS` key, while this stage calls `c_s` its units carrier and records it as distinct from that frozen-wall speed (§4). This records the proposal and evidence only; §1.7(2)/(2b) owns the determination. |
| `sourced.omega` | `baseDims[omega]` | `.wl:247-248` | sourced declaration | `.wl:199-202,225,286-295,462` | PROPOSE EMIT as `sourced_dims.omega` | Record denotes the clean declared dimension of `omega`; both engines spell the identifier `omega`. B1 flag: NOT FLAGGED; §1.7(2) identifies it with represented stages 011/012, whose emitted names are `OmegaDim` and `omega_dim` (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**)) and therefore do not share candidate key `omega`; no unrepresented same-key counterpart is stated. |
| `sourced.M0` | `baseDims[M0]` | `.wl:247-249` | sourced declaration | `.wl:199,201,225,286-295,452,454,461-463` | PROPOSE EMIT as `sourced_dims.M0` | Record denotes the clean declared dimension of `M0`; both engines spell the identifier `M0`. B1 flag: FORWARD GROUP HAZARD; stages 008/009 are not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and this stage identifies the `M0` moment as the same object flowing 008→009/010→022/023 (§4; §1.7(2)); a future `M0` record will meet key `M0`. The physics leg in §1.7 governs the identity decision. |
| `sourced.D1` | `baseDims[D1]` | `.wl:247-249` | sourced declaration | `.wl:200,202,225,286-295,453,455,461-463` | PROPOSE EMIT as `sourced_dims.D1` | Record denotes the clean declared dimension of `D1`; both engines spell the identifier `D1`. B1 flag: FORWARD GROUP HAZARD; stages 008/009 are not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and this stage identifies the `D1` moment as the same object flowing 008→009/010→022/023 (§4; §1.7(2)); a future `D1` record will meet key `D1`. The physics leg in §1.7 governs the identity decision. |
| `sourced.R0` | `baseDims[R0]` | `.wl:247-249` | sourced declaration | NONE_FOUND (searched `.wl:1-816`; symbol `R0` occurs at `.wl:28,249,305` — the `$Assumptions` list, this declaration, and the step-(d) emitter's read of `baseDims[R0]` — while the tokens at `.wl:355,768` are a provenance tag string and print text; none of the three symbol loci enters the seven expressions at `.wl:286-295,461-463`) | PROPOSE EMIT as `sourced_dims.R0` | Record denotes the clean declared dimension of `R0`; both engines spell it `R0`. Emission exposes non-consumption. B1 flag: FORWARD GROUP HAZARD; stage008 is not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and this stage consumes its `R0=−M0` target under the same `R0` name (opening; §4); a future stage008 `R0` record will meet key `R0`. The physics leg in §1.7 governs the identity decision. |
| `sourced.R1` | `baseDims[R1]` | `.wl:247-249` | sourced declaration | NONE_FOUND (searched `.wl:1-816`; symbol `R1` occurs at `.wl:28,249,306` — the `$Assumptions` list, this declaration, and the step-(d) emitter's read of `baseDims[R1]`; textual `R1` at `.wl:355,710,719-723,768,770` denotes probe/provenance labels, not this symbol; none of the three symbol loci enters the seven expressions at `.wl:286-295,461-463`) | PROPOSE EMIT as `sourced_dims.R1` | Record denotes the clean declared dimension of `R1`; both engines spell it `R1`. Emission exposes non-consumption. B1 flag: FORWARD GROUP HAZARD; stage008 is not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and this stage consumes its `R1=−D1` target under the same `R1` name (opening; §4); a future stage008 `R1` record will meet key `R1`. The physics leg in §1.7 governs the identity decision. |
| `sourced.D0` | `baseDims[D0]` | `.wl:247-250` | sourced declaration | `.wl:105-116,120-121,225,286-295,462` | PROPOSE EMIT as `sourced_dims.D0` | Record denotes the clean declared dimension of `D0`; both engines spell the identifier `D0`. B1 flag: FORWARD GROUP HAZARD; stages 021 and 027 are not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and §1.7(2) identifies their `D0` with this one; future stage021/027 `D0` records will meet key `D0`. The physics leg in §1.7 governs the identity decision. |
| `sourced.K0c` | `baseDims[K0c]` | `.wl:247-250` | sourced declaration | `.wl:138-141,199,201,225,286-295,447-462` | PROPOSE EMIT as `sourced_dims.K0c` | Record denotes the clean declared dimension of `K0c`; both engines spell the identifier `K0c`. B1 flag: NOT FLAGGED; no current `K0c` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §1.7(2) explicitly distinguishes it from stage013 `k_shared`, whose candidate key differs; no future same-key counterpart is stated. |
| `sourced.Keta` | `baseDims[Keta]` | `.wl:247-250` | sourced declaration | `.wl:122,139,141,200,202,225,286-295,448-462` | PROPOSE EMIT as `sourced_dims.K_eta` | Record denotes the clean declared dimension of Wolfram `Keta`; the SymPy symbol name is `K_eta`, proposed as join spelling. B1 flag: GROUP HAZARD; it will join stage013 `K_eta` (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:97`) and stage016 `dim_rules.K_eta` (`:121`), whose existing group is `NEEDS_ADJUDICATION` (`:184`, **the `KEta` candidate-key group**); this stage says its scalar is not identified with those raw densities (§4). This records the proposal, group, and evidence only; §1.7(2)/(2b) owns the determination. |
| `sourced.TOmega` | `baseDims[TOmega]` | `.wl:247,250-251` | sourced declaration | `.wl:122,139,141,200,202,225,286-295,448-462` | PROPOSE EMIT as `sourced_dims.T_Omega` | Record denotes the clean declared dimension of Wolfram `TOmega`; the SymPy symbol name is `T_Omega`, proposed as join spelling. B1 flag: GROUP HAZARD; it will join stage016 `dim_rules.T_Omega` (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:124`); this stage says its scalar is not identified with the raw stage017 density (§4). This records the proposal, group, and evidence only; §1.7(2)/(2b) owns the determination. |
| `sourced.Z0ret` | `baseDims[Z0ret]` | `.wl:247,251` | sourced declaration | `.wl:138,140,199,201,225,286-295,447-462` | PROPOSE EMIT as `sourced_dims.Z0_ret` | Record denotes the clean declared dimension of Wolfram `Z0ret`; the SymPy symbol name is `Z0_ret`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `Z0Ret` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §1.7(2) distinguishes this return admittance from stage009 `Z`; no future `Z0_ret` counterpart record is stated. |
| `sourced.Z1ret` | `baseDims[Z1ret]` | `.wl:247,251` | sourced declaration | `.wl:139,141,200,202,225,286-295,448-462` | PROPOSE EMIT as `sourced_dims.Z1_ret` | Record denotes the clean declared dimension of Wolfram `Z1ret`; the SymPy symbol name is `Z1_ret`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `Z1Ret` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §1.7(2) distinguishes this return admittance from stage009 `Z`; no future `Z1_ret` counterpart record is stated. |
| `sourced.OmegaU` | `baseDims[OmegaU]` | `.wl:247,252` | sourced declaration | `.wl:105-120,225,286-295,462` | PROPOSE EMIT as `sourced_dims.Omega_U` | Record denotes the clean declared dimension of Wolfram `OmegaU`; the SymPy symbol name is `Omega_U`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `OmegaU` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and this stage/§1.7 uses it only as an internal port generator without naming a same-key cross-stage record. |
| `sourced.OmegaW` | `baseDims[OmegaW]` | `.wl:247,252` | sourced declaration | `.wl:105-120,225,286-295,462` | PROPOSE EMIT as `sourced_dims.Omega_W` | Record denotes the clean declared dimension of Wolfram `OmegaW`; the SymPy symbol name is `Omega_W`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `OmegaW` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and this stage/§1.7 uses it only as an internal port generator without naming a same-key cross-stage record. |
| `sourced.Rmix` | `baseDims[Rmix]` | `.wl:247,252` | sourced declaration | `.wl:105-120,225,286-295,462` | PROPOSE EMIT as `sourced_dims.R_mix` | Record denotes the clean declared dimension of Wolfram `Rmix`; the SymPy symbol name is `R_mix`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `RMix` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and this stage/§1.7 uses it only as an internal port generator without naming a same-key cross-stage record. |
| `sourced.gU` | `baseDims[gU]` | `.wl:247,253` | sourced declaration | `.wl:105-120,225,286-295,462` | PROPOSE EMIT as `sourced_dims.g_U` | Record denotes the clean declared dimension of Wolfram `gU`; the SymPy symbol name is `g_U`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `gU` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and this stage/§1.7 uses it only as an internal port generator without naming a same-key cross-stage record. |
| `sourced.gW` | `baseDims[gW]` | `.wl:247,253` | sourced declaration | `.wl:105-120,225,286-295,462` | PROPOSE EMIT as `sourced_dims.g_W` | Record denotes the clean declared dimension of Wolfram `gW`; the SymPy symbol name is `g_W`, proposed as join spelling. B1 flag: NOT FLAGGED; no current `gW` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and this stage/§1.7 uses it only as an internal port generator without naming a same-key cross-stage record. |
| `sourced.etaNull` | `baseDims[etaNull]` | `.wl:247,254` | sourced declaration | NONE_FOUND (searched `.wl:1-816`; `etaNull` enters injected physics at `.wl:189-194` but no expression passed to `dimensionAudit` at `.wl:286-295,461-463`; the step-(d) emitter reads `baseDims[etaNull]` at `.wl:318` without walking it) | PROPOSE EMIT as `sourced_dims.eta_null` | Record denotes the clean declared dimension of Wolfram `etaNull`; the SymPy symbol name is `eta_null`. Emission exposes dimensional non-consumption. B1 flag: NOT FLAGGED; no current `etaNull` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §4/§1.7 classify it as a local tracked-not-counted control without a cross-stage identity. |
| `sourced.gain0` | `baseDims[gain0]` | `.wl:247,254` | sourced declaration | `.wl:451,461-463` through `dimOf` at `.wl:225` in mode `decouple` | PROPOSE EMIT as `sourced_dims.gain0` | Record denotes the clean declared dimension of `gain0`; both engines spell the identifier `gain0`. B1 flag: NOT FLAGGED; no current `gain0` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §4/§1.7 classify it as a local tracked-not-counted control without a cross-stage identity. |
| `sourced.gain1` | `baseDims[gain1]` | `.wl:247,254` | sourced declaration | `.wl:451,461-463` through `dimOf` at `.wl:225` in mode `decouple` | PROPOSE EMIT as `sourced_dims.gain1` | Record denotes the clean declared dimension of `gain1`; both engines spell the identifier `gain1`. B1 flag: NOT FLAGGED; no current `gain1` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §4/§1.7 classify it as a local tracked-not-counted control without a cross-stage identity. |
| `sourced.qfree` | `baseDims[qfree]` | `.wl:247,254` | sourced declaration | NONE_FOUND (searched `.wl:1-816`; `.wl:282,668` positively establish absence from all checked expressions; the step-(d) emitter reads `baseDims[qfree]` at `.wl:321` without walking it) | PROPOSE EMIT as `sourced_dims.q_free` | Record denotes the clean declared dimension of Wolfram `qfree`; the SymPy symbol name is `q_free`. Emission exposes the deliberately unconsumed control declaration. B1 flag: NOT FLAGGED; no current `qFree` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §1.7(1) calls it an undetermined free-carrier control without a cross-stage identity. |
| `expected.A0` | `expectedDims["A0"]` | `.wl:257-258` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration and comparison oracle; `.wl:276,655-661` compares the walked `A0` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `expected.A1` | `expectedDims["A1"]` | `.wl:257-258` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration and comparison oracle; `.wl:276,655-661` compares the walked `A1` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `expected.T0` | `expectedDims["T0"]` | `.wl:257-259` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration, comparison oracle, and zero-expression fallback; `.wl:276,655-661` compares the walked `T0` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `expected.T1` | `expectedDims["T1"]` | `.wl:257-259` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration, comparison oracle, and zero-expression fallback; `.wl:276,655-661` compares the walked `T1` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `expected.epsilon0` | `expectedDims["epsilon0"]` | `.wl:257,260` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration, comparison oracle, and zero-expression fallback; `.wl:276,655-661` compares the walked `epsilon0` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `expected.epsilon1` | `expectedDims["epsilon1"]` | `.wl:257,260` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration, comparison oracle, and zero-expression fallback; `.wl:276,655-661` compares the walked `epsilon1` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `expected.P0Physical` | `expectedDims["P0Physical"]` | `.wl:257,261` | expected target | `.wl:272,276,655-661` | PROPOSE NOT EMIT | Hand-typed declaration and comparison oracle; `.wl:276,655-661` compares the walked `P0Physical` against it, so emitting both would duplicate one quantity under target/result names. A symmetric two-engine omission remains invisible to the comparator. |
| `baseline.A0` | `baselineDimAudit["Computed"]["A0"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.A0` | Record denotes the baseline dimension computed by `dimOf` for the live `A0lead` expression. B1 flag: NOT FLAGGED; no current `A0` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §1.3/§1.7 describe a pathA_29 residual continuation but name no future dimension record with candidate key `A0`. |
| `baseline.A1` | `baselineDimAudit["Computed"]["A1"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.A1` | Record denotes the baseline dimension computed by `dimOf` for the live `A1lead` expression. B1 flag: NOT FLAGGED; no current `A1` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and §1.3/§1.7 describe a pathA_29 residual continuation but name no future dimension record with candidate key `A1`. |
| `baseline.T0` | `baselineDimAudit["Computed"]["T0"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.T0` | Record denotes the baseline dimension computed by `dimOf` for the live `T0dc` expression. B1 flag: NOT FLAGGED; no current `T0` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and the stage/§1.7 names no future dimension record with candidate key `T0`. |
| `baseline.T1` | `baselineDimAudit["Computed"]["T1"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.T1` | Record denotes the baseline dimension computed by `dimOf` for the live `T1dc` expression. B1 flag: NOT FLAGGED; no current `T1` key exists in `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:44-167` (**the whole `## Quantities` table**), and the stage/§1.7 names no future dimension record with candidate key `T1`. |
| `baseline.epsilon0` | `baselineDimAudit["Computed"]["epsilon0"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.epsilon0` | Record denotes the baseline dimension computed by `dimOf` for the live `eps0` expression. B1 flag: FORWARD GROUP HAZARD; stage009 is not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and §1.7(2) identifies this record with stage009's same-named `epsilon0` transmission; a future stage009 `epsilon0` record will meet key `epsilon0`. The physics leg in §1.7 governs the identity decision. |
| `baseline.epsilon1` | `baselineDimAudit["Computed"]["epsilon1"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.epsilon1` | Record denotes the baseline dimension computed by `dimOf` for the live `eps1` expression. B1 flag: FORWARD GROUP HAZARD; stage009 is not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and §1.7(2) identifies this record with stage009's same-named `epsilon1` transmission; a future stage009 `epsilon1` record will meet key `epsilon1`. The physics leg in §1.7 governs the identity decision. |
| `baseline.P0Physical` | `baselineDimAudit["Computed"]["P0Physical"]` | `.wl:264-278,286-288` | clean walked result | `.wl:276,655-661` | PROPOSE EMIT as `computed_dims.P0_physical` | Record denotes the baseline walked dimension of Wolfram `P0Physical`; Python's record string key is `P0_physical`, proposed as join spelling. B1 flag: FORWARD GROUP HAZARD; stages 021 and 027 are not represented (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md:13-16`), and §1.7(2) identifies their `P0_physical` with this snake-case join key; future stage021/027 `P0_physical` records normalize to `P0Physical`. The physics leg in §1.7 governs the identity decision. |
| `corruptSourcedDims` | `corruptSourcedDims` | `.wl:289` | mutation-only declaration map | `.wl:290-292,444,461-463` | PROPOSE NOT EMIT | Copy of `baseDims` with the `M0` declaration deliberately replaced for the sourced-dimension failure probe. |
| `corruptFreeDims` | `corruptFreeDims` | `.wl:293` | mutation-only declaration map | `.wl:294-296` | PROPOSE NOT EMIT | Copy of `baseDims` with the unconsumed `qfree` declaration deliberately replaced for the free-carrier control. |
| `corruptSourcedComputed` | `corruptSourcedDimAudit["Computed"]` | `.wl:264-278,290-292` | mutation-only walked result set | `.wl:276,739,745-750` | PROPOSE NOT EMIT | Seven-result set produced under the deliberate `M0` declaration corruption; only mutation evidence, never a clean artifact value. |
| `corruptFreeComputed` | `corruptFreeDimAudit["Computed"]` | `.wl:264-278,294-296` | mutation-only walked result set | `.wl:276,739,745-750` | PROPOSE NOT EMIT | Seven-result set produced under the deliberate `qfree` declaration corruption; only control evidence, never a clean artifact value. |
| `caseComputedFamily` | `caseFor[mode]["Dimension"]["Computed"]` | `.wl:424-494` (construction `.wl:461-463`, returned `.wl:490-492`) | memoized conditional walked result family | `.wl:276,477` | PROPOSE NOT EMIT | Fourteen memoized states are reached: clean-equivalent `baseline`, `neutral`, `neutralized_independent_rerun`; control `selector`; and mutations `asserted_selector`, `wrong_sign`, `perfect`, `no_consistent`, `corrupt_v1`, `corrupt_order`, `corrupt_dimension`, `assert_not_derive`, `emit_epsilon`, `decouple`. They duplicate or deliberately alter the clean seven-result set and an in-body emission would mix states. |
<!-- ENUM_ROW_TABLE_END -->

**Axis order (§3.1).** `DETERMINED: (L,M,T)`. At the **enumerated pre-emission revision** the order was readable only at *prose
strength*: from the executed label `subheading["Exact (L,M,T) dimensional gate and free-carrier control"]` (`.wl:654`) — a print,
not a slot→label data structure — plus the non-zero slot coverage of `baseDims` (`a`, `omega` and `M0` at `.wl:248-249`, with
mixed/fractional entries through `.wl:253`), so no slot was evidenced only by zeros. That revision carried no mechanical binding
comparable to `{{"L",d[[1]]},...}`, and the enumeration said so honestly.
⭐ **The committed file now binds the order MECHANICALLY, so the evidence grade IMPROVES — that improvement is a result of this
conversion.** The step-(d) emitter declares exactly one slot→label structure, `dimensionAxisSlots = {{"L",1},{"M",2},{"T",3}}`
(`.wl:212`); the emitted `axes=` label is `StringRiffle[dimensionAxisSlots[[All,1]], ","]` (`dimensionAxesLabel[]`, `.wl:213`) and
every exponent vector is `(dimension[[#[[2]]]] &) /@ dimensionAxisSlots` (`dimensionComponents[]`, `.wl:214`), the two read
together by `printDimRecord` (`.wl:216-220`). Label and exponents therefore derive from that one structure, so permuting it moves
both together — which is what `manifests/DIMENSION_REWRITE.md` §4-(b) requires.

<!-- JOIN_KEY_TEXT_BEGIN -->
**B2 determination — join-key evidence and object kind.** Only names were read from the Python engine; no Python value, axis
order, or dimension was used. For sourced symbols, the join key is the **SymPy symbol name** written in constructor name strings
at `research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:163-177`, not the Python variable
identifier used as a dictionary key at the same file `:460-482`. For the computed `P0Physical` record, the join key is the Python
record string key `P0_physical` at `:493`. The applied rule is: retain an identical cross-engine spelling; otherwise use the
SymPy symbol name for a sourced symbol, or the Python record string key for a computed record. It reproduces
`cs→c_s`, `Keta→K_eta`, `TOmega→T_Omega`, `Z0ret/Z1ret→Z0_ret/Z1_ret`,
`OmegaU/OmegaW→Omega_U/Omega_W`, `Rmix→R_mix`, `gU/gW→g_U/g_W`,
`etaNull→eta_null`, `qfree→q_free`, and `P0Physical→P0_physical`. The Python variable identifiers `Keta`, `TOmega`,
`Z0ret`, `Z1ret`, `OmegaU`, `OmegaW`, `Rmix`, `gU`, and `gW` are therefore not the evidence for the proposed underscore
spellings. Names remain provisional join keys; the fresh physics leg must adjudicate them before emission.
<!-- JOIN_KEY_TEXT_END -->

**Reachability and provenance (§3.2, §3.4, §3.6).** All 29 proposed values survive to readable clean bindings:
the declarations at `baseDims[symbol]` (`.wl:247-255`) and the walked results at
`baselineDimAudit["Computed"][name]` (`.wl:264-288`). Thus `UNREACHABLE proposed records = {}`. At print time the first 22
hold declarations (18 literal vectors plus four distinct aliases to `zeroDim`); the last seven hold values computed by this
audit's `dimOf` walk. The fallback at `.wl:272` uses a target literal only when `TrueQ[expression == 0]`; no baseline symbolic
expression takes that branch, whereas mode `perfect` sets `A0`, `A1`, `epsilon0`, and `epsilon1` to zero at
`.wl:440,447-463`. Non-emitted objects are reachable too; they are excluded for helper/target/mutation/duplicate-state reasons,
not reachability. `NONE_FOUND` consumers are first-class for `sourced.R0`, `sourced.R1`, `sourced.etaNull`, and
`sourced.qfree`, with the finite full-file scopes in their rows.

**Invocation and one-clean-set context (§3.3).** `result = Catch[runAll[], ...]` calls `runAll` once at `.wl:802-806`;
`runAll` calls `runDimensionalGate[]` once at `.wl:785-798`, specifically `.wl:792`. That routine therefore provides a
single execution context for all 29 readable clean bindings. The three top-level `dimensionAudit` calls occur at
`.wl:286-295`. Memoization at `.wl:424` makes `caseFor` evaluate once for each of 14 reached modes: the direct
`baseline/selector/neutral` calls at `.wl:611-616`, the shared `neutralized_independent_rerun` at `.wl:496-514`, and ten
mutation modes named at `.wl:648-649,669,702-704,712-715`. Total `dimensionAudit` evaluations are therefore 3+14=17.
No A3b terminal condition fires: records are proposed, axes are determined, every proposed value is readable, and an
exactly-once clean emission context exists.

**Fractional values (§3.5).** Among proposed values, only `baseDims[gU]` and `baseDims[gW]` are non-integral, each carrying
exact rational entries `{-1/2,1/2,-2}` at `.wl:253`. No float serialization is permissible. The other proposed declaration
literals and all seven clean walked vectors are integral.

**Reconciliation with the provisional survey and prior stage note.** Agreement below is **not independent corroboration**:
both prior records were read before this comparison.

| crosswalk key | current source-derived item(s) | prior manifest §8 / stage note | verdict |
|---|---|---|---|
| `X.zero` | `zeroDim` | Manifest `.md:522` and stage-note dimensional gate §1.4 (`:96-104`) do not record it as an object. | `CURRENT_ONLY`; no contradiction. |
| `X.numeric` | `v0`, `v1`, and `coeff1` are uniformly ruled out; the walker visits only evaluated/fused numeric atoms, not stable exponent-vector bindings. | Stage note §1.3 (`:77-94`) records `v0/v1` physical use; prior records do not classify these scalar bindings as dimension objects. | `B4_REMEDIATED`; physical use agrees and no numeric binding is a row. |
| `X.sourced.literal` | `sourced.a`, `cs`, `omega`, `M0`, `D1`, `R0`, `R1`, `D0`, `K0c`, `Keta`, `TOmega`, `Z0ret`, `Z1ret`, `OmegaU`, `OmegaW`, `Rmix`, `gU`, `gW` | Manifest `.md:522` says 22 `baseDims` declarations and `.md:557-575` singles out fractional `gU/gW`; stage note §1.4 (`:96-104`) states the gate-driving subset. | `AGREEMENT` for every named source entry and the fractional pair; the stage note is not exhaustive. |
| `X.sourced.alias` | `sourced.etaNull`, `gain0`, `gain1`, `qfree` | Included only in manifest's aggregate 22; stage note §1.4 (`:102-104`) records the `q_free` corruption/control but not the other clean declarations as dimension objects. | `AGREEMENT` where stated; prior object-level coverage is incomplete. |
| `X.expected` | all seven `expected.*` rows | Manifest `.md:522` says seven declared `expectedDims`; stage note §1.4 (`:96-104`) prints the same target/result claims without distinguishing oracle from result. | `AGREEMENT` on presence; current classification adds the target/result distinction. |
| `X.baseline` | all seven `baseline.*` rows | Manifest `.md:522` says seven computed values via global `baselineDimAudit`; stage note §1.4 (`:96-104`) states all seven dimensions. | `AGREEMENT` for all seven values and reachability, not independent corroboration. |
| `X.mutation.maps` | `corruptSourcedDims`, `corruptFreeDims` | Manifest `.md:522` mentions 17 invocations and the zero-substitution hazard; stage note §1.4 (`:102-104`) mentions both corruptions, but neither enumerates the maps as objects. | `CURRENT_ONLY_AS_OBJECTS`; behavior agrees. |
| `X.mutation.results` | `corruptSourcedComputed`, `corruptFreeComputed` | Prior records report only their verdict behavior, not their returned seven-vector sets. | `CURRENT_ONLY`; no prior-only contrary object. |
| `X.case.family` | `caseComputedFamily` with 14 modes | Manifest `.md:522` correctly records 14 memoized modes and the `perfect` four-of-seven fallback; the stage note describes affected probes but does not enumerate the returned dimension set. | `PARTIAL_AGREEMENT`; invocation/hazard agree, object record was omitted. |
| `X.axes` | `(L,M,T)`, now **code-bound** in the committed file: `dimensionAxisSlots` (`.wl:212`) feeds both the `axes=` label and every exponent vector; prose-strength only at the enumerated pre-emission revision (`.wl:654`) | Manifest `.md:522` gives the same order and the same *pre-emission* evidence grade; stage note §1.4's heading (`:96`) names the same tuple machinery. | `AGREEMENT`, not independent corroboration. |
| `X.name` | `P0Physical` versus `P0_physical` | Manifest `.md:576-581` records the same name-pairing debt. | `AGREEMENT`; current proposal chooses `P0_physical` as join key. |
| `X.count` | 42 rows under the declared membership rule; 29 proposed records = 22 declarations + 7 clean results | Manifest census `.md:505` says `py dims=29`, while its reachability prose `.md:522` explicitly lists 7 computed + 22 `baseDims` + 7 `expectedDims` = 36 clean vector objects. | `DISAGREEMENT/SCOPING FINDING`: 29 matches the proposed artifact set, not the survey's own source-object list; source membership must not be reported as 29. |
| `X.prior_only` | none | Every manifest/stage-note dimension object or claim checked above has a source locus; no prior-only item failed to resolve. | `AGREEMENT: NONE_FOUND` (finite prior scope: manifest §8 and stage note §1.4 (`:96-104`), plus its `v0/v1` data flow in §1.3 (`:77-94`)). |

This record establishes keyed correspondence and makes omissions visible; it does not establish that any row is physically right
or that the predicates were exhaustive. Completeness would require a semantics-aware Wolfram analysis that proves every runtime
value carrying an exponent vector is classified across every reachable state, plus independent human review of the membership
definition. Because Wolfram evaluation can synthesize values dynamically and opaque helpers can hide them, absolute completeness
is not realistically provable for arbitrary source; reviewable completeness for this finite file is achievable but remains a
judgement, not what this table proves.

<!-- ENUM_COUNT|source_lines|769 -->
<!-- ENUM_COUNT|raw_hits|179 -->
<!-- ENUM_COUNT|p1_hits|23 -->
<!-- ENUM_COUNT|p2_hits|113 -->
<!-- ENUM_COUNT|p3_hits|25 -->
<!-- ENUM_COUNT|p4_hits|18 -->
<!-- ENUM_COUNT|promoted_hits|52 -->
<!-- ENUM_COUNT|ruled_out_hits|127 -->
<!-- ENUM_COUNT|row_count|42 -->
<!-- ENUM_COUNT|proposed_records|29 -->
<!-- ENUM_COUNT|sourced_proposed|22 -->
<!-- ENUM_COUNT|computed_proposed|7 -->
<!-- ENUM_COUNT|literal_sourced|18 -->
<!-- ENUM_COUNT|alias_sourced|4 -->
<!-- ENUM_COUNT|run_all_invocations|1 -->
<!-- ENUM_COUNT|run_dimensional_gate_invocations|1 -->
<!-- ENUM_COUNT|top_level_dimension_audits|3 -->
<!-- ENUM_COUNT|case_modes|14 -->
<!-- ENUM_COUNT|total_dimension_audits|17 -->
<!-- ENUM_COUNT|fractional_proposed|2 -->
<!-- ENUM_COUNT|b1_flagged|12 -->
<!-- ENUM_COUNT|b1_not_flagged|17 -->

### 1.7 Physics-leg verdict (step (c1)) — derived from the model, independently of the conversion

Run on a fresh leg against this stage's existing declarations, deriving each dimension from the **model** —
`docs/model_map.md` §2, **this stage's own physics, and `notes/parameter_register.md`**, which are the three
sources `manifests/DIMENSION_REWRITE.md:260-262` names — rather than checking a claim, and **blind to the
step-(a) enumeration** (which was being written in parallel; the leg's report was quarantined from that
session and vice versa). Where the two agree below, that is two parties, not one reviewing the other.
⚠ The model map alone would not carry this: `:49` gives the four primitives (`ħ`, `m_GNLS`, `K`, `ρ0`), `:55`
the brane mass density `ρ_br = M L⁻³`, `:62` the EOS speed `c_s² = 5Kρ⁴/m`, and nothing else in §2 reaches the
port kernel, the stiffness sector or the moments. Most of the 34 routes below run through this stage's own
homogeneity or through a register row; the table names which.

**(1) Per-quantity verdict — 34 objects, each with its route.** `manifests/DIMENSION_REWRITE.md:270-271`
requires a verdict **and** its derivation route *per quantity* (*"A verdict without its route is an
assertion"*), so the itemization is seated here, not left as a count. The 34 are the 22 sourced
declarations, the 7 expected/computed records, and **5 dimension-bearing intermediates neither engine
emits** — `P_port`, `Delta_port`, `N0_from_port`, `P0_raw`, `K1` — which the §8 survey's count of 29 does
not carry. ⚠ **All five sit in the same position, and the count is five, not two:** none of them is a key
of the `expressions` map the dimension checker walks (`.py:510-518` / `.wl:266-269`, whose seven keys are
exactly the seven emitted records). Each is reached only as a sub-expression of one of those seven, so its
triple is computed inside the walk, constrained solely by the assertion on that record's *total*, and then
discarded. What differs among the five is not their treatment here but whether any other stage declares a
**dimension triple** for the same object — being *named* outside this stage is **not** the discriminator (all
five are); see **U11**.

⛔ **Tally — it depends on the basis, and `D0`'s uncertainty propagates to `P0_raw` and `P0_physical`.**
- **As identifications against the corpus** (the banked reading, because W1 is open):
  **24 CORRECT / 0 WRONG / 10 UNDETERMINED**.
- **Scored inside this stage's own closure**, where `D0`'s triple is *forced* by the asserted
  `[P0_physical] = 1`: **27 CORRECT / 0 WRONG / 7 UNDETERMINED**.
Apart from rows 8, 29, and 33, every object's verdict is the same on both bases.

⚠ Dimensions are given **labelled** (`M L⁻³ T⁻²`), never as bare tuples: the `(L,M,T)` slot order was
prose-only at the **pre-emission `.wl` revision** (the `subheading["Exact (L,M,T) dimensional gate and
free-carrier control"]` print, now `.wl:654`), and is **code-bound in the committed `.wl`** by
`dimensionAxisSlots` (`.wl:212`), which feeds both the `axes=` label and every exponent vector (§1.6,
*Axis order (§3.1)*). Either way a transposition must not be able to hide inside this table.

| # | object | class | dimension | verdict | derivation route |
|---|---|---|---|---|---|
| 1 | `a` | decl | **L** | CORRECT | `a = ħ/(m_GNLS c_s0)` `CONV` pin (register `:132`) with `model_map:49`'s primitives ⇒ `(M L²T⁻¹)/(M · L T⁻¹) = L`. 2nd route: 016's measure `dV = a²dw dΩ` asserted `L³` (`stage016 sympy:335,353`). |
| 2 | `c_s` | decl | **L T⁻¹** | CORRECT (value); *which* speed UNDETERMINED — (2) | `c_s² = 5Kρ⁴/m` (`model_map:62`) with `[K] = M L¹⁸T⁻²`, `[ρ] = L⁻⁴`, `[m] = M` (`:49`) ⇒ `L²T⁻²`. 2nd route in-stage: `z = aω/c_s` is the DtN argument ⇒ `[c_s] = [a][ω]`. |
| 3 | `omega` | decl | **T⁻¹** | CORRECT | Time convention `e^{−iωt}` ⇒ `[ωt] = 1`. |
| 4 | `M0` | decl | **M T⁻¹** | UNDETERMINED | `M0 = ∫_brane S_leak d³x` (`stage008 sympy:315`). Mass reading `ρ_br = M L⁻³` (`model_map:55`) ⇒ `M T⁻¹` = this declaration; number reading `ρ_B = L⁻⁴` ⇒ `L⁻¹T⁻¹`; 009 declares `T⁻¹` (`sympy:467`). 008 **declines to dimension it** (`sympy:546`) — so **among the four stages that carry this moment (006, 008, 009, 023) none records a chosen convention**: 008 declines and the other three declare different triples. Scope of the negative: those four stages, not the corpus at large. |
| 5 | `D1` | decl | **M L T⁻¹** | UNDETERMINED | `D1_i = ∫x_i S_leak d³x + ∫j_i d³x`; same unchosen anchor as row 4. ⭐ The *ratio* `[D1]/[M0] = L` — the multipole ladder — is CORRECT under all three readings. |
| 6 | `R0` | decl | **M T⁻¹** | UNDETERMINED | `R0 = −M0` (008 target) ⇒ `[R0] = [M0]`; inherits row 4's anchor and is read by no expression (U5). |
| 7 | `R1` | decl | **M L T⁻¹** | UNDETERMINED | `R1 = −D1` ⇒ `[R1] = [D1]`; same, unread (U5). |
| 8 | `D0` | decl | **M L⁻¹ T⁻²** | ⚠ SPLIT — CORRECT in this stage's closure, UNDETERMINED as the 017 identification | Forced in-stage by `[P0_physical] = 1` (`.py:493` / `.wl:261`) with `[N0_from_port] = M L⁻¹`; but 017's exported D-lane gives `M T⁻²`. See the ⭐ note below and **W1**. |
| 9 | `K0c` | decl | **M T⁻²** | CORRECT | ℓ=0 collective stiffness. 013's `K_AB` entries and `k_shared` are `M T⁻²` with `M_AB` at `M` (`CANONICAL_DIMENSIONS.md:98-101,102-105`) ⇒ `K/M = T⁻² = ω²`. ⚠ this stage's gate cannot test it (U3). |
| 10 | `K_eta` | decl | **M T⁻²** | CORRECT | The `L⁰` **reduced radial scalar** `K̃ = ∫[T_w β'² + K_η β²] dV`, registered `M T⁻²` (register `:184`) — not 013's line density `M L⁻¹T⁻²` (`:179`) nor 016's volume density. |
| 11 | `T_Omega` | decl | **M T⁻²** | CORRECT | `T̃_Ω = ∫T_Ω β² dV`, registered `M T⁻²` (register `:184`) — not 017's density `M L⁻³T⁻²` (`:182`). Identified by `K1 = K_eta + 2·T_Omega` being 016's `K̃ + λ_m·T̃_Ω` at `λ_m = ℓ(ℓ+1) = 2`, i.e. ℓ=1. |
| 12 | `Z0_ret` | decl | **M T⁻²** | CORRECT | `ε0 = Z0_ret/K0c` is 009's **dimensionless** DC transmission, aliased at register `:169` ⇒ `[Z0_ret] = [K0c]`. |
| 13 | `Z1_ret` | decl | **M T⁻²** | CORRECT | `ε1 = Z1_ret/K1` ⇒ `[Z1_ret] = [K1] = [K_eta]`. |
| 14 | `Omega_U` | decl | **T⁻¹** | CORRECT | The 2×2 port diagonal is `Ω_U²`, `Ω_W²`; 024's analogue diagonal `ϖ_q2`, `ϖ_Φ2` is `T⁻²` (register `:186`) ⇒ `[Ω] = T⁻¹`. |
| 15 | `Omega_W` | decl | **T⁻¹** | CORRECT | Same slot, same route. |
| 16 | `R_mix` | decl | **T⁻²** | CORRECT | `Δ_port = Ω_U²Ω_W² − R_mix²` homogeneity ⇒ `[R_mix] = [Ω_U] + [Ω_W] = T⁻²`; 024's off-diagonal `λ_c` is `T⁻²` (register `:187`). |
| 17 | `g_U` | decl | **M¹ᐟ² L⁻¹ᐟ² T⁻²** | CORRECT (two routes) | In-stage: `[P_port]² = [N0_from_port]·[Δ_port]² = M L⁻¹T⁻⁸` ⇒ `[P_port] = M¹ᐟ²L⁻¹ᐟ²T⁻⁴`, less `2[Ω_U]`. ⭐ Independently: 024's `g_base = √ρ_eff·c_s²·I25·Ξ_Q/a^{7/2}` (register `:188`) gives the identical value from an unrelated formula. |
| 18 | `g_W` | decl | **M¹ᐟ² L⁻¹ᐟ² T⁻²** | CORRECT | `P_port = Ω_U²g_W + R_mix g_U` is an Add ⇒ `[g_W] = [g_U]` (given `2[Ω_U] = [R_mix]`). |
| 19 | `eta_null` | decl | **1** | CORRECT | Injected as `z0 → z0 + η·K0c` (`.py:299`) — an Add ⇒ `[η] = 1`. ⚠ Correct **and effectively unchecked**: the `.py` computes the injected audit but asserts only its nullity (`.py:837,851`), never its dimensions, and the `.wl`'s inject path (`.wl:189-194`) calls no `dimensionAudit` at all — the three top-level calls are at `.wl:286-295` and the per-mode call at `.wl:461-463`, and `inject` is not a `caseFor` mode. Scope: these two engines. |
| 20 | `gain0` | decl | **1** | CORRECT | `T0 → gain0·T0` with `[T0] = 1` ⇒ `[gain0] = 1`. ⚠ Unchecked in practice: `decoupled` precedes `dimensional` in the ladder (`.py:675,677` / `.wl:402,404`), so in `decouple_knobs` — the only one of the 14 `caseFor` modes in which `gain0` enters an expression at all (`.py:306-308`) — the dimensional rung is never reached. |
| 21 | `gain1` | decl | **1** | CORRECT | Same. |
| 22 | `q_free` | decl | **1** | UNDETERMINED *by construction* | It enters **no** checked expression, and both engines assert exactly that (`.py:917` / `.wl:668`). Every triple is therefore equally admissible; `1` records its role, not a derivation. |
| 23 | `A0` | emit | **M T⁻¹** | UNDETERMINED (correct *as a computation*) | `A0 = i·v₀·(aω/c_s)·M0·(1−T0)` (`.py:342`); `(aω/c_s)` and `(1−T0)` are dimensionless ⇒ `[A0] ≡ [M0]` identically. Inherits row 4's unchosen anchor. ⚠ its target literal *is* `SOURCED_DIMS[M0]` — U2. |
| 24 | `A1` | emit | **M L T⁻¹** | UNDETERMINED (correct *as a computation*) | `A1 = i·v₁·(aω/c_s)³·D1·(1−T1)` (`.py:343`) ⇒ `[A1] ≡ [D1]`. `[A1]/[A0] = L` survives any normalization choice. |
| 25 | `T0` | emit | **1** | CORRECT | `K0c/(K0c+Z0_ret)`: the Add forces `[Z0_ret] = [K0c]`, the quotient is then dimensionless. |
| 26 | `T1` | emit | **1** | CORRECT | `K1/(K1+Z1_ret)` with `K1 = K_eta + 2·T_Omega`: inner Add forces `[T_Omega] = [K_eta]`, outer forces `[Z1_ret] = [K_eta]`. |
| 27 | `epsilon0` | emit | **1** | CORRECT | Its **own** quotient `Z0_ret/K0c` (`.py:304`) — not a consequence of the `T` formula. Matches 009's dimensionless `ε0`. |
| 28 | `epsilon1` | emit | **1** | CORRECT | `Z1_ret/K1` (`.py:305`). |
| 29 | `P0_physical` | emit | **1** in-stage; **L⁻¹** on the corpus `D0` reading | ⚠ SPLIT — CORRECT in this stage's closure, UNDETERMINED as a corpus identification; the corpus reading contradicts this stage's dimensionless target | `(c_s/a)²·P0_raw`: stage-local `[D0]=M L⁻¹T⁻²` gives `T⁻²·T²=1`, while corpus `[D0]=M T⁻²` gives `T⁻²·L⁻¹T²=L⁻¹` (`.py:225,231` / `.wl:110,116`). |
| 30 | `P_port` | interm | **M¹ᐟ² L⁻¹ᐟ² T⁻⁴** | CORRECT | `Ω_U²g_W + R_mix g_U` (`.py:222` / `.wl:107`). |
| 31 | `Delta_port` | interm | **T⁻⁴** | CORRECT | `Ω_U²Ω_W² − R_mix²` (`.py:223` / `.wl:108`). |
| 32 | `N0_from_port` | interm | **M L⁻¹** | CORRECT | `[P_port]²/[Δ_port]²` — ⛔ **both squared** (`.py:224` / `.wl:109`) — `= M L⁻¹T⁻⁸ / T⁻⁸`. Matches 021's `SOURCED_N0_DIM` (`sympy:145`) and 027's registered `[N0_den] = M L⁻¹` (register `:191`), reached here from a different expression. |
| 33 | `P0_raw` | interm | **T²** in-stage; **L⁻¹T²** on the corpus `D0` reading | ⚠ SPLIT — CORRECT in this stage's closure, UNDETERMINED as a corpus identification | `[N0_from_port]/[D0]`: `M L⁻¹/(M L⁻¹T⁻²)=T²` in-stage, but `M L⁻¹/(M T⁻²)=L⁻¹T²` on the corpus reading (`.py:225` / `.wl:110`). |
| 34 | `K1` | interm | **M T⁻²** | CORRECT | `K_eta + 2·T_Omega` (`.py:285` / `.wl:122`) — 016's ℓ-stiffness assembly at ℓ=1. |

⚠ **Six of the ten corpus-basis `UNDETERMINED` share one root cause** — `M0`, `D1`, `R0`, `R1` and hence `A0`, `A1`
are undetermined *not* because this stage is inconsistent but because **none of the four stages that carry
the moment — 006, 008, 009, 023 — chose** the pathA_28 moment convention (number density `ρ_B` `L⁻⁴` vs mass
density `ρ_br` `M L⁻³`; `d³x` vs `d⁴x`); 008 declines to dimension it at `sympy:546` and the other three
declare different triples. That set of four is the scope of the negative.
Both endpoints are internally coherent. **This is a model decision, not a script defect** — recording it
**resolves the convention** that 006, 008, 009 and 023 each read differently; ⚠ it does not by itself edit
those four loci, which still have to be brought into line with whichever reading is recorded. `q_free` is
undetermined by construction: an unread free-carrier control. The final three are `D0` and its propagated
dependents `P0_raw`/`P0_physical`; inside the stage-local closure `D0` is forced and those three leave the
undetermined set, giving the separate 27/0/7 tally above.

⭐ **`D0` — the two readings are not documented normalizations; propagated, they expose a contradiction.**
- **Within this stage's own closure the triple is FORCED.** `P0_physical = (c_s/a)²·N0_from_port/D0` with
  the asserted target `[P0_physical] = 1` (`.py:493` / `.wl:261`) and `[N0_from_port] = M L⁻¹` (row 32)
  gives `[D0] = (c_s/a)²·[N0_from_port] = T⁻²·M L⁻¹ = M L⁻¹T⁻²` — exactly the declaration at `.py:468` /
  `.wl:250`. On that basis: **CORRECT**. ⚠ The anchor is this stage's *own* assertion, so per **U1**/**U2**
  it is forced, not externally pinned.
- **As the identification with stage017's exported D-lane it is `UNDETERMINED`.** 017 builds
  `K₂ = c_self·(K̃ + λ_m·T̃_Ω)` and `D0 = K₂ − (B̃0+Z̃0)`
  (`notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:30-31`) and exports that D-lane to 022/023
  (`:113-115`); with `[K̃] = [T̃_Ω] = M T⁻²` (register `:184`) that is `[D0] = M T⁻²`, which register `:185`
  records. Stages 021 (`sympy:146`), 023 (`sympy:468`) and 027 (`sympy:206`) all declare `M L⁻¹T⁻²` for
  what this note identifies as the same object. Propagating the register reading gives
  `[P0_raw]=L⁻¹T²` and `[P0_physical]=L⁻¹`, in direct conflict with this stage's asserted
  `[P0_physical]=1`. Because **no 017→019 export normalization is recorded in the three places one would
  be** — 017's own export bullet (`notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:113-115`, which
  hands `K₂`, `B̃/Z̃` and the D-lanes onward with no rescaling), the register's `D0` row
  (`notes/parameter_register.md:185`), and 019's audit script (which per
  `manifests/DIMENSION_REWRITE.md:100` has no dimension machinery at all) — the two readings cannot
  presently be defended as alternative normalizations of one object. Which assertion is wrong remains
  undetermined; that does not weaken the contradiction. This is **W1** and **U8**.
- **Settled by:** an explicit 017→019 normalization carrying the missing `L⁻¹`; or a correction that makes
  the D-lane and stage023 `D0` distinct objects; or a correction to one of `[D0]`, `[N0_from_port]`, or the
  dimensionless `[P0_physical]` target. A bare preference between the existing triples settles nothing.

**(2) D4 name determinations.** ⭐ The build proposed its join keys independently, from the two engines'
own identifiers; the two sets agree, including the snake_case direction that *creates* the cross-stage
grouping.
- ✅ `a` **same** as 016/018 (throat/pin radius; 016's measure `a²dwdΩ`, 018/023's `z = aω/c_s`).
  ✅ `omega` **same** as 011/012. ✅ `D0` **same** as 021/027 (⚠ but W1 below). ✅ `epsilon0/1` **same** as
  stage009's same-named DC transmission symbols (`stage009 sympy:113,136-137`), which stage009 feeds into
  the separately derived `alpha(ε)=1/(1+ε)` and `neighbor_fraction(ε)=ε/(1+ε)` (`:190-195,230-235`);
  register `:169` states that stage023's `Z_ret/K` coordinates alias those transmissions. They are therefore
  the two forward same-key hazards §1.6 flags **under the `computed_dims.epsilon0`/`epsilon1` rows** — not
  the only forward hazards it flags (nine of its twelve `B1` flags are forward, and the `M0`/`D1` pair names
  stage009 too) — and not identities with either derived fraction.
  ✅ `P0_raw`/`P0_physical` **same** as 021/027 — and the `.wl` `P0Physical` / `.py` `P0_physical` split is
  aligned to snake_case here.
- ✅ `K_eta`, `T_Omega` vs the **013/016 densities**: a **different reduction level of the same object**
  (`M L⁻¹T⁻²` line, `M L⁻³T⁻²` volume, `M T⁻²` reduced — each off by exactly its measure's L-power).
  Keep the names; the `NEEDS_ADJUDICATION` the canonical table shows is **correct output** and is the
  counted reduction debt. ⛔ Per §7, do not rename the variants apart. *(Contested by the parallel
  enumeration; adjudicated in **(2b)** — the names stand.)*
- ⛔⛔ **THE MERGE TO REFUSE — `K_eta`/`T_Omega` (023) vs `K_tilde`/`T_Omega_tilde` (016).** Both carry
  `M T⁻²`; both are the same construction at the same reduction level. But 016's are the **ℓ=2** radial
  reductions against the profile `β₂`, while 023's appear in `K1 = K̃ + λ_m T̃_Ω` at `λ_m = ℓ(ℓ+1) = 2`,
  i.e. **ℓ=1**. Whether they are the same number turns entirely on whether the frozen wall's radial
  profile is ℓ-independent (`β₁ ≡ β₂`), and **no document states that in the trees searched** —
  `research/pde_ledger_v2/notes/` (incl. `stages/`), `research/pde_ledger_v2/paper/stages/`,
  `research/pde_ledger_v2/manifests/` and `docs/model_map.md`, in which the symbol `β₁` does not occur at
  all outside stage029's unrelated `β₁PN`. Trees outside that list are not covered by this negative.
  ⇒ **UNDETERMINED same-or-different; DO NOT MERGE.** This is the `c_s0`/`c_S` shape, and it is *more*
  attractive here because D4's ✅ direction — and, now, **(2b)**'s adjudication — appears to sanction it — a D4 pass "closing a grouping
  limitation" would merge the ℓ=1 and ℓ=2 reduced stiffnesses where no dimensional, cross-engine or
  re-run check could ever see it. Settled by a statement or derivation of `β₁ ≡ β₂`.
  ⚠ **Inverted convention, specific to `T_Omega`:** 016 uses the bare name for the *density* and the
  tilde name for the *reduced scalar*; 023 uses the bare name for the *reduced scalar*. A reader can
  "fix" the table in either direction and one of those directions is the bad merge.
  ⭐ This is the same hazard 016's own leg recorded from the other side (`ledger_stage016…md` §1.6 (2)),
  reached again independently here and extended from `T_Ω` to `K_η`.
- ⛔ `c_s` — **different** from 012's frozen-wall `c_S` (023's own note §4 says so); **UNDETERMINED**
  against 018's bulk asymptotic `c_s0`, because **nothing in the four places that define this carrier
  states which `ρ` it is evaluated at** — §4 of this note (which gives only `c_s²=5Kρ⁴/m` and the role
  "units carrier in `z = aω/c_s`"), the two engines (which declare `[c_s]` and never a state point),
  `docs/model_map.md:62` (formula plus `d ln c_s/d ln ρ = 2`, no evaluation point), and the register. Scope:
  those four; an evaluation point stated somewhere unsearched would settle it.
  The candidate keys happen not to collide; that is luck, not a guard. *(Also flagged by the
  enumeration; adjudicated in **(2b)** — forward hazard, not a live merge; the name stands.)*
- ⛔ `M0`/`D1` — one intended object flowing 008 → 009/010 → 022/023 with **three** declarations:
  `L⁻¹T⁻¹` (006+008 chain), `T⁻¹` (009 `MOMENT0`), `M T⁻¹` (023). Do **not** rename apart; give it one
  name plus an explicit measure statement. ⚠ 009 spells it `MOMENT0`, so **the table generator's
  `candidate_key` normalization** (`scripts/generate_canonical_dimension_table.py:116-128` — splits on
  separator boundaries, never case-folds) sends `MOMENT0`→`MOMENT0` and `M0`→`M0` and will **never**
  group them; the conflict stays invisible even after both stages convert — renaming `MOMENT0` → `M0`
  *creates* the detection. The "never" is scoped to that generator, not to human review.
- ⛔ `K0c` must **not** merge with 013's `k_shared` (which is the shared dimension tag of the `K_AB`
  entries, not a physical stiffness) — the `m_shared` landmine shape. Which projection of `K_AB` gives
  `K0c` is exactly the `PENDING` scalar-reduction debt at register `:170`.
- `Z0_ret`/`Z1_ret` are return admittances (`M T⁻²`), **not** 009's return amplitude `Z` (`[Z] = [M0]`);
  they differ dimensionally but will never be compared **by the table generator**, whose `candidate_key`
  (`scripts/generate_canonical_dimension_table.py:116-128`) sends `Z0_ret`→`Z0Ret` and `Z`→`Z`, so no group
  forms. Nothing here says a human reader could not compare them.

**(2b) ⚖ ADJUDICATION of a flagged conflict — `K_eta` / `T_Omega` / `c_s`.** *Not a new finding.* This
settles a conflict raised against (2) by the parallel step-(a) enumeration
(`_scratch/stage023_orch/ENUM_REVIEW.md:13-56`, `HIGH`), which read the three proposed names as
"known-wrong merges" — apparently untenable against this stage's own text — because emitting them makes
the canonical table GROUP this stage's `M T⁻²` scalars with the same-named 013/016 densities, while §4 of
this note says they are "**NOT identified with the raw 013/017 densities**" (`:797-804`).
⇒ **The names stand.** `K_eta` and `T_Omega` are proposed as-is, and the `NEEDS_ADJUDICATION` group that
results is **correct output, not a defect.** The two statements the flag sets against each other are both
true; they are not about the same thing.
- **The manifest already assigns this stage the third slot in that group, deliberately.**
  `manifests/DIMENSION_REWRITE.md:477-480`: *"⭐ Most 'conflicts' are REDUCTION LEVELS, not drift. The
  model spans a 4D bulk and a 3D brane, so integrating out a coordinate shifts L-powers legitimately:
  `K_eta` is volume-density (L³, 016) vs line-density (L¹, 013) vs reduced scalar (L⁰, 023), all
  integrating to `M T⁻²`. Verified across four quantities (`K_eta`, `T_w`, `μ_η`, `T_Omega`), each off by
  exactly its measure's L-power."* Then the rule, `:481-483`: *"⛔ **Never resolve one by renaming the
  variants apart** — it destroys the check **and** hides reduction debt … The dimension is the symptom;
  the unperformed reduction is the thing (register `:170` counts it as debt)."* The reduced-scalar slot
  the manifest hands 023 is exactly the `M T⁻²` these two carry, under exactly these names.
- **A group is a review flag, not a merge.** `CANONICAL_DIMENSIONS.md:21-23`: *"Exact emitted names are
  primary; candidate groups below are review flags, never automatic merges."* And `:178-179`: *"A
  differing group is flagged `NEEDS_ADJUDICATION`; the generator never chooses a winner."* Grouping
  therefore asserts **nothing** about identity — only that two names normalize alike and that a human owes
  a determination.
- ⇒ **The collision was never real.** §4's "not identified" (`:797-804`) says that **no earned equality
  exists**: the profile+measure scalar-reduction has not been performed, so the quantities are classed
  `FREE-UNREDUCED`, `PENDING`, and **counted as reduction debt** (`notes/parameter_register.md:170`). The
  group is precisely the mechanism that keeps that debt visible in the generated table; renaming the
  variants apart would delete **the one place a reader of `CANONICAL_DIMENSIONS.md` meets it** — the debt is
  carried independently at `notes/parameter_register.md:170`, so this is not a claim that renaming would
  erase it from the corpus. ⭐ The enumeration's "apparently untenable"
  reading was over-strict in **at least one respect — it treated a group as an identity claim** ("known-wrong
  merges", `ENUM_REVIEW.md:13`). ⚠ **`UNDETERMINED` whether that is its *only* over-strict respect:** B1 also
  objects that the rows carried "a join-spelling argument only … with no same-quantity / different-quantity
  determination" (`:28-30`), citing row text at `:203-204` (B1's own pre-fold numbering for the `sourced.Keta` / `sourced.TOmega` rows, seated here at `:235-236`) that has since been amended to defer explicitly to
  §1.7(2)/(2b). Whether that second objection was over-strict or simply answered by a later edit cannot be
  settled from the seated note, so this is stated as *one* respect, not *the* respect. Its own
  rows already deferred correctly: the `K_eta` and `T_Omega` rows above (the `sourced.Keta` / `sourced.TOmega` rows, `:235-236`) each close *"This
  records the proposal, group, and evidence only; §1.7(2)/(2b) owns the determination."*, and §1.6's B1
  determination paragraph (`:205-218`) says the section *"makes no tenability ruling"* (`:213`).
- ⛔ **Scope — this settles which NAME to emit, and nothing else.** Whether 023's reduced `K_η`/`T_Ω` are
  numerically the reduced 013/016 objects is **not** settled here: that is the PENDING reduction debt and
  it stays open. Per §4 (`:797-804`), an explicit profile+measure scalar-reduction to the calibrated wall
  packet is what would earn `DERIVED` and discharge it.
- ⛔⛔ **The refusal in (2) is UNCHANGED, and it points the other way.** Nothing here licenses merging
  023's `K_eta`/`T_Omega` with 016's `K_tilde` (`CANONICAL_DIMENSIONS.md:122`) / `T_Omega_tilde` (`:125`).
  That is a **different proposal**: same `M T⁻²`, same construction, **same reduction level**, differing
  only ℓ=1 vs ℓ=2 — so no measure's L-power separates them and **none of the three mechanisms this
  workstream owns — the dimensional gate, the cross-engine comparator, or a re-run — could ever catch it.**
  (Scope: those three. A physics derivation of `β₁` vs `β₂` still could.)
  A differing-dimension trio grouped in the table is a *visible flag*; a same-dimension pair merged is a
  *silent identity*. ⚠ Do not read the two cases as one determination — they run in opposite directions.
- **`c_s` — a forward hazard, not a live merge; the name stands.** No group forms today: 023's `c_s`
  normalizes to `cS`, and the only **medium-speed** candidate keys the table currently carries are `csDim`
  (012 `clean_walk.cs_dim`, `CANONICAL_DIMENSIONS.md:80`; also `corrupt_walk.cs_dim` at `:87`, same key) and
  `cS0Dim` (018 `c_s0_dim`, `:134`) — the only other `L T⁻¹` row is stage004's generic `velocity` tag
  (`:64`), which names no speed — and **no `cS` member exists**. The distinctness §4 records (`:792-793`, *"Distinct from the frozen-wall Helmholtz speed
  `c_S`, 011–017"*) stays on the record as the hazard for whenever 011–017's `c_S` is emitted: the group
  would then form as `cS` and require this same adjudication on the same terms. ⛔ With the caveat that
  per `manifests/DIMENSION_REWRITE.md:331-339` this is the *measured* 018 `c_s0`/`c_S` shape, where the
  comparator is blind by construction and two independent parties already made the wrong merge.
⚠ **Locus note:** ENUM_REVIEW B1 cites these §4 passages as `:366-369` and `:359-360`; in the seated note
they are `:797-804` (the `K0c`/`K1` "NOT identified with the raw 013/017 densities" bullet) and `:792-793`
(the `c_s` "Distinct from the frozen-wall Helmholtz speed `c_S`" parenthetical) — B1's numbering predates
the fold, and `:359-360` / `:366-369` now land inside §1.7. ⚠ These loci move whenever §1.6/§1.7 grow; an
earlier revision of this note cited them as `:693-700` and `:689-690`, which had already gone stale. Cite
them by quoted text when in doubt. The enumeration's own rows cite the seated loci correctly.

**(3) Coverage — ⛔ 29 of 29 dimension values are hand-typed literals in BOTH engines: 22 declarations
plus 7 targets typed on both sides, and the 7 computed records are walks over exactly those literals.**
⚠ Read that precisely — the 7 computed records *are* computed, by a live `dimOf` walk; what the finding
says is that the walk consumes no dimensional input from outside the two hand-typed tables, so nothing
in the block is pinned by anything but a literal. All 22 sourced
declarations and all 7 expected targets are literals on both sides; the 7 computed records are derived
from those literals alone. For scale: 016 was 12 of 21 (`ledger_stage016_l2_so3_covariance.md:194`), the 037
spike 8 of 21 and 018 6 of 10 (`manifests/DIMENSION_REWRITE.md:276`, **§4-(c1) item 3, *"stage037: 8 of 21. stage018: 6 of 10."***) — **this is the least independent
dimension block among the four for which a literals-in-both-engines count has been recorded at all: 016,
018, 037 and this stage.** The other converted stages (004, 011, 012, 013) carry no such count, so they are
outside the comparison rather than beaten by it.
- **What the gate does buy.** It is not a token comparison: the checker walks each of the seven
  expressions and compares the walked dimension against the hand-typed target (`.py:510-528` /
  `.wl:266-276`), so a clean result **does** impose algebraic consistency. Concretely it forces the eight
  independent vector relations of **U1**, and it is able to fire — **W3** shows that dropping the `g`
  dimensions moves `[P0_physical]` off `(0,0,0)`, and corrupting `[M0]` reaches `FAIL_DIMENSIONAL` (§2).
- ⛔ **What it does not buy — and this is the coverage finding.** Those eight relations leave a
  **24-parameter family** of declaration sets that pass everything (U1), five of the seven records evaluate
  to `(0,0,0)`, and the other two are the *same literal triples* as `[M0]` and `[D1]` (U2). So the gate
  constrains *relations among* declarations and pins almost none of them individually: **cross-engine
  agreement between two hand-typed tables is most of what a green result buys**, and a wrong declaration
  shared by both engines passes. Under an L↔M relabelling only 1 of the 7 computed records moves — the
  transposition detection this stage gets comes almost entirely from the 22 declarations.

**(4) ⭐⭐ Structurally uncheckable — the deliverable.** No conversion, cross-engine comparison or re-run
can close these:
- **U1** — the gate determines only **8 of its 16 live declarations**; `{[a],[ω],[c_s],[K0c],[K_η],[Ω_U],[Ω_W],[D0]}`
  are free generators, i.e. a **24-parameter family** of declarations passes everything (018's finding, 4× larger).
- **U2** — `EXPECTED_DIMS["A0"]` **is** `SOURCED_DIMS[M0]` and `["A1"]` is `[D1]`: the declaration is
  asserted against its own copy, so a shared wrong pair passes.
- **U3** — the register's load-bearing `K0c/K1 ≠ 013/017` decision rests on a declaration this stage
  provably cannot test: only *differences* are pinned.
- **U4** — the `expr == 0 → substitute the expected triple` shortcut (`.py:525` / `.wl:272`) is
  **identical in both engines**, so cross-engine comparison is blind to it; it is load-bearing in
  `perfect` mode, where 4 of 7 rows then self-match.
- **U5** — the `q_free` control and its companion assertion are one fact stated twice; `R0`/`R1` sit in
  the identical unread slot without being labelled controls.
- **U6** — the dimensions never touch the **rank/nullity headline**, in either direction. `rank0`,
  `native_nullspace_dimension` and `return_moving_nullity` are computed from the symbolic Jacobian with no
  reference to any dimension (`.py:378-426` / `.wl:126-148`), and matrix rank is invariant under the
  diagonal row/column rescalings a change of units induces — so `8`, `2`,
  `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` is **provably unit-robust**, and equally **no dimensional work can
  validate it and no dimensional error can falsify it**. The `dimensional` rung is bolted on beside the
  physics, not inside it. ⚠ **Not** the wider claim that dimensions never touch the verdict at all: a
  dimensional failure is a rung of the verdict ladder in both engines (`.py:677`, fed at `:706`; `.wl:404`,
  fed at `:477`), which is exactly why corrupting `[M0]` reaches `FAIL_DIMENSIONAL` (§2). The Jacobian rows
  are themselves dimensionally inhomogeneous — within row 1, `∂P0_raw/∂Ω_U` and `∂P0_raw/∂D0` differ — so
  the null-basis vectors carry no consistent dimension; only the integer rank/nullity survives, and nothing
  in either engine says so.
- **U7** — the three inconsistent `[M0]` values, with every catching mechanism disabled (008 emits none
  by policy; 009's assertion cannot fail; 023's is U2; the names do not group).
- **U8** — the `[D0]` contradiction: 017 implies `M T⁻²`, while 021/023/027 declare `M L⁻¹T⁻²`; under
  this note's same-object identification the former propagates to `[P0_physical]=L⁻¹`, contradicting the
  asserted dimensionless target. No export normalization is recorded (loci in the `D0` note above), and 017
  is one of the 13 audit scripts with **no dimension machinery** (`manifests/DIMENSION_REWRITE.md:100`), so
  the conflict is unreachable by every mechanism **in the dimension-rewrite workstream** — conversion,
  cross-engine comparison, the comparator, and re-run. It remains reachable by hand, as this leg shows.
- **U9** — the `neutered`/`neutralized` controls are label-level only (`.wl:616` is `X === X`).
- **U10** — `v₀ = 1` is consumed but never ablated; **of the two cited stage022 coefficients, only `v₁`
  carries a mutation flag** — `.py:337` binds `v0 = sp.Integer(1)` unconditionally while `:198,338` gate
  `v1` behind `corrupt_v1`, and `:891-895` is the only coefficient ablation either engine runs. Scope: this
  stage's mutation set, not the corpus's treatment of `v₀` elsewhere.
- **U11** — `N0_from_port` (`M L⁻¹`) and `P0_raw` (`T²` stage-local, `L⁻¹T²` on the corpus `D0` reading)
  are never **separately emitted or separately asserted**. ⚠ They *are* read: `N0_from_port` feeds
  `P0_raw`, which feeds `P0_physical` (`.py:224-225,231` / `.wl:109-110,116`), and `P0_raw` is additionally
  the first rank constraint (`.py:384` / `.wl:124`). But neither enters the `expressions` map the dimension
  checker walks (`.py:510-518` / `.wl:266-269`) — only `P0_physical`'s final `(0,0,0)` is asserted — and no
  record carries their triples out. ⛔ **Being unemitted and unasserted is not what singles these two out —
  per (1) it is shared by all five intermediates.** `P_port`, `Delta_port` and `K1` are absent from that
  same map; they exist only as keys of the build-internal return dicts (`.py:227-228`, `:321` /
  `.wl:112-113,122`), which nothing prints and nothing checks. ⭐ **What singles these two out is that they
  are the only two of the five for which another stage declares a dimension triple:** their triples are what
  would cross-check stage021's `SOURCED_N0_DIM` and computed `[P0_raw]` (`sympy:145,559`) and stage027's
  `[N0_den]` (`sympy:525`) — the same-object identifications D4 records in (2) — whereas **no stage declares
  a triple for `P_port`, `Delta_port` or `K1`**, so no triple of theirs could meet a counterpart even if it
  were emitted. ⚠ **That is a claim about dimension triples, not about names: all three ARE named outside
  this stage** — `P_port` at `notes/stage019_pathA33_prefactor_algebra_source_map.md:31-32` and
  `software/stage1_solver/tools/pathA_33_quadrupole_normalization_sympy.py:190-205` (declared and returned),
  `Delta_port` at `pathA_34_cross_l_unification_sympy.py:294,562`, and `K1` at
  `notes/parameter_register.md:170` and `notes/midway_knob_audit.md:126`. Those loci carry no `(L,M,T)`
  triple for them: stages 019 and 022 are two of the 13 audit scripts with **no dimension machinery at all**
  (`manifests/DIMENSION_REWRITE.md:100`), and the register rows record this stage's own `M T⁻²` rather than
  an independent counterpart. This is stage011's `OmegaDim` shape one step removed, and the uncheckable part
  is the **forgone cross-stage comparison**, not the non-emission on its own.
- **U12** — engine divergence in *how* a failure is expressed: a sum-mismatch gives a graded
  `FAIL_DIMENSIONAL` in the `.py` but aborts the `.wl`; **the only declaration corruption either engine
  probes against a walked expression** — `[M0]`, `{0,1,-1}`→`{1,1,-1}` (`.wl:289-292` / `.py:699,903`) —
  changes a *product*, so the divergence is invisible. (The other probed corruption, `[q_free]`, reaches no
  expression at all, which is what makes it the control.)
- **U13** — **a mis-binding that stays inside an exponent-degeneracy class is invisible to every instrument
  in this workstream.** 25 of the 29 emitted records share their exact triple with at least one sibling —
  only `a`, `c_s`, `D0` and `R_mix` are value-unique; the classes are `(0,0,0)`×9, `(0,1,-2)`×5,
  `(0,0,-1)`×3, `(0,1,-1)`×3, `(1,1,-1)`×3 and `(-1/2,1/2,-2)`×2. A review leg executed **five simultaneous
  same-class rebindings** in `dimension_records()` — `K0c←Z1ret`, `M0←R0`, `g_U←g_W`, `T0←P0_physical`,
  `eta_null←q_free` — obtaining stage exit 0, 111 pass, a **byte-identical** 29-record payload, and
  comparator `status=PASS|mismatches=0`. Scope of the negative: the three mechanisms this workstream owns —
  the stage's own gate, the cross-engine comparator, and a re-run. A source read *does* reach it, and one
  verified all 29 live bindings correct.
- **U14** — **under a joint perturbation the declaration freedom is larger than the one-at-a-time ablation
  shows.** A review leg mutated **16 declarations together** — `a, c_s, omega, D0, K0c, K_eta, T_Omega,
  Z0_ret, Z1_ret, Omega_U, Omega_W, R_mix, g_U, g_W, R0, R1` — and got stage exit 0, 111 pass, and all
  seven `computed_dims` **byte-identical**; with the four inert declarations that is **20 of the 22
  unpinned by anything in the Python**, leaving only `M0` and `D1` constrained, and those only relative to
  the declared `EXPECTED_DIMS`. This is **U1** (8 of 16 live declarations determined, a 24-parameter
  family) measured from the other side, not a separate discovery: a one-at-a-time ablation classifies 14 of
  these 16 as "caught", because moving one breaks a relation while moving them together preserves it.

**(5) Claims this leg believes are wrong elsewhere in the repo** (recorded, not acted on here):
**W1** register `:185` `[D0] = M T⁻²` vs three stages' `M L⁻¹T⁻²` (high that they differ, medium on which
side) — ⭐ **a contradiction, not merely a seam:** under this note's same-object identification the
register reading forces `[P0_physical]=L⁻¹`, contrary to the asserted dimensionless target, and no
normalization supplies the missing `L⁻¹`. This and `D0`'s split verdict in (1) are one finding, not two:
`D0`, `P0_raw`, and `P0_physical` are `UNDETERMINED` on the corpus basis and correct only inside this
stage's forced closure · **W2** stage009 `MOMENT0 = T⁻¹` plus an assertion that cannot fail · **W3** `DIMENSION_REWRITE.md`
§8's `g_U`/`g_W` parenthetical — ⭐ **independently recomputed by the orchestrator and CONFIRMED wrong**:
the `g` content cancels against `D0`, not against itself, and dropping the `g`s gives `L M⁻¹ T⁴`, not
`{0,0,0}` · **W4** `notes/stage023_pathA34_nullspace_underdetermination_source_map.md:250-252` states the 013/017 identification **as performed** — `K0c` ← stage013's Gate-2 collective (a,L) reduction (§8.2), *"`K_eta` = R29-derived `T_wβ²`
(013), `T_Omega` = 017's counted `T_Ω`"* — while `notes/parameter_register.md:170` carries those same objects as *unperformed*: `FREE-UNREDUCED`/`PENDING`, *"NOT identified with the raw 013/017 densities"*; a reachable file supplying the
wrong answer to this leg's question. ⚠ Still **stale pre-build, not a live claim**, but its hedges are thinner than written: `:253`'s *"likely"* qualifies the register **classification** (*"likely DERIVED manifestations, NOT new counted
knobs"*), not the identification — `:250-252` states that flatly; the *"⚠ CONFIRM at the register step"* at `:402` heads the §6 pre-read `:253` routes this into. ⚠ **The filename is load-bearing:** the three loci resolve only in that 478-line
stage023 map; in the *stage019* prefactor map this note names at `:679`, `:250-252` is unrelated `prefactor_ablation` text and `:402` does not exist (330 lines). §1.7(2b) is unaffected: it settles only which NAME to emit · **W5** the note
and both engines call these dimensions "sourced" though stage008 explicitly declines to dimension
`M0`/`D1`/`Q2`; "sourced" here means only "declared in this stage's own dict".

⚠ **Scope:** none of this disturbs `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`. The leg re-derived the `v₁=1/2`
consumption algebra by hand — `A1 − expected_A1 = i(a³ω³/c_s³)·D1·(ε₁/(1+ε₁))·(v₁ − 1/2)`, zero iff
`v₁ = 1/2` — and per **U6** the nullity is unit-invariant. What it disturbs is the belief that this
stage's dimensional gate **pins** the ledger's dimensions: per (3) it constrains *relations among* the
declarations and admits a 24-parameter family of them, so a shared wrong declaration survives it.

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

Both engines are standalone, assert-zero, and exchange nothing: the `.wl` remains print-only with ZERO file I/O, while the
`.py` writes its deterministic labelled-dimension sidecar (no scratch YAML, report, note, card, LaTeX or registration write).
The `.wl` is a **genuinely independent route,
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

### 5.1 Dimension-rewrite steps (g)/(g2) and (h)/(h2) — executed evidence

*Seated here rather than in a new section: §5 is this stage's verification home, and these are verification legs.*

- **(g)/(g2), both engines.** `python3 scripts/compare_dimension_artifacts.py 023` → exit 0,
  `py=29|wl=29|shared=29|py_only=0|wl_only=0`, `mismatches=0`, waivers empty. The orchestrator re-ran `math -script` on the
  committed `.wl`: exit 0, and the transcript reproduces **byte-identically** after `sed -E 's/\$[0-9]+/$N/g'`, carrying 29
  `DIM|` records. The orchestrator also regenerated the Python sidecar: the emitted bytes are identical to the committed artifact.
- **Two orchestrator ablation axes**, both re-run with per-record captures retained and the emitted sidecar deleted before every
  iteration. **Axis 1** — one declaration replaced by a decoy that occurs nowhere in the source: **22 of 22 detected**, 16 of them
  by the stage itself at exit 1 with no sidecar emitted, and **6** (`R0`, `R1`, `eta_null`, `gain0`, `gain1`, `q_free`) only by
  the comparator, the stage staying green at 111 PASS. **Axis 2** — one emitted record repointed at a different-valued source:
  **29 of 29**, every one leaving the stage green and named by the comparator alone. ⚠ **The qualifier that makes axis 2 honest:**
  its decoys are always value-distinct, so it establishes **cross-class binding detection only**; the same-class case is **U13**.
  Evidence in-repo: `research/pde_ledger_v2/notes/stage023_step_h_evidence/` — the method-and-verdict summary
  (`ABLATION_SUMMARY.md`), the 51-row per-target results table (`results.tsv`, 22 axis-1 + 29 axis-2 rows) and the
  orchestrator-owned include list (`include_list.tsv`, the same 51 targets with their old/new text). ⚠ **What did NOT survive:**
  the 51 per-record `.cmp`/`.stdout` captures, the ablation driver, and the `179`-hit candidate ledger
  (`research/pde_ledger_v2/_scratch/stage023/enum/candidate_ledger.md`) remain in gitignored scratch — so the committed evidence is
  the summary and the two tables, **not** the raw captures those tables were read from.
- ⛔ **In none of axis 1's 16 caught rows did the dimensional gate's own assertion fire first.** `base_verdict` ranks
  `dimensional` **third** in the verdict ladder (`.py:674-682`, `dimensional` at `:677`, behind `decoupled` and `tautological`),
  so a corrupted declaration flips the gate verdict before `run_dimensional_gate()` (`.py:900`) is reached; what fires instead is
  the selector control (11 rows) or the baseline earned-verdict check (5 rows), both in `run_native_rank_and_selector()`, which
  `main()` calls first (`.py:1075`). The detection is real; what it is **not** is this stage's dimensional assertions doing it.
- **Two fresh review legs and two remediation rounds.** Round 1 removed an unreachable `fmt` branch, activated `assert_no_float`
  over the 29 records, and regenerated the stale committed transcript. A fresh re-verify leg then measured that the activated
  guard **cannot fire on any input the stage can construct**: every `Dimension` funnels its exponents through the module's
  `_exact` → `sp.Rational`, so `Dim(0.5,0,0)` becomes `1/2` and `from_mapping({L: 0.1})` becomes
  `3602879701896397/36028797018963968`, silently. Round 2 therefore **reverted both edits**, leaving stage023 in its siblings'
  shape. What remains from that round is the ordering fix: the sidecar is emitted **before** the terminal status banner
  (`emit_dimension_sidecar` at `.py:1088`, ahead of `print_verdict_labels()` at `:1092` which prints `AUDIT_STATUS=PASS` at
  `:1067`), so a fault during emission can no longer print `AUDIT_STATUS=PASS` beside `OVERALL FAIL`.
- **(h2) — the sealed step-(e) prediction.** Opened only after both review legs and both remediation rounds had landed, and
  committed unchanged with its adjudication at `research/pde_ledger_v2/notes/stage023_py_rewrite_prediction.md`:
  **5 fully confirmed · 1 falsified (P2) · 1 split (P3: mechanism confirmed, exclusivity falsified).**
  The falsified one is **P2**, which predicted the emitter would key off `symbol.name` so that a
  name defect would be loud; instead `dimension_records()` hand-types each join key beside an independently written source
  expression (`.py:539-571`), which is exactly what makes **U13**'s silent same-class rebinding constructible.
