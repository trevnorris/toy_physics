# The REPORTED-PHYSICAL-RESULT set — stages 016, 023, 043

**Status:** the artifact required by `CENSUS_SCHEMA.md` §7.1.1 step 1 and §5.8 output 17 — the
enumerated result set, one entry per result, each with its locus. It is the **starting set** from which
the census universe is later derived by walking closures (§7.1.1 step 2).

⛔ **This file classifies nothing.** No census row, no tier, no axis value, no derived/declared
judgement, no count of knobs. It answers exactly one question, per stage: *what does this stage present
as a statement about the modelled world?*

⛔ **Scope: stages 016, 023, 043 only.** The set for the remaining stages is not taken here, and this
file is not a partial census.

⛔ **Not amendable after classification begins** (§7.1.1): a result set that can be adjusted once the
rows are visible is a knob on the denominator. Corrections before classification are ordinary; changes
after it are not.

---

## 0. The test applied, and how it was applied

`CENSUS_SCHEMA.md` §7.1.1. An entry is a **reported physical result** iff all three hold:

- **(a)** presented as an **outcome** — in a stage note, a paper section or a script emission — not as
  an intermediate on the way to one;
- **(b)** a claim about **the medium or its phenomenology** — a dimensionless ratio compared or
  comparable to observation, the form / sign / power of a force law, an equation of motion or one of its
  coefficients, a spectrum, a stability statement about the medium, a predicted relation between
  physical quantities;
- **(c)** ⭐ the discriminator — **it would still be a claim about the world if every bookkeeping
  artifact of this ledger were deleted.**

And the spec's own exclusion list: dimension exponent vectors and dimensional CORRECT/WRONG/UNDETERMINED
verdicts, coverage figures and row counts, gate pass/fail flags, ID manifests and reconciliation deltas,
category tallies, and the census's own outputs.

**Method.** Every locus below was opened and read; none is carried from an attribution in another
document (§9.1 rule 1). Loci are full paths (§7.2). Line numbers are as of the reading; where a claim
spans a block, the block is given.

**Two conventions used throughout, stated so they can be contested.**

1. ⭐ **A mathematical identity is not automatically a physical result.** Both `Gram = I₅` and
   `λ_m = 6` are true of the real ℓ=2 spherical harmonics whether or not this model exists. Test (b) is
   what separates them here: the eigenvalue becomes **a coefficient of the medium's ℓ=2 equation of
   motion**, the orthonormality does not. Where this distinction decided an entry, it is said so
   explicitly.
2. ⭐ **A bundle exported under one name may still split.** Stage016 exports "the ℓ=2 SO(3) covariance
   theorem" as a triple; the triple's members are judged individually, and the demotion of one member is
   flagged rather than hidden.

⚠ Each stage's negative section carries a **⚠ NEAR-MISS** marker on the entries a reviewer is most
likely to want to contest. Those are the deliberate seams of this artifact.

---

## 1. stage016 — ℓ=2 SO(3) covariance

Note: `/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage016_l2_so3_covariance.md`
Script: `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py`

### 1.1 Reported physical results — **2**

#### `RES:016:L2-IRREP`

> **Claim.** The frozen throat's ℓ=2 angular response sector is a single 5-dimensional SO(3)
> irreducible representation — all five real ℓ=2 channels share one `−Δ_S²` eigenvalue `λ_m = 6` — so
> the ℓ=2 response is m-degenerate, i.e. angularly isotropic.

| | |
|---|---|
| **loci (note)** | `.../notes/stages/ledger_stage016_l2_so3_covariance.md:43-53` (§1.2; the degeneracy sentence at `:52-53`), restated as the earned export at `:335-336` and as the stage's opening claim at `:24-25` |
| **loci (script)** | `.../scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:287` (Rayleigh quotient computed), `:290` (eigenfunction residual), `:299` (`lambda_all_six`), asserted `:679-682`, emitted `:674`, `:825`, `:864` |
| **(a)** | ✅ presented as an outcome — §1.2 is titled "the crux", and the triple is the stage's **Exported** packet (`:335-336`) |
| **(b)** | ✅ **a spectrum** — the eigenvalue set of the angular operator acting on the medium's ℓ=2 sector — and, via `RES:016:K2-FORM`, an equation-of-motion coefficient |
| **(c)** | ✅ delete every ledger bookkeeping artifact and "the throat's ℓ=2 response is degenerate across the five m channels, so that sector is isotropic" is still a statement about the modelled medium |

⚠ **What is *not* in this entry.** The number 6 is `ℓ(ℓ+1)` and is mathematics. What is claimed about
the world is the **degeneracy** — one eigenvalue shared by all five channels of the medium's ℓ=2
response — and that is what carries into the stiffness below.

#### `RES:016:K2-FORM`

> **Claim.** The ℓ=2 angular stiffness and mass of the frozen throat wall are `K₂ = K̃ + λ_m·T̃_Ω` and
> `M₂ = M̃`, with the `T̃_Ω` coefficient being the computed `−Δ_S²` eigenvalue rather than an
> independently chosen constant.

| | |
|---|---|
| **loci (note)** | `.../ledger_stage016_l2_so3_covariance.md:55-60` (§1.3; the display block at `:57-58`), in the exported packet at `:335-336` |
| **loci (script)** | `.../ledger_stage016_l2_so3_covariance_sympy_audit.py:304-305` (`build_K2`), `:499` (assembled per channel from the live `lambdas`), `:500-504` (extracted-coefficient residual), asserted `:693-697`; emitted `:688-691`, `:864` |
| **(a)** | ✅ presented as an outcome and exported to stage017 (`note:335-337`) |
| **(b)** | ✅ **a coefficient of an equation of motion** — the ℓ=2 stiffness of the wall's reduced dynamics |
| **(c)** | ✅ "the wall's ℓ=2 restoring stiffness is its radial stiffness plus six times its angular stiffness scalar" survives deletion of the ledger's bookkeeping |

⚠ **Scope carried with the claim, from the stage's own §3 (`:311-314`):** the **angular** structure is
what is earned; the radial profile `β₂(w)` and the radial scalars `M̃`, `K̃`, `T̃_Ω` remain frozen
calibration inputs. The claim above is about the angular form and its coefficient, not about the values
of the three scalars.

⚠ **Recorded against this entry, not resolved here:** the stage's own adversarial finding **H8**
(`:259-270`) is that the dimensional block is a *parallel reconstruction* of these expressions and that
dimension-preserving rewrites of the walked expressions — including dropping `λ_m` from the angular
term — leave the audit green. That bears on how well the artifact protects this result; it does not
change whether the result is reported.

### 1.2 NOT reported physical results — **11 listed**

| # | item | locus | why it is not a reported physical result |
|---|---|---|---|
| 1 | ⚠ **NEAR-MISS** — `Gram = I₅`, the orthonormality of the five real ℓ=2 harmonics, and the unit self-overlaps | note `:27-41`; script `:281`, `:513`, asserted `:676-677`, `:742-744` | Passes (a) — it is inside the exported triple (`:335`) — but fails **(b)** and **(c)** *as a claim about the medium*: orthonormality on S² is a property of the sphere and of the basis this artifact writes down, true independently of the medium. Functionally it is a **normalization self-check** on the artifact's own harmonics (and the Gram diagonal feeds the probe machinery, `:66`). ⚠ Contestable: it is exported under the same name as `RES:016:L2-IRREP`. |
| 2 | the `(L,M,T)` dimension records — the 12 `dim_rules.*` declarations and the 9 `baseline_dims.*` walked results | script `:312-326`, `:462-489`, asserted `:716-726`; note §1.4 `:73-83` | §7.1.1 explicitly excludes **dimension exponent vectors**; §11 confirms dimensional verdicts are a different axis from provenance. |
| 3 | the physics-leg per-quantity dimensional verdict, "21 of 21 CORRECT on this stage's own convention" | note `:175-182` | Excluded by name — a **dimensional CORRECT/WRONG/UNDETERMINED verdict**. |
| 4 | the dimension-object enumeration (the `.wl` object table, its membership rule and its emitted/not-emitted dispositions) | note `:85-149` | **Coverage figures and row counts** about the artifact's own emission set. |
| 5 | coverage findings — "12 of 21 records are declared literals in BOTH engines", "0 are computed from any physical input", "21 records but only 12 FREE VALUES … read the census as 12 free / 9 derived" | note `:194-207` | Self-audit / coverage figures about the artifact. Test (c): delete the records and nothing about the world remains. ⚠ These are *evidence for* later census axes; they are not results. |
| 6 | the three able-to-fail probes (`wrong_eigenvalue`, `tautology_hash_collision`, `dimensional_corrupt_T_Omega`), their self-ablations and the aggregate battery | note `:293-305`; script `:516-627`, `:698-703`, `:728-735`, `:746-767` | **Controls and deliberate negatives** (§7.3) — constructions built to prove a check able-to-fail. |
| 7 | the per-tooth ablations (corrupted basis, scaled harmonic, forced coefficient 2, neutered distinctness) | script `:770-816` | Same: controls, run on copies, never asserted about the medium. |
| 8 | the verdict token `ISOTROPY_CALIBRATED`, the earned-label `L2_SO3_COVARIANCE_THEOREM_EARNED`, and the tallies "SymPy 82 PASS / Mathematica 91 PASS" | note `:8-14`; script `:819-830`, `:859-867`, `:899-906` | **Gate pass/fail flags** and run tallies — verdicts about the ledger, not about the medium. |
| 9 | the register accounting — "ZERO new counted knobs", the structural edge **R34**, "Part-II counted CALIB set unchanged … = 4", the deferral of `T_Ω`/`β₂` counting to 017 | note `:330-334` | **Category tallies** and ledger-internal edges; the census's own kind of output. |
| 10 | **H1–H10**, the "structurally uncheckable" deliverable, and the step-(f) 4000-trial differential-harness result | note `:151-161`, `:209-283` | **Self-audit** about what the artifact's own checks can and cannot see. Genuine findings; they are about the instrument, not the medium. |
| 11 | ⚠ **NEAR-MISS** — `T_Ω = T_w/a²` for an isotropic wall, and the identity hazard `T̃_Ω`(016) vs `T_Ω`(023) "different quantity until R42 exists" | note `:180-182`, `:187-192`, `:228-229` | Both are **claims about the world in form** but fail **(a)**: neither is presented as an outcome. The first is presented as an explicitly **UNDETERMINED** possibility ("no dimensional check here can support or refute"); the second is an **identity adjudication** about ledger symbols, resolved only conditionally on an unstated ℓ-independence of the radial profile. ⚠ If a later pass finds either asserted as an outcome elsewhere in the corpus, it belongs to *that* artifact's result set. |

---

## 2. stage023 — native nullspace underdetermination

Note: `/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md`
Script: `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`

### 2.1 Reported physical results — **3**

#### `RES:023:RETURN-UNDERDETERMINED`

> **Claim.** The linearized Gate-5 return law of this medium does not determine the ℓ=0/1 brane↔bulk
> return: over the 11 generator dofs the three collected named constraints have rank 3, and the two
> return admittances `{Z0_ret, Z1_ret}` are untouched by every one of them (each constraint's
> `∂/∂Z_ℓ,ret = 0`) while still moving the DC return transmissions `T0`/`T1` — a genuine 2-parameter
> return freedom (`return_moving_nullity = 2`), so the ℓ=0/1 ↔ ℓ=2 return is not predicted at linear
> order.

| | |
|---|---|
| **loci (note)** | `.../ledger_stage023_nullspace_underdetermination.md:40-61` (§1.1; the constructive witness at `:56-59`), headline verdict at `:5-18`, honest scope at `:771-773`, export at `:821-824` |
| **loci (script)** | `.../ledger_stage023_nullspace_underdetermination_sympy_audit.py:205-217` (the 11 dofs), `:384` (the three collected constraints), `:388-393` (genuine `sp.Matrix(rows).rank()`, nullity, return-moving nullity), `:395-413` (the two constructive `e_{Z_ℓ,ret}` witnesses and their `ΔT0`/`ΔT1`), asserted `:809-812`, `:862`, `:1036`; emitted `:802-808`, `:1068` (`PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE`) |
| **(a)** | ✅ the stage's headline outcome, and the ledger itself separates it from the bookkeeping flag: `AUDIT_STATUS=PASS` (the script ran) vs `PHYSICS_VERDICT=FAIL_…` (the earned physics) — note `:13-18`, script `:1067-1068` |
| **(b)** | ✅ a claim about **the medium's equations of motion** — specifically that the linear return law supplies no equation fixing two of its own dofs, witnessed constructively by directions that preserve every constraint yet change a physical transmission. ⚠ *Interpretive note:* §7.1.1(b)'s list has no explicit entry for a **negative** determination claim; this is read under "an equation of motion or one of its coefficients" (here: the absence of the fixing equation, plus the nonzero `ΔT_ℓ` gradients). Flagged so a reviewer can contest the mapping, not the entry. |
| **(c)** | ✅ delete the ledger's bookkeeping and "in this medium the linear brane↔bulk return law leaves a two-parameter freedom in the ℓ=0/1 return admittances" is still a statement about the modelled world — and it is a **falsifiable** one: a nonlinear closure supplying two independent equations would overturn it |

⚠ **Carried with the claim, from the stage's own §3 (`:778-781`):** the raw native nullity `8` is
**partly bookkeeping** — two of the three collected constraints are self-constraint rows — and the
stage states that the verdict-bearing quantity is `return_moving_nullity = 2`, with the raw nullity a
reported diagnostic. `rank0 = 3` and `native_nullspace_dimension = 8` are therefore listed here as
*components of this result's evidence*, not as separate results.

#### `RES:023:A0-FORM`

> **Claim.** The ℓ=0 (scalar) return-residual amplitude of the medium is
> `A0 = i·v₀·(aω/c_s)·M0·(1−T0)`, which coincides with the independently built pathA_29 residual form
> `i·aω·M0·ε0/(c_s(1+ε0))` and is realized at order `ω¹` — the residual's form, sign and ω-power are
> claimed; its magnitude is not.

| | |
|---|---|
| **loci (note)** | `.../ledger_stage023_nullspace_underdetermination.md:77-94` (§1.3; the forward build at `:82`, the independent comparison form at `:87`, the class/magnitude split at `:94`), export at `:823` |
| **loci (script)** | `.../ledger_stage023_nullspace_underdetermination_sympy_audit.py:342` (`A0`), `:345` (`expected_A0`), `:349` (residual), `:351` (realized ω order), asserted `:882`, `:884`; emitted `:877`, `:879`, `:881` |
| **(a)** | ✅ an outcome — §1.3 is a "what this stage earns" subsection, and the `A0/A1` residual class is in the **Exported** list (`:823`) |
| **(b)** | ✅ **the form, sign and power of a radiative return law** — explicitly a power-of-ω statement about an amplitude |
| **(c)** | ✅ "the scalar return residual goes as `ω¹` and as `ε/(1+ε)` in the return ratio" is a claim about the medium's radiative behaviour independent of any ledger bookkeeping |

⚠ **Conditions the ledger attaches, carried here:** the class is earned **conditional on a positive
bounded branch** (`script:273-280`, `:316`), and the **magnitude is DEFERRED** because the nullspace of
`RES:023:RETURN-UNDERDETERMINED` leaves `ε_eff` free (note `:94`, `:782-784`). This entry is the
form/sign/order claim only.

#### `RES:023:A1-FORM`

> **Claim.** The ℓ=1 (dipole) return-residual amplitude of the medium is
> `A1 = i·v₁·(aω/c_s)³·D1·(1−T1)`, which coincides with the independently built pathA_29 dipole form
> `i·a³ω³·D1·ε1/(2c_s³(1+ε1))` and is realized at order `ω³` — again form, sign and power only.

| | |
|---|---|
| **loci (note)** | `.../ledger_stage023_nullspace_underdetermination.md:77-94` (the dipole half of §1.3; the `v₁ = 1/2` consumption algebra at `:89-92`) |
| **loci (script)** | `.../ledger_stage023_nullspace_underdetermination_sympy_audit.py:343` (`A1`), `:346-348` (`expected_A1`), `:350` (residual), `:352` (realized ω order), asserted `:883`, `:885`; emitted `:878`, `:880-881` |
| **(a)** | ✅ same as above |
| **(b)** | ✅ **the form and power of a radiative return law** — the `ω³` dipole scaling against the `ω¹` monopole |
| **(c)** | ✅ same as above |

⚠ **Why this is listed separately from `RES:023:A0-FORM`.** They are two relations, over different
inputs (`M0`/`v₀`/`T0` vs `D1`/`v₁`/`T1`), asserted by four separate teeth, and their closures differ.
The ledger presents them jointly as "the `A0/A1` residual class"; splitting is at the relation level and
is recorded here so a reviewer can rejoin them if that is judged wrong. ⚠ The `1/2` in `v₁` is the
**checkable** consumption of stage022's coefficient (`note:89-92`; ablated at `script:891-895`) — the
coefficient itself is stage022's result, not this stage's.

### 2.2 NOT reported physical results — **13 listed**

| # | item | locus | why it is not a reported physical result |
|---|---|---|---|
| 1 | ⚠ **NEAR-MISS** — the counterfactual selector `{Z0_ret = K0c, Z1_ret = K_eta + 2·T_Omega}` and the flip to `CROSS_L_RESIDUAL_PREDICTION` | note `:63-75`; script `:235-243`, `:246-266`, `:848-864` | A **control / deliberate negative** (§7.3) — the stage says so itself: `derived_from_named_pde = False`, `control_only = True`, "a COUNTERFACTUAL RANK-COLLAPSE WITNESS … NOT a proven Gate-6 selector" (`:72-75`, `:774-777`). Its role is to prove the gate able-to-fail. ⚠ It **looks** like a physical prediction about the return law; it is explicitly withheld as one. |
| 2 | ⚠ **NEAR-MISS** — the exported "Gate-6 need": two independent equations must fix `{Z0_ret, Z1_ret}` (edge **R42**) | note `:74-75`, `:774-777`, `:821-824` | An **obligation / debt statement** addressed to the ledger's future work. Its world-facing content is already carried by `RES:023:RETURN-UNDERDETERMINED` (the freedom is 2-dimensional); as a separate entry it would double-count one claim. |
| 3 | the raw native nullity `8` and `rank0 = 3` as standalone headline numbers | note `:47-48`, `:778-781`; script `:809-810` | The stage itself de-rates the raw nullity as **partly bookkeeping** and names `return_moving_nullity` as the verdict-bearing quantity. Recorded as evidence inside `RES:023:RETURN-UNDERDETERMINED`, not as its own result. |
| 4 | the ℓ=2 port kernel `P0_raw = (Ω_U²g_W + R_mix g_U)²/(Ω_U²Ω_W² − R_mix²)²/D₀`, the ℓ=0 collective stiffness `K0c`, and `K1 = K_eta + 2·T_Omega` | note `:42-44`, `:800-802`; script `:220-232`, `:236`, `:284-285` | **Consumed inputs** — the three collected named constraints, cited from the handoff §10.2–10.3, Gate-2/§8.2 and §9.4. They are inputs to this stage's rank computation; their physical content belongs to the stages that earned them. |
| 5 | the transfer definitions `T_ℓ = K_ℓ/(K_ℓ + Z_ℓ,ret)`, `ε_ℓ = Z_ℓ,ret/K_ℓ`, and the `forward_relations_ok` identity | note `:50-51`, `:84`, `:765`; script `:302-315`, `:886-889` | **Definitions plus a self-consistency identity.** The stage de-counts the identity itself: a labeled `SELF-CONSISTENCY` check, not a verdict tooth (`:765`, `:870-872`). |
| 6 | pathA_29's `Z_is_premise = True`, `boundary_dof = none` | note `:59-61`, `:796-797`; script `:425` | **Consumed provenance**, and de-counted by the stage to a labeled `PROVENANCE` print (`:764-765`). |
| 7 | the `(L,M,T)` dimensional gate — `[A0]`, `[A1]`, `[T_ℓ]`, `[ε_ℓ]`, `[P0_physical]`, and the 22 sourced + 7 expected declarations | note `:96-104`; script `:429-537`, `:900-917` | Excluded by name: **dimension exponent vectors**. |
| 8 | the 34-object per-quantity dimensional verdict table and its tallies `24/0/10` and `27/0/7` | note `:379-452` | Excluded by name: **dimensional CORRECT/WRONG/UNDETERMINED verdicts**. ⚠ Includes the individually interesting `[D1]/[M0] = L` multipole-ladder remark (`:411`), which is a dimensional-consistency statement, not a reported outcome. |
| 9 | the dimension-object enumeration — 42 rows, 179 predicate hits, 29 proposed records, the `ENUM_COUNT` block | note `:114-365` | **Row counts and coverage figures** about the artifact's own emission set. |
| 10 | the coverage finding "29 of 29 dimension values are hand-typed literals in BOTH engines", and **U1–U14** | note `:603-625`, `:627-719` | **Self-audit** about the instrument: what the dimensional gate can and cannot pin (the 24-parameter family, the same-class rebinding, the joint 16-declaration perturbation). Evidence for later census axes; not claims about the medium. |
| 11 | **W1–W5** — claims this leg believes are wrong elsewhere in the repo, incl. the `[D0]` contradiction, and the D4 name determinations / merge refusals (`K_eta`/`T_Omega` vs 016's tildes, `c_s` vs `c_S`, `M0`/`D1`) | note `:478-534`, `:536-601`, `:721-736` | **Identity adjudications and defect reports about the ledger's own declarations.** Delete the declarations and they say nothing. ⚠ The `β₁ ≡ β₂` question underneath the merge refusal *is* a physical question; nothing here asserts an answer. |
| 12 | the `FAIL_TAUTOLOGICAL` firewall, the transfer/return probes (`wrong_sign_return`, `perfect_return`, `no_consistent_return`, `decouple_knobs`), the R1 port-kernel dependency probe, and the verdict read-set assert | note `:106-112`, `:746-765`; script `:759-794`, `:925-1016` | **Controls and deliberate negatives** (§7.3), each a mutation built to make a check fail. |
| 13 | `AUDIT_STATUS=PASS`, the tallies `111 PASS` / `117 PASS`, the register accounting (ZERO new counted CALIB knobs, `Z0_ret/Z1_ret` as coordinate aliases, `K0c/K1` as counted `FREE-UNREDUCED` debt, edge **R42**), and the step-(g)/(h) ablation-axis results | note `:15`, `:23-24`, `:805-820`, `:877-916`; script `:1067` | **Gate flags, run tallies, category tallies and register edges** — the ledger's statements about itself. |

⚠ **Explicitly withheld, and recorded as withheld:** no `ε_eff` magnitude and no nonzero return
prediction is claimed (`:94`, `:782-784`, `:106-112`). The absence is deliberate and structurally
enforced; it is not an omission for a later pass to fill.

---

## 3. stage043 — irreducible count range

Note: `/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage043_irreducible_count_range.md`
Script: `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage043_irreducible_count_range_sympy_audit.py`

### 3.1 Reported physical results — **ZERO**

⭐ **Stage043 reports no physical result. This is a finding, not a gap, and nothing was manufactured to
fill the slot.**

The stage says so in its own words. Its **KEY FRAMING** block states the scope class outright:

> "Scope class: **STRUCTURAL / dimensionless** — a codimension count carrying no `[L,T,M]` and no new
> knob; it is a **bookkeeping object** that FEEDS Part VII's unified count, it does not itself discharge
> or mint an edge." — `.../ledger_stage043_irreducible_count_range.md:30-34`

Its two reported quantities — the continuous codimension range `[40, 49]` and the discrete
structural-postulate count `11 (= 7+4)` (`:57-62`; script `:748-772`, `:775-787`) — are **category
tallies over the ledger's own register rows and audit buckets**, transcribed into the engines as
self-contained facts (`:411-416`; script `:5-7`, `:137-176`). §7.1.1's exclusion list names exactly this
kind of object: *coverage figures and row counts, ID manifests and reconciliation deltas, category
tallies, and the census's own outputs*. Test **(c)** is decisive: delete the register, the midway audit
buckets and the manifest, and the number `[40, 49]` has nothing left to count — it is a statement about
the ledger's parameter bookkeeping, not about the medium.

⭐ The count is in fact the **predecessor of the census's own headline** (`CENSUS_SCHEMA.md` §1, §10.1
and §5.6 — which cites this very note as the range idiom's source). A census that admitted it as a
physical result would be measuring itself.

⚠ **Stated so it is not mistaken for a verdict on the stage:** "no reported physical result" is a
statement about *what kind of object this stage produces*, not about whether it produces it well. The
stage is dual-engine, 20 teeth, per-tooth-ablated, and its counting rule is explicit.

### 3.2 NOT reported physical results — **10 listed**

| # | item | locus | why it is not a reported physical result |
|---|---|---|---|
| 1 | the headline `IRREDUCIBLE_COUNT_RANGE` object — `continuous_codimension [40, 49]`, `base_continuous [27, 36]`, `E_continuous 13`, `spread 9`, `C1`/`C2` cardinalities 3/6, `convention OPEN` | note `:36-62`, `:140-168`; script `:748-772`, `:775-787`, `:1109-1140` | **Category tally / row count** over the register. Fails (b) — it counts the ledger's knob rows — and fails (c) decisively: with the register deleted there is nothing to count. |
| 2 | the discrete structural-postulate count `11 (= 7+4)` and its itemization D1–D4 (`H` existence + Pöschl–Teller law-form, `±w`-puncture `Z₂` topology, mouth BC-class, `τ_d` time-arrow) plus the 7 Parts-I–II postulates | note `:102-115`, `:170-192` | The **count** is a tally. The **postulates themselves** are real physical content, but they are cited **by register name** as members of a count; each is stated as a postulate in the stage that owns it (030/031/032/…, wall/G0/EOS). This stage adds no claim about them beyond membership and orthogonality of the count. ⚠ Whichever stage *states* a postulate owns the corresponding result entry. |
| 3 | ⚠ **NEAR-MISS** — the `Δr` independence diagnostic: block M `(10, 5, 5)`, block Wχ `(10, 7, 3)`, token `CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS` | note `:228-258`; script `:603-696` (blocks `:606-644`, baselines `:646-653`), asserted `:1024-1034` | The closest candidate in the stage, and it still fails. It computes the Krull dimension of a parameter variety after imposing relations — but the **relations are imported** (`m c_s0² = 5Kρ0⁴`, `a m c_s0 = ħ`, `K_η = T_w β²`, …; script `:609-615`, `:632-636`), each earned in another stage, and the **statement made here** is that the *counted knobs have no hidden multiplicity* — a property of the **count**, i.e. of the census's own §10.1 denominator. The stage itself scopes it as "a **DIAGNOSTIC** that discharges NO knob and does NOT shrink the count" (`:254-258`). Fails (b) as a claim about the medium and (c) as a claim that survives deleting the knob list. ⚠ A reviewer who reads "these model parameters are independent" as phenomenology should contest this entry specifically. |
| 4 | the four Parts-I–II buckets and the raw tally `34` / `43` | note `:84-100`; script `:451-535` | **Category tally** checked against a ratified oracle. |
| 5 | the coherent counting rule `count = nominal − DERIVED-and-EXECUTED − CONV − external-bridge` and its riders (pending debt stays counted, imposed calibrations stay counted, a departure adds no knob) | note `:64-82` | A **counting rule** — method for the ledger's bookkeeping, quoted by `CENSUS_SCHEMA.md` §5.2 as an existing counting convention. |
| 6 | the `REGISTER_TO_COUNT_MANIFEST` (universe ≈ 152; the ten disjoint categories and their counts) | note `:260-285`; script `:137-176`, `:430-449`, `:1096-1106` | An **ID manifest** and category tally — named in §7.1.1's exclusion list, and the object §8.4 tells the census to *reconcile against*, never to inherit. |
| 7 | the reduction-debt sub-count `[18, 21]`, the R49 parallel-track out-of-scope line (8, or 9 with `θ_B`), the departure-no-knob count 4 | note `:194-226`, `:269-281` | **Sub-tallies and scope decisions** about which rows enter the count. |
| 8 | the R35 reconciliation (`DERIVED-in-form` label kept + `C1`-contributor rider) | note `:287-294`; script `:1057-1075` | A **label / provenance reconciliation** about a register row. |
| 9 | the two continuous and two discrete sensitivities (`Q_E`-undeclared total-neutral 3→4; `K_θ`/`κ_phase` as 2 DOF → `[41, 50]`; the MacCullagh law-form and H-closure granularity each `11 → 12`) | note `:296-314`; script `:704-742` | **Sensitivities of the count** to bookkeeping choices — how the tally moves if a row is split. ⚠ Each names a real physical question (is `K_θ` one knob or two?), but the reported quantity is the count's response, and the note states they are "printed, NOT folded into the headline". |
| 10 | the 20 teeth, their per-tooth ablations (`LEDGER_STAGE043_MUTATION`), the dual-engine agreement, and the "what 043 does NOT do" ledger-accounting block | note `:315-390`; script `:32-96`, `:111-124` | **Gate flags, controls and self-audit.** |

---

## 4. Summary of this artifact

| stage | reported physical results | listed non-results |
|---|---:|---:|
| 016 | **2** (`RES:016:L2-IRREP`, `RES:016:K2-FORM`) | 11 |
| 023 | **3** (`RES:023:RETURN-UNDERDETERMINED`, `RES:023:A0-FORM`, `RES:023:A1-FORM`) | 13 |
| 043 | **0** | 10 |
| **total** | **5** | **34** |

⛔ These are counts of *results*, not census counts. No tier, no axis and no derived/declared judgement
appears anywhere above.

---

## 5. Limits of this pass — what could not be applied, and what a reviewer should attack

1. ⭐ **Test (c) has no procedure for a mathematical fact.** "Would it still be a claim about the world"
   is answerable for a physical relation and ambiguous for an identity that is true regardless
   (`Gram = I₅`, `λ_m = ℓ(ℓ+1)`). It was resolved by leaning on (b) — does the quantity enter the
   medium's equations of motion — and every entry decided that way is flagged. A reviewer who disagrees
   with the `Gram` demotion would make stage016 3 results.
2. ⭐ **§7.1.1(b) has no clause for a negative determination.** `RES:023:RETURN-UNDERDETERMINED` is a
   claim that the theory *fails* to fix something. It was admitted under "an equation of motion or one
   of its coefficients"; the mapping is flagged in its entry.
3. ⚠ **(a) "outcome vs intermediate" is not fully decidable inside one stage.** Stage016 exports a
   three-member bundle under one name; stage023 exports an "`A0/A1` residual class" that this file
   splits in two. Both calls are recorded with their reasoning rather than hidden.
4. ⛔ **Cross-stage duplication cannot be resolved from three stages.** §7.1.1 itself says membership is
   not decidable from one artifact and needs one corpus-wide pass. If another stage presents the same
   claim as *its* result (e.g. stage017 riding the ℓ=2 covariance, stage022 owning `{1, 1/2}`, stage009
   owning the `ε/(1+ε)` residual form), the full result set will have to adjudicate whether those are
   one result or two. Nothing here presumes an answer.
5. ⚠ **The consumed/exported distinction was taken from each stage's own §4.** Where a stage cites an
   input as provenance, this file treats that input's physical content as belonging to its source stage.
   That is a judgement, made uniformly, and it is what removes `P0_raw`, `K0c`, `K1`, `v₀`, `v₁` and the
   pathA_29 residual form from stage023's positive list.
