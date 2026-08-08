# S10 · Two transverse photons — and the count comes from the brane's own dimension

**Sector 1 (light), step 2.** Walked side by side with the user, 2026-08-02.
**Engines, harness and verification rebuilt 2026-08-05/06** under the script-rebuild programme.
⚠ **The physics of the walk did not change. The evidence for it did**, and in one place the old evidence
was worth nothing — see `§Q7` below.

---

## What it is

S9 gave light a cone. S10 asks how many independent waves the brane's transverse sector actually
carries, and — the part that matters — **why that number**.

## What it does

Establishes that the mode count is `D_brane − 1`, computed rather than inferred, and registers the
brane's dimension as the input it has silently been all along.

## The walk, move by move

**Move 1 — decompose the amplitude about `k`.** For a plane wave `u = a·exp(i(k·x − ωt))`, the constant
amplitude `a` splits into a piece along `k̂` (one number) and a piece perpendicular to it. In `D`
dimensions the perpendicular subspace has `D − 1` directions:

```
D = 1 + (D − 1)
```

⚠ **The identification made here, flagged before it was used:** `u` is the **in-plane** displacement, a
`D`-vector on the brane — ⛔ **not** a `D+1`-vector including motion out of the brane into `±w`. The
out-of-plane motion is the **`h`-branon**, a separate object carrying charge (`h ≠ u_L`, stated in three
places in the corpus). ⭐ **User confirmed the picture.** ⛔ If that were wrong — if out-of-plane motion
belonged to the same elastic sector — the transverse count at `D=3` would be **3**, not 2.
⇒ **`R-S8-02`** in `SUBSTRATE_REQUIREMENTS.md`: nothing in this ledger yet computes the decoupling.

**Move 2 — the curl sees only the perpendicular part.** `∇×u → i k × a`, so a longitudinal amplitude has
zero curl. ⭐ The clean algebraic form of this (contributed by a review leg, and it is the whole reason):

```
(1/2) Σ_{i,j} (k_i a_j − k_j a_i)²  =  |k|²|a|² − (k·a)²
```

⇒ the stiffness quadratic form **penalises only `a_⊥`**.

**Move 3 — the equation of motion, per piece.** `ρ_br ω² a = μ_R[k²a − k(k·a)]`, giving `ω² = 0` for the
longitudinal and `ω² = (μ_R/ρ_br)k²` for the perpendicular directions.

⭐ **The longitudinal mode is not forbidden — it is FREE.** Zero curl means zero energy cost, so no
restoring force. ⭐ That is precisely what *curl-only* means: the medium resists twisting and is
indifferent to compression along the wave. An ordinary elastic solid resists both, and carries a
longitudinal wave.

⭐ **The two transverse modes are degenerate by SYMMETRY**, not coincidence: nothing in the equation
refers to *which* perpendicular direction `a_⊥` points, so both carry the same `ω`.

## ⭐⭐ The computed result — and the claim it settles

⛔ **The count is not read off a multiplicity in a factorised determinant** (which is how the
orchestrator had it by hand). Both engines compute the **nullity of the dynamical matrix at each root**,
and both compute the **transverse nullity** — the corank of `M` stacked on `kᵀ` — separately, so a mode
that is null for `M` but parallel to `k` is not miscounted as transverse.

| `D` | `ω² = 0` | `ω² = (μ_R/ρ_br)k²` | sum |
|---|---|---|---|
| 2 | nullity 1, transverse nullity 0 | nullity 1, transverse nullity 1 | 2 |
| **3** | nullity 1, transverse nullity 0 | **nullity 2, transverse nullity 2** | 3 |
| 4 | nullity 1, transverse nullity 0 | nullity 3, transverse nullity 3 | 4 |
| 5 | nullity 1, transverse nullity 0 | nullity 4, transverse nullity 4 | 5 |

⭐ Both engines, every cell. `n2_nullity`, `n2_rank`, `n3_stacked_rank`, `n3_transverse_nullity`,
`n4_nullity_difference`, `n7_basis_count` and `n7_basis_count_residual` are **`AGREE` on all 26 declared
per-root cells**; `q3_root_count`, which is declared once per package-and-dimension rather than per root,
is `AGREE` on all **13**. ⇒ this is the load-bearing comparison of the step and it is complete.

⛔⛔ **THE COUNT IS NOT THE CODIMENSION — and the sharp form of that is stronger than "wrong number".**
⭐ **The bulk never enters this computation at all.** There is no codimension anywhere in the
calculation; the nullity depends only on the brane's own dimension. ⇒ Codimension is not merely the wrong
answer (`4−3 = 1`), it is an **absent quantity**.

⇒ ⭐ **Light having exactly two polarisations is a statement that our space is three-dimensional.**

### ⭐⭐ And the count is not merely GENERIC — that was checked, and it fires on a control

A rank computed in symbols is the rank **away from measure-zero loci**. `§Q8` therefore solves, for each
matrix, the locus where the generic rank drops, asks whether that locus meets the allowed region
(`Σk² > 0`, `k` real, control parameters in range), and **recomputes the whole `N1`–`N7` set on every
stratum that does**.

| package | allowed strata | consequence |
|---|---|---|
| **MAIN**, all `D` | **none** | ⭐ `D−1` is not a generic-only statement here. There is **nowhere in the allowed region for it to fail.** |
| `XFORM_FULLGRAD`, `XFORM_DIVONLY`, `XFORM_SIGNFLIP`, `XCOEF_SCALE` | none | |
| **`XFORM_ANISO`** | **one**, `k₂ = k₃ = 0` at `D=3` (`k₂=k₃=k₄=0` at `D=4`) | ⛔ **the generic count is wrong there** |

⭐⭐ **The anisotropic control is the case the machinery exists for, and it fires.** Generically the
anisotropy **splits** the degenerate transverse set into two propagating roots, and on the stratum — `k`
along the anisotropy axis — the degeneracy is **restored**:

| | `D = 3` | `D = 4` |
|---|---|---|
| generic: roots | 3 | 3 |
| generic: transverse nullity, the two nonzero roots | **1** and 0 | **2** and 0 |
| on the stratum: roots | 2 | 2 |
| on the stratum: transverse nullity at the propagating root | **2** | **3** |

⇒ the generic count is short by exactly the piece the anisotropy split off, and it is short by a
different amount in each `D`. ⛔ **Neither figure generalises across `D`** — an earlier draft of this
record quoted the `D=3` number as though it did.

⚠ **This is a hand comparison, not a machine one.** The `§8` tag grammar **these engines were built
from** has no scope token for a stratum (`s10-as-built`; the repaired spec adds one, and neither engine
has been rebuilt to it), so the engines emit `WL_…_D3_STRATUM1_ROOT2_N3_TRANSVERSE_NULLITY` against
`PY_…_D3_Q8_STRATUM1_ROOT2_N3_TRANSVERSE_NULLITY`, and `checks_S10.yaml` declares **no cross-engine pair
for any stratum tag**. Paired by hand on `(package, D, stratum, root, quantity)`: **24 count-bearing keys
pair and every pair matches**; 8 do not pair, all from one suffix spelling (`N7_COUNT_RESIDUAL` against
`N7_BASIS_COUNT_RESIDUAL`, both `0`). ⇒ `reduction/measurements/q8_stratum_manual_comparison.py`.
⛔ **Do not report the stratum agreement as a harness result.** It is the orchestrator's.

## The equations and dimensions

```
L = ½ ρ_br (∂_t u)² − ½ μ_R (curl u)²          an energy density on a D-dimensional brane
(curl u)² ≡ ½ Σ_{i,j} (∂_i u_j − ∂_j u_i)²      the antisymmetric derivative, defined in every D
```

**Dimensions as a closed function of `D`** — new here, and it does retroactive work:

```
[ρ_br] = (−D,   0, 1)          at D=3:  (−3,  0, 1)  ✓ registry
[μ_R]  = (2−D, −2, 1)          at D=3:  (−1, −2, 1)  ✓ registry
```

⭐⭐ **This EXPLAINS an S9 loose end.** An S9 ablation showed that changing the assumed brane dimension
corrupts both constituent dimensions while the speed's dimension check stays green — filed then as a
blind spot in the audit. It is not a blind spot but an **identity**:

```
[μ_R] − [ρ_br] = (2−D) − (−D) = 2   in the length slot, for EVERY D
```

⇒ The speed's dimension **cannot** see the brane dimension. A weak-looking check is a property of the
physics.

### ⭐⭐ Both S10 engines agree with S9's and S11's on these two vectors

⚠ Measured 2026-08-08. Recorded here because the script that measured it served an abandoned track and
was deleted the same day; this paragraph is what survives it.

Five engines across S9, S10 and S11 emit a symbolic dimension for `ρ_br` and `μ_R` — S9-py, S9-wl,
S10-py, S10-wl, S11-py — and all five emit `(−D, 0, 1)` and `(2−D, −2, 1)`. ⭐ **Three carry zero
registry references in their source** (S9-py, S9-wl, **S10-wl**), so they cannot have read the declared
value: `grep -c "registry\|reduction"` returns `0`, `0`, `0` against `14` for S10-py and `18` for S11-py.

⇒ S10's own two engines are on **opposite sides** of that line — `S10-wl` derives it register-free,
`S10-py` does not. ⭐ That is the sharpest form the corroboration takes at this step.

⛔ **Scope.** ⚠ This corroborates the dimension **solve** — the derivative-multi-order extraction and the
linear system — ⛔ not the action it is applied to. The curl-only stiffness is authored identically into
every engine, so a defect in it moves all five together and is caught by none ⇒ the common-mode limit
under **WHAT THIS STEP STILL DOES NOT ESTABLISH** is untouched by this.
⚠ `B_comp` has **one** emitting engine and therefore no corroboration; ⛔ do not group it with these two.

## What's new — the introduction inventory

| item | class | why |
|---|---|---|
| `D_brane = 3` | **postulated** | what light *requires* of the medium; ⛔ not derived from anything ⇒ `R-S1-01` |
| the mode count | **derived** | nullity of the dynamical matrix; `D_brane − 1` transverse, 1 longitudinal |

⛔ **`n_pol` is deliberately NOT a registry row.** Established by test, ⛔ not preference:
`constraint_dimension` **raises** on a relation whose residual involves discrete symbols, and the schema
(`additionalProperties: false`) has no vocabulary for *"this discrete value derives from that one"*. ⇒ A
bare discrete row would make a **derived** count look like a **free choice** — worse than omitting it.

⭐ **`D_brane` closes a real hole.** S9's dimensions depend on it, and it was **declared nowhere**.

## Inputs, by locus

| input | from |
|---|---|
| `L`, `ρ_br`, `μ_R`, the curl-only form | **S9** (`steps/S9_light_requires_shear.md`), which takes them from **S8** |
| `u` is the in-plane displacement | this step's flagged identification, user-confirmed ⇒ `R-S8-02` |
| ⛔ nothing from the bulk | it does not enter |

## Regime

Linear response about an unstrained brane; small oscillations — half one's scope.

## Departure

⛔ **None at this step.** ⚠ But note what S10 does *not* deliver: `ω² = 0` means the longitudinal is
**non-propagating**, ⛔ **not absent**. Maxwell has no longitudinal photon at all. The slot still exists
here, and once the medium is allowed to be **compressible** that zero lifts. → **S11**.

⛔⛔ **And per the user (2026-08-02), that is NOT a defect to be fixed.** Exact Maxwell would be the
**failure**: Maxwell puts charge in by hand, so a model matching it exactly would have **no way to
physically anchor charge**. The extra mode is the anchor — it is what made the drum-head picture click.
→ `CHARTER.md#falsification-standard`, `V3_STEP_PLAN.md` § S11.

---

## ⭐⭐ PRIOR ART — and the honest reading of it

**This sector's algebra is MacCullagh's rotational aether (1839).** A medium whose stored energy depends
only on the *rotation* of the displacement — `L = ½ρ u̇² − ½μ(∇×u)²` — is exactly his construction, and
that it yields `D−1` transverse waves at one speed plus a longitudinal mode with **zero** restoring force
is textbook (Whittaker, *A History of the Theories of Aether and Electricity*, Ch. V). ⇒
`reference_prior_art_maccullagh`.

⭐ **Prior art is an oracle, never a premise.** Our result was computed independently and then checked
against his; ⛔ nothing here assumed it.

⛔⛔ **AND THE OVER-CLAIM TO REFUSE, because it is the most tempting sentence in the sector:**
**S10 linearises about `v₀ = 0`, which *is* MacCullagh's regime.** Agreement with him is therefore
agreement **inside his own domain of validity**, and it is ⛔ **not** evidence that this medium differs
from a nineteenth-century elastic aether. Anything that distinguishes the two has to come from somewhere
S10 does not look.

⭐ **What is genuinely ours here**, and it is very modest:

- **the bulk never enters** — he had no bulk, so *"the count cannot be the codimension"* is not a
  statement he could have made or needed.

⛔⛔ **Three things an earlier draft of this record listed here and should not have.**
- The **in-plane / out-of-plane split** is ⛔ **not ours.** The `h`-branon is a **branon**, and branons
  are filed as **KNOWN prior art** in this project's own register (Cembranos, Dobado & Maroto,
  *Phys. Rev. Lett.* **90** 241301, 2003) ⇒ `docs/medium_requirements_and_prior_art.md:142`. The in-plane-field /
  branon decomposition is the standard brane-world split. ⚠ A leg caught this against our own oracle —
  ⭐ **exactly what an oracle is for**, and a reminder that "I did not find it" is not novelty.
- The **`D`-sweep** is an **evidence exercise** over a standard algebraic identity, ⛔ not a result. It
  is worth having and it is not ours.
- Keeping the longitudinal zero **deliberately as the charge anchor** is recorded only as **NOT-FOUND**
  (`docs/medium_requirements_and_prior_art.md:151-152`), and NOT-FOUND is not originality.

⭐ **What may be residual is the ARCHITECTURE, not any piece of it** — curl-only stiffness on an ordered
thin phase with a shear-free bulk, giving confinement. ⚠ The prior-art register records that
**combination** as NOT-FOUND; ⛔ and `reference_prior_art_maccullagh` is explicit that NOT-FOUND is not
originality.

⚠ The three things that killed the elastic-aether programme — a longitudinal mode, a preferred frame, and
matter–aether coupling — are **not** answered at S10. S11 establishes that all three reduce to the S11b
interface law. ⛔ Until that is rebuilt and closed, this sector has no answer to MacCullagh's critics.

---

## ⭐⭐ VERIFICATION — rebuilt 2026-08-05/06

### Two engines, independently constructed

| | engine | emitted payloads |
|---|---|---|
| Mathematica | `mathematica/S10_brane_mode_spectrum_mathematica_audit.wl` | 2983 |
| SymPy | `scripts/S10_brane_mode_spectrum_sympy_audit.py` | 4227 |

**What is measured** is the inventory: of the emitted payloads, **296 (WL) and 863 (PY) are bare words**,
across 8 and 20 distinct values (`reduction/measurements/bare_word_payload_scan.py`). They fall into:
**CAS booleans** from homogeneity, domain and intersection tests · **symbol names** (`rho_br`, `mu_R`,
`D_brane`) · **route tokens** recording which solver branch was taken (`quadraticFormRoute` / `M_B`,
`factor_returned`, `solve_then_explicit_real_filter`) · **undecidability sentinels** (`Indeterminate`,
`undecided_under_joint_assumptions`, `undefined_zero_root_ratio`, `decidedEmpty`) · **scope
disclosures** (`codingConsistencyOnly` / `same_action_variational_identity`) · **premise echoes**
(`S_curl_supplied`) · **solver-determination words** (`exactlyDetermined` / `exactly_determined`) · and
**skip reasons** (`locus_conflicts_with_positive_wavevector_norm`, `rank_floor_empty_locus`,
`not_skipped_allowed_branch`).

⚠⚠ **That none of them asserts a physics outcome is the ORCHESTRATOR'S READING, ⛔ not an instrument
output** — the script inventories, it does not classify. ⚠ And the last three categories are the ones to
watch: `locus_conflicts_with_positive_wavevector_norm` (×11) is the *"this stratum is not allowed"*
determination **stated as a word**. ⭐ Clause 2 is not violated, because the operands behind each of
those determinations are emitted separately — ⛔ but a category list is not a proof, and an earlier draft
of this record both claimed more than the script measures and enumerated fewer categories than exist.

### ⭐⭐ The five controls, and what each one is actually able to test

All at `D = 3`, both engines agreeing on every cell:

| package | roots | `ω² = 0` | the propagating root | what it establishes |
|---|---|---|---|---|
| **`MAIN`** curl-only | 2 | nullity 1, transverse 0 | `+μ_R k²/ρ_br`, nullity 2, transverse 2 | the baseline |
| **`XFORM_FULLGRAD`** `Σ(∂_iu_j)²` | **1** | ⛔ **gone** | `+μ_R k²/ρ_br`, nullity **3**, transverse **2** | ⭐⭐ **FORM.** The transverse pair is **unchanged** — an ordinary elastic solid carries them at the same speed. What curl-only buys is that the **longitudinal does not propagate**, ⛔ not the transverse modes |
| **`XFORM_DIVONLY`** `(∂_iu_i)²` | 2 | nullity **2**, transverse **2** | nullity 1, transverse **0** | ⭐⭐ **FORM.** The roles **swap** — the transverse pair goes free and the longitudinal propagates |
| **`XFORM_SIGNFLIP`** | 2 | nullity 1, transverse 0 | **`−μ_R k²/ρ_br`**, nullity 2, transverse 2 | ⚠ **COEFFICIENT.** Nullities **identical to `MAIN`**; only the emitted **sign** `−1` distinguishes two growing instabilities from two waves |
| **`XFORM_ANISO`** `ρ₁ → s_ρ ρ_br` | **3** | nullity 1, transverse 0 | pair **splits**: transverse 1 and 0 | ⚠ **COEFFICIENT.** Breaks the degeneracy — and is the only package with an allowed stratum |
| **`XCOEF_SCALE`** `μ_R → s μ_R` | 2 | nullity 1, transverse 0 | `+s μ_R k²/ρ_br`, nullity 2, transverse 2 | ⚠ **COEFFICIENT.** Roots rescale, nullities unmoved |

⛔⛔ **Only `FULLGRAD` and `DIVONLY` are form controls**, despite the `XFORM_` prefix on four packages —
`SIGNFLIP` and `ANISO` vary a coefficient under the `curl` form. The battery classifies them by the
stiffness functional and not by the tag prefix (`ACCEPTANCE 18`:
`XFORM_SIGNFLIP=COEFFICIENT_ONLY`, `XFORM_ANISO=COEFFICIENT_ONLY`). ⇒
`reference_xform_prefix_is_not_a_form_control`.

⭐ **And the sharpened physics is in the first two rows.** S10's claim is often read as *"curl-only gives
two transverse modes"*. The `FULLGRAD` control refutes that reading: **an ordinary elastic solid gives the
same two, at the same speed.** ⇒ what the curl-only form contributes is the **absence of a propagating
longitudinal mode**, which is Maxwell's third demand — and it is a statement about the form of the
stiffness functional, ⇒ **`R-S8-01`**.

### ⭐⭐ The pre-registration, re-checked against the rebuilt engines

`steps/S10_PREREGISTERED_PREDICTION.md` was committed at `bc276485`, **before either script ran**, and
says plainly *"none is verified; if a script disagrees, the script wins."* Every prediction in it holds
against the rebuilt output, **including both flagged uncertainties**:

| predicted, target-blind | measured now |
|---|---|
| longitudinal nullity 1 in every `D`; transverse nullity `D−1`; sum `D` | ✓ all four `D` |
| form control: matrix becomes isotropic ⇒ **one** root, nullity `D`, the `ω²=0` root **disappears** | ✓ `XFORM_FULLGRAD` `D=3`: one root, nullity **3** |
| coefficient control: roots rescale, nullities unchanged | ✓ `XCOEF_SCALE` `D=3`: root `s·μ_R k²/ρ_br`, nullities 1 and 2 |
| *"least sure of"* — the `½` in the antisymmetric form giving exactly `\|∇×u\|²` at `D=3` | ✓ `§Q7` residual `0`, both engines |
| *"least sure of"* — whether `D=2` is degenerate **in kind** | ✓ it is not: nullity 1 / 1, transverse 0 / 1, sum 2 |

⚠ **What this is and is not.** The predictions were the orchestrator's hand derivation from the walk, so
confirming them says the hand algebra and the engines agree. ⛔ It is **not** an independent test of the
`§Q7` factor, because the same `½` convention sits in the shared spec that both engines read.

⛔⛔ **And the fourth row was not honestly testable until this rebuild.** Before `bad20207` the SymPy
engine's `§Q7` computed `curl − curl`, so its residual was `0` **for every package by construction** —
the "least sure of" item was being confirmed by a check that could not have reported anything else. ⇒
⭐ **only one engine tested that prediction until now, and the record said two.**

### What the harness measured

`reduction/engine_output_checks.py --config reduction/checks_S10.yaml` over both committed outputs.
⚠ **The block below is ABRIDGED and re-ordered** — every value in it is literal, but fields are dropped
(`no_comparison`, `not_applicable`, `unassessable`, the coverage formula string) and lines are grouped
two-per-line. ⛔ It is not a transcript; re-run the comparator for one. ⚠ And note the process **exits
`2`** on the committed inputs, printing a long `OPERATIONAL_FAILURE:` line to stderr — that is the
design (operational findings, ⛔ not physics failures), and it is the first thing a re-runner sees:

```
CROSS_ENGINE: agree=517 cardinality_invalid=55 disagree=28 empty=18 naming_mismatch=13 shape_mismatch=59
CROSS_ENGINE_COVERAGE: 545/690
CONTROL_RESPONSE[wl]: compared=6 responsive=6 invariant=0    CONTROL_RESPONSE[py]: 6/6 responsive
TAG_PARITY[wl]: cells=13 gaps=0                              TAG_PARITY[py]: cells=13 gaps=0
DIMENSIONS[wl]: total=2983 compared=430  homogeneous=428  unwalked=3
DIMENSIONS[py]: total=4227 compared=822  homogeneous=802  unwalked=0
LAGRANGIAN_ACTION_COVERAGE 13/13        EULER_LAGRANGE_ACTION_COVERAGE 13/13
```

⭐⭐ **The action is compared, and this is the first time it was.** All 13 Lagrangian and 13
Euler–Lagrange rows `AGREE`. Verified twice independently by review legs: one rebuilt the comparator with
a stricter key, one recovered each package's stiffness density from the emitted Lagrangian and divided it
out. **Single-difference controls** — one derivative order, one variable, one function, one *member* of a
differentiated function's argument list, one coefficient, one form — each still reach `DISAGREE`, so the
agreements are content and not looseness. ⭐ And the discrimination is the right shape: a difference in
argument **order** still reaches `AGREE`, because order is an emission convention while membership is
content. ⚠ That distinction was a live false-agreement channel until `bad20207` — the canonical
derivative atom dropped the argument list entirely, so fields with different coordinate dependence
compared equal. Reverting the fix puts the arity attack back to 26/26 `AGREE`.

### The 28 remaining disagreements — ⭐ every one is now accounted for, and NONE is a computational conflict

| class | rows | what it is |
|---|---|---|
| `*_q5_scale_ratio` | 11 | **Both engines agree.** At `ω²=0` the scaling ratio is `0/0` and both emit an *undefined* sentinel — `undefinedRatio(0)` against `undefined_zero_root_ratio`. Two spellings of the same determination; `§8` never pinned sentinel names |
| `*_n6_nullspace_basis` | 11 | A nullspace **basis** is not canonical. The engines return different bases; the comparator compares the literal list. ⛔ **Nothing computes whether they span the same space.** The *dimension* — which is the physics — agrees on all 26 cells |
| `*_q3_sign` | 3 | ⛔ **Not a sign conflict** — an undischarged premise. See below |
| `*_q6_unknown_coefficient_dimension_count` | 3 | ⛔ **Not a computation** — a gap in the shared spec. See below |

⭐ Plus the `R4` registry residual, nonzero for a known reason: it compares a squared speed against a
speed.

⭐ **In no row that REACHED A VERDICT do the two engines compute different physics.**

⛔⛔ **And that sentence has to carry its denominator.** Of the 690 declared rows, **545 reach a verdict**
and **145 do not**: `shape_mismatch=59`, `cardinality_invalid=55`, `empty=18`, `naming_mismatch=13`. ⇒
**a fifth of the declared comparison never happens**, and in those rows the two engines could be computing
anything at all without this instrument noticing. ⛔ **Never quote the agreement count without the
545/690.**

⚠ **The gaps are not scattered — they are eight families, and every one is a shape or naming difference
rather than a hard case:**

| family | rows | |
|---|---|---|
| `q8_rank_drop_locus` · `q8_rank_drop_minors` · `q8_root_coincidence_loci` | 26 each | the locus algebra: solution **sets** against solution **lists** |
| `q8_allowed_strata` · `q3_root_coincidence_loci` · `q6_dimension_equations` · `premise_period_average` | 13 each | ditto, one per package-and-dimension |
| `q2_downstream_route` | 13 | the bare-symbol route token — see below |
| `run_pairs`, `skipped_pairs` | 1 each | run bookkeeping |

⇒ ⭐ **the whole `§Q8` locus layer — 91 rows — is emitted by both engines and compared by neither.** That
is where `S10`'s generic-versus-stratum reasoning lives, and it is exactly the material the record above
had to pair **by hand**.

#### The three `q6_unknown_coefficient_dimension_count` rows are a SPEC CONVENTION GAP

`py=6` against `wl=9`, on exactly the three packages that carry a scale factor. Measured across all
thirteen package-and-dimension pairs:

- **SymPy emits `6` for every one of the thirteen** — the count is a constant and ⛔ **does not respond
  to the action at all**; it counts the two coefficients that are not declared dimensionless, times three
  components;
- **Mathematica emits `6` for the four packages without a scale factor and `9` for the three with one** —
  it counts the scale factor as an unknown constrained by the supplied declaration.

⇒ ⭐ **both are internally consistent, and both conform to the shared spec**, because the spec never says
whether a symbol declared dimensionless is still an unknown. ⭐ **The cross-engine comparison found a hole
in the artifact both engines read** — which is the one defect class comparison is supposed to be blind to,
caught here only because the two engines resolved the ambiguity differently.

⚠ **And two legs, independently, computed the sharper form.** With **no** dimensionless declaration
made, the coefficient-dimension system has **nullity 0** for `XFORM_ANISO` — so `[s_ρ]` is **determined by
the action** — and **nullity 3** for `XCOEF_SCALE`, so `[s]` genuinely is a premise. ⇒ the as-built spec
supplied as an unfalsifiable premise a fact the action determines, and `§7`'s stated reason for the
declaration was true of one symbol and false of the other.

⭐ **The spec now states the counting rule in both directions**: a symbol **declared** dimensionless is
not a `Q6` unknown; one whose dimension the action **determines** is, and stays one. ⇒ `s` is excluded
(so `XCOEF_SCALE` is **6**) and `s_ρ` is included (so `XFORM_ANISO` is **9**), and each engine's present
answer is non-conforming on one of the two.
⚠⚠ A first attempt at that repair **wrote the solved value into the spec** and left the `XFORM_ANISO`
case decidable both ways; a leg measured both and it was corrected again. ⛔ **The pinned engines predate
all of it**, so the three rows still read `DISAGREE` and close only when the engines are rebuilt.

#### ⛔⛔ The three `q3_sign` rows are an UNDISCHARGED PREMISE, not a disagreement

⚠ In **none** of the three does one engine assert a sign opposite to the other's. Taking each disputed
root as the engines emitted it and asking the CAS what the declared premises entail
(`reduction/measurements/q3_sign_adjudication.py`): every one is a sum of `k_i²` terms whose coefficients are **each
declared positive**, over a **declared-positive** denominator. What differs is only whether an engine
*discharged* the premise `Σk² > 0`.

| row | WL | PY |
|---|---|---|
| `main_d5_root2` | `Sign[k1²+…+k5²]` — declined | `1` — discharged |
| `xform_aniso_d3_root3` | `1` — discharged | `undecided_under_joint_assumptions` — declined |
| `xform_aniso_d4_root3` | `Sign[k2²+k3²+k4²+k1²·sRho]` — declined | `undecided…` — declined |

⭐ WL discharges **24 of its 26** sign rows and leaves exactly the two **largest expressions** symbolic —
a simplifier-depth artefact, ⛔ not a physics boundary.

⛔⛔ **And the deeper finding: `Σk² > 0` is declared by both engines and reliably discharged by
neither.** It is a **relation among symbols**, not a property of one, so SymPy's assumption system cannot
carry it — `ask(Q.positive(root))` returns `None` on all three roots even though each is manifestly a
positive-weighted sum of squares. ⇒ ⛔ **neither the shared spec nor the harness requires an engine to
show that it discharged a declared premise before reporting a sign**, so each sign is whatever that
engine's simplifier happened to reach. ⭐ A new defect-register entry, and it generalises past `q3_sign`
to every emitted decision that depends on a relational premise.

### ⛔⛔ `§Q7` — the check that could not fail, and what fixing it did and did not buy

⚠⚠ **The SymPy engine's `§Q7` computed `curl − curl` for every package.** The action was built from each
package's own stiffness form, but `§Q7` hardcoded `"curl"`, so the emitted residual asked whether the
curl density equals the curl norm — and answered it six times. ⇒ **the previous version of this record's
claim that "the reduction to `|∇×u|²` at `D=3` is COMPUTED, not assumed" was, on that engine, worth
nothing.** Repaired at `bad20207`: the engine now emits the density its own action used, established by
mutation and from the output alone.

⭐ **What `§Q7` now reads**, both engines, `D=3`: residual `0` for the four curl-form packages, nonzero
for `XFORM_FULLGRAD` and `XFORM_DIVONLY`.

⛔⛔ **AND IT COMPARES THE FORM ONLY.** `§Q7` drops the stiffness sign and the coefficient. Recovering
each package's coefficient from its **own emitted Lagrangian** — dropping the time-derivative terms,
substituting the gradient symbols, and dividing by the emitted density (residual `0` in all twelve
engine × package cases, `reduction/measurements/stiffness_coefficient_recovery.py`):

| package | recovered coefficient, both engines | what the engine *reports* as `Q6_STIFFNESS_COEFFICIENTS` |
|---|---|---|
| `MAIN` | `−μ_R/2` | `μ_R` |
| `XFORM_FULLGRAD` | `−μ_R/2` | `μ_R` |
| `XFORM_DIVONLY` | `−μ_R/2` | `μ_R` |
| **`XFORM_SIGNFLIP`** | **`+μ_R/2`** | **`μ_R`** |
| `XFORM_ANISO` | `−μ_R/2` | `μ_R` |
| `XCOEF_SCALE` | `−s·μ_R/2` | `s·μ_R` |

⇒ ⛔ **`XFORM_SIGNFLIP`'s `§Q7` payload is byte-identical to `MAIN`'s** while its action carries the
opposite sign, and the reported coefficient tag drops **both** the sign and the ½ — on **both** engines,
which share the convention. ⭐ **A zero `§Q7` residual is a claim about the FORM, never about the TERM.**

⚠⚠ **And `§Q7`'s reference side is TYPED, not computed, in both engines.** The repaired spec requires
`c_i = Σ_{j,k} ε_ijk g_jk` to be reached **by computation from the Levi-Civita definition**, precisely so
that one side of the comparison is something the engine derived rather than transcribed. ⛔ **Neither
as-built engine does this** — both write the three components out by hand
(`.wl:1315-1318`, `.py:1544-1550`) and take the dot product of the typed list. ⇒ `§Q7` compares an object
the action produced against an object the author typed, and the only thing standing behind the typed side
is that **two authors typed it independently and agree**. ⭐ That is worth something; ⛔ it is not what the
spec now asks for, and a rebuilt engine must close it.

⚠ **The sign control is caught, but by exactly one emitted object.** `XFORM_SIGNFLIP` has the same
nullities as `MAIN` at every root; what distinguishes them is `ROOT2_Q3_SIGN = −1` against `+1`, i.e.
`ω² = −μ_R k²/ρ_br` — two exponentially growing modes, not two waves. ⇒ **the mode count alone cannot
tell a wave from an instability.** ⭐ Emitting the sign per root is what makes the sign control able to
fail, and it is the same lesson S9 learned.

### ⭐⭐ The gradient-symbol declaration — and the exact size of what it rests on

The two engines spell the nine `D=3` gradient symbols differently: `g{r}x{c}` against `g{r}{c}`. Until
2026-08-06 **eight of the twelve** declared `§Q7` rows reached no verdict: the comparator fell back to a
per-row bijection search and reported `NAMING_MISMATCH`, which is ⛔ **not a comparison**. ⚠ The other
four reached `AGREE` — their payload is literally `0` on both sides, so there was nothing to relabel and
no search ran. ⇒ the operand is `ACCEPTANCE 19`'s `ABSENT=[AGREE=4,NAMING_MISMATCH=8]`, printed on every
battery run.

The map is now declared in `reduction/checks_S10.yaml`, **sourced from the two engine definitions**
(`.wl:148-152, 216, 1304-1311`; `.py:309, 1536, 1538-1542`): both index the symbol `(coordinate, field)`,
so `g{r}x{c} = g{r}{c} = ∂u_c/∂x_r`, the identity on the index pair.

⛔⛔ **THE PAYLOADS COULD NOT HAVE TOLD US THIS, and the measurement is exact.** Over all `9! = 362 880`
bijections of the nine symbols, **288 leave every emitted `§Q7` density invariant**
(`reduction/measurements/q7_payload_invariance_group.py`; reproduced independently by a review leg). The
group is

```
S3 on the three diagonal symbols  ×  S3 on the three off-diagonal PAIRS  ×  2³ within-pair swaps  =  288
```

⇒ ⭐ **the payloads fix only which symbols are diagonal and how the off-diagonal ones pair.** The
transpose is one member of that set; so is every permutation of `g11, g22, g33`. Under the transpose the
harness's **entire counter line is byte-identical**, not merely its `§Q7` tally. ⇒ ⛔ no bijection search,
no ablation and no verdict-keyed assertion can narrow the map. ⭐ **That is why it is sourced from the
engine definitions**, and the config says so where the declaration is made.

⭐ **What the declaration buys is not the agreement** — the per-row search had already reconciled those
eight — **but that ONE map is applied to EVERY row**, so a relabelling that is wrong in one row alone is
no longer absorbed by a rename no other row uses.

⚠ **"Per-row bijection" is narrower than it sounds, and the rule below depends on which.** The comparator
reconciles only symbols appearing on **one** side (`engine_output_checks.py:1128-1136`); if no symbol is
left over, no bijection exists and the row reaches `DISAGREE` rather than `NAMING_MISMATCH`.

⛔⛔ **A single-row derangement is governed by stabilisers, not by a status attached to that row.**
A derangement drawn from the **288-element joint stabiliser** produces no signal on **any** row. Outside
it, a nonzero payload moves unless the permutation lies in that payload's own stabiliser: when it
moves, the row reaches `DISAGREE` if the corrupted payload cannot be re-reconciled by a per-row
bijection, and `NAMING_MISMATCH` if it can.

⚠ **Five of the twelve rows are invisible to every single-row derangement, in two distinct ways.**
Exactly one payload-bearing row, **`xform_fullgrad_d3_q7_stiffness`**, has stabiliser **`9!`**: its
nonzero `Σ g_ij²` payload is invariant under every relabelling. The other four are `_q7_difference`
rows whose payload is literally `0`, hence vacuously `9!`-invariant. ⇒ **No payload check can police
`xform_fullgrad_d3_q7_stiffness`, which is why the unit test pins its nine renames by hand.**

⚠⚠ **This paragraph was wrong three times because each draft generalised one status from a single
measured branch.**

**Measured effect, whole declared row set:** exactly **8 rows** change status, all `§Q7`, all
`NAMING_MISMATCH → AGREE`; a leg confirmed *touched-but-did-not-move* and *moved-but-not-touched* are both
empty. `agree 509→517`, `naming_mismatch 21→13`, coverage `537→545`, **`disagree` unchanged at 28** —
and with five injected physical differences (coefficient, sign, dropped term, changed form) every one
still reaches `DISAGREE`. ⛔ Not one row moved from `DISAGREE` to `AGREE`.

⚠⚠ **And `AGREE = 12` is twelve rows but only FIVE distinct informative comparisons.** The twelve are six
stiffness rows and six residual rows. `MAIN`, `XFORM_SIGNFLIP`, `XFORM_ANISO` and `XCOEF_SCALE` emit
**byte-identical** `§Q7` stiffness payloads — `§Q7` drops the sign and the coefficient — so the six
stiffness rows carry only **3 distinct values** (curl, full-gradient, divergence-only). And four of the
six residual rows are literally `0 == 0`, leaving **2** non-vacuous. ⇒ `3 + 2 = 5` distinct informative
comparisons, replicated across twelve rows.
⛔ **Quote the 12 only with that beside it.** Likewise of the 18 emitted `D=3` pairs, **13** carry
information about the map; four are `0` on both sides and `XFORM_FULLGRAD`'s stiffness pair is
`9!`-invariant.

**Able to fail:** `harness_ablation.py` ACCEPTANCE 19 — deranging or removing the map changes the tally;
transposing it does not, and the battery prints that rather than asserting it. ⚠ A leg then broke both
gates and found two mutations confined to the **diagonal** entries that passed both — a wrong map
(`g22 ↔ g33`) and a silently dropped entry. ⭐ Closed in the same session: the unit test now pins the
complete nine-rename set on the one row whose payload carries all nine, and the battery now requires
**all twelve `§Q7` rows to be `AGREE`**; merely forbidding `NAMING_MISMATCH` would admit a wrong map that
produces `DISAGREE`.

---

## ⛔ WHAT THIS STEP STILL DOES NOT ESTABLISH

- ⛔ **That the in-plane and out-of-plane sectors decouple.** The count at `D=3` is 2 rather than 3
  *because* `u` was identified as the in-plane displacement alone. Nothing computes that. ⇒ `R-S8-02`.
- ⛔ **`D_brane = 3`** is postulated. The bulk does not enter, so nothing here can derive it ⇒ `R-S1-01`.
- ⛔ **The absence of a propagating longitudinal wave is ASSUMED, not derived.** The curl-only action sets
  the longitudinal restoring stiffness to **zero by construction**. It removes the restoring *force*,
  ⛔ **not the degree of freedom** — which is why the longitudinal test is emitted at every root, and why
  S11 exists.
- ⛔ **The nullspace bases are never compared for span.** Eleven rows report `DISAGREE` on representation
  alone and nothing computes the question that matters there.
- ⛔⛔ **Thirteen `q2_downstream_route` rows are a check that CANNOT FAIL**, and it is the same family as
  the `§Q7` defect this rebuild removed. Both payloads are **bare symbols** (`quadraticFormRoute` against
  `M_B`), and the comparator's bijection search reconciles any single symbol with any other. Substituting
  the SymPy side and re-running the comparator (`reduction/measurements/route_row_information_content.py`): a token
  naming **the other route**, a token naming **no route** (`banana`), a token asserting **the routes
  disagreed**, and a token asserting **`FAILED`** all return the identical `NAMING_MISMATCH` — the same
  verdict the committed payloads return. ⇒ ⭐ **the row distinguishes only literal equality from
  everything else, and the harness does not count that as a verdict at all.** It is carried here as a
  named non-check, ⛔ not as evidence about the two routes.
- ⛔ **The stratum re-run is not machine-compared** — the `§8` grammar these engines were built from
  (`s10-as-built`) has no stratum scope token; the repaired `directives/S10_SHARED_PHYSICS.md` adds
  `STRATUM<s>`. The 24-key agreement above is a hand pairing.
- ⚠ **Cross-engine agreement did not improve when `§Q7` was repaired.** Across `bad20207`,
  `DISAGREE 32→28` and `NAMING_MISMATCH 17→21`, coverage numerator `541→537`: four rows moved from
  *compared and wrong* to *not compared*. The declaration then moved eight rows the other way. ⛔ **Do not
  write the repair up as improving agreement.** Net across both: `disagree 32→28`,
  `naming_mismatch 17→13`, coverage `541→545`.
  ⚠ **Provenance of that first triple:** it is **cited from `bad20207`'s commit message**, ⛔ not
  re-measured here — the pre-repair SymPy output it describes was overwritten by the repair. The second
  leg is measured, twice, this session. ⇒ a reader wanting the pre-repair numbers independently must
  re-run the pinned engines at `s10-as-built`.
- ⚠⚠ **Agreement rests on declared spellings, and the margin is large.** The harness computes each row's
  verdict a second time with **no** declared naming applied: `NAMING_EFFECT: legacy_before_agree=299
  declared_after_agree=517`. ⇒ **218 of the 517 agreements exist only because a spelling was declared.**
  ⭐ Every declaration is printed on its verdict line and every one carries a reason; ⛔ but the number to
  quote for agreement *between two engines left to their own orthography* is **299**, not 517.

  **Which declaration carries what**, by removing each one and re-running
  (`reduction/measurements/declaration_load_ablation.py`):

  | declaration | rows it moves | of which would be `DISAGREE` without it |
  |---|---|---|
  | `s_rho ↔ sRho` | **30** | **18** |
  | `s ↔ coefficientScale` | **14** | **8** |
  | `D ↔ braneDimension` | 13 | 0 |
  | the nine gradient symbols, jointly | 8 | 0 |
  | `rho_z`, `mu_F`, `mu_G`, `lambda_rho`, `lambda_mu`, `lambda_scale` | **0 each** | — |

  ⛔⛔ **This corrects a figure recorded at `bad20207`**, which attributed the 30-row declaration to
  `lambda_scale`. On this config `lambda_scale` moves **nothing**; the 30 rows belong to `s_rho`. The
  other two figures in that note (`D` 13, `s` 14 with 8 from `DISAGREE`) are confirmed exactly.
  ⚠⚠ ⇒ **`XFORM_ANISO`'s and `XCOEF_SCALE`'s agreement is where the declarations actually bite** — 26
  rows would read `DISAGREE` without those two.
  ⚠ And `s ↔ coefficientScale` is the **weakest** declaration in the config: unlike the others it is
  ⛔ **not** a snake-to-camel transliteration but two engines choosing genuinely different names. What
  justifies it is that both bind the symbol as the dimensionless multiplier of `μ_R` in the same
  package's emitted stiffness coefficient (`{coefficientScale*muR}` against `(mu_R*s,)`) — ⭐ a reason,
  ⛔ not a mechanical rule.
- ⚠ **`XFORM_SIGNFLIP` and `XFORM_ANISO` are not form controls despite the `XFORM_` prefix.** Both vary a
  coefficient under the `curl` form; only `FULLGRAD` and `DIVONLY` change the stiffness functional.
  ⇒ `reference_xform_prefix_is_not_a_form_control`.
- ⚠ **Twenty dimension payloads on the PY side and two on the WL side are non-homogeneous**, and all of
  them sit on the `XFORM_ANISO` stratum, where an explicit numeric wavevector has been substituted — so a
  `k²` factor's dimension is gone by construction. Both engines report the same thing. Three WL payloads
  are **unwalked** (`Q6_DIMENSION_PREMISES`, relation with non-scalar sides).
- ⚠ **Two instrument hazards, both carried deliberately.** (i) The comparator's own bijection-verification
  line is **untested** — deleting it leaves the suite green. (ii) `OpaqueDerivative` subclasses
  `sp.Symbol`, which caches on the rendered name, so rebuilding one **returns the same object and
  re-stamps its attributes in place**. A leg confirmed the hazard is **real** (`rebuilt is original`) and
  **not reached**: instrumenting the branch over the real 690-row run gives **0 firings** across 61
  distinct atoms and 595 occurrences, and none of the 17 declared `wl` renames collides with any
  derivative's function name or differentiation variable. ⚠ **New datum:** reachability is
  **LRU-cache dependent** — rebuilding immediately mutates in place, rebuilding after churning 20 000
  symbols does not. ⇒ whether it fires is a function of unrelated allocation volume, so ⛔ **it must be
  fixed before any field or coordinate spelling is ever declared.**
- ⚠ **The spec both engines read was repaired *after* they were built.** The as-built text is pinned at
  `s10-as-built`; the repaired `directives/S10_SHARED_PHYSICS.md` describes what a rebuild should compute,
  ⛔ not what these engines were told.
- ⚠ Limits taken: linearisation about `v₀ = 0` · sharp zero-width sheet · no dissipation ·
  frequency-independent moduli · continuum limit · amplitude → 0.

## ⚠ A correction to S9's record, found here

S9's record cites its dimension layer as **`1219/1219 homogeneous`**. Re-run under the rebuilt harness
that figure is **312 comparisons** on the WL side (`compared=312 homogeneous=312` of `total_tags=1559`,
with `no_comparison=1228`) and **329** on the PY side. ⛔ The cited figure counted tags reached, not
comparisons made, and overstated the check by about a factor of four. S9 is a signed-off step; the
correction is applied to its record rather than left here.

## Registry additions (executed at this step)

`Q.brane.D_brane` — the brane's spatial dimension, `kind: discrete-choice`, `value: 3`, dimensionless,
**postulated**. ⛔ No `n_pol` row, for the reason above.

⚠ **One provenance dent, found while writing this record.** That row's `source_locus`
(`reduction/quantities.yaml:177`) points at `steps/S9_light_requires_shear.md:68-83` — a passage about the
two-sided shear requirement, ⛔ not about the brane's dimension. `D_brane` is introduced **here**, at S10,
which is what this record's introduction inventory says. ⇒ the registry entry cites the step **before**
the one that introduces it, which is the same provenance inversion the `{ρ_br, μ_R}` ordering correction
was made to prevent. ⛔ **Not fixed here** — a registry edit is out of this step's scope and would move a
gate; recorded so the next registry pass finds it.
