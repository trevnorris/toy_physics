# ledger_stage032_electric_sign_bc_ensembles_landing

## Status

**Part IV — Charge. IV-3 (build-order 032; the R1-LANDING stage of the 4-stage Part-IV split,
user decision 2026-07-22).** Reshape of the §4 landing table + the Q-BC / Q-AMEND / Q-COMBINE / Q-MAG
blocks of the puncture-deflection source build
`software/em_charge_attribute/puncture_deflection_electric_sign_result.md` (+ its `check.{py,wl}` /
`independent_recompute.py`). Stage 031 earned the target-blind far-field STRUCTURE; **stage 032 owns the four
boundary-condition ensembles that fix the force SIGN, shows the committed bare model does not select among
them, and lands the honest terminal token**:

- **`R1_REQUIRED(bc_selection)`** — the far-field two-body force `F_X = s₁s₂ A_X/(4πR²)` is target-blind
  EARNED in FORM (031) but its SIGN is not invariant across the admissible mouth boundary-condition
  ensembles. On a nondegenerate response (`m_gg>0`, `S_gg>0`, nonzero magnitudes): **V (Dirichlet) repels,
  J (fixed-source) attracts, M (fixed-monopole) is indefinite, MIXED spans**. The committed bare model
  (G0 + U2-unresolved `𝔅` + `S_hold` scoped to `r_B−½`) does not select the class ⇒ `outcome_not_invariant`
  ⇒ the terminal R1. Co-blocker **`R1_REQUIRED(magnitude)`** (`c_a`,`c_ξ` core normalizations unbounded at
  tier-A). Resolver = the SIM-DEFERRED nonlinear parent-throat boundary functional / barrier / map `s↦h_A`
  (a sibling of gravity's R10/R30/R33). `internal_inconsistency = none` (charge coexists with gravity;
  owned + verified here via Q-AMEND consistency).

**⚠ CRITICAL TYPING (do not overstate).** The electric SIGN is `R1_REQUIRED(bc_selection)` — **neither EARNED
nor CALIBRATED**. What is **DECIDED** (not R1): the four coefficient formulas `A_V/A_J/A_M/A_MIXED` *given a
class*, and the coupling `g=g_χh` (consumed from 031). What is **R1** (deferred): the class SELECTION among
{V,M,J,MIXED} and the ensemble data `{φ,q,j,λ}`. The strict signs (`A_V>0`, `A_J<0`, the strict MIXED
regimes) are CONDITIONAL on a nondegenerate response (`m_gg>0` at the `z_g≠0` witness, `S_gg>0`, nonzero
magnitudes); weak (`≥0`/`≤0`) in general; at `z_g=0` every `A_X` vanishes. The **magnitude**
(`R1_REQUIRED(magnitude)`, `c_a`/`c_ξ` unbounded at tier-A) and the **variant / MIXED-parameter /
REPLACE-vs-ADD** selections are SEPARATE downstream R1 debts — not part of the sign-class R1, and not
resolved by picking a class.

Ledger-local landing label (surviving-solution framing): the committed bare model **does not determine the
boundary operator**, hence does not pick the electric sign.

Charge = the static ±w throat. Stage 030 paid for the static substrate; stage 031 dressed it with the
target-blind puncture-deflection mechanism and emitted the neutral shell `A=m_gg·C`. Stage 032 turns that
neutral shell into the honest sign landing: it selects the class numerator `C` per BC class, runs the sealed
exhaustive landing ladder, verifies internal consistency, re-homes the charge magnitude, and lands the R1.

## Purpose

Stage 031 stopped at the EARNED, target-blind structure: the field identity `ξ_w=ℓh`, the reconstructed
bare-mouth source (`Q_χ=s`), the exterior `h(r)=h_A a/r`, the full completed-square response matrix `m`, and
the neutral far-field shell `A=m_gg·C` with `[A]=EL`, the `1/R²` falloff, and the `s₁s₂` product — all
target-blind, picking NO sign. It also EXPORTED the named fact `NONZERO_HA_REQUIRES_CORE_HOLDER`. Stage 032
is the R1-LANDING half of the puncture-deflection gate: it (i) selects the class-dependent numerator `C` per
admissible mouth BC class, giving the four DECIDED conditional coefficients `A_V/A_J/A_M/A_MIXED`; (ii)
characterizes the MIXED three-regime admissibility; (iii) verifies `internal_inconsistency=none` (Q-AMEND);
(iv) runs the sealed 23040-cell §4 landing ladder — recomputed from typed upstream facts — to the honest
terminal `R1_REQUIRED(bc_selection)`; (v) re-homes `Q_E`/charge-magnitude as the co-blocker
`R1_REQUIRED(magnitude)`; and (vi) certifies the complete 44-entry source-to-stage predicate manifest. It
picks NO sign: WHICH BC class holds is the R1, and its resolver is sim-deferred.

Consumes (by CITATION, NOT re-derived): from stage 031 (`THROAT_H_SOURCE_1_OVER_R2`) the bundle
`{m, m_gg=B_eff z_g²/D, det m=z_g²/D, S_gg, A=m_gg·C, 1/R², NONZERO_HA_REQUIRES_CORE_HOLDER}` (engines print
`CONSUMES_STAGE031=…`); from stage 030 the block `{D, D*=7/4}` transitively (engines print
`CONSUMES_STAGE030_TRANSITIVELY=…`). None of that mechanism is re-derived here — the ownership firewall
(`PASS_SCOPE_FIREWALL`) rejects any response-matrix / `1/R²` / mouth-mechanism re-derivation into 032.

## 1. The four DECIDED conditional ensembles (consuming 031's `m_gg`, `S_gg`)

The neutral shell `U = s₁s₂ A/(4πR)`, `F = s₁s₂ A/(4πR²)`, `A = m_gg·C` (031) leaves the `[E²]` numerator `C`
class-dependent. Selecting `C` per admissible mouth boundary-condition class gives four coefficients `A_X`,
each **reshaped through 031's two-mouth response without reconstructing the response matrix or its `1/R²`
mechanism**. All four are DECIDED *given a class* (NOT R1).

**V — Dirichlet / imposed value (`PASS_A_V_FORMULA`).** Reshaping the source V-ensemble through the
symmetric self/cross response `[[S_gg, ε m_gg],[ε m_gg, S_gg]]` on the held datum `(s₁φ, s₂φ)` and taking the
`O(ε)` cross-coupling of `Ω_V = −½ yᵀ M⁻¹ y` (divided by `s₁s₂`) yields the pair-coupling amplitude

```text
A_V = m_gg φ² / S_gg²          (V/DIRICHLET_VALUE; [A_V]=EL; repel at the witness).
```

**M — fixed monopole / imposed conormal (`PASS_A_M_FORMULA`).** The same-field source/reservoir work
`m_gg(q+g)² − 2g·m_gg(q+g)` factors to

```text
A_M = m_gg (q² − g²)           (M/FIXED_MONOPOLE; INDEFINITE — sign follows sgn(q²−g²)).
```

**J — fixed source (`PASS_A_J_FORMULA`).** The fixed-source functional `m_gg(j+g)² − 2(j+g)·m_gg(j+g)` factors
to

```text
A_J = −m_gg (j + g)²           (J/FIXED_SOURCE; attract at the witness).
```

**MIXED — imposed Robin relation (`PASS_A_MIXED_FORMULA`).** With the interpolation weight `λ`,
`m_gg(q+g)² − 2(g+λq)·m_gg(q+g)` factors to

```text
A_MIXED = m_gg [ (1−2λ) q² − 2λ q g − g² ]     (MIXED; [λ]=1; spans — see §2).
```

`COEFFICIENTS = {V: A_V, J: A_J, M: A_M, MIXED: A_MIXED}`. All four consume 031's `m_gg` (and `A_V` the
self-response `S_gg`); `g = g_χh` is the DECIDED committed coupling of 031 (`g_χh=J_m/ℓ ≠ 0`), NOT the sign R1.

### Weak vs strict sign typing (CONDITIONAL — `PASS_WEAK_SIGNS_GENERAL`, `PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS`, `PASS_DEGENERATE_Z_G_ZERO`)

031 earns only `m_gg≥0`. Accordingly:

- **Weak, general:** `A_V ≥ 0` (`m_gg≥0`, `S_gg>0`, `φ` real) and `A_J ≤ 0` (`m_gg≥0`) hold as weak
  certificates over the whole admissible range — verified by symbolic `ask`/`Reduce`, NOT as strict signs.
- **Strict, witness-only:** at the labelled NON-DEGENERATE witness (`z_g≠0` ⇒ `m_gg = B_eff z_g²/D > 0`,
  `S_gg>0`, `φ≠0`, `j+g≠0`; concrete witness `m_gg=8/7`, `S_gg=2`, `φ=3`, `j=2`, `g=1`) one has `A_V>0`
  (repel) and `A_J<0` (attract). The `strict_sign_context` predicate must be ACTIVE for these to be asserted.
- **Degenerate:** at `z_g=0` (⇒ `m_gg=0`) EVERY `A_X` vanishes and `strict_sign_context` is INACTIVE — the
  strict asserts are NOT raised (a degenerate response must not be typed strict). This is the tooth that keeps
  the sign honestly conditional rather than unconditional.

Type register (`PASS_DECIDED_CONDITIONAL_TYPING`):
`coefficients = DECIDED_given_class`; `bc_selection = R1`; `magnitude = R1`; `mixed_parameters = R1`;
`variant = R1`. Mis-typing the class selection as `DECIDED` fires the tooth.

## 2. MIXED three-regime admissibility (`PASS_MIXED_ENDPOINTS`, `PASS_MIXED_INTERIOR_ZERO`, `PASS_MIXED_REGIME_*`, `PASS_MIXED_FULL_DOMAIN_SPANS`)

`A_MIXED` interpolates between the two realized endpoints and carries one algebraic zero:

```text
λ = 0  → A_MIXED = m_gg(q²−g²) = A_M          (endpoint = the monopole coefficient),
λ = 1  → A_MIXED = −m_gg(q+g)²                (endpoint = the fixed-source–like attractor),
A_MIXED = 0  ⇔  λ* = (q − g)/(2q)             (the algebraic interior zero).
```

"Spans" means over the FULL combined parameter-and-magnitude domain (all `q,g>0`, `0≤λ≤1`), NOT for every
fixed `q/g`. The three admissibility regimes (`q,g>0`, `0≤λ≤1`) are:

| Regime | `λ*=(q−g)/2q` | Behavior on `0≤λ≤1` |
|---|---|---|
| `q > g` | `∈ (0, ½)` admissible | positive (`λ<λ*`), null (`λ=λ*`), negative (`λ>λ*`) — all three occur |
| `q = g` | `= 0` | zero ONLY at `λ=0`; negative in the interior |
| `q < g` | `< 0` (inadmissible) | endpoints both negative + negative slope ⇒ MIXED stays negative on `[0,1]` |

The full-domain span is witnessed by the cross-regime samples `(q=3,g=1)`, `(q=1,g=1)`, `(q=1,g=2)` at `λ=0`
giving signs `(+, 0, −)`. Claiming every fixed `q/g` spans all signs (rather than the full domain) fires
`PASS_MIXED_FULL_DOMAIN_SPANS`; a positive algebraic root for `q<g` fires `PASS_MIXED_REGIME_Q_LT_G`.

## 3. Q-AMEND consistency → `internal_inconsistency = none` (`PASS_AMEND_REPLACE`, `PASS_AMEND_ADD`, `PASS_ZERO_LEDGER_13`, `PASS_S_HOLD_SCOPE`, `PASS_INTERNAL_INCONSISTENCY_NONE`)

Charge must coexist with the gravity substrate without contradiction. The two admissible amendments to the G0
ledger are:

- **REPLACE:** the core holder RE-TYPES the existing `h`-source BC (`source_row =
  core_holder_retypes_existing_h_source_BC`, **new_rows = 0**).
- **ADD:** the existing `h`-source BC is unchanged and ONE core `h`-holding row is added
  (`new_rows = (core_embedding_h_holding_row:R_w_odd)`, **new_rows = 1**).

Both keep `S_hold` scoped to `r_B−½` only (extending `S_hold` to `h` fires `PASS_S_HOLD_SCOPE`) and both
preserve the **13 `[POSTULATE]` zero rows** (the neutral zero-ledger of bulk / dynamic-drain / mixing /
Berry / cubic-quartic / independent-primitive / phase-drain / wall-collar / return-kernel / drain-source /
geon-derivative / viscosity / mixture terms; dropping one fires `PASS_ZERO_LEDGER_13`). The
`amendment_sectors` map returns the EMPTY tuple over {replace_ledger, add_ledger, S_hold, zero_ledger} ⇒

```text
internal_inconsistency = none        (computed, owned + verified at 032).
```

Injecting a Q-AMEND inconsistency sector fires `PASS_INTERNAL_INCONSISTENCY_NONE`; the same empty-sector fact
is the `internal=False` input to the §4 landing (an inconsistency would route to `NO_GO(sector)` first).

## 4. The typed Q-BC classifier and the four controls (`PASS_BC_ACTUAL_CLASSIFIER`, `PASS_BC_*_CONTROL`, `PASS_ADMISSIBLE_CLASSES`)

The committed model's ACTUAL evidence — admissible variation present, but with the missing parent-throat /
boundary closure (`NONZERO_HA_REQUIRES_CORE_HOLDER`, consumed from 031) — classifies to

```text
Q_BC(actual) = UNDETERMINED_ANALYTICALLY(missing parent-throat/boundary closure).
```

Removing the missing-closure premise fires `PASS_BC_ACTUAL_CLASSIFIER`. The classifier is proved able-to-fail
by the four control evidences that DO determine a class:

| Control | Evidence | Class |
|---|---|---|
| free mouth | admissible variation, no imposed data, no signed stationary / barrier | `FIXED_SOURCE` |
| imposed value | essential value, no admissible variation | `DIRICHLET_VALUE` |
| imposed conormal | fixed conormal + admissible variation | `FIXED_MONOPOLE` |
| imposed relation | mixed relation + admissible variation | `MIXED` |

The admissible class set is exactly `{DIRICHLET_VALUE, FIXED_MONOPOLE, FIXED_SOURCE, MIXED}` (4; dropping MIXED
fires `PASS_ADMISSIBLE_CLASSES`).

## 5. Q-COMBINE: REPLACE/ADD totals, variant, no-double-count (`PASS_REPLACE_TOTALS`, `PASS_ADD_TOTALS`, `PASS_VARIANT_UNRESOLVED`, `PASS_NO_DOUBLE_COUNT`, `PASS_OUTCOME_NOT_INVARIANT`)

The two amendment realizations give DIFFERENT per-class totals:

```text
REPLACE:  DIRICHLET = m_gg φ²/S_gg²,  MONOPOLE = m_gg q²,     FIXED_SOURCE = −m_gg j²,       MIXED = m_gg(1−2λ)q²
ADD:      DIRICHLET = A_V,            MONOPOLE = A_M,          FIXED_SOURCE = A_J,            MIXED = A_MIXED
```

(REPLACE re-types the existing BC with no added `g`-holder ⇒ no `g` contribution; ADD dresses the source with
the core `h`-holding row ⇒ the `g`-dressed `A_X`.) The class×variant outcome table is

| class | replace | add |
|---|---|---|
| DIRICHLET_VALUE | POSITIVE_R2 | POSITIVE_R2 |
| FIXED_MONOPOLE | POSITIVE_R2 | outcome_not_invariant |
| FIXED_SOURCE | NEGATIVE_R2 | NEGATIVE_R2 |
| MIXED | outcome_not_invariant | outcome_not_invariant |

The REPLACE-vs-ADD realization is **`variant_unresolved`** (selecting `replace` fires
`PASS_VARIANT_UNRESOLVED`). The held-`h` monopole work is not double-counted:
`m_gg(q+g)² − 2g·m_gg(q+g) = A_M` (adding a spurious `m_gg g²` fires `PASS_NO_DOUBLE_COUNT`). Because the
class×variant outcomes are not all equal, the flattened outcome is `outcome_not_invariant` (forcing all
outcomes positive fires `PASS_OUTCOME_NOT_INVARIANT`).

## 6. Q-MAG and the `Q_E`/magnitude re-home (`PASS_MAGNITUDE_FREE_FACTOR`, `PASS_QMAG_*_HOOK`, `PASS_MAGNITUDE_COBLOCKER`)

The magnitude carries an unbounded free factor `{c_a, c_ξ}` (`magnitude_free_factor`; hiding `c_ξ` fires
`PASS_MAGNITUDE_FREE_FACTOR`). The four Q-MAG hooks are all negative / undetermined:

```text
density        = NO(no_local_prediction)
radial_monopole = UNDETERMINED(core source/conormal not fixed)
modulus        = NO(continuous core modulus)
close_range    = UNDETERMINED(out_of_scope R comparable to r_e)
```

so `Q_E`/charge-magnitude is **re-homed from pathA_38 into the puncture-deflection build** as the deferred
co-blocker **`R1_REQUIRED(magnitude)`** (replacing it with `R1_REQUIRED(bc_selection)` fires
`PASS_MAGNITUDE_COBLOCKER`). Both the sign (`bc_selection`) and the magnitude discharge only under the shared
sim-deferred throat solve.

## 7. The sealed 23040-cell §4 landing ladder → `R1_REQUIRED(bc_selection)`

The landing is a first-matching-predicate ladder over an EXHAUSTIVE typed grid (`adjudicate`), cross-checked
against an independent `declarative_oracle` on every cell.

**Exhaustive totality (`PASS_TRUTH_TABLE_TOTAL`, `PASS_TRUTH_TABLE_COUNTS`, `PASS_TRUTH_TABLE_DIGEST`,
`PASS_VERDICT_TOTALITY`).** The product of the typed domains (5 BC classes × 6 replace-outcomes ×
6 add-outcomes × 4 variants × 2 magnitudes × 2 tiers × 2³ booleans) is **23040 cells**, every one classified
(no `unclassified`), with `adjudicate == declarative_oracle` on all cells. The grouped counts are

```text
NO_GO(sector)                    11520
R1_REQUIRED(variant_selection)    3840
R1_REQUIRED(sign_and_magnitude)   2856
R1_REQUIRED(bc_selection)         1152
R1_REQUIRED(mixed_bc_parameters)  1152
MECHANISM_FALSIFIED(wrong_range)  1008
MECHANISM_FALSIFIED(wrong_sign)    504
R1_REQUIRED(subleading)            504
SIGN_EARNED                        252
R1_REQUIRED(magnitude)             252     (Σ = 23040)
```

The serialized truth table hashes to the **committed digest**

```text
7627417ace0f99a17187a90efa2523ca9b68df8b64f7960d38be2dc6f553ac49
```

asserted against the committed literal (not a same-path copy); incrementing the total, the `bc_selection`
count, or corrupting the first serialized row each fires its own tooth.

**Precedence + production landing (`PASS_VERDICT_PRECEDENCE`, `PASS_PRODUCTION_LANDING`).** The class gap
fires BEFORE the variant / MIXED-parameter / tier / magnitude gaps. The production tuple (from §§1–6:
`qbc = UNDETERMINED_ANALYTICALLY`, `replace_outcome = add_outcome = outcome_not_invariant`,
`variant = variant_unresolved`, `magnitude = magnitude_free_factor`, `tier = tier_A_conditional`,
`internal = False`, `all_classes_agree = False`, `mixed_range_invariant = False`,
`overall_outcome = outcome_not_invariant`) lands at the FIRST matching predicate

```text
qbc == UNDETERMINED_ANALYTICALLY  ∧  not all_classes_agree   ⇒   R1_REQUIRED(bc_selection).
```

A `SIGN_EARNED` invariance witness (an `UNDETERMINED_ANALYTICALLY` cell with all classes agreeing POSITIVE_R2,
`magnitude_no_free_factor`) confirms the ladder can land elsewhere — the production first match is genuinely
`R1_REQUIRED(bc_selection)`, not a hardcoded token.

**Landing ownership (`PASS_LANDING_OWNERSHIP`, `PASS_TYPED_NEUTRAL_FACTS`, `PASS_TARGET_BLINDNESS`).** The
landing is RECOMPUTED from the typed neutral upstream facts (Q-BC status, the class×variant outcomes, variant,
magnitude, tier, inconsistency sectors) and compared against a FRESH adjudication; an INJECTED landing token
(or any `R1_REQUIRED`/`SIGN_EARNED` string smuggled into the upstream facts) is REJECTED. A named alternative
(`qbc → FIXED_SOURCE`) re-lands at `R1_REQUIRED(sign_and_magnitude)` (`≠` the production token), proving
Q-BC/Q-COMBINE own the tuple rather than the adjudicator being hardcoded. The upstream vocabulary is
target-blind: it contains `outcome_not_invariant` but NO `R1_REQUIRED`/`SIGN_EARNED` token, and no coefficient
mentions a target sign (injecting `SIGN_EARNED` into an upstream token fires `PASS_TARGET_BLINDNESS`).

## 8. Units firewall and range controls (`PASS_UNITS_A`, `PASS_UNITS_U`, `PASS_UNITS_F`, `PASS_RANGE_*`, `PASS_TYPED_NEUTRAL_FACTS`)

The ensemble/landing dimensions are exponent-triple homogeneous, units-restored `[L,T,M]`, able-to-fail:

```text
[A] = E L = M L³ T⁻²   ([φ]=[λ]=1, [q]=[j]=[g]=E),
[U] = E,
[F] = E L⁻¹ = M L² T⁻²·L⁻¹ = M L T⁻².
```

(`[A]=E`, `[U]=EL`, `[F]=E` each fire their own tooth.) Range controls certify the invariance classifier is
able-to-fail: a sign flip `(−1,0,1)` and a zero-touch-without-flip `(1,0,1)` both classify
`outcome_not_invariant`, while a derived subdominant constant `(3/4, 1/2)` classifies `CONSTANT_OUTCOME`
(mis-labelling any of the three fires its tooth).

## Scope, deferrals, and first-class departures (recorded honestly)

1. **The SIGN is `R1_REQUIRED(bc_selection)` — neither EARNED nor CALIBRATED.** The committed bare model
   (G0 + U2-unresolved `𝔅` + `S_hold` scoped to `r_B−½`) does not determine the mouth boundary operator, so it
   does not select among {V,M,J,MIXED} and does not pick the sign. Only the four `A_X` formulas *given a class*
   and `g=g_χh` are DECIDED.
2. **The strict signs are CONDITIONAL.** `A_V>0` (repel), `A_J<0` (attract), and the strict MIXED regimes hold
   only on an admissible nondegenerate response (`m_gg>0` at the `z_g≠0` witness, `S_gg>0`, nonzero
   magnitudes); weak `≥0`/`≤0` in general; at `z_g=0` all `A_X` vanish.
3. **The magnitude is a SEPARATE R1** (`R1_REQUIRED(magnitude)`; `c_a`,`c_ξ` unbounded at tier-A; far-field
   truncation `O(a/R)`), co-blocking with the sign. The **variant / MIXED-parameter / REPLACE-vs-ADD**
   selections are further SEPARATE downstream R1 debts — not part of the sign-class R1 and not resolved by
   picking a class.
4. **The resolver is SIM-DEFERRED.** The nonlinear parent-throat boundary functional / barrier / map `s↦h_A`
   that would SELECT the BC class (and fix the magnitude) is the shared throat solve — a sibling of gravity's
   R10/R30/R33. 032 records it as the R1 resolver and does NOT attempt it. A guaranteed-nonzero far-field
   amplitude needs the same core holder (`NONZERO_HA_REQUIRES_CORE_HOLDER`, R61, consumed from 031).
5. **`internal_inconsistency = none`** is a POSITIVE structural fact (charge coexists with gravity), computed
   from the Q-AMEND sectors — not an absence of checking.

## Source-to-stage predicate manifest (all 36 atomic teeth + Q-block / App-A–E extras = 44 entries)

Completeness certificate: **no silently-dropped source claim**. Every atomic tooth of the source engine
`puncture_deflection_electric_sign_check.py` (the 36-entry `TOOTH_ORDER`) plus the 8 Q-block/App-A–E extras
lands as **PRESERVED** (folded here as-is), **REPLACED_BY_STRONGER** (a stronger reconstruction tooth), or
**SCOPED_OUT** (owned by 031's mechanism / consumed by citation, with a defensible reason). The partition is
disjoint + exhaustive, computed at runtime in BOTH engines:

```text
partition = { PRESERVED: 22, REPLACED_BY_STRONGER: 12, SCOPED_OUT: 10 }   (44 total),
manifest digest = e2cfd11b41cdd3d95111f25215b16e917b1ce0a9623929619338a1a83fa81014.
```

Flipping any one entry's disposition (e.g. `AMEND_REPLACE` PRESERVED→SCOPED_OUT while keeping all 44
identifiers) fires `MANIFEST_MISDISPOSITION`; dropping one entry fires `PASS_SOURCE_PREDICATE_MANIFEST`.

### The 36 atomic source teeth

The Owner column is the EXACT canonical engine token — part of the digest-protected `(id, disposition, owner)`
triple in both engines' `SOURCE_MANIFEST` (verbatim); the Reason column is explanatory only.

| Source tooth | Disposition | Owner (canonical token) | Reason |
|---|---|---|---|
| `FIELD_IDENTITY_UNITS` | SCOPED_OUT | `STAGE031_MECHANISM` | 031 mechanism |
| `FIELD_PARENT_MAP` | SCOPED_OUT | `STAGE031_MECHANISM` | 031 mechanism |
| `FIELD_LIVE_QH` | SCOPED_OUT | `STAGE031_MECHANISM` | 031 mechanism |
| `ACTION_TRANSCRIPTION` | SCOPED_OUT | `STAGE030_031_CITED` | 030/031 cited |
| `AMEND_REPLACE` | PRESERVED | `STAGE032_Q_AMEND` | Q-AMEND (REPLACE, new_rows=0) |
| `AMEND_ADD` | PRESERVED | `STAGE032_Q_AMEND` | Q-AMEND (ADD, one core h-holding row) |
| `SHOLD_SCOPE` | PRESERVED | `STAGE032_Q_AMEND` | Q-AMEND (`S_hold` = `r_B−½` only) |
| `ZERO_LEDGER` | PRESERVED | `STAGE032_Q_AMEND` | Q-AMEND (13 zero rows) |
| `MATRIX_DETERMINANT` | SCOPED_OUT | `STAGE031_RESPONSE` | 031 response (`det m=z_g²/D`) |
| `BC_ACTUAL_GAP` | REPLACED_BY_STRONGER | `STAGE032_TYPED_Q_BC` | typed Q-BC classifier |
| `BC_FREE_CONTROL` | PRESERVED | `STAGE032_TYPED_Q_BC` | typed Q-BC (→ FIXED_SOURCE) |
| `BC_VALUE_CONTROL` | PRESERVED | `STAGE032_TYPED_Q_BC` | typed Q-BC (→ DIRICHLET_VALUE) |
| `BC_MONOPOLE_CONTROL` | PRESERVED | `STAGE032_TYPED_Q_BC` | typed Q-BC (→ FIXED_MONOPOLE) |
| `BC_MIXED_CONTROL` | PRESERVED | `STAGE032_TYPED_Q_BC` | typed Q-BC (→ MIXED) |
| `FORCE_V_FUNCTIONAL` | REPLACED_BY_STRONGER | `STAGE032_FORMULA_SIGN` | `A_V` formula + sign (response-route reshape) |
| `FORCE_M_HWORK` | REPLACED_BY_STRONGER | `STAGE032_FORMULA_SIGN` | `A_M` formula + sign |
| `FORCE_J_FUNCTIONAL` | REPLACED_BY_STRONGER | `STAGE032_FORMULA_SIGN` | `A_J` formula + sign |
| `FORCE_MIXED_FUNCTIONAL` | REPLACED_BY_STRONGER | `STAGE032_FORMULA_SIGN` | `A_MIXED` formula + sign |
| `MIXED_FULL_RANGE` | REPLACED_BY_STRONGER | `STAGE032_THREE_REGIMES` | MIXED three-regime admissibility |
| `FALLOFF` | SCOPED_OUT | `STAGE031_NEUTRAL_SHELL` | 031 neutral shell (`1/R²`) |
| `UNITS_RESTORED` | REPLACED_BY_STRONGER | `STAGE032_UNITS_FIREWALL` | ensemble units firewall `[A]/[U]/[F]` |
| `COMBINE_REPLACE` | PRESERVED | `STAGE032_Q_COMBINE` | Q-COMBINE (REPLACE totals) |
| `COMBINE_ADD` | PRESERVED | `STAGE032_Q_COMBINE` | Q-COMBINE (ADD totals) |
| `NO_DOUBLE_COUNT` | PRESERVED | `STAGE032_Q_COMBINE` | Q-COMBINE (no double-counted held-h work) |
| `RANGE_SIGN_FLIP` | PRESERVED | `STAGE032_RANGE` | range control (sign flip) |
| `RANGE_TOUCH_ZERO` | PRESERVED | `STAGE032_RANGE` | range control (zero-touch) |
| `RANGE_SUBDOMINANT` | PRESERVED | `STAGE032_RANGE` | range control (subdominant constant) |
| `MAG_FREE_FACTOR` | PRESERVED | `STAGE032_Q_MAG` | Q-MAG (`c_a,c_ξ` unbounded) |
| `DENSITY_HOOK` | PRESERVED | `STAGE032_Q_MAG` | Q-MAG (density NO) |
| `MONOPOLE_HOOK` | PRESERVED | `STAGE032_Q_MAG` | Q-MAG (radial monopole UNDETERMINED) |
| `MODULUS_HOOK` | PRESERVED | `STAGE032_Q_MAG` | Q-MAG (modulus NO) |
| `VERDICT_TOTALITY` | REPLACED_BY_STRONGER | `STAGE032_GRID_COUNTS_DIGEST` | 23040-cell grid + counts + digest |
| `VERDICT_PRECEDENCE` | PRESERVED | `STAGE032_LADDER` | class gap fires before variant/tier/mag |
| `LANDING_OWNERSHIP` | REPLACED_BY_STRONGER | `STAGE032_UPSTREAM_RELANDING` | upstream-relanding (injection rejected) |
| `TARGET_BLINDNESS` | PRESERVED | `STAGE032_NEUTRAL_UPSTREAM` | neutral upstream vocabulary |
| `DUAL_ENGINE_TERMS` | REPLACED_BY_STRONGER | `STAGE032_INDEPENDENT_ENGINES` | independent SymPy/Wolfram engines |

### The 8 Q-block / App-A–E extras

| Extra source claim | Disposition | Owner (canonical token) | Reason |
|---|---|---|---|
| `EXTRA_VARIANT_UNRESOLVED` | PRESERVED | `STAGE032_Q_COMBINE` | Q-COMBINE (`variant_unresolved` landing input) |
| `EXTRA_QMAG_CLOSE_RANGE` | PRESERVED | `STAGE032_Q_MAG` | Q-MAG (close-range UNDETERMINED) |
| `EXTRA_SGG_SELF_RESPONSE` | SCOPED_OUT | `STAGE031_CONSUMED` | 031 consumed (`S_gg` definition) |
| `EXTRA_G0_RESPONSE_WITNESS` | SCOPED_OUT | `STAGE031_RESPONSE` | 031 response (star `m*`) |
| `EXTRA_RESPONSE_MATRIX_ESCAPE_FACTORS` | SCOPED_OUT | `STAGE031_RESPONSE` | 031 response (`z_g`,`z_b` opaque) |
| `EXTRA_EXTERIOR_CORE_GAP` | SCOPED_OUT | `STAGE031_CONSUMED_NAMED_FACT` | 031 consumed named fact (`NONZERO_HA_REQUIRES_CORE_HOLDER`) |
| `EXTRA_SECTION4_DIGEST_COUNTS` | REPLACED_BY_STRONGER | `STAGE032_COMMITTED_LITERAL` | committed-literal digest + counts |
| `EXTRA_APP_E_ACCEPTANCE` | REPLACED_BY_STRONGER | `STAGE032_STANDALONE_ABLATIONS` | standalone per-tooth ablations |

## Consumes (by citation, NOT re-derived)

- **From stage 031** (`THROAT_H_SOURCE_1_OVER_R2`): the response matrix `m`; the pair coupling
  `m_gg = B_eff z_g²/D`; `det m = z_g²/D`; the self-response `S_gg`; the neutral shell `A = m_gg·C`
  (`[A]=EL`, `[C]=E²`, `1/R²`, `s₁s₂`); and the named fact `NONZERO_HA_REQUIRES_CORE_HOLDER`. The ownership
  firewall (`PASS_SCOPE_FIREWALL`) rejects re-deriving the response matrix, the `1/R²` structure, or the
  mouth mechanism into 032.
- **From stage 030** (transitively through 031): the kernel determinant `D` and its dimensionless witness
  `D* = 7/4` (`PASS_STAGE030_TRANSITIVE_D_WITNESS`).

`{m_gg, z_g, det m, S_gg}` are CONSUMED (no re-count); `D` is transitive from 030 — none re-derived.

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage032_electric_sign_bc_ensembles_landing_sympy_audit.py` — **SymPy 57 PASS**.
  `mathematica/ledger_stage032_electric_sign_bc_ensembles_landing_mathematica_audit.wl` —
  **Mathematica 57 PASS**. Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]`), no
  argparse harness, no JSON/YAML payload, zero file-I/O between engines; float-free / machine-real-free
  payload throughout. Verdict token `R1_REQUIRED(bc_selection)`; co-blocker `R1_REQUIRED(magnitude)`;
  resolver `SIM_DEFERRED(parent-throat boundary functional/barrier/map s→h_A)`; `INTERNAL_INCONSISTENCY=none`.
- **The `.wl` is a genuinely independent route.** It uses native Wolfram machinery: `Solve`/`Reduce` for the
  ensemble coefficients, the MIXED zero, and the three-regime case split; an independent native construction
  of the exhaustive 23040-cell landing table + SHA-256 digest; an independent landing-ladder evaluation over
  the typed facts. Both engines compute the truth-table digest, the manifest partition digest, and the
  production landing at RUNTIME and agree.
- **58 mutations total, each `FIRED_AT_OWN_ASSERT`** (env switch `LEDGER_STAGE032_MUTATION`): the 57 teeth
  (a primitive mutation for every folded predicate — the four `A_X` formulas, the weak/strict/degenerate sign
  typing, the MIXED endpoints/interior-zero/three-regime teeth, the Q-AMEND/Q-BC/Q-COMBINE/Q-MAG blocks, the
  range controls, the units firewall, the exhaustive total/counts/digest, the precedence + production landing,
  the upstream-relanding ownership tooth, target-blindness, and the 44-entry manifest) plus
  `MANIFEST_MISDISPOSITION`. The landing tooth mutates an UPSTREAM fact (Q-BC / class-variant outcome) and
  re-lands on a NAMED DIFFERENT §4 landing — never the computed summary (the 030 X≡X lesson).
- **Tri-review CLEAN.** Fidelity + adversarial fresh-agent review both CLEAN; per-tooth ablation confirms
  every tooth (incl. the mundane units/anchor/range checks) fires at its own assert; the landing tooth is
  genuinely upstream-driven; no vacuous / subsumed / X≡X tooth. Directive Codex→Grok→Codex bookend + build +
  fresh-agent re-verify all CLEAN.

## Downstream consumers

- **Stage 033 (IV-4):** `NATIVE_P_NO_EMERGENT_GAUSS` (the native-`P` departure) — independent; NOT this stage.
- **Parameter register:** rows + edges R62–R65 added (the class numerator `C`, the four `A_X` decided-
  conditional coefficients, the R1 ensemble-selection data `{φ,q,j,λ}`, `internal_inconsistency=none`, the
  `Q_E`/magnitude re-home); `{m_gg, z_g, det m, S_gg}` consumed from 031 + `D` from 030 by citation (no
  re-count).
- **Part VII:** the sign + magnitude R1 debts and their shared sim-deferred throat resolver enter the
  codimension / knob audit as reduction debt with a named address, alongside gravity's R10/R30/R33.

## Provenance

- **Physics source:** `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` (§4 landing
  table + the Q-BC / Q-AMEND / Q-COMBINE / Q-MAG blocks + App-A–E) + its
  `puncture_deflection_electric_sign_check.{py,wl}` + `..._independent_recompute.py` (the MIXED three-regime
  admissibility at `:109`).
- **Consumes:** stage 031 `{m, m_gg, det m, S_gg, A=m_gg·C, 1/R², NONZERO_HA_REQUIRES_CORE_HOLDER}` + stage
  030 `{D, D*=7/4}` (cite, no re-count).
- **Governing:** `notes/ledger_v2_blueprint.md` §5/§6; `notes/part4_charge_atomic_split.md` (IV-3 row +
  register preview); `docs/model_map.md` §3.4. Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage032_reshape_directive.md`. ⛔ **Not retained** — it lived in gitignored `_scratch/` and no copy survives; this line records that a directive existed, it is not an auditable citation.
