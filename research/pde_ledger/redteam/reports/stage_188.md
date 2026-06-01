---
unit_id: 188
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage188_branch_observables_completion.md]
  paper_appendix: present
---

# Audit unit 188 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_188.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage188_branch_observables_completion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 107 status row; line 265 part overview; lines 881-916 "Branch-observable coordinates" subsection; line 1466 anchor MTDC-T9.4)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The stage card's `\stagefield{Output}` reads: "Compiles direct branch observables to \((\Theta_1,\Xi_1,\mathcal R_1)\) and records equivalent observable packets." The notes elaborate four deliverables: (1) the exact coefficient identity `A_tr,* = B_* C_tr,*`; (2) an exact invertible first-order compiler `C_obs->quot` from the branch-observable packet `(δln R_tr, δln N_*, δln ε_η)` to the Stage-238 tangent quotient packet (matrix `diag(-1/C_tr,*, 1, 1)`, determinant `-1/C_tr,* ≠ 0`); (3) an exact invertible compiler `C_obs->def` from the same observable packet to the defect packet `(Θ_1, Ξ_1, R_1)` (a lower-triangular matrix with determinant `-ε_η,*/(1-ε_η,*) ≠ 0`), which must factor exactly as `C_obs->def = C_quot->def · C_obs->quot`; and (4) the sharp zero-defect theorem `Θ_1=Ξ_1=R_1=0 ⇔ δln R_tr=δln N_*=δln ε_η=0`, plus the complementary-observable drift `δln(1-ε_η) = R_1+Ξ_1 = -(ε_η,*/(1-ε_η,*)) δln ε_η`. The appendix subsection (lines 881-916) confirms the same `R_tr`, `N_* := T^2 R_tr^{B_*}`, `B_*`, and `C_*` definitions verbatim.

## What the script claims to verify

The SymPy script defines the coherent-branch coefficients `Ctr`, `Cstar=1/Ctr`, `Bstar`, `Astar` directly from the closed forms the paper states, then verifies: (I) the coefficient identity `Astar - Bstar*Ctr == 0`; (II) the observable->quotient compiler matrix and its independently computed inverse (`C_quot->obs * C_obs->quot - I == 0` both orders); (III) the quotient->defect compiler reproduces the paper's componentwise `Θ/Ξ/R` identities; (IV) the central factorization claim `C_obs->def = C_quot->def · C_obs->quot` matches the independently typed expected direct compiler, with `det = ε_η/(ε_η-1)`; (V) the independently inverted `C_def->obs` matches the paper's expected inverse and round-trips the observable packet; (VI) the complementary-observable drift identity; and (VII) a "zero-set equivalence" block. This matches the four paper deliverables one-for-one.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `A_tr,* = B_* C_tr,*` | line 60 `Astar - Bstar*Ctr == 0` | match |
| (2) `C_obs->quot` invertible, det `-1/C_tr` | lines 71-85 matrix + det print + inverse-identity checks | match |
| (3a) `C_obs->def` componentwise `Θ,Ξ,R` | lines 136-141 (+ III lines 107-112) | match |
| (3b) factorization `C_obs->def = C_quot->def · C_obs->quot` | line 127 `C_obs_to_def - C_obs_to_def_expected == 0` | match |
| (3c) `det C_obs->def = -ε/(1-ε) ≠ 0`; inverse `C_def->obs` | lines 133, 156-161 | match |
| (4a) zero-defect iff theorem | invertibility (II/IV/V) genuinely establishes it; section VII presents it vacuously | partial |
| (4b) `δln(1-ε_η) = R_1+Ξ_1` | line 167 `(Rcal+Xi) - dln_E == 0` | match |

Dominant pattern: aligned. The single `partial` is a presentation defect in Section VII, not a paper↔script disagreement (see F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(Cstar*Ctr - 1) == 0` | none (sanity) | no — tautological (Cstar := 1/Ctr) |
| A2 | sympy | 60 | `simplify(Astar - Bstar*Ctr) == 0` | deliverable 1 | yes |
| A3 | sympy | 84-85 | `C_quot->obs*C_obs->quot - I == 0` (both orders) | deliverable 2 (invertibility) | yes |
| A4 | sympy | 107-112 | `Theta/Xi/R_from_quot - <paper rhs> == 0` | deliverable 3a (via quotient) | yes |
| A5 | sympy | 127 | `C_obs_to_def - C_obs_to_def_expected == 0` | deliverable 3b (factorization) | yes |
| A6 | sympy | 136-141 | `Theta/Xi/R - <paper rhs> == 0` | deliverable 3a | yes |
| A7 | sympy | 156-158 | `C_def->obs - expected == 0`, `inv - I == 0` (both orders) | deliverable 3c (inverse) | yes |
| A8 | sympy | 161 | `C_def->obs*Delta_def - obs == 0` | deliverable 3c (round-trip) | yes |
| A9 | sympy | 167 | `(Rcal+Xi) - dln_E == 0` | deliverable 4b | yes |
| A10 | sympy | 175-176 | `C_def->obs*0 == 0`, `C_obs->quot*0 == 0` | deliverable 4a (claimed) | no — vacuous (M*0=0 always) |

## Findings

### F1 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:169-176`

**What's wrong:**
Section VII is labeled "Zero-set equivalence" and its comment (lines 170-171) claims to "Verify by substituting zero defects and recovering zero observables exactly," i.e. to exercise the paper's central theorem (deliverable 4a):
> `Θ_1=Ξ_1=R_1=0 ⇔ δln R_tr=δln N_*=δln ε_η=0` (notes lines 405-414).

But the two assertions are:
```
zero_obs_from_zero_def  = simplify(C_def_to_obs * Matrix([0,0,0]))
zero_quot_from_zero_obs = simplify(C_obs_to_quot * zero_obs_from_zero_def)
expect_zero("zero observables from zero defect", zero_obs_from_zero_def)
expect_zero("zero quotient packet from zero observables", zero_quot_from_zero_obs)
```
Any matrix multiplied by the zero vector is the zero vector. These checks hold no matter what `C_def_to_obs` or `C_obs_to_quot` contain — they cannot fail and so verify nothing about the *equivalence* (which is the nontrivial direction: that the maps are bijections, so the zero set is *unique*, i.e. nonzero observables cannot map to zero defects). The "⇔" content of deliverable 4a is the invertibility of the two compilers.

**Why this matters:**
The theorem's real content (the iff / shared zero set) is *already* genuinely established elsewhere in the script — by the four invertibility assertions A3 (lines 84-85), A7 (lines 156-158), the printed nonzero determinants `det C_obs->quot = -1/C_tr` (lines 79-80) and `det C_obs->def = ε/(ε-1)` (line 133), and the round-trip A8 (line 161). So the theorem is not under-verified overall. The defect is narrowly that the block presented as the proof of deliverable 4a is vacuous scaffolding; a reader auditing Section VII in isolation would mistake `M*0=0` for evidence of the equivalence. Left alone, this is a silent-pass trap: the section title overstates what its assertions test.

**Required change:**
Replace the vacuous trivial-direction substitution in lines 169-176 with an assertion that actually exercises injectivity / shared zero set. Concretely, verify that a *generic* (symbolic, nonzero) observable packet maps to defects whose simultaneous vanishing forces the observable packet to vanish — equivalently, assert the determinants are nonzero and that the only solution of `C_obs->def · obs == 0` is `obs == 0`. A minimal, non-vacuous replacement:
```python
subbanner("VII. Zero-set equivalence (shared zero set via invertibility)")
# Equivalence is the nontrivial content: both compilers are bijections,
# so Delta_def == 0 forces Delta_obs == 0 (and conversely), and likewise
# for the quotient packet. Exercise it on the GENERIC packet, not on 0.
expect_zero("nonzero det forces unique zero: C_obs->def^{-1} (C_obs->def obs) - obs",
            sp.simplify(C_def_to_obs * (C_obs_to_def * obs)) - obs)
expect_zero("nonzero det forces unique zero: C_obs->quot^{-1}(C_obs->quot obs) - obs",
            sp.simplify(C_quot_to_obs * (C_obs_to_quot * obs)) - obs)
# Determinants are nonzero on the physical domain (printed in II and IV):
#   det C_obs->quot = -1/C_tr,  det C_obs->def = -eps/(1-eps),  0<eps<1.
det_obs_to_def = sp.simplify(C_obs_to_def.det())
det_obs_to_quot = sp.simplify(C_obs_to_quot.det())
expect_zero("1/det(C_obs->def) is finite & nonzero (reciprocal well-defined)",
            sp.simplify(det_obs_to_def * (1/det_obs_to_def) - 1))
expect_zero("1/det(C_obs->quot) is finite & nonzero (reciprocal well-defined)",
            sp.simplify(det_obs_to_quot * (1/det_obs_to_quot) - 1))
```
(The two round-trip-on-generic-`obs` checks are the load-bearing additions; they are non-vacuous because they would fail if either compiler were singular. The `det*(1/det)-1` lines document nonzero-determinant explicitly. Codex may keep the existing `M*0` prints as informational output but they must not be the only checks in this section.)

**Verification:**
After the fix, Section VII contains at least two assertions whose left operand depends on the generic symbols `dln_Rtr, dln_Nstar, dln_epseta` (not the zero vector), and the script still exits 0. The verifier confirms the new `expect_zero` lines reference `obs` (or `C_obs_to_def * obs`) rather than `Matrix([0,0,0])`.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:46,59`

**What's wrong:**
Line 46 defines `Cstar = sp.simplify(1 / Ctr)`. Line 59 then asserts `expect_zero("C_* C_tr,* - 1", Cstar * Ctr - 1)`. Since `Cstar` is literally `1/Ctr`, the expression `Cstar*Ctr - 1 = (1/Ctr)*Ctr - 1` is identically 0 by construction; the assertion cannot fail regardless of the physics. It does not independently exercise the paper's `C_*` definition (appendix eq:app-part05-Tstar-def gives `C_* := (1+χ₀,*)(1+δU,*)(1+χ₀,*+δU,*)/(χ₀,* δU,*)`), which is just the reciprocal of `C_tr,*` by inspection.

**Why this matters:**
Low impact: the genuinely load-bearing identity of Section I is the adjacent `Astar - Bstar*Ctr == 0` (A2), which is fully substantive. But the tautological line is presented in the same banner as a verified "coefficient identity" and pads the PASS transcript with a check that proves nothing. It is the textbook `tautological_check` shape (`x = expr; assert x == expr`).

**Required change:**
Make the `C_*` check independent: define `C_*` from its own closed form (the appendix definition) and assert it equals `1/C_tr,*`, rather than defining `Cstar := 1/Ctr` and checking the reciprocal of itself. Replace line 46:
```python
Cstar = sp.simplify(1 / Ctr)
```
with an independent definition matching appendix eq:app-part05-Tstar-def, and rewrite the line-59 check to compare the two independent expressions:
```python
Cstar = sp.simplify(
    (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)
    / (chi0s * deltaUs)
)
...
expect_zero("C_* - 1/C_tr,*", Cstar - 1 / Ctr)
```
This turns a tautology into a real cross-check (the closed-form `C_*` versus the reciprocal of the independently-built `C_tr,*`); it would fail if either closed form were mistyped. All downstream uses of `Cstar` (the `diag(-Cstar,1,1)` compiler in Section II) are unaffected because the value is unchanged.

**Verification:**
After the fix, line 46 defines `Cstar` from the four positive symbols (not as `1/Ctr`), the Section I assertion reads `Cstar - 1/Ctr` (or `Cstar*Ctr - 1` with `Cstar` now independent), and the script exits 0 with `C_* - 1/C_tr,* = 0` printed. The `det(C_obs->quot)` print (line 80) is unchanged (`-1/C_tr`).

## Independent-derivation check (Mathematica)

No `.wl` script exists for this unit, so `mathematica_transliteration` does not apply.

## Engine cross-check

Only one engine present; `engines_agree: n/a`.

## On the missing Mathematica engine (line-114 judgment)

I did **not** raise a `missing_mathematica` finding. The stage's entire content is a set of exact rational-function matrix identities over four positive symbols (`χ₀,*, δU,*, ε_η,*` and the derived `dln_*` packet): a coefficient identity, two compiler matrices verified by *independent composition and comparison to independently typed target matrices*, two matrix inverses computed independently by SymPy, the printed nonzero determinants, and a complementary-drift identity. These are decidable, closed-form symbolic equalities that SymPy settles completely and unambiguously (every residual prints exactly `0` / the zero matrix). There is no transcendental, branch-cut, integral, or numerically-delicate step where a second CAS would add confidence, and no claimed result that SymPy fails to genuinely verify. This matches established pipeline precedent for SymPy-only non-status-only stages (e.g., 121/122/123). Per the prompt's line-114 guidance, single-engine verification is acceptable here; the card already states "Mathematica audit: none yet."

## Verdict justification

The script is well-aligned with the paper: every one of the four notes-level deliverables maps to a substantive, non-tautological assertion (A2, A3-A8, A9). The central factorization claim (deliverable 3b, the paper's §5 hinge) is genuinely tested at line 127 by composing `C_quot->def · C_obs->quot` and comparing to an independently typed target matrix — it would fail on any sign or factor error, and indeed the output shows the simplified equality. Invertibility (the real content of the zero-defect theorem) is genuinely established by four `*inv - I` checks plus printed nonzero determinants. I attacked the constants (Ctr, Bstar, Astar all match the paper/appendix closed forms exactly), the `positive=True` assumptions (justified: χ₀,*, δU,*, ε_η,* are physical positive quantities, and `0<ε_η,*<1` is the paper's stated domain for the `-ε/(1-ε)` entries), and the symbol wiring (the `obs` packet symbols are genuinely carried through every compiler). Two genuine defects remain, both low severity and both presentation/scaffolding rather than wrong physics: Section VII's "zero-set equivalence" is verified vacuously (`M*0=0`) when the theorem's real content is invertibility (F1), and the `C_*·C_tr,*−1` check is tautological because `C_* := 1/C_tr` (F2). Neither undermines the paper's claims (which are independently verified elsewhere in the same script), so the verdict is `findings` (not clean), `stop_cold: null`. No downstream propagation: both fixes change only how two already-true facts are *checked*, not any derived value (`Cstar`'s numeric value is unchanged).

## Self-test notes

I checked the self-test traps: (1) Variable independence — there are no `sp.diff` calls in this script, so the differentiate-w.r.t.-unwired-symbol trap does not arise; the proposed F1 fix uses `C_def_to_obs * (C_obs_to_def * obs)` on the generic `obs` packet, which depends on all three `dln_*` symbols, so the round-trip is non-vacuous. (2) No unbounded integrals, so parity is N/A. (3) Trivial-case pre-check — I confirmed the F1 replacement is non-vacuous (round-trip on generic `obs` reduces to `obs - obs = 0` only because the matrices are genuine inverses; a singular matrix would not round-trip) and that F2's independent `C_*` closed form equals `1/C_tr,*` by inspection (numerator/denominator swap). (4) Path: no missing-script finding, so no `.py`/`.wl` directory call needed. (5) Paper round-trip — both fixes leave `Cstar`'s value and all compiler matrices unchanged, introducing no new paper_misalignment.
