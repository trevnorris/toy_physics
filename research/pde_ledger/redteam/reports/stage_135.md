---
unit_id: 135
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage135_outlet_consistent_mouth_closure.md
  paper_appendix: present
---

# Audit unit 135 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_135.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage135_outlet_consistent_mouth_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the include row for stage_135 was relevant; line 1304)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.txt`

## What the paper claims

The stage card (line 16) states verbatim:

> Inheriting the compensated outlet ratio gives \(M_s=4\Sigma_m\), \(M_q=-\Sigma_m\), and \(\Sigma_m^*\approx0.451485\).

The notes elaborate the claim into a set of distinct deliverables: (i) imposing the 4:-1 outlet ratio gives `M_s = 4 Sigma_m`, `M_q = -Sigma_m`; (ii) the Family-1 mouth-layer equation reduces to `Pi = Sigma_m[4 - S_q(Pi)]`; (iii) at the canonical compensation point `Pi_* ≈ 1.50882951349316` the script must confirm `S_q(Pi_*) ≈ 0.6581 < 1`; (iv) `Sigma_m^* ≈ 0.451485277739090`; (v) `M_s^* ≈ 1.80594111095636` and `M_q^* ≈ -0.451485277739090`; (vi) the mixed-lane correction `M_q^* * S_q(Pi_*) ≈ -0.297111597463199`. The card also lists three explicit `\stagefield{Checks}` items: gain pair against outlet consistency; self-matched susceptibility closure; numeric fixed points recorded as numerical.

## What the script claims to verify

The SymPy script (lines 19-49) defines `S_q = S(Pi, pi/2)`, prints the reduced law `Pi_eq = Sigma_m*(4 - S_q)`, hardcodes `Pi_* = 1.50882951349316`, solves `Pi_* = Sigma_var*(4 - S_q(Pi_*))` for `Sigma_var`, prints `Sigma_m^*, M_s^*, M_q^*`, and asserts the closure residual `Pi_* - Sigma_star*(4 - S_star) < 1e-12`. No other assertion is present; `0 < S_q < 1`, the symbolic substitution `M_s + M_q*S_q -> Sigma_m*(4 - S_q)`, the numerical values of `Sigma_m^*, M_s^*, M_q^*`, and the mixed-lane correction are printed but never tested. The Mathematica script (lines 38-79) builds `genericLaw = ms + mq*sQ`, substitutes `{ms -> 4*sigmaM, mq -> -sigmaM}`, derives `reducedLaw`, and `expectZero[reducedLaw - expectedLaw]` checks the substitution; it also `expectApprox`-anchors all numerical deliverables (Sigma_m^*, M_s^*, M_q^*, mixed-lane correction, residual) and `expectTrue`s both `0 < S_q < 1` and `3 Sigma_m^* < Pi_* < 4 Sigma_m^*`.

## Paper ↔ script cross-check

| Paper deliverable | SymPy | Mathematica |
|---|---|---|
| (i) `M_s = 4 Sigma_m`, `M_q = -Sigma_m` (closure choice) | implicit (literal in formula); no test | substitution exercised by `expectZero[reducedLaw - expectedLaw]` (line 56) |
| (ii) Reduced law `Pi = Sigma_m[4 - S_q(Pi)]` | printed; no assertion | `expectZero[reducedLaw - expectedLaw]` (line 56) |
| (iii) `S_q(Pi_*) ≈ 0.6581 < 1` | printed; `bool(...)` evaluated but not asserted (line 41) | `expectTrue["0 < S_q(Pi_star) < 1", ...]` (line 67) |
| (iv) `Sigma_m^* ≈ 0.451485277739090` | computed; no anchored assertion | `expectApprox[..., 0.451485277739089696..., 10^-15]` (line 74) |
| (v) `M_s^*`, `M_q^*` numerical values | computed; no anchored assertion | `expectApprox` (lines 75, 76) |
| (vi) Mixed-lane correction `M_q^* S_q(Pi_*) ≈ -0.297111597463199` | not computed | `expectApprox` (line 77) |
| closure residual identity | `raise AssertionError` if `>1e-12` (line 48) | `expectApprox["closure residual", residual, 0, 10^-14]` (line 78) |

Coverage pattern: Mathematica = match across all six deliverables; SymPy = missing or printed-only for all six. The dominant mismatch is on the SymPy side, hence `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48-49 | `if abs(residual) > 1e-12: raise AssertionError(...)` | none — residual is identically ~0 by construction (see F1) | no |
| A2 | mathematica | 56 | `expectZero["outlet-consistent reduction", reducedLaw - expectedLaw]` | (i)+(ii) substitution → reduced law | partial — algebraically trivial distributive law, but does exercise the `{ms->4σ, mq->-σ}` map |
| A3 | mathematica | 67 | `expectTrue["0 < S_q(Pi_star) < 1", 0 < sStar < 1]` | (iii) | yes |
| A4 | mathematica | 74 | `expectApprox["Sigma_m^*", sigmaStar, 0.451485277739089696..., 10^-15]` | (iv) | yes |
| A5 | mathematica | 75 | `expectApprox["M_s^*", msStar, 1.80594111095635878..., 10^-14]` | (v) | yes |
| A6 | mathematica | 76 | `expectApprox["M_q^*", mqStar, -0.451485277739089696..., 10^-15]` | (v) | yes |
| A7 | mathematica | 77 | `expectApprox["mixed-lane correction", mixedCorrection, -0.297111597463198..., 10^-14]` | (vi) | yes |
| A8 | mathematica | 78 | `expectApprox["closure residual", residual, 0, 10^-14]` | none — identically zero by `Solve` construction (see F1) | no |
| A9 | mathematica | 79 | `expectTrue["3 Sigma_m^* < Pi_* < 4 Sigma_m^*", ...]` | sanity bracket from notes §2 | yes (low value) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:36-49`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl:60-64,78`

**What's wrong:**
In the SymPy script:

```python
Sigma_star = sp.N(sp.solve(sp.Eq(Pi_star, Sigma_var * (4 - S_star)), Sigma_var)[0], 30)
...
residual = sp.N(Pi_star - Sigma_star * (4 - S_star), 30)
...
if abs(residual) > 1e-12:
    raise AssertionError("Outlet-consistent threshold did not close.")
```

`Sigma_star` is defined as the closed-form solution of `Pi_star = Sigma_var*(4 - S_star)`, so the residual `Pi_star - Sigma_star*(4 - S_star)` is identically zero up to floating-point noise no matter what `Pi_star` or `S_star` are. The check cannot fail. The Mathematica script repeats the same construction at lines 60-64 and asserts the residual at line 78 with `expectApprox["closure residual", residual, 0, 10^-14]`.

**Why this matters:**
This is the only assertion in the SymPy script. As a result, the SymPy audit cannot detect any error in `Pi_*`, in the kernel `S(Pi, pi/2)`, in the outlet ratio `4:-1`, or in the carried-forward constant `Sigma_m^* ≈ 0.451485`. The Mathematica residual is similarly empty as a check; the substantive Mathematica content lives in the other `expectApprox`/`expectTrue` lines.

**Required change:**
Replace the tautological residual with an anchored numerical check. Specifically, in the SymPy script, add assertions that anchor `Sigma_star`, `M_s^*`, `M_q^*`, and the mixed-lane correction `M_q^* * S_q(Pi_*)` to the paper-stated values from the stage card (`0.451485`) and the notes file (`0.451485277739090`, `1.80594111095636`, `-0.451485277739090`, `-0.297111597463199`). Also assert `0 < S_q(Pi_*) < 1` (paper claim iii). The Mathematica residual line may be left in place since it is shadowed by the substantive `expectApprox` lines, but should not be the load-bearing check.

**Verification:**
After fix, SymPy `assert` statements for at least four numerical deliverables (Sigma_m^*, M_s^*, M_q^*, mixed-lane correction) and the inequality `0 < S_q(Pi_*) < 1` must appear and the script must exit 0. The verifier should confirm those checks were not present in the prior revision.

### F2 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:19-49` (entire script)

**What's wrong:**
The SymPy script has exactly one `raise AssertionError` (line 49), and that check is tautological per F1. None of the paper's six deliverables (i)-(vi) are asserted on the SymPy side:

- The substitution `Pi = M_s + M_q*S_q(Pi)` with `M_s=4Σ, M_q=-Σ` reducing to `Pi = Σ(4 - S_q)` (deliverables i+ii) is never exercised — the script writes `Pi_eq = Sigma_m*(4 - S_q)` directly (line 29) without constructing the generic two-channel law first.
- The numerical values `Sigma_m^* ≈ 0.451485277739090`, `M_s^* ≈ 1.80594111095636`, `M_q^* ≈ -0.451485277739090` (deliverables iv+v) are printed at lines 42-44 but never compared against the paper-stated values.
- The mixed-lane correction `M_q^* * S_q(Pi_*) ≈ -0.297111597463199` (deliverable vi) is never even computed.
- The inequality `0 < S_q(Pi_*) < 1` (deliverable iii) is evaluated as `bool(0 < S_star < 1)` and printed at line 41, but the boolean result is not asserted.

Compare to the Mathematica script lines 56, 67, 74-77, 79, which cover all six deliverables.

**Why this matters:**
This is a checkpoint-policy violation under the second-engine principle: both engines must independently verify the paper's claims. Currently SymPy is a print-only narrative and contributes nothing to the audit. If the upstream `Pi_*` carry-forward changes, or if `S_q` is miscomputed, the SymPy script will still exit 0.

**Required change:**
Augment the SymPy script with the missing assertions to bring it to parity with Mathematica. Specifically:

1. Construct `generic_law = M_s + M_q * S_q` symbolically with independent symbols `M_s, M_q`, substitute `{M_s: 4*Sigma_m, M_q: -Sigma_m}`, simplify, and assert `simplify(reduced - Sigma_m*(4 - S_q)) == 0`.
2. Assert `bool(0 < S_star < 1)`.
3. Add anchored numerical assertions:
   - `abs(Sigma_star - sp.Float("0.451485277739090", 30)) < 1e-12`
   - `abs(M_s_star - sp.Float("1.80594111095636", 30)) < 1e-11`
   - `abs(M_q_star - sp.Float("-0.451485277739090", 30)) < 1e-12`
4. Compute `mixed_correction = M_q_star * S_star` and assert `abs(mixed_correction - sp.Float("-0.297111597463199", 30)) < 1e-11`.

Tolerances may be loosened if needed; the anchoring digits come from the notes file (lines 86-101, 124-127).

**Verification:**
After the fix, SymPy output must contain at least five new assertion confirmations corresponding to deliverables i+ii, iii, iv, v (two values), and vi. Script must still exit 0. The verifier confirms these checks are present and that the values agree with the Mathematica output to ≥10 significant figures.

## Independent-derivation check (Mathematica)

The Mathematica script is not a line-by-line transliteration of the SymPy script. It introduces an additional construct that the SymPy script lacks: it builds `genericLaw = ms + mq*sQ` (line 51), substitutes `{ms -> 4 sigmaM, mq -> -sigmaM}` (line 52), and then proves `reducedLaw - expectedLaw == 0` via `FullSimplify` (line 56). The SymPy script writes the reduced form directly at line 29 and never constructs the generic two-channel law. Both scripts share the same kernel definition and the same numerical-substitution pattern (`Pi_*` hardcoded → `Sigma_star = Pi_*/(4 - S_q(Pi_*))`), so the two engines are not maximally independent on the numerical side, but the Mathematica script does exercise the substitution claim independently. No `mathematica_transliteration` finding is raised.

## Engine cross-check

Both engines agree on the load-bearing numerical values:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `S_q(Pi_*)` | 0.658075937605428494269581645208 | 0.6580759376054294867405036727 |
| `Sigma_m^*` | 0.451485277739089696513730132210 | 0.4514852777390898089760109878 |
| `M_s^*` | 1.80594111095635878605492052884 | 1.8059411109563592359040439511 |
| `M_q^*` | -0.451485277739089696513730132210 | -0.4514852777390898089760109878 |
| closure residual | -7.22e-17 | 0 (precision-limited) |

Agreement to ≥14 significant figures throughout; differences are at the level of the 14th-15th digit (precision-engineering rounding). No `engine_disagreement` finding.

## Verdict justification

`findings`. The Mathematica script faithfully audits all six paper-side deliverables and is internally consistent. The SymPy script, however, contains a single tautological assertion (F1) and no anchored numerical checks against any of the paper-stated values (F2). The engines agree on the numerics that are computed, and the paper's stated values match what both engines produce, so the unit is in no danger of `paper_misalignment` or `stop_cold`. The script-side fix is straightforward: bring SymPy to parity with Mathematica by adding the missing assertions.

Cosmetic note (not a finding): both scripts print a banner reading `STAGE 118 — OUTLET-CONSISTENT ONE-PARAMETER CLOSURE` (sympy line 27, mathematica line 38), and the notes file's H1 reads `Moving-Throat PDE — Stage 237: Outlet-Consistent Mouth Closure`. The paper card, file paths, and appendix all consistently identify the unit as stage 135. The 118/237 labels are stale prose artifacts and have no effect on the math; they are not flagged as a paper_misalignment because no `\stagefield{Output}` quantity is affected.

## Self-test notes

I checked: (1) the closure residual `Pi_* - Sigma_*(4 - S_*)` — by construction `Sigma_* := Pi_*/(4-S_*)`, so residual is identically 0 (this is F1); (2) symmetry/parity is not relevant since no integrals over unbounded domains appear; (3) trivial-case pre-check on the proposed F2 substitution assertion: `(4σ) + (-σ)*s - σ*(4 - s) = 4σ - σs - 4σ + σs = 0`, so the proposed `assert simplify(generic_law.subs(...) - Sigma_m*(4 - S_q)) == 0` will succeed correctly under any S_q; (4) the proposed numerical-anchor tolerances (1e-11 to 1e-12) are comfortably above the engine precision (~1e-14 to 1e-16 observed); (5) the paper round-trip: the four numerical values requested for SymPy anchoring (Sigma_m^*, M_s^*, M_q^*, mixed correction) all come directly from notes lines 86-101 and 124-127 and match the Mathematica `expectApprox` targets, so no new `paper_misalignment` is introduced by the fix.
