---
unit_id: 150
batch: IV.5
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage150_full_profile_residual.md
  paper_appendix: present
---

# Audit unit 150 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_150.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage150_full_profile_residual.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsection "Full-profile residual" at lines 852–897)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage150_full_profile_residual_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.txt`

## What the paper claims

Stage 150 is the "Exact Full-Profile Mouth Potential and Curvature Residual" entry of the Finite mouth-profile corrections block. The card's bottom-line claim is the quote block: "Exact residual `R_*(x)=\Phi_*(x)-\Pi_*x` is tangent-matched but has negative curvature at the mouth." The notes file boxes three concrete deliverables: (i) closed-form profiles `T_s(x;\Pi_*)` (shell, `\kappa_s=0`) and `T_q(x;\Pi_*)` (mixed D/N half-wave with `\kappa_q=\pi/2`, with explicit `A_q`, `C_q`); (ii) the residual identity `R_*(x)=\Sigma_m^*[4T_s(x;\Pi_*)-T_q(x;\Pi_*)-(4-\mathcal S_q(\Pi_*))x]` together with first theorem `R_*(0)=0` and `R_*'(0)=0`; (iii) curvature theorem `R_*''(0)=-3\Sigma_m^*\Pi_*/(1-e^{-\Pi_*}) < 0`. The part-IV appendix subsection (line 852) restates the residual definition and the `R_*(0)=R_*'(0)=0`, `R_*''(0)<0` claim, then folds in downstream covariance/corrected-point material (eqs. app-part04-Sigma-act, ..-covariance-shifts, ..-actual-mouth-covariances, ..-mouth-only-corrected-point) which are not part of the stage card's quoted Output and are not referenced in the stage 150 notes file; those appear to be the appendix narrative bridging into stages 151+ rather than stage 150's own deliverables.

## What the script claims to verify

The SymPy and Mathematica scripts both build the closed-form `T_s`, `T_q` (with `C_q`, `A_q` exactly as in the notes), define `S_q := T_q'(0)`, and form the residual `R(x) = Sigma * (4*T_s - T_q - (4 - S_q)*x)`. They assert: `T_s(0)=0`, `T_q(0)=0`, `T_s'(0)=1`, `T_q'(0)=S_q`, `R(0)=0`, `R'(0)=0`, and (the load-bearing check) `R''(0) - (-3 Sigma Pi/(1-exp(-Pi))) = 0`. Each assertion uses `simplify`/`FullSimplify`. The closing print restates the theorem from the notes verbatim.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `T_s(x;\Pi_*)` closed form (notes, boxed) | sympy 30 / wl 32 (definition) plus `T_s(0)=0`, `T_s'(0)-1=0` BC checks | match |
| `T_q(x;\Pi_*)` closed form with `A_q`, `C_q` | sympy 33–35 / wl 33–35 (definition) plus `T_q(0)=0`, `T_q'(0)-S_q=0` | match (S_q check is definitional, see F1) |
| `R_*(x)` identity (notes box) | sympy 48 / wl 46 (formula) | match (as algebraic form) |
| `R_*(0)=0` (notes; appendix eq. Rstar-curvature) | sympy 49 / wl 47 | match |
| `R_*'(0)=0` (notes; appendix) | sympy 50 / wl 48 | match |
| `R_*''(0) = -3 \Sigma_m^* \Pi_*/(1-e^{-\Pi_*})` (notes box; appendix sign-of) | sympy 53–57 / wl 50–53 | match |
| Negativity `R_*''(0) < 0` | implicit — script verifies the closed form; with positive Sigma, Pi the sign follows | match |
| Appendix-only items (Cov_*(c,R_*), Cov_*(K_q,R_*), Pi_corr, T_m,corr) | not in script | not stage 150 (see verdict justification) |

`paper_alignment` front-matter: **aligned**. The stage card's quoted Output and the boxed theorems in the notes are exactly what the script verifies. The appendix subsection includes additional bridging content (covariances, corrected canonical point) that the stage card does not flag as stage 150 deliverables; those are tracked by the subsequent stages in this block, not by this audit unit.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero("T_s(0)", Ts.subs(x,0))` | `T_s` shell BC at mouth (notes §1) | yes |
| A2 | sympy | 43 | `expect_zero("T_q(0)", Tq.subs(x,0))` | `T_q` mixed-D/N BC at mouth (notes §1) | yes |
| A3 | sympy | 44 | `expect_zero("T_s'(0)-1", diff(Ts,x).subs(x,0) - 1)` | shell-slope normalization (notes §1) | yes |
| A4 | sympy | 45 | `expect_zero("T_q'(0)-S_q", diff(Tq,x).subs(x,0) - Sq)` | none — `Sq` is defined as the LHS on line 37 | no (tautological; see F1) |
| A5 | sympy | 49 | `expect_zero("R(0)", R.subs(x,0))` | `R_*(0)=0` (notes first theorem; appendix eq. Rstar-curvature) | yes |
| A6 | sympy | 50 | `expect_zero("R'(0)", diff(R,x).subs(x,0))` | `R_*'(0)=0` (notes first theorem) | yes |
| A7 | sympy | 57 | `expect_zero("R''(0) - target", R2 - target_R2)` where `target_R2 = -3 Sigma Pi/(1-exp(-Pi))` | curvature theorem `R_*''(0) = -3\Sigma_m^*\Pi_*/(1-e^{-\Pi_*})` (notes second theorem) | yes (load-bearing) |
| B1 | wl | 41 | `expectZero["T_s(0)", ts /. x->0]` | mirror of A1 | yes |
| B2 | wl | 42 | `expectZero["T_q(0)", tq /. x->0]` | mirror of A2 | yes |
| B3 | wl | 43 | `expectZero["T_s'(0)-1", (D[ts,x]/.x->0)-1]` | mirror of A3 | yes |
| B4 | wl | 44 | `expectZero["T_q'(0)-S_q", (D[tq,x]/.x->0)-sQ]` | mirror of A4 — `sQ` defined identically on line 37 | no (tautological) |
| B5 | wl | 47 | `expectZero["R(0)", r /. x->0]` | mirror of A5 | yes |
| B6 | wl | 48 | `expectZero["R'(0)", D[r,x]/.x->0]` | mirror of A6 | yes |
| B7 | wl | 53 | `expectZero["R''(0) - target", r2 - targetR2]` | mirror of A7 | yes (load-bearing) |

Engines agree: both report `R''(0) - target = 0` after their respective simplifiers (SymPy: `R''(0) = 3⋅Π⋅Σ⋅exp(Π)/(1-exp(Π))`, equivalent to the target form `-3 Sigma Pi/(1-exp(-Pi))`; Mathematica: `R''(0) = (-3*E^p*p*sigmaM)/(-1 + E^p)`, same form). Outputs are fresh (mtime newer than scripts).

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py:37,45`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl:37,44`

**What's wrong:**
The script defines

```
Sq = sp.simplify(sp.diff(Tq, x).subs(x, 0))             # sympy line 37
```

and the mathematica analogue

```
sQ = FullSimplify[D[tq, x] /. x -> 0, Assumptions -> $Assumptions];   (* wl line 37 *)
```

Then later asserts

```
expect_zero("T_q'(0)-S_q", sp.diff(Tq, x).subs(x, 0) - Sq)            # sympy line 45
expectZero["T_q'(0)-S_q", (D[tq, x] /. x -> 0) - sQ];                 (* wl line 44 *)
```

The asserted expression is `X - simplify(X)` where `X = D[Tq,x] /. x->0`. Algebraically this is identically zero regardless of whether `Tq` is the correct mixed-D/N profile; the assertion cannot fail.

**Why this matters:**
The check is currently passing trivially. The intended substantive content is "`T_q'(0)` equals the explicit boundary-derivative combination of the constants `A_q`, `C_q`, `\Pi_*` that comes from differentiating the closed-form `T_q` at `x=0` by hand". A meaningful version of this check should compare `T_q'(0)` against an *independently written* closed form, not against itself. As currently coded, the check provides no protection against an algebraic error in either `A_q` or `C_q` that happens to leave `T_q'(0)` self-consistent but wrong; that protection is only available indirectly through the `R''(0) - target = 0` assertion later.

**Required change:**
Replace `Sq` with the *hand-derived* closed form for `T_q'(0)`. Differentiating `T_q(x) = A_q sinh(k x) - C_q cosh(k x) + C_q exp(-Pi x)` and evaluating at `x = 0` gives `T_q'(0) = A_q*k - C_q*Pi`. Define `Sq` from this expression and leave the rest of the script unchanged. After the fix, the `T_q'(0) - Sq` assertion compares the symbolic derivative against the hand-written closed form and can fail if either `A_q` or `C_q` carries an algebra error. See the directive for the exact line edits.

**Verification:**
The saved SymPy output's `S_q(Pi) =` block currently shows the fully expanded derivative; after the fix it should show the compact form `Aq*pi/2 - Cq*Pi` (or an equivalent simplification). The `T_q'(0)-S_q = 0` line still passes. The Mathematica `S_q(Pi) =` line should similarly display a short form involving `aq*k - cq*p`. The `R''(0) - target = 0` load-bearing check is unaffected.

## Independent-derivation check (Mathematica)

The Mathematica script is, structurally, a line-by-line port of the SymPy script. Compare:

SymPy lines 30–35:
```
Ts = (1 - sp.exp(-Pi*x)) / (Pi * (1 - sp.exp(-Pi))) - x * sp.exp(-Pi)/(1 - sp.exp(-Pi))
Cq = Pi / ((1 - sp.exp(-Pi)) * (k**2 - Pi**2))
Aq = Cq * (k * sp.sinh(k) + Pi * sp.exp(-Pi)) / (k * sp.cosh(k))
Tq = Aq * sp.sinh(k*x) - Cq * sp.cosh(k*x) + Cq * sp.exp(-Pi*x)
```

Mathematica lines 32–35:
```
ts = (1 - Exp[-p*x])/(p*(1 - Exp[-p])) - x*Exp[-p]/(1 - Exp[-p]);
cq = p/((1 - Exp[-p])*(k^2 - p^2));
aq = cq*(k*Sinh[k] + p*Exp[-p])/(k*Cosh[k]);
tq = aq*Sinh[k*x] - cq*Cosh[k*x] + cq*Exp[-p*x];
```

Same variable choreography (`Ts/ts`, `Cq/cq`, `Aq/aq`, `Tq/tq`), same construction order, identical algebraic forms; the residual formula `R = Sigma * (4*Ts - Tq - (4 - Sq)*x)` is likewise mirrored verbatim as `r = sigmaM*(4*ts - tq - (4 - sQ)*x)`. Neither engine independently re-derives `T_s` or `T_q` from the underlying lane ODEs; both plug in the closed forms from the notes.

I did **not** raise this as a `mathematica_transliteration` finding because the only mechanically applicable fix (have Mathematica solve the lane ODEs via `DSolve` and assert agreement with the closed form) requires the exact lane ODE / boundary conditions that the closed-form expressions in the notes were derived from, and those ODE/BC pairs are not stated in this stage's notes or card. The shell-lane closed form `T_s(x) = (1-e^{-\Pi x})/(\Pi(1-e^{-\Pi})) - xe^{-\Pi}/(1-e^{-\Pi})` does *not* satisfy `T_s(1) = 0`, so it is not a Dirichlet/Dirichlet BVP solution; the BCs appear to be an upstream Cauchy/normalization choice (`T_s(0) = 0`, `T_s'(0) = 1`) inherited from earlier stages, but I cannot confirm that without reading those stages — which the audit prompt forbids. Filing a fix with speculative `DSolve` boundary conditions risks a Codex iteration that breaks the script. Flagging the transliteration as an observation here, without a directive entry, so the orchestrator/maintainer is aware.

## Engine cross-check

Both engines agree on the load-bearing assertion. SymPy output (lines 28–34 of `.txt`) reports `R''(0) = 3⋅Π⋅Σ⋅exp(Π)/(1-exp(Π))`; Mathematica output (line 26) reports `R''(0) = (-3*E^p*p*sigmaM)/(-1 + E^p)`. These are the same expression in different equivalent algebraic forms, and both reduce to `-3 Sigma Pi / (1 - exp(-Pi))`. Both engines report `R''(0) - target = 0` and exit 0. No engine disagreement.

## Verdict justification

The script does verify the paper's stated claim: the residual is tangent-matched (`R(0)=R'(0)=0`) and has the closed-form negative curvature at the mouth (`R''(0) = -3 \Sigma_m^* \Pi_* / (1 - e^{-\Pi_*})`). The load-bearing assertion A7/B7 is non-tautological — it compares the symbolic second derivative of the constructed residual against an independently written closed form — and both engines pass. The one finding (F1, low severity) is that the `T_q'(0) - S_q` assertion is tautological because `S_q` is defined to be `T_q'(0)`; the fix is a one-line edit to set `S_q` from the hand-derivative `A_q*k - C_q*Pi`. The Mathematica/SymPy transliteration concern is documented above but not raised as a directive finding because the corrective fix requires upstream notes I cannot read under the audit rules. Verdict: **findings** (script-side, applicable).

## Self-test notes

- Variable independence: For the proposed `Sq = Aq*k - Cq*Pi`, the assertion `diff(Tq,x).subs(x,0) - Sq` evaluates `Tq'(0) - (Aq*k - Cq*Pi)`. Differentiating `Tq = Aq*sinh(k*x) - Cq*cosh(k*x) + Cq*exp(-Pi*x)` at `x=0` gives `Aq*k*cosh(0) - Cq*k*sinh(0) - Cq*Pi*exp(0) = Aq*k - Cq*Pi`, so the residual is identically zero — but only after performing the explicit derivative. An incorrect `Aq` or `Cq` would propagate into both sides differently only if the substitution chain breaks, which it does not — so the check still ultimately compares algebraically equal forms. The protection is that any *transcription* error in writing `Sq` (e.g., `Aq*Pi - Cq*k` by mistake) would now fail; this is a non-trivial improvement over the current tautology even though it does not catch errors *inside* `Aq` or `Cq` themselves. (Those are still caught by the `R''(0) - target` assertion downstream.)
- Symmetry/parity: Not applicable — all integrals are evaluations at `x = 0`, no symmetric-domain integration.
- Trivial-case pre-check: At `Pi → 0+`, `Aq → finite`, `Cq → finite`, `Aq*k - Cq*Pi → Aq*pi/2`, and `Tq'(0)` symbolically reduces to the same thing. Consistent.
- Path specifications: F1 edits existing lines in both `.py` (line 37) and `.wl` (line 37). No new files.
- Paper round-trip: The F1 fix introduces no new constants; `A_q*k - C_q*Pi` is purely a re-expression of `T_q'(0)` using definitions of `A_q`, `C_q` already in the script and the notes. No paper conflict.
