---
unit_id: 169
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage169_no_linear_p2_scalar_slippage.md]
  paper_appendix: present
---

# Audit unit 169 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_169.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage169_no_linear_p2_scalar_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 69, 265, 373-409, 411 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}`: "Pure grouped real \(P_2\) anisotropy cannot linearly source the scalar off-bundle slippages; scalar feed-down starts quadratically." The notes elaborate four supporting deliverables: (1) the exact grouped bilinear invariant `I[x,y] = (1/5) δxᵀ G_grp δy = 4 a_x a_y + (4/5) b_x b_y` with `G_grp = diag(1,2,2)`; (2) the weak-axisymmetric branch signature `b_x = 3 a_x` and `I_axis = (7/10) ε² x⁽¹⁾y⁽¹⁾` (notes §1, eq. boxed); (3) the monopole-selection / harmonic-average fact that `∫_{S²} Y_{2m} = 0`, hence the linear coefficient of any averaged scalar log-observable vanishes (notes §2); and (4) the Stage-168 weighting that transports the quadratic invariants into `ε_⊥` via `Ξ_⊥ = g_* Ξ_T + (g_* + 1/(2√(1+r_*²)))Ξ_v + (2g_* + 3/(4√(1+r_*²)))Ξ_L`, with the numeric Family-1 evaluation `Ξ_⊥ ≈ 0.758035078944663 Ξ_T + 1.00314310113848 Ξ_v + 1.88373219118005 Ξ_L` (notes §5). The card itself is a terse audit card; the notes are authoritative on intent.

## What the script claims to verify

The SymPy docstring enumerates exactly the four deliverables above. The assertions: (1) `I[x,y] − [4 a_x a_y + (4/5) b_x b_y] == 0` by symbolic identity; (2) `b_x − 3 a_x == 0` and `I_axis − (7/10) ε² x1 y1 == 0` under the axisymmetric substitution `(1, 1/2, −1)`; (3) `<Y20> == 0`, `<Y20²> − 1/(4π) == 0`, and the linear-in-`e` coefficient of `<log(X0(1+e·Y20))>` `== 0`; (4) `ε_perp − Ξ_perp·Igrp == 0` (the factorization), plus numeric checks that the three `Ξ_perp` coefficients evaluate to the three paper values within `1e-12`. The Mathematica script asserts the identical set.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| Grouped invariant `I = 4a_xa_y + (4/5)b_xb_y`, `G=diag(1,2,2)` | sympy L54-59 / wl L54-59 | match |
| Weak-axis signature `b=3a`, `I_axis=(7/10)ε²x1y1` | sympy L74-82 / wl L70-78 | match |
| No linear scalar feed-down (`∫Y20=0`, lin coeff of avg log =0) | sympy L96-109 / wl L82-93 | match |
| Stage-168 transport `Ξ_⊥` form + numeric Family-1 coeffs | sympy L121-144 / wl L97-119 | match |

Every stated deliverable has a faithful, non-tautological script-side check. `paper_alignment: aligned`. The one finding (transliteration) is a second-engine-policy defect, not a paper-alignment defect.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(I − [4a a +4/5 b b]) == 0` | deliverable 1 | yes |
| A2 | sympy | 81 | `b_x − 3 a_x == 0` (axis subs) | deliverable 2 | yes |
| A3 | sympy | 82 | `I_axis − (7/10)ε²x1y1 == 0` | deliverable 2 | yes |
| A4 | sympy | 100 | `<Y20> == 0` | deliverable 3 | yes |
| A5 | sympy | 101 | `<Y20²> − 1/(4π) == 0` | deliverable 3 | yes |
| A6 | sympy | 109 | `∂_e <log(X0(1+e Y20))>|_0 == 0` | deliverable 3 | yes |
| A7 | sympy | 125 | `ε_perp − Ξ_perp Igrp == 0` | deliverable 4 (factorization) | partial (algebraic, see note) |
| A8 | sympy | 136-144 | 3× `|coeff − paper| < 1e-12` | deliverable 4 (numeric) | yes |
| A1'-A8' | mathematica | 59,77,78,86,87,93,99,105-118 | same set | same | same as A1-A8 |

Note on A7: `ε_perp` and `Ξ_perp` are both built from the *same* hand-typed weight expressions one line apart (sympy L121-122), so `ε_perp − Ξ_perp·Igrp == 0` is algebraically guaranteed — it verifies that `Igrp` factors out, not that the weights are correct. That is weak in isolation, but the weights are independently pinned by the numeric checks A8 against the three paper values, so deliverable 4 is genuinely exercised overall. Not flagged as a separate `tautological_check` because A8 carries the real anchoring; noted for completeness.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:40-119`
- (mirror of) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:38-144`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The two share identical variable choreography, identical block ordering, identical intermediate definitions, identical assertion targets, and even identical hand-typed numeric literals. Corresponding excerpts:

Invariant block — sympy L43-50:
```
xbar = (xb20 + 2*xb21 + 2*xb22) / 5
ax = (2*xb20 - xb21 - xb22) / 10
bx = (xb21 - xb22) / 2
```
wl L45-50:
```
xBar = (x20 + 2*x21 + 2*x22)/5;
aX = (2*x20 - x21 - x22)/10;
bX = (x21 - x22)/2;
```

Transport block — sympy L121-122:
```
eps_perp = g*XiT*Igrp + (g + sp.Rational(1,2)/s)*Xiv*Igrp + (2*g + sp.Rational(3,4)/s)*XiL*Igrp
Xi_perp  = g*XiT + (g + sp.Rational(1,2)/s)*Xiv + (2*g + sp.Rational(3,4)/s)*XiL
```
wl L97-98:
```
epsPerp = g*xiT*iGrp + (g + 1/(2*s))*xiv*iGrp + (2*g + 3/(4*s))*xiL*iGrp;
xiPerp  = g*xiT + (g + 1/(2*s))*xiv + (2*g + 3/(4*s))*xiL;
```
Same hardcoded `r_num`/`g_num` constants (`1.77799353547498`, `0.758035078944663`) and the same three hand-typed paper-comparison values (sympy L137-139 ↔ wl L110-112). The two engines echo each other's algebra rather than re-deriving the result from the physical premises.

**Why this matters:**
The dual-engine policy requires the second engine to derive each result independently so a transcription error in the first engine cannot be silently mirrored. Here a wrong weight expression (e.g. `3/(4·s)` mistyped) would be copied verbatim into both engines and both would "agree." The most load-bearing quantities — the Stage-168 weights and the three Family-1 coefficients — are precisely the hand-typed literals that are mirrored, so the second engine adds essentially no independent assurance for deliverable 4.

**Required change:**
Make the Mathematica script re-derive the structural pieces rather than transcribe them. Concretely (minimal, non-feature-adding):
1. Derive the invariant from the metric definition only: keep `iXY = (dX . gMat . dY)/5` (already structural) but obtain `aX, bX` etc. and the `4 a a + (4/5) b b` target from `iXY` itself rather than re-typing the closed form, OR derive `<Y20²>=1/(4π)` from the normalization constant rather than asserting the literal `1/(4*Pi)` (already integral-based — acceptable).
2. For the transport coefficients, derive `r_*` (= `Sqrt[4107 - 100 Pi^2]/(10 Pi)`, the canonical Family-1 radius) and `g_*` symbolically from their upstream definitions instead of pasting the decimal literals `1.77799353547498`/`0.758035078944663`, then compute the three coefficients and compare to the paper values. That makes the Mathematica numeric check an independent evaluation rather than a transcription of the SymPy numbers.

If full independence is judged impractical for this audit card (it carries forward Stage-168/Family-1 results), at minimum the `r_*`/`g_*` literals should be replaced by their symbolic upstream forms in the `.wl` so the second engine does not depend on the same hand-typed decimals as the first.

**Verification:**
After the fix, `wl` no longer contains the bare decimals `1.77799353547498` / `0.758035078944663` as the *source* of the coefficients (they may appear only as the paper-side comparison targets), and the three coefficient checks still PASS in the refreshed `mathematica/output/...txt`. The SymPy output is unchanged.

## Independent-derivation check (Mathematica)

Not independent. See F1. The `.wl` mirrors the `.py` block-for-block (banner helpers, `expectZero`, `sphereAvg`, the four numbered blocks, the same substitution rule `subsAxis`, the same hand-typed weights and decimals). Confidence: high that this is a port.

## Engine cross-check

Both engines pass all eight checks and agree. Final symbolic forms differ only in cosmetic grouping (sympy prints the invariant fully expanded, `4*x20*y20/25 - ...`; Mathematica prints the factored form `(-2*(x22*(...)+...))/25` — algebraically identical). Numeric coefficients agree to ~15 digits (sympy diffs `3.97e-15`/`4.46e-15` vs paper; Mathematica prints the 20-digit values, same). `<log...> = -e²/(8π) + log(X0)` in both. No `engine_disagreement`.

## Verdict justification

Every paper deliverable maps to a faithful, non-tautological script check, and both engines pass with agreeing residuals — so `paper_alignment: aligned` and the math holds up. Attacks tried and failed: (a) re-derived `<Y20²> = (5/16π)·(4/5) = 1/(4π)` by hand — matches A5; (b) re-derived the three transport coefficients from `g_*=0.758035…` and `r_*=√(4107−100π²)/(10π)=1.77799…`: `g`, `g+1/(2s)=1.003143`, `2g+3/(4s)=1.883731` — all match the paper values in §5; (c) confirmed `r_num` IS the canonical Family-1 radius, just numeric, and is used only to regenerate the §5 coefficients (so it is internal, not a stray constant); (d) confirmed the axisymmetric subs `(1,1/2,−1)` matches the appendix signature eq. (app-part05-weak-axisymmetric-signature) and `a_x=ε x1/4`, `b_x=3ε x1/4` matches eq. (app-part05-b-equals-3a-general). The single defect is that the second engine is a transliteration of the first, which weakens the dual-engine guarantee on exactly the hand-typed load-bearing constants — hence `verdict: findings` with one `mathematica_transliteration` finding, no stop-cold.

## Self-test notes

Checked: (1) Variable independence — the only derivative is `∂_e` of `log_avg`, which genuinely depends on `e` (the `O(e²)` series has an `e²` term that vanishes under `∂_e|_0`), so A6 is a real zero, not a trivially-zero derivative. (2) Symmetry/parity — `<Y20>=0` is the odd-Legendre/angular-average vanishing, correctly nonzero only at second order; verified `<Y20²>=1/(4π)≠0`. (3) Trivial-case — substituting `x=y=(1,1/2,−1)·ε` reproduces `I_axis=(7/10)ε²` (A3) and the axis `a/b` ratios. (4) Path: F1 targets the existing `.wl` under `mathematica/` (no new-script path needed). (5) Paper round-trip: the prescribed `.wl` change only replaces hand-typed decimals with their symbolic upstream forms and the canonical Family-1 radius `√(4107−100π²)/(10π)`; it introduces no constant absent from the paper/notes.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `I[x,y] = 4 a_x a_y + (4/5) b_x b_y` (with `G_grp=diag(1,2,2)`) | sympy L59 / wl L59; sympy out L10, wl out L10-11 | notes §1 boxed (L65-72), §4 (L209-211); appendix eq. (app-part05-grouped-invariant) L376-383 | MATCH |
| `b_x = 3 a_x` (weak-axis) | sympy L81 / wl L77; out L17/L19 | notes §1 boxed `b_x=3a_x` (L90); appendix eq. (app-part05-b-equals-3a-general) L407 | MATCH |
| `I_axis = (7/10) ε² x⁽¹⁾y⁽¹⁾` | sympy L82 / wl L78; out L18/L20-21 | notes §1 boxed `A_x²=(7/10)ε²(x⁽¹⁾)²` (L92-94), §4 (L249-252) | MATCH |
| `<Y20> = 0` | sympy L100 / wl L86; out L23/L26 | notes §2 `∫ Y_{2m}=0` (L115-118) | MATCH |
| `<Y20²> = 1/(4π)` | sympy L101 / wl L87; out L24/L27 | (normalization scaffolding; not a stated deliverable) | INTERNAL |
| linear coeff of `<log(X0(1+e Y20))>` = 0 | sympy L109 / wl L93; out L28/L33 | notes §2 `δ⁽¹⁾_{P2} S = 0` (L121-125), §3 boxed log-observables (L153-176) | MATCH |
| `Ξ_⊥ = g_* Ξ_T + (g_*+1/(2√(1+r_*²)))Ξ_v + (2g_*+3/(4√(1+r_*²)))Ξ_L` | sympy L122 / wl L98; out L33-34 | notes §5 boxed `Ξ_⊥^{(XY)}` (L279-285) | MATCH |
| coeff on Ξ_T = `0.758035078944663` | sympy L137 / wl L110; out L35/L42 | notes §5 numeric box (L293) | MATCH |
| coeff on Ξ_v = `1.00314310113848` | sympy L138 / wl L111; out L36/L44 | notes §5 numeric box (L294) | MATCH |
| coeff on Ξ_L = `1.88373219118005` | sympy L139 / wl L112; out L37/L46 | notes §5 numeric box (L295) | MATCH |

INTERNAL items (accounted for, no finding): `<Y20²>=1/(4π)` (normalization check); `r_num=1.77799353547498` (the canonical Family-1 radius `√(4107−100π²)/(10π)`, numeric form, used only to regenerate the §5 coefficients — verified by hand to ~6 digits); `g_num=0.758035078944663` (same as the Ξ_T coefficient, doubles as the §5 `g_*`); the symbolic intermediates `a_x,b_x,xbar,δx`; pass/fail flags and the `1e-12` tolerance.

reconciliation: complete; 10 deliverable values checked, 0 misaligned.

No `value_mismatch` or `script_missing_paper_claim` arises from the reconciliation. The only finding remains the `mathematica_transliteration` (F1), which does not route to the user gate.
