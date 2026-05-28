---
unit_id: 169
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
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
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 69, 265, 363-411 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` states verbatim: "Pure grouped real \(P_2\) anisotropy cannot linearly source the scalar off-bundle slippages; scalar feed-down starts quadratically." The notes flesh this out into several concrete deliverables: (1) the exact grouped bilinear invariant \(\mathcal I[x,y]=\frac15\delta x^T G_{\rm grp}\delta y = 4a_xa_y+\frac45 b_xb_y\) with \(G_{\rm grp}=\operatorname{diag}(1,2,2)\) and the anomaly variables \(a_x=(2x_{20}-x_{21}-x_{22})/10\), \(b_x=(x_{21}-x_{22})/2\); (2) the weak axisymmetric signature \((1,1/2,-1)\) giving \(b_x=3a_x\) and \(\mathcal A_x^2=\mathcal I[x,x]=\frac{7}{10}\epsilon^2(x^{(1)})^2\); (3) the monopole-selection theorem that the sphere average of every real \(P_2\) harmonic vanishes, so \(\delta^{(1)}_{P_2}\mathcal S=0\) for every scalar observable and hence \(\varepsilon_L^{(1,P_2)}=\varepsilon_v^{(1,P_2)}=\varepsilon_T^{(1,P_2)}=0\); and (4) the Stage 253 weighting that transports the quadratic invariants into \(\varepsilon_\perp\) with the numeric Family-1 combination \(\Xi_\perp\approx 0.758035078944663\,\Xi_T + 1.00314310113848\,\Xi_v + 1.88373219118005\,\Xi_L\). The bottom-line claim (the Output field) is the linear no-feed-down result; section 5 is a downstream restatement using imported Stage 253/168 weights.

## What the script claims to verify

Both engines (identical structure) verify four blocks. (1) The grouped bilinear identity: build \(\mathcal I[x,y]\) from the matrix definition \(\frac15\delta x^T G\delta y\) and assert it equals \(4a_xa_y+\frac45 b_xb_y\). (2) Substitute the axisymmetric signature and assert \(b_x=3a_x\) and \(\mathcal I_{\rm axis}=\frac{7}{10}\epsilon^2 x_1y_1\). (3) Compute \(\langle Y_{20}\rangle\) and \(\langle Y_{20}^2\rangle\) explicitly over the sphere, assert \(\langle Y_{20}\rangle=0\) and \(\langle Y_{20}^2\rangle=\frac{1}{4\pi}\), then expand \(\langle\log(X_0(1+eY_{20}))\rangle\) to \(O(e^2)\) and assert the linear coefficient is zero. (4) Construct `eps_perp` as a sum of `Xi*Igrp` terms and assert `eps_perp - Xi_perp*Igrp == 0`, then print (not assert) the numeric Family-1 coefficient combination.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Verdict |
|---|---|---|
| (1) Invariant \(\mathcal I=4a_xa_y+\frac45 b_xb_y\), \(G=\mathrm{diag}(1,2,2)\) | block 1, `Ixy - Ixy_expected == 0` | match |
| (2) \(b_x=3a_x\), \(\mathcal I_{\rm axis}=\frac{7}{10}\epsilon^2 x_1y_1\) | block 2, two `expect_zero` | match |
| (3) Linear scalar feed-down vanishes (monopole selection) | block 3, `<Y20>=0`, `<Y20^2>=1/4pi`, lin-coeff `=0` | match (axisymmetric \(Y_{20}\) realization; see note) |
| (4) Stage 253 transport law + numeric combination | block 4, tautological `eps_perp - Xi_perp*Igrp`, numeric only printed | partial |

Block 3 verifies the no-feed-down result on the operative weak-axisymmetric branch (which is built solely on \(Y_{20}\)); the paper's broader "every \(Y_{2m}\) averages to zero" is not exhaustively tested, but for this stage's axisymmetric branch only \(Y_{20}\) is in play, so coverage is appropriate. The dominant pattern is `match`; `paper_alignment: aligned`. Block 4 is the only weak row (covered by F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `expect_zero(Ixy - Ixy_expected)` | claim 1 (invariant) | yes |
| A2 | sympy | 81 | `expect_zero(bx_axis - 3*ax_axis)` | claim 2 (b=3a) | yes |
| A3 | sympy | 82 | `expect_zero(Ixy_axis - 7/10 eps^2 x1 y1)` | claim 2 (7/10 law) | yes |
| A4 | sympy | 100 | `expect_zero(Y20_avg)` | claim 3 (monopole) | yes |
| A5 | sympy | 101 | `expect_zero(Y20_sq_avg - 1/(4pi))` | claim 3 (normalization) | yes |
| A6 | sympy | 109 | `expect_zero(lin_coeff)` | claim 3 (no linear feed-down) | yes |
| A7 | sympy | 125 | `expect_zero(eps_perp - Xi_perp*Igrp)` | claim 4 (transport) | no — tautological |
| A8 | sympy | 130-131 | `print(Xi_perp_num)` (no assert) | claim 4 (numeric) | no — not asserted |
| B1-B7 | mathematica | 59,77,78,86,87,93,99 | `expectZero[...]` | mirror of A1-A7 | mirror |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:121-131`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:96-103`

**What's wrong:**
The block-4 assertion is tautological. The script defines
```
eps_perp = g*XiT*Igrp + (g + 1/(2*s))*Xiv*Igrp + (2*g + 3/(4*s))*XiL*Igrp
Xi_perp  = g*XiT      + (g + 1/(2*s))*Xiv      + (2*g + 3/(4*s))*XiL
```
and then asserts `expect_zero("eps_perp - Xi_perp Igrp", eps_perp - Xi_perp*Igrp)`. By construction `eps_perp` is term-for-term `Xi_perp` with each summand multiplied by the same `Igrp`; factoring `Igrp` is the distributive law and the residual is identically zero regardless of whether the Stage 253 weighting coefficients are correct. The check therefore exercises none of the substantive content of paper section 5 — the load-bearing coefficients are \(\mathfrak g_*\) (on \(\Xi_T\)), \((\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}})\) (on \(\Xi_v\)), and \((2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}})\) (on \(\Xi_L\)). The numeric combination that *would* test those coefficients against the paper's stated values
\(0.758035078944663\,\Xi_T + 1.00314310113848\,\Xi_v + 1.88373219118005\,\Xi_L\)
is only `print`-ed (sympy line 131, mathematica line 103), never asserted.

Note the paper-side numbers are reproduced correctly: with \(g=0.758035078944663\) and \(r=1.77799353547498\), \(s=\sqrt{1+r^2}\approx 2.039917\), the coefficients evaluate to \(g=0.758035\), \(g+\tfrac{1}{2s}=1.003143\), \(2g+\tfrac{3}{4s}=1.883732\), matching the paper to its stated precision. So this is a verification-strength gap, not a value mismatch.

**Why this matters:**
Block 4 is the script's only coverage of paper section 5 (the quadratic transport law / numeric Family-1 combination). As written, that coverage is vacuous: an incorrect weighting coefficient would still pass A7. The numeric forwards \(0.758035\ldots / 1.003143\ldots / 1.883732\ldots\) that downstream Part V citations rely on are displayed but unchecked.

**Required change:**
Keep A7 (harmless), but add a real numeric assertion comparing each computed `Xi_perp` coefficient to the paper's stated Family-1 values. Replace the bare `print` of `Xi_perp_num` with explicit per-coefficient `expect_close`-style checks against `0.758035078944663` (on `XiT`), `1.00314310113848` (on `Xiv`), `1.88373219118005` (on `XiL`), at a tolerance consistent with the paper's quoted digits. Mirror in the `.wl`.

**Verification:**
New assertions appear after sympy line 131 / mathematica line 103; sympy/mathematica outputs show three new PASS lines confirming each numeric coefficient matches the paper value; scripts still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl` (whole file)
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py`

**What's wrong:**
The `.wl` is a line-for-line port of the `.py`, not an independent re-derivation. Corresponding sections:
- Invariant: WL `iXY = FullSimplify[(dX . gMat . dY)/5]; iXYExpected = FullSimplify[4*aX*aY + (4/5)*bX*bY]` (lines 54-55) vs SymPy `Ixy = sp.simplify((dx.T * G * dy)[0] / 5); Ixy_expected = sp.simplify(4*ax*ay + sp.Rational(4,5)*bx*by)` (lines 54-55) — same construction, same intermediate names.
- Axisymmetric subs: WL `subsAxis = {x20 -> x0 + eps*x1, x21 -> x0 + eps*x1/2, x22 -> x0 - eps*x1, ...}` (lines 61-68) vs SymPy `subs_axis = {xb20: x0 + eps*x1, xb21: x0 + eps*x1/2, xb22: x0 - eps*x1, ...}` (lines 65-72) — identical rules.
- Transport block: WL `epsPerp = g*xiT*iGrp + (g + 1/(2*s))*xiv*iGrp + (2*g + 3/(4*s))*xiL*iGrp` and the same hardcoded `rNum = 1.77799353547498`, `gNum = 0.758035078944663` (lines 97, 101-102) vs SymPy lines 121, 127-128 — identical algebra and identical literal constants.

The two-engine policy requires each engine to derive the result independently from the physical premises; here the second engine echoes the first engine's exact variable choreography and even reuses the same hardcoded numeric literals.

**Why this matters:**
A transliterated second engine provides no genuine cross-check — it cannot catch a derivation-level error in the first engine, because it reproduces the same algebra. The "engines agree" status is therefore weaker than it appears.

**Required change:**
This is a known structural pattern across this paper's batches; flag for the same disposition applied to prior `mathematica_transliteration` findings. No mechanical code edit is mandated by this finding alone (the algebra is correct); the disposition (accept-as-mirror vs. re-derive) is a policy call. If the policy is to accept transliteration mirrors, mark resolved; otherwise the `.wl` should re-derive the invariant via an independent route (e.g. direct \(G\)-orthogonal projection / eigen-decomposition rather than re-substituting the same \(a,b\) closed forms).

**Verification:**
If re-derivation is required, the `.wl` invariant block no longer references the same `4*aX*aY + (4/5)*bX*bY` closed form as the construction target; otherwise auditor confirms the policy disposition.

### F3 — stale_output (informational) / banner mislabel

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:31`

**What's wrong:**
The opening banner in both scripts reads `STAGE 152 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE`, and this wrong stage number propagates into both saved transcripts (sympy output line 11, mathematica output line 11). This is a copy-paste artifact (this is stage 169, not 152). Outputs themselves are fresh (sympy output mtime 12:47 > script 11:58; mathematica output 13:20 > script 11:58), so this is not a true `stale_output`; it is a provenance-label defect in the transcript header.

**Why this matters:**
The transcript header misidentifies which stage produced it, which can corrupt provenance auditing when transcripts are inspected out of context. Purely cosmetic for the math; non-zero for audit traceability.

**Required change:**
Change the banner string from `STAGE 152` to `STAGE 169` in both files (sympy line 30, mathematica line 31). Re-run both engines so the transcripts carry the correct header.

**Verification:**
Both transcript line-11 banners read `STAGE 169 — ...` after re-run; scripts exit 0.

## Independent-derivation check (Mathematica)

Transliteration, not independent. See F2. The `.wl` reproduces the SymPy script's exact variable names, the same `(dX . gMat . dY)/5` construction, the same `4*aX*aY + (4/5)*bX*bY` comparison target, the same `subsAxis` rule set, the same `Y20 = Sqrt[5/(16 Pi)](3 Cos[th]^2 - 1)` harmonic, the same `Series[Log[x0s(1+e y20Harm)],{e,0,2}]` expansion, the same tautological `epsPerp - xiPerp*iGrp`, and the same hardcoded `rNum`/`gNum` literals. The only differences are Mathematica syntax and `FullSimplify` vs `simplify`.

## Engine cross-check

Both engines agree on every asserted quantity: `I[x,y] - [...] = 0`; `a_x = eps*x1/4`, `b_x = 3*eps*x1/4`, `b_x - 3 a_x = 0`, `I_axis - 7/10 eps^2 x1 y1 = 0`; `<Y20> = 0`, `<Y20^2> = 1/(4 pi)`, linear coefficient `= 0`; `eps_perp - Xi_perp Igrp = 0`; and the printed numeric combination matches to the precision shown (`1.883732191180.../0.758035078944663/1.003143101138...`). No residual, sign, or factor disagreement. `engines_agree: true`.

## Verdict justification

The bottom-line paper claim — pure grouped real \(P_2\) anisotropy has no linear scalar feed-down, with feed-down starting quadratically — is substantively and non-tautologically verified. I attacked each of A1-A6 and they hold: I hand-expanded the grouped invariant identity over 25 (matches SymPy output term by term), verified the axisymmetric reduction gives \(a_x=\epsilon x_1/4\), \(b_x=3\epsilon x_1/4\) and \(\mathcal I_{\rm axis}=\frac{7}{10}\epsilon^2 x_1y_1\), and confirmed \(\langle Y_{20}\rangle=0\), \(\langle Y_{20}^2\rangle=\frac{1}{4\pi}\), and the vanishing linear coefficient \(-\frac{e^2}{8\pi}+\log X_0\) of the averaged log-observable. These exactly realize the monopole-selection and quadratic-invariant content of the notes/appendix. The paper alignment is `aligned`. The findings are: (F1) block 4's transport assertion is tautological and the load-bearing numeric Family-1 coefficients are printed but not asserted; (F2) the Mathematica engine is a transliteration of the SymPy engine; (F3) both scripts mislabel the banner as STAGE 152. None of these touches the verified core, so verdict is `findings`, not `stop_cold`. No `paper_misalignment`: every numeric the script displays matches the paper's stated values to the quoted precision.

## Self-test notes

Checked variable-independence (block 4 has no derivatives; the only `diff` is `D[logAvg, e]` at A6, and `logAvg` genuinely depends on `e` via the \(O(e^2)\) term, so the derivative is real, not identically zero). Checked parity for the sphere integrals: \((3\cos^2\theta-1)\sin\theta\) integrates to zero over \([0,\pi]\) (odd-about-equator weight kills the \(l=2\to l=0\) overlap), and \((3\cos^2\theta-1)^2\sin\theta\) integrates to \(8/5\) giving \(\langle Y_{20}^2\rangle=1/(4\pi)\) — both confirmed by hand. Trivial-case pre-check on the proposed F1 numeric assertion: \(g=0.758035\), \(g+1/(2s)=1.003143\), \(2g+3/(4s)=1.883732\) reproduce the paper's three stated constants, so the new `expect_close` assertions evaluate to nonzero literals matching the paper and pass at the quoted tolerance.
