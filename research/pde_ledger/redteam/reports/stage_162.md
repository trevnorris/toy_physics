---
unit_id: 162
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage162_parent_compensation_rigidity.md]
  paper_appendix: present
---

# Audit unit 162 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_162.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage162_parent_compensation_rigidity.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (just an `\input{stages/stage_162}` line at L1358)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage162_parent_compensation_rigidity_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage162_parent_compensation_rigidity_mathematica_audit.txt`

## What the paper claims

The card's quoted bottom-line claim is: "Staying on the exact parent compensation family gives automatic D/N similarity preservation and \(\Delta_Q=0\) at first order." The Checks list mandates: (a) deviations about the renormalized canonical point, (b) even-preservation constraint imposed before reading the residual odd defect, (c) tangent motion on the parent compensation family yields \(\delta_\perp=0\). The notes flesh out three loadbearing equations: (1) the similarity identity \(\delta\ln\gamma_0 - 2\,\delta\ln(L_W/a) = 0\) along the exact parent family (i.e., \(\Xi_{\rm slip}=0\)), using \(\gamma_0=(1+\mathfrak r^2)/9\) and \(L_W/a=(\pi/2)\sqrt{(1+\mathfrak r^2)/3}\); (2) the lower-branch differential law \(\delta\mathfrak g = (1 - \mathfrak r/(2\sqrt{1+\mathfrak r^2}))\,\delta\mathfrak r\) from \(\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2}\); (3) the exact positivity decomposition \(1 - \mathfrak r/(2\sqrt{1+\mathfrak r^2}) = (4+3\mathfrak r^2)/(2\sqrt{1+\mathfrak r^2}(2\sqrt{1+\mathfrak r^2}+\mathfrak r))>0\), forcing \(\delta\mathfrak g=0 \Rightarrow \delta\mathfrak r=0\). The notes also state a numerical Family-1 slope \(\approx 0.564199521046343\) at \(\mathfrak r_{F1}\approx 1.77799353547498\).

## What the script claims to verify

Both scripts (SymPy and Mathematica) verify four items: (1) the similarity identity `dlog gamma0 - 2 dlog Lratio == 0` from independently posed `gamma0 = (1+r^2)/9` and `Lratio = pi/2 sqrt((1+r^2)/3)`; (2) the lower-branch differential law `slope - (1 - r/(2 sqrt(1+r^2))) == 0` from differentiating `g_lower = r - sqrt(1+r^2)/2`; (3) the positive-slope decomposition `slope - (4+3 r^2)/(2 sqrt(1+r^2)(2 sqrt(1+r^2)+r)) == 0`; (4) numerical evaluation of the slope and its reciprocal at the carry-forward `r_F1 = 1.77799353547498`. The docstrings on both scripts say exactly these four checks.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Similarity identity \(\delta\ln\gamma_0 - 2\delta\ln(L_W/a)=0\) along family | A1 sympy L53 / A1m wl L43 `expectZero(... dlog_gamma - 2*dlog_L)` | match |
| Lower-branch differential law \(\delta\mathfrak g = (1-r/(2\sqrt{1+r^2}))\,\delta\mathfrak r\) | A2 sympy L58–61 / A2m wl L47 `expectZero("lower-branch differential law", slope - (1 - r/(2 sqrt(1+r^2))))` | match |
| Exact positive-slope decomposition \((4+3r^2)/(2\sqrt{1+r^2}(2\sqrt{1+r^2}+r))\) | A3 sympy L67 / A3m wl L50 `expectZero("positive slope decomposition", slope - slope_pos)` | match |
| Numerical Family-1 slope \(\approx 0.564199521046343\) | sympy L71–76 / wl L53–58 prints `slope_num`, `inv_slope_num` | match (informational) |
| \(\Xi_{\rm slip}=0\), \(\delta\mathfrak r=0\Rightarrow\delta\mathfrak B_W=\delta\gamma_W=\Delta_Q=0\) | Logical consequences once (1)–(3) hold; explicitly summarized in the script's carry-forward print block | match (carry-forward) |

All checked rows are `match`. Set `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero("similarity identity", dlog_gamma - 2*dlog_L)` (residual `simplify`/`expand` then compare `== 0`) | claim 1 (Ξ_slip=0 along family) | yes |
| A2 | sympy | 58–61 | `expect_zero("lower-branch differential law", slope - (1 - r/(2*sqrt(1+r^2))))` | claim 2 (δg in terms of δr) | yes |
| A3 | sympy | 67 | `expect_zero("positive slope decomposition", slope - slope_pos)` | claim 3 (positivity → rigidity) | yes |
| A1m | mathematica | 43 | `expectZero["similarity identity", dlogGamma - 2*dlogL]` (`FullSimplify[Together[Expand[...]]]` to 0) | claim 1 | yes |
| A2m | mathematica | 47 | `expectZero["lower-branch differential law", slope - (1 - r/(2*Sqrt[1+r^2]))]` | claim 2 | yes |
| A3m | mathematica | 50 | `expectZero["positive slope decomposition", slope - slopePos]` | claim 3 | yes |

All six rows are anchored. None of A1/A2/A3 is tautological: each compares the symbolic output of a `diff`/`D[]` on an independently posed parent-family formula against a separately-written target. The "positive slope decomposition" check, in particular, is a real algebraic identity (multiplying numerator and denominator of \(1 - r/(2\sqrt{1+r^2})\) by \(2\sqrt{1+r^2}+r\)) that the script does not pre-bake.

## Findings

None.

I attempted the following attacks:

- **Tautology on A1**: `gamma0` and `Lratio` are independently defined; `dlog_gamma = 2 r/(1+r^2) dr`, `dlog_L = r/(1+r^2) dr`. The identity holds because both reduce to the same `r/(1+r^2)` (up to factor 2), but this is a derived algebraic consequence of the two physical inputs, not a tautology of construction.
- **Tautology on A2**: `slope = diff(g_lower, r)` is computed by SymPy/Mathematica, not declared; comparing to `1 - r/(2 sqrt(1+r^2))` is a real check that the script's hand-written target matches the engine's autodiff output.
- **Tautology on A3**: `slope_pos` is `(4+3 r^2)/(2 sqrt(1+r^2)(2 sqrt(1+r^2)+r))`, an algebraically rewritten form. The check is non-trivial; SymPy must simplify the difference. By hand: \(1-r/(2\sqrt{1+r^2})=(2\sqrt{1+r^2}-r)/(2\sqrt{1+r^2})\); rationalize by \((2\sqrt{1+r^2}+r)\) → numerator \(4(1+r^2)-r^2 = 4+3r^2\), denominator \(2\sqrt{1+r^2}(2\sqrt{1+r^2}+r)\). Identity confirmed by hand; script's `expect_zero` is a real check.
- **Branch / sign attacks**: `r` is declared `real=True` (SymPy) and `Element[{r,dr}, Reals]` (Mathematica). \(\sqrt{1+r^2}>0\) for all real `r`. \(2\sqrt{1+r^2}+r > 0\) for all real `r` (since \(4(1+r^2)>r^2\)). No branch ambiguity in the simplifications.
- **Hardcoded numeric constant `r_F1 = 1.77799353547498`**: used only in informational prints of the slope and its reciprocal — no assertion depends on it. The notes box names this exact Family-1 value as a carry-forward; flagging it as `hardcoded_result` is not warranted because it does not load-bear any assertion.
- **Engine cross-check**: SymPy and Mathematica produce the same simplified residuals (all 0) and matching numerical values to ~28 sig figs (sympy `(dg/dr)|_F1 = 0.564199521046342514375259007973`, mathematica `0.564199521046342510304460143977…`). Agreement is at the printed precision.
- **Banner label "STAGE 145"** (sympy L34, mathematica L26): cosmetic copy-paste leftover; not a math finding. The script docstrings and trailing print statements correctly identify "Stage 162", and the file paths are also correct. Not raised as a finding since it doesn't affect any assertion or claim.
- **Source-stage comment "Stages 99 and 102"** (sympy L39): the notes attribute the formulas to "Stages 99 and 170". This is a comment-only reference; the script's math does not depend on it. Not raised as a finding since the formulas \(\gamma_0=(1+r^2)/9\) and \(L_W/a=(\pi/2)\sqrt{(1+r^2)/3}\) are correct per the notes regardless of which upstream stage label is attached.

## Independent-derivation check (Mathematica)

The two scripts are structurally parallel because the math is short: each defines `gamma0`, `L_W/a`, `g_lower`; differentiates with the engine's primitive (`sp.diff` vs `D[]`); simplifies (`sp.simplify(sp.expand(...))` vs `FullSimplify[Together[Expand[...]]]`); and compares to a target expression. There is no large algorithmic choreography being shadowed — the parallelism is the irreducible minimum any second-engine implementation of "differentiate and simplify" would have. The Mathematica simplifier reaches `1 - r/(2 Sqrt[1 + r^2])` directly; SymPy's prints `-r/(2*sqrt(r**2+1)) + 1` (commuted form). Both arrive at residual `0` for all three checks. This is acceptable engine independence; not `mathematica_transliteration`.

## Engine cross-check

Both engines produce residual `0` on all three assertions. SymPy prints (output L18–21):
- `similarity identity = 0`
- `lower-branch differential law = 0`
- `positive slope decomposition = 0`

Mathematica prints (output L18–24):
- `similarity identity = 0` then `PASS: similarity identity`
- `lower-branch differential law = 0` then `PASS: lower-branch differential law`
- `positive slope decomposition = 0` then `PASS: positive slope decomposition`

Numerical Family-1 slope: SymPy `0.564199521046342514375259007973`, Mathematica `0.564199521046342510304460143977…`. Reciprocal: SymPy `1.77242263188284677523314446134`, Mathematica `1.77242263188284678802148576124…`. Agreement to ~28 sig figs is consistent with both engines using 30-digit precision and rounding differently in the last few digits. Engines agree.

Saved outputs are newer than their source scripts (mtime check: sympy script 1778522211 vs output 1778525232; mathematica 1778522212 vs output 1778527108) — no `stale_output`.

## Verdict justification

The paper's claim (Ξ_slip=0 along the exact parent compensation family, plus rigidity δ𝔤=0 ⇒ δ𝔯=0 on the lower branch) decomposes into three concrete algebraic identities. Both engines verify all three non-tautologically. The numerical Family-1 slope is reported informationally and matches across engines. Outputs are fresh. Comment-level minor blemishes (banner label says "STAGE 145"; source-stage comment names "Stages 99 and 102" vs notes' "Stages 99 and 170") are not math findings — neither affects any assertion. No attack landed. Verdict: `clean`.

## Self-test notes

I checked variable independence (each `diff(EXPR, r)` operates on `EXPR` that genuinely depends on `r`), parity/symmetry (no integrals to worry about — this stage is pure algebra), trivial-case behavior (`r=0`: slope = `1 - 0 = 1`, slope_pos = `4 / (2·1·(2·1+0)) = 1`, identity holds; `r=1`: slope = `1 - 1/(2√2) ≈ 0.6464`, slope_pos = `7/(2√2(2√2+1)) ≈ 0.6464`, agrees), and paper round-trip (the script's load-bearing formulas \(\gamma_0=(1+r^2)/9\), \(L_W/a=(\pi/2)\sqrt{(1+r^2)/3}\), \(g_-=r-\frac12\sqrt{1+r^2}\) match the notes' boxed equations verbatim). No issues uncovered.
