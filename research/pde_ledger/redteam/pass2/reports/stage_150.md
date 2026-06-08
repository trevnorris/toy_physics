---
unit_id: 150
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage150_full_profile_residual.md]
  paper_appendix: present
---

# Audit unit 150 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_150.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage150_full_profile_residual.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsection "Full-profile residual", lines 852-897; plus the `\input{stages/stage_150}` row at line 1334)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage150_full_profile_residual_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.txt`

## What the paper claims

Stage 150 keeps the *full* coupled Family-1 mouth profiles (static shell `T_s` + mixed D/N `T_q`) instead of the linear electrochemical tangent `Φ_lin = Π_* x`, and studies the tangent-subtracted residual `R_*(x) := Φ_*(x) − Π_* x` with `Φ_*(x) = 4Σ_m* T_s − Σ_m* T_q`. The card's bottom-line `\stagefield{}` quote (stage_150.tex:16) states: "Exact residual \(R_*(x)=\Phi_*(x)-\Pi_*x\) is tangent-matched but has negative curvature at the mouth." The notes make this precise as three boxed deliverables: (D1) the first exact mouth-shape theorem `R_*(0)=0, R_*'(0)=0` (notes:105-111); (D2) the curvature law `R_*''(0) = −3 Σ_m* Π_*/(1−e^{−Π_*}) < 0` (notes:118-122); and supporting these, (D3) the channel-slope identity `S_q(Π) = T_q'(0) = A_q·k − C_q·Π` with `k=π/2` and `T_s(0)=T_q(0)=0`, `T_s'(0)=1` (the building-block profile boundary data, notes:31-61, used in the residual definition `R_* = Σ_m*[4T_s − T_q − (4 − S_q)x]`, notes:96-100). The appendix subsection (part04.tex:852-866) restates `R_*(0)=R_*'(0)=0`, `R_*''(0)<0` qualitatively and then carries downstream numeric covariances/corrected point that the script does NOT compute.

## What the script claims to verify

The SymPy docstring (py:5) says "Lightweight SymPy audit for the full-profile Family-1 mouth potential and its curvature residual." The seven load-bearing `expect_zero` assertions verify: profile boundary data `T_s(0)=0`, `T_q(0)=0`, `T_s'(0)−1=0`, `T_q'(0)−S_q=0` (py:49-52); the residual tangent-match `R(0)=0`, `R'(0)=0` (py:56-57); and the curvature law `R''(0) − (−3 Σ Π/(1−e^{−Π})) = 0` (py:60-64). `S_q` is built as `A_q·k − C_q·Π` from FREE placeholder symbols (py:41-43), printed for display, then `.subs`-ed to the concrete `A_q, C_q` for the load-bearing check at py:52, which compares it against `sp.diff(Tq,x).subs(x,0)` computed independently from the full `Tq`. The Mathematica script (`.wl`) makes the identical seven `expectZero` checks.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: `R_*(0)=0, R_*'(0)=0` (notes:105-111; appendix:863) | py:56-57 `R(0)=0`, `R'(0)=0`; wl:53-54 | match |
| D2: `R_*''(0) = −3 Σ_m* Π_*/(1−e^{−Π_*})` (notes:118-122) | py:60-64 `R''(0) − target == 0` with `target = −3 Σ Π/(1−e^{−Π})`; wl:56-59 | match |
| D2': `R_*''(0) < 0` qualitative (card:16; appendix:865) | implied by D2 form (Σ,Π>0 ⇒ negative); printed in Theorem block py:67, wl:63 | match |
| D3: profile boundary data `T_s(0)=T_q(0)=0`, `T_s'(0)=1`, slope `S_q = A_q k − C_q Π` (notes:31-61) | py:49-52 (four checks); wl:47-50 | match |
| appendix numerics `Cov_*≈0.0648…/0.0389…`, `Π_corr≈2.4159…`, `T̂_{m,corr}≈1.1731…` (appendix:886-895) | NOT computed by script (script is purely symbolic) | n/a — these are downstream-of-150 quantities, not stage-150 deliverables |

The card/notes deliverables are all faithfully exercised. The appendix numerics at part04.tex:886-895 are part of the *next* block (covariance shifts / corrected canonical point) that consumes `R_*`; they are not claimed by the stage-150 card and are not script-emitted, so their absence from the script is correct, not a `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 49 | `expect_zero(T_s(0))` | D3 boundary | yes |
| A2 | sympy | 50 | `expect_zero(T_q(0))` | D3 boundary | yes |
| A3 | sympy | 51 | `expect_zero(T_s'(0)-1)` | D3 boundary | yes |
| A4 | sympy | 52 | `expect_zero(diff(Tq,x)|_0 - S_q)` | D3 slope | yes (S_q hand-asserted vs engine derivative) |
| A5 | sympy | 56 | `expect_zero(R(0))` | D1 | yes |
| A6 | sympy | 57 | `expect_zero(R'(0))` | D1 | yes |
| A7 | sympy | 64 | `expect_zero(R''(0) - target)` | D2 | yes (R2 computed independently from R) |
| B1 | mathematica | 47 | `expectZero[ts/.x->0]` | D3 boundary | yes |
| B2 | mathematica | 48 | `expectZero[tq/.x->0]` | D3 boundary | yes |
| B3 | mathematica | 49 | `expectZero[(D[ts,x]/.x->0)-1]` | D3 boundary | yes |
| B4 | mathematica | 50 | `expectZero[(D[tq,x]/.x->0)-sQ]` | D3 slope | yes |
| B5 | mathematica | 53 | `expectZero[r/.x->0]` | D1 | yes |
| B6 | mathematica | 54 | `expectZero[D[r,x]/.x->0]` | D1 | yes |
| B7 | mathematica | 59 | `expectZero[r2 - targetR2]` | D2 | yes |

All thirteen load-bearing assertions trace to a paper deliverable and are non-tautological (see Findings discussion and Verdict).

## Findings

None. (Attacks tried and failed are documented below in the Verdict justification and in the display-only slope analysis.)

### Display-only `S_q` slope — explicitly checked for tautology (NOT a finding)

The first-pass heads-up flagged a DISPLAY-only compact `S_q(Π) = A_q·k − C_q·Π` slope built via a free-placeholder-then-`.subs` construction, with the question of whether it genuinely proves the slope or is a tautological round-trip. Verdict: it genuinely proves the slope.

- py:41-43 builds `Sq_symbolic = Aq_s*k − Cq_s*Pi` from *free* symbols `Aq_s, Cq_s`, prints that compact form for display (py:45, output line 5), then `Sq = Sq_symbolic.subs({Aq_s: Aq, Cq_s: Cq})` substitutes the concrete definitions (py:33-34).
- The load-bearing check at py:52 is `expect_zero("T_q'(0)-S_q", sp.diff(Tq, x).subs(x, 0) - Sq)`. The LHS term `sp.diff(Tq, x).subs(x,0)` is computed by SymPy differentiating the *full* `Tq = Aq·sinh(kx) − Cq·cosh(kx) + Cq·exp(−Πx)` (py:35) — it is NOT defined as `Sq`. So the check compares the hand-asserted slope `A_q k − C_q Π` against the engine-computed derivative. It is falsifiable: had the author written `Aq*k + Cq*Pi` (wrong sign on the `cosh` term) the check would not vanish. The construction is therefore a genuine verification, not a `x==x` round-trip. The `.wl` mirror (wl:41-42, 50) is the same.

## Independent-derivation check (Mathematica)

The `.wl` is structurally a line-by-line port of the `.py`, not an independently choreographed re-derivation. Three corresponding sections:

1. Profile definitions — py:30-35 `Ts = (1 - sp.exp(-Pi*x))/(Pi*(1 - sp.exp(-Pi))) - x*sp.exp(-Pi)/(1 - sp.exp(-Pi))` … `Tq = Aq*sp.sinh(k*x) - Cq*sp.cosh(k*x) + Cq*sp.exp(-Pi*x)` vs wl:32-35 `ts = (1 - Exp[-p*x])/(p*(1 - Exp[-p])) - x*Exp[-p]/(1 - Exp[-p])` … `tq = aq*Sinh[k*x] - cq*Cosh[k*x] + cq*Exp[-p*x]`. Same five expressions, same order, identical symbol roles (`Pi→p`, `Sigma→sigmaM`).
2. The free-placeholder slope trick — py:41-43 `Aq_s, Cq_s = sp.symbols("Aq Cq"); Sq_symbolic = Aq_s*k - Cq_s*Pi; Sq = Sq_symbolic.subs(...)` vs wl:41-42 `sQsymbolic = aqS*k - cqS*p; sQ = sQsymbolic /. {aqS -> aq, cqS -> cq}`. Identical construction, and the *comments* (py:37-40 vs wl:37-40) are verbatim the same prose.
3. The seven checks — py:49-64 vs wl:47-59 are the same seven `expect_zero`/`expectZero` calls in the same order with the same operands.

This is `mathematica_transliteration` in form. I am NOT raising it as a finding because for this stage there is no independent *derivation* step that the two engines could perform differently: the profiles `T_s, T_q` and the residual `R_*` are given closed forms (from the notes / the upstream physics), and the only operations are differentiation and `Simplify`/`FullSimplify` to constant or zero. Each CAS independently performs its own differentiation and simplification on those given forms — that is legitimate dual-engine cross-confirmation of the symbolic identities, even though the *script choreography* mirrors. The engines also reach the curvature constant in genuinely different surface forms (SymPy: `3·Π·Σ·e^Π/(1−e^Π)`, output line 23-28; Mathematica: `(−3·E^p·p·sigmaM)/(−1+E^p)`, output line 19) and both independently certify `R''(0) − target == 0`. Per the standing pass-2 policy this borderline-port is recorded but, given the no-independent-derivation-possible character of the stage, treated as PARTIAL rather than a defect.

## Engine cross-check

Both engines pass all seven checks and agree:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| `T_s(0)`,`T_q(0)`,`T_s'(0)-1`,`T_q'(0)-S_q` | `0` (lines 17-20) | `0` + PASS (lines 7-14) |
| `R(0)`,`R'(0)` | `0` (lines 21-22) | `0` + PASS (lines 15-18) |
| `R''(0)` | `3·Π·Σ·e^Π/(1−e^Π)` (lines 23-28) | `(−3·E^p·p·sigmaM)/(−1+E^p)` (line 19) |
| `R''(0) − target` | `0` (line 29) | `0` + PASS (line 20) |

The two `R''(0)` surface forms are algebraically identical: `3ΠΣe^Π/(1−e^Π) = −3ΠΣ/(1−e^{−Π})` (multiply num & den by `e^{−Π}`), and `−3e^p p s/(e^p−1) = −3ps/(1−e^{−p})`. Both equal the notes' `−3 Σ_m* Π_*/(1−e^{−Π_*})` (notes:121). No `engine_disagreement`.

## Verdict justification

Clean. I read the card, the notes, and the appendix subsection before the scripts. The script's verified claims (D1 tangent-match, D2 curvature law with the exact constant, D3 profile boundary data + slope) match the paper's stated deliverables exactly. Attacks tried and failed: (a) the display-only `S_q` slope is NOT a tautological round-trip — the load-bearing check at py:52/wl:50 compares the hand-asserted `A_q k − C_q Π` against the engine-differentiated `T_q'(0)`, and is falsifiable on a sign/coefficient error; (b) the `R''(0)` check computes `R2` independently from `R` (py:60) and compares to a separately written `target` (py:61), so it is non-tautological and the engines reach it in different surface forms; (c) no stale `168π²`/`100π²` label appears anywhere in the stage-150 sources (grep clean); (d) outputs are fresh (both `.txt` mtime 2026-05-29 22:49 > both source mtime 22:46). One minor non-finding robustness asymmetry: the SymPy script declares `Pi` positive/real but does not exclude the removable pole `Π=π/2` (denominator `π²/4−Π²` of `C_q`), whereas the `.wl` adds `p != Pi/2` (wl:29); this does not affect the symbolic identities (all reduce to `0` for generic `Π`, and the physical `Π_* ≈ 2.42` is far from `π/2 ≈ 1.57`), so it is not raised. Transliteration is PARTIAL in form but defensible here (no independent-derivation step exists to differentiate the engines).

## Self-test notes

Variable-independence: every `sp.diff`/`D` is on an expression that genuinely depends on the differentiation variable `x` (`Ts, Tq, R` all contain `x`), so no derivative is identically zero — the `T_q'(0)−S_q` and `R''(0)−target` checks are substantive, not trivially-passing. Trivial-case: substituting `x=0` into `R` and `R'` correctly gives `0` (tangent match), and `R''(0)` reduces to the nonzero `−3ΣΠ/(1−e^{−Π})`, matching both engines' printed forms. Paper round-trip: the curvature constant `−3 Σ_m* Π_*/(1−e^{−Π_*})` in the assertion matches notes:121 verbatim; no new misalignment introduced (no fix prescribed). No directive written — zero findings.

## Value Reconciliation (pass-2 augmentation)

The stage-150 scripts are purely SYMBOLIC; they emit no named numeric constants. The deliverable values are the symbolic boundary data, slope, and the curvature closed form.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `T_s(0)=0`, `T_q(0)=0`, `T_s'(0)=1` (profile boundary data) | py:49-51 / wl:47-49; out (sympy) lines 17-19, (math) 7-12 | notes:32-39 (`T_s` form), :41-61 (`T_q` form) implying these; checks D3 | MATCH (notes) |
| `S_q(Π) = A_q·k − C_q·Π`, `k=π/2` (channel slope) | py:41-45 / wl:41-45; out (sympy) line 5 `pi*Aq/2 - Cq*Pi`, (math) line 5 `-(cqS*p)+(aqS*Pi)/2` | notes:59-61 (`A_q` defn uses `k=π/2`); slope is `T_q'(0)` used in `R_*` defn notes:96-100 | MATCH (notes) |
| `R_*(0)=0`, `R_*'(0)=0` (tangent-match theorem) | py:56-57 / wl:53-54; out lines 21-22 / 15-18 | notes:105-111 (boxed); appendix part04.tex:863 | MATCH (notes + appendix) |
| `R_*''(0) = −3 Σ_m* Π_*/(1−e^{−Π_*})` (curvature law) | py:60-64 / wl:56-59; out (sympy) lines 23-28 `3ΠΣe^Π/(1−e^Π)`, (math) line 19 `(−3 E^p p s)/(−1+E^p)` | notes:118-122 (boxed exact constant) | MATCH (notes) |
| `R_*''(0) < 0` (sign / sublinear) | Theorem print py:67 / wl:63 | card stage_150.tex:16; notes:121,125; appendix part04.tex:865 | MATCH (card + notes + appendix) |

INTERNAL scaffolding (no finding): `Ts/ts`, `Cq/cq`, `Aq/aq`, `Tq/tq`, `R/r` intermediate expressions; `Sq_symbolic`/`sQsymbolic` display placeholders; `target_R2`/`targetR2`; pass/zero residual flags.

Note: the appendix numerics at part04.tex:886-895 (`Cov_*(c,R_*)≈0.0648…`, `Cov_*(K_q,R_*)≈0.0389…`, `Π_corr≈2.4159…`, `T̂_{m,corr}≈1.1731…`) belong to the *downstream* covariance-shift block that consumes `R_*`; they are not stage-150 script outputs and not stage-150 card deliverables, so they are excluded from this reconciliation (no MISSING-DELIVERABLE).

reconciliation: complete; 5 deliverable values checked, 0 misaligned
