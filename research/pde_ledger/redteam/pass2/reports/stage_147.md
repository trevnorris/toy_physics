---
unit_id: 147
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage147_first_order_rigidity_kernel.md
  paper_appendix: present
---

# Audit unit 147 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_147.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage147_first_order_rigidity_kernel.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (Family-1 mouth-correction block, lines 768-850, plus the `\input{stages/stage_147}` row at 1328)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.txt`

## What the paper claims

Stage 147 is a finite mouth-profile-corrections ledger step. The card's verification target (`stage_147.tex:16`) is:

> Traction shift is one kernel projection with coefficients \(A_T\approx-4.27264\), \(B_T\approx0.134875\).

The notes expand this into the full first-order rigidity law. On the lower compensated canonical branch (`R_q=1/4`, `\Sigma_0=\Pi/(1-S_q/4)`, `\widehat T_m=\sqrt{9\Sigma_0/20}`), the traction shift is `\delta\widehat T_m = \epsilon[A_T(\bar g_\varsigma-\mathfrak g_*)+B_T(\bar S_\varsigma-\mathcal S_*)]`, with the closed forms for `A_T` (notes:48-56) and `B_T` (notes:57-62) in terms of the starred quantities `T_{m,*},\Pi_*,\mathcal S_*,\mathcal S'_*,\mathfrak g'_*`. The deliverable numeric values are (notes:66-69, appendix:846-848):

> `A_T \approx -4.27263956256927`, `B_T \approx 0.134875005736706`,

plus the computed cross-check `|A_T|/B_T \approx 31.6785` (notes:77 — present in notes, not in the card). Section 2 of the notes elevates this to the single-kernel rigidity statement `\delta\widehat T_m = \epsilon\int_0^1 \mathcal W_*(x)[\varsigma(x)-\Sigma_*(x)]\,dx` with the centered weight `\mathcal W_*(x)=A_T(c(x)-\mathfrak g_*)+B_T(K_q(x)-\mathcal S_*)`, valid because `\Sigma_\epsilon-\Sigma_*` integrates to zero. Supporting starred constants quoted in the appendix: `\Pi_*\approx1.50882951349316` (663), `\Sigma_0(\Pi_*)\approx1.80594111095636` (773), `\widehat T_m(\Pi_*)\approx0.901484054174205` (775), `\mathfrak g'_*\approx0.0714453558083195` (831), with `c(x)=\cos(\pi x/2)` and `K_q(x)=\cosh(\tfrac\pi2(1-x))/\cosh(\pi/2)` (821-823).

## What the script claims to verify

Both engines solve `gPi(\Pi)=g_-` for `\Pi_*` (with `g_-` the lower-branch root built from the Family-1 anchor `rF1=\sqrt{12(37/20)^2/\pi^2-1}`), evaluate `g_*,S_*,g'_*,S'_*,\Sigma_*,T_*`, and assemble the closed-form `A_T`, `B_T`. They then assert: (1) `A_T`/`B_T` match the paper literals to `1e-12` and `|A_T|/B_T` matches `31.6785` to `1e-3`; (2) the closed-form `A_T` agrees with an independent autodiff route `A_T=-(dT_m/d\Pi)/(dg/d\Pi)` at `\Pi_*` (residual `<1e-20`); (3) the kernel projection `\int W_*(\varsigma-\Sigma_*)dx` reproduces the two-moment formula for a concrete non-canonical deformation `\varsigma(x)=2x` (residual `<1e-22`); (4) the kernel is source-centered, `\int_0^1\Sigma_*W_*\,dx=0` (residual `<1e-22`); (5) `g_*,S_*` equal their source-moment quadratures `\int\Sigma_* c` and `\int\Sigma_* K_q` (residual `<1e-25`). The SymPy script also prints the symbolic `\delta T_m` and `W_*(x)`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `A_T \approx -4.27263956256927` | py:56-58 / wl:67-68 assert vs literal `<1e-12` | match |
| `B_T \approx 0.134875005736706` | py:59-61 / wl:69-70 assert vs literal `<1e-12` | match |
| `|A_T|/B_T \approx 31.6785` (computed cross-check, notes-only) | py:62-64 / wl:71-72 assert vs `31.6785` `<1e-3` | match |
| Closed forms for `A_T`,`B_T` (in starred quantities) | py:33-42 / wl:49-53 build them; py:66-80 / wl:74-83 cross-check `A_T` against autodiff route | match |
| Single-kernel projection `\delta T_m=\epsilon\int W_*(\varsigma-\Sigma_*)dx` | py:93-116 / wl:98-112 projection identity (kernel quadrature == two-moment formula) | match |
| Centering / `\Sigma_\epsilon-\Sigma_*` integrates to zero | py:117-128 / wl:114-129 orthogonality `\int\Sigma_*W_*=0`; py:104-108/wl:104-106 normalization | match |
| `g_*,S_*` as source moments of `c`,`K_q` | py:130-143 / wl:131-141 quadrature vs closed forms | match |
| `\Pi_*,\Sigma_*,T_*` starred constants (appendix 663/773/775) | computed + printed by both engines (py:45-47/wl:56-58) | match (values agree to quoted precision) |
| `\mathfrak g'_*\approx0.0714...` (appendix:831) | computed internally (`gp_star`/`gPrimeStar`) but never printed/asserted | n/a (internal intermediate, carried by appendix) |

Dominant pattern: every paper-side deliverable is exercised by a non-tautological assertion in both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 56-58 | `abs(N(AT) - AT_paper) < 1e-12` | `A_T` literal | yes |
| A2 | sympy | 59-61 | `abs(N(BT) - BT_paper) < 1e-12` | `B_T` literal | yes |
| A3 | sympy | 62-64 | `abs(N(Abs(AT)/BT) - 31.6785) < 1e-3` | `|A_T|/B_T` cross-check | yes |
| A4 | sympy | 78-79 | `abs(AT_autodiff - AT_30) < 1e-20` | closed-form `A_T` (independent route) | yes |
| A5 | sympy | 107-108 | `abs(norm_s-1)<1e-30 and abs(norm_Sigma-1)<1e-30` | deformation normalization (centering precondition) | yes |
| A6 | sympy | 114-115 | `abs(lhs_proj - rhs_moment) < 1e-22` | single-kernel projection identity | yes |
| A7 | sympy | 126-127 | `abs(center_resid) < 1e-22` (`\int\Sigma_*W_*=0`) | source-centering of kernel | yes |
| A8 | sympy | 139-142 | `abs(g_star_moment-g_star)<1e-25`, `abs(S_star_moment-S_star)<1e-25` | `g_*,S_*` source moments | yes |
| B1 | mathematica | 67-68 | `expectZero[If[Abs[aT-aTPaper]<1e-12,0,...]]` | `A_T` literal | yes |
| B2 | mathematica | 69-70 | `expectZero[If[Abs[bT-bTPaper]<1e-12,0,...]]` | `B_T` literal | yes |
| B3 | mathematica | 71-72 | `expectZero[If[Abs[Abs[aT]/bT-31.6785]<1e-3,0,...]]` | `|A_T|/B_T` cross-check | yes |
| B4 | mathematica | 82-83 | `expectZero[If[Abs[aTAutodiff-aT]<1e-20,0,...]]` | closed-form `A_T` (independent route) | yes |
| B5 | mathematica | 95-96 | `expectZero[Chop[D[wStar-(aT c + bT kq), x]]]` | `W_*` offset is constant in x | yes (structural, weak; see ID-check) |
| B6 | mathematica | 106 | `expectZero[normalization]` | deformation normalization | yes |
| B7 | mathematica | 111-112 | `expectZero[If[Abs[lhsProj-rhsMoment]<1e-20,0,...]]` | single-kernel projection identity | yes |
| B8 | mathematica | 128-129 | `expectZero[If[Abs[centerResid]<1e-20,0,...]]` | source-centering of kernel | yes |
| B9 | mathematica | 138-141 | `expectZero[g_*/S_* vs moment integral <1e-20]` | `g_*,S_*` source moments | yes |

No tautological or orphaned assertions. Every row traces to a paper deliverable.

## Findings

None. The three first-pass findings (F1 SymPy had zero assertions; F2 the lone Mathematica check `R_q(g_minus)-1/4` was tautological; F3 the centering/projection structure was not exercised) have all been remediated in the current scripts: the SymPy script now carries eight substantive assertions, the tautological `R_q` check is removed, and both engines now verify the projection identity, the source-centering orthogonality, and the moment integrals. The stale `STAGE 130` banner flagged at first-pass `wl:188` is also fixed — `wl:26` now reads `STAGE 147`.

## Independent-derivation check (Mathematica)

A `.wl` exists; this is a dual-engine stage. Verdict: **PARTIAL** — same derivation *route* in both engines, but each load-bearing value is independently *recomputed* by a different engine's primitives, and within each engine the closed form is cross-checked against a genuinely different (autodiff) route. This is the established and acceptable dual-engine pattern for this ledger; it is NOT a `mathematica_transliteration` finding.

Corresponding sections:

- Closed-form assembly is identical choreography. SymPy py:33-42:
  > `AT = sp.N(-(sp.Rational(9,1)/(40*T_star)) * (1/(gp_star*(1-S_star/4)) + Pi_star*Sp_star/(4*gp_star*(1-S_star/4)**2)), 30)`

  Mathematica wl:49-52:
  > `aT = N[-(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)), 30];`

  Same hand-written closed form in both — but it is not left unchecked: each engine independently re-derives `A_T` by autodiff and asserts agreement.

- Autodiff route is computed independently in each CAS. SymPy py:73-79:
  > `Tm_of_Pi = sp.sqrt(sp.Rational(9, 20) * (Pi / (1 - Sformula/4)))` … `AT_autodiff = sp.N(-(dTm_dPi.subs(Pi, Pi_star)) / (dg_dPi.subs(Pi, Pi_star)), 30)`

  Mathematica wl:78-81:
  > `tmOfP = Sqrt[(9/20)*(p/(1 - sFormula/4))]` … `aTAutodiff = N[-(dTmDp /. p -> pStar)/(dgDp /. p -> pStar), 30];`

  I verified by hand that `-(dT_m/d\Pi)/(dg/d\Pi)` expands to exactly the closed form (the `1/(2\sqrt u)` prefactor reduces to `9/(40 T_{m,*})` since `T_{m,*}=\sqrt{9/20}\sqrt u`, and `du/d\Pi = 1/(1-S/4) + \Pi S'/(4(1-S/4)^2)`). So A4/B4 are valid, non-tautological independence checks: a hand sign/power slip in the closed form would fail them.

- The integral checks use genuinely different integrators: SymPy `sp.integrate` (symbolic→numeric) vs Mathematica `NIntegrate` (numerical quadrature). E.g. source-centering SymPy py:125 `sp.integrate(Sigma_star_x * Wstar_x, (x, 0, 1))` vs Mathematica wl:120-127 `NIntegrate[Evaluate[SetPrecision[sigmaStarX*wStar, 60]], ..., WorkingPrecision -> 60]`. These do not echo each other's algebra.

B5 (`D[wStar-(aT c+bT kq), x]==0`) is structurally guaranteed (the offset is a constant by construction) and therefore weak on its own, but the script's own comment (wl:115-119) acknowledges this and adds the substantive orthogonality check B8 that actually tests the centering constants. SymPy reaches the same conclusion via B7/A7 directly without the weak intermediate check. Not a finding.

## Engine cross-check

Both engines agree on all printed values to ~17-30 significant digits (precision differences are from the differing `WorkingPrecision`/`nsolve` settings, not disagreement):

| Quantity | SymPy (output) | Mathematica (output) |
|---|---|---|
| `Pi_*` | 1.508829513493155527… | 1.5088295134931555830… |
| `Sigma_*` | 1.805941110956353804… | 1.8059411109563538723… |
| `T_*` | 0.901484054174204022… | 0.9014840541742040389… |
| `A_T` | -4.27263956256927466681875206233 | -4.2726395625692746412223747285… |
| `B_T` | 0.134875005736705540316968999027 | 0.1348750057367055429617462273… |
| `|A_T|/B_T` | 31.678512554876561144 | 31.67851255487656033274… |

All assertions in both transcripts print `PASS`. The Mathematica `NIntegrate::precw` precision warnings (output lines 31, 38) are advisory only; each enclosing `expectZero` still reports `= 0` and PASSes, and the same residuals are confirmed symbolically by SymPy. `engines_agree: true`.

## Verdict justification

`clean`. I read the card, notes (both sections), and the part-04 appendix block (768-850) first, then attacked the scripts. Attacks tried and failed: (1) checked the autodiff route is mathematically equivalent to the closed form rather than a restatement of it — it is a genuine second route that would catch a sign/power error; (2) checked the source-centering `\int\Sigma_*W_*=0` is not tautological by construction — `g_*,S_*` come from the closed forms while the integral is independent quadrature, so a closed-form error fails it; (3) checked the moment-integral check catches a transcription error in `gPi`/`Sformula`; (4) verified the Family-1 anchor uses `100\pi^2` (`12(37/20)^2=4107/100`), confirming the stale `168\pi^2` warning does NOT apply (output line 120 shows `-100\pi^2`); (5) confirmed `\Sigma_*,T_*,\Pi_*` match the appendix literals to quoted precision; (6) confirmed `|A_T|/B_T=31.6785` is correctly treated as a notes-side computed cross-check (asserted at `1e-3`), not flagged MISSING merely because the card omits it; (7) confirmed all stage self-labels read 147 (no stale-renumber labels). Every paper deliverable is exercised by a non-tautological, paper-anchored assertion in both engines, and the two engines agree. Outputs are fresh (both `.txt` newer than their scripts).

## Self-test notes

Variable independence: the autodiff `D[..., p]`/`sp.diff(..., Pi)` act on `Tm_of_Pi`/`gPi`, which genuinely depend on `Pi`/`p` (via `Sformula`/`gFormula`); derivatives are not identically zero. Symmetry/parity: the projection and centering integrals are over the bounded `[0,1]` (not symmetric/unbounded), so no parity-cancellation trap — the residuals are genuine. Trivial-case: the source-centering integral is nonzero unless the `-g_*`/`-S_*` constants are exactly right, so the `assert_zero` is substantive. No directive written (no script-side findings; nothing for Codex to apply).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 8 values checked, 0 misaligned

Deliverable-level reconciliation (every RESULT/labeled value the scripts emit):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `A_T = -4.27263956256927…` | py:48 / wl:59; sympy out:8, math out:12 | notes:67; appendix:846; card `\approx-4.27264` (stage_147.tex:16) | MATCH |
| `B_T = 0.134875005736706…` | py:49 / wl:60; sympy out:9, math out:13 | notes:69; appendix:848; card `\approx0.134875` (stage_147.tex:16) | MATCH |
| `|A_T|/B_T = 31.6785…` (computed cross-check) | py:50 / wl:61; sympy out:10, math out:14 | notes:77 (`\approx 31.6785`) | MATCH (notes-only; correctly not a card literal) |
| `Pi_* = 1.50882951349316…` | py:45 / wl:56; sympy out:5, math out:9 | appendix:663 (`\Pi_*\approx1.50882951349316`) | MATCH |
| `Sigma_* = 1.80594111095636…` | py:46 / wl:57; sympy out:6, math out:10 | appendix:773 (`\Sigma_0(\Pi_*)\approx1.80594111095636`) | MATCH |
| `T_* = 0.901484054174205…` | py:47 / wl:58; sympy out:7, math out:11 | appendix:775 (`\widehat T_m(\Pi_*)\approx0.901484054174205`) | MATCH |
| `delta T_m = eps*(A_T(gbar-g_*)+B_T(Sbar-S_*))` (symbolic) | py:84-86 / wl:85-86; sympy out:15-62, math out:23 | notes:36-45; appendix:835-842 | MATCH (symbolic form agrees) |
| `W_*(x) = A_T(c-g_*)+B_T(K_q-S_*)` (centered kernel) | py:89-91 / wl:88-89; sympy out:63-125, math out:24 | notes:98-106; appendix:850 (description) | MATCH (symbolic form agrees) |

INTERNAL items (genuine scaffolding, raise no finding):
`g_*`, `S_*`, `g'_*`/`gp_star`/`gPrimeStar` (=0.0714453558083195, carried in appendix:831 but emitted by the script only as an internal intermediate, never printed/asserted), `S'_*`/`Sp_star`, `gminus`/`gMinus`, `rF1`, the centering-offset constant `3.15005267…` (=`-(A_T g_* + B_T S_*)`, the constant term of `delta T_m`/`W_*`), the test deformation `\varsigma(x)=2x`, normalization residuals, projection/centering/moment residuals, and all PASS flags.

All eight emitted deliverable values reconcile to the `.tex` card and/or `.md` notes at the quoted precision. The `g'_*` appendix literal (0.0714453558083195) is internally consistent with the computation but is not a script-emitted labeled result, so it is INTERNAL, not a missing deliverable.
