---
unit_id: 148
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage148_representative_positive_families.md]
  paper_appendix: present
---

# Audit unit 148 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_148.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage148_representative_positive_families.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read the Family-1 mouth-correction section, lines ~540-960, that carries this stage's deliverables)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage148_representative_positive_families_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_148.tex:7,15-17`) is terse: it is a "finite mouth-profile corrections ledger step" whose audit target is the quoted block "Uniform and self-matched derivative families bracket the first-order correction; interpolation reproduces the earlier compensation fraction." The notes (`notes/stages/...148...md`) are authoritative on intent and enumerate the deliverables: (1) the uniform-broadening family `ς_u(x)=1` with `ḡ_u=2/π`, `S̄_u=2tanh(π/2)/π` and boxed shifts `δΠ_u/ε≈+1.699414961314297`, `δT̂_{m,u}/ε≈+0.508756302215084` (notes:42-46); (2) the self-matched derivative family `ς_d(x)=(π/2)cos(πx/2)` with `ḡ_d=π/4`, `S̄_d≈0.657844575502831` and boxed shifts `δΠ_d/ε≈-0.382993186095928`, `δT̂_{m,d}/ε≈-0.116943802151811` (notes:74-79); (3) the convex interpolation giving affine `δΠ_λ/ε=1.699414961314297-2.082408147410224λ`, `δT̂_{m,λ}/ε=0.508756302215084-0.625700104366895λ` (notes:98-113), the bias-neutral point `λ_{Π,0}≈0.816081594488460` (notes:119), traction-neutral `λ_{T,0}≈0.813099276577333` (notes:141), and the key consistency claim `1-λ_{Π,0}≈0.183918405511540` agreeing with the Stage-126 closed form `ξ_*` (notes:124-135). The first-order shift machinery (`δΠ=-ε(ḡ_ς-g_*)/g'_*`, `δT̂_m=ε[A_T(ḡ_ς-g_*)+B_T(S̄_ς-S_*)]`) and the published literals `A_T≈-4.27263956256927`, `B_T≈0.134875005736706` live in appendix part04 (`stage_appendix_part04.tex:826-849`).

## What the script claims to verify

The SymPy script derives the starred canonical point `Pi_star` from `nsolve(gPi - gminus)` (py:24), computes `g_*`, `S_*`, `g'_*`, `S'_*`, `T_*` from the explicit exponential-family formulas, builds `A_T`/`B_T` from a hand-written closed form that includes the `Pi_star*Sp_star` S-follows-Π chain term (py:32-41), anchors `A_T`/`B_T` to the published appendix literals via `assert abs(...)<1e-11` (py:44-45), then computes and prints both families' `δΠ/ε`, `δT/ε`, the affine interpolation slopes, the neutral points, and asserts the EXACT symbolic residual `(1-λ_{Π,0}) - ξ_*_closed == 0` (py:109). The Mathematica script does the same physical setup but derives `aT` via its OWN `D[Tm[p],p]` autodiff along the S=sFormula(p) curve (wl:47-49), independently anchors `aT`/`bT` to the same appendix literals (wl:57-60), and independently asserts `(1-λ_{Π,0})-ξ_*_closed === 0` via `FullSimplify` (wl:98).

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `ḡ_u=2/π`, `S̄_u=2tanh(π/2)/π` (notes:42-46) | `g_u`, `S_u` printed (py:50-51, out:5-6) | match |
| `δΠ_u/ε≈+1.699414961314297` (notes:44) | `dPi_u` printed = 1.69941496131429... (py:53, out:7) | match |
| `δT̂_{m,u}/ε≈+0.508756302215084` (notes:45) | `dT_u` printed = 0.50875630221508... (py:54, out:8) | match |
| `δΠ_d/ε≈-0.382993186095928` (notes:76) | `dPi_d` printed = -0.38299318609592... (py:65, out:11) | match |
| `δT̂_{m,d}/ε≈-0.116943802151811` (notes:78) | `dT_d` printed = -0.11694380215181... (py:66, out:12) | match |
| affine slopes -2.082.../-0.6257... (notes:103,112) | `dPi_lam`,`dT_lam` expanded (py:80,82, out:14,18) | match |
| `λ_{Π,0}≈0.816081594488460` (notes:119) | `lam_Pi_zero` (py:84, out:21) | match |
| `λ_{T,0}≈0.813099276577333` (notes:141) | `lam_T_zero` (py:85, out:22) | match |
| `1-λ_{Π,0}=ξ_*` exact consistency (notes:124-135) | `assert exact_resid==0` (py:109, out:25) | match |
| `A_T≈-4.27263956256927` (app:846) | `assert abs(AT - (-4.27...))<1e-11` (py:44) | match |
| `B_T≈0.134875005736706` (app:848) | `assert abs(BT - 0.1348...)<1e-11` (py:45) | match |

Every paper-side deliverable has a faithful, load-bearing script-side counterpart. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `assert abs(AT - (-4.27263956256927)) < 1e-11` | A_T literal (app:846) | yes |
| A2 | sympy | 45 | `assert abs(BT - 0.134875005736706) < 1e-11` | B_T literal (app:848) | yes |
| A3 | sympy | 109 | `assert exact_resid == 0` (symbolic `(1-λ_{Π,0})-ξ_*`) | claim 3 consistency (notes:124-135) | yes |
| A4 | mathematica | 57-58 | `expectZero["aT vs paper literal A_T", If[Abs[aT-(-4.27...)]<1e-11,0,...]]` | A_T literal (app:846) | yes |
| A5 | mathematica | 59-60 | `expectZero["bT vs paper literal B_T", If[Abs[bT-0.1348...]<1e-11,0,...]]` | B_T literal (app:848) | yes |
| A6 | mathematica | 98 | `expectZero["(1-lambda_(Pi,0)) - xi_*", FullSimplify[(1-lamPiZero)-xiStarClosed]]` | claim 3 consistency | yes |

The printed family/interpolation lines (py:56-88, wl:67-94) are not assertions; they are reconciled value-for-value against the notes in the Value Reconciliation section below. The three assertion pairs are non-tautological: A1/A2/A4/A5 anchor independently-derived `A_T`/`B_T` (different routes per engine) to EXTERNAL appendix literals (a dropped chain term would push `A_T` off `-4.2726...` and fail); A3/A6 assert an exact symbolic collapse of nested radicals that would be nonzero if the convex-family algebra or the `r_F1` constant were wrong.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **INDEPENDENT** derivation of the load-bearing quantity (`A_T`), not a transliteration. The first-pass bug was that the Mathematica `dT` route dropped the S-follows-Π chain term; the fix re-routes `A_T` through Mathematica's own automatic differentiation. Three corresponding sections:

1. **A_T derivation (the divergence point).** SymPy hand-writes the closed form including the chain term:
   - py:32-37 `AT = sp.N(-(9/(40*T_star)) * ( 1/(gp_star*(1-S_star/4)) + Pi_star*Sp_star/(4*gp_star*(1-S_star/4)**2) ), 30)` — the `Pi_star*Sp_star/...` term IS the explicit S-follows-Π chain term, written by the author.
   - wl:47-49 `Tm[pp_] := Sqrt[(9/20)*(pp/(1 - sFormula/4))] /. p -> pp;  aT = N[-(D[Tm[p], p] /. p -> pStar)/gPrimeStar, 30];` — Mathematica substitutes `S = sFormula(p)` then takes `D[Tm[p], p]`, so the chain term is produced by autodiff, never hand-transcribed from the SymPy formula. The comment wl:43-46 documents this is "Mathematica's OWN symbolic differentiation (NOT a port of SymPy AT)" and that it "restores the sPrimeStar chain term the dSigmaOfDeltas route dropped." This is a genuinely different algebraic route to the same number, and the outputs confirm both land on `A_T = -4.2726...` (both pass A1/A4).

2. **B_T derivation.** Both engines use the same explicit fixed-Π S-sensitivity form (py:38-41 `(9/(40*T_star))*Pi_star/(4*(1-S_star/4)**2)`; wl:51 `(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2)`). This single quantity is structurally the same expression in both — acceptable, since it is the trivial closed-form S-derivative and both anchor it to the external appendix literal `0.1348...`; the independence that matters is on `A_T`, which carried the bug.

3. **The ξ_* consistency assertion.** SymPy: py:102-105 builds `one_minus_lam_exact = (π/4 - gminus_exact)/(π/4 - 2/π)` from the EXACT `gminus`, then `exact_resid = sp.simplify(one_minus_lam_exact - xi_star_closed)` and `assert exact_resid == 0`. Mathematica: wl:90 solves `lamPiZero` from `Solve[gLam == gMinus, lam]` directly and wl:98 asserts `FullSimplify[(1 - lamPiZero) - xiStarClosed] === 0`. Different construction (SymPy substitutes the convex-family closed form for `1-λ`; Mathematica solves the linear equation symbolically) reaching the same exact-zero collapse.

Net: PARTIAL-leaning-INDEPENDENT, but the load-bearing/previously-buggy quantity (`A_T`) is now genuinely independently derived (autodiff vs hand closed form), and the other two checks use distinct constructions. I record this as INDEPENDENT for the cross-engine call, since the only shared expression is the trivial `B_T` closed form and the two engines disagree in route everywhere it matters. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit the family shifts to high precision; they AGREE at the precision claimed (~16+ sig figs; remaining digits are precision artifacts of `N[...,30]` vs `sp.N(...,30)`):

| quantity | SymPy (out line) | Mathematica (out line) | agree? |
|---|---|---|---|
| uniform dPi/eps | `1.69941496131429664915...` (7) | `1.6994149613142967238...` (11) | yes (~16 sf) |
| uniform dT/eps | `0.508756302215083911371...` (8) | `0.5087563022150839246579...` (12) | yes (~16 sf) |
| derivative dPi/eps | `-0.382993186095927685585...` (11) | `-0.3829931860959276343730...` (15) | yes (~16 sf) |
| **derivative dT/eps** | `-0.116943802151810766015...` (12) | `-0.1169438021518107487257...` (16) | **yes (~16 sf)** |
| lambda_(Pi,0) | `0.816081594488460293490...` (21) | `0.8160815944884603716916...` (19) | yes (~16 sf) |
| lambda_(T,0) | `0.813099276577333163196...` (22) | `0.8130992765773331742159...` (20) | yes (~16 sf) |
| (1-lam_Pi,0)-xi_* | exact `0` (25) + numeric 7.8e-17 (26) | exact `0` (23) | yes |

Critically, `dTU` and `dTD` (the quantities the first-pass bug corrupted on the Mathematica side) now AGREE cross-engine. The Mathematica `dT/eps` for the uniform family (`0.50875630221508392...`) matches SymPy's (`0.50875630221508391...`); the derivative-family `dT/eps` (`-0.11694380215181074...` both) likewise matches. No silent disagreement remains. No `engine_disagreement` finding.

## Verdict justification

CLEAN. I read the stage card, the full notes file, and the relevant appendix section (Family-1 mouth correction, lines ~540-960) before the scripts. Every paper/notes deliverable maps to a faithful, non-tautological script-side check, and the script constants (`r_F1` radical `4107-100π²`, `A_T`, `B_T`, the convex-family moments) all match the paper. Attacks tried that failed: (1) I checked whether the Mathematica `A_T` still drops the S-follows-Π chain term — it does not; it is now derived via `D[Tm[p],p]` autodiff with `S=sFormula(p)` already substituted, a route independent of SymPy's hand-written closed form, and both land on the published `-4.2726...`. (2) I compared the emitted `dTU`/`dTD` across both committed `.txt` outputs digit-by-digit — they agree to ~16 sig figs (precision-artifact tail only), so the prior silent cross-engine disagreement is resolved. (3) I grepped for the stale `168` class across all six stage-148 files plus the appendix — it is absent; the radical is correctly `100π²` (`stage_appendix_part04.tex:562`), and the lone `168` digit-string in the appendix (line 956, `4.651033550168876`) is an interior mantissa digit, not a `168π²` value. (4) I checked the A1-A6 assertions for tautology — A1/A2/A4/A5 anchor independently-derived quantities to EXTERNAL literals (would fail on a wrong derivation), and A3/A6 assert an exact symbolic radical collapse that would be nonzero under a wrong constant. Outputs are fresh (both `.txt` mtimes 1780117328 > both script mtimes 1780116764/1780116805).

## Self-test notes

Checked: (1) variable-independence — `D[Tm[p], p]` (wl:49) acts on `Tm[pp_]:=Sqrt[(9/20)*(pp/(1-sFormula/4))]/.p->pp`, which genuinely depends on `p` both directly (the `pp`) and through `sFormula(p)`, so the derivative is non-trivial and the chain term is present (the whole point of the fix); the SymPy `sp.diff(gPi,Pi)` and `sp.diff(Sformula,Pi)` (py:27-28) likewise act on Π-dependent expressions. (2) Trivial-case — the two exact-zero assertions (py:109/wl:98) hinge on the nested-radical identity `sqrt(1+r_F1^2)=37√3/(10π)` and `sqrt(4107)=37√3`; the SymPy guard comment (py:96-101) correctly tracks the collapse, and both engines emit exact 0. (3) Anchors — `A_T`/`B_T` literals match appendix part04:846/848 exactly; `100π²` (not `168π²`) confirmed. No directive needed.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 14 values checked, 0 misaligned

All deliverable values the scripts emit/assert reconcile with the notes and/or appendix. Output `.txt` files are fresh (not stale), so the reconciliation is based on script source + committed outputs.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `ḡ_u = 2/π = 0.636619772367581` | py:50 / out:5; wl:62 / out:9 | notes:35 `0.636619772367581` | MATCH |
| `S̄_u = 0.583877311158896` | py:51 / out:6; wl:63 / out:10 | notes:37 `0.583877311158896` | MATCH |
| `δΠ_u/ε = 1.699414961314297` | py:53 / out:7; wl:64 / out:11 | notes:44 `+1.699414961314297` | MATCH |
| `δT̂_{m,u}/ε = 0.508756302215084` | py:54 / out:8; wl:65 / out:12 | notes:45 `+0.508756302215084` | MATCH |
| `ḡ_d = π/4 = 0.785398163397448` | py:62 / out:9; wl:72 / out:13 | notes:63 `ḡ_d=π/4` (symbolic) | MATCH |
| `S̄_d = 0.657844575502831` | py:63 / out:10; wl:73 / out:13 | notes:70 `0.657844575502831` | MATCH |
| `δΠ_d/ε = -0.382993186095928` | py:65 / out:11; wl:74 / out:15 | notes:76 `-0.382993186095928` | MATCH |
| `δT̂_{m,d}/ε = -0.116943802151811` | py:66 / out:12; wl:75 / out:16 | notes:78 `-0.116943802151811` | MATCH |
| `δΠ_λ slope = -2.082408147410224` | py:80 / out:14; wl:87 / out:17 | notes:103 `2.082408147410224 λ` | MATCH |
| `δT̂_{m,λ} slope = -0.625700104366895` | py:82 / out:18; wl:88 / out:18 | notes:112 `0.625700104366895 λ` | MATCH |
| `λ_{Π,0} = 0.816081594488460` | py:84 / out:21; wl:90 / out:19 | notes:119 `0.816081594488460` | MATCH |
| `λ_{T,0} = 0.813099276577333` | py:85 / out:22; wl:91 / out:20 | notes:141 `0.813099276577333` | MATCH |
| `1-λ_{Π,0} = 0.183918405511540` | py:88 / out:23; wl:94 / out:21 | notes:124 `0.183918405511540` | MATCH |
| `A_T = -4.27263956256927` | py:44 (asserted); wl:57-58 / out:5 | appendix:846 `-4.27263956256927` | MATCH |
| `B_T = 0.134875005736706` | py:45 (asserted); wl:59-60 / out:7 | appendix:848 `0.134875005736706` | MATCH |
| `ξ_* (Stage 126 closed form) = 0.183918405511540` | py:94 / out:24; wl:97 / out:22 | notes:131-134 closed form `≈0.183918405511538` | MATCH (exact resid 0 asserted) |

INTERNAL (scaffolding, not expected in prose; no finding): `Pi_star`, `g_star`, `S_star`, `gp_star`/`g'_*` (does appear at appendix:831 `≈0.0714453558083195`, but is an intermediate, not a script-printed deliverable), `Sp_star`, `Sigma_star`, `T_star`/`tStar`, and the per-engine exact/numeric consistency residuals (py:105-107, out:25-26; wl:98).

Note: the notes give `ξ_*` numeric `0.183918405511538` (notes:134) vs the script's `1-λ_{Π,0}` `0.183918405511540` — they differ at the ~15th digit, but this is the expected numerical noise the notes themselves call out ("agrees to numerical precision", notes:126); both engines independently assert the EXACT symbolic residual is 0 (py:109/wl:98, out:25/out:23). MATCH, not a value_mismatch.
