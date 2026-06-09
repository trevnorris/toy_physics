---
unit_id: 175
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage175_wall_normalized_load_shape.md]
  paper_appendix: present
---

# Audit unit 175 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_175.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage175_wall_normalized_load_shape.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 81 status row; lines 610-626 outgoing-load theorem subsection; line 1464 result-anchor catalog; line 1527 stage input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.txt`

## What the paper claims

Stage 175 (`\stagefield{Output}`: "Factors the remaining defect into wall-normalized shape/load variables and proves conservative-shape preservation leaves only outgoing load slippage.") factors the Stage 174 microscopic static bundles by the wall baseline `K`. The notes give the load-bearing deliverables as boxed identities:
(1) homogeneity factorization `B0 = K χ²`, `Z0 = K Υ`, `N0 = Λ²` (notes §2.1-2.2, with intermediate scaling laws `Δ = K²Δ̂`, `Q = K³Q̂`, `P = K²P̂`);
(2) the differential defect rewrites `Σ_B = δln χ²`, `Σ_Z = δln Υ`, `Σ_N = 2 δln Λ − δ_K = δln(Λ²/K)` (notes §2.3, with `δ_K := K₁/K`);
(3) the conservative-shape reduction: with `δln χ² = δln Υ = 0`, `Ξ_load = Σ_r ρ_r^(N)(2 δln Λ_r − δ_K)` (notes §4);
(4) the naive common-self-similarity no-go: freezing ALL shapes gives `Ξ_load = −δ_K` because `Σ_r ρ_r^(N) = 1` (notes §5).
There are NO numeric figures-of-merit; every deliverable is a symbolic identity. The appendix (lines 610-626) re-states `Λ_r := P_r/Δ_r` as the controlling factor and notes Stage 176 (not 175) further factors `Λ²/K` into mixed-leg/interference/hybridization ratios.

## What the script claims to verify

The SymPy docstring enumerates four checks matching the notes exactly: (1) wall-normalized factorization `B0=Kχ²`, `Z0=KΥ`, `N0=Λ²`; (2) defect rewrites `Σ_B=δln χ²`, `Σ_Z=δln Υ`, `Σ_N=δln(Λ²/K)`; (3) conservative-shape reduction `Ξ_load = ⟨2 δln Λ − δK⟩_N`; (4) common-shape branch `Ξ_load = −δK`. Each `expect_zero`/`expectZero` asserts the residual of one boxed identity is exactly 0. The homogeneity laws are exercised by building `Δ, Q, P` from physical primitives, applying the wall scaling `subs_hat` (`ou→K·ouhat`, …), and confirming the K-powers. The differential identities introduce an `eps`-flow (`K = K0·exp(κ·eps)`, etc.), take the eps-slope of the log, and confirm `2 δln(P/Δ) − κ = δln(Λ²/K)` where `κ` plays `δ_K`.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| `B0 = K χ²` | `B0 - K*chi^2` (py 50 / wl 40) | match |
| `Δ = K²Δ̂`, `Q = K³Q̂`, `P = K²P̂` | three K-power checks (py 67-78 / wl 53-58) | match |
| `Z0 = K Υ` | `Z0 - K*Upsilon` (py 82 / wl 62) | match |
| `N0 = Λ²` | `N0 - Lambda^2` (py 83 / wl 63) | match |
| `Σ_B = δln χ²` | `Sigma_B - dln(chi^2)` (py 118 / wl 85) | match |
| `Σ_Z = δln Υ` | `Sigma_Z - dln(Upsilon)` (py 125 / wl 91) | match |
| `Σ_N = 2 δlnΛ − δ_K = δln(Λ²/K)` | `Sigma_N - dln(Lambda^2/K)` (py 135 / wl 100) + series route (wl 108) | match |
| conservative-shape: `Σ_B=Σ_Z=0` on frozen shapes | `Conservative-shape branch Sigma_B/Sigma_Z` (py 154-155 / wl 121-122) | match |
| no-go: `Σ_N = −δ_K` frozen, `Ξ_load = −δ_K` | `Common-shape branch Sigma_N + dK` + `Xi_load (all shapes frozen) + dK` (py 156/162 / wl 123/128) | match |

`paper_alignment: aligned`. The appendix's `Λ²/K = M²(1+I)²/(1−H)²` (line 622) is explicitly a Stage 176 result; absence of a 175-side check is correct.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `simplify(B0 - K*chi^2)==0` | claim 1 (B0) | yes |
| A2 | sympy | 67-78 | three `expand(...)−K^n(...)==0` | claim 1 (Δ/Q/P scaling) | yes |
| A3 | sympy | 82-83 | `Z0−K*Υ==0`, `N0−Λ²==0` | claim 1 (Z0,N0) | yes |
| A4 | sympy | 118 | `Σ_B_direct − dln(chi²)==0` | claim 2 (Σ_B) | yes |
| A5 | sympy | 125 | `Σ_Z_direct − dln(Υ)==0` | claim 2 (Σ_Z) | yes |
| A6 | sympy | 135 | `2 dln(P/Δ)−κ − dln(Λ²/K)==0` | claim 2 (Σ_N) | yes |
| A7 | sympy | 154-156 | frozen-branch `Σ_B=Σ_Z=0`, `Σ_N+κ=0` | claims 3,4 | yes |
| A8 | sympy | 162 | `(ρ1+ρ2)Σ_N|_{sum=1}+κ==0` | claim 4 (Ξ_load=−δ_K) | yes |
| B1-B8 | math | 40,53-63,85,91,100,121-128 | `expectZero[...]` mirrors A1-A8 | same claims | yes |
| B9 | math | 108 | `2 dlogSeries(P/Δ)−κ − dlogSeries(Λ²/K)==0` | claim 2 (Σ_N) via Series route | yes |

All rows non-tautological: each builds LHS and RHS by independent symbolic routes (e.g. A6 builds `2 dln(P/Δ)` from physical primitives via `subs_hat` BEFORE the cached `Lambda`, so the `−κ` subtraction is load-bearing — a wrong `κ` fails the check). The frozen-branch checks substitute concrete shape slopes (`schi=0`, `su=sw=sr=sgu=sgw=0`) and confirm the residual collapses to the predicted `0` or `−κ`.

## Findings

None.

## Independent-derivation check (Mathematica)

This is the unit at the END of the 105-175 transliteration watch and was a first-pass batch-8 finding. I re-audited fresh.

**Assessment: the `.wl` is PREDOMINANTLY A TRANSLITERATION of the `.py`, with one genuinely independent route bolted on for the `Sigma_N` slope.** Confidence: high.

Evidence of porting (the bulk of the script):
- Identical variable choreography and naming: `b0`, `delta`, `q`, `p`, `subsHat`, `subsEps`, `upsilon`, `lambda`, `sigmaBDirect`, `sigmaZDirect`, `sigmaNDirect`, `xiLoadFrozen` — one-to-one with the `.py` names.
- Identical 14-check sequence in identical order with byte-identical check labels (`"B0 - K*chi^2"`, `"Delta - K^2*Delta_hat"`, … `"Xi_load (all shapes frozen) + dK"`). The `.wl` adds exactly one extra (`"Sigma_N - dln(Lambda^2/K) [series route]"`), total 15 vs 14.
- The primary `dlog` helper is a line-for-line port of the SymPy one:
  - `.py:100-101` — `def dlog(expr): return sp.simplify(sp.diff(sp.log(sp.simplify(expr)), eps).subs(eps, 0))`
  - `.wl:26-29` — `dlog[expr_] := FullSimplify[D[Log[FullSimplify[expr, ...]], eps] /. eps -> 0, ...]`
  Same physical premise, same algebra (differentiate `Log`, evaluate at `eps→0`), same outer `simplify`/`FullSimplify` wrapper.
- The Sigma_N comment block (`.wl:109-115`) is a verbatim transliteration of the `.py:136-142` "red-team F1 resolution" comment.

The ONE independent element (the batch-8 remediation):
- `.wl:31` — `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps];` — extracts the first-order log-slope via Taylor-series coefficient (`Series`+`Coefficient`), a mechanically distinct procedure from symbolic `D[Log[...]]`.
- `.wl:106-108` — re-derives `Sigma_N` by the Series route and confirms `2 dlogSeries(P/Δ) − κ = dlogSeries(Λ²/K)` (B9). Saved output line 27-28 shows it PASSing at residual 0.

So the remediation is genuinely present and the `Sigma_N` slope identity is now confirmed by TWO independent slope mechanisms inside Mathematica. However, this does NOT make the whole `.wl` an independent re-derivation — `B0/Z0/N0` homogeneity, `Σ_B`, `Σ_Z`, and the conservative/common-shape reductions still ride the ported `dlog` route only. Per the dual-engine policy this is acceptable here because (a) the homogeneity checks (A1-A3/B1-B3) are pure polynomial/`Expand` identities with no engine-specific algebra to "port" — they are the same identity any CAS must produce, and (b) the slope identity that was the actual first-pass concern now has an engine-native second route. I record the residual port nature as an observation, not a finding: re-flagging `mathematica_transliteration` would contradict the user-accepted batch-8 resolution that ADDED the independent route rather than rewriting every line. The independent route is real and load-bearing; the structural porting elsewhere is on identities where independence is not meaningfully definable.

## Engine cross-check

Both engines present; outputs agree exactly. Every shared check (`B0 - K*chi^2`, …, `Xi_load (all shapes frozen) + dK`) emits residual `0` in both transcripts (sympy txt lines 5-25; mathematica txt lines 5-40, each followed by `PASS`). The Mathematica-only series route (`Sigma_N - dln(Lambda^2/K) [series route] = 0`, mathematica txt 27-28) PASSes, corroborating the symbolic-`D` route. No sign, factor, or residual disagreement.

## Verdict justification

`clean`. Attacks tried and failed: (1) tautology hunt on `Sigma_N` — the `2 dln(P/Δ) − κ` route is built from physical primitives via `subs_hat` BEFORE the cached `Lambda`, so the `−κ` term is load-bearing (a wrong `κ` breaks the identity), not a `P/Δ ≡ Lambda` simplify-commute; the prior simplify-commute lines were correctly removed per the F1-resolution comments. (2) Symbol-domain check — all symbols are `positive, real` matching the physical setup; the `eps`-flow uses `exp(slope·eps)` so logs are well-defined; no aggressive assumption hides a branch. (3) Naming-convention check — both engines use `Δ = ou·ow − r²` with `ou ≡ Ω_U²` (not `Ω_U`), consistent with `subs_hat: ou→K·ouhat` giving `Δ = K²Δ̂`, matching the notes' `Δ = K²Δ̂`; convention is shared by both engines and self-consistent. (4) No-go check — `Ξ_load = (ρ1+ρ2)·Σ_N` with `ρ1+ρ2→1` correctly encodes `Σ_r ρ_r^(N) = 1` and yields `−κ`. (5) Stale-value watch — searched all four 175 files for `168`/`100π²`/`4107`/`10π`; none present (they belong to a different stage; not a 175 concern). I read the paper card, notes, and appendix first; the script's verified claims match the paper's stated deliverables one-to-one. Outputs are fresh (both `.txt` mtimes 2026-05-29 23:50, newer than `.py` 2026-05-28 16:16 and `.wl` 2026-05-29 23:44).

## Value Reconciliation (pass-2 augmentation)

This stage emits NO numeric constants or figures-of-merit. Every deliverable is a symbolic closed-form identity, asserted to residual 0. The reconciliation maps each script-emitted symbolic result to its boxed home in the notes (the `.tex` card is terse by design and carries the Output sentence + `Λ_r` reference; the notes are the natural carrier).

| value (symbolic deliverable) | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `B0 = K χ²` | py:50, wl:40; sympy txt 5 / math txt 5-6 | notes:117 (boxed) | MATCH |
| `Δ = K²Δ̂` | py:67-70, wl:53; sympy txt 6 / math txt 7-8 | notes:173 (boxed) | MATCH |
| `Q = K³Q̂` | py:71-74, wl:54-57; sympy txt 7 / math txt 9-10 | notes:174 (boxed) | MATCH |
| `P = K²P̂` | py:75-78, wl:58; sympy txt 8 / math txt 11-12 | notes:175 (boxed) | MATCH |
| `Z0 = K Υ`, `Υ = Q/(KΔ)` | py:80,82, wl:60,62; sympy txt 9 / math txt 13-14 | notes:184,186 (boxed) | MATCH |
| `N0 = Λ²`, `Λ = P/Δ` | py:81,83, wl:61,63; sympy txt 10 / math txt 15-16 | notes:191,193 (boxed); .tex:`eq:app-part05-load-factor-Lambda` appendix:614 | MATCH |
| `Σ_B = δln χ²` | py:118, wl:85; sympy txt 15 / math txt 21-22 | notes:126 (boxed) | MATCH |
| `Σ_Z = δln Υ` | py:125, wl:91; sympy txt 16 / math txt 23-24 | notes:218 (boxed) | MATCH |
| `Σ_N = δln(Λ²/K) = 2 δlnΛ − δ_K` | py:135, wl:100,108; sympy txt 17 / math txt 25-28 | notes:223-224 (boxed) | MATCH |
| conservative branch: `Σ_B=Σ_Z=0` | py:154-155, wl:121-122; sympy txt 22-23 / math txt 33-36 | notes:288,290 (Θ_B=Θ_Z=0) | MATCH |
| `Ξ_load = −δ_K` (all shapes frozen) | py:156,162, wl:123,128; sympy txt 24-25 / math txt 37-40 | notes:333 (boxed); appendix:608 (`Ξ_load`) | MATCH |

INTERNAL scaffolding (accounted for, no finding): `expect_zero`/`expectZero`/`pass`/`fail`/`banner`/`fmt` helpers; `dlog`/`dlogSeries` helpers; the `eps`-flow drivers (`K=K0 exp(κ eps)`, etc.); intermediate `Sigma_*_direct`/`Sigma_*_shape` expressions; `subs_hat`/`subs_eps` substitution maps; the conclusions `print` block.

Stale-value watch (per task): searched for `168`, `100π²`, `4107`, `10π`, Family-1 radius `√(4107−100π²)/(10π)` in all four 175 files — NONE present. Those values belong to a different stage; no 175-side mismatch.

reconciliation: complete; 11 deliverable values checked, 0 misaligned

## Self-test notes

Checked: (1) variable independence — every `diff(log(EXPR), eps)`/`D[Log[EXPR], eps]` and `dlogSeries` operates on expressions genuinely carrying `eps` (via `subs_eps`/`subsEps` exponentials), so no slope is identically zero by construction; the frozen-branch substitutions deliberately zero specific slopes and the residuals correctly collapse to `0` or `−κ`. (2) Symmetry/parity — n/a, no integrals over unbounded domains. (3) Trivial-case — substituting the frozen profile (`schi=0`, all `s*=0`) into `Σ_N_direct` gives `−κ` and into `Σ_B/Σ_Z` gives `0`, matching the asserted residuals. No directive written (zero findings).
