---
unit_id: 013
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 013 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_013.tex`
- notes: `(none)` — no `notes/stages/moving_throat_pde_stage013_*.md` exists (confirmed by directory listing; prompt declared `(none)`).
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 48 + `\input{stages/stage_013}` at line 105)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.txt`

## What the paper claims

Stage 013 is the "Mouth-Taylor master map." The card's `\stagefield{Output}` states verbatim: "Stage~013 exports the mouth-Taylor primitive map \eqref{eq:stage007-z0n0}--\eqref{eq:stage007-z4}." The deliverables are the one-sided-Taylor primitive slippage coefficients evaluated along the mouth coordinate `s`: the leading slots `z_0=∂_s(Q/Δ)`, `n_0=∂_s(P²/Δ²)` (eq:stage007-z0n0); the next conservative slot `z_2=∂_s((Q S₂ − H_port Δ)/Δ²)` (eq:stage007-z2); and `z_4=∂_s((Q(S₂²−Δ) − H_port S₂ Δ)/Δ³)` (eq:stage007-z4). The gate-feed paragraph adds the first-order prefactor-slope relation `Ξ_load = n_0/N_0 + z_0/D_0` (eq:stage013-xi), and states the audit "checks the one-sided Taylor derivatives and their dependence on the primitive bottlenecks." `N_0` and `D_0` are carried as abstract forward symbols in the card (not given closed forms here). The appendix row marks the stage `\StatusExactClosure{}`: "Mouth-local Taylor map for projected source and flux corrections." There is no stage-013 notes file; the .tex card is authoritative.

## What the script claims to verify

The SymPy script verifies three things. (1) The one-sided exponential-weight Taylor projection reproduces the first moment: `∫₀^∞ e^{-u}X du` expanded to order `ell` equals `X0 + ell X1` (line 26), with a deliberately mutated control that must fail (line 29). (2) The `W2 = u e^{-u}` weight has first moment `μ₁=∫u² e^{-u}=2` (line 34) and reproduces `X0 + ell·μ₁·X1` (line 35), again with a mutated control (line 38). (3) The paper round-trip: it builds `Xi` from the prefactor-slope combination `2p1/P − 2d1/Δ + q1/(D0 Δ) − Q d1/(D0 Δ²)` and independently builds `Xi_paper = n0_form/N0_form + z0_form/D0` with `N0 = P²/Δ²`, asserting they are equal (line 107), plus `dXi/dPprime = 2/P` (line 108) with a nonzero-dependence control (line 110). The closed-form literals `z0,z2,z4,n0,n2,n4` (lines 58–98) are assigned in SymPy but only `z0_form`/`n0_form` are exercised there; the comments (lines 49–57, 71–79) explicitly delegate independent verification of all six literals to the Mathematica script (claims M3/M4). The Mathematica script independently derives `z0..z4` and `n0..n4` from generating sources `Zsource[x]=(q(t)−h(t)x²)/(d(t)−s(t)x²+x⁴)` and `Nsource[x]=(p(t)−g(t)x²)²/(d(t)−s(t)x²+x⁴)²` via `Series` to order 4, time-derivative at `t=0`, and asserts each equals the same literal (M3/M4), plus mirrors M1/M2/M5.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `z_0=∂_s(Q/Δ)` = `(Δ q1 − Q d1)/Δ²` | sympy `z0` line 58 / `z0_form` line 103; mathematica `z0Expected` line 72 asserted M3 (line 81) | match |
| `n_0=∂_s(P²/Δ²)` = `2P(Δ p1 − P d1)/Δ³` | sympy `n0` line 80 / `n0_form` line 104; mathematica `n0Expected` line 91 asserted M4 (line 104) | match |
| `z_2=∂_s((Q S₂ − H Δ)/Δ²)` | sympy `z2` line 59; mathematica `z2Expected` line 74 asserted M3 (line 82) | match (hand-expansion confirms) |
| `z_4=∂_s((Q(S₂²−Δ) − H S₂ Δ)/Δ³)` | sympy `z4` lines 62–69; mathematica `z4Expected` line 77 asserted M3 (line 83) | match (full hand-expansion confirms, all 8 terms) |
| Gate feed `Ξ_load = n_0/N_0 + z_0/D_0` | sympy `Xi == Xi_paper` line 107; mathematica M5 line 128 | match (with `N_0 = P²/Δ²` natural identification, paper leaves `N_0` abstract) |
| "checks one-sided Taylor derivatives" (projection mechanics) | sympy M1 line 26, W2 lines 34–35; mathematica M1 line 41, M2 lines 47–53 | match |
| `n_2,n_4` (flux-correction higher slots, implied by appendix "flux corrections") | sympy literals lines 81–98; mathematica `n2Expected`/`n4Expected` asserted M4 (lines 105–106) | match (cross-verified by generating source) |

No paper deliverable is unmatched; no script assertion is orphaned from a paper deliverable. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 26 | `series(Xproj,ell)−(X0+ell X1)==0` | "one-sided Taylor derivatives" mechanics | yes |
| A2 | sympy | 29 | `assert_nonzero` mutated first moment | control for A1 | yes |
| A3 | sympy | 34 | `μ₁(W2)−2==0` | Taylor projection (W2 weight) | yes |
| A4 | sympy | 35 | `series(Xproj_W2)−(X0+ell μ₁ X1)==0` | Taylor projection | yes |
| A5 | sympy | 38 | `assert_nonzero` mutated | control for A4 | yes |
| A6 | sympy | 107 | `Xi − Xi_paper == 0` | gate feed `Ξ_load=n0/N0+z0/D0` (eq:stage013-xi), embeds z0/n0 | yes |
| A7 | sympy | 108 | `diff(Xi,Px) − 2/P == 0` | "dependence on primitive bottlenecks" (P′) | yes |
| A8 | sympy | 110 | `assert_nonzero diff(Xi,Px)` | control for A7 | yes |
| M1 | mathematica | 41 | `Series(Xproj)−(X0+ell X1)==0` | Taylor projection | yes |
| M2a | mathematica | 47 | `∫u²e^{-u}−2==0` | Taylor projection (W2) | yes |
| M2b | mathematica | 51 | `Series(Xproj2)−(X0+2 ell X1)==0` | Taylor projection (W2) | yes |
| M3 | mathematica | 81–83 | `z0,z2,z4 − Expected == 0` | z_0, z_2, z_4 deliverables | yes |
| M4 | mathematica | 104–106 | `n0,n2,n4 − Expected == 0` | n_0 (+ n_2,n_4) deliverables | yes |
| M5a | mathematica | 128 | `Xi − XiPaper == 0` | gate feed Ξ_load | yes |
| M5b | mathematica | 129 | `D[Xi,Px] − 2/P == 0` | P′ dependence | yes |

Every row is non-tautological and traces to a stated paper deliverable.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent derivation, not a transliteration**. For the load-bearing z/n coefficients the two engines use structurally different methods:
- SymPy hardcodes the closed-form literals `z0..z4`, `n0..n4` (lines 58–98) and exercises only `z0_form`/`n0_form` in the Xi round-trip. It does NOT expand the generating ratio at all.
- Mathematica builds the generating sources `Zsource[x]=(qFun[t]−hFun[t]x²)/(dFun[t]−sFun[t]x²+x⁴)` (line 66) and `Nsource[x]=(pFun[t]−gFun[t]x²)²/(...)²` (line 85), `Series`-expands in `ell` to order 4 (lines 67, 86), extracts the ell⁰/ell²/ell⁴ coefficients, takes the time-derivative at `t=0` (`timeDerivativeAtZero`, line 63), and asserts the result equals each literal (M3 lines 81–83, M4 lines 104–106). This is exactly the derivation the SymPy comments (lines 51–57, 73–79) delegate to it — Mathematica is the derivation engine for the coefficients, SymPy is the literal-bearer.

The overlap on M1/M2/M5 (both engines compute the exponential-weight Taylor moments and the Xi round-trip by direct integration/simplification) is benign and expected; the load-bearing coefficient derivation is genuinely independent. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass and agree. SymPy output: `STATUS: PASS` (3 print lines, all in-file asserts passed by not raising). Mathematica output: 11 `OK … residual = 0` lines (M1, M2 second moment, M2 first moment, M3 z0/z2/z4, M4 n0/n2/n4, M5 Xi round-trip, M5 dXi/dPprime) then `STAGE 013 MATHEMATICA AUDIT: PASS`. The Mathematica run additionally certifies (residual=0) the four higher coefficients `z2,z4,n2,n4` that SymPy carries as unexercised literals, so the two engines jointly cover all six coefficient deliverables with zero residual. No `engine_disagreement`.

## Verdict justification

`clean`. I read the paper card, confirmed there is no stage-013 notes file, and read the appendix row before opening the scripts. Attacks tried and failed: (1) hand-expanded the paper's explicit closed forms for `z_0`, `z_2`, and `z_4` (the latter all 8 numerator terms) and confirmed each matches the script literal exactly — no factor/sign error. (2) Checked the Xi round-trip is non-tautological: `Xi` and `Xi_paper` are built from two distinct algebraic expressions, and `n0_form/N0_form` genuinely reduces to `2p1/P − 2d1/Δ` while `z0_form/D0` reduces to `q1/(D0Δ) − Qd1/(D0Δ²)`, so A6/M5a are real identities, not constructions equal by definition. (3) Checked `dXi/dPx=2/P`: `Xi` does depend on `Px` (avoids the zero-derivative trap), and the W2-control / mutated-moment `assert_nonzero` guards (A2/A5/A8) make the positive checks meaningful. (4) Verified the `μ₁(W2)=2` literal is the correct `Γ(3)=2`. The only convention assumption — identifying `N_0=P²/Δ²` in the Xi round-trip — is explicitly flagged in the script comment as "the natural identification," and the paper card leaves `N_0` as an abstract forward symbol, so it is not a misalignment. Symbol domains (`positive` for `u,ell`; `nonzero` for `Δ,P,D0,μ1`; `Reals` in Mathematica) are consistent with the physical setup and across engines. Outputs are fresh (txt mtimes 2026-05-25 21:54 newer than script mtimes 21:52).

## Value Reconciliation (pass-2 augmentation)

Every deliverable value the scripts emit is a **symbolic** closed form (this stage produces no named numeric figures of merit). I enumerate each boxed/asserted symbolic deliverable and locate it in the `.tex` card (no `.md` notes exist for this stage).

| value | source (py/wl + output line) | .tex location | status |
|---|---|---|---|
| `z_0 = (Δ q1 − Q d1)/Δ²` = `∂_s(Q/Δ)` | py:58, wl:72 (M3 OK, mathematica out line 5) | stage_013.tex:16 (`z_0=∂_s(Q/Δ)`) | MATCH |
| `n_0 = 2P(Δ p1 − P d1)/Δ³` = `∂_s(P²/Δ²)` | py:80, wl:91 (M4 OK, mathematica out line 8) | stage_013.tex:18 (`n_0=∂_s(P²/Δ²)`) | MATCH |
| `z_2 = (−Δ²h1 + Δ(H d1+Q s1+S₂ q1) − 2Q S₂ d1)/Δ³` = `∂_s((Q S₂ − H Δ)/Δ²)` | py:59, wl:74 (M3 OK, mathematica out line 6) | stage_013.tex:23–25 | MATCH (hand-expansion) |
| `z_4` (8-term Δ⁴ numerator) = `∂_s((Q(S₂²−Δ) − H S₂ Δ)/Δ³)` | py:62–69, wl:77 (M3 OK, mathematica out line 7) | stage_013.tex:27–31 | MATCH (full hand-expansion) |
| `Ξ_load = n_0/N_0 + z_0/D_0` (verified `Xi==Xi_paper`) | py:107, wl:128 (M5 OK, mathematica out line 11) | stage_013.tex:35–38 (eq:stage013-xi) | MATCH (paper leaves N_0,D_0 abstract; script uses natural N_0=P²/Δ²) |
| `μ₁(W2)=2` (Taylor-projection moment) | py:34, wl:47 (mathematica out line 3) | stage_013.tex:39–41 ("checks the one-sided Taylor derivatives") | MATCH (mechanics, prose-level; not a boxed deliverable) |

INTERNAL items (genuine scaffolding, no prose expected, no finding):
`Xproj`, `Xproj_W2`, `Xproj2`, `W`, `W2`, `X` (projection intermediates); `n_2`, `n_4` literals and `z2Expected`/`z4Expected`/`n2Expected`/`n4Expected` (higher flux/source slots, only generically referenced by the appendix "source and flux corrections" wording — not individually boxed in the card; cross-verified zero-residual by M3/M4 so not MISSING); the mutated-control residuals; `subs_der`/`subsDer`; `dXi/dPx=2/P` (a dependence check, not a reported value); all PASS/OK flags and `residual=0` strings.

reconciliation: complete; 6 deliverable values checked, 0 misaligned. Every emitted closed-form deliverable reconciles with the stage-013 card. The card legitimately states `z_0,n_0,z_2,z_4` as its boxed map and leaves the higher `n_2,n_4` and the `N_0,D_0` symbols abstract; nothing the scripts box as a stage deliverable is absent from the card.

## Self-test notes

Checked: (1) Variable independence — `diff(Xi,Px)` is non-trivial because `Xi` genuinely contains `Px` (via `p1→μ1 Px`); avoids the identically-zero-derivative trap. (2) Symmetry/parity — the integrals are over `[0,∞)` with `e^{-u}` / `u e^{-u}` weights (Gamma integrals), all convergent and correctly evaluated (`∫u²e^{-u}=2`); no symmetric-domain cancellation claim is made. (3) Trivial-case — the z₀/z₂/z₄ literals were each reconstructed by independent hand chain-rule from the paper's Δ-ratio closed forms and matched term-by-term, and the Xi round-trip reduces algebraically as claimed. No trap fired; verdict stands as clean.
