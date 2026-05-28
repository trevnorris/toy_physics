---
unit_id: 170
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage170_linear_grouped_outlet_map.md]
  paper_appendix: present
---

# Audit unit 170 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_170.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 71, 416-459, 1463 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.txt`

## What the paper claims

Stage 170 solves the linear grouped-`P2` direct outlet map. Card `\stagefield{Output}`: "Reduces the direct outlet deformations to \(\mathcal K_A=\delta D_{A,2}+\delta D_{A,0}/9\) and \(\mathcal G_A=\delta N_{A,0}-P_0\delta D_{A,0}\)." The notes enumerate the full deliverable set: (1) linear grouped transport laws `δu2 = -(δD2 + δD0/9)/D0`, `δu4 = -(δD4 + 2δD2/9 + 5δD0/81)/D0`, `δP0 = (δN0 - P0 δD0)/D0` (Sec. 1); (2) the exact outlet map `δκ_W = 3(1-σ*)/(σ* D0)·(δD2+δD0/9)` and `δγ_W = -(1-σ*)/(9σ* N0)·(δN0 - P0 δD0)` (Sec. 2); (3) the even one-parameter consistency relation `δD4 = (2/3)δD2 + (1/27)δD0` together with `δu4 = (8/9)δu2` (Sec. 3); (4) the grouped trace/anomaly transport for `a_κ,b_κ,a_γ,b_γ` and `a_{D,4},b_{D,4}` (Sec. 4); and (5) the weak-axisymmetric branch result that the maps inherit the grouped signature `(λ_20,λ_21,λ_22)=(1,1/2,-1)`, collapsing the problem to two scalar amplitudes `κ1`, `γ1` with explicit forms (Sec. 5). The card's `\stagefield{Checks}` explicitly lists "Check the weak-axisymmetric signature \((1,1/2,-1)\) before reducing grouped defects to a scalar."

## What the script claims to verify

The SymPy docstring lists four checks: (1) the linear transport formulas for δu2, δu4, δP0; (2) the exact map into δκ_W and δγ_W; (3) the even-consistency relation `δD4 = (2/3)δD2 + (1/27)δD0`; (4) the grouped trace/anomaly transport formulas. The assertions linearize the rational definitions `u2=-D2/D0`, `u4=(D2²-D0 D4)/D0²`, `P0=N0/D0` about the canonical branch (`u2=1/9`, `u4=4/81`, so `D2=-D0/9`, `D4=-D0/27`) via a first-order series in `eps`, invert the Stage-244 hybrid relations by `solve`, and check each reduced form against the paper's boxed formula by `expect_zero`/`expectZero`. The Mathematica script mirrors this exactly. Neither script touches the Section-5 weak-axisymmetric signature/amplitude deliverable.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Sec.1 δu2 transport | py L57 / wl L60 | match |
| Sec.1 δu4 transport | py L58 / wl L61 | match |
| Sec.1 δP0 transport | py L59 / wl L62 | match |
| Sec.2 δκ_W map | py L73-76 / wl L77-80 | match |
| Sec.2 δγ_W map | py L77-80 / wl L81-84 | match |
| Sec.3 δu4 = (8/9)δu2 | py L89 / wl L90 | match |
| Sec.3 δD4 = (2/3)δD2+(1/27)δD0 | py L90-91 / wl L91-95 | match |
| Sec.4 trace/anomaly a_κ,b_κ,a_γ,b_γ | py L120-128 / wl L118-126 | match |
| Sec.5 weak-axisym signature (1,1/2,-1), κ1, γ1 | (none) | missing |

Section 5 is the only paper-side deliverable with no script-side counterpart; the card's `\stagefield{Checks}` calls it out explicitly. Sets `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 57 | `expect_zero(du2 + (dD2+dD0/9)/D0)` | Sec.1 δu2 | yes |
| A2 | sympy | 58 | `expect_zero(du4 + (dD4+2dD2/9+5dD0/81)/D0)` | Sec.1 δu4 | yes |
| A3 | sympy | 59 | `expect_zero(dP0 - (dN0-P0 dD0)/D0)` | Sec.1 δP0 | yes |
| A4 | sympy | 73-76 | `expect_zero(dkappa_from_du2 - 3(1-σ)(dD2+dD0/9)/(σ D0))` | Sec.2 δκ_W | yes |
| A5 | sympy | 77-80 | `expect_zero(dgamma_from_dP0.subs(...) + (1-σ)(...)/(9σ N0))` | Sec.2 δγ_W | yes |
| A6 | sympy | 89 | `expect_zero(du4_from_kappa - (8/9)du2)` | Sec.3 δE4/δE2=8/9 | yes (du2 cancels; checks the two prefactors 1/3, 8/27) |
| A7 | sympy | 90-91 | `expect_zero(relation - (2/3 dD2 + dD0/27))` | Sec.3 δD4 relation | yes |
| A8 | sympy | 120-121 | `expect_zero(a/b_kappa_from_map - a/b_kappa)` | Sec.4 κ trace/anomaly | yes (derived map vs posited formula) |
| A9 | sympy | 122-128 | `expect_zero(a/b_gamma_from_map - a/b_gamma)` under P0→N0/D0 | Sec.4 γ trace/anomaly | yes |
| A10 | mathematica | 60-127 | `expectZero[...]` mirror of A1-A9 | same claims | yes (but transliterated — see F2) |
| — | — | — | (no check for Sec.5 signature/κ1/γ1) | Sec.5 | missing |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** script_missing_paper_claim
**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_170.tex:21` (Checks) and notes Sec.5
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md:370-432`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py` (no corresponding block)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl` (no corresponding block)

**What's wrong:**
The notes (Sec. 5) and the stage card both enumerate the weak-axisymmetric branch as a deliverable, but neither script verifies it. Notes lines 385-424 box four results: the inherited signature
`δκ_W^(20)=ε κ1, δκ_W^(21)=(ε/2)κ1, δκ_W^(22)=-ε κ1`, `κ1 = 3(1-σ*)/(σ* D0)·(D2^(1)+D0^(1)/9)`, and the matching `δγ_W^(A)`/`γ1` set. The card explicitly requires it:
`stage_170.tex:21` — "Check the weak-axisymmetric signature \((1,1/2,-1)\) before reducing grouped defects to a scalar."
The scripts contain no occurrence of `λ`, `(1,1/2,-1)`, `κ1`, `γ1`, or `eps·λ_A` substitution (confirmed by grep — only the `eps` series variable appears, used solely for the Sec. 1 linearization). The grouped trace/anomaly checks (Sec. 4, `a_*`/`b_*`) cover the real/imaginary decomposition but NOT the 20/21/22 lane signature.

**Why this matters:**
The card lists a check the verification scripts do not perform. Section 5 is a near-trivial corollary of the linearity already established (the maps are linear, so scaling each lane's input bundle by `λ_A` scales the output by `λ_A`), which is why a reasonable user may decide no separate script check is warranted. But the direction (add a script check for the `(1,1/2,-1)` inheritance and the κ1/γ1 closed forms, OR accept that Sec. 5 needs no independent verification and leave the card's Checks list as documentation of a downstream convention) is the user's call, not Codex's.

**Required change:**
See `## Resolve before fix_loop` in the directive. Codex must not auto-resolve this.

**Verification:**
After user resolution: if direction (a), a new `expectZero`/`expect_zero` block verifies `δκ_W^(21) = (1/2)δκ_W^(20)`, `δκ_W^(22) = -δκ_W^(20)` under `δD_{A,n}=ε λ_A D_n^(1)` with `λ=(1,1/2,-1)`, and that `κ1` matches `3(1-σ)/(σ D0)(D2^(1)+D0^(1)/9)`; both engines exit 0. If direction (b), no script change.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:43-127`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:39-129`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The two share identical variable choreography, identical intermediate steps, identical ordering, and even identical SymPy-style idioms transliterated into Mathematica. Three corresponding excerpts:

1. Canonical-branch setup, identical:
   - py L43-46: `u2 = sp.Rational(1,9); u4 = sp.Rational(4,81); D2 = -u2*D0; D4 = -D0/27`
   - wl L43-46: `u2 = 1/9; u4 = 4/81; D2 = -u2*D0; D4 = -D0/27`

2. Series-coefficient extraction, identical:
   - py L52: `du2 = sp.expand(sp.series(u2_full, eps, 0, 2).removeO()).coeff(eps, 1)`
   - wl L52: `du2 = FullSimplify[Coefficient[Normal[Series[u2Full, {eps, 0, 1}]], eps, 1], ...]`

3. The SymPy dummy-symbol Solve indirection, transliterated verbatim into Mathematica (where it is unnecessary — Mathematica `Solve` does not need a `du2sym` placeholder):
   - py L69: `sp.solve(sp.Eq(sp.Symbol('du2sym'), du2_hyb), dkappa)[0].subs(sp.Symbol('du2sym'), du2)`
   - wl L67-70: `dkappa /. First[Solve[du2sym == du2Hyb, dkappa, Reals]] /. du2sym -> du2`

The output `Print` strings are also identical between the two files. This violates the second-engine independence policy: both engines run the same algebra, so an algebra-level error in the chosen derivation path would be reproduced by both rather than caught.

**Why this matters:**
The point of a second engine is to catch a derivation/algebra mistake the first engine's path makes. A transliteration cannot do that; it only catches transcription typos. The stage's identities are all symbolically true (verified by hand above), so no live bug is hidden here, but the cross-check provides far less assurance than the two PASS transcripts imply.

**Required change:**
Re-derive at least the two load-bearing outlet maps (δκ_W, δγ_W) and the even-consistency relation in the `.wl` by an independent route rather than mirroring the `.py` line-for-line: e.g., construct `u2(eps),u4(eps),P0(eps)` and obtain first-order pieces via `D[..., eps] /. eps->0` (a different mechanism than `Series`+`Coefficient`), and invert the hybrid relations algebraically with direct `Solve[du2 == du2Hyb, dkappa]` without the `du2sym` placeholder. The final `expectZero` targets stay identical (they are the paper's boxed formulas); only the derivation path must differ.

**Verification:**
The `.wl` no longer contains the `du2sym`/`dP0sym` placeholder idiom and uses `D[...,eps]` (not `Series`/`Coefficient`) for the linearization; both `expectZero` final checks still print `= 0` and the script exits 0.

## Independent-derivation check (Mathematica)

Not independent. As shown in F2, the `.wl` mirrors the `.py` step-for-step, including the SymPy-specific dummy-symbol `Solve` indirection and identical output strings. The Mathematica file is a transliteration.

## Engine cross-check

Both engines agree. SymPy output and Mathematica output both print `= 0` (or `PASS`) for every `expect_zero`/`expectZero`, and the displayed intermediate forms agree up to surface formatting:
- `a_kappa`: py `-(aD0 + 9*aD2)*(sigma - 1)/(3*D0*sigma)` ≡ wl `-1/3*((aD0 + 9*aD2)*(-1 + sigma))/(D0*sigma)`.
- `a_gamma from map`: py `(-P0*aD0*sigma + P0*aD0 + aN0*sigma - aN0)/(9*D0*P0*sigma)` ≡ wl `((aN0 - aD0*P0)*(-1 + sigma))/(9*D0*P0*sigma)`.
Both exit 0. (Agreement is unsurprising given F2, but the residuals do match.)

## Verdict justification

The math of the eight implemented identities holds up against every attack I tried. I re-derived by hand the canonical branch (`D2=-D0/9`, `D4=-D0/27` giving `u4=4/81`), the first-order `δu2`, `δu4`, `δP0` (all match the boxed notes formulas exactly, including the `5/81` coefficient and the `δu4` sign), the inversion to `δκ_W`/`δγ_W` (the `P0→N0/D0` substitution in the γ check is consistent on both sides and the check is genuine, not tautological), the `δu4=(8/9)δu2` consistency (reduces to `8·3/27 = 8/9`, a real check on the Stage-244 prefactors 1/3 and 8/27), the `δD4 = (2/3)δD2 + (1/27)δD0` solve, and the trace/anomaly cross-checks (derived map vs posited formula — non-tautological). No symbol-domain error invalidates any simplify; no branch is hidden; no hardcoded answer. Verdict is `findings` (not clean) for two reasons: (F1) the paper card and notes enumerate the Section-5 weak-axisymmetric signature `(1,1/2,-1)` and amplitudes κ1/γ1 as a deliverable, and the card's Checks list explicitly requires the signature check, but neither script tests it — a `script_missing_paper_claim` paper_misalignment that needs user resolution before any fix; and (F2) the Mathematica script is a line-by-line transliteration of the SymPy script, including a SymPy-only dummy-symbol Solve idiom, so the second engine does not independently re-derive the result. No stop-cold: the implemented identities are all correct, nothing downstream is invalidated, and F1 is a paper_misalignment which is never auto-CRITICAL_DOWNSTREAM. (Observed but not raised as a categorized finding: both scripts' `banner` and saved transcripts mislabel the stage as "STAGE 153" rather than 170 — a cosmetic label error with no effect on the assertions; flagged as a low-priority relabel in the directive.)

## Self-test notes

Checked variable-independence: the only derivative-like operations are first-order series coefficients in `eps`, and `u2_full/u4_full/P0_full` genuinely depend on `eps` through every `δ`-term, so no identically-zero-derivative trap. No unbounded integrals, so no parity trap. Trivial-case pre-check: substituting the canonical branch (`D2=-D0/9`, `D4=-D0/27`) reduces each residual to the literal paper formula, confirmed by hand to vanish. Paper round-trip: F2's prescribed re-derivation keeps the identical `expectZero` targets (the paper's boxed formulas), so it introduces no new misalignment; F1 is routed to the user, not auto-fixed.
