---
unit_id: 174
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage174_static_self_similarity.md]
  paper_appendix: present
---

# Audit unit 174 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_174.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage174_static_self_similarity.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 79, 559-609, 715, 1464, 1478)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage174_static_self_similarity_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage174_static_self_similarity_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.txt`

## What the paper claims

Stage 174 rewrites the residual linear grouped `2.5`PN loading defect `\Xi_{\rm load} = N_{01}/N_0 - D_{01}/D_0` (the Stage-241/173 scalar) as a weighted failure of static self-similarity relative to the wall-baseline slope. `\stagefield{Output}` reads verbatim: "Expresses \(\Xi_{\rm load}\) as weighted failure of static self-similarity between wall, BdG, conservative port, and outgoing port slopes." The deliverables, drawn from the notes and appendix eqs (`app-part05-static-slope-defs`, `weights-static`, `Xiload-selfsimilarity`, `Sigma-self-fields`, `Xiload-Sigma-fields`): (D1) weight identity `\omega_K-\omega_B-\omega_Z=1` with `D_0=K-B_0-Z_0`; (D2) `\delta_D=\omega_K\delta_K-\omega_B\delta_B-\omega_Z\delta_Z`; (D3) the wall-referenced form `\Xi_{\rm load}=(\delta_N-\delta_K)+\omega_B(\delta_B-\delta_K)+\omega_Z(\delta_Z-\delta_K)`; (D4-D6) each sector slope is an exact weighted average of microscopic log deformations — `\delta_B=\sum\rho^{(B)}2\,\delta\ln(c/\varpi)`, `\delta_Z=\sum\rho^{(Z)}\delta\ln(Q/\Delta)`, `\delta_N=\sum\rho^{(N)}2\,\delta\ln(P/\Delta)`; (D7) static self-similarity (all defect fields vanish) `\Rightarrow\Xi_{\rm load}=0`; (D8) the scalar mismatch-field collapse `\Xi_{\rm load}=\sum\rho^{(N)}\Sigma^{(N)}+\omega_B\sum\rho^{(B)}\Sigma^{(B)}+\omega_Z\sum\rho^{(Z)}\Sigma^{(Z)}`. All are exact-closure (algebraic) claims inside the weak-axisymmetric grouped bundle.

## What the script claims to verify

The SymPy docstring lists four checks: weighted decomposition of `\delta_D`; wall-referenced decomposition of `\Xi_{\rm load}`; weighted-average formulas for the BdG, conservative Maxwell/mixed, and outgoing-transfer slopes; and the self-similarity theorem `\Xi_{\rm load}=0`. The assertions verify exactly these as symbolic rational-function identities: the weight identity (line 56), the `\delta_D` decomposition (57-60), the wall-referenced form (64), the three two-mode/two-port weighted-average slope identities (84, 103, 122), the self-similarity collapse via `K1,B01,Z01,N01 -> delta_star*(K,B0,Z0,N0)` (130-131), and the scalar mismatch-field form via `B01,Z01,N01 -> (deltaK + Sigma)*(...)` (135-136). The Mathematica `.wl` mirrors all eight checks.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 weight identity `\omega_K-\omega_B-\omega_Z=1` | sympy L56 / math L45 | match |
| D2 `\delta_D` weighted decomposition | sympy L57-60 / math L46-49 | match |
| D3 wall-referenced `\Xi_{\rm load}` form | sympy L62-64 / math L51-56 | match |
| D4 `\delta_B` weighted average (2-mode) | sympy L84 / math L79 | match |
| D5 `\delta_Z` weighted average (2-port) | sympy L103 / math L101 | match |
| D6 `\delta_N` weighted average (2-port) | sympy L122 / math L123 | match |
| D7 self-similarity `\Rightarrow\Xi_{\rm load}=0` | sympy L130-131 / math L130-134 | match |
| D8 scalar mismatch-field collapse | sympy L135-136 / math L136-140 | match |

All eight paper deliverables map to a non-tautological script check, in both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 56 | `expect_zero(omegaK - omegaB - omegaZ - 1)` | D1 | yes |
| A2 | sympy | 57-60 | `expect_zero(deltaD - (omegaK deltaK - omegaB deltaB - omegaZ deltaZ))` | D2 | yes |
| A3 | sympy | 64 | `expect_zero(Xi_load - Xi_wall_ref)` | D3 | yes |
| A4 | sympy | 84 | `expect_zero(B01_two/B0_two - deltaB_weighted)` | D4 | yes |
| A5 | sympy | 103 | `expect_zero(Z01_two/Z0_two - deltaZ_weighted)` | D5 | yes |
| A6 | sympy | 122 | `expect_zero(N01_two/N0_two - deltaN_weighted)` | D6 | yes |
| A7 | sympy | 131 | `expect_zero(Xi_self_similar)` | D7 | yes |
| A8 | sympy | 136 | `expect_zero(Xi_sigma - (SigmaN + omegaB SigmaB + omegaZ SigmaZ))` | D8 | yes |
| A9 | math | 45 | `expectZero[omegaK - omegaB - omegaZ - 1]` | D1 | yes |
| A10 | math | 46-49 | `expectZero[deltaD - (...)]` | D2 | yes |
| A11 | math | 56 | `expectZero[xiLoad - xiWallRef]` | D3 | yes |
| A12 | math | 79 | `expectZero[b01Two/b0Two - deltaBWeighted]` | D4 | yes |
| A13 | math | 101 | `expectZero[z01Two/z0Two - deltaZWeighted]` | D5 | yes |
| A14 | math | 123 | `expectZero[n01Two/n0Two - deltaNWeighted]` | D6 | yes |
| A15 | math | 134 | `expectZero[xiSelfSimilar]` | D7 | yes |
| A16 | math | 140 | `expectZero[xiSigma - (...)]` | D8 | yes |

Every assertion traces to a specific paper deliverable. None is tautological: each defines an object (`D0`, `deltaD`, `B01_two`, etc.) from primitive symbols and then asserts a non-trivial rational-function identity that depends on the cancellation structure (e.g. A1 hinges on `D0 = K-B0-Z0`; A4 hinges on `B01_two` being the exact differential of `B0` so that `B01_two/B0_two` collapses to the weight-averaged log slope). I hand-verified each identity algebraically and all hold.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl:64-79` (and 86-123, 130-140)

**What's wrong:**
The `.wl` is a line-for-line port of the `.py`, not an independent re-derivation. The structural choreography is identical (same banners, same `expect_zero`/`expectZero` helper, same object order, same carry-forward print block), and — more tellingly — the load-bearing intermediate expressions are hand-typed identically rather than re-derived. The "differential of B0" is copied verbatim:

SymPy L78:
`B01_two = ... 2 * c1 * dc1 / w1**2 - 2 * c1**2 * dw1 / w1**3 + 2 * c2 * dc2 / w2**2 - 2 * c2**2 * dw2 / w2**3`

Mathematica L67-71:
`b01Two = ... 2*c1*dc1/w1^2 - 2*c1^2*dw1/w1^3 + 2*c2*dc2/w2^2 - 2*c2^2*dw2/w2^3`

The same verbatim transcription holds for `z01Two` (L89-93 ↔ SymPy L97), `n01Two` (L111-115 ↔ SymPy L116), the weighted-slope candidates (`deltaBWeighted`, `deltaZWeighted`, `deltaNWeighted`), and the two theorem substitutions (L131, L137 ↔ SymPy L130, L135). Even the stale banner text `"STAGE 157 — STATIC SELF-SIMILARITY DECOMPOSITION"` (math L26 ↔ sympy L34) is copied. Because the Mathematica engine re-types the same closed-form differentials as the SymPy engine, both engines confirm the *same hand-entered algebra*; a transcription error in the differential (e.g. a wrong power of `w1`) would be invisible because both files would carry the identical mistake.

**Why this matters:**
The second-engine policy exists so a mistake in one engine's algebra is caught by the other deriving the same fact differently. Here the differentials `B01_two`, `Z01_two`, `N01_two` are the only non-trivial inputs (everything else is a closed-form identity), and they are entered identically in both engines, so the cross-engine check provides no independent confirmation that those differentials are correct — it only confirms that the same typed string equals itself in two CAS engines.

**Required change:**
In the Mathematica script, derive the sector differentials by actual symbolic differentiation instead of re-typing the closed form, so the second engine independently reconstructs the load-bearing inputs. Concretely (perturb each primitive by an infinitesimal and read off the first-order coefficient, or use `D[..., eps]` at `eps -> 0`):
- For BdG: build `b0Eps = (c1 + eps*dc1)^2/(w1 + eps*dw1)^2 + (c2 + eps*dc2)^2/(w2 + eps*dw2)^2` and set `b01Two = D[b0Eps, eps] /. eps -> 0` (replacing the hand-typed L67-71).
- For Z: `z0Eps = (q1 + eps*q1p)/(delta1 + eps*delta1p) + (q2 + eps*q2p)/(delta2 + eps*delta2p)`, `z01Two = D[z0Eps, eps] /. eps -> 0` (replacing L89-93).
- For N: `n0Eps = (p1 + eps*p1p)^2/(delta1 + eps*delta1p)^2 + (p2 + eps*p2p)^2/(delta2 + eps*delta2p)^2`, `n01Two = D[n0Eps, eps] /. eps -> 0` (replacing L111-115).
The existing `expectZero` checks at L79, L101, L123 then remain and should still pass, but now `b01Two`/`z01Two`/`n01Two` are independently derived rather than transcribed. Also fix the stale banner text at L26 from `"STAGE 157 ..."` to `"STAGE 174 ..."` to match the unit.

**Verification:**
After the edit, the Mathematica transcript at `mathematica/output/...stage174...txt` should still print `BdG weighted-average formula = 0`, `Z weighted-average formula = 0`, `N weighted-average formula = 0` with `PASS`, the banner line should read `STAGE 174 ...`, and the script exits 0. The verifier confirms the `b01Two`/`z01Two`/`n01Two` assignments now use `D[..., eps]` rather than the hand-typed four-term forms.

## Independent-derivation check (Mathematica)

The `.wl` is not an independent derivation. It mirrors the `.py` step-for-step: identical helper functions (`expectZero` vs `expect_zero`), identical object-definition order, identical hand-typed differentials (`b01Two` ↔ `B01_two`, `z01Two` ↔ `Z01_two`, `n01Two` ↔ `N01_two`), identical theorem substitutions, identical carry-forward print block, and the same copied "STAGE 157" banner. See F1 for quoted corresponding sections. This is the basis for the `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass all eight checks with zero residuals. The final `Xi_load` forms are displayed differently but are algebraically equal:
- SymPy: `Xi_load = (N0*(-B01 + K1 - Z01) + N01*(B0 - K + Z0))/(N0*(B0 - K + Z0))`
- Mathematica: `Xi_load = n01/n0 - (b01 - k1 + z01)/(b0 - k + z0)`
Expanding the Mathematica form: `n01/n0 - (b01 - k1 + z01)/(b0 - k + z0) = N01/N0 - (B01 - K1 + Z01)/(B0 - K + Z0)`, and `-(B01-K1+Z01)/(B0-K+Z0) = (-B01+K1-Z01)/(B0-K+Z0) = (K1-B01-Z01)/(-(K-B0-Z0)) = D01/(-D0)`... wait: `B0 - K + Z0 = -(K - B0 - Z0) = -D0`, so `-(B01-K1+Z01)/(-D0) = (B01-K1+Z01)/D0 = -(K1-B01-Z01)/D0 = -D01/D0`. So the Mathematica form is `N01/N0 - D01/D0 = deltaN - deltaD`, exactly `Xi_load`. SymPy form expands to the same. The mismatch-field prototypes likewise agree (`(-B0*SigmaB + SigmaN*(B0-K+Z0) - SigmaZ*Z0)/(B0-K+Z0)` = `SigmaN - (B0*SigmaB + SigmaZ*Z0)/(B0-K+Z0)` = `SigmaN + omegaB*SigmaB + omegaZ*SigmaZ`). Engines agree.

## Verdict justification

The math holds up against every attack I tried: the weight identity, `\delta_D` decomposition, wall-referenced `\Xi_{\rm load}` form, the three two-mode/two-port weighted-average slope identities, and both theorem forms (common-slope self-similarity and the scalar mismatch-field collapse) are all genuine, non-tautological rational-function identities that I hand-verified, and they match the paper's eqs `app-part05-static-slope-defs/weights-static/Xiload-selfsimilarity/Sigma-self-fields/Xiload-Sigma-fields` exactly. Sign, factor-of-2 (the BdG and N slopes correctly carry the factor 2 from squared quantities; Z does not, matching the paper), and the `\omega_K = 1+\omega_B+\omega_Z` cancellation that powers D3 all check out. Symbol domains (real, denominators nonzero) are appropriate; `D0` is not declared nonzero but the identities are rational-function-true regardless. Outputs are fresh and the two engines agree. The one defect is that the Mathematica script is a transliteration of the SymPy script — it re-types the same closed-form differentials rather than re-deriving them — so the cross-engine check does not independently confirm the load-bearing `B0`/`Z0`/`N0` differentials. Hence `verdict: findings`, one finding, no stop-cold (the fix is mechanical and changes no assertion or constant, so there is no downstream propagation).

## Self-test notes

Variable-independence: F1's prescribed `D[b0Eps, eps] /. eps->0` — `b0Eps` depends on `eps` through every `(c + eps*dc)` and `(w + eps*dw)` factor, so the derivative is non-trivially nonzero; substituting `eps->0` reproduces the four-term `B01_two` exactly (checked by hand for the `c1`/`w1` pair: `D[(c1+eps*dc1)^2/(w1+eps*dw1)^2, eps]|_{eps=0} = 2 c1 dc1/w1^2 - 2 c1^2 dw1/w1^3`), so A12 still passes. Parity/integral traps: none — this unit has no integrals, only rational identities. Trivial-case pre-check: setting all primed (perturbation) symbols to 0 sends every `expect_zero` residual to `0-0=0` consistently, and the self-similarity substitution drives `Xi_load` to `delta_star - delta_star = 0`; both are genuine, not vacuous. Paper round-trip: the F1 fix introduces no new constant and leaves all eight assertions and the paper-stated forms unchanged, so it cannot create a new paper_misalignment.
