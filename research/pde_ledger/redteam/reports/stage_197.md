---
unit_id: 197
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage197_conditional_packetA_closure_theorem.md]
  paper_appendix: present
---

# Audit unit 197 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_197.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage197_conditional_packetA_closure_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/blocks at lines 125, 1216–1326, 1467 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The card's `\stagefield{Output}` states verbatim: "Proves \(\Delta_{\rm branch}=0\iff \chi_Q=1\) and hence \(\Delta_Q=\chi_Q-1\) is the sole Packet-A scalar." The card is explicitly a "claim boundary, not a second independent proof" carry-forward stage (line 17): it imports the isotropic grouped-`P2` one-pole front end, the outgoing `l=2` DtN fingerprint, and the natural point-particle source-map reduction. The notes file (which numbers the stage internally as "Stage 248") enumerates five deliverables: (1) collapse of the full Packet-A residual to the single normalization slot once the isotropic front end is imposed; (2) the source-map reduction `Delta_norm = P0^target(chi_Q^{-1}-1)` with `P0^target = 54 G c_s^5/(5 a^5 c^5)`; (3) the equivalence `Delta_branch=0 ⟺ chi_Q=1 ⟺ N_Q=1 ⟺ Delta_Q:=chi_Q-1=0`; (4) the exact deformation-algebra form `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)` and `chi_Q=1 ⟺ 3S(beta^5-1)+Sigma_0+27Sigma_5=0`; (5) higher-odd stability (any tail beginning at `O(omega^7)` leaves `chi_Q` unchanged). The appendix block (lines 1216–1326) restates the same chain and the same closed forms (`eq:app-part05-chiQ-deformation`, `eq:app-part05-packetA-finish-line`). The card states "Mathematica audit: none yet."

## What the script claims to verify

The SymPy script verifies, in seven sections: (I) the grouped anisotropy projectors annihilate an isotropic outgoing input, so `a_{P0}=b_{P0}=0`; (II) the residual vector reduces to `(0,...,0,Delta_norm)`; (III) it constructs `chi_from_series` by series-expanding the carried one-pole/DtN response `L0/(L0 + L2 z^2 + L4 z^4 + i L5 z^5 + i L7 z^7)` (with the matched even moments substituted), extracts the `z^5` coefficient, scales by `-27 i`, and asserts this equals the deformation-algebra closed form `chi_from_def = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)`; (IV) several algebraic identities encoding the zero-set equivalence, the `chi_Q = 1 + Delta_Q` reparametrization, and an explicit anti-vacuity guard that closure FAILS at `chi_Q = 6/5`; (V) the deformation-algebra closure-numerator identity `(3S - Sigma_0)(chi_Q - 1) = 3S(beta^5-1)+Sigma_0+27Sigma_5` and the matching `Delta_norm` form; (VI) the first-order linearized map with coefficients `5, 1/(3S), 9/S`; (VII) `d chi_from_series/dL7 = 0` and `d Delta_norm_pt/dL7 = 0`, encoding higher-odd irrelevance. The output transcript shows exit code 0, all checks zero.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Residual collapse to `(0,...,0,Delta_norm)`; `a_{P0}=b_{P0}=0` | Sec I (`a_{P0}`,`b_{P0}` via projectors) + Sec II (matrix slots) | match (Sec I substantive; Sec II decorative — see notes) |
| (2) `Delta_norm = P0^target(chi_Q^{-1}-1)`, `P0^target=54Gc_s^5/(5a^5c^5)` | Sec III `Delta_norm_pt` defined as `P0_target*(1/chi_from_series - 1)`; checked against itself | match (carry-forward hypotheses; check is definitional — see notes) |
| (3) `Delta_branch=0 ⟺ chi_Q=1 ⟺ N_Q=1 ⟺ Delta_Q=0` | Sec IV check-3 (`chi_Q=1+Delta_Q` reparametrization) + `Delta_norm_bad` anti-vacuity guard at `chi_Q=6/5` | match (genuine) |
| (4) `chi_Q=3(S beta^5+9 Sigma_5)/(3S-Sigma_0)`; `chi_Q=1 ⟺ 3S(beta^5-1)+Sigma_0+27Sigma_5=0` | Sec III identity-1 (`chi_from_series - chi_from_def`) + Sec V both checks | match (genuine, independent-route derivation) |
| (5) Higher-odd `O(omega^7)` irrelevance | Sec VII (`d/dL7` of `chi_from_series` and `Delta_norm_pt`) | match (genuine — `L7` is wired into the series) |
| (5b) Linearized map (notes §5) | Sec VI (`chi_linear`, `Delta_norm_linear`) | match (genuine Taylor check) |

`paper_alignment: aligned`. Every paper-side deliverable maps to a script-side check; the genuinely load-bearing NEW content of the stage (the equivalence and the deformation-algebra forms) is exercised by non-tautological checks.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `expect_zero(a_{P0})` | claim 1 (isotropy ⟹ anisotropy vanishes) | yes |
| A2 | sympy | 56 | `expect_zero(b_{P0})` | claim 1 | yes |
| A3 | sympy | 67 | `expect_zero(Delta_branch[:7] - zeros)` | claim 1 (restatement) | no (defined-zero minus zero; decorative) |
| A4 | sympy | 108 | `chi_from_series - chi_from_def == 0` | claim 4 (key derivation) | yes (independent route) |
| A5 | sympy | 109–112 | `Delta_norm_pt - P0_target(1/chi-1) == 0` | claim 2 | no (X−X; `Delta_norm_pt` defined as exactly this) |
| A6 | sympy | 124–127 | `chi*(Delta_norm/P0_target+1)-1 == 0` | claim 3 | partial (forced by A5's definition) |
| A7 | sympy | 128–131 | `Delta_norm/P0_target+(chi-1)/chi == 0` | claim 3 | partial (forced by A5's definition) |
| A8 | sympy | 135–138 | `Delta_norm(chi=1+Delta_Q)+P0_target Delta_Q/(1+Delta_Q) == 0` | claim 3 (reparametrization) | yes |
| A9 | sympy | 141–142 | `assert Delta_norm(chi=6/5) != 0` | claim 3 (⟸ direction, anti-vacuity) | yes |
| A10 | sympy | 155–158 | `(3S-Sigma0)(chi-1) - closure_num == 0` | claim 4 | yes |
| A11 | sympy | 163–166 | `Delta_norm_def + P0_target closure_num/[3(...)] == 0` | claim 4 | yes |
| A12 | sympy | 185–188 | `chi_linear - [1+eps(5 eps_beta + ...)] == 0` | claim 5b (linearized) | yes |
| A13 | sympy | 199–202 | `Delta_norm_linear + eps P0_target(...) == 0` | claim 5b | yes |
| A14 | sympy | 209 | `d chi_from_series/dL7 == 0` | claim 5 (higher-odd) | yes (L7 wired into series) |
| A15 | sympy | 210 | `d Delta_norm_pt/dL7 == 0` | claim 5 | yes |

## Findings

None. (The decorative/definitional checks A3, A5, A6, A7 are documented below for transparency but do not rise to a defect, because they restate explicitly-carried upstream hypotheses and the stage's actual new content is independently verified by A4, A8, A9, A10, A11, A14, A15.)

### Non-defect notes on the definitional checks

- **A3 (line 67):** `Delta_branch[:7,:] - sp.zeros(7,1)` is `0 - 0 == 0` because the first seven entries of `Delta_branch` were literally constructed as `0` (line 64). This is a tautological restatement of the carried front-end conditions (paper notes §1.1–1.2; appendix Stage 193 freezes `a_2=b_2=a_4=b_4=Delta_pole=0` and Stage 194 isotropy gives `a_{P0}=b_{P0}=0`). The substantive content for the `P0` slots is actually exercised by A1/A2 in Section I. Not a defect: it restates an input the paper explicitly carries, and hides no wrong result.
- **A5 (lines 109–112):** `Delta_norm_pt` is defined (line 98) as `simplify(P0_target*(N_Q - 1))` with `N_Q = 1/chi_from_series` (line 97), and the check compares it to `P0_target*(1/chi_from_series - 1)` — algebraically the same object, so this is `X - X == 0`. The two ingredients (`Delta_norm = P0^target(N_Q-1)` and `N_Q = 1/chi_Q`) are carried-forward hypotheses from upstream Stages 246/195, not new derivations here; the card itself states this stage is a "claim boundary, not a second independent proof." Not a defect.
- **A6/A7 (lines 124–131):** Both reduce to identities forced once `Delta_norm_pt = P0_target(1/chi_from_series - 1)` (A5's definition): A6 collapses to `chi·(1/chi) - 1 = 0`, A7 to `(1/chi - 1) + (1 - 1/chi) = 0`. They confirm clean cancellation but cannot fail given A5. Genuinely tautological, but harmless presentational restatements of claim 3; the real ⟸ test is A9.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit; `mathematica_transliteration` does not apply.

## Engine cross-check

Only one engine present; `engine_disagreement` does not apply.

## Missing-Mathematica judgment (line-114 vs line-118)

I did NOT flag `missing_mathematica`. Per the prompt's line-114 standard (this unit is `is_status_only_candidate: False` but SymPy-only, matching pipeline precedent for stages 121/122/123 and 56 SymPy-only units), a `missing_mathematica` finding is valid only if the SymPy verification is genuinely insufficient as a standalone check AND a second engine is the only route to confidence. Here the stage's actual new content is a set of **pure rational-function and power-series algebraic identities**: the `z^5`-coefficient extraction (A4), the closure-numerator identity (A10–A11), the reparametrization (A8), the linearization (A12–A13), and the `L7`-independence (A14–A15). These are exactly the class of symbolic identities SymPy settles definitively — there is no transcendental special function evaluated numerically, no branch-cut ambiguity, no floating-point tolerance, and no quadrature where a second engine would catch an error. The key derivation (A4) is computed by a genuinely independent route (series expansion of the response function) and reproduces the closed form, and the anti-vacuity guard (A9) confirms the equivalence is not trivially true. SymPy fully and genuinely settles every identity the stage claims; single-engine is acceptable.

## Verdict justification

Verdict `clean`. I read the paper card, the notes file, and the appendix block (lines 1216–1326), built the model of the five deliverables, then attacked the script. Attacks tried and failed: (1) I checked whether the `z^5`-coefficient extraction in Section III is tautological — it is NOT: `chi_from_series` is built from a different representation (a truncated series of the rational one-pole/DtN response with matched even moments) and the asserted equality to `chi_from_def` is a genuine constraint that a wrong factor (e.g. `9` vs `27` in the numerator) would break; hand-derivation confirms `z^5` coeff `= -i L5/L0`, giving `chi_from_series = (3S beta^5 + 27 Sigma_5)/(3S - Sigma_0) = chi_from_def`. (2) I checked the Section VII `d/dL7` checks against self-test trap #1 (differentiating w.r.t. an unwired symbol): `L7` IS wired into `Y_stage194_hi` (line 86, `i*L7*z**7`) and captured by `sp.series(..., 0, 8)`, so the zero derivative is a genuine result (the `z^5` slot does not see the `z^7` tail), not a vacuous one. (3) I checked the equivalence is not vacuously true: the `Delta_norm_bad` guard (A9) explicitly asserts closure FAILS at `chi_Q = 6/5`, output `-9 G c_s^5/(5 a^5 c^5) != 0`, securing the `⟸` direction. (4) I verified all load-bearing constants against the paper: `P0^target = 54 G c_s^5/(5 a^5 c^5)`, `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)`, closure numerator `3S(beta^5-1)+Sigma_0+27Sigma_5`, and linearization coefficients `5, 1/(3S), 9/S` all match the notes and appendix exactly — no value or target mismatch. (5) Output mtime (12:48) is newer than script mtime (11:58): not stale. The only weak checks (A3, A5, A6, A7) are definitional restatements of explicitly-carried upstream hypotheses and hide no wrong result; the stage's genuine new content is fully and non-tautologically verified. The script's claim matches the paper's claim.

## Self-test notes

I walked the required traps. Variable-independence (trap #1): the `d/dL7` checks in Section VII are NOT vacuous — `L7` is genuinely wired into the denominator of `Y_stage194_hi` (line 86) and survives into the order-8 series, so the zero derivative is a real higher-odd-irrelevance result, not a differentiate-an-absent-symbol artifact. Trivial-case / anti-vacuity (trap #3): the `assert Delta_norm_bad != 0` at `chi_Q=6/5` confirms the `assert_zero` cluster is not vacuously satisfiable — the closure genuinely fails off the canonical point, output confirms `-9 G c_s^5/(5 a^5 c^5)`. Paper round-trip (trap #5): all load-bearing constants and closed forms in the script match the notes/appendix verbatim, so no fix is prescribed and no new misalignment is introduced. No directive is written (zero findings).
