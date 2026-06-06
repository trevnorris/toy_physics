---
unit_id: 102
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage102_higher_odd_irrelevance.md]
  paper_appendix: present
---

# Audit unit 102 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_102.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage102_higher_odd_irrelevance.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows at lines 263–295 cover this stage's retarded factorization block)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.txt`

## What the paper claims

Stage 102 is a retarded 2.5PN factorization ledger step. The card's boxed claim (stage_102.tex:16) is: "Retarded terms beginning at \(O(\omega^7)\) do not affect the point-particle \(2.5\)PN theorem." The notes (lines 18–32) make this concrete: starting from the most general one-pole grouped-`P2` retarded denominator with one extra higher-odd tail, `Yhat_Q^ret(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5 - i tau_Q omega^7) + O(omega^8)`, the omega^5 expansion is `1 + omega^2/(4 Omega^2) + omega^4/(4 Omega^4) + i chi_Q 9/(32 Omega^5) omega^5 + O(omega^6)`, and the `tau_Q` term first appears only at `O(omega^7)`. The appendix (part04, eqs at lines 275/281/295) pins the carried constants: `sigma_Q^can = 9/(8 Omega_Q^5)`, the odd coefficient `Gamma_5 = chi_Q · 9 K_0/(32 Omega_Q^5)`, and the verbatim irrelevance statement at line 295. Distinct deliverables: (D1) the omega^5 (2.5PN) coefficient is independent of `tau_Q`; (D2) `tau_Q` first surfaces at omega^7; (D3) the canonical omega^5 odd coefficient equals `chi_Q · 9/(32 Omega^5)`.

## What the script claims to verify

The SymPy docstring (lines 2–14) states three checks, mirrored by the Mathematica script. It builds `Y = 3/4 + (1/4)/(1 - omega^2/Omega^2 - I·chiQ·sigma_can·omega^5 - I·tauQ·omega^7)` with `sigma_can = (9/8)/Omega^5`, series-expands to O(omega^5) and O(omega^7), then asserts: D1 — `d/dtauQ Im[coeff(omega^5)] = 0` (tauQ cannot back-propagate to omega^5); D2 — `d/dtauQ Im[coeff(omega^7)] = 1/4` (tauQ's leading imaginary appearance); D3 — `Im[coeff(omega^5)] - chiQ·9/(32 Omega^5) = 0`. The Mathematica `.wl` does the identical three via `Coefficient[..]/I` instead of `sp.im`. All three are can-fail `expect_zero`/`expectZero` gates, not prints.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| D1: tauQ irrelevant at omega^5 (card line 16, appendix line 295, notes line 30) | sympy D1 line 51 / wl line 46 | match |
| D2: tauQ first surfaces at omega^7 (notes line 30) | sympy D2 line 52 / wl line 47 (coefficient pinned to 1/4 = the outer 1/4 prefactor) | match |
| D3: omega^5 odd coeff = chiQ·9/(32 Omega^5) (appendix eq line 281, notes line 28) | sympy D3 lines 53–54 / wl line 48 | match |
| sigma_Q^can = 9/(8 Omega^5) (appendix eq line 275, notes line 22) | script literal `sigma_can = (9/8)/Omega**5` (sympy line 38 / wl line 33) | match |

Every paper-side deliverable maps to a substantive script check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 51 | `expect_zero(tau5)` where `tau5 = diff(im(coeff(omega^5)), tauQ)` | D1 | yes |
| A2 | sympy | 52 | `expect_zero(tau7 - 1/4)` | D2 | yes |
| A3 | sympy | 53–54 | `expect_zero(im(coeff(omega^5)) - chiQ·9/(32 Omega^5))` | D3 | yes |
| A4 | mathematica | 46 | `expectZero[tauCoeff5]` (= `D[Coeff(omega^5)/I, tauQ]`) | D1 | yes |
| A5 | mathematica | 47 | `expectZero[tauCoeff7 - 1/4]` | D2 | yes |
| A6 | mathematica | 48 | `expectZero[Coeff(omega^5)/I - chiQ·9/(32 omegaQ^5)]` | D3 | yes |

All six trace to a specific deliverable and none is tautological: each asserts a property of the freshly series-expanded `Y`, so a wrong expansion would make them fail rather than pass trivially.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses Mathematica-native machinery throughout: `Series[yRet, {omega, 0, 5}]` + `Normal` + `Expand` (wl lines 35–36), and extracts the imaginary coefficient via `Coefficient[ySer5, omega, 5]/I` (wl lines 41, 48) — a genuinely different route from SymPy's `sp.series(...).removeO()` + `sp.im(...coeff(omega,5))` (py lines 40–41, 46, 54). The verification harness also differs: `expectZero` uses `FullSimplify[Together[Expand[...]]]` under `$Assumptions` (wl lines 20–24) versus the Python `expect_zero` `sp.simplify` + `!= 0` (py lines 20–32). Both build the same physical premise (the one-pole denominator with the `tauQ omega^7` tail) — that is required by the paper, not a sign of copying — but the algebra to extract and check the coefficients is engine-native on both sides. The `/I` versus `sp.im` divergence is the cleanest evidence of independence. No `mathematica_transliteration` finding.

## Engine cross-check

The two saved outputs agree exactly at the level they claim:

- omega^5 series: sympy `1 + omega**2/(4*Omega**2) + omega**4/(4*Omega**4) + 9*I*chiQ*omega**5/(32*Omega**5)` (sympy txt line 1); mathematica `1 + ((9*I)/32)*chiQ*omega^5/omegaQ^5 + omega^4/(4*omegaQ^4) + omega^2/(4*omegaQ^2)` (wl txt line 5) — identical up to term ordering and symbol name `Omega`↔`omegaQ`.
- omega^7 series: both carry the extra `I*omega^7*tauQ/4` and `9*I*chiQ*omega^7/(16*Omega^7)` terms (sympy txt line 2, wl txt line 6).
- `tauQ` coeff at omega^5 = 0 in both (txt line 3 / line 7); `tauQ` coeff at omega^7 = 1/4 in both (txt line 4 / line 8).
- All three gates PASS with residual 0 in both engines (sympy txt lines 5–7, wl txt lines 9–14).

`engines_agree: true`.

## Verdict justification

Clean. I read the card, notes, and appendix block (part04 lines 263–295) before the scripts and confirmed the script's verified identities are exactly the stage's deliverables. Attacks tried and failed: (1) tautology — the three checks assert properties of a freshly series-expanded `Y`, so they can fail if the expansion is wrong; the D2 "1/4" is the legitimate outer prefactor of `Y = 3/4 + (1/4)/D`, not an arbitrary pin. (2) hardcoded coefficient — `sigma_can = (9/8)/Omega^5` matches appendix eq line 275 verbatim, and the resulting `9/32` matches eq line 281 and notes line 28. (3) symbol-domain abuse — `Omega, chiQ, tauQ` are declared positive/real (py line 36, wl lines 29–31), which is exactly the physical setup and is what makes `sp.im` / `/I` extract the purely-imaginary odd coefficients correctly; no positivity is smuggled in to force a pass. (4) coverage gap — omega^6 is even and outside the odd-term claim, so its absence of an explicit check is not a defect; the load-bearing omega^5 (D1) and omega^7-surfacing (D2) checks are present. (5) stale output — both `.txt` mtimes (14:28 / 14:31) are newer than their script mtimes (11:12 / 11:12); outputs are fresh and match what the current scripts would produce. Both engines required (checkpoint: False, status-only: False) are present, substantive, and independent.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `sigma_Q^can = 9/(8 Omega^5)` | py line 38, wl line 33 | appendix eq line 275; notes line 22 | MATCH |
| omega^5 odd coeff `= chi_Q·9/(32 Omega^5)` | py txt line 1, wl txt line 5; py lines 53–54, wl line 48 | appendix eq line 281 (`Gamma_5 = chi_Q·9 K_0/(32 Omega_Q^5)`); notes line 28 | MATCH |
| D1: tauQ irrelevant at omega^5 (= 0) | py txt lines 3,5; wl txt lines 7,10 | card line 16; appendix line 295; notes line 30 | MATCH |
| D2: tauQ first surfaces at omega^7, coeff 1/4 | py txt lines 4,6; wl txt lines 8,11 | notes line 30 (surfacing order; 1/4 itself is script scaffolding for the outer 1/4 prefactor) | MATCH (claim = surfacing order; value reconciles) |
| `Delta_Q = chi_Q - 1` | not script-emitted (conceptual finish) | notes line 46; appendix eq line 293 | MATCH (doc-only, consistent) |

INTERNAL (scaffolding, no finding expected in prose): omega^7 chiQ coefficient `9 chiQ/(16 Omega^7)` (printed intermediate term in the O(omega^7) series, not a stated deliverable); the even omega^2/omega^4/omega^6 conservative coefficients `1/(4 Omega^2)`, `1/(4 Omega^4)`, `1/(4 Omega^6)` (intermediate conservative-branch terms, not this stage's deliverables); pass/PASS flags and residual=0 verification values.

reconciliation: complete; 5 deliverable values checked, 0 misaligned.

## Self-test notes

I checked variable-independence (D1/D2 differentiate w.r.t. `tauQ`, which genuinely appears in `Y` via the `i·tauQ·omega^7` denominator term, so the derivatives are not trivially zero — D1 = 0 is a real cancellation from `tauQ` being absent at omega^5, D2 = 1/4 is the real surfacing coefficient). Parity: the odd-power (omega^5, omega^7) coefficients are purely imaginary and the even-power ones purely real, which is exactly why `sp.im` / `/I` cleanly isolate them; no symmetric-domain integral is involved. Trivial-case: substituting the canonical premise, the omega^5 coefficient reduces to `9 i chiQ/(32 Omega^5)` so D3's residual is 0 and D2's tauQ-coefficient is the literal 1/4 — both confirmed by the fresh outputs. No directive written (zero findings).
