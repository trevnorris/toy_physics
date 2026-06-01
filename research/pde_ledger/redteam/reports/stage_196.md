---
unit_id: 196
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage196_higher_odd_irrelevance_theorem.md]
  paper_appendix: present
---

# Audit unit 196 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_196.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage196_higher_odd_irrelevance_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at line 123 and lines 1318-1326)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card `\stagefield{Output}` states verbatim: "Proves all extra isotropic odd tails beginning at \(O(\omega^7)\) are invisible to the point-particle \(2.5\)PN theorem." The appendix row (line 123) repeats this and line 1318 adds that such tails "do not alter any coefficient through the point-particle \(2.5\)PN order." The notes enumerate five distinct deliverables: (1) an exact response-side higher-odd difference identity `Y7 - Y5 = H_Q/[4(1-X_Q)(1-X_Q-H_Q)] = O(omega^7)`; (2) an exact DtN-side difference identity `Ydef7 - Ydef5 = -L0 L_{>=7}/(D5 D7) = O(z^7)`; (3) the theorem that the outgoing-normalization factor `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)` is unchanged by all higher odd DtN data beginning at `O(z^7)`; (4) the consequence that the Stage-246 source-map-reduced Packet-A residual `Delta_norm = P0_target(chi_Q^{-1} - 1)` (with `P0_target = 54 G c_s^5/(5 a^5 c^5)`) is likewise unchanged; and (5) the final reduced statement that the only live retarded obstruction at 2.5PN is `Delta_Q := chi_Q - 1`. The notes also fix the Stage-245 canonical-even matching `-L2/L0 = 1/9`, `L2^2/L0^2 - L4/L0 = 4/81` and the resulting `Sigma_2`, `Sigma_4` laws.

(Provenance note on stage labels: the notes refer to the upstream stages as 242/244/245/246 and to "this stage" as 247; the script's banner and comments say "STAGE 179 / Stage 194 / Stage 195"; the paper card is "Stage 196." These are three different numbering schemes for the same physics content. See the cosmetic observation in Verdict justification — it does not affect any verified identity.)

## What the script claims to verify

The SymPy script is organized into five sections that map one-to-one onto the five notes deliverables. Section I builds the closed forms `Y5 = 3/4 + (1/4)/(1-X_Q)` and `Y7 = 3/4 + (1/4)/(1-X_Q-H_Q)` with `X_Q = omega^2/Omega_Q^2 + i chi_Q sigma_can omega^5`, `sigma_can = 9/(8 Omega_Q^5)`, `H_Q = i tau_Q omega^7`, and asserts the exact difference identity, the vanishing of the difference through `O(omega^5)`, the canonical one-pole form through `O(omega^5)`, and that the first correction starts at `omega^7`. Section II does the same on the DtN side with `D5 = L0 + L2 z^2 + L4 z^4 + i L5 z^5` and `D7 = D5 + i L7 z^7`. Section III substitutes the Stage-245 parameterization, applies the canonical-even matching laws, extracts `chi_Q` from the matched `z^5` coefficient (computed from a denominator that still contains `i L7 z^7`), confirms it equals the deformation-algebra closed form, and confirms it is `L7`-independent. Section IV builds `N_Q = 1/chi_Q` and `Delta_norm` and confirms their `L7`-independence and the `P0_target` relation. Section V confirms `L7` enters the normalized DtN response first at `z^7` with coefficient `-i/L0`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Response-side difference identity `= H_Q/[4(1-X_Q)(1-X_Q-H_Q)]`, `= O(omega^7)` | lines 59, 72; series form lines 74, 78 | match |
| Burke-Thorne `omega^5` coeff `i chi_Q (9/32) omega^5/Omega^5` | line 73-77 (`expected_Y5_through5`) | match |
| (2) DtN difference identity `= -L0 L_{>=7}/(D5 D7)`, `= O(z^7)` | lines 103, 104; series form line 109 | match |
| (3) `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)` unchanged by `L7` | lines 157, 158, 159 | match |
| Canonical-even matching `1/9`, `4/81`; `Sigma_2`,`Sigma_4` laws | line 157 (consequence check) | match |
| (4) `Delta_norm` unchanged by `L7`; `P0_target = 54 G c_s^5/(5 a^5 c^5)` | lines 167, 175, 176, 177 | match |
| (5) `Delta_Q := chi_Q - 1` only obstruction | definitional corollary of (3)+(4); ledger lines 207-209 | match (no separate identity needed) |
| `L7` first enters at `z^7`, coeff `-i/L0` | lines 184-187 | match |

`paper_alignment: aligned`. Every numbered deliverable has a corresponding non-tautological script-side check; the load-bearing `chi_Q` closed form, `P0_target`, and matching constants all match the notes verbatim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(Y7-Y5) - Ydiff_expected == 0` | claim 1 (response identity) | yes |
| A2 | sympy | 72 | `series(Y7-Y5,0,7) == 0` | claim 1 (vanishes through omega^5) | yes |
| A3 | sympy | 74 | `series(Y7,0,6) - expected_Y5_through5 == 0` | claim 1 (canonical one-pole form + BT slot) | yes |
| A4 | sympy | 78 | `Ydiff_series - i tau_Q omega^7/4 == 0` | claim 1 (first correction at omega^7) | yes |
| A5 | sympy | 103 | `simplify(Ydef7-Ydef5) - Ydef_diff_expected == 0` | claim 2 (DtN identity) | yes |
| A6 | sympy | 104 | `series(Ydef7-Ydef5,0,7) == 0` | claim 2 (vanishes through z^5) | yes |
| A7 | sympy | 109 | `series(Ydef7,0,6) - (1 - L2 z^2/L0 + ... ) == 0` | claim 2 (normalized DtN form) | yes |
| A8 | sympy | 157 | `Y_stage194_through5 - expected_through5 == 0` | claim 3 (1/9, 4/81 matching + chi_Q slot) | yes |
| A9 | sympy | 158 | `chi_from_series - chi_from_def == 0` | claim 3 (chi_Q closed form from matched series) | yes |
| A10 | sympy | 159 | `diff(chi_from_series, L7) == 0` | claim 3 (chi_Q is L7-independent) | yes (L7 wired in at line 129) |
| A11 | sympy | 175 | `diff(N_Q_pt, L7) == 0` | claim 4 (N_Q L7-independent) | partial (redundant; see below) |
| A12 | sympy | 176 | `diff(Delta_norm_pt, L7) == 0` | claim 4 (Delta_norm L7-independent) | partial (redundant; see below) |
| A13 | sympy | 177 | `Delta_norm_pt - P0_target*(1/chi - 1) == 0` | claim 4 (Delta_norm form / P0_target) | yes |
| A14 | sympy | 187 | `L7_coeff_in_Y + i/L0 == 0` | L7 first enters at z^7, coeff -i/L0 | yes |

A11/A12 are marked "partial" only because, by the time `N_Q_pt`/`Delta_norm_pt` are built (lines 168-169) from the already-`L7`-free `chi_from_series`, the differentiation target literally contains no `L7` symbol, so each derivative is `0` for a trivial reason. They are not the load-bearing evidence: the genuine `L7`-independence is established upstream by the exact DtN difference identity A5/A6 (`Ydef7-Ydef5 = O(z^7)`) and by A9/A10, where `chi_Q` is extracted from a series whose denominator (line 129) still carries `i L7 z^7` and is shown to be `L7`-free only because `L7` first appears at `z^7` (confirmed by A14 and output lines 103-106, 122-128). A11/A12 are therefore harmless redundant confirmations, not insufficient verification.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit (paper card line 11: "Mathematica audit: none yet"). Not applicable.

## Engine cross-check

Only one engine present; not applicable.

## Verdict justification

I attacked every assertion and could not break the stage. Specific attacks tried and failed:

- **Tautology on the difference identities (A1, A5):** `Ydiff_expected` (line 51) and `Ydef_diff_expected` (line 97) are constructed independently from the closed forms, not by re-arranging `Y7-Y5`; an error in the boxed formula would fail the assert. Not tautological.
- **Burke-Thorne factor (A3):** I hand-derived `Y5 = 1 + X/4 + X^2/4 + ...` with `X` containing `i chi_Q (9/(8 Omega^5)) omega^5`, giving the `omega^5` coefficient `i chi_Q 9/(32 Omega^5)`, matching `expected_Y5_through5` exactly. Real check.
- **chi_Q closed form (A9):** I hand-derived `chi_from_series = -27 i * z5_coeff = -27 i * i(S beta^5/9 + Sigma_5)/(3S - Sigma_0) = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)`, matching `chi_from_def` and the notes formula. Real check.
- **The variable-independence trap (A10):** `diff(chi_from_series, L7)` is NOT a vacuous "differentiate w.r.t. an unwired symbol" — `L7` is genuinely wired into `Y_stage194_hi` (line 129, `+ i L7 z^7` in the denominator). The theorem content is precisely that the `z^5` coefficient is `L7`-free, which the series extraction demonstrates (output: `L7` first appears at `z^7`). So `diff = 0` is the correct, physically meaningful outcome of a check that could have failed if `L7` had leaked into a lower-order coefficient.
- **Hardcoded `Sigma_2`/`Sigma_4` (lines 133-134):** these are taken from the notes rather than solved from the matching conditions, but A8 checks the *consequence* (`z^2` coeff `= 1/9`, `z^4` coeff `= 4/81`), which would fail if the laws were wrong. Anchored by downstream consistency, not a free pass.
- **Matching of all constants:** `P0_target = 54 G c_s^5/(5 a^5 c^5)` (line 167) matches notes line 297; `sigma_can = 9/(8 Omega_Q^5)` matches notes; `chi_Q` form matches notes line 274.

On the missing-Mathematica question, I applied the prompt's line-114 judgment rather than a reflexive line-118 demand. Every claim this stage makes is an *exact symbolic rational-function identity* or a *series coefficient of a rational function* — operations for which SymPy's `simplify`/`series`/`coeff` are deterministic and complete. There are no integrals, branch cuts, transcendental evaluations, or numerical-conditioning hazards where a second CAS would plausibly catch a SymPy-specific defect; a `.wl` would almost certainly be a transliteration of the same rational algebra. SymPy fully and genuinely settles the identities the stage claims, so single-engine is acceptable here (matching the cited precedent of SymPy-only verified stages 121/122/123). I therefore did NOT flag `missing_mathematica`.

One cosmetic, non-blocking observation: the script banner prints "STAGE 179" (lines 35, 189) and its comments reference "Stage 194/195" (lines 116, 118, 162, 164), while the notes use "Stage 242/244/245/246/247" and the paper card is "Stage 196." This is stale internal labeling across three numbering schemes. It does not change a single verified identity — every equation in the script faithfully matches the corresponding notes equation — so per the no-stylistic-changes / no-scope-expansion constraints it is recorded here as an observation only and is NOT filed as a directive finding or a `paper_misalignment` (the verified claim matches the paper claim exactly; only the printed labels differ).

Output freshness: script mtime `May 11 11:58`, output mtime `May 11 12:48` — output is newer than the script, exit code 0, all in-file `expect_zero` checks print `= 0`. Not stale.

Verdict: **clean**. The script's five sections faithfully and non-tautologically exercise all five notes deliverables and the paper card's `\stagefield{Output}`, with all load-bearing constants matching the paper/notes.

## Self-test notes

I checked the variable-independence trap on all three `diff(..., L7)` checks: A10 (line 159) is genuine because `L7` is wired into the line-129 denominator and the theorem is that it drops out by `z^5`; A11/A12 (lines 175-176) are trivially zero because their targets are built from the already-`L7`-free `chi_from_series`, but they are redundant confirmations, not the sole evidence (the genuine proof is the exact DtN difference identity plus the `z^5`-coefficient extraction). No integrals are present, so the parity/symmetry trap is not applicable. Trivial-case pre-check: every `expect_zero` target is an exact symbolic identity confirmed to simplify to literal 0 in the committed output, and the hardcoded `Sigma_2`/`Sigma_4` are guarded by the downstream `1/9`, `4/81` consequence check.
