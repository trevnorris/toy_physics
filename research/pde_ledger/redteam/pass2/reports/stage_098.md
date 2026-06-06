---
unit_id: 098
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md]
  paper_appendix: present
---

# Audit unit 098 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md` (only stage-098 notes file present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (relevant rows: line 1230 `\input{stages/stage_098}`; lines 251-258 the blocked support-demand equation; lines 1175, 1186 the result-anchor/conditioning rows)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Purpose}` declares: "Stage 098 is a geometry-lane firewall ledger step. Its audit target is the verification output quoted below." The quoted verdict block (card lines 15-17) is the bottom-line claim: "With `\rho_\alpha=4/3`, any branch with `\zeta_{\max}>1` passes the support test; Family--1 does." The notes elaborate the derivation: with the actual isotropic loading `rho_alpha = 4/3`, the Stage-68 blocked support demand `zeta_req = (rho_alpha-1)/(1 - eps_blk(2-rho_alpha))` collapses to `zeta_req^act(eps_blk) = 1/(3 - 2 eps_blk)` (matched verbatim by appendix eq:app-part04-blocked-zeta-req, lines 252-256). This demand is monotone increasing in `eps_blk`, so on the admissible window `0 <= eps_blk < 1/zeta_max` the inequality `1/(3-2/zeta_max) < zeta_max` reduces to `3 zeta_max(zeta_max-1)/(3 zeta_max-2) > 0`, true for all `zeta_max > 1` (notes §2). Specializing to Family-1's external ceiling `zeta_max^F1 ≈ 2.46752922945601 > 1` gives an edge demand `zeta_req^act < 0.456730991107963 < 2.46752922945601`, so the Family-1 support side is automatic (notes §3, line 96). The card's `\stagefield{Checks}` list (static-limit `c_pole=1/4`, `l=0/l=2` orthogonality, minimal-module hypothesis) is verbatim shared boilerplate across all nine geometry-lane firewall cards 091-099 and is not a stage-098-specific deliverable — those derivations live in their dedicated stages (091/092/...); this stage's own deliverable is solely the support-test automaticity.

## What the script claims to verify

Both scripts verify exactly the support-test theorem. They (1) substitute `rho = 4/3` into the Stage-68 demand and assert it reduces to `1/(3-2 eps)`; (2) compute `d zeta_req/d eps` and assert it equals `2/(3-2 eps)^2` (positive ⇒ monotone increasing); (3) evaluate the worst-case edge demand at `eps = 1/zmax`, form `gap = zmax - zeta_edge`, and assert it factors as `3 zmax(zmax-1)/(3 zmax-2)`; the Mathematica additionally asserts symbolically (`expectTrue`) that this gap is positive for all `zMax > 1` under `$Assumptions = zMax>1 && 0<=epsBlk<1/zMax`. (4) Both then plug the external Family-1 ceiling `zmax_F1 = 2.46752922945601` (flagged in-script as an external carry-forward anchored to the notes) and assert the edge demand `0.45673...` and margin `2.0108...` match pinned targets, with SymPy also asserting `gap_F1 > 0`. The docstring's stated purpose ("actual isotropic branch support demand is automatic for any explicit family with zeta_max > 1 on the admissible blocked interval") matches the card.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| `zeta_req^act = 1/(3-2 eps_blk)` (appendix eq, notes §1) | sympy L11-12 / wl L45,56 assert reduction | match |
| demand monotone increasing in eps_blk (notes §1) | sympy L14-15 / wl L46,57 derivative `= 2/(3-2eps)^2` | match |
| automatic-support inequality `3 zmax(zmax-1)/(3 zmax-2) > 0` for zmax>1 (notes §2) | sympy L18-21 gap factorization + wl L48-49,58-59 factorization AND symbolic `gap>0` | match (wl is stronger: proves positivity for all zmax>1) |
| Family-1 `zeta_max^F1≈2.46752922945601>1`, edge demand `<0.456730991107963` (notes §3, line 96) | sympy L26-34 / wl L64-71 numeric pins + `gap_F1>0` | match |
| card `Checks` block (c_pole=1/4, l=0/l=2 orthog, minimal-module) | none here | n/a — shared boilerplate across 091-099, owned by sibling stages, not a 098 deliverable |

Dominant pattern: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 12 | `simplify(zeta_req - 1/(3-2 eps)) == 0` | demand closed form (appendix eq) | yes |
| A2 | sympy | 15 | `simplify(dz - 2/(3-2 eps)**2) == 0` | monotonicity (notes §1) | yes |
| A3 | sympy | 21 | `simplify(gap_factored - 3*zmax*(zmax-1)/(3*zmax-2)) == 0` | automatic-support theorem (notes §2) | yes |
| A4 | sympy | 29 | `assert gap_F1 > 0` | Family-1 passes (notes §3) | yes |
| A5 | sympy | 33-34 | `abs(zeta_edge_F1 - 0.45673...) < 1e-15`; `abs(gap_F1 - 2.0108...) < 1e-15` | Family-1 numeric (notes line 96) | yes (edge literal in notes; margin derived) |
| A6 | math | 56 | `expectZero[zetaReq - 1/(3-2 eps)]` | demand closed form | yes |
| A7 | math | 57 | `expectZero[dZ - 2/(3-2 eps)^2]` | monotonicity | yes |
| A8 | math | 58 | `expectZero[gap - gapExpected]` | automatic-support theorem | yes |
| A9 | math | 59 | `expectTrue[gap > 0]` (symbolic, all zMax>1) | automatic-support positivity | yes |
| A10 | math | 70-71 | `expectApprox` edge `0.45673...` & margin `2.0108...` (tol 1e-15) | Family-1 numeric | yes |

All rows trace to a specific paper deliverable; no orphaned assertions. No tautological rows: each assertion compares an independently constructed/substituted expression against a separately-stated target form or literal, none of which is the assertion's own LHS by construction.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses Mathematica-native machinery (`$Assumptions` with the admissible window `zMax>1 && 0<=epsBlk<1/zMax`, `FullSimplify`, `expectZero`/`expectTrue`/`expectApprox` harness) rather than echoing SymPy's `sp.simplify`/`sp.diff` call sequence. Crucially it adds a result the SymPy lacks: A9 (`expectTrue[gap > 0]`, wl L59) *symbolically* proves the support gap positive for every `zMax > 1` under the stated assumptions, whereas SymPy only checks `gap_F1 > 0` numerically at the single Family-1 point (sympy L29). Corresponding sections show parallel-but-not-identical algebra: wl L45 builds `zetaReq` directly from `(4/3-1)/(1-epsBlk*(2-4/3))` under assumptions, vs sympy L10-11 builds it from a `Rational(4,3)` symbol; wl L47 substitutes `epsBlk -> 1/zMax`, sympy L18 `subs(eps, 1/zmax)`. The two reach the same four symbolic results (zeta_req form, derivative, edge, gap factorization) by their own routes. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit identical symbolic results (modulo sign-normalized denominators):
- `zeta_req`: sympy `-1/(2*eps-3)` = math `(3-2*epsBlk)^-1` (same).
- derivative: sympy `2/(2*eps-3)**2` = math `2/(3-2*epsBlk)^2` (same).
- gap: sympy `3*zmax*(zmax-1)/(3*zmax-2)` = math `(3*(-1+zMax)*zMax)/(-2+3*zMax)` (same).
- Family-1 numerics agree to ~18 sig figs: edge `0.4567309911079631...` (sympy `...69017835980412` vs math `...649053851860788`, diff ≈ 4.1e-18); margin `2.01079823834804688...` (diff ≈ 4.1e-18). Both within the 1e-15 tolerance; the sub-1e-17 difference is the expected consequence of entering `zmax_F1` at machine-double precision in both engines. All four `expectZero`/`expectTrue` Mathematica checks PASS and the SymPy script prints "STAGE 098 AUDIT PASSED". No `engine_disagreement`.

## Verdict justification

Clean. I read the card, the single notes file, and the relevant appendix rows (251-258, 1175, 1186) before opening the scripts, and the scripts' claim matches the paper's stated audit target exactly. Attacks tried and failed: (a) tautology — each assertion compares an independently-built expression to a separately stated form/literal, not to itself; A1/A3/A8 could fail if `rho_alpha != 4/3` or the factorization were wrong. (b) variable-independence trap on the derivative (A2/A7) — `eps`/`epsBlk` genuinely appears in `zeta_req`, so `dz` is nonzero, not a vacuous zero-derivative. (c) symbol-domain — positivity assumptions (`zMax>1`, `0<=epsBlk<1/zMax`) match the paper's admissible blocked window verbatim, so the `gap>0` simplification is justified. (d) hardcoded-result — the sole literal, `zmax_F1 = 2.46752922945601`, is explicitly flagged in both scripts as an external Family-1 carry-forward anchored to the notes (line 22/78); the edge target `0.45673...` is the notes' own value (line 96) compared against a formula-derived quantity, not self-confirmed. (e) missing card "Checks" — the c_pole=1/4 / l=0,l=2 orthogonality / minimal-module items are verbatim boilerplate shared across all nine 091-099 firewall cards, owned by the dedicated derivation stages, not a 098 deliverable. Outputs are fresh (both `.txt` mtimes ~3h after their scripts). Both engines required (checkpoint:false but pass-2 dual-engine policy) and both present, independent, and agreeing.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_req^act = 1/(3-2 eps_blk)` | sympy out L1; wl out L5 | appendix eq:app-part04-blocked-zeta-req lines 252-256; notes §1 line 40 | MATCH |
| `d zeta_req/d eps = 2/(3-2 eps)^2` | sympy out L2; wl out L6 | notes §1 "monotone increasing" line 42 (form not numerically printed but implied) | MATCH (qualitative carrier) |
| `gap = 3 zmax(zmax-1)/(3 zmax-2)` | sympy out L4; wl out L8 | notes §2 line 62 | MATCH |
| `rho_alpha = 4/3` (input constant) | sympy L10; wl L45 | card line 16, notes lines 18/37, appendix line 251 | MATCH |
| `zeta_max^F1 = 2.46752922945601` | sympy L26; wl L64 | notes lines 22/78/96 | MATCH (declared external carry-forward) |
| Family-1 edge demand `0.456730991107963...` | sympy out L5; wl out L17 | notes line 96 `< 0.456730991107963` | MATCH |
| Family-1 margin `2.01079823834804688...` | sympy out L6; wl out L18 | (absent from .tex and .md) | INTERNAL — derived margin `zmax_F1 - zeta_edge_F1`; cross-engine sync pin, not a stated deliverable |
| `zeta_req^act(0) = 1/3` (zero-blocking) | implied by formula (not separately printed) | notes line 84; card downstream line 27 (`zeta_req=1/3`) | MATCH |

INTERNAL scaffolding (no finding): pass/fail flags, `expectZero`/`expectApprox`/`expectTrue` residual prints, the `1e-15` tolerances, the cross-engine numeric-pin targets `zeta_edge_F1_target`/`gap_F1_target`.

The only non-doc value is the Family-1 margin `2.0108...`, which is a genuine internal derived quantity (the support gap at the F1 ceiling) used to keep the two engines numerically locked; it is not one of the stage's stated deliverables, so it is INTERNAL, not MISSING. All stated deliverable values appear and agree.

reconciliation: complete; 8 deliverable/input values checked, 0 misaligned

## Self-test notes

Checked: (1) variable-independence on A2/A7 — `eps` genuinely appears in `zeta_req`, so the derivative is non-trivially nonzero (not the earlier-unit zero-derivative trap). (2) Trivial-case pre-check on A3/A8 — at `zmax=2`, `gap = 3*2*1/4 = 3/2 > 0` and `zeta_edge = 2/4 = 1/2`, matching `2 - 1/2 = 3/2`; consistent. (3) Positivity-assumption legitimacy — the `gap>0`/`expectTrue` rests on `zMax>1`, exactly the paper's admissible-branch hypothesis, so no hidden over-assumption. No trap fired; no directive written (zero findings).
