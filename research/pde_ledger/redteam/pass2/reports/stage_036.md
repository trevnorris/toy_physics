---
unit_id: 036
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage036_support_feasibility_frontier.md]
  paper_appendix: present
---

# Audit unit 036 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_036.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage036_support_feasibility_frontier.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 62)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`

## What the paper claims

Stage 036 (Part II, checkpoint anchor MTDC-T5) isolates the support-feasibility condition that remains after the Stage-035 normalization locus has selected `xi_req`. `\stagefield{Output}` reads verbatim: "Stage~036 outputs the dimensionless mixed baseline \eqref{eq:app-stage036-Mmix}, the support-feasibility function \eqref{eq:app-stage036-G}, the required support loading \eqref{eq:app-stage036-gBreq}, the monotonicity derivative \eqref{eq:app-stage036-G-derivative}, and the final admissibility test \eqref{eq:app-stage036-final-test}." The five boxed deliverables are: (1) `M_mix := 8 alpha_mix/(pi^2 A) = 8 Chi^2/(pi^2 A Omega_U^2 Delta_0)`; (2) `G(xi,delta) := 8 alpha_req/(pi^2 A) = 9 xi(xi+delta)/(9 delta + 11 xi)`; (3) `g_{B,req}^2/varpi^2 = (pi^2 A/8)[G(xi_req,delta) - M_mix]`; (4) `dG/dxi = 9(9 delta^2 + 18 delta xi + 11 xi^2)/(9 delta + 11 xi)^2 > 0`; (5) the final test `R_target >= 1`, `F(xi_req,delta) = R_target`, `M_mix <= G(xi_req,delta)`. The card and notes also state endpoints `G(0,delta)=0`, `G_max = 9(1+delta)/(9 delta + 11)`, the parametric frontier `xi -> (F,G)`, the near-onset expansion `G = xi - 2 xi^2/(9 delta) + O(xi^3)`, and a corollary `M_crit ~= (R_target-1)/(1+8/(9 delta))` (not boxed, not in the Output list). The notes mirror all of this.

## What the script claims to verify

Both engines verify: (1) the closed form `G = 9 xi(xi+delta)/(9 delta + 11 xi)` equals `8 alpha_req/(pi^2 A)`; (2) the factorization `g_{B,req}^2/varpi^2 - (pi^2 A/8)(G - M_mix) = 0`; (3) `dG/dxi` equals the manifestly-positive closed form (SymPy) / the positivity numerator polynomial `11 xi^2 + 18 delta xi + 9 delta^2` with discriminant `-72 delta^2 <= 0` (Mathematica); (4) endpoints `G(0,delta)=0`, `G_max = 9(1+delta)/(9 delta + 11)`; (5) the final-test support-inequality ↔ nonnegativity identity at `xi_req`; (6) admissible/inadmissible numeric witnesses for `R_target >= 1` and `M_mix <=> G`; (7) the genuine non-tautological anchor — a symbolic reconstruction of `R_target_sym` from the Stage-18 microscopic kappa expansion (`kappa0^2 = 8/pi^2`, `kappa1^2 = 16/(9 pi^2)`), confirming `F - R_target_sym = 0`; and (8) the near-onset series `G = xi - 2 xi^2/(9 delta)`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `M_mix = 8 Chi^2/(pi^2 A Omega_U^2 Delta_0)` | sympy L56 / wl L49 `Mmix_expr`, printed | match |
| (2) `G = 9 xi(xi+delta)/(9 delta+11 xi) = 8 alpha_req/(pi^2 A)` | sympy L72-75 `expect_zero` on factorization; wl L57-58 `G - 8 alpha_req/(Pi^2 A)`, `G - closed form` | match |
| (3) `g_{B,req}^2/varpi^2 = (pi^2 A/8)(G - M_mix)` | sympy L72-75 / wl L63-66 | match |
| (4) `dG/dxi = 9(9 d^2+18 d xi+11 xi^2)/(9 d+11 xi)^2 > 0` | sympy L78-81; wl L72-85 (polynomial + discriminant) | match |
| (5) final test `R_target>=1`, `F(xi_req,d)=R_target`, `M_mix<=G` | sympy L92-95,104,122-125,146-155; wl L104-107,116-154 | match |
| Endpoints `G(0,d)=0`, `G_max=9(1+d)/(9d+11)` | sympy L82-86; wl L86-96 | match |
| Parametric frontier `xi->(F,G)` | sympy L88-95 (printed + identity); wl L98-107 | match |
| Near-onset `G = xi - 2 xi^2/(9 d)` | sympy L158-161; wl L156-166 | match |
| Corollary `M_crit ~= (R_target-1)/(1+8/(9 d))` | (not boxed, not in Output list; derived from G≈xi + Stage-035 onset) | not required (not a stated deliverable) |

Dominant pattern: every boxed/Output deliverable has a faithful, non-tautological script-side check in both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 72-75 | `expect_zero(gBreq - (pi^2 A/8)(G-Mmix))` | (2),(3) | partial — definitional given hardcoded a_req,G (author notes this; anchor is A8) |
| A2 | sympy | 78-81 | `expect_zero(dG - dG_target)` | (4) | yes — diff(G) computed, not declared |
| A3 | sympy | 82 | `expect_zero(G.subs(xi,0))` | endpoint | yes |
| A4 | sympy | 83-86 | `expect_zero(Gmax - closed form)` | endpoint | yes — limit computed |
| A5 | sympy | 92-95 | `expect_zero` final-test identity at xi_req | (5) | partial — definitional |
| A6 | sympy | 104,146-155 | `expect_true` admissible/inadmissible numeric | (5) | yes — numeric witnesses |
| A7 | sympy | 122-125 | `expect_zero(F(host) - R_target_host)` | (5) middle conjunct | yes — host built from kappa, not F |
| A8 | sympy | 142-145 | `expect_zero(F - R_target_sym)` | (2)/(5) genuine anchor | yes — symbolic kappa rebuild |
| A9 | sympy | 158-161 | `expect_zero(G_series - target)` | near-onset | yes — series computed |
| B1 | mathematica | 57-58 | `expectZero(G - 8 alphaReq/(Pi^2 A))`, `G - closed form` | (2) | yes |
| B2 | mathematica | 63-66 | `expectZero(gBreq - (Pi^2 A/8)(G-Mmix))` | (3) | partial — definitional |
| B3 | mathematica | 76-79 | `expectZero(dG polynomial - (11 xi^2+18 d xi+9 d^2))` | (4) | yes — D[G] then polynomial reduction |
| B4 | mathematica | 82-85 | `expectZero(disc + 72 d^2)` | (4) positivity | yes — Discriminant computed |
| B5 | mathematica | 86,96 | `expectZero(G(0))`, `expectZero(Gmax - closed form)` | endpoints | yes — Limit computed |
| B6 | mathematica | 104-107 | `expectZero` final-test identity at xiReq | (5) | partial — definitional |
| B7 | mathematica | 116,145-154 | `expectTrue` admissible/inadmissible numeric | (5) | yes — witnesses |
| B8 | mathematica | 126-129 | `expectZero(F(host) - R_target_host)` | (5) middle conjunct | yes — host from kappa |
| B9 | mathematica | 141-144 | `expectZero(F - R_target_sym)` | (2)/(5) genuine anchor | yes — symbolic kappa rebuild |
| B10 | mathematica | 159-166 | `expectZero` series coeffs c0=0,c1=1,c2=-2/(9 d) | near-onset | yes — coeffs read from Series |

## Findings

### F1 — stale_output

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt` (mtime 2026-05-26 00:44:48)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt` (mtime 2026-05-26 00:44:54)

**What's wrong:**
Both saved-output `.txt` files have mtimes (2026-05-26) older than their source scripts (both 2026-06-03 15:59:11). The staleness is visible in content, not just timestamps: the SymPy script now prints banners `STAGE 36.1 … STAGE 36.4` (`.py` L63,77,88,157) and a final line `All Stage 19 checks passed.` (L163), while the committed SymPy output still shows the OLD banners `STAGE 19.1 … STAGE 19.4` (`.txt` L3,12,21,32) and `All Stage 19 checks passed.` (L36). The Mathematica script prints `STAGE 036 — …`, `STAGE 036.3`, `STAGE 036.4` (`.wl` L31,98,156) and `Stage 036 Mathematica audit passed.` (L169), but the committed output shows `STAGE 019 …`, `STAGE 019.3`, `STAGE 019.4` (`.txt` L3,29,45) — the banner labels were relabeled in the scripts after the outputs were captured. The mathematical content of every assertion in the saved outputs still agrees with what the current scripts assert (all checks `= 0` / PASS), so this is purely a label-relabel staleness, not a content disagreement.

**Why this matters:**
The committed transcript no longer matches the banner labels the current scripts emit, which is a provenance inconsistency. It does not affect any verified result.

**Required change:**
Re-run both scripts to refresh the saved outputs (`python3 scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py > scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`; `math -script mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl > mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`). No script logic change is required. (Per project policy the orchestrator's independent exec re-run refreshes these; this finding is informational.)

**Verification:**
Refreshed SymPy output should show `STAGE 36.1 … 36.4` banners and the same `= 0` residual lines; refreshed Mathematica output should show `STAGE 036 …` banners and the same `PASS:` lines. All assertions still pass.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses materially different verification moves than the `.py`:
- dG/dxi positivity: SymPy compares `diff(G,xi)` to a declared `dG_target` closed form (`.py` L78-81); Mathematica instead multiplies `D[gTarget,xi]` through by `(9 delta+11 xi)^2/9` to expose the polynomial `11 xi^2+18 delta xi+9 delta^2` (`.wl` L72-79) AND independently computes `Discriminant[...] = -72 delta^2 <= 0` (`.wl` L82-85) to prove non-negativity — a check absent from the SymPy script.
- Near-onset: SymPy asserts `series(G) - (xi - 2 xi^2/(9 delta)) = 0` against a declared target (`.py` L158-161); Mathematica reads `c0,c1,c2` out of the `Series` via `Coefficient` and asserts each separately (`.wl` L159-166), declaring no target series.
- The genuine anchor (A8/B9) is parallel in intent but each engine rebuilds `R_target_sym` from the kappa literals in its own CAS (`.py` L129-145, `.wl` L136-144).

These are independent strategies for the same claims. No `mathematica_transliteration` finding.

## Engine cross-check

The engines agree on every shared result. `G = 9 xi(delta+xi)/(9 delta+11 xi)` (py `.txt` L5, wl `.txt` L5); `M_mix = 8 Chi^2/(pi^2 A Delta_0 Omega_U^2)` (py L7, wl L7); `R_target = pi^2 A NQ/(8 beta0)` (py L8, wl L8); `G_max = 9(1+delta)/(9 delta+11)` (py L17, wl L24); near-onset `xi - 2 xi^2/(9 delta)` (py L34, wl L47). Numeric witnesses match: admissible `M_mix=0.02795…, G=0.46551…`, inadmissible `M_mix=0.81056…` (py L28-29, wl L39,41). Sample `R_target = 1414562/558009 ≈ 2.535 >= 1` agrees in both. No `engine_disagreement`.

## Verdict justification

The stage holds up against the paper. Every boxed Output deliverable (M_mix, G, gBreq factorization, dG/dxi positivity, final admissibility test) plus the endpoints, parametric frontier, and near-onset series is checked in both engines, and the constants/forms match the card and notes exactly (F from Stage 035 L47-48; alpha_req from Stage 035 L104-105). I attacked the obvious tautology risk — F, G, and a_req are hardcoded closed forms, so the factorization identities (A1/A5/B2/B6) are definitional — but the authors flagged this explicitly and provided a genuine non-tautological anchor (A7/A8, B8/B9) that rebuilds `R_target_sym` from the Stage-18 microscopic kappa literals (`8/pi^2`, `16/(9 pi^2)`) and confirms it equals F; this is the load-bearing check and it is not circular. The dG/dxi positivity is independently witnessed by the discriminant `-72 delta^2 <= 0` in Mathematica. The only finding is `stale_output`: the committed transcripts carry the old `STAGE 19/019` banner labels while the current scripts emit `STAGE 36/036`, with no mathematical content disagreement. Verdict is `findings` solely on that informational staleness; no script logic is wrong.

## Self-test notes

I checked: (1) variable-independence of every `diff`/`D` — `diff(G,xi)` and `D[gTarget,xi]` genuinely depend on `xi` (G is a non-trivial rational in xi), so dG is not identically zero; the discriminant is taken in `xi` of a genuine quadratic. (2) The `xi -> 1^-` limit is one-sided (`dir='-'` / `Direction->"FromBelow"`), correctly matching the open interval `xi in [0,1)` and the `F -> +infty` pole at xi=1 (handled by taking the limit of G, which is finite there). (3) Trivial-case: `G(0,delta)=0` substitution and the near-onset coefficients (c0=0,c1=1,c2=-2/(9 delta)) reduce as claimed; numeric admissible/inadmissible witnesses straddle G=0.4655 as the saved outputs confirm. (4) No missing-script directive (both engines present). (5) Paper round-trip: the only fix (refresh outputs) introduces no new constant and no paper-side change, so no new paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

Every deliverable value the scripts emit was traced to the `.tex` card and/or `.md` notes. All reconcile.

| value | source (py/wl + output line) | .tex / .md location | status |
|---|---|---|---|
| `G(xi,delta) = 9 xi(xi+delta)/(9 delta+11 xi)` | py L54,64 / wl L43,53; out py L5, wl L5 | tex L34-37 (eq G, boxed); md L21-22,61-62 | MATCH |
| `M_mix = 8 Chi^2/(pi^2 A Omega_U^2 Delta_0)` | py L56,66 / wl L49,55; out py L7, wl L7 | tex L26-29 (eq Mmix, boxed); md L29-30,56-57 | MATCH |
| `g_{B,req}^2/varpi^2 = (pi^2 A/8)(G - M_mix)` | py L73-74 / wl L64-65 | tex L42-44 (boxed); md L25,66 | MATCH |
| `dG/dxi = 9(9 d^2+18 d xi+11 xi^2)/(9 d+11 xi)^2` | py L79 / wl L72-79 (poly form); out py L15-16, wl L14 | tex L58-60 (boxed); md L79 | MATCH |
| `G(0,delta) = 0` | py L82 / wl L86; out py L22, wl L16 | tex L65; md L84 | MATCH |
| `G_max(delta) = 9(1+delta)/(9 delta+11)` | py L84 / wl L96; out py L17, wl L24 | tex L67-68; md L86-87 | MATCH |
| `R_target = pi^2 A NQ/(8 beta0)` (script form `NQ A/(beta0·(8/pi^2))`) | py L57 / wl L50; out py L8, wl L8 | (intermediate; R_target=F is the carrier, tex L84/L107) — value `pi^2 A NQ/(8 beta0)` not separately boxed in 036 card; consistent with Stage-035 normalization product | INTERNAL (intermediate; defined via Stage 035) |
| `F(xi,delta) = (9d+11xi)^4/[81(1-xi)(9d^2+18d xi+11xi^2)^2]` | py L53 / wl L44-47; out py L6, wl L6 | carried from Stage 035 tex L47-48; referenced in 036 tex L79,L107 frontier/final-test | MATCH (carry-forward, Stage 035) |
| near-onset `G = xi - 2 xi^2/(9 delta)` | py L159 / wl L158-166; out py L34, wl L47-48 | tex L91-92; md L128 | MATCH |
| final admissibility test (`R_target>=1`, `F(xi_req,d)=R_target`, `M_mix<=G`) | py L104,122-125,146-155 / wl L116,126-154 | tex L104-109 (boxed); md L117-118 | MATCH |
| parametric frontier `xi -> (F(xi,d), G(xi,d))` | py L90 (printed) / wl L98 | tex L77-79 (boxed); md L108 | MATCH |

INTERNAL items (scaffolding/witnesses, no prose expected): kappa0Sq=`8/pi^2`, kappa1Sq=`16/(9 pi^2)` (microscopic anchor literals); a_req=`9 pi^2 A xi(xi+delta)/[8(9d+11xi)]` (Stage-035 carry-forward, matches Stage 035 tex L104-105); host sample values `A_host=3, beta0_host=5, delta=1, xi=1/2`; numeric witnesses `M_mix=0.02795…/0.81056…`, `G=0.46551…`, `R_target=1414562/558009`; discriminant `-72 delta^2`; alpha_mix=`Chi^2/(Omega_U^2 Delta_0)`.

reconciliation: complete; 11 deliverable values checked, 0 misaligned (R_target line is an intermediate defined via Stage 035, treated INTERNAL; F is a verified Stage-035 carry-forward, MATCH).
