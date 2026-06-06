---
unit_id: 093
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: present
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage093_grouped_p2_status_update.md]
  paper_appendix: present
---

# Audit unit 093 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_093.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage093_grouped_p2_status_update.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows 1175, 1183–1186, 1220)
- sympy: (missing) — card states "SymPy audit: none yet"; `is_status_only_candidate: True`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl`
- sympy output: (missing)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.txt`

## What the paper claims

Stage 093 is a geometry-lane firewall ledger/status step (Part IV, anchor MTDC-T8.1). Its `\stagefield{Purpose}` says "Its audit target is the verification output quoted below," and the `\stagefield{Derivation ledger}` states the computation isolates the forced conservative carrier `\widehat Y_Q^{cons} = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` or the obstruction variables `(eps_2, eps_4)`. The notes (lines 18–55) make the intent precise: under (i) one isotropic grouped-P2 pole and (ii) a geometry lane static through O(omega^4), the conservative quadrupole module is forced to `3/4 + 1/4`, hence the deliverables `rho_alpha = 4/3` and `zeta_req = 1/3`; otherwise the deviation is controlled by the obstruction formula `c_pole = (1 + eps_4) / [4 (1 + eps_2)^2]`. The card's `\stagefield{Checks}` enumerate three: (1) the static limit `eps_2 = eps_4 = 0` returns `c_pole = 1/4`; (2) `l=0`/`l=2` orthogonality before the firewall; (3) a hypothesis-carrying caveat on support/source success. Distinct deliverable values: `c_pole = 1/4`, the `3/4 + 1/4` (i.e. `c_geom = 3/4`) module, `rho_alpha = 4/3`, `zeta_req = 1/3`, plus the obstruction-formula structure as the open-case fallback.

## What the script claims to verify

The `.wl` hardcodes the static limit `eps2 = 0`, `eps4 = 0` (lines 28–29), forms `cPole = (1 + eps4)/(4*(1 + eps2)^2)` (the obstruction formula at the static point, line 30), then derives `cGeom = 1 - cPole`, `rhoAlpha = 1/cGeom`, `zetaReq = cPole/cGeom` (lines 31–33). It prints all four (lines 35–38) and asserts via `expectZero` that each equals its literal deliverable target: `cPole - 1/4`, `cGeom - 3/4`, `rhoAlpha - 4/3`, `zetaReq - 1/3` (lines 40–43). The output transcript confirms all four print as `0` and PASS. So the script verifies: in the static limit, the obstruction formula evaluates to `c_pole = 1/4`, the carrier is `3/4`, and the derived `rho_alpha = 4/3`, `zeta_req = 1/3`.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| Static limit eps2=eps4=0 ⇒ c_pole = 1/4 (Check 1) | `expectZero["c_pole - 1/4", ...]` (line 40), eps hardcoded 0 | match |
| 3/4+1/4 module ⇒ c_geom = 3/4 | `expectZero["c_geom - 3/4", ...]` (line 41) | match |
| rho_alpha = 4/3 | `expectZero["rho_alpha - 4/3", ...]` (line 42) | match |
| zeta_req = 1/3 | `expectZero["zeta_req - 1/3", ...]` (line 43) | match |
| Obstruction formula c_pole = (1+eps4)/[4(1+eps2)^2] structure (open-case) | formula appears at line 30 but only evaluated at eps=0; eps-dependence never exercised | partial |
| Check 2: l=0/l=2 orthogonality before firewall | (none) — belongs to upstream isotropic/decoupling modules | missing (acceptable for status-only consolidation) |
| Check 3: hypothesis-carrying caveat | (prose status; not a numeric check) | n/a |

Dominant pattern: all numeric deliverables MATCH; `paper_alignment: aligned`. The two non-matched rows are the obstruction-formula eps-structure (partially exercised — see F1) and an orthogonality precondition that is an upstream-module concern legitimately absent from this consolidation stage.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | mathematica | 40 | `expectZero[cPole - 1/4]` | Check 1 (static-limit c_pole) | partial — eps hardcoded to 0, so the obstruction formula collapses to arithmetic 1/4−1/4 |
| A2 | mathematica | 41 | `expectZero[cGeom - 3/4]` | 3/4 carrier | partial — `cGeom = 1 - cPole`, fully determined once A1 fixed |
| A3 | mathematica | 42 | `expectZero[rhoAlpha - 4/3]` | rho_alpha deliverable | partial — `rhoAlpha = 1/cGeom`, determined by A2 |
| A4 | mathematica | 43 | `expectZero[zetaReq - 1/3]` | zeta_req deliverable | partial — `zetaReq = cPole/cGeom`, determined by A1/A2 |

All four checks reduce to the same single static-point evaluation: with `eps2 = eps4 = 0`, `cPole` is the literal `1/4` and the remaining three are deterministic functions of it. None of the four can fail for any reason other than a typo in the literal targets; the eps-dependence of the obstruction formula — the genuine content of the firewall block — is never varied.

## Findings

### F1 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl:28-43`

**What's wrong:**
The script hardcodes `eps2 = 0; eps4 = 0;` (lines 28–29) and then forms the obstruction formula `cPole = (1 + eps4)/(4*(1 + eps2)^2)` (line 30). With both perturbations fixed to zero, this expression is the literal constant `1/4`, and `cGeom`, `rhoAlpha`, `zetaReq` are deterministic arithmetic functions of it. The four `expectZero` assertions (lines 40–43) therefore reduce to confirming `1/4 - 1/4`, `3/4 - 3/4`, `4/3 - 4/3`, `1/3 - 1/3` — they can only fail if the literal targets on the right-hand side are mistyped, never if the obstruction formula itself is wrong. The genuine deliverable of the geometry-lane firewall block (MTDC-T8.1, appendix line 1175: "the obstruction variables `(eps_2, eps_4)`, and the forced `3/4+1/4` module") is the *structure* of `c_pole = (1+eps_4)/[4(1+eps_2)^2]` and the fact that the `3/4+1/4` split is recovered **specifically in the static limit** (notes lines 50–55). The script never exercises the formula at nonzero eps, so it does not actually confirm that (a) the formula reduces to `1/4` *only* in the static limit (and not for some other eps pair) or (b) the symbolic obstruction form matches what the notes state. This is the card's stated static-limit target (Check 1) evaluated at one trivial point, which is necessary but weaker than the firewall block's stated content.

Note: this is NOT a `tautological_check` in the strict "defines x=expr then asserts x==expr" sense — the right-hand-side targets (1/4, 3/4, 4/3, 1/3) are independent literals, not echoes of the left-hand-side definitions, so a wrong target literal *would* be caught. It is `insufficient_verification`: the assertions exercise only a single substituted point of a formula whose eps-dependence is the actual claim.

**Why this matters:**
A reader citing MTDC-T8.1 expects the script to anchor the obstruction formula, not merely to evaluate it at the trivial point. As written, the audit would pass identically if the obstruction formula were `(1 + eps4)/(4*(1 + eps2)^3)`, `(1 - eps4)/(4*(1 - eps2)^2)`, or any other expression that happens to equal `1/4` at `eps2 = eps4 = 0`. The static-limit reduction and the symbolic obstruction structure are exactly the firewall block's deliverable.

**Required change:**
Add a symbolic anchor for the obstruction formula that does NOT depend on the static-limit substitution, alongside the existing static-limit checks (keep all four current `expectZero` lines — they correctly report the deliverable values). Specifically:
1. Declare `eps2`, `eps4` as free symbolic parameters (not 0) in a separate block and define `cPoleGen = (1 + eps4)/(4*(1 + eps2)^2)`.
2. Add `expectZero["c_pole static-limit", (cPoleGen /. {eps2 -> 0, eps4 -> 0}) - 1/4]` to confirm the *general* formula reduces to `1/4` exactly at the static point (this is a non-trivial substitution of a symbolic expression, not the pre-collapsed constant).
3. Add a can-fail directional check confirming the formula genuinely deviates off the static point, e.g. assert that `cPoleGen` is NOT identically `1/4`: `If[TrueQ[FullSimplify[cPoleGen - 1/4] === 0], fail["c_pole has spurious eps-independence", ...], pass["c_pole eps-dependent"]]`. (This guards against a degenerate formula that is constant in eps.)
4. Optionally add `expectZero["c_pole at probe eps", cPoleGen /. {eps2 -> 1, eps4 -> 3} - 1]` — at eps2=1, eps4=3 the formula gives `(1+3)/(4*4) = 1/4`... (DO NOT use this probe; see Self-test) — instead use a probe whose value differs from 1/4, e.g. `eps2 -> 1, eps4 -> 0` gives `1/(4*4) = 1/16`, so `expectZero["c_pole probe", (cPoleGen /. {eps2 -> 1, eps4 -> 0}) - 1/16]` confirms the eps-dependence carries the correct algebraic structure.

**Verification:**
The refreshed transcript should show new PASS lines for "c_pole static-limit", "c_pole eps-dependent", and "c_pole probe" in addition to the existing four. The existing four deliverable values (1/4, 3/4, 4/3, 1/3) must remain unchanged. The new checks must be able to fail if the obstruction formula's exponents/signs are altered.

## Independent-derivation check (Mathematica)

Only the Mathematica engine exists for this unit; there is no SymPy script to transliterate from (`missing_sympy` does not apply — see status-only note). The `.wl` derives `c_pole` from the obstruction formula and the three downstream quantities from it; it is not a port of any peer script. No `mathematica_transliteration` concern.

## Engine cross-check

Not applicable — single engine. The card explicitly states "SymPy audit: none yet" and the manifest carries `is_status_only_candidate: True`. This is a status/consolidation stage that recovers the forced `3/4+1/4` module values; the Mathematica-only design is legitimate-by-design here, and the carry-forward values (`rho_alpha = 4/3`, `zeta_req = 1/3`, `c_pole = 1/4`) are the ones cited downstream (card `\stagefield{Downstream use}`: Stages 100–106). No upstream script is referenced for a value this stage's own script fails to verify, so no `missing_sympy` finding is raised.

## Verdict justification

The numeric deliverables all reconcile exactly between script, output, card, and notes (`c_pole = 1/4`, `c_geom = 3/4`, `rho_alpha = 4/3`, `zeta_req = 1/3`), so `paper_alignment: aligned`. The single finding is a low-severity `insufficient_verification`: with `eps2 = eps4 = 0` hardcoded, all four assertions reduce to a single trivial static-point evaluation and never exercise the eps-dependence of the obstruction formula `c_pole = (1+eps_4)/[4(1+eps_2)^2]`, which is the actual stated content of the geometry-lane firewall block (MTDC-T8.1). The right-hand-side targets are independent literals (so it is not strictly tautological), but the formula's symbolic structure is unanchored — the audit would pass for any wrong formula that happens to equal `1/4` at the static point. Attacks tried and failed: (a) checked for value mismatch across card/notes/script/output — none; (b) checked the derivation chain `c_geom = 1 - c_pole`, `rho_alpha = 1/c_geom`, `zeta_req = c_pole/c_geom` — internally consistent and reproduces 3/4, 4/3, 1/3; (c) checked output freshness — `.txt` (14:29) newer than `.wl` (11:12), fresh; (d) confirmed Mathematica-only is by-design per card and status-only manifest flag. Verdict `findings` with one fixable script-side finding.

## Value Reconciliation (pass-2 augmentation)

| value | source (wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `c_pole = 1/4` | wl:30,40; output:5,9 | card Check 1 line 22 ("returns `c_{\rm pole}=1/4`"); notes obstruction formula line 53 (static limit) | MATCH |
| `c_geom = 3/4` (the `3/4` of the `3/4+1/4` module) | wl:31,41; output:6,11 | card line 13 (`3/4 + (1/4)/...`); notes line 25 / appendix line 1175 ("forced `3/4+1/4` module") | MATCH |
| `rho_alpha = 4/3` | wl:32,42; output:7,13 | card `\stagefield{Downstream use}` line 27 (`\rho_\alpha=4/3`); notes line 29 | MATCH |
| `zeta_req = 1/3` | wl:33,43; output:8,15 | card `\stagefield{Downstream use}` line 27 (`\zeta_{\rm req}=1/3`); notes line 31 | MATCH |
| obstruction formula `c_pole = (1+eps4)/[4(1+eps2)^2]` (symbolic) | wl:30 (eps hardcoded 0) | notes line 53; card `\stagefield{Derivation ledger}` line 13 / appendix line 1175 | MATCH (form present in docs; eps-structure under-exercised — see F1, an `insufficient_verification` not a `value_mismatch`) |

INTERNAL (scaffolding, no finding): `eps2 = 0`, `eps4 = 0` (static-limit substitution inputs); intermediate relations `cGeom = 1 - cPole`, `rhoAlpha = 1/cGeom`, `zetaReq = cPole/cGeom` (script's own bridge definitions, not literal in docs but standard for the 3/4+1/4 carrier/pole decomposition); `expectZero`/`pass`/`fail`/`banner`/`fmt` helpers; PASS flags; banner strings.

reconciliation: complete; 5 deliverable values checked, 0 misaligned

## Self-test notes

Variable-independence: no `D[...]` derivatives in the proposed fix; the new checks are substitutions and a non-vanishing assert. Symmetry/parity: no integrals involved. Trivial-case pre-check: verified `cPoleGen /. {eps2->0, eps4->0} = (1+0)/(4*(1+0)^2) = 1/4` (proposed static-limit check residual = 0, correct); verified `cPoleGen /. {eps2->1, eps4->0} = (1+0)/(4*(1+1)^2) = 1/(4*4) = 1/16` so the probe residual `1/16 - 1/16 = 0` is correct and the value 1/16 ≠ 1/4 makes it a genuine off-static probe — I explicitly rejected the eps2->1,eps4->3 probe in the directive because `(1+3)/(4*4) = 1/4` would have re-hit the static value and silently re-introduced the same weakness. The `cPoleGen - 1/4 =!= 0` directional assert can fail iff the formula is eps-independent, which it is not, so it passes as a genuine can-fail guard. Paper round-trip: the fix adds checks only; it leaves the four existing deliverable values (1/4, 3/4, 4/3, 1/3) untouched, so it introduces no new paper_misalignment.
