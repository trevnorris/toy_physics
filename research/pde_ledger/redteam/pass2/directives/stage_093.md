---
unit_id: 093
batch: IV.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T05:10:34Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 093

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond the named edit. Do NOT touch paper.tex, notes/, or any prose documents — the red-team only modifies scripts.

After editing, RUN the script (`math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl`) and iterate until it exits 0 with all in-file checks passing.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl:28-43`

**Issue:** The script hardcodes `eps2 = 0; eps4 = 0;` then forms the obstruction formula `cPole = (1 + eps4)/(4*(1 + eps2)^2)`. With both perturbations zero, this is the literal `1/4`, and `cGeom`, `rhoAlpha`, `zetaReq` are deterministic functions of it. The four `expectZero` checks (lines 40–43) therefore reduce to a single trivial static-point evaluation and never exercise the eps-dependence of the obstruction formula — which is the actual deliverable of the geometry-lane firewall block (MTDC-T8.1, appendix line 1175: "the obstruction variables `(eps_2, eps_4)`, and the forced `3/4+1/4` module"; notes line 53: `c_pole = (1 + eps_4)/[4 (1 + eps_2)^2]`). The audit would pass identically for any wrong formula that happens to equal `1/4` at the static point.

**Required change:**
KEEP all existing lines unchanged (lines 28–46 stay; the four deliverable checks correctly report 1/4, 3/4, 4/3, 1/3). ADD a new symbolic-anchor block AFTER line 43 (before the `Print[""]` at line 45). The new block must:

1. Define the general (symbolic-eps) obstruction formula with free symbols:
```mathematica
ClearAll[e2, e4];
cPoleGen = (1 + e4)/(4*(1 + e2)^2);
```
   (Use fresh symbol names `e2`, `e4` so the earlier `eps2 = 0; eps4 = 0;` bindings do not leak in.)

2. Confirm the GENERAL formula reduces to the static-limit value `1/4` exactly at the static point (a non-trivial substitution of a symbolic expression):
```mathematica
expectZero["c_pole static-limit (symbolic)", (cPoleGen /. {e2 -> 0, e4 -> 0}) - 1/4];
```

3. Add a can-fail directional guard that the formula is genuinely eps-dependent (NOT identically `1/4`):
```mathematica
If[TrueQ[FullSimplify[cPoleGen - 1/4] === 0],
  fail["c_pole spurious eps-independence", FullSimplify[cPoleGen - 1/4]],
  pass["c_pole eps-dependent"]];
```

4. Add an off-static probe whose target value is NOT `1/4`, anchoring the formula's algebraic structure:
```mathematica
expectZero["c_pole off-static probe (e2=1,e4=0)", (cPoleGen /. {e2 -> 1, e4 -> 0}) - 1/16];
```
   Rationale: at `e2 -> 1, e4 -> 0` the formula gives `(1+0)/(4*(1+1)^2) = 1/(4*4) = 1/16 ≠ 1/4`, so this probe fails if the exponent (the `^2`) or the numerator/denominator structure is wrong. DO NOT pick a probe such as `e2 -> 1, e4 -> 3` — that re-hits `1/4` and would silently re-introduce the same weakness.

**Self-test confirmation (already done by auditor):**
- `(cPoleGen /. {e2->0,e4->0}) = (1+0)/(4*(1+0)^2) = 1/4` → check 2 residual = 0. Correct.
- `(cPoleGen /. {e2->1,e4->0}) = (1+0)/(4*(1+1)^2) = 1/16` → probe residual `1/16 - 1/16 = 0`, and `1/16 ≠ 1/4`. Correct, genuine off-static point.
- `FullSimplify[cPoleGen - 1/4]` is not identically 0 (it depends on e2,e4), so the directional guard takes the `pass` branch — a real can-fail check.

**Verification command:**
After Codex applies, the verifier runs the Mathematica script and confirms the new PASS lines "c_pole static-limit (symbolic)", "c_pole eps-dependent", "c_pole off-static probe (e2=1,e4=0)" appear in the transcript AND the four original deliverable PASS lines (c_pole − 1/4, c_geom − 3/4, rho_alpha − 4/3, zeta_req − 1/3) are unchanged AND the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl`
- summary: Added symbolic-eps anchor checks for the obstruction formula, including static-limit, eps-dependence, and off-static probe verification.
- deviation: none
