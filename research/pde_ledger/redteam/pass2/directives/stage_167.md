---
unit_id: 167
batch: V.1
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T22:07:40Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 167

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit only the file:line range named. Do NOT touch `paper.tex`, `notes/`, or the SymPy script — only the Mathematica `.wl` may change.

After editing, RUN `math -script mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl` and iterate until it exits 0 with all in-file PASS lines printing.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl:47-61`

**Issue:** The `.wl` reaches `dv, dT, dgq, dgs, dIsq, dlam` through the exact same algebraic assembly as the SymPy script (`scripts/...stage167...py:53-65`) — identical variable names, identical RHS decompositions (`(dZ-drho)/2 + 3*dcsw/2 ± dcs - (5|3)*da/2`), and the identical internal name `dIsq`/`I_sq`. This is a structural port: a shared sign/coefficient error in the assembly would be reproduced in both engines and would escape the engine cross-check. The current `.wl` does add genuine extra closure assertions (the carry-forward `expectZero`s at lines 53–54, 55–56, 68–70), so this is low-severity, not a blocking echo — but the second engine should reach the invariants/channels by a route that does NOT replay the SymPy assembly.

**Required change:**
Add an independent NUMERIC second route in the `.wl` that does not depend on the symbolic assembly path. Concretely, after the existing symbolic checks (after line 87, before the print-summary block at line 89), insert a numeric-substitution closure that picks several concrete drift tuples and confirms every invariant/channel vanishes when computed from the ground-truth primitive definitions:

1. Build a substitution-rule generator from the primitive definitions directly (not from the already-assembled `drc/dr/dg/chan1/chan2`). For example, define a pure function that takes numeric `(dTheta, dKs, dKq, dP)` and recomputes `drhoN, daN, dcsN, dZN, dcswN, dellN, dLWN, dvN, dTN, dgqN, dgsN, dIsqN, dlamN` numerically via the SAME physical formulas, then forms `r_c, frak r, frak g, chan1, chan2, deltaPerp` numerically.
2. For at least three independent tuples (e.g. `{1,0,0,0}`, `{0,1,0,0}`, `{0,0,1,0}`, `{0,0,0,1}`, and one mixed `{2,-3,5,-7}` with `gstar->Rational[1,3]`, `rstar->2`), assert each of `{r_c, frak r, frak g, chan1, chan2, deltaPerp}` evaluates to exactly `0`. Use `expectZero["delta ln r_c numeric (tuple k)", value]` etc. so the verifier sees new PASS lines.

This gives the Mathematica engine a verification route (numeric evaluation over multiple sampled drifts) that is genuinely independent of the symbolic-assembly echo, while leaving the existing symbolic block intact so nothing already-passing regresses.

(If a full numeric block is judged too large, the minimal acceptable change is: for the single mixed tuple `{dTheta->2, dKs->-3, dKq->5, dP->-7, gstar->1/3, rstar->2}`, evaluate `drc, dr, dg, chan1, chan2, deltaPerp` and `expectZero` each — five-plus new PASS lines that exercise the invariants through concrete substitution rather than symbolic cancellation.)

**Self-test (auditor-verified):** With `{dTheta=2, dKs=-3, dKq=5, dP=-7}`: `drho=1`, `da=-3/2-1/2=-2`, `dcs=-3/2-1/2-7/5=-17/5`, `dZ=5+14/5=39/5`, `dcsw=2`, `dell=-2`, `dLW=-2`, `dIsq=2(-2)+(-2)+(-1)= -7`, `dv=(39/5-1)/2+3+(-17/5)-5(-2)/2 = (34/5)/2+3-17/5+5 = 17/5+8-17/5 = 8`, `dlam=dv+dIsq=8-7=1`; then `r_c = 2(1)-(-3)-5 = 2+3-5 = 0` ✓, `frak r = 1-((-3)+5)/2 = 1-1 = 0` ✓, `chan2 = (-3)+5-2(1) = 0` ✓. The invariants vanish for this nonzero mixed drift, confirming the numeric route is non-trivial (the residual would be nonzero if any upstream coefficient were wrong) and not tautological.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 167` and confirms the new numeric `expectZero` PASS lines appear in `mathematica/output/...stage167...txt` AND the script exits 0. The assembly block at `wl:47-61` may remain, but the invariants are now also confirmed through a non-port numeric route.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl`
- summary: Added an independent exact numeric closure block that recomputes primitive drifts for five tuples and checks all invariant/channel residuals with `expectZero`.
- deviation: none
