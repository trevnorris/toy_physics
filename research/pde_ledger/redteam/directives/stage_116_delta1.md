---
unit_id: 116
batch: IV.3
iteration: 2
created_at: 2026-05-28T20:50:00-06:00
parent_directive: stage_116.md
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-28T20:58:38-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex delta-directive — unit 116 (iteration 2)

Iteration 1 applied F1 and F2. F2 is fine. **F1 regressed on re-run**: `exec-mathematica 116` now exits 1 with

```
D/N eigenvalue solves Cos[kW lW]==0 = Cos[kSym*lW]
FAIL: D/N eigenvalue solves Cos[kW lW]==0
  residual -> Cos[kSym*lW]
```

i.e. `kWValue` was left unbound (still the symbol `kSym`). Root cause: `Solve[Cos[kSym*lW] == 0 && 0 < kSym*lW < Pi, kSym]` solves *for* `kSym` while `lW` is symbolic, so Solve must divide by `lW` and returns an empty / `ConditionalExpression` result non-deterministically across `math -script` runs. (Your iteration-1 run happened to bind it; the re-run did not.)

## F1-fix — make the eigenvalue solve deterministic

Replace the eigenvalue-extraction lines (the `kWValue = ...` Solve/Reduce line, and keep the two `expectZero` checks) so the solver acts on the **product** `u = kSym*lW` over the bounded open interval `(0, Pi)` — there is exactly one root `u == Pi/2`, and no division by symbolic `lW` is needed. Then recover `kWValue = u/lW`. Suggested form (adapt syntax so the script runs and exits 0):

```
(* The Neumann BC q'(lW)=0 on the q(0)=0 mode gives the characteristic
   equation Cos[kW lW] == 0. Solve for the product u = kW*lW on (0, Pi):
   exactly one root u == Pi/2 (deterministic; no division by symbolic lW). *)
uRoot = u /. First[Solve[Cos[u] == 0 && 0 < u < Pi, u]];
kWValue = uRoot / lW;
expectZero["D/N eigenvalue solves Cos[kW lW]==0", Cos[kWValue*lW]];
expectZero["D/N eigenvalue kW = Pi/(2 lW)", kWValue - Pi/(2*lW)];
```

If `Solve` with the inequality is itself unreliable in this environment, the acceptable robust fallback is `Reduce[Cos[u] == 0 && 0 < u < Pi, u]` (which returns `u == Pi/2`), then extract that root and set `kWValue = (Pi/2)/lW` from the *solved* `u` (NOT typed directly as `Pi/(2 lW)`). The DSolve general-solution comment block from iteration 1 may stay or be trimmed — the load-bearing requirement is unchanged: `kWValue` must come from a solver acting on the characteristic equation `Cos[u]==0`, and the script must exit 0 deterministically.

RUN `math -script` on the file (and re-run it a second time to confirm determinism) until it exits 0 with both eigenvalue checks PASS. Do not touch the SymPy script. Do not alter the F2 renormalization block (it is correct). Append an `## Applied: F1-fix` block when done.

## Applied: F1-fix

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl`
- summary: Replaced symbolic-k root extraction with a deterministic bounded solve for u = kW*lW and recovered kWValue from the solved product.
- deviation: none
