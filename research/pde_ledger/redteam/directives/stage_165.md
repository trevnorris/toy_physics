---
unit_id: 165
batch: V.1
created_at: 2026-05-28T15:20:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 165

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:39-43`

**Issue:** The headline deliverable δlnL_W = δln a (the stage's `\stagefield{Output}` claim, used immediately as `subs(dLW, da)`) is not actually verified by the SymPy script. Line 41 short-circuits with `if False else 1`, always assigning the literal `1` and never differentiating `LW_law`; line 43 prints the claim as text only. The Mathematica script verifies it genuinely at its line 33. SymPy must do the same so the headline result is exercised by both engines.

**Required change:**
In `scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py`, the current block is:
```python
# 1. D/N law
LW_law = sp.pi * a * sp.sqrt(1 + r**2) / (2 * sp.sqrt(3))
dlog_LW = sp.simplify(sp.diff(sp.log(LW_law), sp.log(a)) if False else 1)  # documentary only
print("D/N law: L_W =", LW_law)
print("At fixed r_* : d ln L_W = d ln a")
```
Replace it with:
```python
# 1. D/N law
LW_law = sp.pi * a * sp.sqrt(1 + r**2) / (2 * sp.sqrt(3))
dlog_LW = sp.simplify(a * sp.diff(sp.log(LW_law), a))
print("D/N law: L_W =", LW_law)
print("At fixed r_* : d ln L_W = d ln a")
expect_zero("d ln L_W - d ln a at fixed r_*", dlog_LW - 1)
```
Notes:
- Differentiate w.r.t. `a` directly and multiply by `a` (logarithmic derivative); do NOT differentiate w.r.t. `sp.log(a)` (SymPy cannot differentiate w.r.t. a composite expression and would error or return 0).
- Use the existing `expect_zero(name, expr)` helper (defined at line 24) — do not invent a new one.
- Keep the two existing `print` lines exactly as-is; add the `expect_zero` line after them.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 165` and confirm the new transcript line `d ln L_W - d ln a at fixed r_* = 0` appears AND the script exits 0.

## Applied: F1
- files_changed: `scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py`
- summary: Replaced the documentary `if False else 1` short-circuit with a genuine logarithmic derivative `dlog_LW = sp.simplify(a * sp.diff(sp.log(LW_law), a))` and added `expect_zero("d ln L_W - d ln a at fixed r_*", dlog_LW - 1)` after the two existing print lines.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:80-84`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:67-71`

**Issue:** The "Stage 164 channel closure" check (D5, δ⊥=0) is tautological. `channel_g`/`channel_r` are the exact left-hand sides of `eq_g`/`eq_r`; substituting `dv_sol`, `dT_sol` (the solution of `solve([eq_r, eq_g], [dv, dT])`) back into them yields 0 by definition of "solution." The `expect_zero`/`expectZero` PASS lines therefore verify only that the solver works, not that the off-family channels vanish for any independent reason. The Mathematica mirror (`/. sol` at lines 68-69) is identically tautological.

**Required change (minimal, safe variant — downgrade to print):**
The fully substantive variant (constructing the Stage-249 δ⊥ normal coordinate as an explicit linear combination of the two channels and asserting it vanishes on-branch while being nonzero off-branch) requires the δ⊥ coefficients, which are not stated in this stage's notes. Unless you can source those exact coefficients from this unit's notes (`notes/stages/moving_throat_pde_stage165_exact_branch_drifts.md` — section 8 does not give them), do NOT fabricate them. Instead, demote the tautological assertions to clearly-labelled prints so the transcript no longer reports a tautology as a verified PASS:

SymPy — replace lines 81-84:
```python
channel_g = sp.simplify((dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW).subs({dv: dv_sol, dT: dT_sol}))
channel_r = sp.simplify((dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW).subs({dv: dv_sol, dT: dT_sol}))
expect_zero("fixed-g channel", channel_g)
expect_zero("fixed-r channel", channel_r)
```
with:
```python
# NOTE: channel_g/channel_r are the LHS of eq_g/eq_r evaluated at the solved (dv,dT).
# They vanish by construction (substituting a linear solution into its own system),
# so this is a consistency print of the solver, not an independent verification of
# delta_perp = 0. Reported, not asserted.
channel_g = sp.simplify((dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW).subs({dv: dv_sol, dT: dT_sol}))
channel_r = sp.simplify((dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW).subs({dv: dv_sol, dT: dT_sol}))
print("fixed-g channel (solver consistency) =", channel_g)
print("fixed-r channel (solver consistency) =", channel_r)
```

Mathematica — replace lines 67-71:
```
banner["Stage 164 channel closure"];
channelG = FullSimplify[(dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW) /. sol, Assumptions -> $Assumptions];
channelR = FullSimplify[(dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW) /. sol, Assumptions -> $Assumptions];
expectZero["fixed-g channel", channelG];
expectZero["fixed-r channel", channelR];
```
with:
```
banner["Stage 164 channel closure"];
(* channelG/channelR are the LHS of eqG/eqR evaluated at the solved (dv,dT); *)
(* they vanish by construction, so this is a solver-consistency print, not an *)
(* independent verification of delta_perp = 0. Reported, not asserted. *)
channelG = FullSimplify[(dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW) /. sol, Assumptions -> $Assumptions];
channelR = FullSimplify[(dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW) /. sol, Assumptions -> $Assumptions];
Print["fixed-g channel (solver consistency) = ", fmt[channelG]];
Print["fixed-r channel (solver consistency) = ", fmt[channelR]];
```
Do NOT change `eqR`, `eqG`, `sol`, `dvSol`, `dTSol`, or any of the substantive `expectZero` calls for the ratio/product/n=5 laws.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 165` and `redteam exec-mathematica 165`; confirm both scripts still exit 0, the substantive ratio/product/n=5 PASS lines are unchanged, and the two former `fixed-g channel`/`fixed-r channel` PASS lines are now plain "(solver consistency)" prints (no PASS:).

## Applied: F2
- files_changed: `scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py`, `mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl`
- summary: Demoted the tautological `fixed-g`/`fixed-r` channel `expect_zero`/`expectZero` assertions to clearly-labelled "(solver consistency)" prints in both engines and added the explanatory comment block; channel_g/channel_r and channelG/channelR definitions and the substantive ratio/product/n=5 checks left unchanged.
- deviation: none
