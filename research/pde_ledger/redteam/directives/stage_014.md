---
unit_id: 014
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 014

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_014.tex:30-33` quote: "Stage~014 exports the mechanism sieve: projected EM mouth data can tune the even gates only away from the degeneracy \eqref{eq:stage007-compensation-denominator}." (The card body, lines 12-28, exhibits only `K_1`, `H_{even}`, and the determinant `Q S_2 - \Delta H_{port}`. No mention of `Xi_load`, `\delta P_2`, `\delta P_4`, or any "5PN bottleneck" / "constant-prefactor transport" content.)
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:50` quote: "014 & Mouth-Taylor gate bridge & \StatusExactClosure{} & Gate conditions for carrying mouth-local projected data into the grouped bundle." (Broad summary; does not name the extra transport content.)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:67` quote: `Xi = sp.simplify((2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)).subs(subs_der)/mu1)`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:74-75` quote: `deltaP2 = sp.simplify(((D0**2*n2 - 2*D0*D2*n0 + 2*D0*N0*z2 - 2*D2*N0*z0)/D0**3).subs(subs_der)/mu1)` and the analogous five-term `deltaP4` definition.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:126-128` quote: load-bearing assertions on `coef_Xi_Px`, `coef_dP2_Gx`, `coef_dP4_Gx`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl:113,117-126,138-140` quote: parallel definitions and `expectZero` / `expectNonzero` assertions M2, M3, M4 for `XiLoad`, `deltaP2`, `deltaP4`.

## Resolve before fix_loop

The paper card's `Output` paragraph and the appendix row both describe Stage 014 narrowly as "the mechanism sieve" for the even gates `K_1, H_{even}` and the compensation degeneracy `Q S_2 - \Delta H_{port}`. The scripts additionally test the 5PN bottleneck `Xi_load = 2 P'/P - ...` and the constant-prefactor transports `\delta P_2`, `\delta P_4`, including specific derivative-coefficient identities (`d(Xi_load)/dP' = 2/P`, `d(\delta P_2)/dG_W' = -2P/(D_0 \Delta^2)`). Which side is canonical for Stage 014?

Possible directions (the user picks one):

- (a) **Paper is incomplete; expand it.** The scripts correctly cover the stage's full content; the paper card was written before the transport ledger was added. Expand `stage_014.tex`'s `Output` paragraph (or add a new paragraph) to enumerate `Xi_load`, `\delta P_2`, `\delta P_4` as additional deliverables of Stage 014, with the constants `2/P` and `-2P/(D_0 \Delta^2)` displayed. No script change.

- (b) **Scripts over-cover; trim them.** Stage 014 really is only the even-gate mechanism sieve, and the `Xi_load`/`\delta P_2`/`\delta P_4` content belongs to a later stage (e.g., a transport-ledger stage). Remove lines 67, 71-75, 93-95, 126-128 from `scripts/moving_throat_pde_stage014_...sympy_audit.py` and the analogous mathematica lines 113, 117-126, 138-140, and migrate them to whichever stage actually claims them. (The user must identify that stage before Codex can act.)

- (c) **Light acknowledgement.** Leave the scripts as-is and add a single sentence to `stage_014.tex` (e.g., between current lines 33 and 34) stating: "The audit script also exercises the constant-prefactor transport coefficients $d \Xi_{\text{load}}/dP' = 2/P$, $d(\delta P_2)/dG_W' = -2P/(D_0 \Delta^2)$, and the dependence of $\delta P_4$ on $G_W'$, which are carried into the downstream grouped-response ledger." This is the lightest-touch resolution.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:139-143`

**Issue:** Lines 142-143 (`assert_zero("compensation K1", K1.subs(comp_surface))` and the H_even analogue) are tautological by `sp.solve`'s contract: `comp_surface` is defined on line 80 as the output of `sp.solve([sp.Eq(K1, 0), sp.Eq(He, 0)], [Hx, Sx], dict=True)[0]`, so substituting it back into the equations Solve was applied to yields zero by construction. The author has self-flagged this on lines 139-141 ("Note: the next two assertions are tautological by sp.solve's contract... kept here for visual symmetry"). The substantive compensation-surface verification is fully captured by the explicit denominator factorization on lines 131-132 (and the Z2-slot identification on lines 145-146, and the sign-flip mutation on line 133). The tautological pair is redundant scaffolding that obscures the assertion inventory.

**Required change:**

Step 1: Delete sympy lines 139-143 inclusive (the three-line comment block on lines 139-141 followed by the two tautological `assert_zero` calls on lines 142-143). Do not modify any surrounding lines.

Concretely, the current text to remove is:
```
# Note: the next two assertions are tautological by sp.solve's contract
# (comp_surface is defined as the solve output), kept here for visual
# symmetry with the explicit denominator factorization assertions above.
assert_zero("compensation K1", K1.subs(comp_surface))
assert_zero("compensation H_even", He.subs(comp_surface))
```

After deletion, line 138 (`raise AssertionError(f"Unexpected pure spectral solve: {sh_only}")`) should be followed directly by what is currently line 144 (`Z2_slot = (Q*S2 - Hport*Delta)/Delta**2`).

Step 2: Verify by inspection that no later line references the removed assertions (they have no side effects; `assert_zero` raises on failure but returns nothing on success). The variables `comp_surface`, `K1`, `He` remain in scope and are not used after the deletion point.

**Self-test for the change:**

- The deletion removes only `assert_zero` calls (no assignments, no variable definitions). `comp_surface`, `K1`, `He` are referenced earlier (lines 80, 68-69) and not used after line 143 in the current file. So the deletion is local.
- The `assert_zero` function (defined lines 24-27) is silent on success, so the saved-output transcript does not contain any "compensation K1" / "compensation H_even" lines today. No transcript line changes after the deletion.
- The substantive content on lines 144-146 (`Z2_slot = ...`, two `assert_zero` calls tying `Hx_den` / `Sx_den` to `-9*Delta^4*Z2_slot` / `-9*Delta^3*Z2_slot`) is untouched.
- No other script in the repo imports from this script, so the deletion has no cross-file impact.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 014` and confirm the script exits 0 and the saved-output transcript at `scripts/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.txt` is byte-identical to the current one (since `assert_zero` is silent on success, the removed assertions emit no transcript content today). The script file should be shorter by exactly 5 lines.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`
- summary: Removed the self-flagged tautological compensation solve round-trip assertions and left the Z2 denominator checks directly after the spectral solve guard.
- deviation: none
