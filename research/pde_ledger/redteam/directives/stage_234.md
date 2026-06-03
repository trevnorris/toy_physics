---
unit_id: 234
batch: VII.2
created_at: 2026-06-02T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T22:06:35-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 234

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the targets named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_mathematica_audit.wl` (NEW file — create it in `mathematica/`, the `_mathematica_audit.wl` suffix is mandatory)

**Issue:** Unit 234 has only a SymPy engine; it is not a checkpoint and not status-only, so the dual-engine rule applies. Mathematica can independently verify every deliverable of this stage. You must write a NEW Mathematica script that derives and checks the claims below independently — using native Mathematica primitives via a DIFFERENT decomposition than the SymPy script, NOT a line-by-line transliteration. In particular: do NOT echo the SymPy `D[expr, eps] /. eps -> 0` linearization choreography (use `Series[expr, {eps, 0, 1}]` + `Normal`/`Coefficient` instead) and do NOT echo the explicit Lagrange-multiplier `Solve` for the balanced family (use `Minimize` instead). Use `FullSimplify`/`Simplify` with `Assumptions` for the log-chart identities. Strip any `ConditionalExpression[0, ...]` wrappers from simplifier output before the zero test, and verify each in-file check with an explicit `If[FullSimplify[lhs - rhs] =!= 0, Print["FAIL ..."]; Exit[1]]` pattern (or equivalent) so the script exits nonzero on any failure.

Declare symbols with the same physical domains the SymPy script uses: `C_star, B_star, Rtr_ref, Rtarget_ref, eps_eta_ref, eps_eta_star, Rtr, Rtarget, eps_eta` positive; `q_tr, q_nt, q_eta, r1, R1, E1, eps, xi` real; `vareps` real and nonzero. Use `c_eta = eps_eta_star/(1 - eps_eta_star)`.

**Claim manifest** (each must be an independent in-file check):

- **M1 — exact finite quotient chart + inverse roundtrip.** With
  `q_tr = -C_star*Log[Rtr/Rtr_ref]`,
  `q_nt = B_star*Log[Rtr/Rtr_ref] + Log[(1-eps_eta)/(1-eps_eta_ref)] - Log[Rtarget/Rtarget_ref]`,
  `q_eta = Log[eps_eta/eps_eta_ref]`,
  and the inverse map
  `Rtr_inv = Rtr_ref*Exp[-q_tr/C_star]`,
  `eps_eta_inv = eps_eta_ref*Exp[q_eta]`,
  `Rtarget_inv = Rtarget_ref*Exp[-q_nt - (B_star/C_star)*q_tr]*(1-eps_eta_inv)/(1-eps_eta_ref)`,
  substituting the inverse map back into the forward chart returns `q_tr`, `q_nt`, `q_eta` respectively (residual `== 0` under positive assumptions).

- **M2 — first-order linearization.** With `Rtr = Rtr_ref*Exp[eps*r1]`, `Rtarget = Rtarget_ref*Exp[eps*R1]`, `eps_eta = eps_eta_star*Exp[eps*E1]` and `eps_eta_ref -> eps_eta_star`, the order-`eps` coefficients (via `Series`/`Coefficient`, NOT `D[...]/.eps->0`) satisfy
  `q_tr^(1) = -C_star*r1`,
  `q_nt^(1) = B_star*r1 - c_eta*E1 - R1`,
  `q_eta^(1) = E1`.

- **M3 — triangular compiler + inverse drift map.** With `Theta1 = -q_tr^(1)/C_star`, `Xi1 = q_nt^(1) + (B_star/C_star)*q_tr^(1)`, `Rcal1 = -Xi1 - c_eta*q_eta^(1)`:
  `Theta1 == r1`, `Xi1 == -R1 - c_eta*E1`, `Rcal1 == R1`, and the inverse `-(1 - eps_eta_star)/eps_eta_star * (Rcal1 + Xi1) == E1`.

- **M4 — exact R_tr cancellation.** `D[Xi1, r1] == 0` (equivalently, `FullSimplify[Xi1]` contains no `r1`). This is the load-bearing physics claim; confirm it independently (form `Xi1 = q_nt^(1) + (B_star/C_star) q_tr^(1)` directly from M2's coefficients, do not reuse a pre-simplified `-R1 - c_eta E1`).

- **M5 — rigid-mouth reduction.** At `r1 -> 0`: `Theta1 == 0` and `Xi1 == q_nt^(1)|_{r1=0}` (i.e. `q_nt = Xi1` on the rigid-mouth slice). Also confirm the finite rigid-mouth form `q_nt|_{Rtr=Rtr_ref} == Log[(1-eps_eta)/(1-eps_eta_ref)] - Log[Rtarget/Rtarget_ref]`.

- **M6 — two-observable strip form + ceiling ordering.** With `robust = 0.367930328492646`, `nonempty = 0.737619063660757`, and `Xi1 = -R1 - c_eta*E1`: substituting `R1 = -c_eta*E1 - robust/Abs[vareps]` gives `Xi1 == robust/Abs[vareps]` (and the `+` edge gives `-robust/Abs[vareps]`); same for `nonempty`. ADDITIONALLY assert `0 < robust < nonempty` so a corrupted/ swapped ceiling fails the check (see F2 — mirror that ordering guard here).

- **M7 — canonical direct-branch families.** Pure-target `(R1,E1)=(-xi,0)` gives `Xi1 == xi`; pure-dressing `(R1,E1)=(0,-xi/c_eta)` gives `Xi1 == xi`; and the balanced minimal-norm family obtained via `Minimize[{R1^2 + E1^2, -R1 - c_eta*E1 == xi}, {R1, E1}]` yields the minimizer `(R1,E1) = (-xi/(1+c_eta^2), -c_eta*xi/(1+c_eta^2))`, which also satisfies `Xi1 == xi`.

**Self-test confirmation (already done by auditor):**
- M2's `Series` route: each linearized expression genuinely depends on `eps` via the `Exp[eps*...]` substitutions, so the order-1 coefficient is non-trivially nonzero — no zero-coefficient trap.
- M4's `D[Xi1, r1]`: `Xi1` formed from M2 coefficients genuinely contains `r1` in both `q_nt^(1)` (`+B_star*r1`) and `q_tr^(1)` (`-C_star*r1`); the cancellation is real and contingent on the exact `B_star/C_star` weight, so the zero result is load-bearing, not a vacuous derivative.
- No integrals / parity concerns in this unit.

**Verification command:**
The verifier will run `redteam exec-mathematica 234` and confirm the new `.wl` appears, exits 0, and independently reproduces M1-M7 (especially M4 cancellation and M7 balanced minimal-norm closed form).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_mathematica_audit.wl`
- summary: Created the Mathematica audit witness for M1-M7 using Series coefficients and Minimize for the balanced family.
- deviation: none

## F2 — tautological_check

**Target:** `scripts/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.py:146-163`

**Issue:** The Section 5 strip-endpoint assertions (lines 160-163) are algebraically guaranteed by construction. Each substitutes a strip edge `R1 = -c_eta*E1 ∓ const/Abs(vareps)` into `Xi1_twoobs = -R1 - c_eta*E1`; the `const` cancels identically, so the residual is zero for ANY value of `robust`/`nonempty`. The literals `0.367930328492646` (line 146) and `0.737619063660757` (line 147) are never actually exercised — a swap or mistranscription would still pass green. The constants are legitimate carried imports (do NOT derive them; do NOT change their values), but the check must fail if the values or their ordering are corrupted.

**Required change:**
Keep the existing strip prints (lines 165-168) and the `Xi1 = -R1 - c_eta*E1` rearrangement. Replace/augment the four tautological endpoint `assert_zero` calls at lines 160-163 with assertions that pin the ceiling values and their ordering. Specifically, after defining `robust` and `nonempty` (lines 146-147), add:

1. An ordering/positivity guard: assert `0 < float(robust) < float(nonempty)` (the nonempty corridor must be the strictly looser bound, and both strictly positive). Raise `AssertionError` if violated.
2. A half-width identity that still ties `robust`/`nonempty` to the strip but is NOT invariant to their value: e.g. assert the strip total width equals `2*robust/Abs(vareps)` for the robust gate and `2*nonempty/Abs(vareps)` for the nonempty gate by computing `simplify(robust_upper - robust_lower)` and comparing to `2*robust/sp.Abs(vareps)` (and likewise nonempty). This stays true to the strip definition while making a digit-drop in `robust`/`nonempty` detectable, since the literal now appears on both sides only through the SAME symbol (a corrupted value still matches itself) — therefore ALSO add (3).
3. A distinctness/sentinel guard so a swap is caught: assert `simplify(nonempty - robust) != 0` AND assert the documented gap, `abs(float(nonempty) - float(robust) - 0.369688735168111) < 1e-12` (this is `0.737619063660757 - 0.367930328492646`), so swapping the two literals or zeroing one fails.

Do not introduce a derivation for the constants. Do not touch the print statements' values. Edit only within the lines 146-163 region (you may add lines).

**Self-test confirmation (already done by auditor):**
- Guard (1) and (3) are genuinely value-sensitive: `0.367930328492646 < 0.737619063660757` holds, and `0.737619063660757 - 0.367930328492646 = 0.369688735168111`, so the assertions pass on the correct values and fail on swap/corruption.
- No new constant or paper claim is introduced (both ceilings are already stated verbatim in the notes and appendix), so no `paper_misalignment` is created by this fix.

**Verification command:**
The verifier will run `redteam exec-sympy 234` and confirm: the script exits 0 on the correct constants; the new ordering/gap assertions appear near lines 146-163; and (manual spot-check) swapping the two literals would now raise `AssertionError`.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.py`
- summary: Added positivity, ordering, distinctness, documented-gap, and strip half-width assertions for the two-observable ceiling constants.
- deviation: none
