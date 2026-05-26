---
unit_id: 026
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-25T23:36:02-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 026

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:197-202`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:157-162`

**Issue:**
The paper card (`/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_026.tex`) lists eq:app-stage026-Kreq as one of the four boxed `\stagefield{Output}` deliverables and the `\stagefield{Checks}` block explicitly says "Solving (eq:app-stage026-normalization) for K gives (eq:app-stage026-Kreq)." The paper's closed form is `K_req = kappa^2*lambda_B^2/varpi^2 + Q/Delta + mhat_rad^2 * kappa^2 * (Omega_U^2*lambda_W + lambda_R*lambda_U)^2 / (N_Q^target * Delta^2)`.

The script currently verifies only `residual @ K_req == 0` after solving for K_req via `sp.solve` / `Solve`. That back-substitution is algebraically guaranteed by the solver and cannot detect a divergence between the script's K_req and the paper's three-term structural form. A future change that altered the residual but happened to produce a Solve output back-substituting to zero would not be caught.

Add an explicit assertion that the script's K_req matches the paper's three-term decomposition.

**Required change:**

Edit `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`. In the `normalization_test` function, locate the existing block (currently lines 195-202):

```python
    # Solve exactly for the required wall stiffness and verify by back-substitution
    # that the residual vanishes when K is set to the solver's output.
    K_req = sp.solve(sp.Eq(residual, 0), K)[0]

    subbanner("IV.2 — Exact required wall stiffness")
    print("K_req =")
    sp.pprint(sp.simplify(K_req))
    expect_zero("residual @ K_req", residual.subs(K, K_req))
```

Replace it with:

```python
    # Solve exactly for the required wall stiffness and verify by back-substitution
    # that the residual vanishes when K is set to the solver's output.
    K_req = sp.solve(sp.Eq(residual, 0), K)[0]

    subbanner("IV.2 — Exact required wall stiffness")
    print("K_req =")
    sp.pprint(sp.simplify(K_req))

    # Independent structural check: the paper's eq:app-stage026-Kreq states the
    # three-term decomposition K_req = B0 + Q/Delta + mhat^2 * kappa^2 * (...)^2
    # / (target * Delta^2). Verify the solver's output matches that form.
    K_req_paper = (
        B0
        + Q / Delta
        + mhat**2 * kappa**2 * (Omega_U**2 * lambda_W + lambda_R * lambda_U)**2
          / (target * Delta**2)
    )
    expect_zero("K_req - K_req_paper", sp.simplify(K_req - K_req_paper))
    expect_zero("residual @ K_req", residual.subs(K, K_req))
```

All referenced symbols (`B0`, `Q`, `Delta`, `mhat`, `kappa`, `Omega_U`, `lambda_W`, `lambda_R`, `lambda_U`, `target`) are already in scope at this point in `normalization_test`. No imports change.

Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`. In `normalizationTest[]`, locate the existing block (currently lines 157-162):

```mathematica
  kReq = k /. First[Solve[residual == 0, k]];
  kReq = FullSimplify[kReq, Assumptions -> $Assumptions];

  subbanner["IV.2 — Exact required wall stiffness"];
  Print["K_req = ", fmt[kReq]];
  expectZero["residual @ K_req", residual /. k -> kReq];
```

Replace it with:

```mathematica
  kReq = k /. First[Solve[residual == 0, k]];
  kReq = FullSimplify[kReq, Assumptions -> $Assumptions];

  subbanner["IV.2 — Exact required wall stiffness"];
  Print["K_req = ", fmt[kReq]];

  (* Independent structural check: the paper's eq:app-stage026-Kreq states the
     three-term decomposition K_req = B0 + Q/Delta + mhat^2 * kappa^2 * (...)^2
     / (target * Delta^2). Verify the solver's output matches that form. *)
  kReqPaper = FullSimplify[
    b0
    + q/delta
    + mhat^2 * kappa^2 * (omegaU^2*lambdaW + lambdaR*lambdaU)^2
      / (target * delta^2),
    Assumptions -> $Assumptions
  ];
  expectZero["K_req - K_req_paper", kReq - kReqPaper];
  expectZero["residual @ K_req", residual /. k -> kReq];
```

All referenced symbols (`b0`, `q`, `delta`, `mhat`, `kappa`, `omegaU`, `lambdaW`, `lambdaR`, `lambdaU`, `target`) are already in scope inside `normalizationTest[]` (passed in via the tuple from `branchSubstitution[]` and the local `target` definition). The local Module variable list at line 148 already includes `kappa` (it should also include the new `kReqPaper` if Mathematica complains about an undeclared local — if so, add `kReqPaper` to the `Module[{...}]` declaration at line 148).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 026` and `redteam exec-mathematica 026`. Confirm:
- SymPy output contains a new line `K_req - K_req_paper = 0` (immediately before `residual @ K_req = 0` in section IV.2).
- Mathematica output contains a new `PASS: K_req - K_req_paper` line (immediately before `PASS: residual @ K_req`).
- Both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- summary: Added explicit checks that the solved K_req matches the paper's three-term structural decomposition before back-substitution.
- deviation: none
