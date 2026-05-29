---
unit_id: 005
batch: I.1
created_at: 2026-05-20T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T07:13:40Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 005

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl` (file does not exist; create it)

**Issue:** Unit 005 has `is_checkpoint: false` and `is_status_only_candidate: false` in `redteam/MANIFEST.yaml`, so both engines are required. No Mathematica audit script exists for this unit (subtype `missing_mathematica`). The SymPy script `scripts/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py` derives new projected-Maxwell identities and verifies them on a Gaussian kernel; without an independent Mathematica derivation, the SymPy PASS is internal consistency rather than two-engine verification.

**Required change:**

Create a new file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl` that independently verifies the claim manifest M1..M5 below. The script must:

1. Begin with a docstring block (Mathematica `(* ... *)`) stating that it independently verifies the projected inhomogeneous Maxwell law and projected charge continuity. Do NOT copy text from the SymPy docstring verbatim.

2. Define the projection kernel `Wg[w_] := Exp[-w^2]/Sqrt[Pi]` and its derivative `Wgp[w_] := D[Wg[w], w]`. Use `Integrate[Wg[w] expr, {w, -Infinity, Infinity}, Assumptions -> {t, x, y, z} \[Element] Reals]` for `Pg[...]` and analogously for `Pgp[...]`. Do **not** mirror the SymPy helpers' Python function names; use natural Mathematica naming such as `projW`, `projWPrime`, `boundaryValue`.

3. **Choose test profiles that differ from the SymPy choices**, so the two engines are exercising different algebra:
   - For the commutation check M1, use `Phi = (Sin[t] + x^2) (w^2 + 3)` (the SymPy script uses `(t^2 + x) (w^2 + 1)`).
   - For the IBP check M2, use `Q = w Exp[-w^2/4]` (a profile whose `Wg Q` boundary value decays differently than the SymPy polynomial `w^3 + w`; the IBP identity must still hold).
   - For the inhomogeneous-law check M3, use `G0 = Cos[t](w^2 + 2)`, `Gx = x^2 w^2`, `Gy = y(w^4 + 1)`, `Gz = z w^2`, `Gw = w + w^3/3`, `Gamma = (Sin[t] + x)(w^2 + 1)`.
   - For the leakage value M4, use `Gw = w` (the same as SymPy, since M4 asserts the universal Gaussian moment `-Pgp[w] == 1`; this overlap is intentional — it's the cross-check point).
   - For the continuity check M5, use `J0 = Cos[t](w^2 + 2)`, `Jx = x^2 w^2`, `Jy = y(w^4 + 1)`, `Jz = z w^2`, `Jw = w + w^3/3`.

4. Assert each item M1..M5 with the pattern
   ```
   If[Simplify[lhsExpr - rhsExpr] =!= 0, Print["FAIL: <label>"]; Exit[1], Print["PASS: <label>"]];
   ```
   and emit `Print["STATUS: PASS"]` as the last line.

5. Include at least one mutant guard per assertion, paralleling the SymPy `assert_nonzero` checks: e.g. assert that flipping the sign of the leakage term produces a nonzero residual, so a sign error in the derivation would be caught. Use:
   ```
   If[Simplify[<mutated_residual>] === 0, Print["MUTANT FAIL: <label>"]; Exit[1]];
   ```

**Claim manifest** (for missing-script findings only):

Reconstructed from the SymPy script's assertions A1-A14 (sympy script lines 226-303). Each `Proj_W[...]` denotes `Integrate[Wg[w] (...), {w, -Infinity, Infinity}]`; `Proj_{W'}[...]` denotes the same with `Wgp[w]`.

- **M1 (commutation with brane derivatives).** For any smooth `Phi(t, x, y, z, w)` with sufficient decay in `w`,
  `Proj_W[partial_t Phi] - partial_t Proj_W[Phi] == 0`,
  and the same with `partial_t` replaced by `partial_x, partial_y, partial_z`. Mathematica must check at least one of these (the script in SymPy checks `partial_t`).

- **M2 (projected IBP identity, boundary-decay form).** For any `Q(w)` with `Wg Q -> 0` at `w -> +/- infinity`,
  `Proj_W[partial_w Q] + Proj_{W'}[Q] == 0`.
  Additionally verify `boundaryValue[Wg Q] == 0` for the chosen `Q`, and that flipping the sign in the IBP residual (i.e. `Proj_W[partial_w Q] - Proj_{W'}[Q]`) is nonzero.

- **M3 (projected inhomogeneous Maxwell law, boundary-decay form).** For test profiles `G^a` (`a in {0, x, y, z}`), `G^w`, `Gamma` and source `J_nu`,
  ```
  Proj_W[ partial_t G^0 + partial_x G^x + partial_y G^y + partial_z G^z + partial_w G^w + Gamma ]
    ==
  partial_t Proj_W[G^0] + partial_x Proj_W[G^x] + partial_y Proj_W[G^y] + partial_z Proj_W[G^z]
    - Proj_{W'}[G^w] + Proj_W[Gamma]
  ```
  under the assumption `boundaryValue[Wg G^w] == 0` (verify this separately for the chosen `G^w`). Also assert that flipping the sign of the `-Proj_{W'}[G^w]` term produces a nonzero residual.

- **M4 (analytic transverse-leakage value, Gaussian kernel).** With `Wg[w] = Exp[-w^2]/Sqrt[Pi]`,
  `-Proj_{W'}[w] == 1` exactly.
  Also assert that `-Proj_{W'}[w]` is nonzero (i.e. the leakage does not vanish trivially under the kernel choice). This is the engine cross-check point with SymPy line 278.

- **M5 (projected charge continuity, open-system form).** For test current components `J^0, J^x, J^y, J^z, J^w` with `boundaryValue[Wg J^w] == 0` (verify separately),
  ```
  Proj_W[ partial_t J^0 + partial_x J^x + partial_y J^y + partial_z J^z + partial_w J^w ]
    ==
  partial_t Proj_W[J^0] + partial_x Proj_W[J^x] + partial_y Proj_W[J^y] + partial_z Proj_W[J^z]
    - Proj_{W'}[J^w]
  ```
  i.e. `partial_t Proj_W[J^0] + partial_a Proj_W[J^a] == Proj_{W'}[J^w]`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 005` and confirm:
- the file `mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl` exists,
- it exits 0,
- the saved output `mathematica/output/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.txt` contains a `STATUS: PASS` line and `PASS:` lines for each of M1..M5 plus their mutant guards,
- the leakage cross-check (M4) reports `1` (matching SymPy line 278).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl`
- summary: Created the independent Mathematica projected-Maxwell audit covering M1 through M5 with Gaussian projection, boundary checks, leakage reporting, and mutant guards.
- deviation: none
