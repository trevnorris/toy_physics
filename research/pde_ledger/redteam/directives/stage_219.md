---
unit_id: 219
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T13:58:53-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 219

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose document. The red-team only modifies scripts.

After editing, RUN the new Mathematica script (`math -script <path>`) and iterate until it exits 0 with every in-file check printing PASS. Use `timeout 600`. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.wl`

**Issue:** Stage 219 has only a SymPy audit; `stage_219.tex:11` records "Mathematica audit: none yet." The stage is non-status, non-checkpoint, and every deliverable is a finite closed-form symbolic computation (3x3 determinant, Schur complement, matrix inverse, scalar quadratic form, two rational identities, a product-kernel family decomposition, and a numeric positive-definiteness witness). Independent Mathematica verification is fully possible, so the dual-engine rule requires a second-engine `.wl`. Author a NEW Mathematica audit that re-derives the claims independently.

**Required change:**
Create the file at the Target path. Define the bundle symbolically (Mathematica symbols, NOT carried over from any Python intermediate). Verify the claim manifest M1..M7 below. End the script with a clear PASS/Exit[1] gate per claim so a failure is non-silent.

Recommended structure (Codex chooses exact code; the constraint is the manifest + the anti-transliteration guard):
- Build `Kred = {{Kstar,-GU,-GW},{-GU,OmegaU^2,-R},{-GW,-R,OmegaW^2}}`.
- Define `Delta = OmegaU^2 OmegaW^2 - R^2`, `Q = GU^2 OmegaW^2 + 2 GU GW R + GW^2 OmegaU^2`, `P = OmegaU^2 GW + R GU`, `PU = GU OmegaW^2 + R GW`, `D0 = Kstar - Q/Delta`.
- Use a single helper `zeroQ[e_] := TrueQ[Simplify[Together[e]] == 0]` and gate each claim with `If[!zeroQ[...], Print["FAIL Mn"]; Exit[1], Print["PASS Mn"]]`.

**Claim manifest** (each must be an independent Mathematica check; M-numbering matches the report):

- **M1** (determinant identity): `Det[Kred] - Delta*D0 == 0`. Native primitive: `Det`.
- **M2** (Schur complement): with internal block `Kint = {{OmegaU^2,-R},{-R,OmegaW^2}}` and coupling row `cpl = {-GU,-GW}`, the Schur complement `Kstar - cpl . Inverse[Kint] . cpl == D0`. Native primitive: `Inverse` on the 2x2 block (NOT a hand-substituted `Q/Delta`).
- **M3** (six inverse entries): with `Kinv = Inverse[Kred]`, verify
  `Kinv[[1,1]] == 1/D0`,
  `Kinv[[1,2]] == PU/(Delta*D0)`,
  `Kinv[[1,3]] == P/(Delta*D0)`,
  `Kinv[[2,2]] == (Kstar*OmegaW^2 - GW^2)/(Delta*D0)`,
  `Kinv[[2,3]] == (Kstar*R + GU*GW)/(Delta*D0)`,
  `Kinv[[3,3]] == (Kstar*OmegaU^2 - GU^2)/(Delta*D0)`.
  Native primitive: `Inverse` (the full 3x3, so this is a DIFFERENT decomposition from M2's block inverse).
- **M4** (collinear-source factorization): with `J = S*{sq,sU,sW}`, define `dV = -1/2 * (J . Kinv . J)`. Define
  `Ns = Delta*sq^2 + 2 PU sq sU + 2 P sq sW + (Kstar OmegaW^2 - GW^2) sU^2 + 2 (Kstar R + GU GW) sU sW + (Kstar OmegaU^2 - GU^2) sW^2`,
  `chiS = Ns/(Delta*D0)`. Verify `dV - (-1/2 chiS S^2) == 0`. Native primitive: dot products on `Inverse[Kred]`.
- **M5** (outgoing-load identities): with `Lambda = P/Delta`, `N0 = Lambda^2`, `P0 = N0/D0`, and `chiqW = Kinv[[1,3]]`, verify BOTH `chiqW - Lambda/D0 == 0` AND `chiqW^2 - P0/D0 == 0`.
- **M6** (primitive product-kernel family decomposition — the suppression theorem): set `SQ = x^(-3)`, `SY = Exp[-2 kappa x]/x`, `Jp = {betaQ*SQ, betaU*SY, betaW*SY}`. Define `dVp = -1/2 * (Jp . Kinv . Jp)`. Define `C6 = (1/D0) betaQ^2`, `C4 = (PU/(Delta*D0)) betaQ betaU + (P/(Delta*D0)) betaQ betaW`, `C2 = ((Kstar OmegaW^2 - GW^2)/(Delta*D0)) betaU^2 + 2 ((Kstar R + GU GW)/(Delta*D0)) betaU betaW + ((Kstar OmegaU^2 - GU^2)/(Delta*D0)) betaW^2`. Verify `dVp - (-1/2 (C6/x^6 + 2 C4 Exp[-2 kappa x]/x^4 + C2 Exp[-4 kappa x]/x^2)) == 0`.
  IMPORTANT — make the family-membership explicit (this is the load-bearing suppression claim, not just an equality): independently confirm there is NO `1/x` or `Exp[-2 kappa x]/x` term, e.g. multiply `dVp` by `x^6` and `Exp[4 kappa x]` and confirm the result is a polynomial in `x` of degree <= 4 with the exponential structure exactly `{1, Exp[2 kappa x], Exp[4 kappa x]}` — use `Exponent`/`CoefficientList` (in `x`) or `Collect[..., {x, Exp[...]}]` to demonstrate that the longest-range power is `x^(-2)` (paired with `Exp[-4 kappa x]`) and the slowest exponential decay is `x^(-4) Exp[-2 kappa x]`, with no surviving `x^(-1)` family. This must be a structural extraction, NOT a re-statement of M6's equality.
- **M7** (admissible positive-definiteness witness): substitute `Kstar->11, OmegaU->3, OmegaW->4, R->2, GU->1, GW->2`. Confirm `Delta == 140`, `D0 == 74/7`, `Det[Kred] == 1480`, all leading principal minors `> 0` (`11`, `98`, `1480`), and all `Eigenvalues[Kred]` are real and positive (equivalently `PositiveDefiniteMatrixQ`). Use a DIFFERENT route from the `.py` where possible (e.g. `PositiveDefiniteMatrixQ[Kred /. subs]` in addition to / instead of the explicit minor list).

**Anti-transliteration guard (MANDATORY):**
The `.wl` must be an INDEPENDENT re-derivation, not a line-by-line port of the `.py`. Concretely:
- Use Mathematica-native `Det`, `Inverse`, `Eigenvalues`/`PositiveDefiniteMatrixQ`, `Coefficient`/`CoefficientList`/`Exponent`/`Collect` — do not mirror the SymPy `together`/`simplify` choreography step-for-step.
- For M6, the suppression must be demonstrated by STRUCTURAL family extraction (`Collect`/`CoefficientList`/`Exponent`), a decomposition the SymPy script does NOT perform (the `.py` only asserts the equality at line 164). This different decomposition is the whole point of the second engine.
- Do not import, read, or reproduce any printed SymPy intermediate value as an input; derive every closed form inside Mathematica from the matrix.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 219` and confirms the new `.wl` exists, exits 0, and prints `PASS M1` .. `PASS M7`.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.wl`
- summary: Created the independent Mathematica audit for M1-M7, including native determinant/inverse checks, structural product-kernel family extraction, and the positive-definiteness witness.
- deviation: none
