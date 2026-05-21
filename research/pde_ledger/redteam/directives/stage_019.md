---
unit_id: 019
batch: I.2
created_at: 2026-05-21T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T19:35:28Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 019

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl`

**Issue:** The unit has no Mathematica auditor. The full `(P0, P2, P4)` bundle bridge — closed-form `KSigma` on the one-pole and normalization branches, `N2_const`/`N4_const` constant-prefactor closures, Jacobian determinant `D0^3`, M-root factorization, and the concrete Gaussian wall-integral example — rests on a single algebra engine. Add an independent Mathematica auditor that derives each claim from the bundle definition rather than transliterating the SymPy variable choreography.

**Required change:**
Create the file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl` containing an independent Mathematica derivation of M1–M12 below. The script MUST:

1. Define symbolic variables: `KSigma, MSigma, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, mhat0, G, cs, a, c, eps`.
2. Define `D0 = KSigma - B0 - Z0`, `D2 = -(MSigma + B2 + Z2)`, `D4 = -(B4 + Z4)`.
3. Define `u2`, `u4` from the small-x expansion of `1/(D0 + D2*x + D4*x^2)` truncated at quadratic order (i.e., derive them, do not copy from sympy). Mathematica-native way: `Series[1/(D0 + D2*x + D4*x^2), {x, 0, 4}]` and read off the coefficients with `Coefficient[Normal[%], x, k]/Coefficient[Normal[1/D0], x, 0]^(-1)` — or compute `1 + (D2/D0)*x + ...` by hand from `1/(1 + (D2 x + D4 x^2)/D0)` using `Series`.
4. Define `P0 = N0/D0`, and derive `P2, P4` by expanding `Series[(N0 + N2 x^2 + N4 x^4)/(D0 + D2 x + D4 x^2), {x, 0, 4}]` and reading off coefficients of `x^0, x^2, x^4`. Do NOT hand-copy the SymPy form `(D0 N2 - 2 D2 N0)/D0^2`; let `Series` produce it.
5. Use `FullSimplify` (not `Simplify`) on each residual. On any mismatch, print `Print["FAIL: <label> residual = ", residual]; Exit[1]`. On full pass, print `Print["STATUS: PASS"]` and `Exit[0]`.
6. At the end print a one-line per claim summary (`Print["M1 OK"]`, etc.) so the saved output shows the checks executed.

**Claim manifest** (for missing-script findings only):

- **M1**: `u4 - 4 u2^2 == (D0 (B4 + Z4) - 3 (MSigma + B2 + Z2)^2)/D0^2`. Verify with `FullSimplify[(u4 - 4 u2^2) - (D0 (B4 + Z4) - 3 (MSigma + B2 + Z2)^2)/D0^2] == 0`.

- **M2** (one-pole `KSigma`): `Solve[u4 - 4 u2^2 == 0, KSigma]` returns exactly `KSigma -> B0 + Z0 + 3 (MSigma + B2 + Z2)^2/(B4 + Z4)`. Verify by `FullSimplify[KSigma /. First[Solve[u4 - 4 u2^2 == 0, KSigma]] - (B0 + Z0 + 3 (MSigma + B2 + Z2)^2/(B4 + Z4))] == 0`.

- **M3** (normalization `KSigma`): With `P0target = 54 G cs^5/(5 a^5 c^5 mhat0^2)`, `Solve[P0 == P0target, KSigma]` returns `KSigma -> B0 + Z0 + N0/P0target`. Verify analogously.

- **M4** (`N2_const`): `Solve[P2 == 0, N2]` returns `N2 -> 2 N0 (B2 + MSigma + Z2)/(B0 - KSigma + Z0)`. Verify `FullSimplify[(N2 /. First[Solve[P2 == 0, N2]]) - 2 N0 (B2 + MSigma + Z2)/(B0 - KSigma + Z0)] == 0`.

- **M5** (`N4_const`): substituting the M4 solution into `P4` and solving `P4 == 0` for `N4` yields
  `N4 -> N0 (2 B0 B4 + 2 B0 Z4 + B2^2 + 2 B2 MSigma + 2 B2 Z2 - 2 B4 KSigma + 2 B4 Z0 - 2 KSigma Z4 + MSigma^2 + 2 MSigma Z2 + 2 Z0 Z4 + Z2^2)/(B0 - KSigma + Z0)^2`. Verify with `FullSimplify`.

- **M6** (Jacobian determinant): with `P2zeroEq = Expand[D0^2 P2]` and `P4zeroEq = Expand[D0^3 P4]`, the matrix `{{D[P2zeroEq, N2], D[P2zeroEq, N4]}, {D[P4zeroEq, N2], D[P4zeroEq, N4]}}` has determinant `D0^3`. Verify `FullSimplify[Det[...] - D0^3] == 0`. NOTE: `P2zeroEq` does not depend on `N4`, so `D[P2zeroEq, N4]` is identically 0 by construction — the determinant simplifies to `D[P2zeroEq, N2] * D[P4zeroEq, N4]`. This is correct, not a bug; just verify the resulting product equals `D0^3`.

- **M7** (P2 factorization): `FullSimplify[P2 - (N2 - N2closed)/D0] == 0` where `N2closed = 2 N0 (B2 + MSigma + Z2)/(B0 - KSigma + Z0)`.

- **M8** (P4 factorization on N2 = N2closed): `FullSimplify[(P4 /. N2 -> N2closed) - (N4 - N4closed)/D0] == 0`.

- **M9** (mutation guards): with `eps` declared nonzero, `FullSimplify[P2 - (N2 - (N2closed + eps))/D0]` and `FullSimplify[(P4 /. N2 -> N2closed) - (N4 - (N4closed + eps))/D0]` MUST be nonzero (both equal `-eps/D0`). Test with `FullSimplify[residual] === 0` returning False (or use `If[FullSimplify[residual] === 0, Print["FAIL: mutation guard"]; Exit[1]]`).

- **M10** (M-root factorization): the polynomial `D0 (B4 + Z4) - 3 (MSigma + B2 + Z2)^2` has roots `MSigma -> Sqrt[D0 (B4 + Z4)/3] - (B2 + Z2)` and `MSigma -> -Sqrt[D0 (B4 + Z4)/3] - (B2 + Z2)`. Verify the factorization `D0 (B4 + Z4) - 3 (MSigma + B2 + Z2)^2 == -3 (MSigma - Mplus)(MSigma - Mminus)` via `FullSimplify[lhs - rhs] == 0` (use `PowerExpand` if needed to handle the sqrt). Verify Vieta: `Mplus + Mminus == -2 (B2 + Z2)` and `Mplus Mminus == (B2 + Z2)^2 - D0 (B4 + Z4)/3`.

- **M11** (u2 on M-roots): `FullSimplify[(u2 /. MSigma -> Mplus) - Sqrt[D0 (B4 + Z4)/3]/D0]` and `FullSimplify[(u2 /. MSigma -> Mminus) + Sqrt[D0 (B4 + Z4)/3]/D0]` both equal 0. Use `Assumptions -> {D0 > 0, B4 + Z4 > 0}` or `PowerExpand` so the sqrt doesn't pick the wrong branch.

- **M12** (concrete Gaussian wall integrals): with `beta[w_] := Exp[-w^2/2]`, `muEta = 1`, `Tw = 1`, `Tomega = 1/6`, `Keta = 0`:
  - `MSigmaExample = Integrate[muEta beta[w]^2, {w, -Infinity, Infinity}]` equals `Sqrt[Pi]`.
  - `KSigmaExample = Integrate[Tw D[beta[w], w]^2 + (Keta + 6 Tomega) beta[w]^2, {w, -Infinity, Infinity}]` equals `3 Sqrt[Pi]/2`.

  Verify `MSigmaExample - Sqrt[Pi] == 0` and `KSigmaExample - 3 Sqrt[Pi]/2 == 0`.

After defining `P0target`, `Mplus`, `Mminus`, `N2closed`, `N4closed` as standalone Mathematica quantities (not copied from SymPy output text), run the assertions M1–M12 in order and abort with `Exit[1]` on any failure.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 019` and confirm the new `.wl` produces saved output containing the literal strings `M1 OK`, `M2 OK`, …, `M12 OK`, `STATUS: PASS`, and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl`
- summary: Added the independent Mathematica auditor for M1-M12 covering the isotropic bundle bridge, constant-prefactor closures, M-root checks, and Gaussian wall integrals.
- deviation: none
