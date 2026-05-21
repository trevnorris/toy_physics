---
unit_id: 007
batch: I.1
created_at: 2026-05-20T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T11:02:18-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 007

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl` (new file)

**Issue:**
Unit 007 currently has only a SymPy verification script. Per the second-engine policy (unit is not status-only and not a checkpoint carve-out), an independent Mathematica derivation is required. The Mathematica script must not transliterate the SymPy algebra — it must perform `Integrate` calls natively, name its own intermediates, and arrive at the same closed forms by independent route.

**Required change:**
Create the file `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl`. It must:

1. Declare symbols with explicit positivity assumptions: `lambda > 0`, `mu0 > 0`, `sigma > 0`, `tau > 0`, `epsilon > 0`. Use Mathematica's native `$Assumptions` or `Assuming[...]` rather than copying SymPy's `symbols(..., positive=True)` style.
2. Define `Z[w_] := Exp[-w^2/lambda^2]` and compute `Zint = Integrate[Z[w], {w, -Infinity, Infinity}]` and `Z2int = Integrate[Z[w]^2, {w, -Infinity, Infinity}]` directly. Verify symbolically with `FullSimplify[Zint - Sqrt[Pi]*lambda] === 0` and `FullSimplify[Z2int - Sqrt[Pi/2]*lambda] === 0`. On failure, `Print` the residual and `Exit[1]`.
3. Define smooth normalized Gaussians `Wsmooth[w_] := Exp[-w^2/sigma^2]/(Sqrt[Pi]*sigma)` and `Ssmooth[w_] := Exp[-w^2/tau^2]/(Sqrt[Pi]*tau)`. Compute `IWZsmooth = Integrate[Wsmooth[w] Z[w], {w, -Infinity, Infinity}]` and `IWSsmooth = Integrate[Wsmooth[w] Ssmooth[w], {w, -Infinity, Infinity}]`. Assert `FullSimplify[IWZsmooth - lambda/Sqrt[lambda^2 + sigma^2]] === 0` and `FullSimplify[IWSsmooth - 1/(Sqrt[Pi]*Sqrt[sigma^2 + tau^2])] === 0`. Each failure -> `Exit[1]`.
4. Verify the field-mutation Gaussian-moment identity: with `etaSym` a real nonzero parameter and `Fmut[x_, w_] := Sin[k x] + x^2 + etaSym*x*w^2`, define `fieldMutationLHS = Integrate[Wsmooth[w] Z[w] D[Fmut[x, w], x], {w, -Infinity, Infinity}]` and `fieldMutationDelta = FullSimplify[fieldMutationLHS - IWZsmooth*D[Sin[k x] + x^2, x]]`. Assert `FullSimplify[fieldMutationDelta - etaSym*lambda^3*sigma^2/(2*(lambda^2 + sigma^2)^(3/2))] === 0`. Failure -> `Exit[1]`.
5. Verify the source-mutation Gaussian-moment identity: with `Jmut[x_, w_] := Ssmooth[w]*(Cos[k x] + x) + etaSym*x*w^2*Ssmooth[w]`, define `sourceMutationRHS = Integrate[Wsmooth[w] mu0 Jmut[x, w], {w, -Infinity, Infinity}]` and `sourceMutationDelta = FullSimplify[sourceMutationRHS - mu0*IWSsmooth*(Cos[k x] + x)]`. Assert `FullSimplify[sourceMutationDelta - etaSym*mu0*x*sigma^2*tau^2/(2*Sqrt[Pi]*(sigma^2 + tau^2)^(3/2))] === 0`. Failure -> `Exit[1]`.
6. Matched-observer case. Define `Wmatch[w_] := Z[w]/Zint`. Compute `IWZmatch = Integrate[Wmatch[w] Z[w], {w, -Infinity, Infinity}]`. Assert `FullSimplify[IWZmatch - 1/Sqrt[2]] === 0` (this is the substantive check; do not just restate `Z2int/Zint`). For the delta-source side, evaluate `IWSmatch = Wmatch[0]` to install the delta-sampling identity explicitly, then `mu0ProjMatch = mu0*IWSmatch/IWZmatch` and `mu0Red = mu0/Zint`. Assert `FullSimplify[mu0ProjMatch/mu0Red - Sqrt[2]] === 0`. Failure -> `Exit[1]`.
7. Regulator case. Define `Weps[w_] := Exp[-w^2/epsilon^2]/(Sqrt[Pi]*epsilon)`. Compute `IWZeps = Integrate[Weps[w] Z[w], {w, -Infinity, Infinity}]` and `IWSeps = Integrate[Weps[w]^2, {w, -Infinity, Infinity}]`. Assert `FullSimplify[IWZeps - lambda/Sqrt[epsilon^2 + lambda^2]] === 0` and `FullSimplify[IWSeps - Sqrt[2]/(2*Sqrt[Pi]*epsilon)] === 0`. Then verify `Limit[IWZeps, epsilon -> 0, Direction -> "FromAbove"] === 1`. Failure -> `Exit[1]`.
8. After all checks pass, `Print["STATUS: PASS"]` and `Exit[0]`.

Use the `wl` style consistent with the repository's other `mathematica_audit.wl` files. Do not import or read the SymPy file's variable names directly (no `W_smooth`, `df_test`, `I_WZ_smooth` — use the Mathematica-native names above).

**Claim manifest:**
The new Mathematica script must independently verify the following physical/algebraic results:

- M1: `int_{-infty}^{+infty} exp(-w^2/lambda^2) dw = sqrt(pi)*lambda` (i.e. `Z_int = sqrt(pi)*lambda`).
- M2: `int_{-infty}^{+infty} exp(-2 w^2/lambda^2) dw = sqrt(pi/2)*lambda` (i.e. `Z2_int = sqrt(pi/2)*lambda`).
- M3: Smooth Gaussian overlap `int W(w) Z(w) dw = lambda/sqrt(lambda^2 + sigma^2)` for `W(w) = exp(-w^2/sigma^2)/(sqrt(pi)*sigma)`.
- M4: Smooth Gaussian source overlap `int W(w) S(w) dw = 1/(sqrt(pi)*sqrt(sigma^2 + tau^2))` for `S(w) = exp(-w^2/tau^2)/(sqrt(pi)*tau)`.
- M5: Field-mutation moment: `int W(w) Z(w) w^2 dw = sigma^2*lambda^3/(2*(sigma^2 + lambda^2)^{3/2})`, equivalently the delta from a `eta*x*w^2` admixture equals `eta*lambda^3*sigma^2/(2*(lambda^2 + sigma^2)^{3/2})`.
- M6: Source-mutation moment: `int W(w) S(w) w^2 dw = sigma^2*tau^2/(2*sqrt(pi)*(sigma^2 + tau^2)^{3/2})`, equivalently the source delta is `eta*mu0*x*sigma^2*tau^2/(2*sqrt(pi)*(sigma^2 + tau^2)^{3/2})`.
- M7: Matched-observer overlap `I_WZ_match := int (Z/Z_int) Z dw = 1/sqrt(2)`.
- M8: Matched-observer / delta-source projection-vs-reduction ratio: `mu0_proj_match / mu0_red = sqrt(2)`, where `mu0_proj_match := mu0 * (Z(0)/Z_int) / I_WZ_match` and `mu0_red := mu0 / Z_int`.
- M9: Regulator overlap `int W_eps(w) Z(w) dw = lambda/sqrt(epsilon^2 + lambda^2)`, with `W_eps(w) = exp(-w^2/epsilon^2)/(sqrt(pi)*epsilon)`.
- M10: Regulator self-overlap `int W_eps(w)^2 dw = sqrt(2)/(2*sqrt(pi)*epsilon)`, demonstrating divergence as `epsilon -> 0`.
- M11: Sharp-sampling limit `lim_{epsilon -> 0+} I_WZ_eps = 1`.

Every claim must be asserted with an explicit `FullSimplify[... ] === 0` (or `=== 1` / `=== Sqrt[2]` as appropriate) and a guard that calls `Exit[1]` on failure.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 007` and confirm the new script exists at the path above, that all eleven claims M1-M11 are explicitly asserted, and that the script exits 0 with output `STATUS: PASS`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl`
- summary: Created the independent Mathematica audit script with native Gaussian integrations and explicit M1-M11 residual checks.
- deviation: none
