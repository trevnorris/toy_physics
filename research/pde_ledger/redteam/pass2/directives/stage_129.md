---
unit_id: 129
batch: IV.4
created_at: 2026-06-06T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-06T17:02:21-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 129

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl:33`

**Issue:** The `.wl` hard-codes the same closed-form profile `sigma = piM*Exp[-piM*z/lM]/(lM*(1 - Exp[-piM]))` (l.33) that the SymPy script defines (py l.13), then runs the identical normalization/zero-flux/residual checks with the identical `v1 -> piM*thetaSigma/lM` substitution. It is a line-by-line port, not an independent derivation. The notes (l.98-119) DERIVE the profile by solving the zero-flux ODE `Θσ' + V₁σ = 0` and imposing normalization; the second engine should reproduce that derivation rather than copy the answer.

**Required change:**
Make the `.wl` independently derive the profile, then confirm it equals the boxed form. Concretely, BEFORE the existing checks (i.e., add after l.34 where `jSigma` is defined, keeping `sigma` as currently defined so the downstream checks are unchanged), insert an independent derivation block:

```
(* Independent re-derivation: solve the stationary zero-flux ODE and normalize *)
sol = DSolve[thetaSigma*sigmaFn'[z] + v1*sigmaFn[z] == 0, sigmaFn, z];
sigmaGen = sigmaFn[z] /. First[sol];                 (* = C[1]*Exp[-v1 z/thetaSigma] *)
cVal = C[1] /. First[Solve[Integrate[sigmaGen, {z, 0, lM}] == 1, C[1]]];
sigmaDerived = FullSimplify[(sigmaGen /. C[1] -> cVal) /. v1 -> piM*thetaSigma/lM,
  Assumptions -> $Assumptions];
Print["Independently derived sigma = ", fmt[sigmaDerived]];
expectZero["derived profile matches boxed sigma_Pi", sigmaDerived - sigma];
```

Place this block after the `jSigma = ...` definition (current l.34) and before the `Print["sigma_Pi(z) = "...]` line (current l.36), so the existing prints and the existing three checks remain intact and now serve as confirmations of an independently-derived result. Do NOT delete or alter the existing `sigma` definition (l.33) or the three existing `expectZero` checks.

**Self-test (verified by auditor):**
- `DSolve` variable `z` is a genuine dependency of `sigmaFn` → derivative non-trivial; solution is `C[1]*Exp[-v1 z/thetaSigma]`.
- `Integrate[C[1]*Exp[-v1 z/thetaSigma], {z,0,lM}] = C[1]*thetaSigma*(1 - Exp[-v1 lM/thetaSigma])/v1`; setting `=1` and substituting `v1 -> piM thetaSigma/lM` gives `C[1] = piM/(lM(1 - Exp[-piM]))`, so `sigmaDerived = piM Exp[-piM z/lM]/(lM(1-Exp[-piM]))` = `sigma` exactly → `sigmaDerived - sigma` reduces to 0. Non-tautological: a wrong ODE coefficient or wrong normalization would make the residual nonzero.

**Verification command:**
After Codex applies, the verifier runs `math -script` on the `.wl` and confirms: a `DSolve[...]` call and a new `expectZero["derived profile matches boxed sigma_Pi", ...]` appear, the new check PASSes (residual 0), and the script exits 0. The SymPy script is unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl`
- summary: Added an independent Mathematica DSolve-and-normalize derivation and checked it against the boxed sigma profile.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py:22`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl:34`

**Issue:** The chemical-potential → Onsager-current link (notes l.64-89), which is the physical motivation of the whole stage, is never asserted. The SymPy script computes `mu = Theta*sp.log(sigma/sigma_star) + V1*z` (l.22) but never uses it; the `.wl` omits `μ` entirely. The identity `−M σ ∂_z μ = −M(Θσ' + V₁σ)` (deliverable D4) must be checked in both engines.

**Required change:**

SymPy (`..._sympy_audit.py`): the `mu` definition already exists at l.22. After l.24 (after the `res` computation and its print), add:
```
J_from_mu = -M*sigma*sp.diff(mu, z)
mu_link_res = sp.simplify(J_from_mu - J)
print("Onsager current from mu identity residual =", mu_link_res)
if mu_link_res != 0:
    raise AssertionError("Onsager current does not match -M*sigma*d(mu)/dz.")
```
(`J` is defined at l.14 as `-M*(Theta*sp.diff(sigma, z) + V1*sigma)`; this uses the existing `mu` and `J`, no new symbols.)

Mathematica (`..._mathematica_audit.wl`): after the `jSigma` definition (current l.34), add the `μ` definition and after the existing residual check (current l.46) add the identity check:
```
mu = thetaSigma*Log[sigma/sigmaStar] + v1*z;
```
(place near l.34) and
```
expectZero["Onsager current from mu identity", (-mobility*sigma*D[mu, z]) - jSigma];
```
(place after the existing `expectZero["boundary-layer ODE residual", residual];` at l.46). Use `jSigma` as defined at l.34, NOT `jSub` (the identity must hold before the `v1` substitution).

**Self-test (verified by auditor):**
- `mu` depends on `z` via both `Log[sigma/...]` and the `v1 z` term → `D[mu, z]` non-trivial.
- `∂_z μ = Θ σ'/σ + V₁`, so `−M σ ∂_z μ = −M(Θσ' + V₁σ) = J`; residual reduces to exactly 0. Non-tautological: a wrong coefficient on the entropic (`Θ`) or potential (`V₁`) term in `μ` would make it nonzero.
- This identity holds for ALL `v1` (it does not require the `v1 -> piM thetaSigma/lM` substitution), which is why the `.wl` check uses `jSigma` not `jSub`.

**Verification command:**
After Codex applies, the verifier runs `python3` on the `.py` and `math -script` on the `.wl`, confirms the new μ→J residual check appears and PASSes (residual 0) in BOTH engines, and both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl`
- summary: Added SymPy and Mathematica checks that the chemical-potential gradient reproduces the Onsager current.
- deviation: none
