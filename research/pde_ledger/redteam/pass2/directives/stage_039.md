---
unit_id: 039
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-05T13:58:24Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 039

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:114-117`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:95-98`

**Issue:** The `z1/z0`-ratio assertion is the only check intended to exercise the boxed direction factor `R_U` (paper `\stagefield{Output}`, eq:app-stage039-RU). It is tautological: `z0 = kappa0*g_W*(1+rho0)` and `z1 = kappa1*g_W*(1+rho0/(1+deltaU))` are defined with `kappa0`/`kappa1` prefactors, so `z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))` is identically 0 by construction and never references `R_U`. A wrong `R_U` closed form would not be caught.

**Required change:**
Replace the tautological assertion with a falsifiable one that ties the loading-vector component ratio to the independently-defined `R_U`.

SymPy — replace the `expect_zero(...)` block at lines 114–117:
```python
expect_zero(
    "z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))",
    sp.simplify(z1 * (1 + rho0) - (kappa1 / kappa0) * z0 * (1 + rho0 / (1 + deltaU))),
)
```
with:
```python
expect_zero(
    "z1/z0 - (kappa1/kappa0)*R_U",
    sp.simplify(z1 / z0 - (kappa1 / kappa0) * R_U),
)
```
`R_U` is already defined at line 109 (`R_U = sp.simplify((1 + rho0 / (1 + deltaU)) / (1 + rho0))`); keep it and its print (line 113) unchanged.

Mathematica — replace the `expectZero[...]` block at lines 95–98:
```wolfram
expectZero[
  "z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))",
  z1*(1 + rho0) - (kappa1/kappa0)*z0*(1 + rho0/(1 + deltaU))
];
```
with:
```wolfram
expectZero[
  "z1/z0 - (kappa1/kappa0)*R_U",
  z1/z0 - (kappa1/kappa0)*rU
];
```
`rU` is already defined at line 90; keep it and its print (line 94) unchanged.

**Self-test (pre-verified):** substituting the definitions, `z1/z0 = (kappa1/kappa0)*(1+rho0/(1+deltaU))/(1+rho0) = (kappa1/kappa0)*R_U`, so the residual reduces to 0 only because `R_U` equals the boxed form; a deliberately-wrong `R_U` (e.g. dropping the `(1+deltaU)`) makes it nonzero. Non-tautological.

**Verification command:** after Codex applies, the verifier runs `redteam exec-sympy 039` and `exec-mathematica 039`; a check named `z1/z0 - (kappa1/kappa0)*R_U` must appear, evaluate to `0`/`PASS`, and the scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Replaced the tautological loading-vector ratio residual with a direct check against the independently defined `R_U`/`rU`.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:147-151`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:127-131`

**Issue:** The surviving exact product law is a stated deliverable (notes §5: `R_target^(splitU) M_mix^(splitU) = 8 Lambda (1 - eps_W_split)/pi^2`, the headline "factorization survives" claim). The scripts compute `product = M_mix_split * R_target_split` and print it but never assert it equals the closed form, so a regression passes silently.

**Required change:**
Add a falsifiable assertion right after the existing `product` print, using the already-defined `eps_W_split`/`epsWSplit`.

SymPy — after line 151 (`print("product =", product)`), add:
```python
expect_zero(
    "product law survives",
    product - 8 * Lambda * (1 - eps_W_split) / pi**2,
)
```
(`product` is defined at line 147, `eps_W_split` at line 94, `Lambda`/`pi` already in scope.)

Mathematica — after line 131 (`Print["product = ", fmt[product]];`), add:
```wolfram
expectZero["product law survives", product - 8 lambda (1 - epsWSplit)/Pi^2];
```
(`product` is defined at line 127, `epsWSplit` at line 79, `lambda` is the script's symbol for `Lambda`.)

**Self-test (pre-verified):** at `deltaU=0` both sides give `8 Lambda(1-eps_W)/pi^2` (residual 0); for `deltaU≠0` the script's printed `product` (sympy output l.45) equals `8 Lambda(11 deltaU - eps_W(9 deltaU+11)+11)/(11 pi^2(deltaU+1))`, which is exactly `8 Lambda(1-eps_W_split)/pi^2` with `eps_W_split = eps_W(9 deltaU+11)/(11(1+deltaU))`. So the assert_zero is satisfiable and non-trivial: `product` is built from independently-`simplify`d `M_mix_split`/`R_target_split`, the RHS from the notes' closed form. No new constant introduced (no paper round-trip violation).

**Verification command:** the verifier runs `redteam exec-sympy 039` and `exec-mathematica 039`; a check `product law survives` must appear, evaluate to `0`/`PASS`, scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Added the surviving exact product-law residual checks immediately after the existing product prints.
- deviation: none

## F3 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`

**Issue:** the SymPy source carries stale pre-renumber SELF-labels (pre-renumber self-number `22`; canonical = `039`) in the docstring header + the `22.k` subbanner indices, matching the file's canonical banner `STAGE 39`. The main banner was already fixed; these self-labels were missed and print INTO the transcript, so a plain re-run leaves them stale. (The committed `.txt` are also stale — refreshed by the F1/F2 re-run. The `.wl` source labels are already canonical.)

**Required change (label-only — change ONLY the indicated numeric token, preserve the rest of each string verbatim):**
- line 3: `Moving-throat PDE — Stage 22 SymPy audit.` → `Moving-throat PDE — Stage 39 SymPy audit.`
- line 68: subbanner `22.1` → `39.1`
- line 89: subbanner `22.2` → `39.2`
- line 103: subbanner `22.3` → `39.3`
- line 140: subbanner `22.4` → `39.4`
- line 153: subbanner `22.5` → `39.5`

**DO NOT TOUCH (deferred CROSS-refs — dedicated pass owns these):**
- lines 15, 142, 178 (`Stage-21 continuum placement map` / `Stage-21 factorization` cross-refs to upstream stage 038) — LEAVE.
- the already-canonical `STAGE 39` banners (lines 62, 169) — LEAVE.

**Verification command:** the `.py` docstring + subbanners read `39` / `39.1`–`39.5`; after the re-run both transcript banners read `STAGE 39`/`STAGE 039`, all PASS lines remain, and the strip-the-number diff is byte-identical.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage039_split_u_sector_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage039_split_u_sector_mathematica_audit.txt`
- summary: Updated the SymPy self-label numeric tokens to Stage 39 and regenerated the saved audit transcripts.
- deviation: none
