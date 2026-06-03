---
unit_id: 232
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T21:53:56-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 232

Apply each non-`paper_misalignment` finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 is a `paper_misalignment` the USER HAS RESOLVED — see `## RESOLVED — F1` below: prefactor = **100** (the script is correct). Apply F1 by making ONLY the one authorized notes edit specified there (correct notes `168`→`100` on lines 153/157). Do NOT change the script. **F2 (the new `.wl`) now PROCEEDS** against the current SymPy literals (the 100-based `Pe_*`/`zeta_phys`/margin decimals already listed in the manifest below).

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

The ONLY authorized prose/notes edit is the F1 notes correction in `## RESOLVED — F1`; do not touch `paper.tex` or any other prose document.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper/notes side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md:153` quote: `\Xi_\chi = 168\,\Theta_w^{(\chi)}\Lambda_\ell^2 \approx 5.5548332017764099\times 10^5`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md:157` quote: `\Xi_J = 168\,\Theta_w^{(J)}\Lambda_\ell^2 \approx 1.2663707072528143\times 10^5`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.py:149` quote: `Xi_chi = 100 * Theta_w_chi * Lambda_ell**2`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.py:150` quote: `Xi_J = 100 * Theta_w_J * Lambda_ell**2`

## RESOLVED — F1 (user direction 2026-06-02: correct notes to script)

Direction (a): prefactor = **100** (script canonical). The notes' written `168` is a typo; only `100·Θ_w·Λ²` reproduces the notes' OWN quoted decimals (`Xi_chi ≈ 5.5548332017764099×10^5`, `Xi_J ≈ 1.2663707072528143×10^5`); `168·Θ_w·Λ²` gives `≈ 933212 / 212750`, contradicting them. No script change; no downstream regeneration. **Authorized notes edit (Codex APPLIES; Claude reviews):**

In `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md`, correct the prefactor `168` → `100` (leave the decimals — they already match 100):
- line 153: `\Xi_\chi = 168\,\Theta_w^{(\chi)}\Lambda_\ell^2` → `\Xi_\chi = 100\,\Theta_w^{(\chi)}\Lambda_\ell^2`
- line 157: `\Xi_J = 168\,\Theta_w^{(J)}\Lambda_\ell^2` → `\Xi_J = 100\,\Theta_w^{(J)}\Lambda_\ell^2`

The new F2 `.wl` uses `c = 100` (M4) and cross-engine-corroborates Xi_chi/Xi_J and the downstream Pe_*/margins.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md`
- summary: Corrected the two notes-side Xi prefactors from 168 to the resolved canonical value 100 while leaving the matching decimals unchanged.
- deviation: none

## F2 — missing_verification_script (subtype missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.wl`

**Issue:** Unit 232 has a SymPy audit and no Mathematica second engine (the stage card records "Mathematica audit: none yet"). The dual-engine rule requires a `.wl` wherever Mathematica CAN independently verify the stage; here it clearly can. The new script must be an INDEPENDENT route, not a transliteration of the SymPy `Delta_closed` algebra. In particular: obtain `Delta(Pe;kappa,eta)` by NATIVE symbolic/numeric integration of the kernel against the source profile (`Integrate`/`NIntegrate` of `K(x)·Sigma(x)` over `x∈[0,1]`), NOT by re-coding the hand-derived `Delta_closed` closed form. This gives a genuine cross-check of the SymPy script's analytic primitive. Use Mathematica's arbitrary-precision arithmetic; do not import the SymPy decimals as the source of truth — derive them.

**F1 RESOLVED (prefactor = 100, script unchanged):** the `.wl` targets the current SymPy literals — the 100-based `Pe_*`/`zeta_phys`/margin decimals already listed in the manifest below. Use prefactor `c = 100` in M4.

**Required change:**
Create the `.wl` at the Target path. It must independently derive and assert (with `If[Abs[...] > tol, Print["FAIL ..."]; Exit[1]]` style checks) the claims in the manifest below. Use Mathematica primitives end-to-end. Match the (post-F1-resolution) SymPy values to the tolerance each claim states. The integration route for M3/M5 must be `Integrate`/`NIntegrate` of the kernel·source product, NOT a port of `Delta_closed`.

**Claim manifest** (each must be independently verified in the new `.wl`):

- **M1 — Refreshed geometry.** With `x01 = BesselJZero[0,1]`, `Lambda_ell = 20 Sqrt[2] Pi / x01`, `chi_s = Lambda_ell/2`, `eta = Lambda_ell`, `kappa = (9/5) Lambda_ell^2`. Assert `Lambda_ell ≈ 36.94973154240256`, `kappa ≈ 2457.508789900114` (tol ≤ 5e-12).
- **M2 — Robin support ceiling.** Solve `y Tan[y] = eta` for the lowest root `y ∈ (0, Pi/2)` via `FindRoot`/`Reduce` (independent of bisection). Form `A_K = (kappa + Pi^2/4)/(kappa + y^2)` and `zeta_max = A_K Pi^2/4`. Assert `y ≈ 1.5294278190457656`, `A_K ≈ 1.0000521380385143`, `zeta_max ≈ 2.4675297457259358` (tol ≤ 5e-13).
- **M3 — Support-drop endpoints.** With `alpha = Sqrt[kappa]` and `denom = alpha Sinh[alpha] + eta Cosh[alpha]`, derive `Delta_0` and `Delta_inf`. Cross-check the endpoint formulas against the `Pe→0` and `Pe→∞` limits of the integral `Delta(Pe) = Integrate[K(x) Sigma_Pe(x), {x,0,1}]` (use the kernel `K(x) = (Cosh[alpha x] + (eta/alpha) Sinh[alpha x] - Cosh[alpha(1-x)])/denom` and source `Sigma_Pe(x) = Pe E^{Pe x}/(E^{Pe}-1)`, with the numerically stable rewrite for large Pe). Assert `Delta_0 ≈ 1.7377393923469950e-4`, `Delta_inf ≈ 2.0172162594593645e-2`.
- **M4 — Figure-of-merit.** `Xi_chi = c·Theta_w_chi·Lambda_ell^2`, `Xi_J = c·Theta_w_J·Lambda_ell^2` with `Theta_w_chi = 4.06863235008162`, `Theta_w_J = 0.927552032539308` and prefactor `c = 100` (F1 resolved). Assert `Xi_chi`, `Xi_J` decimals.
- **M5 — Fixed-point roots.** Solve `Pe = Xi · Delta(Pe)` via `FindRoot` (independent of the SymPy bisection), where `Delta(Pe)` is the NATIVE integral of M3, evaluated by Mathematica (`NIntegrate` or exact `Integrate` then evaluate). Assert the roots `Pe_*^(chi)`, `Pe_*^(J)` lie in `[Xi·Delta_0, Xi·Delta_inf]` and satisfy `Abs[Pe - Xi·Delta(Pe)] ≈ 0`. Assert their decimals (resolved values).
- **M6 — Physical support ratios.** `Omega_Pe = Pi Pe (2 Pe + Pi E^{-Pe})/((4 Pe^2 + Pi^2)(1 - E^{-Pe}))`; `zeta_phys = A_K Omega_Pe^2`; `rho_alpha,max = 1 + zeta_phys`. Assert `zeta_phys^(chi) ≈ 2.4675296478814376`, `zeta_phys^(J) ≈ 2.467527805167508` (resolved values; tol ≤ 5e-13).
- **M7 — Injected margins.** With `zeta_req = 1/3`, `rho_alpha^req = 4/3`: assert each `margin = zeta_phys - zeta_req > 0`, each `rho` margin `> 0`, each ceiling gap `zeta_max - zeta_phys > 0`, and the ratio values `zeta_phys/zeta_req ≈ 7.4026`, `rho/(4/3) ≈ 2.6006` (resolved values).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 232` (i.e. `math -script` on the new `.wl`) and confirms the script exits 0 with all in-file checks passing, AND that M3/M5 used a native integration route (not a `Delta_closed` transliteration), AND that the final `Delta`, `Pe_*`, `zeta_phys`, and margin values agree with the resolved SymPy values.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.wl`
- summary: Added the missing Mathematica audit using a native `Integrate`-derived support-drop function and independent `FindRoot` checks for the resolved Stage 232 constants, roots, support ratios, and margins.
- deviation: none
