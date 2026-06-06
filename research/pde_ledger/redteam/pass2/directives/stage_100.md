---
unit_id: 100
batch: IV.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T23:39:08-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 100

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond the named work. Do NOT touch paper.tex, notes/, or any prose documents — the red-team only modifies scripts.

Do NOT modify the SymPy script `scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py` — it is the reference engine and stays as-is. ONLY the Mathematica script is rewritten.

After editing, RUN the script (`math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration

**Subtype:** independent-route required (dual-engine policy).

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl` (whole derivation body, lines ~40–72)

**Issue:**
The `.wl` is a line-by-line transliteration of the SymPy `.py`, not an independent re-derivation. Corresponding sections:

- `.py:32-34` builds `sigma_can = (9/8)/Omega**5`, `Y = 3/4 + (1/4)/(1 - omega**2/Omega**2 - I*chiQ*sigma_can*omega**5)`, then `Yser = expand(series(Y, omega, 0, 6).removeO())`.
- `.wl:40-42` builds the SAME `sigmaCan = (9/8)/omegaQ^5`, the SAME `yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5)`, then `ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]]`.
- `.py:36-38` extracts `K2 = K0*Yser.coeff(omega,2)`, `K4 = K0*Yser.coeff(omega,4)`, `Gamma5 = im(Yser.coeff(omega,5))*K0`.
- `.wl:44-46` extracts the SAME `k2 = k0*Coefficient[ySeries, omega, 2]`, `k4 = k0*Coefficient[ySeries, omega, 4]`, `gamma5 = (Coefficient[ySeries, omega, 5]/I)*k0`.
- The target definitions (`.py:40-44` ↔ `.wl:48-52`), the three ratio checks (`.py:51-53` ↔ `.wl:60-62`), and the closure construction (`.py:65-67` ↔ `.wl:68-72`) are identical algebra.

Both engines build the identical rational `Y`, hand it to the same symbolic-series black box, and read the same coefficients. The only differences are surface syntax (`sp.im(...)` vs `.../I` is the same operation). This violates the second-engine policy: a `.wl` is required wherever Mathematica CAN independently verify the result, and here it can. The user authorized the independent-route rewrite on 2026-06-05.

**Required change (state of the requirement — YOU design the route):**
Re-author the Mathematica derivation so the low-frequency coefficients `K2`, `K4`, `Gamma5` are obtained by a method **structurally different** from the SymPy script's `Series[]`-on-the-full-rational-then-`Coefficient` choreography. The new route must reach the coefficients by genuinely different algebra, not by echoing the `.py`'s steps.

One acceptable independent method (you may use this or any equally-independent route of your choosing): expand the pole denominator analytically via the geometric series `1/(1 - u) = 1 + u + u^2 + O(u^3)` with `u = omega^2/omegaQ^2 + I*chiQ*sigmaCan*omega^5`, collect by powers of `omega` by hand (note `u` starts at order `omega^2`, so only `1 + u + u^2` contributes through order `omega^5`; `u^2 = omega^4/omegaQ^4 + O(omega^7)`), and read `K2`, `K4`, `Gamma5` from the hand-collected coefficients. This cross-checks the `.py`'s `Series[]` extraction by a method that does not call `Series` on the full `Y`. Other acceptable routes: residue/`SeriesCoefficient`-with-explicit-order-bookkeeping that does not reconstruct the `.py` choreography, or closed-form coefficient formulas derived from the geometric expansion. (Do NOT simply rename the existing `Series[yRet,...]` call — that is the transliteration being removed.)

**Acceptance criteria (the verifier will check ALL of these):**
1. The coefficient extraction does NOT use `Series[yRet, {omega, 0, n}]` on the full rational `yRet` as the source of `K2`/`K4`/`Gamma5` (that is the transliterated step). The new route is auditable as a different derivation method.
2. All original deliverables are still verified, by genuine (can-fail, non-tautological) checks, with the SAME values:
   - M1: `K2 = K0/(4*Omega^2)` form — verified via `K2/K2_target - NQ == 0`.
   - M2: `K4 = K0/(4*Omega^4)` form — verified via `K4/K4_target - NQ == 0`.
   - M3: `Gamma5 = chiQ * NQ * Gamma5_target` — verified via `Gamma5/Gamma5_target - chiQ*NQ == 0`.
   - M4: the closure `closure_ratio - (mhat0^2*chiQ*NQ - 1) == 0`, recovering the headline `mhat_0^2 chi_Q N_Q = 1`.
   where `K0_target = 64 G Omega^5/(45 c^5)`, `K2_target = K0_target/(4 Omega^2)`, `K4_target = K0_target/(4 Omega^4)`, `Gamma5_target = 2 G/(5 c^5)`, `NQ = K0/K0_target`.
3. `chiQ` remains a FREE real symbol (its pin to 1 is upstream at stage 097); do NOT constrain it positive or substitute 1.
4. No emitted/derived value changes versus the current committed Mathematica output (the headline forms and all PASS lines stay the same); only the *derivation method* becomes independent.
5. Script exits 0 with all checks passing.

**Self-test (perform before finalizing):**
- Geometric route arithmetic: `u = omega^2/omegaQ^2 + I*chiQ*sigmaCan*omega^5`; `1 + u + u^2` through `omega^5` gives `omega^2`-coeff `= 1/omegaQ^2`, `omega^4`-coeff `= 1/omegaQ^4`, `omega^5`-coeff `= I*chiQ*sigmaCan`. Then `Y = 3/4 + (1/4)*(...)` gives `K2 = K0/(4 omegaQ^2)`, `K4 = K0/(4 omegaQ^4)`, `Gamma5 = K0*chiQ*sigmaCan/4 = K0*chiQ*(9/8)/(4 omegaQ^5)`. Confirm these match the existing `.py`/`.wl` values before wiring the asserts.
- Confirm `K2/K2_target = K4/K4_target = K0/K0_target = NQ`, `Gamma5/Gamma5_target = chiQ*NQ`, and the closure residual over `Gamma5_target` reduces to `mhat0^2*chiQ*NQ - 1` — all as in the reference `.py`.

**Verification command:**
After Codex applies, the orchestrator runs `redteam exec-mathematica 100` (and `redteam exec-sympy 100` to confirm the untouched reference engine still passes) and confirms: (a) no `Series[yRet, ...]` coefficient-extraction remains as the source of the K-coefficients; (b) all four deliverable checks (M1–M4) still report `= 0`/PASS; (c) the script exits 0; (d) the new route is auditable as structurally distinct from the `.py`.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl`
- summary: Replaced the full-rational `Series` coefficient extraction with a hand-collected geometric expansion through `omega^5` while preserving the original checks and emitted forms.
- deviation: none
