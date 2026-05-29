---
unit_id: 052
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-22T23:19:50Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 052

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py:85-92`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl:58-59`

**Issue:**
The two "threshold equality" assertions test residuals of the form `zeta_0^(phys).subs(Omega0**2, zeta_req*Kphi0/KW) - zeta_req` (and the analogous `Kphi0` form). Because `Omega0_req_sq` and `Kphi0_req` are *defined* as the rearrangements that make `zeta_0^(phys) = zeta_req` hold, the substitution is an algebraic identity by construction and the residual is zero regardless of physics. The check cannot fail and therefore proves nothing about the rescue thresholds.

**Required change:**

1. In the SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`, replace the block at lines 85-92:

   Before (lines 85-92):
   ```python
   expect_zero(
       "threshold equality at fixed stiffness",
       zeta0_phys.subs(Omega0**2, Omega0_req_sq) - zeta_req,
   )
   expect_zero(
       "threshold equality at fixed overlap",
       zeta0_phys.subs(Kphi0, Kphi0_req) - zeta_req,
   )
   ```

   After:
   ```python
   Omega0_sq_sym = sp.symbols("Omega0_sq", positive=True, real=True)
   sol_Omega0_sq = sp.solve(
       (KW * Omega0_sq_sym / Kphi0) - zeta_req, Omega0_sq_sym
   )
   assert len(sol_Omega0_sq) == 1, "expected unique Omega0^2 solution"
   expect_zero(
       "solve(zeta_phys = zeta_req) for Omega0^2 - expected",
       sp.simplify(sol_Omega0_sq[0] - Omega0_req_sq),
   )

   sol_Kphi0 = sp.solve(
       (KW * Omega0**2 / Kphi0) - zeta_req, Kphi0
   )
   assert len(sol_Kphi0) == 1, "expected unique Kphi0 solution"
   expect_zero(
       "solve(zeta_phys = zeta_req) for Kphi0 - expected",
       sp.simplify(sol_Kphi0[0] - Kphi0_req),
   )
   ```

   The new block solves `zeta_phys = zeta_req` for the rescue variable from the equation itself, then compares the engine-derived solution against the pre-existing `Omega0_req_sq` / `Kphi0_req` definitions. The check is now substantive: a wrong `Omega0_req_sq` formula would no longer pass.

2. In the Mathematica script `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`, replace lines 58-59:

   Before:
   ```mathematica
   expectZero["threshold equality at fixed stiffness", (zetaPhys /. omega0^2 -> omega0ReqSq) - zetaReq];
   expectZero["threshold equality at fixed overlap", (zetaPhys /. kPhi0 -> kPhi0Req) - zetaReq];
   ```

   After:
   ```mathematica
   omegaSqSol = First[omega0Sq /. Solve[(kW omega0Sq/kPhi0) - zetaReq == 0, omega0Sq]];
   expectZero["solve(zeta_phys = zeta_req) for Omega0^2 - expected",
              FullSimplify[omegaSqSol - omega0ReqSq, Assumptions -> $Assumptions]];

   kPhi0Sol = First[kPhi0 /. Solve[(kW omega0^2/kPhi0) - zetaReq == 0, kPhi0]];
   expectZero["solve(zeta_phys = zeta_req) for Kphi0 - expected",
              FullSimplify[kPhi0Sol - kPhi0Req, Assumptions -> $Assumptions]];
   ```

   Note: introduce `omega0Sq` as a fresh symbol used only inside the `Solve` call; it is not added to `$Assumptions` because `Solve` works without an explicit positivity declaration on the dummy variable.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 052` and `redteam exec-mathematica 052` (single-seat). The new output transcripts must contain:

- SymPy: `solve(zeta_phys = zeta_req) for Omega0^2 - expected = 0` and `solve(zeta_phys = zeta_req) for Kphi0 - expected = 0`.
- Mathematica: `PASS: solve(zeta_phys = zeta_req) for Omega0^2 - expected` and `PASS: solve(zeta_phys = zeta_req) for Kphi0 - expected`.
- The old `threshold equality at fixed stiffness` / `threshold equality at fixed overlap` lines must be absent.
- Both scripts must exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`
- summary: Replaced the tautological threshold substitutions with solver-derived Omega0^2 and Kphi0 rescue checks.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl:33-71`

**Issue:**
The Mathematica script is a near-line-by-line port of the SymPy script. The "expected" closed forms `dZExpected`, `deltaExpected`, and `softExpected` are typed in verbatim from the SymPy file (with only camelCase variable substitutions). The engine never derives these forms independently; it merely subtracts them from the same engine's `D[zetaReq, piTr]`, `zetaReq - 1`, and `1 - kPhi0Req/kW`, simplifies the difference, and confirms zero. Any error in the hand-typed closed form would survive both engines.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`, introduce independent derivations for the closed forms before each `expectZero` consumes them. Make the following insertions; line numbers refer to the file *before* any other edits.

1. Replace line 34
   ```mathematica
   zetaReq = FullSimplify[(sReq - 1)/(1 + eps (sReq - 2)), Assumptions -> $Assumptions];
   ```
   with:
   ```mathematica
   (* Derive zetaReq by solving the lowest-lane support equation directly. *)
   zetaSym = zetaSym;  (* fresh symbol *)
   zetaSolList = Solve[(sReq - 1) - zetaSym (1 + eps (sReq - 2)) == 0, zetaSym];
   If[Length[zetaSolList] =!= 1, fail["unique zeta solution"]];
   zetaReq = FullSimplify[zetaSym /. First[zetaSolList], Assumptions -> $Assumptions];
   ```

2. Replace lines 41-44 (the `dZdPi` / `dZExpected` block)
   ```mathematica
   dZdPi = FullSimplify[D[zetaReq, piTr], Assumptions -> $Assumptions];
   dZExpected = FullSimplify[cMix (1 - eps)/(cMix - eps (2 cMix - piTr))^2, Assumptions -> $Assumptions];
   Print["d zeta_req / d Pi_tr = ", fmt[dZdPi]];
   expectZero["dzeta_req/dPi - expected", dZdPi - dZExpected];
   ```
   with:
   ```mathematica
   dZdPi = FullSimplify[D[zetaReq, piTr], Assumptions -> $Assumptions];
   (* Alternative path: logarithmic differentiation of the rational form. *)
   zetaTogether = Together[zetaReq];
   numZ = Numerator[zetaTogether];
   denZ = Denominator[zetaTogether];
   dZdPiAlt = FullSimplify[
     zetaReq (D[numZ, piTr]/numZ - D[denZ, piTr]/denZ),
     Assumptions -> $Assumptions];
   expectZero["dZdPi vs dZdPiAlt (independent path)", dZdPi - dZdPiAlt];
   dZExpected = FullSimplify[cMix (1 - eps)/(cMix - eps (2 cMix - piTr))^2, Assumptions -> $Assumptions];
   Print["d zeta_req / d Pi_tr = ", fmt[dZdPi]];
   expectZero["dzeta_req/dPi - expected", dZdPi - dZExpected];
   ```

3. Replace lines 46-49 (the `deltaZeta` / `deltaExpected` block)
   ```mathematica
   deltaZeta = FullSimplify[zetaReq - 1, Assumptions -> $Assumptions];
   deltaExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(cMix - eps (2 cMix - piTr)), Assumptions -> $Assumptions];
   Print["Delta_zeta = ", fmt[deltaZeta]];
   expectZero["Delta_zeta - expected", deltaZeta - deltaExpected];
   ```
   with:
   ```mathematica
   (* Build deltaZeta via Together-based algebra, independent of the SymPy form. *)
   deltaZetaDerived = FullSimplify[Together[zetaReq - 1], Assumptions -> $Assumptions];
   deltaExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(cMix - eps (2 cMix - piTr)), Assumptions -> $Assumptions];
   Print["Delta_zeta = ", fmt[deltaZetaDerived]];
   expectZero["Delta_zeta - expected", deltaZetaDerived - deltaExpected];
   ```

4. Replace line 65 (the `softExpected` definition)
   ```mathematica
   softExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(piTr - cMix), Assumptions -> $Assumptions];
   ```
   with:
   ```mathematica
   (* Independent derivation of the softening fraction via Together[1 - 1/zetaReq]. *)
   softFracDerived = FullSimplify[Together[1 - 1/zetaReq], Assumptions -> $Assumptions];
   softExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(piTr - cMix), Assumptions -> $Assumptions];
   expectZero["softFrac vs Together[1 - 1/zetaReq] (independent path)", softFrac - softFracDerived];
   ```

   This adds an engine-independent cross-check between the loop-driven `softFrac` and the `Together`-derived form, *before* both are compared to the hand-typed `softExpected`.

Keep all other lines unchanged. Do not modify the SymPy script for F2.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 052`. The new output transcript must contain:

- `PASS: dZdPi vs dZdPiAlt (independent path)`
- `PASS: softFrac vs Together[1 - 1/zetaReq] (independent path)`
- The pre-existing `PASS: dzeta_req/dPi - expected`, `PASS: Delta_zeta - expected`, and `PASS: softening fraction - expected` lines must remain.
- Script must exit 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`
- summary: Added Mathematica-side independent derivations for zetaReq, dZdPi, Delta_zeta, and the softening fraction before the expected-form checks.
- deviation: none
