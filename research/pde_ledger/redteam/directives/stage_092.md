---
unit_id: 092
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 092

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding requires. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl:33-68`

**Issue:**
The Mathematica derivation is structurally identical to the SymPy derivation: same series-expand, same `Solve[branch == 0, kg0]`, same intermediate variables, same print labels. The second-engine policy requires an independent algebraic path. Restructure the body so the derivation works directly in the dimensionless `(eps2, eps4)` variables and verifies the branch identity there, rather than going through `Solve[branch == 0, kg0]`.

**Required change:**

Replace lines 33–68 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl` with the following body (preserve lines 1–32 and 70–73 unchanged). The new body works in the dimensionless variables from the start and verifies the same identities by an algebraically different route — it derives `c_pole` from `K_0 = 4 K_pole (1+eps_2)^2/(1+eps_4)` directly (notes Section 3) rather than from `Solve[K_0 K_4 - 4 K_2^2 == 0, kg0]`:

```mathematica
(* Work in dimensionless (eps2, eps4) variables from the outset.
   Notes Section 3 gives K_0 = 4 K_pole (1+eps_2)^2 / (1+eps_4) on the branch.
   We verify this is consistent with the original symbolic K_0, K_2, K_4
   coefficients of the conservative-carrier expansion, then read off c_pole. *)

(* Step A: define K_0, K_2, K_4 from the conservative-carrier expansion,
   then substitute eps2, eps4 directly without first solving for kg0. *)
k0Full = kg0 + kPole;
k2Full = kg2 + kPole/omegaQ^2;
k4Full = kg4 + kPole/omegaQ^4;

(* Dimensionless rewrite: kg2 = eps2 kPole/omegaQ^2, kg4 = eps4 kPole/omegaQ^4. *)
k2Eps = FullSimplify[k2Full /. {kg2 -> eps2*kPole/omegaQ^2}, Assumptions -> $Assumptions];
k4Eps = FullSimplify[k4Full /. {kg4 -> eps4*kPole/omegaQ^4}, Assumptions -> $Assumptions];
Print["K_2 in eps variables = ", fmt[k2Eps]];
Print["K_4 in eps variables = ", fmt[k4Eps]];

(* The branch identity K_0 K_4 = 4 K_2^2 gives K_0 = 4 K_2^2/K_4.
   Compute this directly without Solve[]. *)
k0FromBranch = FullSimplify[4*k2Eps^2/k4Eps, Assumptions -> $Assumptions];
Print["K_0 from branch (eps form) = ", fmt[k0FromBranch]];

(* Notes Section 3 prediction: K_0 = 4 K_pole (1+eps2)^2 / (1+eps4). *)
k0Predicted = FullSimplify[4*kPole*(1 + eps2)^2/(1 + eps4), Assumptions -> $Assumptions];
expectZero["K_0 closed form matches 4 K_pole (1+eps2)^2/(1+eps4)", k0FromBranch - k0Predicted];

(* c_pole = K_pole / K_0. *)
cPole = FullSimplify[kPole/k0FromBranch, Assumptions -> $Assumptions];
cPoleExpected = FullSimplify[(1 + eps4)/(4*(1 + eps2)^2), Assumptions -> $Assumptions];
Print["c_pole = ", fmt[Factor[cPole]]];
expectZero["c_pole - (1+eps4)/(4(1+eps2)^2)", cPole - cPoleExpected];

(* Static-geometry limit: eps2 = eps4 = 0 should give c_pole = 1/4. *)
expectZero["static-geometry limit c_pole = 1/4", (cPole /. {eps2 -> 0, eps4 -> 0}) - 1/4];

cGeom = FullSimplify[1 - cPole, Assumptions -> $Assumptions];
Print["c_geom in (eps2,eps4) variables = ", fmt[Factor[cGeom]]];
expectZero["static-geometry limit c_geom = 3/4", (cGeom /. {eps2 -> 0, eps4 -> 0}) - 3/4];

(* Small-contamination first-order expansion (notes Section 4). *)
smallSeries = Expand[Normal[Series[Normal[Series[cPoleExpected, {eps2, 0, 1}]], {eps4, 0, 1}]]];
linearPart = Expand[(1/4)*(1 + eps4 - 2*eps2)];
Print["First-order expansion of c_pole = ", fmt[smallSeries]];
Print["Linear part = ", fmt[linearPart]];
Print["Dropped higher-order tail = ", fmt[Expand[smallSeries - linearPart]]];
```

After this replacement, line counts shift; line numbers below F2 refer to the post-F1 state. (See F2 for the additional first-order assertions to append to this body.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 092` and confirm: (a) the script exits 0; (b) the transcript contains `PASS: K_0 closed form matches 4 K_pole (1+eps2)^2/(1+eps4)`; (c) the transcript contains `PASS: c_pole - (1+eps4)/(4(1+eps2)^2)`; (d) the transcript contains `PASS: static-geometry limit c_pole = 1/4` and `PASS: static-geometry limit c_geom = 3/4`; (e) the script no longer contains the substring `Solve[branch`.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py:78` (append after this line, before the FINAL LEDGER banner at line 80)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl` (after the small-series printout block introduced/preserved by F1; before the final `Print[""]` at line 70 of the original file)

**Issue:**
The first-order small-`(eps_2, eps_4)` expansion of `c_pole` is computed but no assertion enforces it. Add explicit coefficient checks so a regression in the predicted first-order behavior fails the script.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py`, insert the following three lines immediately after line 78 (`print("Dropped higher-order tail  =", remainder)`), before the `banner("FINAL LEDGER")` call at line 80:

```python
expect_zero("first-order eps^0 coefficient", small_series.coeff(eps2, 0).coeff(eps4, 0) - sp.Rational(1, 4))
expect_zero("first-order eps2 coefficient", small_series.coeff(eps2, 1).coeff(eps4, 0) - sp.Rational(-1, 2))
expect_zero("first-order eps4 coefficient", small_series.coeff(eps4, 1).coeff(eps2, 0) - sp.Rational(1, 4))
```

In the Mathematica file (`moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl`), after the `"Dropped higher-order tail = "` print (which is the last statement of the F1 replacement body), insert these three lines before the final `Print[""]` / `Print["Stage 092 Mathematica audit passed."]` block:

```mathematica
expectZero["first-order eps^0 coefficient", Coefficient[Coefficient[smallSeries, eps2, 0], eps4, 0] - 1/4];
expectZero["first-order eps2 coefficient", Coefficient[Coefficient[smallSeries, eps2, 1], eps4, 0] - (-1/2)];
expectZero["first-order eps4 coefficient", Coefficient[Coefficient[smallSeries, eps4, 1], eps2, 0] - 1/4];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 092` and `redteam exec-mathematica 092`. The sympy transcript should now contain three new lines `first-order eps^0 coefficient = 0`, `first-order eps2 coefficient = 0`, `first-order eps4 coefficient = 0`. The Mathematica transcript should contain three new `PASS:` lines with the same names. Both scripts must exit 0.

---

## Applied: F1+F2 (orchestrator-direct)

- files_changed: scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py, mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl
- summary: Restructured Mathematica derivation to work in dimensionless (eps2, eps4) variables from the start, deriving K_0 from 4 K_2^2/K_4 directly (no Solve[branch == 0, kg0]); added 3 first-order series coefficient asserts on both engines (1/4, -1/2, 1/4); plus banner sweep STAGE 75/075 → STAGE 092.
- deviation: none
