---
unit_id: 031
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T05:43:14Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 031

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. (None in this directive.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl` (full rewrite of body, lines 26-83; preserve the existing helpers/banner block at lines 1-24).

**Issue:**
The Mathematica script is a structural line-by-line port of the SymPy script. Both files use the same `(2 A + DK - alpha sigma ± R)/2` parameterization of the eigenvalues, the same Hellmann-Feynman shortcut `s_- := -D[lam_-, alpha]`, the same five-part decomposition with identical banner strings (`PART I` – `PART V`, including the typo `PREFATOR`), the same intermediate variable choreography (`sigma`, `deltaKappa`, `kappaProd`, `r`, `lamMinus`, `lamPlus`, `sMinus`, `dsExact`, `dsExpected`, `p0Sel`, `t0`, `alphaCrit`, `radCritDerived`, `p0Factored`), and the same closed-form expressions typed in the same order. The two engines are not independently re-deriving the result from physical premises — the `.wl` is transliterating the `.py`'s algebraic choreography. The second-engine policy requires each engine to derive results independently from the physical premises (here, from `Eigenvalues`/`Eigenvectors` of the loaded 2x2 wall matrix), so the `.wl` must be restructured.

**Required change:**

Replace the body of `moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl` from line 26 through line 82 (everything between the helpers/banner block and the final `Exit[0];`) with the following *independent-derivation* sequence. The helpers `banner`, `pass`, `fmt`, `fail`, `expectZero` at lines 4-24 are unchanged. The final `Exit[0];` at line 83 is unchanged.

Replacement body (insert in place of current lines 26-82, keeping lines 1-24 and line 83 intact):

```
banner["PART I — INDEPENDENT EIGEN-DECOMPOSITION AND HF IDENTITY"];

Clear[a, dK, alpha, kappa0, kappa1, beta0];
$Assumptions = Element[{a, dK, alpha, kappa0, kappa1, beta0}, Reals] &&
  a > 0 && dK > 0 && alpha >= 0 && kappa0 > 0 && kappa1 > 0 && beta0 > 0;

(* Loaded 2x2 wall matrix, built directly from the physical setup:
   diagonal stiffnesses (a, a+dK) minus rank-one loading alpha * v v^T
   with loading vector v = {kappa0, kappa1}. This is the rank-one loaded
   wall problem of Stage 028 expressed without the (sigma, deltaKappa)
   choreography of the SymPy mirror. *)
wallMatrix = {{a, 0}, {0, a + dK}} -
  alpha*Outer[Times, {kappa0, kappa1}, {kappa0, kappa1}];

(* Independent eigenvalue extraction via Eigenvalues, then identify
   the lower branch as the one that reduces to a at alpha=0. *)
eigvals = Sort[Eigenvalues[wallMatrix]];
lamMinusIndep = eigvals[[1]];
lamPlusIndep  = eigvals[[2]];
expectZero["lam_-(0) initial value", (lamMinusIndep /. alpha -> 0) - a];
expectZero["lam_+(0) initial value", (lamPlusIndep  /. alpha -> 0) - (a + dK)];

(* Independent eigenvector extraction, then loading-vector overlap.
   The selected overlap s_- = (v . e_-)^2 / (e_- . e_-) is the
   physical definition used in the paper notes. *)
eigvecs = Eigenvectors[wallMatrix];
(* Pair eigenvectors with eigenvalues so we always pick the lower one *)
pairs = Sort[Transpose[{Eigenvalues[wallMatrix], eigvecs}], #1[[1]] < #2[[1]] &];
eMinus = pairs[[1, 2]];
loadingVec = {kappa0, kappa1};
sMinusOverlap = FullSimplify[
  (loadingVec . eMinus)^2 / (eMinus . eMinus),
  Assumptions -> $Assumptions
];

(* HF identity: d lam_- / d alpha = - s_-. Independent verification,
   not a definition. *)
hfResidual = FullSimplify[
  D[lamMinusIndep, alpha] + sMinusOverlap,
  Assumptions -> $Assumptions
];
expectZero["Hellmann-Feynman d lam_-/d alpha + s_- = 0", hfResidual];

banner["PART II — EXACT OVERLAP DERIVATIVE FROM OVERLAP DEFINITION"];

dsOverlap = FullSimplify[D[sMinusOverlap, alpha], Assumptions -> $Assumptions];
rSquared = FullSimplify[(dK + alpha*(kappa0^2 - kappa1^2))^2 + 4*alpha^2*kappa0^2*kappa1^2,
                        Assumptions -> $Assumptions];
dsExpected = FullSimplify[2*dK^2*kappa0^2*kappa1^2 / rSquared^(3/2),
                          Assumptions -> $Assumptions];
expectZero["ds_-/dalpha closed form", FullSimplify[dsOverlap - dsExpected, Assumptions -> $Assumptions]];
Print["ds_-/dalpha = ", fmt[dsExpected]];

banner["PART III — EXACT PREFACTOR DERIVATIVE"];

p0SelIndep = FullSimplify[beta0*sMinusOverlap/lamMinusIndep, Assumptions -> $Assumptions];
dPDirect = FullSimplify[D[p0SelIndep, alpha], Assumptions -> $Assumptions];
dPClosed = FullSimplify[
  beta0*(dsExpected*lamMinusIndep + sMinusOverlap^2)/lamMinusIndep^2,
  Assumptions -> $Assumptions
];
expectZero["dP0_-/dalpha closed form", FullSimplify[dPDirect - dPClosed, Assumptions -> $Assumptions]];
Print["dP0_-/dalpha = ", fmt[dPClosed]];

banner["PART IV — INITIAL VALUES AT alpha=0"];

expectZero["s_-(0) = kappa0^2", (sMinusOverlap /. alpha -> 0) - kappa0^2];
expectZero["P0_-(0) = beta0 kappa0^2 / a", (p0SelIndep /. alpha -> 0) - beta0*kappa0^2/a];

banner["PART V — SOFTENING THRESHOLD FROM Det[M]=0"];

(* alpha_crit is the root of Det[M] = 0, derived from the matrix directly. *)
detM = FullSimplify[Det[wallMatrix], Assumptions -> $Assumptions];
critSolutions = Solve[detM == 0, alpha];
(* The stable-side threshold is the smaller positive root. Both solutions are
   non-negative; pick the one that equals the closed form. *)
alphaCritClosed = FullSimplify[
  a*(a + dK) / ((a + dK)*kappa0^2 + a*kappa1^2),
  Assumptions -> $Assumptions
];
alphaCritFromDet = FullSimplify[alpha /. First[critSolutions], Assumptions -> $Assumptions];
expectZero["alpha_crit from Det[M]=0", FullSimplify[alphaCritFromDet - alphaCritClosed, Assumptions -> $Assumptions]];
Print["alpha_crit = ", fmt[alphaCritClosed]];

expectZero["lam_-(alpha_crit) = 0",
  FullSimplify[lamMinusIndep /. alpha -> alphaCritClosed, Assumptions -> $Assumptions]
];

banner["PART VI — STABLE-SIDE POLE STRUCTURE"];

(* lam_- lam_+ = Det[M] (basic linear algebra), so factorization through
   alpha_crit is automatic from Det[M]. *)
detExpanded = FullSimplify[detM, Assumptions -> $Assumptions];
t0Closed = FullSimplify[(a + dK)*kappa0^2 + a*kappa1^2, Assumptions -> $Assumptions];
expectZero["Det[M] - t0(alpha_crit - alpha)",
  FullSimplify[detExpanded - t0Closed*(alphaCritClosed - alpha), Assumptions -> $Assumptions]
];
expectZero["lam_- lam_+ - Det[M]",
  FullSimplify[lamMinusIndep*lamPlusIndep - detExpanded, Assumptions -> $Assumptions]
];

p0Factored = FullSimplify[
  beta0*sMinusOverlap*lamPlusIndep / (t0Closed*(alphaCritClosed - alpha)),
  Assumptions -> $Assumptions
];
expectZero["P0_- pole factorization",
  FullSimplify[p0SelIndep - p0Factored, Assumptions -> $Assumptions]
];
Print["P0_- factored = ", fmt[p0Factored]];

banner["STAGE 031 INDEPENDENT MATHEMATICA AUDIT COMPLETE"];
Print["Verified (independently from Eigenvalues/Eigenvectors of the loaded 2x2 wall matrix):"];
Print["  HF identity d lam_-/d alpha + s_- = 0 (derived, not assumed)"];
Print["  exact overlap derivative ds_-/dalpha = 2 dK^2 kappa0^2 kappa1^2 / R^3"];
Print["  exact prefactor derivative formula"];
Print["  initial values at alpha=0"];
Print["  alpha_crit from Det[M] = 0"];
Print["  pole factorization of P0_- at alpha_crit"];
```

Notes for Codex:
- The replacement text uses `kappa0, kappa1` (not `x0, x1`) and builds the loaded matrix `wallMatrix` directly. This breaks the verbatim variable correspondence with the SymPy file.
- `Eigenvalues` and `Eigenvectors` calls are the substantive independent step — they let Mathematica derive `lam_-` and `e_-` itself rather than reading off the SymPy parametrization.
- The HF identity now appears as a *verified theorem* (`hfResidual = 0`), not as a definition of `s_-`. This converts what was a tacit shortcut in the SymPy file into an explicit cross-check.
- All assertions previously present (`ds_-/dalpha closed form`, `P0_- factorization`, etc.) are preserved with substantively equivalent names; the verifier should see all PASS lines.

**Claim manifest:** The replacement script must independently verify the following physical results, numbered M1-M7:
- M1 (HF identity): `D[lamMinusIndep, alpha] + sMinusOverlap == 0`, where `sMinusOverlap = (v . e_-)^2 / (e_- . e_-)` from `Eigenvectors[wallMatrix]`.
- M2 (overlap derivative closed form): `D[sMinusOverlap, alpha] == 2 dK^2 kappa0^2 kappa1^2 / R^3`.
- M3 (prefactor derivative closed form): `D[beta0 sMinusOverlap / lamMinusIndep, alpha] == beta0 (dsExpected lamMinusIndep + sMinusOverlap^2) / lamMinusIndep^2`.
- M4 (initial values): `sMinusOverlap(0) == kappa0^2` and `P0_-(0) == beta0 kappa0^2 / a`.
- M5 (alpha_crit from Det[M]): `alpha_crit` extracted from `Solve[Det[wallMatrix] == 0, alpha]` matches the closed form `a(a + dK) / ((a+dK) kappa0^2 + a kappa1^2)`.
- M6 (eigenvalue collapse): `lamMinusIndep /. alpha -> alpha_crit == 0`.
- M7 (pole factorization): `P0_- == beta0 sMinusOverlap lamPlusIndep / (t0 (alpha_crit - alpha))`.

**Self-test check (Codex should not skip):**
- `Eigenvalues[{{a, 0}, {0, a + dK}}]` (alpha=0 case) returns `{a, a + dK}` — `Sort` picks `a` as the lower branch. This guarantees `lamMinusIndep(0) = a`, satisfying the existing initial-value claim.
- `Det[wallMatrix] = a(a+dK) - alpha((a+dK) kappa0^2 + a kappa1^2)`. `Solve[detM == 0, alpha]` returns a single root `alpha = a(a+dK)/((a+dK) kappa0^2 + a kappa1^2)`, which matches `alphaCritClosed` and verifies M5 non-tautologically.
- `Eigenvectors[wallMatrix]` for the rank-one perturbation is well-defined and Mathematica returns closed forms for both eigenvectors; the overlap `(v . e_-)^2 / (e_- . e_-)` simplifies to the same expression the SymPy file gets via `-d lam_-/d alpha`, so M1 will pass non-vacuously.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 031` and confirm:
- The saved output shows `PASS` for at least the M1-M7 assertions (named in the replacement code).
- The file no longer contains the substrings `sigma = x0 + x1`, `deltaKappa = x0 - x1`, `kappaProd = x0*x1`, `radCritDerived`, or the `PART I — EXACT SELECTED OVERLAP DERIVATIVE` banner string (these were the smoking-gun markers of the transliteration).
- The script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl`
- summary: Replaced the Mathematica audit body with an independent Eigenvalues/Eigenvectors derivation from the loaded 2x2 wall matrix.
- deviation: none

---

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:60-69`

**Issue:**
The "generic quotient/HF identity" assertion at sympy lines 60-69 is algebraically forced by SymPy's quotient rule plus the substitution `dL/dalpha -> -Ssym`. It tests SymPy's symbolic-differentiation machinery, not the moving-throat physics. The companion direct identity at line 73 (`expect_zero("dP0_-/dalpha direct identity", sp.simplify(dP_direct - dP_physical))`) is the load-bearing check; the generic identity is redundant scaffolding. The Mathematica script (current line 47-49) already skips this generic check, confirming it is not needed for the second-engine cross-check.

**Required change:**

In `moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py`, delete lines 60-69 inclusive. The file currently reads:

```
P0_sel = sp.simplify(beta0 * s_minus / lam_minus)

# Verify the quotient-rule identity generically, using the exact Hellmann–Feynman law
# d lambda_-/d alpha = - s_- together with the explicit ds_-/dalpha formula from Part I.
Lsym, Ssym, DSsym = sp.symbols("Lsym Ssym DSsym", nonzero=True, real=True)
dP_generic = sp.diff(beta0 * sp.Function("S")(alpha) / sp.Function("L")(alpha), alpha)
dP_generic = dP_generic.subs({
    sp.diff(sp.Function("S")(alpha), alpha): DSsym,
    sp.Function("S")(alpha): Ssym,
    sp.diff(sp.Function("L")(alpha), alpha): -Ssym,
    sp.Function("L")(alpha): Lsym,
})
dP_expected = sp.simplify(beta0 * (DSsym * Lsym + Ssym**2) / Lsym**2)
expect_zero("generic quotient/HF identity", sp.simplify(dP_generic - dP_expected))

dP_direct = sp.diff(P0_sel, alpha)
```

Change it to:

```
P0_sel = sp.simplify(beta0 * s_minus / lam_minus)

dP_direct = sp.diff(P0_sel, alpha)
```

Concretely: keep line 56 (`P0_sel = ...`) and line 71 (`dP_direct = sp.diff(P0_sel, alpha)`) onward. Remove lines 58-69 inclusive (the comment block, the `Lsym, Ssym, DSsym` declaration, the `dP_generic` construction, the substitution dict, the `dP_expected` line, and the `expect_zero("generic quotient/HF identity", ...)` line). Also remove the blank line between line 69 and line 71 if any. Do not modify line 73 onward.

**Self-test check (Codex should not skip):**
- After deletion, `dP_direct = sp.diff(P0_sel, alpha)` is still well-defined because `P0_sel` (defined at line 56) and `alpha` (defined at line 35) both remain in scope. The remaining direct identity check at the now-line-65 (formerly line 73) is unchanged structurally.
- The symbols `Lsym, Ssym, DSsym` are not used anywhere else in the script (grep confirms only the deleted block references them), so removing the declaration is safe.
- No PART II output line beyond the direct identity is affected.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 031` and confirm:
- The saved sympy output's PART II section contains exactly `dP0_-/dalpha direct identity = 0` and `dP0_-/dalpha = ...` — no `generic quotient/HF identity = 0` line.
- The script exits 0.
- All other assertions in PART I, III, IV, V continue to pass.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py`
- summary: Removed the generic quotient/HF identity scaffolding so PART II relies on the direct physical prefactor derivative check.
- deviation: none
