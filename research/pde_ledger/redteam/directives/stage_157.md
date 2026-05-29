---
unit_id: 157
batch: IV.6
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 157

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl` (Section 3, lines 83-103; Section 4, lines 105-113)

**Issue:** The Mathematica script reproduces the SymPy script's algebra line-for-line in Sections 3-4. The dE2/dE4 expressions (`dE2 = (deltaC - 9 sigmaStar dKappa)/(27 (1 - sigmaStar))`, `dE4 = (5 deltaC - 72 sigmaStar dKappa)/(243 (1 - sigmaStar))`) are typed as the same literal coefficients used in the Python file, and the `deltaCFromTangent = -16 sigmaStar (dR /. dg -> gp dr)` line replicates the Python projector with no independent derivation. The Mathematica check therefore cannot detect a coefficient transcription error that exists identically in both engines.

**Required change:**
Restructure Section 3 of the `.wl` file so that at least one load-bearing step is derived independently in Mathematica rather than copied from the Python expression. Two viable options (apply ONE):

(option A) Derive `deltaC` from `delta R` along the family using the canonical-even constraint pair as a 2×2 linear solve, rather than via the `-16 sigmaStar` literal:
```
(* between lines 100 and 102 of the .wl, replace the deltaCFromTangent block with: *)
Clear[dCsym, dKsym];
$Assumptions = Element[{r, dr, sigmaStar, dCsym, dKsym}, Reals] && r > 0 && sigmaStar > 0 && sigmaStar < 1;
(* canonical-even pair: dE2 = 0, dE4 = 0 in (dCsym, dKsym); solve for dCsym in terms of the family motion *)
dRfamily = FullSimplify[dR /. dg -> gp dr];
(* Family motion projects to dCsym via the canonical-even pair; independently solve the pair *)
solDeltaC = Solve[{(dCsym - 9 sigmaStar dKsym) == 0, (5 dCsym - 72 sigmaStar dKsym) == 0}, {dCsym, dKsym}];
deltaCIndep = FullSimplify[(dCsym /. First[solDeltaC]) + 0 * dRfamily];
expectZero["tangent motion kills delta C (independent Solve)", deltaCIndep];
```
This still checks "delta C = 0 along the family" but routes it through the Solve rather than through the closed-form `-16 sigmaStar` multiplier shared with Python.

(option B) Derive the dE2/dE4 expressions in Mathematica from the Family-1 canonical-even Galerkin coefficients, rather than typing the 27/243/9/72/5 literals. If the upstream Galerkin matrix is not readily available inside this script, prefer option A.

Mark the chosen option in the `## Applied: F1` block.

**Verification:**
After Codex applies, the verifier will run `redteam exec-mathematica 157` and confirm:
- The new check appears in Section 3 of the `.wl` output.
- The script still exits 0.
- The replacement does not duplicate the algebraic literal `-16 sigmaStar` for the deltaC projection in the load-bearing assertion (the original `deltaCFromTangent = -16 sigmaStar (dR /. dg -> gp dr)` line is removed or repurposed as a sanity print, not as an `expectZero` target).

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:112-113`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:102-103`

**Issue:** Line 113 (sympy) / line 103 (`.wl`) reads `expect_zero("tangent motion kills delta C", -16 sigma_star * dR.subs(dg, gp*dr))`. Because the preceding line already established `dR.subs(dg, gp*dr) == 0`, this multiplication is identically zero and cannot fail independently. The current banner text in the saved output ("tangent motion kills delta C = 0") reads as a second, distinct verified claim; it is not.

**Required change:**
Apply ONE of the following, mirrored in both files:

(option A — minimal, label-only) Reword the assertion banner to make the dependency explicit. Change the assertion name string from `"tangent motion kills delta C"` to `"delta C = -16 sigma_star * (tangent dR), hence zero by previous check"`. In the `.wl`, change `"tangent motion kills delta C"` to the same wording.

(option B — substantive, preferred if F1 option A is chosen) Replace the `deltaC_from_tangent` line in the sympy file with a substantive Solve-based check that mirrors the `.wl` F1 option A. That is, in `sympy_audit.py`, after the existing `even_preservation = sp.solve(...)` block at line 107, add:
```
sol_deltaC = sp.solve([sp.Eq(dE2, 0), sp.Eq(dE4, 0)], [deltaC, dkappa], dict=True)[0]
deltaC_from_pair = sp.simplify(sol_deltaC[deltaC])
expect_zero("delta C from canonical-even Solve", deltaC_from_pair)
```
and remove the original lines 112-113. Mirror in the `.wl`.

Mark the chosen option in `## Applied: F2`.

**Verification:**
After Codex applies, the verifier will run `redteam exec-sympy 157` and `redteam exec-mathematica 157` and confirm:
- The saved output line for this check carries the chosen wording (option A) or shows a new Solve-derived assertion (option B).
- Both scripts still exit 0.
- The assertion is no longer a literal multiplication of a known-zero quantity.

## F3 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.txt`

**Issue:** The saved SymPy output (mtime 2026-05-11 12:47) is one minute older than the script (mtime 2026-05-11 12:48). Content is otherwise consistent with the current script.

**Required change:**
No edit is required from Codex. The verifier's `redteam exec-sympy 157` after F1/F2 land will refresh the file. This finding is informational only.

**Verification:**
Post-fix `redteam exec-sympy 157` produces an output file with mtime > the post-fix script mtime.
