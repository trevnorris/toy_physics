---
unit_id: 176
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-30T01:12:20-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 176

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py:100-103`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl:92-94`

**Issue:** The rigidity-corollary check (paper deliverable D4 — the square-root mixed-leg law) is verified against the *constructed* `Sigma_factored`, so it is true by construction. Setting `dlnI -> 0`, `dlnH -> 0` inside `Sigma_factored = 2*dlnM + wI*dlnI + wH*dlnH` and asserting the result equals `2*dlnM` is the vacuous statement `a + b*0 + c*0 == a`; it cannot fail regardless of the physics. The fix re-derives the rigidity reduction from the *independent* `Sigma_exact` (the series/derivative of the actual `Lambda` definition), which is non-tautological.

**Required change (SymPy):**
Replace lines 100-103:
```python
banner("Rigidity corollary")

Sigma_rigid = sp.simplify(Sigma_factored.subs({dlnI: 0, dlnH: 0}))
expect_zero("rigidity reduction to 2 d ln M", Sigma_rigid - 2 * dlnM)
```
with a reduction taken from `Sigma_exact` under the primitive rigidity constraints
`dlnI = 0` (i.e. `dR + dGU - dOU - dGW = 0`) and `dlnH = 0` (i.e. `2*dR - dOU - dOW = 0`):
```python
banner("Rigidity corollary")

# Rigidity conditions delta ln I = 0 and delta ln H = 0 expressed on the
# primitive log-drifts (solve dlnI, dlnH for two of them):
#   dlnI = dR + dGU - dOU - dGW = 0  ->  dGU = dOU + dGW - dR
#   dlnH = 2 dR - dOU - dOW   = 0  ->  dOW = 2 dR - dOU
rigid = {dGU: dOU + dGW - dR, dOW: 2 * dR - dOU}

# Reduce the INDEPENDENT exact drift (not the constructed factored form)
# under rigidity and confirm it equals 2 d ln M reduced the same way.
Sigma_exact_rigid = sp.simplify(Sigma_exact.subs(rigid))
dlnM_rigid = sp.simplify((2 * dlnM).subs(rigid))
expect_zero("rigidity reduction of Sigma_exact to 2 d ln M", Sigma_exact_rigid - dlnM_rigid)
```
(`Sigma_exact` is already defined at lines 62-72; `dlnM`, `dGW`, `dOW`, `dOU`, `dR`, `dGU`, `dK` are already in scope.)

**Required change (Mathematica):**
Replace lines 92-94:
```mathematica
banner["Rigidity corollary"];
sigmaRigid = FullSimplify[sigmaFactoredForm /. {dlnI -> 0, dlnH -> 0}, Assumptions -> $Assumptions];
expectZero["rigidity reduction to 2 d ln M", sigmaRigid - 2*dlnM];
```
with the analogous reduction of the independent `sigmaExact`:
```mathematica
banner["Rigidity corollary"];
(* Rigidity: dlnI = dR + dGU - dOU - dGW = 0 -> dGU -> dOU + dGW - dR;
            dlnH = 2 dR - dOU - dOW = 0     -> dOW -> 2 dR - dOU.        *)
rigid = {dGU -> dOU + dGW - dR, dOW -> 2*dR - dOU};
sigmaExactRigid = FullSimplify[sigmaExact /. rigid, Assumptions -> $Assumptions];
dlnMRigid = FullSimplify[(2*dlnMExpr) /. rigid, Assumptions -> $Assumptions];
expectZero["rigidity reduction of Sigma_exact to 2 d ln M", sigmaExactRigid - dlnMRigid];
```
(`sigmaExact` is defined at lines 60-63; `dlnMExpr` at line 66. Note: use `dlnMExpr` here, not the bare symbol `dlnM`, since in the `.wl` `dlnM` is an unassigned symbol and only `dlnMExpr` carries the expression.)

After editing, also update the conclusion `Print` text only if it references the old check name; the existing conclusion prose (sympy:105-109, math:96-101) is still correct and may stay.

**Verification command:**
The verifier runs `redteam exec-sympy 176` and `redteam exec-mathematica 176` and confirms: (1) a check named `rigidity reduction of Sigma_exact to 2 d ln M` appears in both transcripts with residual `= 0`; (2) the residual is computed from `Sigma_exact`/`sigmaExact`, not from `Sigma_factored`/`sigmaFactoredForm`; (3) both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl`
- summary: Replaced the rigidity corollary with a reduction of the independently derived exact drift under primitive rigidity constraints.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl:59-63`

**Issue:** The `.wl` is a near line-by-line port of the `.py`. The one genuine independence is that the first-order coefficient is extracted by symbolic differentiation (`D[Log[...], eps] /. eps -> 0`) in Mathematica vs. Taylor series (`sp.series(...).removeO()/eps`) in SymPy. This divergence is currently undocumented, so the independence is invisible to a reader. Severity is low (the deliverables are pure algebraic identities); the only required remediation is to make the existing divergence explicit. Do NOT rewrite the algebra or add new physics.

**Required change:**
Add a one-line comment immediately above the `sigmaExact = ...` block (currently line 60) noting the engines use different first-order extraction methods, e.g.:
```mathematica
(* Independent first-order extraction: Mathematica uses symbolic D[Log,eps]
   at eps->0; the SymPy audit instead Taylor-expands and reads the eps^1 term. *)
sigmaExact = FullSimplify[
  D[Log[((lambdaP^2/kP)/(lambda^2/k))], eps] /. eps -> 0,
  Assumptions -> $Assumptions
];
```
Make no other change for F2.

**Verification command:**
The verifier confirms the `.wl` retains `D[Log[...], eps]` (distinct from the SymPy `series` route), the clarifying comment is present, and the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl`
- summary: Added the requested comment documenting Mathematica's derivative-based first-order extraction as distinct from the SymPy Taylor-series route.
- deviation: none

## F3 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py:32`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl:26`

**Issue:** Both scripts' top banner mislabels the stage as "STAGE 159"; this is Stage 176. The mislabel propagates into both saved transcripts (line 11 of each).

**Required change:**
- SymPy line 32: change `banner("STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION")` to `banner("STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION")`.
- Mathematica line 26: change `banner["STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION"];` to `banner["STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION"];`.
Re-run both scripts so the transcripts pick up the corrected banner.

**Verification command:**
The verifier confirms both transcripts' line 11 reads `STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION` and both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl`
  - `scripts/output/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.txt`
- summary: Corrected both top banners to Stage 176 and regenerated the saved transcripts with the corrected banner.
- deviation: none
