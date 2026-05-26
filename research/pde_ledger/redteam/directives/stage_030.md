---
unit_id: 030
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-25T23:40:09-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 030

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:73` (insertion point — between the `lamPlus = ...` block and the `Print["lambda_- = ..."]` calls)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:80` (insertion point — between the `lam_plus = ...` definition and the `s_minus_hf = ...` line)

**Issue:**
The paper stage card's `\stagefield{Checks}` item 1 (paper/stages/stage_030.tex line 101) states the Hellmann-Feynman identity in full: `d lambda_-/d alpha_0 = -e_-^T v v^T e_- = -s_-`. The boxed output equation eq:app-stage030-s-minus (paper/stages/stage_030.tex line 45-50) likewise asserts `s_- = (v.e_-)^2 = (1/2)[sigma + ...]`. The scripts verify only the rightmost equality (closed form vs derivative), not the leftmost (eigenvector projection vs derivative). Neither engine constructs `e_-` or computes `(v.e_-)^2` from the explicit eigenvector. The Mathematica script already builds the wall matrix `mMat` at lines 55-56 and calls `Eigenvalues[mMat]` at line 63 — adding a parallel `Eigenvectors[mMat]` call to verify `(v.e_-)^2 = sMinusClosed` is one extra block. The SymPy script does not yet construct the matrix; it should construct the same matrix (basis ordering matching the Mathematica file) and compute the eigenvector.

**Required change:**

The HF identity has three sides — `(v.e_-)^2 = -d lambda_-/d alpha = closed-form` — and the existing checks only span the two rightmost. Add a new assertion to each engine that verifies the leftmost equality `(v.e_-)^2 = closed-form` by constructing the eigenvector of the loaded wall block explicitly.

The Mathematica script's `mMat` basis (math:55-56) places the kappa_1-mode in row 1 (diagonal `B - alpha*kappa_1^2 = a + dK - alpha*x1`) and the kappa_0-mode in row 2 (diagonal `A - alpha*kappa_0^2 = a - alpha*x0`); the off-diagonal is `-alpha*kappa_0*kappa_1 = -alpha*Sqrt[x0*x1]`. In this basis the loading direction `v = (kappa_0, kappa_1)^T` of the paper becomes `v = {Sqrt[x1], Sqrt[x0]}` (the first component carries the kappa_1-mode coefficient, the second carries the kappa_0-mode coefficient). The same matrix and basis ordering must be used in SymPy for the cross-check to match the Mathematica result.

*Mathematica* (`mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`)

Insert the following block immediately after the existing `lamPlus = FullSimplify[...]` block (current line 73) and before the existing `Print["lambda_- = ", fmt[lamMinus]];` (current line 75):

```mathematica
(* HF eigenvector check: compute the lower-eigenvalue eigenvector of mMat
   directly and verify (v.e_-)^2 = sMinusClosed. This closes the HF chain
   (v.e_-)^2 = -d lambda_-/d alpha = closed form. In the mMat basis (row 1
   = kappa_1-mode, row 2 = kappa_0-mode), the loading direction is
   vVec = {Sqrt[x1], Sqrt[x0]}. *)
eigPairs = Eigensystem[mMat];
eMinusRaw = First[Pick[Transpose[eigPairs][[All, 2]],
  Map[FullSimplify[# - lamMinus, Assumptions -> $Assumptions] === 0 &,
      First[eigPairs]]]];
eMinusNorm = FullSimplify[eMinusRaw/Sqrt[eMinusRaw.eMinusRaw],
  Assumptions -> $Assumptions];
vVec = {Sqrt[x1], Sqrt[x0]};
sMinusEig = FullSimplify[(vVec.eMinusNorm)^2, Assumptions -> $Assumptions];
```

Then, immediately after the existing `expectZero["selected overlap: HF - closed form", sMinusHF - sMinusClosed];` (current line 83), insert a new assertion:

```mathematica
expectZero["HF eigenvector check", sMinusEig - sMinusClosed];
```

*SymPy* (`scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`)

Insert the following block between the existing `lam_plus = sp.simplify(...)` line (current line 74) and the `print("lambda_- =")` line (current line 76):

```python
# HF eigenvector check: build the same loaded wall block used in the
# Mathematica script (basis row 1 = kappa_1-mode, row 2 = kappa_0-mode),
# extract the eigenvector at the lower eigenvalue, and verify (v.e_-)^2
# equals s_minus_closed. The matrix entries reproduce
# trace = 2 A + DK - alpha*(x0+x1) and
# det = A*(A+DK) - alpha*((A+DK)*x0 + A*x1).
M_eff = sp.Matrix([
    [A + DK - alpha * x1, -alpha * sp.sqrt(x0 * x1)],
    [-alpha * sp.sqrt(x0 * x1), A - alpha * x0],
])
v_vec = sp.Matrix([sp.sqrt(x1), sp.sqrt(x0)])
e_minus_raw = None
for val, _mult, vects in M_eff.eigenvects():
    if sp.simplify(val - lam_minus) == 0:
        e_minus_raw = vects[0]
        break
if e_minus_raw is None:
    raise AssertionError("HF eigenvector check: lower-branch eigenvector not found")
e_minus_norm = e_minus_raw / sp.sqrt((e_minus_raw.T * e_minus_raw)[0, 0])
s_minus_eig = sp.simplify(((v_vec.T * e_minus_norm)[0, 0]) ** 2)
```

Then, immediately after the existing `expect_zero("selected overlap: HF - closed form", s_minus_hf - s_minus_closed)` (current line 87), insert:

```python
expect_zero("HF eigenvector check", s_minus_eig - s_minus_closed)
```

The new code in both engines must be inserted, not overwritten — none of the existing assertions are removed.

**Self-test sketch** (done by auditor; reproduce mentally before applying):

- At `alpha = 0`, `M_eff -> diag(A + DK, A) = diag(B, A)`. Lower eigenvalue is `A`, with eigenvector `{0, 1}`. Then `v_vec.e_- = (Sqrt[x1], Sqrt[x0]).(0, 1) = Sqrt[x0]`, so `(v.e_-)^2 = x0`. The closed form `sMinusClosed` at `alpha = 0` likewise equals `x0` (already verified by the existing weak-loading check). So the new check passes at this limit, as required.
- For general alpha, the lower-eigenvalue eigenvector of `M_eff` is proportional to `{-Sqrt[x0*x1]*alpha, (A+DK-alpha*x1) - lam_minus}` (or equivalent up to scaling, depending on which row Mathematica/SymPy chooses). Normalizing and projecting against `v_vec` yields an alpha-dependent expression; the assertion `(v.e_-)^2 - sMinusClosed = 0` then reduces to an algebraic identity that FullSimplify / simplify will dispatch. This is genuinely independent of the `-d lam_-/d alpha` chain, so the assertion can fail (and would fail) under a sign or normalization mistake.
- No new symbol is introduced; no paper constant is added or changed; no existing check is touched.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 030` then `redteam exec-mathematica 030`. Expected outcomes:

1. The SymPy output transcript contains a new line `HF eigenvector check = 0` (positioned just after `selected overlap: HF - closed form = 0`).
2. The Mathematica output transcript contains a new line `HF eigenvector check = 0` followed by `PASS: HF eigenvector check`.
3. All existing assertions still pass with residual 0; no other captured output line is changed except for the new insertion.
4. Both scripts exit 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
  - `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`
- summary: Added explicit lower-eigenvector projection checks in both engines to verify `(v.e_-)^2 = s_-` independently of the derivative check.
- deviation: none
