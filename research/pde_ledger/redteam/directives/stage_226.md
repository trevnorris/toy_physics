---
unit_id: 226
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T16:49:38-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 226

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond creating the required Mathematica script. Do NOT touch `paper.tex`, `notes/`, or any prose document, and do NOT modify the existing SymPy script.

After creating the script, RUN it (`timeout 600 math -script <path>`) and iterate until it exits 0 with every in-file check (`expectZero` / `expect-close`) passing. Getting it to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl`

(Exact path above — `.wl` lives in `mathematica/`; do NOT drop the `_mathematica_audit` suffix.)

**Issue:** Stage 226 (`checkpoint: False`, not status-only) has only a SymPy witness. Every load-bearing claim — the load-defect bridge, the even-gate closed forms, the one-pole compatibility base data, the two mixed-sector corridor matrices with their ranks/nullities and null bases, the same-charge functional projections, and the two QR-projector corridor norms — is mechanizable independently in Mathematica. The dual-engine rule requires a second-engine re-derivation. The card already notes "Mathematica audit: none yet" (`paper/stages/stage_226.tex:11`); this directive supplies it.

**Required change:**
Create a NEW standalone Mathematica audit at the Target path that independently re-derives and checks the claim manifest M1–M7 below. The script must be self-contained (define all primitives from scratch), print a clear per-section log, and end by `Print`-ing a success line then `Exit[0]`; on any failed check it must `Print` the offending residual and `Exit[1]`.

Use a helper of the form:
```
expectZero[expr_, tag_] := If[PossibleZeroQ[Chop[N[expr, 30]]] || Abs[N[expr, 30]] < 10^-12,
  Print["OK  ", tag], Print["FAIL ", tag, " -> ", N[expr, 30]]; Exit[1]];
expectClose[a_, b_, tol_:10^-11, tag_] := If[Abs[N[a, 30] - N[b, 30]] <= tol,
  Print["OK  ", tag], Print["FAIL ", tag, " -> ", N[a, 30], " vs ", b]; Exit[1]];
```
(Strip any `ConditionalExpression[0, ...]` wrapper before testing zero, per the project Mathematica idioms; prefer `1/x == 0` style pole tests over `=!= Infinity` if a pole check ever arises — it does not here.)

**Anti-transliteration guard (MANDATORY — this must NOT be a line-by-line port of the `.py`):**
- M1 bridge: compute `\Xi_1` by forming `P0A = N0A/D0A` with `N0A=N0+eps lam N01`, `D0A=D0+eps lam D01`, then `Normal[Series[P0A, {eps,0,1}]]` and read off the `O(eps)` coefficient divided by `lam P0`; do NOT replicate the SymPy `diff(...).subs(eps,0)/lam` choreography. Confirm it equals `N01/N0 - D01/D0` via `Simplify[... == 0]` / `Together`.
- Nullspaces (M5, M6): extract the coefficient rows with `CoefficientArrays[expr, mixedVars]` (or `Coefficient[expr, v]` per variable assembled into a matrix) on the SAME mixed-restricted symbolic `K1`/`H_even` and `D01/D21/D41` expressions, then use Mathematica's native `NullSpace` and `MatrixRank`. Do NOT re-implement `coeff_vector` step-by-step.
- Projector norms (M5, M6): build the orthogonal projector via `Orthogonalize[nullBasis]` (Gram–Schmidt) → `Transpose[Q].Q`-style projector, or equivalently `B.PseudoInverse[Transpose[B].B].Transpose[B]` with `B` the matrix whose COLUMNS are the null basis; then `sigma = Sqrt[XiCoeff . proj . XiCoeff]`. Do NOT mirror SymPy's `QRdecomposition` call sequence.
- Build the one-pole/hidden-even bundle (B/Z/N/P primitives) and its `eps`-dressed copy from the physical definitions, differentiate with `D[expr, eps] /. eps -> 0` (a derivative is fine — it is a different surface form than Series), and evaluate at the sample point. Independent assembly, not a transcription of the `.py` variable names.

Because the null bases / `\Xi_1(w_i)` / `\Xi_1(t_i)` are returned only up to basis choice and scaling, do NOT hard-compare the raw basis vectors entry-by-entry against the `.py` `expected_*` arrays. Instead verify the BASIS-INVARIANT facts: ranks, nullities, the constraint residuals (M4/M6c), and the two corridor norms `\sigma_{\rm even}` / `\sigma_{\rm transfer}` (which ARE basis-invariant). This both protects against transliteration AND is the correct mathematical test.

**Claim manifest** (each must be independently verified in the `.wl`):

- **M1 — load-defect bridge.** With `D_A=D_0+\eps\lambda D_{01}`, `N_A=N_0+\eps\lambda N_{01}`, `P=N/D`, the first-order load defect satisfies
  `\Xi_1 = (1/\lambda)\,\partial_\eps(N_A/D_A)|_{0}\,/\,(N_0/D_0) = N_{01}/N_0 - D_{01}/D_0 =: \Xi_{\rm load}.`
  Verify `Simplify[\Xi_1 - (N01/N0 - D01/D0)] == 0`.

- **M2 — even-gate closed forms on the compensation surface.** With `u_2=-D_2/D_0`, gates `K_1=D_{21}+D_{01}/9`, `H_{\rm even}=D_{41}-\tfrac23 D_{21}-D_{01}/27`, and the compensation surface `D_{21}=-u_2 D_{01}`, `D_{41}=(D_4/D_0)D_{01}`:
  `K_1|_{\rm comp} = (1/9 - u_2)\,D_{01}`,  `H_{\rm even}|_{\rm comp} = (D_4/D_0 + \tfrac23 u_2 - 1/27)\,D_{01}.`
  Using `D_4/D_0=u_2^2-u_4` and the one-pole identity `u_4=4u_2^2`, also:
  `H_{\rm even}|_{\rm one-pole} = (-3u_2^2 + \tfrac23 u_2 - 1/27)\,D_{01}.`
  Verify all three as `Simplify[lhs - rhs] == 0`.

- **M3 — explicit compatibility base data.** With sample primitives
  `(\kappa,\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M) = (2\sqrt2/\pi, 1/2, 3/10, 2/5, 1/4, 1, 7/5, 2, 1)`
  and `K` fixed by the compatibility condition `D_0 = 3(1+B_2+Z_2)^2/(B_4+Z_4)` with `K=B_0+Z_0+D_0`, verify (tol `1e-12` on `N[...,16]`):
  `D_0\approx 24.2373099886223`, `D_2\approx -1.18562046858190`, `D_4\approx -0.173991572849491`,
  `u_2\approx 0.0489171640391802`, `u_4\approx 0.00957155575054425`, `D_4/D_0\approx -0.00717866681290820`,
  `P_0=N_0/D_0\approx 0.00206979231806289`, and the one-pole identity `u_4-4u_2^2=0` exactly.
  (Bundle definitions — reproduce from the physical premises: `C=\kappa\lambda_B`, `G_U=\lambda_U`, `G_W=\kappa\lambda_W`, `R=\kappa\lambda_R`, `\Delta=\Omega_U^2\Omega_W^2-R^2`, `S_2=\Omega_U^2+\Omega_W^2`, `H=G_U^2+G_W^2`, `Q=G_U^2\Omega_W^2+2G_UG_WR+G_W^2\Omega_U^2`, `P=\Omega_U^2 G_W+R G_U`; `B_0=C^2/\varpi^2`, `B_2=C^2/\varpi^4`, `B_4=C^2/\varpi^6`; `Z_0=Q/\Delta`, `Z_2=(QS_2-H\Delta)/\Delta^2`, `Z_4=(Q(S_2^2-\Delta)-S_2H\Delta)/\Delta^3`; `N_0=P^2/\Delta^2`; and `D_0=K-B_0-Z_0`, `D_2=-(M+B_2+Z_2)`, `D_4=-(B_4+Z_4)`.)

- **M4 — concrete even-gate coefficients on the branch.** `K_1`-coefficient of `D_{01}` `= 1/9 - u_2 \approx 0.0621939470719309`; `H_{\rm even}`-coefficient `= D_4/D_0 + \tfrac23 u_2 - 1/27 \approx -0.0116042611571584` (tol `1e-12`).

- **M5 — strict mixed-sector even-gate corridor.** Restrict the primitive compiler to the mixed family `(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W})` (set the four non-mixed slopes `x_K,x_M,x_{\lambda_B},x_{\varpi}=0`). Form the `2x5` matrix whose rows are the mixed-variable coefficients of `K_1` and `H_{\rm even}`. Verify `MatrixRank == 2` and `dim NullSpace == 3`. Compute the same-charge functional `\Xi_1` coefficient row on the mixed vars and the ambient-Euclidean orthogonal projector onto the nullspace; verify the induced norm
  `\sigma_{\rm even} = \|\Pi_{\rm even}\Xi_1\|_2 \approx 2.67386816837173` (tol `1e-11`).

- **M6 — pure-transfer subcorridor.** Form the `3x5` matrix whose rows are the mixed-variable coefficients of `D_{01}`, `D_{21}`, `D_{41}`. Verify (M6a) `MatrixRank == 3` and `dim NullSpace == 2`; (M6c) on each null vector `t`, `D_{01}\!\cdot\! t = D_{21}\!\cdot\! t = D_{41}\!\cdot\! t = 0` (residual `< 1e-12`); and the induced same-charge norm
  `\sigma_{\rm transfer} = \|\Pi_{\rm transfer}\Xi_1\|_2 \approx 2.31561904386057` (tol `1e-11`).

- **M7 — transported ceiling budgets.** With carried base budgets (from the upstream same-charge ceiling stage — same literals the SymPy script carries)
  `0.367930328492646`, `0.737619063660757`, `2.94889585703134`, `4.63505472371892`,
  verify `budget/\sigma_{\rm even}` reproduces `[0.137602269567650, 0.275862165676603, 1.10285760977778, 1.73346419189450]` and `budget/\sigma_{\rm transfer}` reproduces `[0.158890698998242, 0.318540765855427, 1.27348056877049, 2.00164821411704]` (tol `1e-11`).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 226` and confirms the new `.wl` appears at the Target path AND exits 0 with all checks passing, reproducing the ranks/nullities and both corridor norms above.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl`
- summary: Created the standalone Mathematica audit that independently verifies M1-M7, including native nullspace/projector corridor norms and transported budgets.
- deviation: none

## Authorized notes renumber (USER-AUTHORIZED 2026-06-02)

The user authorized renumbering stale VII.1 notes-prose stage labels to canonical. The audit logged a notes-side renumbering drift here (notes reference pre-renumber stage numbers; the card/appendix/scripts use canonical). Notes-only cleanup in THIS fix loop (Codex applies notes/*.md; Claude reviews):
- In `notes/stages/moving_throat_pde_stage226_..._sympy_audit.md`, renumber every stale stage-number reference to match the canonical numbering used in this stage's SymPy script comments + the paper card (self-reference → Stage 226; cited upstream stages → the numbers the .py comments use). Math/content unchanged.
- Do NOT touch scripts, paper.tex, or appendix. Acceptance: notes stage labels match the .py comments + card. Append `## Applied: notes-renumber`.

## Applied: notes-renumber

- files_changed:
  - `notes/stages/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.md`
- summary: Renumbered stale notes-stage labels from 240-243 to canonical 223-226 references while leaving the mathematical content unchanged.
- deviation: none
