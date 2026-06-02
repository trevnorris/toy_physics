---
unit_id: 214
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-02T11:42:02-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 214

Apply EVERY finding below in order (F1, F2, F3). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 was a `paper_misalignment` (notes arithmetic typo) that the USER HAS NOW RESOLVED with direction (a): the script is correct (150), the notes carry a typo (218). You are EXPLICITLY AUTHORIZED to apply the F1 notes edit specified in F1's `## RESOLVED` block. Codex applies the notes edit; Claude reviews. The SymPy script is correct and must NOT change for F1 (only the F3 tautology replacement touches the `.py`).

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md:435` quote: `\boxed{5\cdot 5\cdot 6 = 218.}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:206-209` quote: `bezout_projected = deg_Crs * deg_Crt * deg_Sr` … `if bezout_projected != 150: raise AssertionError("Unexpected projected Bezout bound")`

## RESOLVED (user direction: (a) — notes typo `218` → `150`; AUTHORIZED)

The user confirmed direction **(a)**: the script's `150` is arithmetically correct (5·5·6 = 150) and the notes' `218` is a typo. Codex is AUTHORIZED to apply this single notes-prose edit; Claude reviews. No script change for F1.

**Authorized notes edit:** in `notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md` (line ~435; content-anchor on the boxed expression), change `\boxed{5\cdot 5\cdot 6 = 218.}` → `\boxed{5\cdot 5\cdot 6 = 150.}`. (5·5·6 = 150 is forced and matches the script's `bezout_projected` assertion.) Do not change any other notes content for F1.

After the notes edit, the scripts are re-run for F2/F3 as usual; the notes edit does not affect any script result.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md`
- summary: Corrected the projected one-chart arithmetic typo from `218` to `150`.
- deviation: none

## F2 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`

**Issue:** Stage 214 is non-status-only and non-checkpoint; the project dual-engine contract requires a Mathematica audit alongside the SymPy one, and `\stagefield{Verification}` records "Mathematica audit: none yet." Every claim of this stage is independently verifiable in Mathematica. Create the `.wl` at the Target path above (note: `.wl` files live under `mathematica/`, not `scripts/`).

**Required change:**
Author a new Mathematica audit script that independently establishes the claim manifest below. It MUST derive each result independently using native Mathematica primitives (`D[]`, `Series`+`Coefficient` or `Exponent`/`MonomialList`/`Total`, `Solve`/`Reduce`, `Resultant`/`GroebnerBasis`, `FindRoot`, `MatrixForm`/matrix ops) via a DIFFERENT decomposition than the `.py` (for example: obtain the derivative laws by differentiating `Phi` directly and `Together`/`Cancel`-ing rather than asserting against a pre-formed `N/(...)` quotient; obtain total degrees via `Exponent`/`MonomialList` rather than `Poly(...).total_degree()`; obtain the elimination eliminants via `Resultant[..., y]` rather than hand-coded `Ms*Lr - Mr*Ls`). A line-by-line port of the SymPy algebra (same variable choreography, same intermediate `M`/`L`/`N` hand-expansions, same function names rewritten in Mathematica syntax) is REJECTED as transliteration. Each manifest check must `Print` its residual and `Exit[1]` (or `Abort[]`) on a nonzero/false result so the script's exit code reflects failure.

**Claim manifest** (the new script must independently verify each):

- **M1 (exact derivative laws / stationary numerators).** With
  `Delta = A + B r + C s + D t + E r^2 + F r s + G r t + H s^2 + I s t + J t^2`,
  `Phi = (k_i + k_j r + k_k s + k_l t + Sqrt[Delta]) / Sqrt[1 + r^2 + s^2 + t^2]`,
  `M_r = (1+r^2+s^2+t^2) k_j - r (k_i + k_j r + k_k s + k_l t)` (and `M_s`, `M_t` by cyclic `r↔s↔t`, `k_j↔k_k↔k_l`),
  `L_r = (1+r^2+s^2+t^2) ∂_r Delta - 2 r Delta` (and `L_s`, `L_t` likewise),
  `N_r = 2 M_r Sqrt[Delta] + L_r` (and `N_s`, `N_t`):
  verify `∂_r Phi == N_r / (2 (1+r^2+s^2+t^2)^(3/2) Sqrt[Delta])`, and the analogues for `s` and `t`. (Residual of LHS−RHS must `FullSimplify`/`Together` to 0.)

- **M2 (lifted polynomial degree ledger).** With `F_r = 2 M_r y + L_r`, `F_s = 2 M_s y + L_s`, `F_t = 2 M_t y + L_t`, `F_Delta = y^2 - Delta`, verify the total degrees (jointly in `r,s,t,y`) are `deg F_r = deg F_s = deg F_t = 3` and `deg F_Delta = 2`.

- **M3 (preferred lifted Bézout bound).** Verify `(deg F_r)(deg F_s)(deg F_t)(deg F_Delta) = 3·3·3·2 = 54`. (This must reproduce the paper card `\stagefield{Output}` value `54` and the appendix `eq:app-part06-four-bezout`.)

- **M4 (projected square-root-free eliminant degrees).** With `C_rs = M_s L_r - M_r L_s`, `C_rt = M_t L_r - M_r L_t`, `C_st = M_t L_s - M_s L_t`, and `S_r = L_r^2 - 4 M_r^2 Delta` (and `S_s`, `S_t` likewise), verify each `C_*` has total degree 5 (quintic) in `(r,s,t)` and each `S_*` has total degree 6 (sextic) in `(r,s,t)`. Do NOT hardcode the projected numeric product — leave M5 to carry the numeric bound so the F1 dispute is not baked in.

- **M5 (projected one-chart Bézout product).** Verify the product of the one-chart eliminant degrees `(deg C_rs)(deg C_rt)(deg S_r) = 5·5·6 = 150`. (This is the script's internally consistent value; it is the subject of the F1 notes-vs-script dispute. Use the computed degree product, not a hardcoded `150` or `218`.)

- **M6 (diagonal-isotropic reduction → gradient-optimal ray).** Under the substitution `A=k_i^2-2 H_0 u, B=2 k_i k_j, C=2 k_i k_k, D=2 k_i k_l, E=k_j^2-2 H_0 u, F=2 k_j k_k, G=2 k_j k_l, H=k_k^2-2 H_0 u, I=2 k_k k_l, J=k_l^2-2 H_0 u`, verify `Delta = (k_i + k_j r + k_k s + k_l t)^2 - 2 H_0 u (1 + r^2 + s^2 + t^2)`, hence `tau = 2 H_0 / (k_rst + Sqrt[k_rst^2 - 2 H_0 u])` with `k_rst = (k_i + k_j r + k_k s + k_l t)/Sqrt[1+r^2+s^2+t^2]`; and verify the gradient-optimal ray `a_grad = (k_i,k_j,k_k,k_l)/Sqrt[k_i^2+k_j^2+k_k^2+k_l^2]` is unit-norm and gives slope `a_grad·(k_i,k_j,k_k,k_l) = Sqrt[k_i^2+k_j^2+k_k^2+k_l^2]`.

- **M7 (full-symmetry → equal-mix barycenter stationary).** Under `k_i=k_j=k_k=k_l=k` and the permutation-symmetric envelope `A=E=H=J=k^2-2 H_0 u_d`, `B=C=D=F=G=I=2 k^2-4 H_0 u_x`, verify `N_r(1,1,1)=N_s(1,1,1)=N_t(1,1,1)=0`, and that the equal-mix barycenter `a_eq=(1,1,1,1)/2` is unit-norm.

(The winner/non-improvement integer ordering of section VI is a pure-logic transitivity model and is left to SymPy; a Mathematica re-derivation of M1-M7 above is the required independent second engine.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 214` (or `math -script <Target>`) and confirm the manifest checks (M1-M7) appear and the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering manifest checks M1-M7 with derivative, degree, resultant, Bezout, and symmetry-reduction verifications.
- deviation: none

## F3 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:179-185`

**Issue:** The six section-III elimination identities are zero by algebraic construction and cannot fail:
`Ms*Fr - Mr*Fs - Crs ≡ 0` (the `2 Mr Ms y` terms cancel and the remainder is the definition of `Crs = Ms*Lr - Mr*Ls`), and `Fr*(Fr - 4 Mr y) + 4 Mr^2 FDelta - Sr ≡ 0` (it reduces to `Sr - Sr` once `FDelta = y^2 - Delta` is substituted). They confirm bookkeeping, not that the eliminants capture the stationary conditions.

**Required change:**
Replace the six `expect_zero("cross-elimination identity ...")` / `expect_zero("square elimination identity ...")` asserts at lines 179-185 with a non-tautological verification that holds the eliminants accountable to the actual stationary geometry. The replacement must establish BOTH of the following for each of the six eliminants `C_rs, C_rt, C_st, S_r, S_s, S_t`:

1. **Non-vacuity:** assert each eliminant is NOT the identically-zero polynomial in `(r,s,t)` (e.g., confirm it has at least one nonzero monomial / nonzero total degree). This guarantees the vanishing in step 2 is meaningful and not a trivial `0 == 0`.

2. **Vanishing on a genuine stationary point:** construct at least one concrete numeric stationary point of the lifted system and confirm each eliminant evaluates to `0` there. Acceptance: choose fixed rational, non-degenerate numeric values for `A..J, k_i, k_j, k_k, k_l, H0` (the diagonal-isotropic substitution already in section IV is a convenient source of a closed-form stationary point — the gradient-optimal ray — but any independently-obtained stationary `(r*,s*,t*,y*)` satisfying `Fr=Fs=Ft=FDelta=0` is acceptable), then `assert` each of the six eliminants evaluates to exactly `0` at that point and that the six lifted polynomials `Fr,Fs,Ft,FDelta` also vanish there (so the point is genuinely stationary, not arbitrary).

State both the non-vacuity and the stationary-vanishing as `expect_zero`/`expect_true` checks so a regression (e.g., a sign error introduced upstream into `M`/`L`/`Delta`) would actually fail them. Do NOT delete the degree checks at lines 201-204 — those remain. Do not change any other section.

**Verification:**
The verifier confirms lines ~179-185 no longer contain the definitional `Ms*Fr - Mr*Fs - Crs`-style asserts, that the script now asserts (i) each eliminant is a nonzero polynomial and (ii) each eliminant and each lifted polynomial vanishes at a constructed stationary point, and that `python3 <path>` still exits 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- summary: Replaced the tautological elimination identities with non-vacuity checks and exact stationary-point vanishing checks for the lifted and projected polynomials.
- deviation: none
