---
unit_id: 222
batch: VII.1
created_at: 2026-06-02T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T20:29:17Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 222

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F1) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If F2's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F2` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

After editing, RUN the affected script (`math -script <path>` for the new Mathematica file) and iterate until it exits 0 with all in-file checks passing.

---

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper/notes side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md:395` quote: `| 0.2 | 0.00594740531769 | 2.82723442158450 | 2.04402272302752 | 213.483858657863 |`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.py:253` quote: `(0.2, 0.005947405317693074, 2.8272344215844973, 2.04402272302752, 145.4838586578641),`
- saved output line 56 confirms the script value: `... 2.04402272302752, 145.4838586578641)`

## Resolve before fix_loop

In the lambda_W scan table, the lambda_W=0.2 upper-wall residue/linewidth figure R_{Q,*} is **213.483858657863** in the notes but **145.4838586578641** in the script and saved output. All other columns of that row (P0, D0, omega_*) match, and every other scan row (0.4/0.6/0.8/1.0) matches between notes and script — only this one R_{Q,*} cell disagrees. Which value is correct?

Possible directions (the user picks one):
- (a) Script is correct (145.483858657864) → update notes line 395 R_{Q,*} cell to 145.483858657864 (paper-side edit; Codex applies only after this follow-up authorization). No script change. This is the likely case: the script recomputes R_{Q,*} symbolically from F(y) and N(omega), and all four 0.4-row poles already reproduce the notes exactly, so a notes-side typo (e.g. 213 vs 145) is the probable cause.
- (b) Notes value (213.483858657863) is correct → the script's symbolic R_{Q,*} computation for the 0.2 slice must be re-examined and the hardcoded expected_scan literal at line 253 corrected to match; re-run sympy.
- (c) Both trace to a different intended definition of N at the pole → flag for deeper review.

F1 is RESOLVED by the user (see the ## RESOLVED block at the end) — apply the authorized notes edit.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md`
- summary: Corrected the lambda_W=0.2 upper-wall R_{Q,*} notes table cell from 213.483858657863 to 145.483858657863.
- deviation: none

---

## F2 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_mathematica_audit.wl`

**Issue:** No Mathematica engine exists for unit 222. The dual-engine rule requires a `.wl` wherever Mathematica CAN independently verify the stage; every deliverable here is natively within Mathematica's reach (trig integral, polynomial algebra, symbolic ratio simplification, quartic root-finding, numeric substitution). Create an INDEPENDENT re-derivation — NOT a transliteration of the SymPy script.

**Anti-transliteration guard (MANDATORY — do not mirror the .py choreography):**
- Build `D(omega)` symbolically, then obtain the quartic via `Together`/`Numerator` and read the degree off `Exponent[..., y]` and the coefficients off `CoefficientList`. Do NOT hand-expand `F(y)` term-by-term the way the SymPy script does at its lines 92-95.
- Verify the quartic identity by `FullSimplify[Numerator[Together[D]] / ((varpi^2 - y) Delta_y) - Together-form] == 0` form, i.e. compare against the built-up `D`, not against a pre-typed polynomial.
- Find the sample-slice poles with `NSolve[F[y] == 0, y]` (work in `y = omega^2`, take positive real `y`, then `Sqrt`), `WorkingPrecision -> 30`. Do NOT reuse the SymPy `nroots` literal list as input.
- Compute `R_{Q,*}` from `27 c_s^5 / (a^5 omega^5 N[omega])` with `N[omega] = P[omega]^2 / Delta[omega]^2`, `P[omega] = A[omega] G_W + R G_U` — derive `N` from these primitives, not from a stored value.

**Claim manifest** — the new `.wl` must independently verify each of:
- **M1** — overlap constant: with `u0 = 1/Sqrt[L]`, `f0 = Sqrt[2/L] Sin[Pi s/(2 L)]`, `Integrate[u0 f0, {s, 0, L}]` simplifies to `2 Sqrt[2]/Pi` (assume `L > 0`).
- **M2** — quartic pole polynomial: with `C = kap lam_B`, `G_U = lam_U`, `G_W = kap lam_W`, `R = kap lam_R`, `K_B = K - M w^2 - C^2/(varpi^2 - w^2)`, `A = Omega_U^2 - w^2`, `W = Omega_W^2 - w^2`, `Delta = A W - R^2`, `Q = G_U^2 W + 2 G_U G_W R + G_W^2 A`, `D = K_B - Q/Delta`: the numerator of `Together[D]` equals `F(y)` with `y = w^2` up to the factor `(varpi^2 - w^2) Delta`, i.e. `FullSimplify[D - F[w^2]/((varpi^2 - w^2) Delta)] == 0`, and `Exponent[F[y], y] == 4`. Here `F(y) = ((K - M y)(varpi^2 - y) - C^2)((Omega_U^2 - y)(Omega_W^2 - y) - R^2) - (varpi^2 - y)(G_U^2(Omega_W^2 - y) + 2 G_U G_W R + G_W^2(Omega_U^2 - y))`.
- **M3** — residue/linewidth cancellation: with `Gamma5 = a^5/(27 c_s^5)`, `Aqq = 1/D0p`, `gamma = Gamma5 omega_*^5 N_* / D0p`, the ratio `Aqq/gamma` `FullSimplify`s to `27 c_s^5/(a^5 omega_*^5 N_*)` (the symbolic pole-derivative `D0p` must cancel). Also confirm the low-loss survival threshold `2 DeltaVreq (1+eta^2)/eta x^6` and peak `(1/2) RQ eta/(1+eta^2)/x^6` reproduce the stage forms.
- **M4** — static sample-slice data on `(lam_B,lam_U,lam_W,lam_R,Omega_U,Omega_W,varpi,K,M) = (1/2,3/10,2/5,1/4,1,7/5,2,3,1)`, `kap = 2 Sqrt[2]/Pi`, `a=c_s=1`: `Delta0 ≈ 1.9093394081788311`, `D0 ≈ 2.7635551093312736`, `N0 ≈ 0.05016619802495911`, `P0 ≈ 0.018152776420332848` (agree to 1e-12).
- **M5** — sample-slice pole census and residue figures: the four positive `omega` roots of `F(y)=0` are `≈ 0.9382727417467537, 1.3914108765380409, 1.7204537104800286, 2.045399487836587`, and the corresponding `R_{Q,*}` are `≈ 18.7069287828307, 0.380740659074003, 16.0250330226177, 32.0025481088465` (agree to 1e-10). Uncoupled wall/BdG roots `≈ 1.6814318259147836, 2.0427400751933362`; uncoupled internal U/W roots `≈ 0.9746017237463136, 1.417798109771174`.
- **M6** — illustrative low-loss thresholds at x=1 with `V_known = 1.181909222592`, `epsilon = 0.1`, `DeltaVreq = V_known - epsilon`: `2 DeltaVreq (1+eta^2)/eta` gives `≈ 21.854566296358396` at `eta=0.1` and `≈ 7.8618736841685335` at `eta=0.3`. AND the lambda_W scan `lambda_W ∈ {1/5,2/5,3/5,4/5,1}` reproduces P0 monotonically increasing and the upper-wall R_{Q,*} monotonically decreasing; report the lambda_W=0.2 upper-wall `R_{Q,*}` value explicitly (this is the F1-disputed cell: the independent engine determines whether it is ≈145.48 or ≈213.48).

For each M1-M6, use `If[<check fails>, Print[...]; Exit[1]]` (or equivalent `expectZero`/`expectClose`) so the script exits nonzero on any mismatch. Strip `ConditionalExpression[0, ...]` wrappers from any Solve/FullSimplify result before the zero test, and use `1/expr == 0` rather than `=!= Infinity` for any pole check.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 222`, confirms the new `.wl` at the Target path exists, exits 0, and its M4/M5/M6 numerics agree with the SymPy figures (and reports the lambda_W=0.2 R_{Q,*} value for F1 cross-reference).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_mathematica_audit.wl`
- summary: Added an independent Mathematica audit that verifies M1-M6, including the native quartic reconstruction, sample pole census, residue values, thresholds, and lambda_W scan.
- deviation: none

## RESOLVED — F1 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: correct the notes to match the verified SymPy script (script is authoritative; re-derived symbolically in-script and committed in `scripts/output/...stage222..._sympy_audit.txt`). Notes-only; Codex applies, Claude reviews. Do NOT change the script, paper.tex, or appendix.
- In `notes/stages/moving_throat_pde_stage222_..._sympy_audit.md`, the λ_W=0.2 upper-wall R_{Q,*} table cell reads `213.483858657863`; correct it to `145.483858657863` (a +68 integer-part slip; fractional tail unchanged).
- Acceptance: `213.483858657863` no longer appears in the notes; `145.483858657863` is present in the corrected cell. Append `## Applied: F1`.
