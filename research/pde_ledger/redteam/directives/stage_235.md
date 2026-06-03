---
unit_id: 235
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T22:13:40-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 235

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose document — the red-team only modifies scripts.

After editing, RUN the new script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.wl`

**Issue:** Stage 235 has a SymPy script but no Mathematica second engine, and the unit is neither status-only nor a checkpoint carve-out. Every Stage-235 claim is finite linear algebra over a 2×2 matrix with one symbolic parameter `c_eta`, which Mathematica can verify independently. Write a NEW Mathematica audit script that re-derives the Stage-235 results from the physical premises using native Mathematica primitives. This must be an INDEPENDENT derivation, NOT a transliteration of the `.py`: do not mirror the SymPy variable names/choreography line-by-line. Use a different decomposition where natural (e.g. native `MatrixPower`, `IdentityMatrix`, `Solve`/`Reduce` for the lock conditions, `Simplify` of explicit inner products for the norm). You design the route; the manifest below states only WHAT must be verified.

Use a hard-failing check helper (e.g. `expectZero[expr_] := If[Simplify[expr] =!= 0, Print["FAIL: ", expr]; Exit[1]]`, and a matrix variant that fails on any nonzero entry). Per project idiom, strip `ConditionalExpression[0, ...]` from any `Solve`/`Reduce` output before comparison, and test poles via `1/expr == 0` rather than `=!= Infinity`. Print each verified result and finish with an all-passed line; exit 0 only if every check passes.

**Required change:**
Create the `.wl` at the Target path. Declare `c_eta` as a positive real (it equals `eps/(1-eps)` with `0 < eps < 1`); `R1, E1, q_nt, q_eta, L` real (`L` positive). Verify the claim manifest items M1–M6. Then run it and iterate to a clean exit-0.

**Claim manifest** (the new `.wl` must independently verify each; symbolic forms are the paper's stated targets — eqs `app-part07-mrm`, `app-part07-direct-projectors`, `app-part07-rm-codim-two`, `app-part07-static-blind-line`, notes §1–§5):

- **M1 — involutive compiler.** With `M_rm = {{-1, -c_eta}, {0, 1}}`: verify `M_rm . {R1, E1} == {-R1 - c_eta*E1, E1} == {Xi1, E1}` where `Xi1 = -R1 - c_eta*E1`, AND `MatrixPower[M_rm, 2] == IdentityMatrix[2]`, AND the inverse map `M_rm . {q_nt, q_eta} == {-q_nt - c_eta*q_eta, q_eta}`. (`Det[M_rm] == -1`.)

- **M2 — direct-space projectors.** Build them by similarity, `P_nt = M_rm . DiagonalMatrix[{1,0}] . M_rm` and `P_eta = M_rm . DiagonalMatrix[{0,1}] . M_rm` (valid since `M_rm = M_rm^{-1}`), and verify they equal the paper closed forms `P_nt == {{1, c_eta}, {0, 0}}`, `P_eta == {{0, -c_eta}, {0, 1}}`. To make the route genuinely independent, ALSO confirm the same two matrices satisfy the defining projector properties directly: `P_nt.P_nt == P_nt`, `P_eta.P_eta == P_eta`, `P_nt.P_eta == 0`, `P_eta.P_nt == 0`, `P_nt + P_eta == IdentityMatrix[2]`.

- **M3 — direct-space decomposition.** With `x_rm = {R1, E1}`: verify `P_nt . x_rm == {R1 + c_eta*E1, 0} == {-Xi1, 0}`, `P_eta . x_rm == {-c_eta*E1, E1}`, and `x_rm == P_nt.x_rm + P_eta.x_rm`.

- **M4 — codimension-two orbit lock.** Verify the equivalence `{q_nt == 0, q_eta == 0} <=> {R1 == 0, E1 == 0}` on the rigid-mouth slice via the invertibility of `M_rm` (`Det != 0`), AND independently that `Solve[{Xi1 == 0, E1 == 0}, {R1, E1}]` and `Solve[{Xi1 == 0, R1 == 0}, {R1, E1}]` each yield the unique solution `{R1 -> 0, E1 -> 0}` (strip any `ConditionalExpression`). The static strip `Xi1 == 0` is codimension one; the lock point is codimension two.

- **M5 — static-blind dressing line and exact norm.** On the line `x_blind = {-c_eta*q_eta, q_eta}`: verify `Xi1` evaluated on it is `0` (i.e. `-x_blind[[1]] - c_eta*x_blind[[2]] == 0`), and `x_blind . x_blind == (1 + c_eta^2)*q_eta^2`. Then with `q_eta = L/Sqrt[1 + c_eta^2]`, verify the squared norm equals `L^2` exactly (so `Xi1 = 0` admits points of arbitrary direct-space size `L`).

- **M6 — correction compilers.** With `q = {q_nt, q_eta}` and `x_q = M_rm . q`: verify `Delta_x_static = -M_rm.DiagonalMatrix[{1,0}].q == {q_nt, 0}`, `Delta_x_orbit = -x_q == {q_nt + c_eta*q_eta, -q_eta}`, the additive relation `Delta_x_orbit == Delta_x_static + {c_eta*q_eta, -q_eta}`, and the full-lock identity `x_q + Delta_x_orbit == 0`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 235` and confirm the new `.wl` appears, exercises M1–M6, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.wl`
- summary: Created the Stage 235 Mathematica audit script verifying M1-M6 with hard-failing scalar and matrix checks.
- deviation: none
