# Directive pathA_10 — Chunk B2a remediation (post-adversarial-review fixes)

**Owner:** Codex (codes + runs + iterates). Claude reviews afterward. **Status:** B2a build DONE (report
`reports/patha_b2a_bdg_derivation_report.md`); two fidelity audits returned FAITHFUL (no operator divergence);
adversarial review found NO fatal flaw but THREE real gate-quality / convergence problems. This directive closes
them before B2a is committed and before B2c consumes the `B_n`. These are math/gate-quality fixes only — they do NOT
change the conceptual nature of B2a (still: derive the BdG bundle on the self-consistent background, target-blind).

Do NOT re-architect anything that passed. Keep all FAITHFUL operator code. No `git add`/commit (orchestrator
commits). Constraints unchanged: `timeout 600` per script (exit 124 = reformulate), ≤2 `math -script` seats, CPU,
firewall (no writes/imports under `research/pde_audit/simulation/`, no touching `physical_export_permitted`), NO
`R_norm`/`R_pole`/`P2`/`P4`/root-find/Maxwell in this chunk.

## R1 (PRIMARY — modal-truncation convergence; the substantive gap)
The eigensolve returns ~100 positive modes; B2a exported only 3 and summed `B_n = Σ c_j²/ϖ_j^{2(n+1)}` over those 3.
The grid-convergence gate refined the spatial mesh but held the mode count fixed, so the modal-truncation error in
`B_n` is unmeasured. Fix:
- **Sweep `B_n` over the number of modes included in the moment sum** `M ∈ {3, 5, 8, 15, 30, all-positive}` at the
  final spatial grid (and confirm the picture is stable at one finer spatial grid). Report `B0,B2,B4` vs `M` and the
  relative change between successive `M`. (This is cheap — the modes are already computed; you are only extending the
  SUM, not re-solving.)
- **Choose the exported mode count `K` by a principled tolerance** (e.g. `B_n` stable to ≤ a stated % between `K`
  and all-positive), NOT a hard-wired 3. The bundle's authoritative `B0,B2,B4` must be the **converged** sum (enough
  modes), and `bdg_modes[]` must contain enough modes that `patha_extraction.bdg_moments(bdg_modes)` reproduces the
  authoritative `B_n` to tolerance (so the consumer stays consistent — see R2).
- **Add a `modal_truncation` gate that CAN FAIL**: assert `|B_n(K) − B_n(all)| / |B_n(all)| ≤ tol` for n=0,2,4.
  State the wrong answer it catches (an under-truncated moment sum with a fat high-mode tail).
- **Record the residual modal-truncation error** for each `B_n` in the bundle metadata (R4 consumes it).
- Both engines (MMA + Python) must agree on the converged `B_n` (extend the dual-engine comparison to the converged
  sum, not the 3-mode sum).

## R2 (de-tautologize the consumer gate)
`_b1_moment_check` currently recomputes `bdg_moments` from the bundle's own `bdg_modes` and compares to the bundle's
own stored `B_n` — the same formula on the same numbers (a self-comparison; the "3e-16" is float round-trip). Replace
it with a GENUINE, falsifiable check:
- Assert that `patha_extraction.bdg_moments` applied to the **MMA-engine** modes equals the **Python-engine**
  converged `B_n` (a cross-engine consumer check), to a stated tolerance. This tests that the two independently
  assembled spectra produce the same moments through the B1 consumer — something that can actually fail.
- (If you prefer, keep a separate cheap "bundle shape / key" structural check, but it must be labelled as a shape
  check, not a physics/consistency gate.) State the wrong answer the new gate catches.

## R3 (tighten the dual-engine gate)
The dual-engine pass uses `abs ≤ tol OR rel ≤ tol`; for the tiny-magnitude `B_n` the abs branch passes almost for
free. Change to **AND** (`abs ≤ tol_abs AND rel ≤ tol_rel`) with tolerances appropriate to each quantity's
magnitude (`ϖ`, `c`, `B_n`). Confirm it still passes on the real numbers and state the looser-branch wrong answer it
now rejects.

## R4 (honest convergence reporting + error budget for B2c)
- Stop reporting the successive-difference ratio as an "observed order ≈3.4" — that is not a Richardson order on the
  non-uniform 6→8→10 ladder. Either compute a proper Richardson order on a **uniform** refinement ladder (e.g.
  N, 2N) or simply report the relative changes per refinement WITHOUT the misleading order claim.
- **Carry BOTH error contributions as explicit numbers into the bundle metadata**: (a) the spatial-discretization
  error (finest-vs-previous relative change in `ϖ_j` and `B_n`), and (b) the modal-truncation error from R1. These
  are the `B_n` (and `ϖ_j`) uncertainties B2c's §J error budget will consume. Do not ship a ~2% change as
  "converged" without the error bar attached.

## R5 (honest gate scope + the τ finding — report text)
In `reports/patha_b2a_bdg_derivation_report.md`:
- Re-describe the **dual-engine** gate to its ACTUAL scope: it independently re-derives the BdG matrix assembly +
  eigensolve + overlap quadrature; the **background** and the **wall mode χ** are shared single-engine (Python)
  inputs (this sharing is required by directive pathA_09), and those common-mode pieces were independently validated
  earlier — the closed background in chunks 1b/1c and χ against the analytic oracle in B1. Name those validations.
- Record the **τ-sensitivity physics finding** as a B2c design input (NOT a defect): doubling τ moves the matter
  background sub-percent (τ enters only the wall sector `T_w,T_Ω,U` in the frozen Hooke family), so B2c's τ-leverage
  on `R_norm` is dominated by `K = τκ̂` (exact) and the wall/Maxwell sectors, not by `B_n`. State the measured
  background/`B_n` movement for τ=1→2.

## Acceptance criteria (iterate until all pass, exit 0)
1. Modal-truncation sweep present; exported `B_n` are the converged sum; `bdg_modes[]` reproduces them via
   `bdg_moments`; `modal_truncation` gate present and passing with documented `K` + tolerance.
2. The consumer gate is now a genuine cross-engine check (not a self-comparison) and CAN fail.
3. Dual-engine gate uses AND; still passes; reports per-quantity abs+rel diffs.
4. No "observed order 3.4" claim; spatial + modal error bars recorded in the bundle metadata.
5. Report re-describes the dual-engine + τ gates honestly and records the τ finding for B2c.
6. Full `pytest` for `test_patha*.py` green; no `R_norm`/Maxwell/root-find anywhere; firewall untouched; no commit.

## Report back
The `B_n`-vs-`M` modal sweep table + chosen `K` + truncation error; the new converged `B0,B2,B4` (and how much they
moved vs the old 3-mode values); the cross-engine consumer-check result; the recorded spatial+modal error bars; and
a one-line "what wrong answer it catches" for every gate that changed (modal_truncation, consumer, dual-engine).
List files modified.
