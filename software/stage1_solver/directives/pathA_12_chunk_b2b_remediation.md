# Directive pathA_12 — Chunk B2b remediation (post-adversarial-review fixes)

**Owner:** Codex (codes + DESIGNS the route + runs + iterates). Claude reviews afterward. **Status:** B2b build DONE
(report `reports/patha_b2b_maxwell_transfer_derivation_report.md`); two transliteration-fidelity audits returned the
operator/transfer FORMS **FAITHFUL** to canonical (no sign/factor error; `Γ_port=4a⁵/27c_s⁵` correctly resolves the
decision-05 factor-4 flag; target-blind CLEAN). An adversarial pass + the second fidelity audit found the **validation
layer substantially over-claims and the numbers are not converged**. This directive closes those before B2b is
committed and before B2c consumes `{Z_n,N_n}`.

**User methodology call (2026-06-18):** "Full rigor + true 2nd engine." So convergence must be driven properly AND a
**genuinely independent second discretization** must be built (not a transliteration of the same scheme). Do NOT
silently ship coarse to beat the clock; if a step would exceed `timeout 600`, **reformulate** the math (sparse /
modal / better-conditioned extraction), never raise the cap. If a requirement is genuinely infeasible, HALT and
surface the specific obstruction to the user — do not downgrade unilaterally
([[feedback-never-alter-calibrated-process]], [[feedback-dual-engine-required]]).

Keep all FAITHFUL operator/transfer code (the VSH operator, DtN, `Γ_port`, `Z=gᵀ·inv(K)·g`, `N=−Im Σ_ret/(Γ_port
ω⁵)`, the decision-05 D3 Fréchet forward source). Do NOT re-architect what passed fidelity. No `git add`/commit
(orchestrator commits). Constraints unchanged: `timeout 600` per script (exit 124 = reformulate, never raise), ≤2
`math -script` seats, CPU, firewall (no writes/imports under `research/pde_audit/simulation/`, don't touch
`physical_export_permitted`; `research/pde_audit/{scripts,notes}/` are READ-ONLY), and **target-blind**: NO
`R_norm`/`R_pole`/`P2`/`P4`/`D0`/`P0`/root-find/`mt15_05` preview.

## R1 (PRIMARY — genuine convergence; the substantive gap)
The transfer is **not converged**: the radial-truncation sweep moves `N0` by 45% then 32% per step (a monotone march,
not a settling tail), the sweep only ever *shrinks* the domain (`radial_scale 0.8→0.9→1.0` — the open boundary is
never probed in the converging direction), the grid is a hardcoded 7×7, and `N4/Z4` oscillate across the grid sweep.
The recorded 31.5% `N0` error bar was stamped PASS by a `≤0.75` (75%) tolerance — a can't-fail gate. Fix:
- **Probe the open boundary OUTWARD.** Extend the radial/axial domain in the converging direction (e.g. domain
  factor `≥1.0` and growing: `1.0, 1.5, 2.0, …`) with a genuine outgoing/sponge boundary (see
  `docs/boundary_and_noise_methods.md` — real sponge ≠ absorber; the outgoing-impedance treatment applies to the
  wave-like radiative tangent that carries `N_n`). Show `{Z_n,N_n}` settle as the domain grows.
- **Refine the grid past 7×7** until the per-refinement change *shrinks* toward a stated few-percent tolerance.
  Provide a ≥3-level ladder so a Cauchy/Richardson (decreasing-increment) statement is possible — a single
  finest-vs-previous diff is not a convergence demonstration. If the dense `O((5·nr·nw)³)` solve cannot reach a
  converged grid under `timeout 600`, **reformulate** (CPU sparse-direct and/or the smaller modal/port formulation
  per [[project-gpu-disabled-machine]]) — do NOT cap the grid for speed and ship the coarse number.
- **Re-derive the CONVERGED `{Z0,Z2,Z4,N0,N2,N4}`** and carry the honest residual error per coefficient into the
  bundle `error_budget` (B2c §J consumes it). The authoritative bundle values must be the converged ones.
- **Rewrite the convergence gate to actually fail.** Replace the `≤0.75` rubber stamp with a real shrinking-increment
  / few-percent tolerance across the ≥3-level ladder (grid, outward-domain, ω-window). Add a negative-control TEST
  that feeds a deliberately under-resolved/unconverged sequence and asserts the gate FAILS. State the wrong answer it
  catches (an unconverged transfer with a fat truncation tail).

## R2 (build a GENUINELY independent second engine — the dual-engine claim was false)
The current "dual-engine" MMA worker is a **line-by-line transliteration** of the Python FD engine (same stencils,
same grid, same quadrature, same anchor, same ω-window, same fit; intermediates match to 17 digits). Their 1e-15
agreement certifies transcription only — it cannot catch a faithful-but-wrong *discretization* (both would be wrong
identically). Fix:
- Build a **genuinely different second discretization** for `{Z_n,N_n}` on the same Path-A `A0` background — a
  different scheme, not a re-typed copy. Acceptable forms include: a higher-order or staggered FD stencil, a spectral
  / Chebyshev collocation, or a finite-element / weak-form assembly. The two engines must (a) NOT call each other,
  (b) NOT share the discretization (different stencil/mesh/quadrature/basis), and (c) **both converge to the same
  `{Z_n,N_n}`** as resolution increases (agreement at the converged value is the gate, to a tolerance reflecting the
  combined discretization error — abs AND rel).
- It is fine to keep the existing Python↔MMA transliteration pair additionally as a transcription check, but it must
  be **labelled** as transcription-only; the correctness-certifying dual-engine comparison must be between the two
  genuinely-different discretizations.
- Report honestly what makes the two engines independent (the concrete differences in stencil/mesh/basis), and the
  converged-value agreement. If a genuinely independent second discretization cannot be made to converge under the
  constraints, HALT and surface to the user — do not relabel the transliteration as "independent."

## R3 (replace the tautological basis-invariance gate with a genuine test)
`basis_invariance_check` rotates operator + source + boundary by the SAME dense random orthogonal `Q`, so
`⟨Qs,(QAQᵀ)⁻¹Qs⟩ = ⟨s,A⁻¹s⟩` is an algebraic identity that holds for ANY operator (including a wrong one); the 4e-16
result tests nothing physical, and the `broken=True` control (rotate `A` but not `s`) is a strawman. Fix:
- Make the check a **genuine lane-rebasis / gauge-deflation transform** within the physical port/gauge subspace — the
  kind of basis change a *posited-U/W-port* extraction would FAIL but the basis-invariant `Z=gᵀ·inv(K)·g` passes
  (this is the exact property decision-05 D4 was designed to guarantee). State the wrong answer it catches (a transfer
  that secretly depends on a posited U/W eigenbasis).
- Replace the strawman negative control with a **realistic port-leak** failure mode (e.g. inject a basis-dependent
  U/W-port term and assert the gate FAILS). Add it as a TEST.

## R4 (gate honesty + negative controls + N-channel robustness)
- The `no_cant_fail_gates` gate is hardcoded `passed: True` (it asserts the property it should verify). Remove the
  hardcode; make it (or the test suite) actually demonstrate, per gate, a negative control that fails.
- Add negative-control TESTS for the gates that currently have none: convergence (R1), self_consistency (feed a
  smoke/stale background → must fail), τ-sensitivity (feed a frozen τ=1 bundle as τ=2 → must fail).
- **N-channel robustness:** `−Im Σ_ret` is ~1e-9 of `Re Σ_cons` and `N_n` is extracted through a `cond≈5.5e6` fit;
  `N4/Z4` oscillate. Demonstrate the radiated-power signal is **above** the solve/fit noise floor and that `N0,N2,N4`
  are stable once converged (R1); report the fit conditioning honestly and, if the polynomial-in-ω² extraction is the
  conditioning culprit, use a better-conditioned extraction (orthogonal basis / analytic low-ω expansion of `Σ`).
  A `N_n` that cannot be shown above the noise floor must not be shipped as a derived value.

## R5 (resolve the forward-source provenance + honest report text)
The fidelity audit found the live forward source is the **canonical decision-05 D3 step-8c Fréchet source**
(`δρ=2(ψ_R0δψ_R+ψ_I0δψ_I)`, `δj_i=(ℏ/m)(…)−(q/m)(δA_i0ρ0+A_i0δρ)`), driven by the **B2a BdG-mode response**
(`coupling·profile/(ϖ²−ω²)`) — NOT a B1-χ(w) overlap via the v2_22a adapter. The CODE is faithful to D3; the
mismatch is in the words. Fix the docs to match the code (do NOT rewrite a faithful source):
- Correct the term-map / diagnostics that assert "B1 χ orientation": the source is the D3 Fréchet source over the
  B2a BdG-mode response; state that accurately and remove the unexercised "B1 χ" provenance claim.
- Add a dated correction note to `directives/pathA_11_chunk_b2b_maxwell_transfer.md` req-3: the prescribed
  forward-source construction is the canonical decision-05 D3 Fréchet source over the B2a BdG-mode response (the χ /
  v2_22a wording was imprecise; the BdG bundle already encodes the χ coupling via the B2a overlaps `c_j`).
- **Honest `A0` framing:** the report overstates how much of `A0` enters. State plainly that the B2a *rest* background
  has `A_r0 ≡ A_w0 ≡ 0` (physically expected for a rest defect — no currents → no vector potential), so the transfer
  runs through the scalar (`A_00`) + kinematic (`j_E`) + `q·δρ` channels; the full gauge-current coupling lights up
  only for non-rest (moving/excited) defects, a later held-out regime. This is correct for the calibration point, not
  a defect — say so.

## Acceptance criteria (Codex iterates until ALL pass, exit 0)
1. Convergence driven properly: outward-domain + grid + ω-window ladders (≥3 levels each) with *shrinking* increments
   to a stated few-percent tolerance; authoritative `{Z_n,N_n}` are the converged values; honest per-coefficient
   error bars in the bundle `error_budget`. The convergence gate now FAILS on an unconverged sequence (negative-control
   test present).
2. A genuinely independent second discretization exists (concretely different stencil/mesh/basis, no shared kernel,
   no cross-calls); both engines converge to the same `{Z_n,N_n}` (abs AND rel agreement at the converged value
   reported). Any retained transliteration pair is labelled transcription-only. (Or a documented HALT-to-user if
   infeasible.)
3. Basis-invariance gate is a genuine lane-rebasis/gauge-deflation test with a realistic port-leak negative control
   (test present and failing on the leak).
4. `no_cant_fail_gates` no longer hardcoded; negative-control tests added for convergence, self_consistency, τ.
   N-channel shown above the noise/fit floor with honest conditioning reporting (or a better-conditioned extraction).
5. Forward-source provenance corrected in the term-map/diagnostics + a dated note in pathA_11; honest `A0`/rest-defect
   framing in the report.
6. Full `pytest` for `test_patha*.py` green (no regression); NO `R_norm`/`R_pole`/`P2`/`P4`/`D0`/`P0`/root-find/
   `mt15_05` anywhere; firewall untouched; no commit.

## Report back
The converged `{Z0,Z2,Z4,N0,N2,N4}` at `τ=1` (and how much they moved vs the unconverged 7×7 values); the
outward-domain + grid + ω-window convergence ladders with the shrinking increments + chosen resolution + honest error
bars; the description of the genuinely-independent second discretization and the two engines' converged-value
agreement (abs+rel); the genuine basis-invariance result + its port-leak negative control; the N-channel
noise-floor/conditioning demonstration; the corrected forward-source provenance + honest `A0` framing; for EACH gate
the specific wrong answer it now catches and the negative-control test that proves it can fail; and every file
created/modified.
