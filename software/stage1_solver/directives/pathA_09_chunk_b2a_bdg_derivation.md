# Directive pathA_09 — Chunk B2a: DERIVE the BdG spectrum on the Path-A self-consistent background

**Owner:** Codex (codes + DESIGNS the route: the Path-A background exporter, the `mt15_02` adaptation, and an
independent Python cross-check). Claude reviews afterward (standing transliteration-fidelity audit per math module +
adversarial). **Status:** post-freeze (GATE-A `frozen: YES`, commit `1703f4c`, `candidate_freeze_hash ed3585…`);
B1 extraction DONE (`398ba27`); architecture = decision-12 Reading A (derive, don't inherit).

## Scope of THIS chunk (B2a only)
DERIVE the matter-sector BdG spectrum and its wall coupling on the **Path-A self-consistent closed background**
(the frozen `homogeneous_isotropic_hooke_v1` family, at a neutral reference `τ`), producing the BdG bundle
`{ϖ_j, φ_j, c_j, B0, B2, B4}` that B1's extraction module consumes.

**OUT of scope (do NOT do here):**
- NO Maxwell transfer `{Z_n, N_n}` — that is **B2b** (Spike-2 Green/self-energy chain).
- NO `R_norm(τ)=0` root-find, NO `R_norm`/`R_pole`/`P2`/`P4` assembly — that is **B2c**.
- Do NOT compute `R_norm` or any held-out residual in this chunk. B2a stays target-blind: it emits the bundle and
  its validation only. (Same discipline B1 followed: prove the piece, don't run the physical verdict.)

## Why this chunk exists (the gap decision-12 surfaced)
The earlier M1c effective-closure run did **not** solve the BdG spectrum on a self-consistent background: `mt15_02`
ran the operator around a **smoke / non-equilibrium** profile (stationary residual norm `243.39…`, hardcoded inline
`rho0=0.035+0.18 Exp[…]`), with M1a open-throat geometry (`L=2.0, R_exit=1.65`) and **placeholder** wall
constitutive coefficients. Those `ϖ_j, c_j, B_n` were labelled "MACHINERY, NOT PHYSICS." Reading A requires the
spectrum DERIVED on the genuine Path-A closed solve under the frozen family. The recurring project sin is a
faithful-but-wrong operator (MMS/arbiter cannot catch it) and a tautological / can't-fail gate — every validation
gate below must be able to FAIL, and the operator must be fidelity-audited term-by-term against the canonical source.

## Build (requirements — Codex designs the implementation)
1. **Path-A background exporter (NEW Python).** No exporter exists for the Path-A closed solve
   (`patha_closed_newton.py` performs none; the existing `m1c_background_export.py` exports the *coupled-branch
   continuation*, not the Path-A closed Newton). Create an exporter that:
   - Runs the **closed** Path-A Newton solve (`patha_closed_newton.py` + `coupled_branch.py`) on the frozen
     `homogeneous_isotropic_hooke_v1` provider (`a=1.0`, `L=37/20=1.85`, `w_min=0`) at a **neutral reference
     `τ`** (use `τ=1.0`; this is a convention, NOT a fit — it must NOT be chosen by peeking at any R-target).
   - Confirms the solve is genuinely self-consistent (converged Newton residual ≪ the smoke `243.39`; report the
     actual residual norm) — this background must NOT be the smoke profile.
   - Exports a **versioned JSON bundle** (machine-to-machine handoff; JSON is correct here) holding the converged
     `{ψ0 (psi_real, psi_imag), A0 (a0, ar, aw), R0(w), μ}`, the grid (`r_centers, w_centers, dr, dw, nr, nw`), the
     geometry (`a, L`), and the solver constants (EOS `K`, `ħ`, particle mass, gauge charge, `τ`). Mirror the field
     layout/convention of `m1c_background_export.py` so downstream tooling is consistent. Record a content hash.
2. **Adapt the BdG engine to the Path-A background (Mathematica).** Re-point `mt15_02_bdg_wall_derivation.wls` so it:
   - READS the exported Path-A background bundle instead of the hardcoded inline smoke profile (currently
     `mt15_02` ~lines 227–232) and uses the **frozen geometry** (`a=1, L=1.85`) instead of the hardcoded
     `L=2.0, R_exit=1.65` (~lines 98–100). Resample/interpolate the background onto the BdG mesh as needed.
   - Keeps the BdG operator **FORM identical** to the canonical operator (kinetic with `l(l+1)=6` angular,
     confinement, quintic enthalpy `h=(5K/4)ρ⁴` and its linearization `dh/dρ=5Kρ³`, `−μ` diagonal, gauge `qA0`,
     wall drive `δV_conf=−4 V_radial r⁴/R0⁵ η`). Only the **background fields + geometry** change, not the operator
     terms. Source of truth = `notes/moving_throat_pde_program_compact.md` (BdG L1406–1428, EOS L578–582, wall
     drive L1080–1085 / L1424–1428).
   - Solves the BdG eigensystem → positive real modes `ϖ_j` and profiles `φ_j`; symplectically normalizes;
     reports per-mode eigen-residual `‖L_BdG v_j − ϖ_j v_j‖/‖v_j‖`.
3. **Wall coupling + moments (reuse the validated B1 χ + the v2_22a overlap algebra).**
   - Use the **B1-derived** wall mode `χ(w)` (from `patha_extraction.solve_l2_wall_eigenproblem` on this SAME
     background; `χᵀWχ=1`, `∫χ>0`) — NOT a re-posited wall shape. B1 already validated `χ` against the analytic
     `sin(πw/(2L))` oracle for this family.
   - Form the overlaps `I_{η,φ_j} = ∫ χ · [drive-weighted BdG density response] dw` via the SAME algebra as
     `research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py`; convert to couplings
     `c_j = λ_B · I_{η,φ_j}`; then `B_n = Σ_j c_j² / ϖ_j^{2(n+1)}` for `n=0,2,4` (i.e. `Σ c_j²/ϖ_j^2`,
     `Σ c_j²/ϖ_j^4`, `Σ c_j²/ϖ_j^6`) — the exact form `patha_extraction.bdg_moments` consumes.
   - Emit the bundle in the shape B1 ingests: a `bdg_modes[]` list of
     `{name, lambda_B, varpi, profile, overlap_I_eta_phi}` per mode, exportable straight into the B1 lane contract.
4. **Independent Python cross-check (dual-engine — REQUIRED).** Independently assemble the SAME BdG operator on the
   SAME exported background in Python (Codex designs; may reuse faithful pieces of the existing tangent operator
   `p2_tangent.py`, which already carries the `l(l+1)=6` angular + the `δV_conf` wall drive) and solve it with an
   independent linear-algebra stack (scipy/torch). The two engines must NOT call each other — genuinely independent
   assembly. Compute `ϖ_j, c_j, B_n` in both and compare.

## Independent validation (the heart of B2a — every gate must be able to FAIL)
1. **Self-consistency gate (catches: running on the smoke/wrong background).** Assert the background actually used
   is the converged Path-A closed solve: converged Newton residual ≪ smoke `243.39` (report the number), geometry
   = `a=1, L=1.85`. A run that silently fell back to the smoke profile or M1a geometry must FAIL here.
2. **Dual-engine agreement (catches: an engine-specific assembly/solve bug).** MMA vs independent Python agree on
   the lowest few `ϖ_j`, the `c_j`, and `B0/B2/B4` to a stated tolerance (report max abs/rel diff per quantity).
   Because the two assemble independently, disagreement flags a real bug — not a copy compared to itself.
3. **Eigen-residual (catches: a wrong eigenpair).** Each exported mode satisfies
   `‖L_BdG v_j − ϖ_j v_j‖/‖v_j‖ ≤ tol`. (Checks the eigensolve, NOT the operator's correctness — gate 6 does that.)
4. **Grid convergence (catches: an unconverged spectrum/moments).** Refine the BdG mesh ≥2 levels; show `ϖ_j` and
   `B_n` settle (report the values + observed trend/order) and pick a resolution where `B_n` is converged. A
   resolution-dependent `B_n` must be caught here, not silently shipped.
5. **Structural sanity (catches: a coding bug only — explicitly NOT a physics gate).** `B0·B4 − B2² ≥ 0`
   (Cauchy–Schwarz for any positive `(c,ϖ)`). State that this catches sign/indexing bugs, not wrong physics.
6. **τ-sensitivity (catches: a background NOT re-solved per τ, i.e. a frozen/stale bundle).** Derive the bundle at
   ≥2 reference values (`τ=1` and one other, e.g. `τ=2`); show `{ϖ_j, B_n}` genuinely MOVE with τ through the
   re-solved background (decision-12 "confirm the τ-scaling early"). Equal bundles across τ would mean the
   background isn't actually re-solved — this gate must catch that. (This is machinery confirmation for B2c's sweep,
   NOT a calibration — no `R_norm` is computed.)
7. **No can't-fail gates.** For EVERY gate above, state in the report the specific wrong answer it would catch. No
   verbatim-copy "independent" recompute, no number compared to itself, no tolerance loose enough to pass a wrong
   operator.

## Acceptance criteria (Codex iterates until ALL pass, exit 0)
1. The Path-A background exporter runs the **closed** solve on the frozen family at `τ=1`, confirms a small
   converged residual (reported), and writes the versioned JSON bundle (content hash recorded).
2. The adapted `mt15_02` reads that bundle (not smoke), uses frozen geometry, solves the BdG eigensystem; all
   exported `ϖ_j>0` with eigen-residuals ≤ tol.
3. Dual-engine (MMA vs independent Python) agreement on `ϖ_j, c_j, B_n` to the stated tolerance (report diffs).
4. Grid-convergence of `ϖ_j` and `B_n` demonstrated (≥2 levels; values + trend reported).
5. The bundle is emitted in the B1-ingestible shape (`bdg_modes[]` with `{name, lambda_B, varpi, profile,
   overlap_I_eta_phi}`); a focused test confirms `patha_extraction.bdg_moments` reproduces the reported `B_n` from
   it. Full `pytest` for the `patha_*` suite stays green (no regression).
6. NO Maxwell transfer, NO `R_norm`/root-find computed anywhere in this chunk.
7. Constraints honored: every script wrapped `timeout 600` (exit 124 = reformulate, never raise the cap); **≤2
   concurrent `math -script` seats**; CPU sparse-direct (GPU off); reasoning-only steps consume 0 seats. Firewall
   untouched — do NOT write under or import from `research/pde_audit/simulation/`, do NOT touch
   `physical_export_permitted`. No `git add` / no commit (orchestrator commits after review).

## Report back
The derived `{ϖ_j, c_j, B0, B2, B4}` at `τ=1`; the converged background residual norm (vs the smoke `243.39`); the
dual-engine max abs/rel diffs per quantity; the grid-convergence table; the τ-sensitivity result (`τ=1` vs `τ=2`);
for EACH validation gate the specific wrong answer it would catch; the files created/modified; and any place the
`mt15_02` operator or the v2_22a overlap algebra had to be adapted (and why) rather than reused verbatim.
