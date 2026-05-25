---
unit_id: 006
batch: I.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 006 (v2 paper-grounded re-audit)

This directive contains only `paper_misalignment` findings that require user resolution. The orchestrator will halt before invoking Codex on this unit until the user has chosen a direction.

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_006.tex:38-43` quote:
  ```
  \nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux}
  =\mu_0\mathbf J_{\rm proj}+\mathbf L_{\rm mix}.
  ```
  followed by "Here `\mathbf L_{\rm mix}` is the vector form of the transverse leakage term in eq:stage005-projected-maxwell-expanded."
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_005.tex:28-35` (upstream definition of the only open-system inhomogeneity, the transverse leakage `−[WZF^{wμ}]_∂ + ∫(∂_w W) Z F^{wμ} dw`). No "gauge-driver" symbol appears anywhere in either card.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:66-69` declares `G0, G1, G2, G3 = sp.Function("Gauge0/1/2/3")(t,x,y,z)` as "gauge-driver" placeholders.
- Same file, lines 135 (Gauss `lhs0 = ... + L0 + G0`), 142–144 (Ampere `lhs_i = ... + L_i + G_i`), 165–167 (Section 4 compact summary: `div D + Leak0 + Gauge0 = mu0 rho_proj` and `curl H − partial_t D + Leak_vec + Gauge_vec = mu0 J_proj`), 193–195 (Section 4 repeated print).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl:104` declares `gauge = Table[Symbol["gauge" <> ToString[i]][t, x, y, z], {i, 0, 3}]`, and line 117 adds `+ gauge[[nu + 1]]` to `projectedInhom[nu]`.

In both scripts the `Gauge_μ` / `gauge[nu]` symbols cancel between LHS and RHS of every assertion they appear in (e.g. SymPy `lhs1 = ... + L1 + G1` versus `amp1_target = ... + L1 + G1`; Mathematica's M2 rearrangement residual `projectedInhom[i] - (...)` has the same `gauge[[i+1]]` on both sides). They are never given a definition or exercised by any concrete check (Section 5 SymPy, M4 Mathematica). They are pure scaffolding.

## Resolve before fix_loop

The stage 006 paper card and the upstream stage 005 card name only `L_mix` (transverse leakage) as the open-system inhomogeneity of the projected Maxwell laws. The SymPy and Mathematica scripts additionally carry symbolic `Gauge_μ` / `gauge` placeholder terms in both Gauss and Ampere equations, advertised in the SymPy docstring (lines 7–12) and Section 4 compact summary as "gauge-driver" contributions. These placeholders are never exercised by any assertion (they cancel trivially between LHS and RHS). Which side is correct?

Possible directions (the user picks one):

- **(a) Paper is correct — drop the `Gauge_μ` placeholders from the scripts.**
  - Edit SymPy file `scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py`:
    - Delete lines 66–69 (the `G0..G3` declarations).
    - Update line 57 comment from "Projected source, leakage, and gauge-driver terms" to "Projected source and leakage terms".
    - Update line 135: remove `+ G0` from `lhs0`.
    - Update lines 142–144: remove `+ G1`, `+ G2`, `+ G3` from `lhs1`, `lhs2`, `lhs3`.
    - Update lines 150–152: remove `+ G1`, `+ G2`, `+ G3` from `amp1_target`, `amp2_target`, `amp3_target`.
    - Update lines 139 and 166–167: change print strings to drop `+ Gauge0` / `+ Gauge_vec`.
    - Update lines 193–195 (Section 4 compact summary) likewise.
    - Update the docstring (lines 7–12) to remove "gauge-driver" mention.
    - Update line 182 ("The Gauge_mu terms are the projected gauge-driver contributions.") — delete this print line.
  - Edit Mathematica file `mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl`:
    - Delete line 104 (`gauge = ...`).
    - Update line 117: remove `+ gauge[[nu + 1]]` from `projectedInhom[nu]`.
    - Update line 119 (`gaussRearranged = projectedInhom[0] - (divergence3[Dflux] + leak[[1]] + gauge[[1]]);`): remove `+ gauge[[1]]`.
    - Update line 121–122 (`ampereRearranged` table): remove `+ gauge[[i + 1]]`.
  - Re-run both engines via the verifier; outputs should still pass.

- **(b) Scripts are correct — paper card needs to introduce the gauge-driver term.**
  - Edit `paper/stages/stage_006.tex`: introduce `\mathbf G` (or analogous) in the inhomogeneous equations, with a brief paragraph explaining where it originates (e.g., from gauge-fixing contributions to the parent Maxwell action that survive projection).
  - Update `paper/stages/stage_005.tex` if the gauge-driver term needs to appear in the covariant law there as well.
  - No script change.

- **(c) The gauge-driver placeholder is intended as a future hook (not active in stage 006) but should be exercised by at least one concrete assertion to justify its presence in the script.**
  - Either: paper introduces the symbol AND the scripts add a concrete check (a specific bulk-coordinate gauge-fixing contribution that survives projection in Section 5 / M4), OR: scripts drop the symbol per (a) and the paper introduces it later when a concrete realization exists.
  - The orchestrator should treat (c) as a deferred (a): drop now, add back when the concrete derivation is ready.

Recommended: **(a)**, since the scripts never exercise the `Gauge_μ` placeholders and the paper card does not derive them. Direction (a) is the minimum-edit, minimum-risk choice. If the user picks (a), Codex can apply the listed edits in a follow-up directive (after explicit user authorization).

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
