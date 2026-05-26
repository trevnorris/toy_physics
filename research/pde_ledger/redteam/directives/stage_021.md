---
unit_id: 021
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 021

The only finding for this unit is a `paper_misalignment`. **Codex must not edit anything for this unit until the user picks a direction below.** Do not touch paper.tex, notes/, or the scripts. The orchestrator is holding pending user resolution.

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_021.tex:77-81` quote:
  > Stage~021 exports the reduced one-port self-energy \eqref{eq:app-stage021-self-energy}, the transfer factor \eqref{eq:app-stage021-transfer-factor}, and the wall-level outgoing quadrupole coefficient \eqref{eq:app-stage021-wall-odd}.
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_021.tex:71-75` quote:
  > δ D_2^{odd}(ω) = -i N_2(0) a^5/(27 c_s^5) ω^5 + O(ω^7)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py:238-244` quote:
  > Gamma_port = sp.symbols("Gamma_port", positive=True, real=True)
  > Dcorr = sp.simplify(-I * Gamma_port * omega**5 * N0)
  > print("If Pi_out = + i Gamma_port omega^5 + O(omega^7), then")
  > print("delta D_wall^(odd) =")
  > sp.pprint(Dcorr)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:135, 140` quote:
  > dCorr = FullSimplify[-I gammaPort omega^5 n0, Assumptions -> $Assumptions];
  > Print["delta D_wall^(odd) = ", fmt[dCorr]];

`Dcorr`/`dCorr` is printed but never asserted, and `Gamma_port`/`gammaPort` is never specialized to the paper's value `a^5/(27 c_s^5)`. The paper's third Output deliverable (the composed wall-level odd quadrupole coefficient) therefore has no script-side check.

## Resolve before fix_loop

The paper's Output paragraph lists three deliverables; the scripts assert two of them and only narrate the third in prose. Which direction should we take?

Possible directions (the user picks one):

- (a) **Script is incomplete — add the composed assertion.** Add a new check in both the sympy and mathematica scripts that substitutes `Gamma_port → a^5/(27 c_s^5)` into the wall-operator-convention odd term and asserts the result equals `-i * (Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² − R²)² * a^5/(27 c_s^5) * ω^5`. Concretely for sympy (insert after Section III's `Dcorr` is defined, or as a new Section VI that imports the constants):
  ```python
  expect_zero(
      "delta D_2^(odd) composed (paper eq:app-stage021-wall-odd)",
      Dcorr.subs(Gamma_port, a**5 / (27 * c_s**5))
      - (-I * (OA**2 * gW + R * gA)**2 / (OA**2 * OW**2 - R**2)**2
         * a**5 / (27 * c_s**5) * omega**5),
  )
  ```
  Note: `a` and `c_s` are currently declared inside `outgoing_l2_fingerprint_audit()` (sympy line 258), not in Section III. The cleanest restructure is a new Section VI that re-declares the symbols and re-imports the formulas. Mirror in the .wl at end of §IV or in a new §VI.

- (b) **Paper is overstating — trim the Output paragraph.** Remove the third item from the Output paragraph (lines 77–81) so it lists only the reduced one-port self-energy and the transfer factor. Also drop `\eqref{eq:app-stage021-wall-odd}` and the surrounding paragraph "Compact outgoing fingerprint" (lines 60–75) if the wall-level composition is genuinely out of scope for this stage. Likely wrong direction given the notes (§7) treat this as the headline result.

- (c) **Both sides are fine; the composition is intended to be implicit.** Document explicitly in the script's final ledger comment or in the paper card that the composition is asserted by inspection from the two separately-verified pieces, not by an end-to-end `expect_zero`. (This is the weakest direction; the auditor flags it because a future edit could silently break it.)

Recommended: (a). The fix is mechanical and adds <10 lines per script.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. If the user picks (a), a follow-up directive should be issued naming the exact insertion points and assertion shapes for both scripts.

## Applied: F1 (iter2)

- files_changed:
  - `scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
- summary: Reworked the composed odd wall-operator assertion so it substitutes the explicit Section III `N(0)` closed form and the Section IV `Gamma5_port = a^5/(27 c_s^5)` coefficient.
- deviation: none
