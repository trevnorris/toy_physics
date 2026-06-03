---
unit_id: 221
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-02T14:29:39-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 221

Apply each non-`paper_misalignment` finding below (F1, F2, F3) in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F4) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl:154-164`

**Issue:** The two survival-window `expectZero` checks are algebraic round-trips of their own definitions and cannot fail. `survivalLeft` (L155) is defined as exactly `2*UdispLowLossMax` (L154), so `survivalLeft - 2 UdispLowLossMax` (L160) is identically 0. The check at L161-164, `residueRequirement * eta/(1 + eta^2) * Sfam^2 - 2 DeltaVreq`, simply inverts the L156-158 definition of `residueRequirement` and is identically 0. Neither exercises the survival inequality against the line shape derived earlier in the script.

**Required change:**
Replace the two tautological `expectZero` calls (L160 and L161-164) with two substantive checks that route through the line-shape `Udisp`/`reExpected` already defined at L116/L135:

1. Low-loss bound recovered from the line shape at the boundary `r = 1/eta`:
```
expectZero[
  "low-loss |U_disp|_max recovered from line shape at r=1/eta",
  ((-Udisp) /. r -> 1/eta) - UdispLowLossMax
];
```
This is non-trivial: it requires `-Udisp = (1/2) Aabs r/(gamma(1+r^2)) Sfam^2` evaluated at `r=1/eta` to equal `(1/2) Aabs/gamma * eta/(1+eta^2) Sfam^2`; it FAILS if the line shape is wrong.

2. Survival inequality saturated through the line-shape left side (not through `residueRequirement`'s own definition):
```
expectZero[
  "residue requirement saturates survival window via line shape",
  (survivalLeft /. Aabs -> residueRequirement gamma) - 2 DeltaVreq
];
```
Substituting `Aabs -> residueRequirement*gamma` into `survivalLeft = Aabs/gamma * eta/(1+eta^2) Sfam^2` gives `residueRequirement * eta/(1+eta^2) Sfam^2`, which must equal `2 DeltaVreq` only because `residueRequirement = 2 DeltaVreq/Sfam^2 (1+eta^2)/eta`; this routes the saturation through `survivalLeft` (the physical left side) rather than re-stating `residueRequirement`. Keep the two `Print` lines at L166-167.

**Verification command:** verifier runs `redteam exec-mathematica 221`; confirms the two new `expectZero` names reference `Udisp`/`survivalLeft` and the `.wl` exits 0 with both new checks PASS.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
- summary: Replaced the tautological survival-window checks with line-shape-routed checks through `Udisp` and `survivalLeft`.
- deviation: none

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py:136-146`

**Issue:** The SymPy survival-window section computes `Udisp_lowloss_max` (L137-139) and `residue_requirement` (L140-142) but only prints them — no assertion exists for paper deliverable #9 in SymPy. Combined with F1, the headline survival window is unverified in both engines.

**Required change:**
After L142 (and before the `print` block at L144), add two substantive asserts mirroring F1, routing through the existing line-shape symbols `U_disp` (L106) and `re_expected` (L75):

```python
    assert sp.simplify((-U_disp).subs(r, 1 / eta) - Udisp_lowloss_max) == 0
    assert sp.simplify(
        (Aabs / gamma * eta / (1 + eta**2) * Sfam**2).subs(Aabs, residue_requirement * gamma)
        - 2 * DeltaV_req
    ) == 0
```
The first asserts the general line-shape `-U_disp = (1/2) Aabs r/(gamma(1+r^2)) Sfam^2` evaluated at the low-loss boundary `r=1/eta` reproduces `Udisp_lowloss_max` (fails if the line shape is wrong). The second asserts the survival left side saturates `2 DeltaV_req` when `Aabs/gamma` equals `residue_requirement`. Do not alter the probe-only numeric slice (L148-174) or its "probe-only" comment.

**Verification command:** verifier runs `redteam exec-sympy 221`; confirms two new asserts appear after L142 referencing `U_disp`/`Udisp_lowloss_max` and the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- summary: Added the two survival-window assertions that recover the low-loss bound from `U_disp` and saturate the line-shape left side with the residue requirement.
- deviation: none

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl` (Sections I, II, IV)

**Issue:** The `.wl` is a line-by-line port of the `.py`: identical intermediate definitions (`chiBW`, the pre-written `Nfun = (Afun GW + R GU)^2/DeltaPi^2`, the hand-collapsed `chiR = Aabs/(r gamma - I gamma)`, the hand-written factorizations), identical ordering, differing only in syntax. It provides no independent confirmation. For a checkpoint the second engine must reach the load-bearing results by a different route.

**Required change:**
Re-derive the load-bearing results in the `.wl` using native Mathematica primitives and a DIFFERENT decomposition than the `.py`, then keep the existing equality cross-checks. Do not merely rename variables. Apply the claim manifest below. Where a native derivation reproduces a quantity the script later reuses, derive it first and assert it equals the previously hand-written form (that equality is the new independent check).

**Claim manifest** (each must be derived independently in the `.wl`, not transliterated):

- **M1** — derivative-identity residue form. Do NOT pre-write `Nfun = (Afun GW + R GU)^2/DeltaPi^2`. Instead compute `NfunDerived = Together[-D[QPi/DeltaPi, portPi]]` from the raw `QPi`, `DeltaPi`, then verify it is the claimed perfect square via `expectZero["N is a perfect square", NfunDerived - (Afun GW + R GU)^2/DeltaPi^2]` and use `NfunDerived` (not the pre-written form) in the `dD_Pi/dPi + N` check. Symbolic claim: `-∂_Π(Q_Π/Δ_Π) = (A·G_W + R·G_U)²/Δ_Π²`.

- **M2** — Breit-Wigner reduction by series, not by re-factoring. Obtain the pole form independently via `chiSeries = Normal[Series[Numstar/(Fprime delta - portPi Zstar) /. portPi -> I GammaOut, {delta, 0, 0}]]` is NOT appropriate (the leading term is the pole); instead derive `chiPassive` by `Together` of the linearized denominator and verify against the Breit-Wigner target `Astar/(delta - I gammaStar)` exactly as a residue identity: confirm `Residue` of the simple pole, e.g. assert `(delta - I gammaStar) chiPassive /. delta -> I gammaStar` equals `Astar` — i.e. `expectZero["A_* is the residue at delta=I gamma_*", Limit[(delta - I gammaStar) chiPassive, delta -> I gammaStar] - Astar]`. Symbolic claim: the residue of `chi_passive` at `δ = iγ_*` is `A_*`.

- **M3** — line shape from the un-collapsed Breit-Wigner. Do NOT start from the pre-collapsed `chiR = Aabs/(r gamma - I gamma)`. Start from `chiBWgeneric = Aabs/(delta - I gamma)`, take `ComplexExpand[Re[...]]`/`ComplexExpand[Im[...]]` in terms of `delta`, THEN substitute `delta -> r gamma`, and verify the results equal `reExpected`/`imExpected`. Symbolic claim: `Re[A/(δ−iγ)] = Aδ/(δ²+γ²)`, `Im = Aγ/(δ²+γ²)`, reducing to `Ar/(γ(1+r²))`, `A/(γ(1+r²))` at `δ=rγ`.

Sections III, V, VI may retain their current equality checks (they are not the headline independence concern), but if the F1 fix is applied, the survival checks there become the substantive line-shape-routed checks.

**Verification command:** verifier runs `redteam exec-mathematica 221`; confirms (a) `Nfun` is derived via `D[QPi/DeltaPi, portPi]` rather than pre-written, (b) a `Residue`/`Limit`-based check for `A_*` appears, (c) the line shape is obtained from `Aabs/(delta - I gamma)` before the `delta -> r gamma` substitution, and (d) the `.wl` exits 0 with all checks PASS.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
- summary: Re-derived the Mathematica load-bearing checks via the raw derivative identity, native residue extraction, and uncollapsed Breit-Wigner line shape before preserving the equality checks.
- deviation: Used `D[QPi/DeltaPi, portPi]` instead of the directive's leading-minus form because the referenced Stage 220 identity is `D_Pi = K_B - Q_Pi/Delta_Pi` with `\partial_\Pi D_\Pi = -N`.

## Resolve before fix_loop

**Finding F4 — paper_misalignment (notes_contradicts_script).** The notes file (`notes/stages/moving_throat_pde_stage221_..._sympy_audit.md`) names this stage "Stage 238" throughout (title, §3, §9) and attributes the wall derivative identity to "Stage 237," while the card (`paper/stages/stage_221.tex:9`), the Part VII appendix, and both scripts name it Stage 221 and attribute the identity to "Stage 220." The math is identical on both sides; only the stage numbering differs (renumbering-era drift in the notes prose).

Possible directions (the user picks one):
- (a) The card/appendix/scripts numbering (221, Stage-220 identity) is canonical → update the notes prose to use Stage 221 / Stage 220 (a notes-only edit; per the RESOLVED block below, Codex now applies it and Claude reviews).
- (b) The notes numbering is canonical → the card/appendix would need renumbering (large prose impact; unlikely).
- (c) Leave the notes as a historical artifact and add a one-line "renumbered from 238/237" annotation.

No script changes are implied either way — the assertions and formulas are correct. F4 is RESOLVED (see the ## RESOLVED block at the end): apply the authorized notes renumber. No script changes are implied — the assertions/formulas are correct either way.

## RESOLVED — F4 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: renumber the notes prose to canonical numbering. Notes-only; Codex applies, Claude reviews. Do NOT touch scripts, paper.tex, or the appendix.
- In `notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md`, renumber every stale stage-number reference to the canonical numbering used by the card + scripts: "Stage 238" → "Stage 221" (title, §3, §9, and anywhere else it appears as this unit's number), and the wall-derivative-identity attribution "Stage 237" → "Stage 220".
- The math/content is identical; ONLY the stage numbers change.
- Acceptance: after editing, no stale `Stage 238`/`Stage 237` stage-number reference remains in that notes file. Append `## Applied: F4` with files_changed + summary + deviation.

## Applied: F4

- files_changed:
  - `notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md`
- summary: Renumbered stale Stage 238 and Stage 237 references, including script-path numerals, to the canonical Stage 221 and Stage 220 numbering.
- deviation: none
