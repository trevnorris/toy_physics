# pathA_38 — Charge as a 4D throat-body interaction: does it give a brane-localized `1/r²` (Coulomb)?

**Status:** ⭐ DRAFT — setup sub-questions to resolve FIRST (next session), THEN gauntlet (Codex design-review → GLM → dual-engine →
tri-review). This is the **carry-over spec** written while full context is fresh (2026-07-03, pre-`/compact`). **Type:** throat/bulk
sector (4D), linearized inter-throat interaction. **Supersedes the mechanism of** `pathA_37` (the flow/counterflow gate — RETIRED;
kept as a documented waypoint). **Author:** orchestrator (scaffolding).

---

## §0. Why this gate — the crux in one paragraph

The v7 conceptual pivot (`docs/conceptual_foundation.md` §3/§4 ⭐ v7; memory `feedback-native-em-mechanisms`): **electric charge is not
a brane flow or deformation** (all three flow/deformation routes were explored and set aside — deformation → wrong energy `1/r⁶`;
one-fluid flow → charge=mass; counterflow → needs a 2nd component the `T=0` condensate lacks in flat 3D). Instead, **charge is the
interaction of two throats' 4D BODIES beyond the mouth** — a geometric interaction of the extended throat structures in the bulk, with
`±w` = the sign. **THE MAKE-OR-BREAK:** an electric field is long-range (`1/r²` force), so the throat-body interaction must come out
**`1/r²` in the 3D brane.** A source acting freely through the **4D bulk** gives the WRONG falloff (4-space → `1/R²` potential →
`1/R³` force); it is correct 3D Coulomb **only if the throat's influence is brane-localized.** This gate tests exactly that.

---

## §1. The leverage — pathA_29 already did this for GRAVITY (reuse the machinery)

**Do not reinvent the localization analysis — pathA_29 is the template.** `pathA_29` (`RETURN_RESIDUAL_PREDICTION`, tri-reviewed;
`reports/pathA_29_*`) proved that **gravity's field survives the finite brane slab as `1/r²`, not `1/r³`**, because the mediating field
has a **normalizable transverse (`w`-direction) zero mode localized to the slab** → the in-brane field solves a real 3D-radial
equation → `p=2` (`1/r²`). A *delocalizing* warp gives `p=3` (`1/r³`) and was the counterfactual that got rejected. **This is the
Randall–Sundrum-style brane-localization mechanism**, already demonstrated in-repo for the gravity sector.

**The pathA_38 question is the same test applied to the CHARGE (throat-body) sector:** does the field that mediates the inter-throat
electric interaction inherit a **normalizable localized transverse zero mode → `1/r²`**, or does it **delocalize into the bulk →
`1/r³` (FAIL)**? If the charge sector localizes the *same way* gravity does, the forces share one localization mechanism — a strong
"connected" result. If it does not, that is a first-class no-go (charge can't be long-range Coulomb in this medium).

---

## §2. ⭐ SETUP SUB-QUESTIONS TO RESOLVE FIRST (next session, before writing the executable gate)

These are genuinely open and must be pinned down (Claude + Codex + GLM) BEFORE the dual-engine derivation. Do not skip.

1. **What field mediates the throat-body electric interaction?** (The load-bearing open piece.) Candidates to adjudicate:
   - a linearized perturbation of the **bulk medium** sourced by the throat's 4D body (a bulk field `δΦ(x,w)` the two throats exchange);
   - the brane's **extrinsic-curvature / embedding** response at the throat (how the brane bends into `w`);
   - a **tension/geometric** field of the throat funnel.
   Pin down: what is the linearized field, and what is its **wave operator** in the (brane-slab + bulk) geometry — i.e. the analog of
   pathA_29's transverse operator whose zero-mode normalizability decides `p=2` vs `p=3`.
2. **The throat-body geometry model.** A tractable linearized model of the throat as a funnel (area > mouth; extrinsic curvature into
   `w`). **Leverage** `[[project-light-4d-throat-hypothesis]]` (the 4D-throat/geon structure) and pathA_29's finite-slab setup. Keep
   it as reduced as pathA_29's — do NOT require the full nonlinear throat interior (that is sim-deferred).
3. **How `±w` enters (the sign).** The charge sign from the throat orientation. The interaction must give **like-repel / unlike-attract**
   — GLM already confirmed (Green's identity + Bernoulli) that a positive-energy interaction gives the right sign IF the field is
   mass-neutral; carry that result forward and check it survives in the localized-mode calculation.

---

## §3. The able-to-fail test (once §2 is resolved)

Compute the **inter-throat interaction energy `U_int(R)`** as a function of 3D brane separation `R`, from the §2 mediating field on
the finite-slab + bulk geometry:
- **PASS (`THROAT_ELECTRIC_LOCALIZED_COULOMB`):** `U_int ∼ 1/R` (Coulomb potential → `1/R²` force), arising from a **normalizable
  localized transverse zero mode** (the pathA_29 mechanism), with `±w` giving like-repel/unlike-attract.
- **FAIL (`FAIL_DELOCALIZED_BULK_1_OVER_R3`):** `U_int ∼ 1/R²` (→ `1/R³` force) — the interaction leaks into the 4D bulk; no
  localized mode; wrong power. (This is the pathA_29 `p=3` counterfactual, now as a *reachable verdict*.)
- **Counterfactual guard (able-to-fail, à la pathA_29):** a non-localizing warp / non-normalizable mode MUST flip the exponent to
  `1/R³` — the test can genuinely fail.
- **Sign guard:** an all-attract (gravity-like) result → `FAIL_SIGN_STRUCTURE`.
- Dual-engine (SymPy + Mathematica), `ENGINE_AGREE` on the exponent, the mode normalizability, and the sign — per house rules.

**If PASS:** the electric sector shares gravity's brane-localization → `1/r²` Coulomb with `±w` sign, from the throat's 4D body →
**the forces are connected**, and the charge sector becomes a calibrated entry in the central `pde_ledger` alongside pathA_29's
gravity localization ([[project-calibrated-pde-goal]]).

---

## §4. Honest scope (the algebra/sim boundary)

- **In-scope for algebra:** the *linearized* inter-throat interaction + the *transverse-mode normalizability* analysis. pathA_29 did
  exactly this for gravity with reduced (algebra-tractable) machinery — the charge version is the same class of calculation.
- **Sim-deferred:** the full nonlinear throat interior, the self-consistent throat-body shape, dynamics. Those stay posed as
  sim-dependent open items ([[project-simulation-deferred-complete-pde-strategy]]) — completing this gate completes the *spec* of the
  charge sector's `1/r²`, it does NOT prove the full nonlinear theory.
- **Harder than the flat-brane far-field:** this is throat-interior / 4D / bulk territory (Gate-T neighborhood), so expect the setup
  (§2) to take real work before the gate is executable.

---

## §5. References (read first, next session)
- `docs/conceptual_foundation.md` §3 + §4 ⭐ v7 (the charge = 4D-throat-body mechanism + this `1/r²` crux).
- `reports/pathA_29_*` + `directives/pathA_29_brane_bulk_return.md` (the gravity brane-localization template — the machinery to reuse).
- `directives/pathA_37_c5_throat_electrostatics.md` (the RETIRED flow-gate — the documented exploration that led here; the GLM reviews
  `_scratch/pathA_37_v5_glm_review.md` established the sign result + the flat-3D counterflow no-go).
- Memories: `[[project-brane-existence-defect-structure]]` (UPDATE 2026-07-03b), `[[feedback-native-em-mechanisms]]` (v7),
  `[[project-light-4d-throat-hypothesis]]` (4D throat structure), `[[project-pn-gravity-ladder]]` / `[[project-calibrated-pde-goal]]`.

---

## §6. Changelog
- v0 DRAFT (2026-07-03) — carry-over spec written pre-`/compact` while context is fresh. Nails the crux (`1/r²` via brane-localization,
  the pathA_29 template), lists the §2 setup sub-questions to resolve first, and the §3 able-to-fail test. NOT yet gauntleted — next
  session: resolve §2 → Codex design-review → GLM → dual-engine → tri-review.
