# pathA_39 — Magnetism = the brane swirl: does a moving charge close Maxwell (B, Lorentz force, curl equations)?

**Status:** ⭐ DRAFT — carry-over spec (2026-07-03, pre-`/compact`), the running start for the NEXT front after `pathA_38` (charge
`1/r²` EARNED). **Type:** brane/throat EM sector; the fourth and last force-sector. **Downstream of** `pathA_38` (electric charge = the
throat-body Goldstone → Coulomb `E`) and `pathA_36` (light = the brane shear wave). **First: resolve §2 setup (Claude↔Codex→GLM), then
gauntlet** (Codex design-review → GLM → dual-engine → tri-review). **Author:** orchestrator (scaffolding).

---

## §0. Why this gate — the crux in one paragraph
The four-way taxonomy (`docs/conceptual_foundation.md` ⭐ FOUR-SECTOR CHAIN / §3 v7): **magnetism = the swirl / vorticity of the
brane.** A **moving** charge (a moving throat) drags/twists the surrounding brane → the magnetic field `B`. EM is ONE thing: `pathA_38`
gave the static charge's Coulomb `E`; `pathA_36` gave the propagating brane-shear wave (the photon). **THE CRUX:** does the *native*
brane-swirl `B` sourced by a moving throat **close the full Maxwell set** — the two curl equations (Ampère + Faraday, with displacement
current), the **Lorentz force** on a second moving charge, and the right relative `E↔B` normalization (so the EM wave built from `E,B`
travels at `c_γ`)? I.e. is the native magnetic mechanism CONSISTENT with the electric sector under motion, or does it depart? A
departure would be a **falsifiable** prediction (and possibly a feature); a clean closure would complete the fourth sector.

---

## §1. The leverage — largely DOWNSTREAM of pathA_38 + pathA_36 (reuse the machinery)
- **`pathA_38`** froze the throat-body charge → Coulomb `E = ∇×(nothing)`-type `1/r²` from the transverse-embedding Goldstone `h`, with
  `±w` sign and the `w→−w` parity split. A moving throat is the same object with a velocity `V`.
- **`pathA_36`** established the brane carries a transverse shear wave at `c_γ²=μ_R/ρ_br` — the photon. The full EM field `(E,B)` must
  reduce to that wave in vacuum.
- **Standard-EM anchor (target-blind):** in ordinary EM, `B` is not independent — it is `E` in a boosted frame, and a moving point
  charge gives the Biot–Savart / Liénard–Wiechert `B = (V×E)/c²` to leading order. The gate is whether the medium's **native swirl**
  reproduces this (with `c→c_γ`), i.e. whether the brane's vorticity response to a moving throat MATCHES the covariant requirement. If
  the medium is (approximately) Lorentz-covariant — an inherited/conceded wall — `B` should fall out; the test is whether the native
  mechanism AGREES, and where it departs.

---

## §2. ⭐ SETUP SUB-QUESTIONS TO RESOLVE FIRST (Claude↔Codex→GLM, before the executable gate)
1. **What is the swirl field, natively, and how does a moving throat source it?** Is `B` the brane vorticity `ω = ∇×v_brane` (a flow
   swirl), or a curl of the embedding/orientation field, or the magnetic partner of `h` under motion? Pin the linearized field + its
   source law for a throat moving at velocity `V` (the analog of pathA_38's `q_eff`, now with a current `q_eff·V`). Keep it on the SAME
   brane+bulk geometry (χ_B wall + finite slab) as pathA_38.
2. **Does `E` (pathA_38) + `B` (swirl) close the full Maxwell set?** The two curl equations with displacement current, from the medium's
   own dynamics — DERIVED, not imposed. What is the computed relative `E↔B` normalization, and does the resulting wave equation give
   `c_γ` (tie to pathA_36)?
3. **The Lorentz force.** Does a second moving charge in the swirl feel `F = q(E + V×B)` with the right sign/magnitude, from the same
   throat-body interaction machinery?
4. **`λγ` / covariance consistency.** Is the light speed from the `(E,B)` wave the same `c_γ` as pathA_36's shear wave, in ONE parameter
   set? (This is the hook into the consistency knit — `λγ=c_γ/c_s=1` cone-lock.)

---

## §3. The able-to-fail test (once §2 is resolved)
- **PASS (`MAGNETISM_CLOSES_MAXWELL`):** a moving throat's native brane-swirl `B` closes the two curl equations + gives the Lorentz
  force `q(E+V×B)` + the right `E↔B` normalization (EM wave → `c_γ`), all DERIVED from the medium on the pathA_38/36 geometry.
- **FAIL modes (each reachable, each physical):** `FAIL_MAXWELL_NOT_CLOSED` (curl equations don't close / need an imposed term);
  `FAIL_WRONG_EB_NORMALIZATION` (`E` and `B` don't share `c_γ` → the EM wave speed disagrees with pathA_36); `FAIL_NO_LORENTZ_FORCE`
  (no `V×B` force, or wrong sign); `FAIL_SWIRL_NOT_SOURCED` (a moving throat doesn't source the swirl / `B=0`); `FAIL_COVARIANCE_DEPART`
  (a computed departure from Maxwell — report it precisely; may be a falsifiable feature, cf. pathA_29's longitudinal-radiation hook).
- **Anti-circularity:** `B`, the closure, the Lorentz force, and the normalization must be OUTPUTS of solves; NO hand-inserted `V×E/c²`
  read back. **Empirically able-to-fail** (the pathA_38 lesson): force a non-closing/wrong-normalization input → the classifier must
  EMIT the `FAIL_*`, not raise. Dual-engine (SymPy + Mathematica), `ENGINE_AGREE`; dimensional firewall with able-to-fail ablations.

---

## §4. Honest scope
- **In-scope for algebra:** the linearized moving-throat swirl + the Maxwell-closure/normalization/Lorentz-force algebra, reusing the
  pathA_38 (charge) + pathA_36 (light) + pathA_29 (localization) machinery on the χ_B brane+bulk geometry.
- **Sim-deferred:** the full nonlinear moving-throat interior + self-consistent dynamics (inherits pathA_38's SIM-DEFERRED risks).
- **This is the FOURTH sector** — a PASS makes all four forces earned from one brane+bulk; it does NOT by itself do the consistency
  knit (`λγ` + NG gauntlet + `pde_ledger` assembly), which is the separate closing step ([[project-brane-existence-defect-structure]]).

---

## §5. References (read first)
- `docs/conceptual_foundation.md` — ⭐ THE FOUR-SECTOR CHAIN + §3 v7 (magnetism = brane swirl).
- `directives/pathA_38_throat_body_electric_localization.md` + `reports/pathA_38_*` (the static charge `E` this builds on).
- `reports/pathA_36_c5_phase_potential.md` (light / the photon at `c_γ`).
- `reports/pathA_29_results.yaml` (the localization machinery).
- Memories: `[[project-brane-existence-defect-structure]]` (UPDATE 2026-07-03c), `[[feedback-native-em-mechanisms]]`,
  `[[feedback-negative-verdict-short-circuit]]` (the pass-by-construction / able-to-fail lesson), `[[project-calibrated-pde-goal]]`.

---

## §6. Changelog
- v0 DRAFT (2026-07-03) — carry-over spec written pre-`/compact` after `pathA_38` (charge) EARNED. Nails the magnetism crux (moving
  throat → swirl → close Maxwell), the leverage (downstream of pathA_38 + pathA_36), the §2 setup sub-questions to resolve first, and
  the §3 able-to-fail test. NOT yet gauntleted — next: resolve §2 (Claude↔Codex → GLM) → Codex design-review → GLM → dual-engine →
  tri-review.
