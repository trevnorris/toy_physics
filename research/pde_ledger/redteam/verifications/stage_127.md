---
unit_id: 127
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 127

## Per-finding outcomes

### F1 — paper_misalignment (banner / H1 mismatch)

**Classification:** resolved

**What changed:**
Cluster A mass-fix landed across the three label sites identified in the original report:
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py:12` now reads `banner("STAGE 127 — GEOMETRIC MOUTH-PENETRATION FAMILIES")` (was `STAGE 110`).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:26` now reads `banner["STAGE 127 — GEOMETRIC MOUTH-PENETRATION FAMILIES"];` (was `STAGE 110`).
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md:2` now reads `# Moving-Throat PDE — Stage 127: Geometric Mouth-Penetration Families` (was `Stage 229`).
- The paper card title was already `Stage 127` and remains unchanged.

The saved `.txt` transcripts also reflect the new banner (`STAGE 127 — GEOMETRIC MOUTH-PENETRATION FAMILIES` heading in both outputs). The Mathematica end-message `Stage 127 Mathematica audit passed.` still agrees.

**Assessment:**
The chosen direction matches the orchestrator's Cluster A resolution (recommended direction (a)/(c)). All four label sites now agree on "Stage 127". No numerical change. No collateral edits beyond the label strings on the three lines named.

### F2 — mathematica_transliteration (independent integral derivation)

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:28–49` now performs an independent symbolic integration before the existing closed-form usage and root-find:
- Line 28: `Clear[x, z];` (added `z`).
- Line 29: `$Assumptions = Element[x, Reals] && x > 0 && Element[z, Reals];` (added `z` real).
- Lines 36–41: `gSlabFromIntegral = FullSimplify[Integrate[(1/x)*Cos[Pi*z/2], {z, 0, x}], Assumptions -> x > 0];` followed by `slabClosedFormResidual = FullSimplify[gSlabFromIntegral - gSlab, ...]` and a `pass`/`fail` gate.
- Lines 44–49: analogous `gExpFromIntegral = FullSimplify[Integrate[(Exp[-z/x]/(x*(1 - Exp[-1/x])))*Cos[Pi*z/2], {z, 0, 1}], ...]` with a comparable residual gate.
- The downstream `FindRoot` + `expectApprox` blocks (lines 54–62) are unchanged in shape.

**Assessment:**
Both new gates exercise `Integrate` rather than re-typing the closed forms, so a typo or sign error in the algebraic forms used for the root-find would now be caught by the symbolic residual. The exec log shows the two new `PASS:` lines: `PASS: slab closed-form matches source integral` and `PASS: exp closed-form matches source integral`, and the existing `PASS: slab compensation root` / `PASS: exponential compensation root` lines are still present. The SymPy script is unchanged, as the directive required. No regressions visible.

The two new assertions are non-tautological — they compare `Integrate[...]` (Mathematica's own evaluator) against the algebraically-rearranged closed form via `FullSimplify[...] === 0`, so a discrepancy would show up as a non-zero residual passed to `fail[]`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `g_-^F1 = 0.758035078944662826919680890414`
- `x_*^slab = 0.79783936090456440508490881279972941789601018958900271741142733512870884256665312`
- `residual slab = -2.66696092767403190071722842020e-83`
- `residual exp  = -9.67012424716328641076781154243e-82`

**Mathematica:** exit=0. Notable lines:
- `PASS: slab closed-form matches source integral`
- `PASS: exp closed-form matches source integral`
- `PASS: slab compensation root`
- `PASS: exponential compensation root`
- `Stage 127 Mathematica audit passed.`

**Output freshness:**
- `scripts/output/...sympy_audit.txt` mtime 2026-05-27 17:45:37; corresponding `.py` mtime 2026-05-27 17:20:55 — output newer.
- `mathematica/output/...mathematica_audit.txt` mtime 2026-05-27 17:46:25; corresponding `.wl` mtime 2026-05-27 17:22:57 — output newer.
Both outputs were re-generated after the edits landed.

## Material-change assessment

`material_change`: false.

No derived numerical value changed. `g_-^{F1}`, `x_*^slab`, and `x_*^exp` are bit-identical pre- and post-fix to all displayed digits. The Mathematica edit added a verification step (`Integrate` vs. closed form) that confirms the same closed forms previously taken as given. The label edits are cosmetic. No downstream unit needs re-audit on numerical grounds.

## Side observations (non-blocking)

- The Mathematica output prints `g_exp(x) = (2*(2*E^x^(-1) + Pi*x))/((-1 + E^x^(-1))*(4 + Pi^2*x^2))` rather than the form `2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x]))` defined in the script. This is `FullSimplify`'s preferred rearrangement (multiply num/denom by `E^{1/x}`); algebraically equivalent, and the new closed-form-vs-integral PASS confirms the equivalence. Not a defect.
- SymPy prints `g_exp(x) = (2*pi*x + 4*exp(1/x))/((pi**2*x**2 + 4)*(exp(1/x) - 1))` for the same reason (sp.simplify's rearrangement). Same algebraic form, same value.

## Verdict justification

Both findings landed as the directive (with Cluster A resolution for F1) required. SymPy and Mathematica both exit 0; both new `PASS:` lines for the independent integral check appear in the Mathematica transcript, the two existing root-residual PASSes remain, and all four label sites now read "Stage 127". Outputs were refreshed post-edit. Numerical results are unchanged, so no downstream invalidation. Verdict: `verified`.
