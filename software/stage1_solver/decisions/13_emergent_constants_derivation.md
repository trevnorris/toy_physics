# Decision 13 — B2c verdict is UNDETERMINED; the real next chunk is the EMERGENT-CONSTANTS derivation

**Date:** 2026-06-19
**Status:** DECIDED direction (user-driven, 2026-06-19). This is the **resume-here record after /compact** for the
Path-A build. Supersedes the "rigorous MISS" reading of B2c (decision-12 B2c STATUS block — now flagged superseded).
**Mechanism:** B2c 3-round build/review → two audits (`pathA_17` validity, `pathA_18` dimensional) → two independent
verification agents → user methodology call (derive the emergent constants before `m̂0²·S_port`).

---

## 1. The headline: B2c does NOT give a validated hit-or-miss
B2c computed, on the self-consistent Path-A background, that the model's `P0` is **6.7–9.6 orders of magnitude below**
the GR quadrupole target `54/5 = 10.8` (measured at every converged τ: `P0=2.795e-9` at τ=1 … `1.22e-6` at
τ=0.029). The numbers are real. **But the comparison is `R_norm = m̂0²·S_port·P0 − 54Gc_s⁵/(5a⁵c⁵)`, and B2c pinned
`m̂0²·S_port = 1`.** Two audits + two independent verification agents established:

- **NOT a dimensional/units bug.** With units restored, the dimensional gaps the audit found are **order-neutral**
  (they equal `×1` under the natural-unit pins `a=c_s=ħ=m=1`); restoring them shifts NOTHING numerically. (This
  *corrected* `pathA_18`, which had overstated "8 inconsistencies" → really 2, and wrongly assumed a 4-spatial `G`.)
- **NOT a validated "missing physics" miss either.** `m̂0²·S_port` is **dimensionful** — it is the conversion factor
  `T_target/[P0]`, NOT a free dimensionless number (both verification agents, independently). Pinning it to `1` is a
  **scale choice** that *forces* the `D0→0` knife-edge. The pre-reg bakes this in: `T := 54Gc_s⁵/(5a⁵c⁵)/(m̂0²S_port)`,
  `D0 = N0/T` (`docs/pathA_preregistration.md` §253). So `m̂0²·S_port` directly sets whether calibration lands at a
  knife-edge (miss) or naturally (match).
- **The ~9-order "miss" is exactly the magnitude `m̂0²·S_port` would supply** (`10.8/P0 = 3.86e9` at τ=1). Whether the
  *derived* `m̂0²·S_port` is ≈1 (→ real miss) or carries large factors (→ match) is **plausible either way** (the
  healing length `a` is microscopic while the GR comparison is at the orbital/radiation scale, so a large
  `(scale-ratio)^n` is on the table) and **cannot be settled by dimensional analysis** — only by deriving it.

**Verdict: UNDETERMINED, hinging on the un-derived dimensionful normalization `m̂0²·S_port`.** B2c is NOT committed as
a miss. The B2c machinery (assembly/dual-engine/verdict-logic/warm-start/negative-controls) was fidelity-clean and is
retained; only the *interpretation* changed.

## 2. The foundational reframe (user-driven): the constants are EMERGENT, derive them first
`c`, `c_s`, `G` are NOT fundamental inputs in this model — they are **outputs of the 4D superfluid PDE**, and their
familiar "3+1" dimensions are the **brane (observed) dimensions**, while the **bulk PDE works with 4D-dimensioned
quantities**. The bulk→brane reduction (integrating out the transverse `w`-direction, width ~`a`) is where the
scale/dimension factors live — and `m̂0²·S_port` is almost certainly a piece of *that same bulk→brane normalization*.
So **deriving the emergent constants is logically PRIOR to `m̂0²·S_port`**; once the constants + reduction are in
closed form, `m̂0²·S_port` should largely fall out.

What the verification already established about the constants:
- **`c_s` (sound speed):** emergent from the EOS, `c_s² = 5Kρ⁴/m`; `[c_s]=L/T`; **value drifts with the background
  `ρ`** (not constant).
- **`G`:** emergent from defect back-reaction. Bulk 4-spatial Laplacian → `1/r²` (`[G_bulk]=L⁴T⁻²M⁻¹`); but the
  observed sector is **brane-projected to 3 spatial dims** → `1/r`, standard Burke-Thorne/Peters, `[G]=L³T⁻²M⁻¹`
  (confirmed from the parent action: defects live on the 3D brane, `w` transverse). **User hypothesis (likely right,
  must DERIVE): `G` reduces to the standard 3D form after the brane reduction, with the transverse-scale (`a`)
  factors being the content of the emergent `G`.**
- The GR target `54Gc_s⁵/(5a⁵c⁵)` is NOT dimensionless on its own (it carries `M⁻¹L⁻²T⁻²` with 3D `G`); it is used as
  the Peters **pure number `54/5`** with constants pinned to 1, and `m̂0²·S_port` silently absorbs the dimension.

## 3. User hypothesis to test: `c` ≠ `c_s` (wave speed vs terminal velocity)
**`c_s` = the speed waves travel *through* the medium (sound).** **`c` = possibly the terminal/drag-limited speed of
an *object* moving through the medium** (analogy: speed of sound vs an object falling through the fluid) — a
*different* physical quantity with a different origin, not just `c_s` relabeled. This would explain why `(c/c_s)³`
appears as a non-trivial bundle ratio (`R_tail`). **Let the math decide** — the constants derivation must determine
what `c` actually is (= `c_s`? a drag/terminal speed? the localized-Maxwell-sector speed?) and its relation to `c_s`.

## 4. NEXT CHUNK — the emergent-constants derivation (the agreed next step)
A full **reasoned + dimensional derivation** of the model's constants, from the parent action outward (NOT assumed —
derived; much already exists scattered across `research/pde_ledger/paper/parts/part01_parent_geometry.tex`, the
sound-speed ledger, and `research/4d_2_5pn/paper/4d_2_5pn.tex` — so partly consolidate + make rigorous + dimension-
check, partly fill gaps):
1. **`c_s`** from the EOS (`c_s²=5Kρ⁴/m`); its dimension + state-dependence.
2. **`c`** — what it physically IS (test the `c≠c_s` / terminal-velocity hypothesis), its definition, and its
   relation to `c_s` (the `c/c_s` ratio).
3. **`G`** — DERIVE from the defect back-reaction + the 4D→3D brane reduction: closed-form in superfluid params
   (`K,a,ρ,ħ,m`), bulk vs brane dimensions, the transverse-scale factors; TEST the "reduces to 3D `L³T⁻²M⁻¹`"
   hypothesis explicitly.
4. **The natural-unit → physical-scale map:** what the pins `a=c_s=ħ=m=1` actually stand in for — which determines
   whether `m̂0²·S_port` is genuinely ≈1 or carries orders.
Run everything through the committed harness `software/stage1_solver/src/stage1_solver/dimensional_check.py`
(machine-verify); produce a consolidated "constants & dimensions" reference doc (4D-bulk vs 3D-brane frame tags +
reduction maps). Discipline: Codex derives/codes, Claude reviews; dimensional-check + fidelity + adversarial.

**THEN:** derive `m̂0²·S_port` (should largely follow from #3/#4) → re-do the B2c `R_norm(τ)=0` calibration with the
derived normalization → real hit-or-miss verdict.

## 5. GATE-A freeze implication (methodology decision, do NOT do unilaterally)
`m̂0²·S_port=1` is currently part of the **GATE-A freeze** (decision-11 §2/§3, hash `ed3585…`). Deriving it means
**un-pinning a frozen convention → a documented freeze amendment (new hash)**. This is a methodology call to settle
with the user when the constants derivation is in hand (it may show the pin was fine, or that it must be replaced by
the derived value).

## 6. Standing step locked in (user request)
**Dimensional-consistency check (units restored, SymPy harness) BEFORE trusting any number** is now a STANDING gate:
`docs/pathA_preregistration.md` "STANDING VERIFICATION STEPS"; memory `[[feedback-dimensional-consistency-check]]`;
harness `dimensional_check.py`. Pairs with `[[feedback-transliteration-fidelity-audit]]`.

## 7. Audit trail (distilled; raw outputs regenerable, in gitignored `_scratch/`)
- `pathA_17` (validity): the "miss" = the pinned `m̂0²·S_port`; `54/5`←Burke-Thorne is correctly derived; no dropped
  4π/Jacobian. Verdict: ARTIFACT for the full claim, REAL only inside the pinned convention.
- `pathA_18` (dimensional, units restored): built the reusable harness; flagged dimensional gaps — BUT overstated
  (8→2 distinct; 4-spatial-`G` wrong). Honest caveat: dimensional analysis fixes no decimal order (factors =1 in
  natural units).
- Verification agent A (adversarial): 8→2 real issues, both order-neutral natural-units conventions, not code bugs;
  "if anything STRENGTHENS the real-missing-physics reading" (i.e. dimensions don't rescue the miss).
- Verification agent B (independent re-derivation): `[G]=L³T⁻²M⁻¹` effective-3D (brane projection, from the parent
  action); `m̂0²·S_port` is the dimensionful conversion factor `T_target/[P0]`, not dimensionless; `D0=K−B0−Z0` is a
  legal subtraction (`K`=modal stiffness eigenvalue) within the frozen normalization conventions.
