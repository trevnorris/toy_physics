# S11c-b SymPy engine — repair round-1 review, two legs (DISAGREED on B1), consolidated

**Artifact:** `scripts/S11c_b_brane_operator_sympy_audit.py` (Codex-repaired round 1, functions
`admissibility_operator`/`build_kernel`/`task_independence`/`uniform_coefficient`). **Legs:** fresh Claude
agent + Grok, prompt `_legs/S11c_b_sympy_repair_review.md`. **Raw:** Grok
`~/.s11_build/S11c_b_sympy_repair_grok.txt`; Claude-agent `…/scratchpad/s11cb_leg2/`; Grok
`/tmp/s11cb_repair_review_grok/`. Verdicts are the orchestrator's after adjudication (rule 13).

## B2/B3/B4 — FIXED (both legs agree, ablation-verified)
- **B2** kernel now from `build_operator` via `weak_operator_blocks`; deleting the operator BREAKS the kernel;
  the off-diagonal block collapses to exactly 0 at jet→0 (§1d decoupled zero, no gradient-independent
  survivor). The old direct routes (`paired_kernel_from_density`/`mixed_variation`/`bulk_kernel_from_density`)
  are now dead code (harmless; cleanup candidate).
- **B3** adjointness emits `NO_INDEPENDENT_SECOND_ROUTE` (the two operator blocks), no Clairaut
  `PAIRING_RESIDUAL`; the independence corruption reaches SLAB_OPERATOR + COUPLING_KERNEL + ADMISSIBILITY
  (no A−A).
- **B4** the `uniform_coefficient` fallback registers a dimension (`dim_div(DIM_ENERGY, invariant_dimension)`);
  no KeyError. No regression: N15 basis (count 26, 16 new invariants), operator, gradient-driven coupling all
  still ablation-clean.

## B1 — the legs DISAGREED; adjudicated NOT fully fixed (Grok right)
- **Claude agent:** "REPAIRED" — verified the operand carries the correct bending content
  `−κ_W W_bg³ W_bg″` (matching its own independent derivation, computed not typed — source ablation of the
  second-jet map removes exactly those terms) and moves under a jet ablation.
- **Grok:** partially fixed — a spurious **jet-independent** force survives the uniform limit:
  `THETA = B_ρ⁽³⁾ + C W₀`, `E_W = C W₀ + k_W W₀²`, from the repair promoting the SCALAR perturbation fields
  `btheta→(1+θ)`, `be→(1+e)` inside the perturbation-bilinear energy terms (`½Bθ²`, `CW₀θe_W`, …).
- **Adjudication (rule 13):** both true simultaneously — the operand = (correct bending force in the jets,
  Claude-agent verified) + (a spurious jet-independent piece from the scalar over-promotion, Grok found;
  Claude-agent checked "moves under jets?" but did not probe the jet=0 survivor). Verified the arithmetic:
  `∂/∂θ[½B(1+θ)²]|₀ = B` vs `∂/∂θ[½Bθ²]|₀ = 0`. Physics: a UNIFORM unloaded slab must source ZERO background
  force; the over-promotion injects one. **B1 is NOT fully fixed.**

## Spec-vs-engine tiebreaker → SymPy-ONLY (no spec change)
Read the sibling WL engine's `constructFullFieldBackgroundEnergy`: it full-fields ONLY the gradient slot
(`gradient[fullWidth] = ∇(W_bg+δW)`, invariant #7) and keeps the scalars (`thetaVariation`, `localEw`) as
perturbations (∝ background order, vanishing at first variation) — the CORRECT construction, with no scalar
promotion. So the blind WL engine read §3d correctly; §3d is adequate; SymPy (Codex) over-read it. **B1 fix is
a SymPy-only round-2 engine patch — no §3d change.** (Rule 1: the two-SymPy-leg disagreement sent me to the
sibling engine, whose independent-correct construction both vindicates the spec and localizes the bug.)

**Repair (SymPy B1 round-2):** in the admissibility full-field substitution, keep `btheta→θ`,
`be→(W₀/W_bg)e_W` (perturbation DOFs / the `e_W,bg` map); put full-field content ONLY in the gradient slots
(`∇(W_bg+δW)`), mirroring WL's `fullWidth`. The correct bending content (κ_W jets) is already present and must
be preserved.

## Non-blocking (step record)
- Admissibility per-face traction is `±e_force/W₀` — a rescaling of the E_W body-force, not an independently
  derived boundary-flux construction (data-dependent, non-vacuous, so B1's core is met). Consider on the WL
  side too.
- Dead code (`paired_kernel_from_density`, `mixed_variation`, `combined_sector_substitution`,
  `bulk_kernel_from_density`) — harmless cleanup.
