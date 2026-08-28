# S11c-b SymPy engine — B1 round-2 fix (admissibility scalar over-promotion)

## Task and authority
Patch `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (committed `10eaa75c`) to fix ONE
residual defect in the admissibility operand. ⭐ **Fix only this; keep everything else unchanged** — B2/B3/B4
and the correct bending content of B1 are ablation-clean and must be preserved. Re-run and regenerate
`scripts/S11c_b_exports.py`. Those two files are the only writes.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §3d is the physics authority
(unchanged — the spec is correct; this is a SymPy implementation over-read). Add no expected value (rule 5).

## The defect (B1 round-2) — a scalar over-promotion injects a spurious jet-independent force
The B1 repair correctly made the admissibility operand the background-order first variation of the full-field
energy and correctly recovers the κ_W bending content (the `W_bg` second-jet terms — KEEP these). BUT it
**over-promoted the SCALAR perturbation fields**: it substituted `btheta→(1+θ)` and `be→(1+e)` into the §3a
invariants that are **perturbation bilinears** (`½B_ρ⁽³⁾θ²`, `C W₀ θ e_W`, `½k_W W₀² e_W²`, …). Those terms
are quadratic in the perturbation, so at 𝔅⁰ they and their first variation vanish by construction; the
`(1+θ)/(1+e)` promotion instead makes them contribute a spurious **jet-independent** body force
(`∂/∂θ[½B(1+θ)²]|₀ = B` — a nonzero linear piece where the real `½Bθ²` gives 0), which is the defect to
remove. The remedy is structural, below; ⛔ its verification is ours, not an emitted check.

Evidence: `directives/_measurements/S11c_b_sympy_repair_review.md` (B1). The sibling blind Wolfram engine
implements §3d correctly and is the reference: its `constructFullFieldBackgroundEnergy` keeps the scalar
perturbation fields as perturbations (∝ background order, vanishing at first variation) and full-fields **only
the gradient slot** (`gradient[fullWidth] = ∇(W_bg+δW)`).

## The fix
In the admissibility full-field construction, **full-field only the gradient content**; keep the scalar DOFs
as perturbations:
- keep `btheta → θ` (the densification perturbation), `be → (W₀/W_bg) e_W` (the exact `e_W,bg` map);
  ⛔ do NOT promote them to `(1+θ)` / `(1+e)`.
- put the full-field promotion ONLY in the gradient slots that carry a background profile — the thickness
  gradient `∇(δW) → ∇(W_bg+δW)` and any coefficient gradient `∇(coeff)` wherever a coefficient varies — so the
  profile's own gradients (and the second spatial derivative the variation generates) enter at background
  order. This is exactly the WL `fullWidth`/`gradient[fullWidth]` construction; mirror it.

The fix is purely structural — it is the gradient-only promotion above, mirroring WL. ⛔ Do not add any guard,
acceptance check, or emitted assertion about the operand's value (in the uniform limit or anywhere): the
engine computes and prints, and the diff happens on our side (rule 5). The correct behaviour follows from the
construction itself (the perturbation-bilinear terms are quadratic in the perturbation, so they and their
first variation vanish at 𝔅⁰ once the scalars are not promoted; only the gradient/profile content survives).
⛔ Do not touch the B2 kernel extraction, the B3 controls, the B4 dimension fallback, the §3a energy basis, or
the §3b operator; keep the κ_W background-jet bending content of the admissibility operand unchanged.

## Report (≤10 lines)
The function(s) changed, tasks run, runtime, confirmation that B2/B3/B4 and the κ_W bending content were not
regressed, and that the scalar over-promotion is removed. No computed value.
