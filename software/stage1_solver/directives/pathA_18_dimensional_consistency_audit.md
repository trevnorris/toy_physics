# Directive pathA_18 — SymPy-enforced dimensional-consistency audit of the Path-A → R_norm chain (units restored)

**Date:** 2026-06-19
**Owner:** Codex (verification/audit; writes a reusable read-only SymPy harness). Claude reviews afterward.
**Trigger:** User methodological call (correct, and a standing gap): we have NOT been doing tool-enforced dimensional
checks. We work in **natural units (`a=c_s=ħ=m=1`)**, which is exactly where a dropped factor (e.g. an `a⁵`, a `c_s⁵`,
a 4D→reduced volume/area Jacobian) hides — it reads as "×1" until the scale isn't 1, then resurfaces as orders of
magnitude. The B2c "miss" is ~9 orders; a dimensional artifact is a prime suspect. **Restore units symbolically and
machine-check every equation is dimensionally homogeneous BEFORE trusting any number.**

## Scope & stance
VERIFICATION ONLY — do not change the model/formulas. Build a **reusable** SymPy dimensional-check harness and run it
over the Path-A chain that feeds the B2c verdict. Deliverable = per-equation CONSISTENT/INCONSISTENT, with any
mismatch's missing factor named + its dimension/order. Read-only; harness + findings under `_scratch/` or `runs/`
(gitignored is fine for run artifacts, but the **harness script itself should live under
`software/stage1_solver/src/` or `tools/` so it persists as a standing check** — confirm location with the existing
tree; if unsure, put it in `software/stage1_solver/src/stage1_solver/dimensional_check.py`). Do NOT use `$RT
exec-sympy` (MANIFEST race); use a standalone `python3` script. Never touch `research/pde_audit/simulation/` or
`physical_export_permitted`. Each script ≤10 min.

## Step 0 — Establish the dimensional system (the foundation everything rests on)
Determine, FROM THE ACTION/derivation (not by guessing), the base dimensions and every quantity's dimensions for the
**4-spatial-dimension** superfluid analog. Pin down:
- Base set (e.g. length `L`, time `T`, mass `M`; note any `ħ=m=1` choices and what they identify).
- `ρ` (condensate density in 4D), `ψ` (field), the **EOS** `P=Kρ⁵`, `h=(5K/4)ρ⁴`, `U=(K/4)ρ⁵` — and hence `[K]`.
- `c_s` (sound speed `L/T`), `a` (core/healing length `L`), the gauge field `A` + gauge charge `q`, throat radius `R`,
  wall coordinate `w`, wall stiffness `τ`, frequency `ω`/`ϖ`.
- **`G` in 4 spatial dimensions** — Newton's constant has D-dependent dimensions; determine `[G]` for this analog and
  hence `[gamma_GR]=[2G/(5c⁵)]`. **This is a prime suspect:** if `gamma_GR` was written with 3D-GR conventions but
  matched to a 4D-analog quantity, it is a dimensional mismatch that would read as a multi-order "miss."
Record this as a single dimensional dictionary the harness consumes.

## Audit items
### D1 — Forward sector equations
Dimensional homogeneity of: the gauged-GPE / static balance, the EOS trio, the wall/throat static-balance operator
(`patha_static_balance.py`), the return-source `Σ_i ΔV_i^r·(−k1·ρ0)`, the BdG operator, and the Maxwell operator.
Every additive term in each equation must share dimensions. Flag any term that doesn't.

### D2 — The REDUCTION steps (highest suspicion — the user's sphere-vs-circle point)
For every reduction/projection from the 4D fields to the reduced coefficients `{K,M,B_n,Z_n,N_n}`:
- What is the **integration measure** at each stage (4D volume `d⁴x`? a 3D shell? the 1D wall line `dw`? shell volumes
  `ΔV_i^r=(4π/3)(r³)` — note the `4π/3` is a **3-ball** volume; is the ambient space 4D, and if so is the correct
  4-ball/3-sphere measure used)? Track the dimension the measure contributes at each step.
- The mode-normalization conventions (`χᵀWχ=1`, the BdG `Re(ψ̄u+ψv)` normalization, the trapz weight): what
  dimensions do `χ`, `W`, the weights carry, and is the normalization dimensionally self-consistent?
- The overlap integrals (`I_eta_phi`, etc.), the couplings `c_j=λ_B·I`, the moments `B_n=Σc_j²/ϖ^{2(n+1)}`, the
  Maxwell `Z_n,N_n` (with `Γ_port=4a⁵/27c_s⁵`, the `ω⁵`).
Confirm each reduced coefficient (`K,M,B_n,Z_n,N_n`) ends up with the dimension the §5 algebra and the `R_norm`
formula ASSUME. **A reduction that changes the measure's dimension without a compensating factor is exactly the
sought artifact — quantify it.**

### D3 — The matching dictionary / R_norm
- `D0=K−B0−Z0`, `D2,D4`, `P0=N0/D0`, `P2,P4`, `R_pole`: every additive term dimensionally equal? Is `P0` the
  dimension the target assumes?
- `R_norm = m̂0²·S_port·P0 − 54·G·c_s⁵/(5·a⁵·c⁵)`: do the two terms share dimensions? Is `P0` dimensionless, and is
  `54G c_s⁵/(5a⁵c⁵)` dimensionless with the 4D `[G]` from Step 0? If NOT, the comparison is dimensionally invalid and
  the "miss" is (at least partly) a units artifact — quantify the offending power (of `a`, `c_s`, `c`, `G`).
- Check whether the natural-unit pins (`a=c_s=ħ=m=1`, `m̂0=S_port=1`) hide a factor that, restored, shifts the target
  by orders.

### D4 — Verdict
State: is the Path-A → `R_norm` chain **dimensionally consistent end-to-end**? If YES, the ~9-order miss is NOT a
dimensional artifact (strengthens "real missing physics"). If NO, name the exact step, the missing factor, its
dimension, and the order of magnitude it contributes in natural units — i.e. how much of the 9 orders it explains.

## Acceptance
A reusable SymPy dimensional-check harness exists and runs; Step 0 dictionary recorded; D1–D3 each report per-equation
CONSISTENT/INCONSISTENT with mismatches quantified; D4 gives a clear consistent-vs-artifact verdict tied to specific
factors. No model/formula changes.

## Review (orchestrator, after Codex)
Claude independently re-derives the key dimensions (esp. `[G]` in 4D, the reduction measures, `[P0]`) and runs a clean
adversarial pass on the harness + findings before we accept or revise the B2c verdict. This audit + `pathA_17`
together decide: commit B2c as a real miss, or revise (artifact). **Adopt the dimensional-check harness as a STANDING
pre-numbers step going forward.**
