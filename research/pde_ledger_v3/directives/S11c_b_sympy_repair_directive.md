# S11c-b SymPy engine — repair directive (round 1)

## Task and authority
Patch `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (committed `0c0e9a8a`) to fix four
defects the two build legs found by FORM ablation. ⭐ **Fix only these; keep everything that passed unchanged**
(the §3a energy-basis spurion construction, the §3b operator, the N15 new invariants, the gradient-driven
coupling — all ablation-clean). Re-run the engine and regenerate `scripts/S11c_b_exports.py`. Those two files
are the only writes; every other file is read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (HEAD `0c0e9a8a`, with the
round-2/round-3-repaired §3a/§3c/§3d/§5a) is the sole physics authority and wins every conflict. The defects
and their ablation evidence are in `directives/_measurements/S11c_b_sympy_build_review.md` (B1–B4); read it.
Add no expected value or acceptance criterion (rule 5): the spec withholds the operand value, the coupling
grade/sign, the new invariants, the admissibility residual, and the basis count — the engine computes and
prints, and the diff happens on our side.

## The four fixes (implement the corrected spec sections; the finding record is the evidence)

**B1 — Admissibility operand is identically zero (BLOCKING).** `admissibility_operator` (L1944–1987) builds the
operand as `zero_wave(euler_derivative(wave_density, …))` — the ε→0 limit of the bilinear §3b wave operator,
≡ 0. Rebuild it per the corrected **§3d**: the **background-order (ε⁰) first variation of the §3a
variable-coefficient energy-and-geometry functional written in FULL fields** (`W=W_bg+δW`, `ρ_4D=ρ_4D,bg⁰(1+θ)`,
brane displacement), with the thickness/coefficient **gradient content of the full fields** (`∇(W_bg+δW)`, …)
so the profile's own gradients are present at background order — ⛔ `|∇(δW)|²` (perturbation-only gradient) is
not an acceptable background-order representative. Take the first variation wrt the brane configuration and
THEN evaluate at `𝔅⁰` (`θ=e_W=δW=0`, `u=0`; profile and its jets retained). ⭐ Retain every spatial derivative
the variation/divergence generates at the background-bookkeeper order of §3a (a second spatial derivative of
`W_bg` is still first order in background amplitude — do NOT discard it). ⛔ Do not use
`zero_wave(euler_derivative(wave_density,…))`, and do not insert `W_bg−W_0`. Emit operator operand, support
operand, residual in the same ordered pairing as `𝒮_hold⁰` (bulk-DOF body force + per-face traction).

**B2 — Coupling kernel not extracted from the operator (§3c).** `build_kernel` (L1834) uses the parallel
`paired_kernel_from_density` (direct mixed variation); `bulk_kernel_from_density`/`build_operator` are unused
for the bulk. Rebuild kernel extraction per the corrected **§3c**: a **weak variational restriction** of the
§3b operator under the §1c pairing, using independent transverse (solenoidal, `∇·u_T=0`) and longitudinal
(irrotational, `∇×u_L=0`) trial and test displacements plus θ/e_W test fields — the transverse→thickness block
is the transverse-trial × {θ,e_W,u_L}-test pairing, the reciprocal the other. This attributes both the
gradient terms AND the undifferentiated-`u` spurion couplings (`g·u`) without a global projector (N5), and it
must reduce to §1d's decoupled zero in the uniform limit. ⛔ Do not use the parallel direct-variation route,
and ⛔ do not implement the split by zeroing only undifferentiated field occurrences (inert on gradient
content). Take trial/test fields of compact in-plane support so the IBP boundary term is zero.

**B3 — Tautological controls (§3c adjointness + §5a independence).** (a) The adjointness residual must be the
pairing-based residual between the two extracted operator blocks (per §3c), ⛔ not `∂²U/∂u_T∂e_W −
∂²U/∂e_W∂u_T` (Clairaut ≡ 0); if adjoint by construction, emit the two blocks and state there is no independent
second route rather than a structural zero. (b) `task_independence` (L2140–2151) uses identical calls for
`COUPLING_KERNEL`/`ADMISSIBILITY` (residual `A−A`). Per corrected **§5a**, the one-sided source mutation must
propagate into every construction where the emitted dependence shows the source structurally present (the
kernel — now extracted from the mutated operator, B2; the admissibility — varied from the mutated background
functional, B1); an object is omitted only when the mutated source is structurally absent (computed, not
asserted). ⛔ Never emit an `A−A` residual as a control.

**B4 — Dormant dimension crash.** The `uniform_coefficient` fallback (L1150–1159) mints `a_s11cb_uniform_XX`
without a registered dimension; `dimension_of` (L816) then raises `KeyError` if the quotient selects an
uncatalogued uniform invariant. Register a dimension (or emit a graceful `_LOCAL_` report) so a future basis
selection does not kill the primary emit.

## Chain, F9, freeze — unchanged from the build directive
Import `LEDGER` from `S11c_a_exports.py`; write `S11c_b_exports.py` (previous LEDGER + this step's primary
§4 objects, merged flat); the §5 controls and §6 homogeneity are ablation evidence, NOT exports; F9 whole,
F9c prefix `s11c_b_`; `BUILD_INPUT_DIGESTS` = {this audit source, `S11c_a_exports.py`, `S11c_b_SHARED_PHYSICS.md`};
F6 publish only if every primary §4 task completed; D3 round-trip; `_RELATIONALS` reviver. The comparator is a
separate downstream artifact — not built here. Do not restate these; they are the build directive's, unchanged.

## Report (≤15 lines)
Files written, the four functions changed, tasks run, runtime, all emitted tag names, any F9c write, any
non-computable object. Confirm the untouched-passing objects were not regressed. Report no computed value.
