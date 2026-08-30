# S11c-b STEP 3 — WL engine repair directive (v2: one item)

**Target engine (the only file this build modifies):**
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
**Regenerated deliverable:** `research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out`

The WL engine is **blind Mathematica**: it re-derives the S11c-b variable-coefficient brane operator from the
shared physics from scratch and **imports nothing** from the SymPy engine or any prior ledger. That property
is load-bearing and must survive this repair unchanged — the whole method is that the two engines can
genuinely disagree, and the disagreement is the measurement. **Do not read, import, or reference the SymPy
engine, its output, or any other engine's result while performing this repair.** Re-derive every changed
object from the shared physics (`directives/S11c_b_SHARED_PHYSICS.md`) alone.

> **v2 scope note.** A v1 of this directive proposed a second repair to the thickness kinetic inertia
> coefficient (`kineticEw = muW WZero^2 …`). Two independent decision-list derivations withdrew it: §1a
> defines `e_W ≡ δW/W_0` with constant `W_0`, so `μ_W W_0²` is the faithful inertia in the named `e_W`
> coordinate — that coefficient is **correct and is not to be touched**. This v2 has **one** item.

## The three standing clauses (bind every emit in this engine)
1. The engine may **PRINT** computed CAS objects. It may **not** state conclusions.
2. **PRINT** the operand; do **not** assert it equals anything. There is no `expected` value in this build.
3. Interpretation belongs to the step record.

⛔ **Withheld deliberately:** this directive does **not** state what the repaired object equals, its value,
coefficient, sign, order, or its residual against any sibling engine, and does **not** tell you whether the
object is currently zero or nonzero. Those are outputs. Re-derive from the shared physics; emit whatever the
derivation yields. A repair that iterates toward a value it was handed is the failure this method exists to
prevent — so no such value is handed to you.

---

## The item — the §3d full-field thickness-gradient content of the background-order admissibility operand

### The object (what §3d defines)
`directives/S11c_b_SHARED_PHYSICS.md` §2b and §3d define `S11CB_ADMISSIBILITY_OPERATOR_OPERAND` as the
**background-order (ε⁰) first variation, with respect to the full brane configuration, of the §3a
variable-coefficient energy-and-geometry functional written in full fields, evaluated at 𝔅⁰**, in the ordered
generalized-force pairing (bulk-DOF body force + per-face traction).

§3d fixes the gradient content precisely (quote): *"Its thickness and variable-coefficient gradient content
must be the gradient content of the **full fields**, including `∇(W_bg+δW)` and the corresponding full-field
gradients wherever a coefficient varies."* Read this as a requirement on the functional's **uniform-basis**
invariants: each §3a **uniform-basis** invariant that carries a **thickness gradient** must carry the gradient
of the **full thickness field** `∇(W_bg+δW)` — not the gradient of the thickness perturbation alone. (The
`N15` gradient-of-background invariants are a separate matter — see the repair below: they already carry a
background jet in their spurion factor and are not lifted.) §3a independently mandates that a
spatial derivative of a background profile is retained at the background-bookkeeper order of its originating
factor (a second spatial derivative of `W_bg` is first order in background amplitude and is not dropped), and
that the mixed thickness/temperature gradient invariant `∇θ·∇e_W` is one of the independent §3a invariants.

§3d additionally forbids three shortcuts, all of which still bind:
- ⛔ do **not** take the first variation of the perturbation (wave) energy and set perturbations to zero;
- ⛔ do **not** reduce it to the `ε→0` limit of the §3b wave operator;
- ⛔ do **not** insert `W_bg−W_0` into the uniform perturbation equations.
And §3a forbids obtaining the variable-coefficient energy by `W₀→W_bg(y)` substitution into the uniform `U`,
and forbids a `ρ_4D` density multiplier that §1c/§3a do not name — the density representatives enter as the
varying background fields of §2a, not as a hand-inserted weight on the energy.

### The measured non-uniformity in WL's construction (a structural fact about the engine's own code)
`constructFullFieldBackgroundEnergy` (L528–577) applies the §3d full-thickness-field gradient to the
**pure-thickness** gradient-energy invariant only: that invariant is built on `fullWidth = anchoredWidth +
WZero ewVariation` (L541) and enters as `anchoredWidth^(-2) Dot[gradient[fullWidth], gradient[fullWidth]]`
(L554). The **mixed** thickness-gradient invariant `Dot[gradTheta, gradLocalEw]` (L555) is instead built on
`gradLocalEw = gradient[localEw]` (L543) with `localEw = (WZero/anchoredWidth) ewVariation` (L540) — the
gradient of the thickness **perturbation**, not of the full thickness field. The §3d full-field requirement is
therefore applied to some thickness-gradient invariants and not to others.

### The repair (name the object; do not prescribe a coefficient or a value)
Reconstruct WL's background-order admissibility construction (`constructFullFieldBackgroundEnergy`, consumed by
`backgroundBalanceFromModel`, L1328–1351) so that the §3a **uniform-basis** invariants carrying a thickness
gradient carry the §3d **full thickness-field** gradient `∇(W_bg+δW)` — the same full-field treatment already
given to the pure-thickness invariant (L554), extended to the mixed `∇θ·∇e_W` invariant (L555), which is
currently on the perturbation-only `gradLocalEw`. ⛔ The `N15` gradient-of-background (spurion) invariants
(`newInvariantExpressions`) already carry the background jet explicitly in their spurion factor; their DOF
thickness-gradient factors are **not** lifted — a lift there introduces a **second** background jet (second
order in the background bookkeeper, outside the requested truncation) and double-counts. Re-derive from
§3a/§3d/§1c independently. Keep both anchorings and both density representatives; retain the full spatial
dependence of `W_bg`, `μ_R,bg`, `ρ_4D,bg⁰`, `ρ_br,bg⁰`, and the `Σ_E⁰` map, and every spatial derivative the
variation generates, at the §3a background-bookkeeper order. Emit the full
`ADMISSIBILITY_OPERATOR_OPERAND` (all components — U, THETA, E_W — both per-face tractions, both anchorings,
both densities), multigraded in `(ε,η,σ_W)` and dimensioned, exactly as now. Do not assert, guard, or target
any value; emit whatever the derivation yields.

### Scope — what must stay byte-identical
This item touches **only** the background-order admissibility construction (`constructFullFieldBackgroundEnergy`
and, if the fix demands it, `backgroundBalanceFromModel`). The §3b slab operator and its origins, the §3c
coupling kernel, the wave-energy construction `constructEnergyData` and everything it feeds, and the two
`kineticEw` definitions must remain byte-identical: the wave operator's mixed invariant is legitimately
perturbation-only (that operator is first order in the perturbations), so the full-field lift belongs to the
background-order admissibility functional alone. If your reconstruction moves any emitted object outside the
admissibility operand, that is a regression — stop and report it rather than absorbing it. (Whether the
pure-thickness E_W admissibility component moves is itself an output — do not force it either way; report what
moved.)

## Method and obligations
- **Independence.** Every changed object is re-derived from `directives/S11c_b_SHARED_PHYSICS.md` alone; the
  engine imports nothing; keep it so.
- **Emit discipline.** Print computed CAS objects; no `expected`/target; preserve emit tags, keys, multigrade,
  and dimension metadata.
- **Regenerate the transcript** to `mathematica/out/S11c_b_brane_operator_mathematica_audit.out`; confirm it is
  complete and non-empty (all emit tags present, no kernel error). Wrap the kernel run in a timeout; a timeout
  is a failed run to report, not to raise. The Mathematica licence has two seats — use one kernel.
- **Report:** the WL functions/lines changed, the diff, the re-derivation you performed (in Mathematica, with
  the emitted objects), and confirmation of what did and did not move elsewhere. State plainly if you could not
  complete the item or if the reconstruction forced a change outside the admissibility construction.

## Definition of done
- The §3a uniform-basis invariants carrying a thickness gradient (pure-thickness and mixed `∇θ·∇e_W`) in WL's
  background-order admissibility functional now carry the §3d full thickness-field gradient `∇(W_bg+δW)`,
  computed independently; the N15 spurion invariants are unchanged (already background-jet-carrying).
- `ADMISSIBILITY_OPERATOR_OPERAND` regenerates with all components/anchorings/densities and its metadata intact.
- The §3b operator, §3c coupling kernel, `constructEnergyData`, and both `kineticEw` definitions are
  byte-identical.
- The transcript regenerates cleanly. No object was tuned to, or checked against, any supplied value.
