# Directive — RE-AUTHOR the S11c-a spec (background & geometry)

**You (Codex) are the AUTHOR.** Rule 15: the orchestrator's fold of this spec bred defects, so a different
author writes the corrected version. Your output is reviewed by a fresh Claude agent + Grok — ⛔ you do not
review it. This is a **physics specification** (an obligation-to-compute document), ⛔ not a script.

## Deliverable
**Overwrite** `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` with the corrected spec. It is the
physics authority two blind engines (a SymPy engine that imports the S11b LEDGER, and a Wolfram engine that
imports nothing and re-derives) will read for **S11c-a**, the first sub-step of the S11c family.

## Read first (the authority)
- `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — the current (defective) spec at HEAD
  `82d14079`. It is your BASE: **keep everything it got right** (below), fix everything the punch-list names.
- `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` — the closed uniform spec S11c-a builds on.
- `research/pde_ledger_v3/directives/S11c_decisions.md` — the folded decision list (`N1`–`N15`); S11c-a must
  honor N3 (full tilted-face shape-derivative), N4 (anchoring is the physics), N6/N14/N15, N11b, N12.
- `research/pde_ledger_v3/V3_STEP_PLAN.md:1179` — the coupling is `∝ k·a`, `∇μ_R≠0` mixes; FORM now,
  MAGNITUDE needs the interior (`R1`).

## ⭐⭐⭐ THE META-RULE THAT BRED THE LAST DEFECTS — obligation-to-compute, ⛔ NOT typed results
The fold's new defects were **typed derived results** (e.g. "the virtual work now carries `t_∥·δ_v u`",
"`n_∥ = O(η)`"). ⛔⛔ **Do NOT state a derived term, its order, whether a channel survives, or a coefficient.**
⭐ **Name the object and supply the setup (the geometry, the DOFs, the laws); every derived form, term, order,
and surviving-or-cancelling channel is REACHED BY COMPUTATION** by the engines and checked by the reviewers.
This is CLAUDE.md rule 3 and `S11b_SHARED_PHYSICS.md:499` (the structural rule). The one place physical
symbols are combined by hand is the §1 setup and the §2 background ansatz; ⛔ everywhere else names an object.

## KEEP (both 2nd-pass legs confirmed these are correct — ⛔ do not regress)
- The **outward-orientation law** `s(n̂_s·ŵ) > 0` fixing both faces outward (⛔ not the bare gradient sign).
- **T-g**: linearise the **material** `δ_v Σ_mat = 0` with the **virtual** `δ_v u` (⛔ not physical `u`),
  dimensionally consistent (⛔ no `Σ`-dimension addend to the dimensionless constraint), preserving the `e_W`
  redefinition that survives even where `∇Σ_E⁰=0`; ⛔ form not pre-stated.
- **μ_R(x) varies** at order of the background inhomogeneity, with anchoring; the operator is S11c-b's (N15).
- **Admissibility seam**: S11c-a names the background state + declared support; the can-fail residual is
  S11c-b's (⛔ not an identically-zero residual of S11b's uniform equations).
- **Two anchoring branches** named as computable maps; **F9 name reservations**; **no VERDICT / physics
  disagreement exits 0 / emit-both-operands-plus-residual**; the scope deferrals to S11c-b…e.

## THE PUNCH-LIST — 12 corrections (both review legs, verified against source)

**1 · Encode contrast vs slope (both legs).** §2e declares contrast `η` and slope `|∇W_bg|` independent, but
the ansatz `W_bg=W̄₀(1+η w₁(x))` forces `∇W_bg=η W̄₀∇w₁`, tying them. ⇒ supply an **independent profile
scale** `L_W`: `w₁ = w₁(x/L_W)`, with slope order `σ_W ≡ η W̄₀/L_W` (or a formal derivative-order map). ⛔
Then **do not assign** tilt/geometric objects (the in-plane normal, tilted traction, `∂_tΩ`) an order by
hand — the engine computes each object's `(ε, η, σ_W)` order. (The `e_W` redefinition and `μ_s=μ_θ/ρ_br⁰(x)`
are contrast-`η`; the tilt is slope-`σ_W`; ⛔ but state that as the parameter *definitions*, not as
pre-assigned per-object orders.)

**2 · T-d as a recomputed object, ⛔ not a typed pairing (both legs).** Restate T-d: "compute the
shape-derivative of `δ_v𝒲_bulk = Σ_s t_s·δ_v x_s`, with `δ_v x_s` taken from the §3 face map **per
anchoring**." ⛔ Do **not** assert an in-plane `t_∥·δ_v u` term exists — it is anchoring-dependent (present
in the natural material map, cancels in the natural lab-held map). ⭐ Supply the parametric face maps so
`δ_v x_s` is defined (see item 3); which pairings survive is computed.

**3 · Give the parametric face maps `R_s(X,t)`, not just level sets (both legs).** Level sets `F_s=0` fix the
geometric normal but ⛔ not the material virtual displacement `δ_v x_s` that virtual work needs. Supply, per
anchoring branch, the face position map `R_s(X,t)` (lab-held vs material-advected) from which both `V_{n,s}`
(T-b) and `δ_v x_s` (T-d) are computed. ⛔ "or `R_s`" wording that lets two engines pick inequivalent
parameterisations must be closed to a single supplied definition each.

**4 · Import BOTH independent face DOFs — ⛔ do not freeze `ζ_c` (Codex).** S11b's conventions carry two
independent combinations: `δW ≡ ζ₊−ζ₋` **and** the centre shift `ζ_c ≡ (ζ₊+ζ₋)/2` (`S11b_SHARED_PHYSICS.md:89`;
B0b reports `Z` separately for both). ⛔ Importing only `ζ_±=±δW/2` hard-wires `ζ_c=0` and removes the
opposite-parity channel S11c-c's two-face DtN needs. ⇒ carry both DOFs into every T-item's shape-derivative.

**5 · Supply T-h's velocity `v=∂_t u` (Codex).** S11b's physical mass balance defines the slab material
velocity `v=∂_t u` (`S11b_SHARED_PHYSICS.md:641`; the balance is at `:644-646`). ⛔ The current spec imports
the balance with `v` undefined, while using unqualified `v` for the bulk velocity in T-c. Supply `v=∂_t u`
for T-h; the linearisation's gradient term `∇Σ_E⁰·∂_t u` is then computed, not dropped.

**6 · Surface measure — `J_s` is per FACE area (Codex).** `J_s` is flux per unit **true face area**
(`S11b_SHARED_PHYSICS.md:210`); a tilted face carries the measure `a_s=√(1+|∇h_s|²)`, so the physical balance
(T-h) and the virtual work (T-d) carry `a_s J_s` / `a_s t_s·δ_v x_s`, whose shape-derivative activates when a
**background** face flux/pressure `J_s⁰, p_s⁰` is nonzero (which §2d permits). ⇒ **either** state
`p_s⁰=J_s⁰=0` and exclude background face gradients explicitly, **or** supply the surface-measure setup and
require its shape-derivative. ⛔ Do not leave the background face fields permitted but their measure terms
uncomputed. Also require the shifted-trace `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰` when a background face
field is nonzero.

**7 · Affinity import is not verbatim (Codex).** S11b supplies `μ_θ ≡ (δU/δθ)|_{u, e_W, all other fields
fixed}` (`S11b_SHARED_PHYSICS.md:224`). ⛔ The current spec abbreviates to `μ_θ ≡ δU/δθ`, dropping the
held-fixed qualifier — an engine could eliminate `θ` via the constraint first and get a different chemical
potential (double-counting in T-i). Supply the **full** qualifier in §1.

**8 · §10 mislabels a premise "tested" (Codex).** §10 classifies "the admissibility premise/state" as
**tested** in S11c-a; it is **supplied** (the can-fail residual is S11c-b's). Reclassify as supplied, and
reserve the S11c-b residual explicitly (name its tag and operands).

**9 · Bind `J`/`V` across the split T-items (Grok).** S11b used one `J` and one `V` everywhere; split across
T-b/T-c/T-c′/T-h/T-i, an engine can keep the flat Cartesian `J_±=ρ_m(v_w−∂_tζ_±)(±1)` on the mass-balance
RHS while T-c carries the tilted flux — silently dropping the coupling image (`:821-823` forbids dropping the
transverse channel by a divergence-free argument; the flat-`J` bookkeeping does it too). ⇒ bind T-h's/T-i's
`J_s`, `V_s` explicitly to T-c's tilted flux and T-b's `V_{n,s}`.

**10 · Bind the ansatz constants to the imported keys (Grok).** `W̄₀`, `μ̄_R` (the constant references in
the §2 ansätze) must **be** the imported constants `W_0`, `mu_R`, so `(η,σ_W)→0` lands on S11b's objects.
(N14 forbids reusing those keys for the *varying* profiles — that ⛔ stays; only the constant references bind.)

**11 · C-1 must ablate separately and police T-d/T-i (both legs).** C-1 is a FORM control that must re-enter
at the **supplied level set / law** (`:502-504`), ⛔ not at an already-emitted `n̂`. Ablate the supplied
`F_s` slope components **separately** (one diff per in-plane direction, `D_brane=3`), ⛔ not one simultaneous
drop (which only proves sensitivity to some slope term). And the control set must cover the channel T-d adds
(the two-route §4 list currently omits T-d).

**12 · T-i "assembled-closure" wording (Grok).** T-i's law is `J_s − Λ_A𝒜_s − Λ_V V_s = 0`; the word
"assembled" must not be read as B0c's `δp = Z·v_bulk` (which S11c-c owns). Clarify T-i is the face-closure
shape-derivative only, ⛔ not the bulk-response assembly.

## Authoring discipline (inherit, and state in the spec)
- Supply the full S11b setup in §1 (don't blind the inputs); flag §1 as **supplied/unfalsifiable within this
  build**. There is **no acceptance value to withhold** for S11c-a (the bench-optics bound is S11c-e).
- Balance laws (⛔ not an action principle); the virtual-displacement rule; variational (⛔ not
  ordinary-partial) derivatives; sign conventions; the three script clauses; **no terminal VERDICT**; a
  physics disagreement **exits 0** (`:518,:543`); controls **emit both operands and the residual**, ⛔ never
  assert it zero.
- Parallel tag sets, one tag per named object; ⛔ do not let the two engines pick two names for one object.
- ⛔ Do not defer the geometric boundary operators (normal, conormal `n̂·∇₄`, tilted kinematic BC, face
  shift) — those are S11c-a's exports; S11c-c *solves* the curved DtN. Do not request a global dispersion
  relation (N5).

## Report back (⛔ under 30 lines)
The corrected file path; a short list of which punch-list items you applied and how; and any place the
punch-list itself looked wrong (you are the author — flag it, ⛔ do not silently deviate). ⛔ No conclusions
about the physics; the spec computes, it does not conclude.
