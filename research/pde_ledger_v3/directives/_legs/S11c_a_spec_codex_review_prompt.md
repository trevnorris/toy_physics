# Independent physics review — the S11c-a spec (Codex-authored replacement)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`

⚠ This spec was **authored by Codex** to replace an orchestrator draft whose own revision kept breeding
defects (CLAUDE.md rule 15). ⛔ Do not assume it is correct because a new author wrote it. Derive your own
view and let it disagree with the artifact.

## What it is, and why this review is load-bearing
S11c is the non-uniform transverse coupling at a brane–bulk interface — *is light's confinement
unconditional?* The prior step **S11b** (uniform) is CLOSED and proved the uniform transverse coupling is
**identically zero**; S11c computes the gradient-driven channel. S11c is a five-sub-step family; **S11c-a**
is the first: it un-freezes the in-plane background (`W₀→W₀(x)`, `μ_R→μ_R(x)`) and computes the first-order
**shape-derivative of every interface law** S11b wrote for a flat wall. Two blind engines (SymPy + Wolfram)
read this spec; an error here corrupts both. ⛔ This is a **physics review**, not a copy-edit.

## ⭐⭐ THE PRIMARY TARGET — FIDELITY OF THE INLINED S11b SETUP (§1)
Unlike a citing spec, this one **inlines** the S11b setup as supplied equations (§1a–§1c). ⚠ **A single
transcription error in an inlined equation makes both engines agree on the same wrong thing** — the worst
failure mode. ⭐ **Verify EVERY supplied equation in §1 term-by-term against
`research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`:** the relative flux `J_s` (outward sign once, per
face area), the closure `J_s=Λ_A𝒜_s+Λ_V V_s`, the affinity `𝒜_s=μ_s−δp_s/ρ_m`, `μ_s=μ_θ/ρ_br⁰`, the
**held-fixed** `μ_θ` qualifier, the traction `t_s`, the kinematic balance `n̂_s·v_bulk=V_s+J_s/ρ_m`, the
material virtual constraint (`δ_vΣ_mat=0`, `δ_vx=δ_vu`, `δ_vΣ_E=0` forbidden), the physical evolution law
`∂_tΣ+∇·(Σv)=−(J₊+J₋)` with `v≡∂_tu`, the energy `U`, the sign conventions, the two face DOFs `δW` and
`ζ_c`, and the drain. Report any sign, factor, held-fixed-qualifier, or normalization discrepancy.
⚠ The rule-2 twin `directives/_measurements/S11c_a_SHARED_PHYSICS.md` carries the S11b sources — use it to
locate each, but verify against the S11b spec itself.

## Derive your own view FIRST, then read the spec
Before opening the spec, form your own answer to *what a correct "background & geometry" sub-step must set
up*, from:
- `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` — the closed uniform spec (§1 geometry `:63-160`,
  §3 flux `:183-207`, §4 closure/affinity `:208-254`, §5 energy `:255-300`, §6 method/constraint/signs/exit
  `:301-380,:499-546`, B0/B1 bulk response + mass balance `:593-646`, B6 `:819-831`).
- `research/pde_ledger_v3/directives/S11c_decisions.md` — the folded decision list; S11c-a must honor N3
  (FULL shape-derivative), N4 (anchoring is the physics), N6/N11b/N12/N14/N15.
- `research/pde_ledger_v3/V3_STEP_PLAN.md:1179` — coupling `∝ k·a`, `∇μ_R≠0` mixes; FORM now, MAGNITUDE
  needs the interior `R1`.
- `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` — the closed uniform result.

## Required method — DOCUMENT branch
Read the sources first, form your own view, then read the spec. Quote **both sides** for every finding.
⚠ Reading order is a method instruction, not a blindness control.

## What to check (beyond §1 fidelity)
1. **Completeness of the shape-derivatives.** Is any interface law whose shape-derivative changes at
   `O(η)`/`O(slope)`/`O(ε·slope)` **missing** from the computed objects (§4)? The set should cover: outward
   normal, conormal `n̂·∇₄`, face velocity, relative flux, tilted kinematic mass balance, traction + its FULL
   virtual work (`t_s·δ_vx_s` from the face map, both anchoring branches), the physical evolution mass
   balance, the assembled face-closure shape-derivative (incl. the affinity-normalization from a varying
   `ρ_br⁰(x)`), the shifted trace, the dynamic-window projection. The **true-face-area measure** must appear
   wherever a background face flux/pressure is nonzero.
2. **⛔ Pre-stated results (rule 3).** Does any task **state** a derived term, its order, whether a channel
   survives/cancels, or a coefficient — instead of naming the object and letting it be computed? (The prior
   draft's defects were exactly typed results: "carries `t_∥·δ_v u`", "`n_∥=O(η)`".) The anchoring-dependent
   in-plane work, and every `(ε,η,slope)` order, must be **computed**, not asserted.
3. **Anchoring & face maps.** Are the two branches (lab-held, material-advected) each a single, unambiguous
   parametric face map `R_s(X,t)` from which `V_{n,s}` and `δ_vx_s` are both computed? Is `δW`/`ζ_c` carried
   (⛔ not frozen to `ζ_c=0`)?
4. **Slope vs contrast.** Are the contrast `η` and an independent slope scale (`L_W`/`σ_W`) genuinely
   distinct and internally consistent, or does the ansatz secretly tie them?
5. **Admissibility (2d).** Is it a supplied premise + named support with the can-fail residual reserved to
   S11c-b (⛔ not an identically-zero residual here, ⛔ not mislabeled "tested")?
6. **Binding across the split.** Are `J_s`/`V_s` bound to one definition across the flux, kinematic balance,
   evolution law, and closure (⛔ so no task can use the flat Cartesian `J` where the tilted one is meant)?
7. **Under-specification.** Any binding relation in prose that must be an equation (the `∇·u`-fell-out-four-
   times failure class)? Any object named but its defining map/scale not actually supplied?
8. **Scope & controls.** Does it defer S11c-b…e correctly (export geometric boundary operators, ⛔ not solve
   the DtN, ⛔ no global dispersion relation)? Are the form control (source-level, per-direction ablation)
   and the two-route/one-sided-corruption control genuine and able to fail? Emit-both-operands-plus-residual
   (⛔ never assert zero)? No terminal VERDICT; physics disagreement exits 0?

## Physics filter
Report a finding only if it catches a way the physics or the set-up could be wrong — not stylistic
preference, not "wrong on a different input." ⛔ A leg returning "nothing survives" is weak evidence: state
what you checked, what you derived, and where you looked hardest (esp. which §1 equations you verified
term-by-term).

## Output
1. Findings, most-serious first — each with the source quote (S11b or decision list), the spec line, and the
   concrete way the physics/set-up goes wrong if uncorrected. ⭐ Put any **§1 fidelity** discrepancy first.
2. Is any interface-law shape-derivative **missing**?
3. Anything **under-specified** (prose where an equation is needed) or **over-specified** (a result form or
   order pre-stated)?
4. Does the scope correctly defer S11c-b…e, and is the geometric-boundary seam to S11c-c right?
