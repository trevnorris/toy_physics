# S11c-a — WL engine fix: ground the T-e shifted density trace in the supplied background

## Target
`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (the blind Wolfram
engine). Physics authority: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`. This is a
single-object correction; do not touch any other emitted family's physics.

## The defect to correct (name the object, not the recipe)
The engine's **T-e shifted-face trace** (`S11CA_FACE_SHIFT`) obtains the **bulk-density background** from
an object that is **not a member of the supplied background state `𝔅⁰`**. Concretely, the density entry
of the trace inventory resolves the background to a bare, **undefined abstract function of the normal
coordinate `w`** (it has no definition anywhere in the engine), and the shifted-trace law then carries
that function's normal derivative(s) as live symbols.

This violates **§3c**, which requires: *"Every background face value or normal derivative appearing in
this law is obtained by differentiating a member of the supplied background state `𝔅⁰` (§2d); none may be
introduced as a free premise,"* and which states the supplied density background *"depends on the
in-plane anchor, not on `w`."* The two supplied density representatives are the members to differentiate
(**§2b**: `ρ_4D,bg⁰` is `ρ_4D,ref⁰` or `rho_br/W_bg(y)`; both functions of the in-plane anchor `y`), and
they are named members of `𝔅⁰` (**§2d**).

Every **other** density-bearing family in this same engine already grounds its background correctly by
differentiating the §2b representative (the engine already carries the grounded, representative-aware
in-plane density profile and uses it throughout the projection, evolution, traction, virtual-constraint,
and background-density-map objects). The shifted trace is the **one** place that reaches for the
ungrounded free premise instead. The pressure, velocity, and current backgrounds in the same trace are
already grounded (their background is zero by §2d/§3b) and must be left exactly as they are.

## What the corrected object must satisfy
1. **The T-e density trace's background — its face value and its normal derivative — must be obtained by
   differentiating the supplied §2b density representative** (the same grounded, in-plane, w-independent,
   representative-aware background the engine already uses for its other density-bearing objects), **not
   from any free-premise or w-argumented abstract function.** After the correction, no `S11CA_FACE_SHIFT`
   operand may contain a normal derivative of an undefined background function; every background normal
   derivative that appears in the shifted-trace law must be the derivative of a supplied §2b member.
2. **Because that background is parameterised by the density representative, the `S11CA_FACE_SHIFT` family
   must be keyed by the density representative** (the DENSITY axis), the same way the other
   density-bearing families in this engine are keyed. Emit one case per representative for the traced
   fields, so the family key is uniform across fields. (The representative-independent fields — pressure,
   velocity, current — are simply emitted once per representative; this is a keying change for them, not
   a physics change.)
3. The **perturbation** part of every traced field is unchanged: the perturbation is evaluated at the
   background face `h_s⁰ = s·W_bg/2` carrying `η` through `W_bg`, retaining first-shape-order dependence
   (§3c) — exactly as it is today. This fix touches the **background** of the density trace only.
4. **This applies to EVERY place the engine builds a shifted-face trace, not only the primary
   `S11CA_FACE_SHIFT` emission.** The same trace object is rebuilt independently for the source-level
   form-ablation control and for the uniform-limit regression. **All** of these — the primary shifted
   trace, the form-control shifted-trace slice, and the uniform-limit shifted-trace slice — must (a) draw
   the density background from the supplied §2b representative and (b) carry the density-representative
   axis. They are the expected co-changers of the same T-e object; leaving any of them on the ungrounded
   or density-agnostic path leaves the defect alive in that package and mismatches the sibling engine's
   keying. Do not fix one site and leave the others.

## Script obligations (every one is mandatory)
1. **The engine may PRINT computed objects; it may NOT state conclusions.** Every emitted payload is a CAS
   object (an expression, a jet, a boolean from a symbolic test) — never a prose sentence describing a
   result.
2. **PRINT operands; do not assert them.** Do not add any `Assert`/`If[...==0]` gate on a measured
   density-trace payload. The residual against the sibling engine is computed on the reviewer's side, not
   here; a genuine disagreement is a finding, not a build failure.
3. **The ONLY place the physical symbols may be combined by hand is in constructing the background state
   and the ansatz. Every other expression involving them must be reached by computation; every control
   re-enters the chain at the supplied background, never at a result.** In particular, the density
   background normal derivative must be produced by an actual `D[…, w]` on the supplied representative,
   never written in by hand; whatever it evaluates to must *emerge* from the differentiation of the
   supplied member, not be asserted.
4. Keep the engine **blind**: it imports nothing and re-derives from the spec. This fix changes only which
   supplied object the density background is drawn from; it does not import or transcribe anything.

## Acceptance (compute-and-print; the diff is on our side)
- Re-run the engine to regenerate its committed transcript
  `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`. The shifted
  density trace in **every** package that carries it — the primary `S11CA_FACE_SHIFT`, the form-ablation
  control's shifted-trace slice, and the uniform-limit regression's shifted-trace slice — must now carry a
  density-representative axis in its keys, and the bulk-density operands must be expressed in terms of the
  supplied §2b background (and the density perturbation), with no normal derivative of an undefined
  free-premise function anywhere in any of them.
- **Form control (must bite):** in a throwaway copy, replace the grounded density background with an
  explicitly `w`-dependent test function and re-run; the emitted density-trace shift term must change.
  Restore. This shows the shift term is genuinely wired to the background you grounded it in — report the
  literal before/after of one density operand. Do **not** report or assert any expected value or residual.
- The expected co-changers are exactly the shifted-trace-bearing packages named above (primary
  FACE_SHIFT, its form-ablation control slice, its uniform-limit slice). Every other emitted family
  (normal, measure, velocity, flux, traction, projection, evolution, closure, kinematic, conormal,
  virtual work, and any control that does not rebuild the shifted trace) must be **unchanged** except for
  any that legitimately consume the density trace; confirm by diffing the regenerated transcript against
  the committed one and reporting which tag families changed and why each change was expected.

## Out of scope
Do not alter the spec, the PY engine, the comparator, or any other family's physics. Do not add an
acceptance gate that references an expected value. Do not "simplify" the density representatives into one
another (§2b forbids it before T-g/T-h/T-i).
