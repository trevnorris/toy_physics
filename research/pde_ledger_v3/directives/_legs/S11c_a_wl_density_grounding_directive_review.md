# Independent review — WL density-grounding fix DIRECTIVE (decision-list gate)

You are reviewing a **build directive** before any builder runs it. The directive proposes a
single-object correction to the blind Wolfram engine of S11c-a. Your job is to find defects in the
**directive** — not to perform the fix.

## Artifact under review
`research/pde_ledger_v3/directives/S11c_a_wl_density_grounding_fix_directive.md`

## Sources of truth — read these and form your own view BEFORE judging the directive
- Spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — especially §1b, §2b, §2d, §3c,
  and the T-e line (~§430).
- The WL engine: `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` —
  especially how it builds the shifted-face trace for the bulk density vs. how it builds density
  backgrounds for its other families (projection, evolution, traction, virtual constraint, background
  density map). Grep the density-background function the trace uses and check whether it has a
  definition anywhere in the file.
- The PY engine for context only: `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`.

## What to check (answer each, with file:line and quotes)
1. **Is the diagnosis correct?** Derive the §3c shifted density trace yourself from the supplied
   background state. Is the WL density background actually an ungrounded free premise (an undefined,
   w-argumented abstract function), and does §3c actually forbid that? Is the value of the background
   normal derivative forced by the spec — and if so, is the directive right that it must *emerge* from
   differentiating a §2b member rather than be asserted? If you can compute the forced value, show a
   runnable snippet and its literal stdout (do not just assert it).
2. **Object vs recipe.** Does the directive name the object/property to achieve, or does it over-specify
   *how* to edit the engine (which would be a rule-3 defect)? Flag either failure: over-specification, or
   under-specification that leaves the builder guessing.
3. **Leakage.** Does the directive leak an expected value, an expected residual, or the sibling engine's
   answer that a builder could iterate toward? The acceptance must be value-free and the diff must be
   deferred to the reviewer's side. Quote anything that leaks.
4. **Able-to-fail acceptance.** Is the form control ("replace the grounded background with a w-dependent
   test function; the shift term must change") a genuine, decisive test — would it actually fail if the
   builder hardcoded the result or left the background ungrounded? Is there any acceptance clause whose
   only honest outcome is a fabricated one?
5. **Keying change.** The directive requires the FACE_SHIFT family to gain a density-representative axis.
   Is that correct and consistent with how the engine keys its other density-bearing families? Could it
   silently corrupt the representation-independent fields (pressure/velocity/current) or any downstream
   consumer of the trace?
6. **Scope / collateral.** Could following this directive change the physics of any *other* emitted family
   (normal, measure, velocity, flux, traction, projection, evolution, closure, kinematic, conormal,
   virtual work, controls, uniform-limit)? Is the "out of scope" section sufficient, or does it miss a way
   the change could leak into another object?
7. **Blindness.** Does the fix preserve the engine's blindness (imports nothing, re-derives from spec)?
   Does it anywhere instruct transcribing the sibling engine?

## Method
- A prose re-derivation is worth nothing on its own: where you assert a computed value, paste a runnable
  snippet (SymPy or Wolfram) and its literal stdout with an absolute path.
- Report a finding only if it changes what the directive should say or catches a way the resulting build
  could be wrong. "It could be clearer" is not a finding unless the ambiguity could produce a wrong build.
- If the directive is correct and buildable as written, say so plainly and name the one or two things you
  checked hardest.

## Return
A short list of directive defects (or "buildable as written"), each with file:line, the quoted text, the
spec/source it contradicts or the ambiguity it creates, and the concrete correction.
