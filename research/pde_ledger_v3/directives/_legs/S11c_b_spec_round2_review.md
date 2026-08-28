# Independent review — S11c-b spec round-2 fixes (§3c, §3d, §5a)

## What this is
The S11c-b SHARED PHYSICS spec was built, then BOTH blind engines (SymPy + Wolfram) were built and reviewed.
The build's four review legs found that three requirements were **under-specified** — both engines independently
implemented the WRONG thing from the same ambiguous spec text (a rule-7 shared-blind-spot). The spec's §3c,
§3d, and §5a were just re-edited (round 2) to remove the ambiguity. **Your job: confirm the re-edited text now
specifies the correct object unambiguously enough that a blind engine will implement it correctly — and that
it leaks no computed value (rule 5) and over-specifies no recipe (rule 3).**

## Artifact (read the changed sections)
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` — sections **§3c** (the
off-diagonal coupling kernel), **§3d** (the background-order balance / admissibility operand), and **§5a**
(the one-sided independence control). Read the whole spec for context but focus your verdict on these three.

## What the build legs found (the failure the round-2 text must prevent) — you are handed these
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_b_sympy_build_review.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_b_wl_build_review.md`

The three defects (all verified by FORM ablation with literal output across four independent legs):
1. **Admissibility vacuous (§3d).** Both engines built the operand as `first-variation(perturbation energy
   with fields scaled by ε) |_{ε→0}` — a bilinear energy's first variation is linear in the perturbations and
   vanishes identically at the background, so the operand ≡ 0 and the can-fail N12 test can never fail. The
   correct object (all legs that computed it agreed) is the first variation of the FULL-field background
   energy at 𝔅⁰ with the profile's own gradients retained — a genuine functional of the background.
2. **Adjointness residual tautological (§3c).** Both engines emitted `∂²U/∂u_T∂e_W − ∂²U/∂e_W∂u_T ≡ 0`
   (Clairaut — zero for any energy, tests nothing). The intended object is the adjoint relation between the two
   extracted OPERATOR blocks under the variational pairing.
3. **Off-diagonal block contaminated (§3c).** One engine built the kernel by a parallel non-operator route;
   the other applied the curl/div sector split only to UNDIFFERENTIATED field occurrences, inert on the
   operator's gradient content, leaving diagonal thickness dynamics in the "coupling" block.

## Check (report a finding only where a blind engine could still get it wrong, or a rule is violated)
1. **§3d — does the round-2 text force the full-field background functional and forbid the vacuous route
   unambiguously?** Could a blind engine still read it as "vary the perturbation energy, then set perturbations
   to zero"? Is the contrast (bilinear-vanishes-at-zero) a correct generic-math reason, and is it stated
   WITHOUT leaking the operand's value (rule 5 — the correct operand, e.g. a `−∇²W_bg`-type bending force, must
   NOT appear)? Is "genuine functional of the background / data dependence" a structural statement, not an
   asserted nonzero value?
2. **§3c adjointness — is the operator-block-under-pairing object well-defined**, and does the text correctly
   forbid the scalar-Hessian Clairaut version? Is the rule-2-corollary-3 escape ("if adjoint by construction,
   emit the blocks and say there is no independent second route, don't dress a structural zero as a check")
   correct and not itself a loophole a builder could abuse to skip the check?
3. **§3c block extraction — does "apply the sector split to the GRADIENT content" + "extract from the operator,
   not a parallel route" unambiguously prevent BOTH engines' bugs** (the inert-undifferentiated-projection and
   the parallel-route)? Or does it over-specify a recipe (rule 3) / assume a decomposition a blind engine can't
   form on a variable background (recall N5 forbids a global spectral projector — is the LOCAL curl/div split
   still well-posed here)?
4. **§5a — does "the mutation must propagate to every derived object (operator, kernel, admissibility)" fix the
   A−A tautology** without authorizing an engine to silently drop an object from the control?
5. **Any NEW leak (rule 5), new ambiguity, or contradiction with the inherited §§1–2 or the unchanged sections**
   introduced by the round-2 edits.

## Method
This is a DOCUMENT review. Form your own view of what a correct spec must say (using the finding records as the
ground truth of the target object), then read the round-2 text and compare. Quote the spec text and the source
for every finding. No CAS run required, but if you claim the §3d method still yields zero (or the §3c object is
ill-posed), show the reasoning concretely.

## Output
A short list of findings (spec quote + why a blind engine could still misimplement it or which rule it
violates + a one-line repair), or "the round-2 §3c/§3d/§5a specify the fixes unambiguously and leak no value"
with the specific ambiguities you checked and ruled out. Read-only; do not modify the tree.
