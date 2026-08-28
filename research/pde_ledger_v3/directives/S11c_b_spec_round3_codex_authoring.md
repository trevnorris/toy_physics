# S11c-b spec round-3 — Codex authoring directive (§3a jet-clause, §3c, §3d, §5a)

## Task and authority
Edit `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` in place, rewriting ONLY four spots to fix
three verified review findings. This is a rule-15 author change: the orchestrator has revised §3c/§3d twice
and they remain ambiguous, so you (the implementer) rewrite them precisely. `CLAUDE.md` binds.

⛔ **Edit only these four spots; leave every other byte of the file unchanged** (do not reflow, renumber, or
touch §0–§2, §3b, §4, §5b/§5c, §6–§8, or any tag name). Report the `git diff` of the file when done.

## The findings to fix (full detail + exact repairs)
Read `research/pde_ledger_v3/directives/_measurements/S11c_b_spec_round2_review.md` (R1/R2/R3) and the two
build-review records it cites (`_measurements/S11c_b_{sympy,wl}_build_review.md`, B1/W1 etc.). Those are the
ground truth for the target objects. Implement their repairs:

**A. §3d (the admissibility operator operand) — fix R1.** Rewrite so a blind engine cannot produce the
identically-zero operand both engines shipped:
- The operator operand is the first variation of the **§3a variable-coefficient (full-field) energy** at `𝔅⁰`
  — NOT the uniform §1c form and NOT the perturbation (wave) energy. The thickness and coefficient **gradient
  content in that functional must be that of the FULL fields** (`∇(W_bg+δW)`, etc.), so the background
  profile's own gradients are present at background order; ⛔ `|∇(δW)|²` (perturbation-only gradient) is not an
  acceptable background-order representative.
- The first variation is taken wrt the brane configuration and evaluated at `𝔅⁰` (`θ=e_W=δW=0`, `u=0`; profile
  and its jets retained). Keep the existing correct ⛔s (do not vary the perturbation/wave energy then set
  perturbations to zero — bilinear ⇒ vanishes identically; do not reduce to the `ε→0` limit of the §3b wave
  operator; do not insert `W_bg−W_0`). Keep the same pairing as `𝒮_hold⁰` and the "value is not stated here /
  is the computed result" clause.

**B. §3a jet clause — fix R1 flip side (Codex F2).** The clause "the divergence-form operator needs only first
coefficient jets; no higher jet is introduced" must not make an engine discard the higher spatial derivatives
the variation legitimately generates. State value-freely: **"first-background-jet order" bounds the number of
independent background AMPLITUDE factors (powers of `η`/`σ_W`), NOT the spatial-derivative order**; a
divergence or a variation may generate a higher spatial derivative of a background profile (e.g. from
`∇·(a∇·)`), and that term is retained at the background-bookkeeper order of its originating factor (a second
spatial derivative of `W_bg` is still first order in the background amplitude, `O(η)` / `O(σ_W)`), not dropped.

**C. §3c (the off-diagonal block + adjointness) — fix R2.** Keep the correct forbids (⛔ the scalar-Hessian
Clairaut adjointness `∂²U/∂u_T∂e_W − ∂²U/∂e_W∂u_T ≡ 0`; ⛔ the parallel non-operator extraction route; ⛔ the
undifferentiated-only projection that is inert on gradient content). Add the missing well-posedness:
- Define the off-diagonal block extraction AND the adjoint residual as a **weak variational restriction using
  independent transverse (solenoidal, `∇·u_T=0`) and longitudinal (irrotational, `∇×u_L=0`) TRIAL and TEST
  displacements** paired under the §1c variational form. This attributes BOTH the gradient-structured terms
  AND the undifferentiated-`u` spurion couplings (`g·u`, the §1a/N15 channels) to a sector without a global
  projector (N5-safe).
- State the in-plane **domain / boundary convention** the pairing uses (compact support, decay at infinity,
  periodicity, or explicit in-plane boundary terms) so the pairing's boundary term is defined and the adjoint
  residual is well-posed. Keep the rule-2-corollary-3 escape (if the two blocks are adjoint by construction,
  emit the blocks and state there is no independent second route rather than dressing a structural zero).

**D. §5a (the one-sided independence control) — fix R3.** The omit-escape must be gated on structure, not
prose: an object is omitted from this control **only when the mutated source is STRUCTURALLY ABSENT from that
object's construction (shown by the computed/emitted dependence, not an asserted reason)**; otherwise the
mutation must propagate into and move that object's operand. Keep the A−A ban.

## Rules that bind the rewrite (non-negotiable)
- ⛔ **Rule 5 — leak no computed value.** Do NOT write the correct operand (`−∇²W_bg`, `−κW_bg″`, any
  coefficient/sign/grade), the coupling grade/sign, the new-invariant forms, the admissibility residual value,
  or the basis count. Name the OBJECT and the METHOD; the engines compute the value. Stating the generic-math
  fact that the *forbidden* construction vanishes identically is allowed (it leaks no S11c-b value).
- ⛔ **Rule 3 — name the object, not a unique CAS recipe.** The weak-restriction and full-field-functional
  language must define WHAT to compute, not prescribe one implementation.
- Cohere with the rest of the spec: the `(ε,η,σ_W)` multigrade, both anchorings, both density
  representatives, the `S11CB_*` tag names, and the §2b reference to §3d/§3a all stay consistent.
- ⛔ No `VERDICT`/expected value/acceptance criterion anywhere.

## Report
The `git diff` of `S11c_b_SHARED_PHYSICS.md`, and a ≤10-line note stating which finding each edit fixes and
confirming no other section changed.
