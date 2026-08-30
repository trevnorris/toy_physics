# Decision-list review — S11c-b row-residual instrument FIX directive (coupling verdict + determinism)

You are reviewing a **fix directive** (orchestrator-written) BEFORE the builder runs it. The fix concerns a
subtle operator ordering for the coupling kernel's weak (modulo-total-in-plane-divergence) verdict at a
requested truncation. Find every way the directive would make the builder compute the **wrong** coupling
verdict, leak the answer, or specify an incoherent operator ordering. Report a numbered list (defect, the
spec/code evidence with file:line + quoted text, correction). If a point is sound, say so briefly, but a
review that finds nothing is weak evidence — probe the ordering hard. Change nothing on disk.

## Artifact under review
`research/pde_ledger_v3/directives/S11c_b_row_residual_fix_directive.md`

## Sources of truth (read and form your own view FIRST)
- The reviewed layer's classifier: `scripts/S11c_b_adjudicated_comparison.py` — `classify_total_divergence`
  (L726-767), `_euler_signature`, `transform` (L373-389, note `_bridge_d` is the LAST step and introduces the
  `(η,σ_W)` bookkeepers), and the layer's own weak-route call (L985-992, `apply_bridge_d=True` on the
  pre-bridge anchored residual).
- The spec: `directives/S11c_b_SHARED_PHYSICS.md` §3c (weak variational restriction, IBP boundary term fixed
  to zero) and §1d (variable-coefficient IBP generates first-jet terms that are physics — the reason a strong
  slab row must NOT be quotiented, but a weak coupling residual IS modulo total in-plane divergence).
- The current instrument: `scripts/S11c_b_row_residual.py` (the defect is at L304 + L406-409: it classifies
  the POST-`_bridge_d` residual with `apply_bridge_d=False`).
- Requested truncation: retain `ε^c η^a σ_W^b` iff `c≤1 ∧ a≤1 ∧ b≤1`.

## The questions (answer with quotes + file:line; reason the ordering yourself)
1. **Is the ordering coherent and correct?** The directive claims: (a) the Euler classification is only correct
   on the **pre-bridge** residual (independent background fields), (b) the requested truncation can only be
   applied to the **bridged** form (`η,σ_W` exist only post-bridge), and (c) therefore one must classify
   pre-bridge, then bridge-and-truncate the classifier's **output** to read the in-scope verdict. Verify each
   of (a),(b),(c) against the layer code and the spec. Is there a soundness gap — e.g., does bridging the
   Euler signature after computing it commute correctly (is the Euler operator's output a well-defined object
   to substitute the profile into and truncate)?
2. **The `RESIDUAL_BULK` branch.** When `_euler_signature` is nonzero, `classify_total_divergence` returns
   route `RESIDUAL_BULK` and does **not** compute a potential `V` or `anchored_remainder` (L735-739). The
   directive says to read the in-scope residual from "the bridged `anchored_remainder`, or the bridged
   Euler-bulk when the route is `RESIDUAL_BULK`." Is using the bridged-and-truncated **Euler signature** as the
   in-scope genuine-bulk representative correct and well-defined in that branch? If not, what is the correct
   in-scope object when the full residual is `RESIDUAL_BULK`?
3. **Does the in-scope truncation change the route classification incorrectly?** The route
   (`RESIDUAL_BULK` / `REPRESENTATIONAL_DIVERGENCE` / `DIVERGENCE_INCOMPLETE`) is computed on the FULL
   pre-bridge residual. The in-scope verdict is a truncation of the output. Could a case be genuine bulk at
   full order but a total divergence in-scope (or vice versa), and does the directive's object capture that
   distinction correctly, or does it conflate the full-order route with the in-scope verdict?
4. **Leak (rule 5).** Does the directive state, hint, or let the builder infer any expected route, any
   in-scope verdict (representational vs genuine), or which engine is spec-correct? Quote any leak.
5. **Determinism fix.** Is canonicalizing the emitted witness sign inside the instrument (not the layer)
   genuinely presentation-only — cannot it change a route or a residual? Is there a risk the sign convention
   interacts with a downstream zero-test?
6. **Scope discipline.** Does the directive correctly forbid routing the strong slab rows / mass row /
   admissibility through the divergence quotient (§1d), and correctly forbid modifying the committed layer?
7. **Anything else physics-bearing** that would make the coupling verdict wrong.

## Method
Document/directive review — no script to run. Ground every claim in quoted spec text and cited code lines.
Physics filter: report a finding only if it catches a way the coupling verdict would be computed wrong or the
answer leaked; not style.
