# Independent review — the S11b SymPy X-1 repair DIRECTIVE (decision list)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_sympy_x1_repair_directive.md`

This is an ORCHESTRATOR-WRITTEN repair directive. It will be handed to a builder (Codex) to repair the SymPy
engine `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py` IN PLACE. Your review is
the gate BEFORE the builder runs (CLAUDE.md rule 7 — no builder until the decision list has two legs). It is
physics-bearing: an error here makes the build wrong.

## What to check — five items, each with computation where the claim is physical

1. **Is the diagnosis correct — or is this a spec question?** The directive claims the engine over-counts the
   §5 stored-energy basis because its independence test (`independent_columns`,
   `scripts/S11b_interface_coupling_law_sympy_audit.py:406`) judges only POINTWISE polynomial independence
   over the raw field-component symbols and OMITS §5's bullet *"equivalence modulo total divergences — two
   densities differing by a total in-plane divergence are the same term; do not count both"*
   (`directives/S11b_SHARED_PHYSICS.md:286-287`).
   Read §5 of the spec (the source of truth) FIRST — the basis-construction task (~L279-281), the symmetry
   group stated in full (~L282-291, **including** both the total-divergence bullet AND the L290 clause
   *"judge independence as field bilinears, with B1's constraint NOT applied"*). Form your OWN view of what
   the symmetry-allowed quadratic basis is. Then decide: is the directive's diagnosis right (the pointwise
   test over-counts; the group's total-divergence equivalence is mandatory and reduces the basis), or is the
   spec genuinely AMBIGUOUS between "carry every field-bilinear-independent invariant" (L290) and "do not
   count two that differ by a total divergence" (L286) — in which case this is a SPEC question, not a
   one-engine repair? ⛔ Do not assume the directive is right.
   ⭐⭐ MANDATORY CAS: write your OWN script that enumerates the symmetry-allowed quadratic invariants of
   `(u, ∇u, θ, ∇θ, e_W, ∇e_W)` under the §5 group and computes the basis DIMENSION both WITH and WITHOUT the
   total-divergence quotient (the practical witness for "X is a total divergence" is that its Euler–Lagrange
   derivative vanishes identically). Save the script AND its literal stdout to named absolute paths and
   report them. ⛔ A prose derivation is discarded. Report BOTH counts; state which the spec mandates and why,
   quoting §5's exact words. ⚠ Do NOT take the directive's or anyone's asserted count as given — derive it.

2. **Does the directive LEAK a withheld value (rule 5)?** The directive is meant to direct the METHOD and
   withhold every value. Flag any sentence that states — or lets a builder iterating to a green exit infer —
   the resulting basis count, WHICH invariant is redundant, the coefficient-fold relation, or an identity
   among the invariants. Report a leak even if it looks harmless.

3. **Is the physics-preservation obligation correct and decisive?** The directive requires the redundant
   invariant be eliminated by REWRITING it (folding its coefficient into the retained invariants), NOT by
   deleting its term, and requires an `ENERGY_REEXPRESSION_RESIDUAL` whose Euler–Lagrange derivative is
   EMITTED (not asserted zero) and must be able to FAIL under a wrong fold. Check with your own CAS where you
   can: (a) is "fold, don't delete" the physically correct operation — does simply dropping the redundant
   term's coefficient lose real stiffness / change the equations of motion? (b) is the EL-derivative a sound,
   decisive witness that two densities differ by a total divergence? (c) is the two-route residual genuinely
   non-tautological (independent routes)?

4. **Is the export-chain preservation list complete?** The engine writes `S11b_exports.py`, importing S11's
   1663-row LEDGER. Is anything a basis-count edit could break missing from the carry-forward / three-valued
   F9 / digest / D3 / `_RELATIONALS` / freeze / F6 obligations the directive lists?

5. **Over- or under-scoped?** The directive deliberately does NOT demand byte-identity of the EOM/energy
   coordinates (they legitimately shift when a coefficient folds) while pinning the PHYSICS invariants
   (dispersion locus/roots, stability class, transverse decoupling, the B2c/G13 slice relation, the breathing
   form). Flag any place it demands too much (byte-preservation of a coordinate that must shift) or too little
   (a physics invariant left unpinned, so a "deleted physics" build could pass).

## What you are handed
- The directive (path above) and its rule-2 measurements twin
  `directives/_measurements/S11b_sympy_x1_repair_directive.md`.
- The spec `directives/S11b_SHARED_PHYSICS.md` (§5 ≈ L255-300).
- The engine baseline `scripts/S11b_interface_coupling_law_sympy_audit.py` (commit `864d6f41`): basis
  construction L454-535, emits L616-640, coefficient/moduli/dimension sites L498-510 / L1533 / ~L1795-1830.
- `CLAUDE.md`. You have the full repo; ⛔ there is no do-not-read list.

## Required method
This is a DOCUMENT review of a directive whose core claim is a computable physics fact. Read §5 first, derive
the basis yourself with a CAS script (item 1) BEFORE judging the directive. Copy anything you run to your own
scratch; ⛔ do not modify the working tree. Save every script and its literal stdout to named absolute paths
and report them — a physics claim without a script behind it is discarded. Physics filter: report a finding
only if it catches a way the REPAIR could produce wrong physics, leak a value, or break the export chain —
not "a builder could misread this."

## Output
For each of the 5 items: verdict + the script path and literal stdout for any physics claim + exact quotes
from the directive and the spec. End with a one-line bottom line: is the directive SAFE TO BUILD FROM as
written, or the specific edits it needs first.
