# Independent spec review — is c2's N6 rep-invariance control correctly DEFINED?

You are ONE of two independent legs reviewing an **orchestrator-written physics spec clause** for correctness
against its own parent specs. Derive your own view from the sources first; quote both sides. This is a
**definition/physics** question, ⛔ NOT a script review and ⛔ NOT a request to run anything. A real
inconsistency is a finding; there is no quota. ⛔ Do not assume the clause is correct because it is "closed."

## The object to judge
`S11c-c2 §5c` defines the **N6 "independent-route / representation-invariance" control**: which two constructed
objects it compares, and whether their residual must vanish. Determine whether that definition is **physically
correct and consistent with its parents**. If it is not, derive what the correct N6 rep-invariance test for c2
should be.

## Sources — read all, quote the load-bearing lines
- **c2 §5c** (the clause under review): `directives/S11c_c2_SHARED_PHYSICS.md` ≈ lines 303–312.
- **c2 §3c** (what the compared object, the self-energy increment, is): same file ≈ lines 203–222.
- **c2 §3d.1 / §5d** (the two density representatives `RHO4_CONSTANT` live / `RHOBR_CONSTANT` frozen): same file
  ≈ lines 228–237 and 314–319.
- **parent N6** (what N6's real control is): `directives/S11c_decisions.md` ≈ lines 88–100.
- **S11c-a §2c** (what the two anchorings LAB_HELD / MATERIAL_ADVECTED physically ARE): `directives/S11c_a_SHARED_PHYSICS.md` ≈ lines 224–244.
- The Eulerian↔material perturbation map `Δρ = δρ_E + u·∇ρ⁰` that c2 §5c invokes — find where it is defined
  (S11c-a / S11b) and confirm what it relates.

## What to determine (derive independently — do NOT accept §5c's assertion, and I am giving you no conclusion)
1. **What does N6 require in general?** Per parent N6, what two independent constructions must agree for a
   representation-invariance test, and what is the map between them?
2. **What are LAB_HELD and MATERIAL_ADVECTED, physically?** Per S11c-a §2c: are they two *representations of one
   physical state*, or two *distinct physical setups*? Quote the defining sentence.
3. **What does the field redefinition `Δρ = δρ_E + u·∇ρ⁰` relate** — two descriptions of the *same* perturbation
   (Eulerian vs material coordinates), or two *different backgrounds*? Is that the same axis as the LAB_HELD /
   MATERIAL_ADVECTED anchoring choice, or a different axis?
4. **Is c2 §5c correct?** It states "the two anchorings ARE the representation-invariance pair … map
   Eulerian↔material … their residual **must vanish** (the same operator in two representations)." Is that
   identification consistent with your answers to 1–3? Specifically: is it legitimate to treat the two
   anchorings as "the same operator in two representations," or are they distinct physics whose increments need
   not coincide?
5. **If §5c is wrong or inconsistent, what is the correct N6 rep-invariance test for c2?** Name the two objects
   that should be compared (within what — a fixed anchoring? across anchorings?), and state whether their
   residual should vanish. If a *cross-anchoring* comparison is still worth emitting, say whether it is N6 or a
   *separate* contract.
6. **Implication for a NONZERO residual:** if the current `S11CC2_REP_INVARIANCE_RESIDUAL` is comparing the
   wrong pair, does a nonzero value indicate a defect, expected physics (distinct anchorings genuinely differ),
   or an as-yet-untested question? ⛔ Do not conclude "it must vanish" merely because §5c says so.

## Method / filter
Read the sources, form your own view, quote both sides for every claim. Report a finding only if it catches a
real way c2 §5c's N6 definition is physically wrong or inconsistent with its parents. ⛔ Not wording nits.

## Output
A verdict: is c2 §5c's N6 rep-invariance pair **correctly defined** (yes / no, with the quoted parent lines that
settle it)? If no, the **correct** N6 test for c2 (the two objects, the axis, whether the residual vanishes),
and what the spec must be changed to say.
