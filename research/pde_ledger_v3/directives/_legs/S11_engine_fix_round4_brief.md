# S11 SymPy engine — fix round 4 brief (radical-entry decidability and cost)

Target: `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` at HEAD (`60dc266f`).
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` wins every conflict. Change nothing beyond
what the four items require: the tag surface — names, order, record shapes, argv — is untouched.

Context, measured (two independent analyses plus two review legs on this brief; the twin carries
every command): on packages whose roots carry an explicit radical, the per-root chain's zero test is
undecided everywhere, its pivot solve uses an effectively non-terminating inverse, and its
determinant route swells — the first full-sweep run ground 8+ hours inside one cell before the user
stopped it. The spec's §Q4 confines all of it to engine territory: the counts are ranks
("independent of which basis anything returns"), the basis is "arbitrary … retained for display
only".

## Item 1 — the zero decision must decide on the engine's real minor class

`algebraic_zero_test` (`:944-979`) returns undecided for the objects the null-space scan and the §Q8
minor machinery actually feed it: rational expressions over one quadratic radical extension —
after clearing a denominator that must itself be proven nonzero, a numerator `a + b·√Q` with
`a, b, Q` polynomial in the declared symbols. Measured consequence: every candidate minor is
undecided, the early exit never fires, the scan computes all `C(·)²` minors, and the pivot loop runs
on undecided residuals — with the latent hazard that an exhausted candidate list emits an empty
basis on a regular matrix.

What must become TRUE: on that class the test decides, exactly — `a + b·√Q` is nonzero whenever
`a² − b²·Q ≠ 0` over the polynomial ring, and a vanishing norm is decided by sign under the spec's
§3 assumptions where they decide it. An expression outside the class — more than one distinct
radical base (the engine's own streams contain such objects among the coincidence-locus solutions),
or anything else undecidable — still returns undecided, honestly. The fix lands in the shared
machinery: the §Q4 scan and the §Q8 minors inherit it together.

## Item 2 — the pivot solve must complete on the measured pivot class

`build_nullspace_from_pivot` (`:1055-1081`) solves the pivot system via `selected_minor.inv()`
(`:1063`); the fallback triggers only on a raised exception, and a grinding inverse raises nothing —
it was never observed to complete on the D4 3×3 radical block.

What must become TRUE: the pivot solve completes on that measured class — the route is the
builder's choice — and produces a valid solution: its defining residual (pivot matrix times solution
minus right-hand side) reduces exactly to zero.

## Item 3 — the determinant route on radical matrices must avoid the measured swell stage

`exact_determinant` (`:1002-1011`) sends every matrix through `reduced_matrix` and then
`DomainMatrix.det()`; on radical (`EX`-domain) blocks that determinant stage is the measured swell
source — one 3×3 block took ~44 s and produced a ~5 MB object carrying hundreds of radical
occurrences, which the subsequent reduction then shrank. Raw fraction-free elimination on the same
blocks runs in seconds, and its value agrees exactly (residual zero). This route serves both
decision-only minors and the **emitted** §Q8 minors, so preserving it anywhere preserves the cost.

What must become TRUE: for matrices whose entries carry radicals, the determinant is produced by a
route that avoids the measured swell stage; the value is preserved exactly (defining residual
against the prior route reduces to zero on the measured blocks); the post-determinant reduction
remains. This covers decision-only and emitted minors alike. Emitted normal forms may move; values
may not — that consequence is judged downstream, and is deliberately not a builder-facing criterion.

## Item 4 — the §Q4 emission pipeline, at both sites, does not pay simplifier-class cost

The per-root chain runs `sp.simplify` on the m·k products, dot products and residual objects at
`:1399-1402` (`emit_q4`) **and** at `:1918-1923` (`emit_scoped_q4` — the §Q8 component-scoped Q4,
the same named-object pipeline). This step was measured at simplifier-class cost growing ~×60 per
dimension step. §Q4's residual objects are exact values and `N6` is display-only.

What must become TRUE: at both sites these objects are produced by an exact value-preserving
reduction instead of `sp.simplify`. The scope is exactly these §Q4 objects; no other emission
pipeline changes (in particular the spec's §Q5 fully-simplified ratio object is untouched).

## Acceptance — executable; operands come from the engine's real objects

The run discipline of prior rounds binds: arm `/tmp/s11_watchdog.sh` for any engine-cell run,
stream observably, never run the full sweep, experiments on `/tmp` copies.

1. Zero-test probes whose operands are real minors, extracted on a `/tmp` copy from the
   `XKIN_ANISO` D3 and D4 second-root matrices — the D4 `(1,2,3)×(1,2,3)` block among them —
   printing each operand and outcome: decided on the single-radical class, and still undecided on a
   multi-radical stream object (a coincidence-locus solution scalar). Perturbing one declared atom
   moves the printed operands.
2. Pivot probe on the D4 3×3 radical pivot block from the real cell: the solve completes within the
   run discipline and the defining residual is printed and reduces exactly to zero.
3. Determinant-route probe on the two D4 blocks the analyses measured: wall time printed; value
   residual against the prior route printed and exactly zero.
4. The `XKIN_ANISO` D3 and D4 cells complete end-to-end, wall times printed, streams flushed as
   computed; in both cells every emitted `N7` residual is zero and every `N6` residual object
   reduces exactly to zero — the identity-class validity of the bases, not a physics target. A
   driven scoped-Q4 demonstration (the established decidable-true premise drive at `:616`, on a
   `/tmp` copy) exercises the `:1918` site and prints the same identity-class checks. If D4 cannot
   complete under the observable-progress rule, that is a recorded failure of this fix — report it;
   do not narrow the objects.
5. `git diff` confined to the functions the four items name plus any helper they introduce; no
   change to tag names, emission order, record shapes, or argv handling.
