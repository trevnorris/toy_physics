# S11 SymPy engine — fix round 5 brief (radical determinants at the size cliff)

Target: `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` at HEAD (`ae105530`).
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` wins every conflict. One item; change
nothing beyond it: the tag surface — names, order, record shapes — is untouched.

Context, measured (two analysts and two review legs on this brief, all convergent; the twin carries
every command): after the round-4 fix cleared §Q4, the rerun sweep ground for hours — longer than
the entire rest of the sweep — inside `emit_q8`'s minor computation for a radical root at D4:
`all_minors` routes radical matrices through `det(method="bareiss")` (`exact_determinant`,
`:1084-1087`), whose default intermediate simplifier expands radical products that roughly square
per elimination level. Minor size 3 completes in seconds; size 4 exceeded every probe budget and
the live run's hours — the measured cliff. The class predates round 4: the pre-fix DomainMatrix
route dies on the same blocks.

## The item — radical determinants must complete at the measured size cliff, value-verified

What must become TRUE: for matrices whose entries, after the engine's reduction, are polynomial in
a single quadratic radical — the measured class of the engine's own §Q8 minor matrices — the
determinant completes at the sizes the D4 cells produce, and the value is exact: verified against
an independent exact route, not asserted. Entries outside that class — radical denominators, nested
radicals, more than one distinct radical — keep the current route through an explicit guard, and
that guard's behavior is demonstrated, not assumed. The post-determinant reduction (`:1096`-class)
remains.

⚠ The analysts measured one achievable route (adjoin `t` with `t² = R`, determinant over the
polynomial ring, reduce mod `t² − R`, substitute back) — but its measured implementation is
value-wrong when a denominator carries the radical (a review leg demonstrated this with a
counterexample), which is why the class boundary above and the acceptance controls below exist. The
route is the builder's choice; the bound, the class boundary, and independent value verification
are not.

## Acceptance — executable; operands discovered from the engine's own objects, no stated values

The run discipline of prior rounds binds: arm `/tmp/s11_watchdog.sh` for any engine-cell run,
stream observably, never run the full sweep, experiments on `/tmp` copies. ⚠ A full-sweep process
may be running (PID in the session record); do not touch it or its output file beyond read-only.

1. Determinant probe with operands discovered dynamically: reconstruct, on a `/tmp` copy, every
   radical-root minor matrix family the current stream's D3 and D4 cells define (derive the ranks
   and minor index sets at runtime from the engine's own functions — state no counts and no
   values). Every determinant completes with wall time printed, and every value — including every
   size-4 block — is verified by BOTH of two independent checks, residuals printed: (a) exact
   rational specialization at several sample points satisfying the declared assumptions, evaluated
   through the new route and through an exact reference evaluation of the original matrix; (b) a
   structurally different exact symbolic route (e.g. cofactor expansion whose cofactors go through
   the prior route at the size below the cliff). A disagreement anywhere is the finding.
2. Class-boundary controls, outcomes printed: constructed entries with a radical denominator, a
   negative odd half-power, a nested radical, and two distinct radicals — each either handled with
   value verified by the specialization check, or routed to the prior route by the guard.
3. The same dynamic probe at D3: values agree with the prior route on every minor — residual zero —
   so the change is invisible below the cliff.
4. The `XKIN_ANISO` D2 cell completes end-to-end in-tree, wall time printed. ⚠ Do not use D3/D4
   cell completion as acceptance: those cells contain fresh-process-fragile locus-solve stages
   outside this item's scope, measured and recorded in round 4.
5. `git diff` confined to `exact_determinant` and any helper it introduces; no change to tag names,
   emission order, or record shapes.

A green run claims exactly one thing: the determinant wall is repaired with values verified. It
does not claim D4 cell completion — the radical locus-solve class after this stage is measured,
deferred, and recorded in the twin.
