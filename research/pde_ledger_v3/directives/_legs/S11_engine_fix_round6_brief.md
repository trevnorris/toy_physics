# S11 SymPy engine — fix round 6 brief (KW zero-locus solve at the radical class)

Target: `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` — the working-tree
copy, content-identical to commit `9e392206` (the repo tip is newer; the engine file is unchanged).
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` wins every conflict. One item; change
nothing beyond it: the tag surface — names, order, record shapes — is untouched, and `PACKAGE_DIMS`
/ driver order are not the builder's to touch.

Context, measured (two analysts, convergent, re-verified by the orchestrator, then two review legs
on this brief; the twin carries the commands): a full sweep on this engine has been silent 4+ hours
inside `emit_locus`'s non-canonical solve (`:840`, `sp.solve(residuals, list(variables),
dict=True, simplify=False, manual=True)`), reached when `canonical_locus` returns NOT_APPLICABLE —
which it does when a residual's radical involves the solve variables (a constant radical still
takes the canonical branch). The stuck residual's radical-free denominator is free of the radical;
after clearing it, the numerator is polynomial in the engine's symbols and a single quadratic
radical whose radicand excludes two of the unknowns. A 120 s / 840 s profile puts ~99% of the time
in `checksol`/`expand` under `_vsolve`: the system loop solves its first symbol in under a second,
then spends the hours on per-symbol solves of the radicand symbols — every one of which its own
discard rule then eliminates as dependent on the already-solved symbol. The same hang reproduces
from a cold process at D3 (whose cell the sweep DID complete warm), so the wall is the route, not
the operand. ⚠ `emit_locus` serves every locus family in the stream, not only KW — the edit
surface is shared, and the acceptance below instruments that.

## The item — the KW-class locus solve must complete at the radical class, value-verified

What must become TRUE: for locus equations routed to the non-canonical branch whose residuals,
after clearing a denominator free of the radical, are polynomial in the engine's symbols and a
single quadratic radical — the measured class of the engine's own KW zero-locus residuals — the
solve completes at the operand sizes the D4 cells produce, and the emitted solution payload is the
exact solution set: sound (every branch annihilates its residual), complete (no branch silently
dropped — the shared spec forbids it), and verified, not asserted. On every cell the current
stream already completed, the repaired path's payload is byte-identical (srepr) to the stream's
own. Residuals outside that class — no radical, a radical not involving the solve variables (the
canonical branch), more than one distinct radical, a nested radical, a radical-bearing
denominator — keep their current route through an explicit guard whose behavior is demonstrated,
not assumed. Everything downstream of the solve (`:847` onward) is untouched.

⚠ The analysts measured two achievable routes, both control-validated to byte-identical payloads:
truncating the system loop at its first successfully-solved symbol while keeping the engine's own
inversion, candidate filtering, and `checksol`/`checkdens` tail; and the round-4 single-radical
lift (adjoin `t` with `t² = R`, solve the lifted polynomial pair, map back, keep branches that
annihilate the residual — the pattern already living in `algebraic_zero_test` /
`_single_quadratic_radical_determinant`). A bare single-symbol `solve` without the engine's
inversion step is measured DEAD: it diverges from the stream on a completed control. The route is
the builder's choice; the class boundary, the byte-identical-on-completed-cells obligation, and
the independent soundness-and-completeness verification are not.

## Acceptance — executable; operands discovered from the engine's own objects, no stated values

The run discipline of prior rounds binds: arm `/tmp/s11_watchdog.sh` for any engine-cell run,
stream observably, never run the full sweep, experiments on `/tmp` copies. ⚠ A full-sweep process
is running (its PID is in the session record; its stream is
`research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`): read the stream
freely, ⛔ never write to it, never signal the process.

1. Locus probe with operands discovered dynamically — the full `emit_locus` surface, not only KW:
   parse, from the live stream, every emitted `*_EQUATIONS` locus payload (every family, package,
   dimension, root — enumerate them from the stream, state no counts). Print the implemented
   guard's classification for every one. Then:
   - For every tag the guard routes to the repaired path that has a stream `*_SOLUTION`: run the
     repaired path; print the srepr equality against the stream payload — byte-identical is the
     bar — with wall time.
   - For every in-class tag with no stream solution: the solve completes within the probe budget;
     if the engine's own `algebraic_zero_test` does not prove the residual identically zero, the
     returned branch set must be non-empty; every returned branch passes BOTH independent checks
     below; and an independently-constructed exact enumeration of the solution set (a structurally
     different route built in the probe, not the engine) must agree with the returned set in BOTH
     directions — set equivalence printed, no counts stated.
     The two per-branch checks, each shown able to fail on a deliberately corrupted branch derived
     from every returned branch: (a) the engine's own `algebraic_zero_test` on the residual after
     substituting the branch; (b) exact rational specialization at several assumption-satisfying
     points — the residual vanishes at the branch and does not vanish at the corrupted branch.
   - The same-class operands the stream has not yet emitted but already determines: reconstruct,
     from the stream's emitted root values at D4 and the engine's own construction, the locus
     residuals of the roots whose `*_EQUATIONS` tags are absent, and put each through the
     no-stream verification above. A disagreement anywhere is the finding.
2. Class-boundary controls, outcomes printed: (a) constructed residuals with no radical, with two
   distinct radicals, with a nested radical, and with a radical-bearing denominator each
   demonstrably keep the current route; (b) a constructed residual whose radical does not involve
   the solve variables — nonempty canonical payload — demonstrably retains the canonical Poly
   route; (c) a constructed in-class residual demonstrably takes the repaired route after first
   printing that its canonical payload is NOT_APPLICABLE; (d) from the engine's own stream: a
   multi-equation locus system sharing one radical (the RANK_DROP/STACKED families contain them)
   is classified by the guard, its classification printed, and at least one such completed case is
   re-solved through whatever route the guard assigns it with payload byte-identical to the
   stream.
3. The `XKIN_ANISO` D2 cell completes end-to-end in-tree, wall time printed, and its stream
   compared to the same cell at the pre-fix engine: integers, statuses, and normal forms must not
   move.
4. The `XKIN_ANISO` D3 cell completes end-to-end from a fresh process — this cell contains
   completed radical KW loci, and its fresh-process completion is what a scoped rerun rests on —
   with its emitted stream byte-compared against the live stream's D3 section: any deviation is
   printed, not hidden. ⚠ This cell's wall is dominated by stages outside this item (measured
   ~13 min cold); give this control its own 1800 s budget and do not optimize anything outside
   the item to make it faster.
5. `git diff` confined to `emit_locus` and any helper it introduces; no change to tag names,
   emission order, record shapes, `PACKAGE_DIMS`, or driver order.

A green run claims exactly one thing: the KW-class locus solve wall is repaired with the solution
sets verified sound and complete. It does not claim D4 or D5 cell completion, and the decision of
whether and how the full sweep reruns is not the builder's.
