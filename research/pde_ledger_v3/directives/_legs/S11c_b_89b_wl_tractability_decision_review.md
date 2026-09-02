# Decision-list review — S11c-b #89b WL directive, the TRACTABILITY amendment (BEFORE the rebuild)

You are one of two independent decision legs. The #89b WL operator-unfreeze directive already passed a full
decision review (its §1–§7 construction is vetted). It has now been **amended** with a new **§3bis
(Tractability)** and a new **§5 control 8 (tractability equivalence controls)** because the naive live-operator
construction did not finish in bounded time (a kernel ran 2h+ at 100% CPU with no output). Your job: review the
AMENDMENT only — does the tractability approach stay algebraically correct (drop NO jet terms / not re-freeze),
and are its equivalence controls decisive? Do not re-litigate the already-vetted construction; flag an
interaction with it only if the amendment breaks it. Report defects; do not fix them.

## The artifact under review
`research/pde_ledger_v3/directives/S11c_b_89b_wl_operator_directive.md` — read the whole thing for context, but
your findings must target **§3bis** and **§5 control 8**.

## What the amendment claims (and must be true)
Two blow-ups were diagnosed and two fixes proposed, both claimed "algebraically equivalent, no jets dropped":
1. **The hang:** `activateSpatialDivergences` (top-down `FixedPoint[ReplaceAll[Inactive[Div][v,c] :> Sum[D[v[[i]],
   c[[i]]],…]]]`) does not terminate on nested `Inactive[Div]` — top-down activates the OUTER Div first, so `D`
   differentiates a vector still holding an inner `Inactive[Div]`, spawning fresh Divs (+3/pass, unbounded).
   **Fix 1:** activate **innermost-first** (only activate a Div whose vector is already Div-free).
2. **The slow reduction:** a global nested `Series[·,{etaBg,0,1}]`/`{sigmaW,0,1}` over the whole operator fully
   distributes grade-≥2 just to discard it (~49 s/row). **Fix 2:** apply the grade truncation **per top-level
   `Plus` summand** (exploit `Series` linearity).

The independent diagnosis evidence (scripts + literal stdout) is under `~/.s11_build/S11c_b_89b_perf_diag/`
(`run6.out`, `run7.out`, `instr6.wl`, `instr7.wl`, `full.diff`, `harness_base.wl`). You may read it, but ⛔ do
not treat it as proof — verify the claims yourself.

## The failure modes to hunt (each would silently corrupt the operator or waste the expensive rebuild)
1. **Is innermost-first activation GENERALLY equal to top-down (confluence), or only on the tested rows?** The
   whole physics rests on `D[Inactive[Div][v,c], x]` (top-down) reaching the SAME fully-activated normal form
   as activating the inner Div first then differentiating. The diagnosis showed `PossibleZeroQ = True` on
   `EW_INTERNAL` and residual-Div 0 on the bundle — but the bundle's top-down HANGS, so there is **no top-down
   reference for the full operator**. Is that a real gap? Could bottom-up produce a term (or miss one) that a
   (hypothetical, non-hanging) top-down would not — undetected because no full reference exists? **Settle
   whether `Div` activation is order-independent as a rewrite system** (is `D[Inactive[Div][v,c],x]` defined so
   that activating first vs differentiating first commute?). If you can, compute it: build a small nested-Div
   expression where the two orders *could* differ and check `PossibleZeroQ`; **save your script + literal
   stdout to a named absolute path** (a prose argument is discarded).
2. **Are the §5-control-8 equivalence checks DECISIVE for the FULL operator, or only spot-checks?** Control (a)
   compares innermost-vs-top-down only on a case where top-down TERMINATES (small); control (b) is residual-Div
   = 0 on the full operator (proves *fully activated*, NOT that it equals the intended result); control (c) is
   per-summand-vs-global Series on a bounded case. Does {(a) small-case equivalence + (b) full activation + (c)
   bounded Series equivalence} actually guarantee the full operator is correct, or is there a hole a
   term-dropping bug slips through? If there is a hole, name the stronger control that closes it (e.g. a
   full-operator reference computed a third, independent way; or a per-row equivalence across ALL rows since
   individual rows terminate under top-down).
3. **Does per-summand `Series` ever differ from the global `Series`?** `Series` is linear over `Plus`, so
   splitting should be exact — but are there summands where it is NOT safe: a summand with a pole /
   `ConditionalExpression` / `O[]`-order interaction, a shared cross-term between summands that the grade-≤1
   truncation needs, or a `thetaRule` rational denominator that behaves differently per-summand vs globally?
   Probe a summand with a rational (`1/…`) coefficient. Compute, don't argue.
4. **Does the bounded-time TARGET pressure the builder to drop terms?** §3bis says "a few minutes per case."
   Could a builder hit that by truncating retained-grade content (re-freezing) and still pass controls 8(a–c)?
   Is "tractable means cheaper never smaller" enforced by the controls, or only stated?
5. **Naming vs recipe / over-specification (rule 3).** Does §3bis supply the *mechanism* (innermost-first;
   per-summand) without dictating exact code the builder can't adapt? Or conversely, is it vague enough that a
   builder produces a still-intractable or a term-dropping version and passes? (Under- and over-specification
   both cost a multi-hour rebuild.)
6. **`activateSpatialDivergences` has other consumers** (e.g. `elRowVector`). Does the amendment's requirement
   ("reproduce top-down wherever top-down terminates") actually protect them, or could the reorder change an
   existing (currently-correct) result elsewhere?

## Method
- Read §3bis + §5 control 8 against the diagnosis evidence and the engine
  (`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`, current = pristine; the
  functions are `activateSpatialDivergences` ~L1197, `truncateScalar`/`backgroundSeries` ~L99–125,
  `variationalSource`/`linearVirtualVariation` ~L275–281, `constrainedRows` ~L1022).
- Where you dispute equivalence or control adequacy, **DEMAND COMPUTATION OF YOURSELF**: write a small WL (or
  sympy) script, run it, save script + literal stdout to a named absolute path you report. ⛔ Wrap any kernel
  in `timeout 600`, one kernel at a time, kill leftover kernels by exact pid (never `pkill -f`).
- ⛔ Do not edit the directive or engine. ⛔ Do not run the full engine.
- Report a numbered list: `severity — directive:line — problem — why it costs a round — suggested correction`.
  Severity ∈ {BLOCKER (rebuild must not launch), SHOULD-FIX, NIT}. Then `VERDICT: N findings (B blockers)` or
  `VERDICT: CLEAR`.
- A leg that finds nothing is weak evidence. If you clear it, name the two or three things you actively tried
  to break (and the computation you ran) so the orchestrator can weigh the clearance.
