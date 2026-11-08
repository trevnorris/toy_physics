Note from the orchestrator. Your repair has been reviewed. The review found item 4 still not delivered, and the user has chosen a further repair, to be done by you in this session. The user's standing approval for repair rounds covers it.

## Where the build stands

Your engines, transcripts and run record are committed as reviewed at 90441b89. The two transcripts are git-annex links again, with read-only targets: before a run writes a transcript, delete the link itself, as the brief says. Everything else in the brief stands, including its contract, rules, running instructions and finishing steps.

Review confirmed that every condition your engines print is true, that no condition restates the question, that your run record's claims are accurate, and that items 1, 2, 3, 5 and 6, and item 4 where `d` is absent, are unchanged from the previous version. It confirmed two defects, both in item 4 for the 48 objects whose item-3 closed form depends on `d`.

**Y1. The conditions are incomplete.**
- Wolfram prints 210 numbered conditions on the 48 objects, and none of them names an exponent (`aa[·]` or `bb[·]`). Each is an amplitude, rate or velocity condition.
- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre transport and radius transport. The other 30 objects have none: the mixed drain-1/drain-2 groups at both surfaces, and the region momentum and rate of both groups. For region momentum, drain 1, component 1, it prints `A_1 = 0` and `B_1 = 0`, each with `Q` equal to 0.
- Independent review reached exponent conditions, by exact routes checked by a second route, for representative objects of these classes, including the region momentum and rate and the outer surface groups of each engine. The printed lists do not contain them. Your engines' own printed item-3 objects vanish under them.

**Y2. The labels stand beside attempts that cannot return an exponent condition.**
- Wolfram: the identity attempt is `SolveAlways[Implies[dist > …, …], dist]` (`mathematica/S15_separation_departure.wl:519–522`). All 48 returns are the unevaluated `Inactive[SolveAlways][Implies[…`. `SolveAlways` leaves an `Implies` form unevaluated even for `a x^2 + b x == 0`. The other construction, `Reduce` over `polynomialSlots` (`:530`, `:546`), excludes the exponents and both radii from its variables. On copies restricted to three objects, changing the family's exponents, scaling the item-3 object by `d/R`, or dropping every special-exponent branch leaves their item-4 conditions byte-identical.
- SymPy: when a value holds `appellf1` or `polylog`, or more than one exponent form, the construction goes to `direct_attempt(coefficients,'all_parameters')` (`scripts/S15_separation_departure.py:306–331`). That attempt treats each `d`-dependent coefficient as a formal unknown; for region momentum it returns `[{A_1(t): 0}, {B_1(t): 0}]` (transcript line 111743). 30 objects carry this entry.

## What must be true

For each of the 48 objects, in each engine:
- the conditions under which the item-3 object vanishes for every `D` on the domain, computed exactly from the printed closed form, complete, with every special exponent value treated, and printed as equations and inequalities in the named quantities: values of, or relations among, the exponents, and equations or inequalities in the amplitudes, the radii and their rates, and the velocities;
- `Q` on the family under each condition: one per condition, computed on that condition;
- each engine's attempt on an object is one that can return conditions on the exponents. A not-obtained label stands only where such an attempt did not return them, and says what was attempted;
- no printed line or run-record sentence states a property that no computation in that engine produced;
- each engine reaches its conditions by its own construction;
- no object of items 1, 2, 3, 5 or 6, and no item-4 object where `d` is absent, changes by value. Every condition printed now that is true stays printed.

## Going on

Which construction to use in each engine is your choice. A construction that does not reach the conditions is a reason to change the construction and go on, here and for the rest of this repair; it does not need a new go-ahead. Repair both engines, each by its own construction, without comparing their outputs.

Stop and report only if:
- the guard refuses admission or ends a run you did not stop; or
- the repair needs a change outside item 4: to items 1, 2, 3, 5 or 6, or to the contract.

Wait for each guarded job by polling it until it ends: check its process and its outcome.json at short intervals. Do not rely on queued notices or detached watchers, and do not use long sleeps. Keep each run's progress lines flowing as the brief says; a run that goes silent for a long stretch has failed, and you may stop your own run through its live guard as the brief describes.

Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
