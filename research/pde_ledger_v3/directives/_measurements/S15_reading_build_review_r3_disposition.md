# S15 reading build, review round 3 (round 4 of the build): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_build_review_r3_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_reading_r3_lookups.sh`). "M-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run. The legs' closed forms and conditions stay in
their scratch files. This file and the lookups cite them only by residual lines and counts.

## The round

- **Artifact.** Repair 3 of the S15 reading build, uncommitted on `0e51c6de`, which changes no build file since the
  round 2 baseline `4c4e2558`. A new author (a fresh Codex session, G4) repaired both engines and the run record for
  X1 of `S15_reading_build_review_r2_disposition.md`. The build is:
  - `scripts/S15_separation_departure.py` (`472cf243…`) and `mathematica/S15_separation_departure.wl` (`952923ed…`);
  - their transcripts, about 0.43 GB (`cf4267b5…`) and 0.47 GB (`1dbe97d8…`), written by the guarded runs
    `repair3_sympy_final_02` (exit 0, 5263 s) and `repair3_wolfram_final_01` (exit 0, 653 s), each started after its
    source's last change (M0);
  - `directives/_measurements/S15_separation_departure_run_record.md` (`35edcc89…`).

  All five matched the r3 baseline after the round, and the working tree was unchanged apart from the lookups (M0).
- **Prompt and legs.** `_scratch/s15/reading_review/s15_reading_build_review_prompt_r3.md`, sha256 `f0f1bc36…` (M0).
  It is the round 2 prompt with item 0 rewritten for X1 (item 4 for the 48 objects, restatement, labels, claims,
  construction, dependence, and every other object unchanged by value against round 2) and the carried run rules: a
  leg may stop its own stalled run, sizes its memory to the job, never wraps the guard in a `timeout`, leaves no run
  behind, writes nothing in the tree, and waits by polling. Both legs received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s15-reading-r3-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-reading-r3-grok/REVIEW.md`, exit 0.

  The orchestrator recorded each report's hash before reading either, and both still match (M0).
- **Verdicts.** Both NEEDS REVISION. Each leg reached Y1 on its own; Y2 is Claude's, and Grok's report bears on it.

## Findings

**Y1. Item 4 is incomplete for the 48 separation-dependent objects, and what it omits includes the exponent relations
under which the object vanishes for every `D` while `Q` is nonzero.** (Claude F1, Grok 1.) **CONFIRMED.** X1 is not
repaired.
- **The contract** (M1). Part B item 4 asks for the conditions under which each closed-form object vanishes for every
  `D`, exact, with every special exponent value treated, "on the exponents and, where they enter, on the amplitudes,
  the radii and their rates, and the velocities", and `Q` under each. X1's "What must be true" adds "complete".
- **The engines** (M2, M3):
  - **Wolfram.** Its 210 numbered item-4 conditions all belong to the 48 objects, and 0 of the 210 name an exponent
    (`aa[·]` or `bb[·]`). Every one of the 48 carries a not-obtained label. For region momentum and outer carried
    current, drain 1, component 1, the conditions are an amplitude set to zero, each with `Q` equal to 0.
  - **SymPy.** Its 782 conditions on the 48 objects include 596 that name an exponent, in 18 objects: the inner and
    outer drain-1 groups of carried current, centre transport and radius transport. The other 30 objects (the 18
    mixed drain-1/drain-2 surface groups and the 12 region momentum and rate components) have no exponent condition.
    For region momentum, drain 1, component 1, it prints the two amplitude zeros, each with `Q` equal to 0. 39 objects
    carry a not-obtained label; the 9 outer drain-1 objects carry an empty unresolved list and no completeness claim.
- **What is omitted** (M4). Each leg reached, by exact routes of its own, exponent relations that the engines do not
  print, with `Q` nonzero and falling off as the inverse square on several of them:
  - **Claude.** Route 1 (closed forms): 18 of 19 residual lines are zero, and the nonzero one is a stratum it reports
    as not vanishing. The region rate: 11 of 15 residual lines are zero; the other 4 are near-miss strata. Route 3
    (the series in the separation): 6 of 6 structure residuals are zero, which carries its completeness argument for
    its representatives. Its closed forms match SymPy's printed `Q` at 12 of 12 points to 30 digits.
  - **Grok.** 12 exact tags are zero (the carried-current limits, radius transport, region momentum and the rate's
    amplitude parts). On SymPy's own printed item-3 form for region momentum, the departure is `0E-129` on the
    relation at two separations and `−3036.03…` off it, while `Q` matches the inverse-square prediction to `0E-130`.
  - **On Wolfram's own objects.** Claude's probes find Wolfram's printed item-3 objects vanish under the relations
    its item 4 omits: 4 exact or numeric zeros, with nonzero controls at a generic exponent.
  - **SymPy's printed exponent relations** for the outer drain-1 groups agree with both legs' derivations (Grok,
    component 1; Claude, all branches spot-checked). Claude finds SymPy's inner drain-1 sets equal to its complete
    sets, while the engine labels them unresolved. Grok did not re-derive the 30 objects with no exponent condition.
- **Effect.** The question Part B asks — for which profiles each object vanishes, and so where `Q` falls off as the
  inverse square — is answered in neither engine for the region objects and the mixed groups, and in no object in
  Wolfram. The conditions printed there are the trivial ones, with `Q` equal to 0. On these objects the comparator
  would compare two sets of amplitude zeros.

**Y2. The not-obtained labels rest on attempts that cannot return an exponent condition.** (Claude F2; Grok's
report bears on it.) **CONFIRMED.**
- **The contract** (M1). X1: "A not-obtained label appears only where an exact attempt that could have returned a
  condition did not."
- **Wolfram** (M2):
  - Its identity attempt is `SolveAlways[Implies[dist > R, …], dist]` (`.wl:519–522`). All 48 returns are `Inactive[SolveAlways][Implies[…`. Claude's probe shows `SolveAlways`
    leaves every `Implies` form unevaluated, including `a x² + b x` (M4, `p01`).
  - Its other construction is `Reduce` over `polynomialSlots` (`.wl:530`, `:546`), which `Pick`s out the exponents
    and both radii.
  - Under three ablations — the family's exponent shifted, `Q` scaled by `d/R` (FORM), and every special-exponent
    branch dropped (FORM) — only item-3 exponent inventories move; Wolfram's item-4 conditions are byte-identical (M7).
    Grok's corruption of the item-3 object by `+1` does move them (M7): the item-4 return reads its input, but its
    conditions are amplitude zeros only.
  - The status line says the native attempt "was attempted and its exact return is printed" (l.15649, M2).
- **SymPy** (M3):
  - Whenever the values hold `appellf1` or `polylog`, or more than one exponent form, the construction goes to
    `direct_attempt(coefficients, 'all_parameters')` (`.py:306–331`). 30 objects carry an `all_parameters` unresolved
    entry.
  - That attempt treats each separation-dependent coefficient as a formal unknown. For region momentum it returns
    `[{A_1(t): 0}, {B_1(t): 0}]` (l.111743), which cannot include an exponent relation.
- **Effect.** The labels are honest about what was tried, but what was tried cannot return the conditions Part B asks
  for, so the labels stand where no capable attempt was made.

**What must be true** (for the repair). For each of the 48 objects, in each engine:
- the conditions under which the item-3 object vanishes for every `D` on the domain, computed exactly from the printed
  closed form, complete, with every special exponent value treated, and printed as equations and inequalities in the
  named quantities: values of, or relations among, the exponents, and equations or inequalities in the amplitudes,
  the radii and their rates, and the velocities;
- `Q` on the family under each condition: one per condition, computed on that condition;
- each engine's attempt on an object is one that can return conditions on the exponents. A not-obtained label stands
  only where such an attempt did not return them, and says what was attempted;
- no printed line or run-record sentence states a property that no computation in that engine produced;
- each engine reaches its conditions by its own construction;
- no object of items 1, 2, 3, 5 or 6, and no item-4 object where `d` is absent, changes by value; every condition
  printed now that is true stays true.

## Examined and not filed

- **The printed conditions are true.** SymPy's identity verification is `True` throughout, and Wolfram prints 210 of
  210 `substitution_verification: True` (Claude). The `nan` residual nodes sit only on `support: False` branches.
- **No restatement.** 0 of 890 SymPy and 0 of 210 Wolfram conditions contain the object, a `Piecewise`, an Appell
  function, `d`, a direction or a logarithm (Claude). X1's restatement is gone.
- **The run record** now states that complete conditions were not obtained, and its counts and hashes match the
  transcripts (M5; Claude). Its line 5 retracts the previous record's claims.
- **Claude's path.** Claude's report cites `derive/d01c2_rate.stdout.partial`; the file is `derive/d01c2_rate.stdout`
  (M4). The lines it quotes are in that file.

## Process

- **Runs the legs stopped.** Claude stopped its own `claude-d01c` (silent about 1.7 h) with `systemctl --user stop`;
  `claude-d01b` and `claude-d01b2` ended with exit 1 on their own (M0). Grok reports stopping three of its runs; three Grok
  runs end with a guard outcome and no child outcome (`grok-derive-02`, `grok-check-closed`, `grok-check-closed-2`),
  and one Wolfram launch has no outcome and was relaunched as `grok-wolfram-form2` (M0). Each leg
  reports the runs it stopped. Both legs left no reservation; the registry is `{}` with no active unit (M0).
- **No orphan, no stray file.** After Grok's exit no Grok run remained, and the tree holds no file a leg wrote (M0).
  The carried run rules held this round.

## What the legs checked and found sound

- **Nothing else changed** (M6):
  - Items 1, 2, 3, 5 and 6 against round 2: Wolfram 20438 lines identical (both legs); SymPy 18581 identical plus 6
    item-2 pressure `whole_minus_expansion` lines that differ in bytes with `expand(r2 − r3) = 0` (both legs).
  - Item 4 where `d` is absent: SymPy 4542 identical; Wolfram 4530 identical plus 270 identical after renumbering
    generated names (Claude).
  - The code diff is confined to item 4 (Claude).
- **Reruns** (M6, Claude). Wolfram reproduces all 33736 tagged lines byte for byte. SymPy (7272 s) reproduces every
  item-4 condition, `Q` restriction, count, unresolved list and status for all 336 objects; its other differing lines
  are written-form variation (1521 intermediate item-4 operands, and 12 item-2/3 lines with expanded difference 0).
- **SymPy dependence** (M7). Under the family shift and the special-branch drop, SymPy's exponent conditions move with
  the family; under `Q·d/R` they move or vanish (83 and 52 differing condition lines; both legs).
- **Construction and independence.** Each engine's item 4 is its own: SymPy's primitive-support recursion, Wolfram's
  `SolveAlways` and `Reduce`. The `.wl` reads nothing from the SymPy side (both legs).
- **X1** states Part B item 4, no more and no less (both legs).

## Next

Item 4 for these 48 objects has now been filed in four rounds:
- round 0, V2: the object `== 0`;
- round 1, W1: solver returns discarded;
- round 2, X1: the object's jet;
- this round, Y1/Y2: true sufficient conditions, but incomplete, from attempts that cannot return exponent conditions.

Repair 3 bred no new defect: every printed condition is true and every label is accurate about what was tried. What
remains is the construction. Both legs reached the omitted conditions by exact routes, so the question is not
infeasible. Because the step's reading turns on these conditions, and this is the fourth round on one item, the
choice of how to proceed goes to the user before any repair:
1. Preserve this reviewed baseline, with the transcripts annexed to GIN.
2. The user decides: a repair that meets the "What must be true" above; accepting the build with Y1/Y2 recorded as its
   limits; or narrowing item 4 by a decision-list change.
