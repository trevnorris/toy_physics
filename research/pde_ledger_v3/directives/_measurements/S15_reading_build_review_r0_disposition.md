# S15 reading build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_build_review_r0_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_reading_r0_lookups.sh`). "M-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run. The legs' closed forms and conditions stay in
their scratch files; this file and the lookups cite only their residuals and existence counts.

## The round

- **Artifact.** The Codex build of the S15 reading, uncommitted on `b6e6f025`. Its contract is Part A of
  `S15_decision_list.md` at `91248b67` and Part B of `S15_reading_decision_list.md` at `7f52b4d2`. The build is:
  - `scripts/S15_separation_departure.py` and `mathematica/S15_separation_departure.wl`;
  - their transcripts, about 190 MB and 324 MB;
  - `directives/_measurements/S15_separation_departure_run_record.md`.

  All five matched the baseline after the round. The working tree was unchanged apart from this file and the lookups,
  and no file in the repository root is newer than the prompt (M0).
- **Prompt and legs.** `_scratch/s15/reading_review/s15_reading_build_review_prompt_r0.md`, sha256 `e82ea291…` (M0).
  Both legs received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s15-reading-r0-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-reading-r0-grok/REVIEW.md`, exit 0.

  The orchestrator recorded both reports' hashes when it read them, and they still match (M0).
- **Verdicts.** Both NEEDS REVISION. Each leg reached the same two findings on its own (M1).

## Findings

**V1. Item 3 is not computed in closed form for any object that depends on the separation.** (Claude F1, Grok 1.)
**CONFIRMED.**
- **The contract** (Part B item 3; M1a): "Where the integrand holds no general field besides the family's amplitudes,
  the region's quantities and the velocities, the integrals are evaluated in closed form."
- **The engines** (M2):
  - SymPy tests whether the angular integrand is a polynomial in the sines and cosines. If it is not, it returns the
    integral unevaluated (`scripts/S15_separation_departure.py:433–434`), with the status "NOT OBTAINED IN CLOSED FORM:
    angular integrand is outside this engine's exact polynomial integration route." (`:486`).
  - Wolfram makes the same kind of test (`PolynomialQ`, `mathematica/S15_separation_departure.wl:266–269`), with the
    status "NOT OBTAINED IN CLOSED FORM: non-polynomial angular factors remain outside this engine's exact integration
    route." (`:310`).
- **The transcripts** (M3):
  - Item-3 statuses:

    | Engine | Not obtained | Closed form | Left unevaluated (general field) |
    |---|---|---|---|
    | SymPy | 51 | 93 | 192 |
    | Wolfram | 48 | 96 | 192 |

  - The objects not obtained are every carried, centre-transport and radius-transport group that contains drain 1, at
    each surface and component, and the region momentum and rate of those groups. SymPy also lists the region rate of
    the drain-2 group, which Wolfram prints as `0`, `Closed form`.
  - Every closed-form object is `0` except the centre transport of the drain-2 group, which holds no `X_1` (Claude's
    tally, M3). So no printed closed form depends on `d`.
  - The item-4 branch for a closed form that depends on the separation was never reached:
    `item4/equation_before_quantifier_elimination` occurs 0 times in each transcript.
- **The closed forms are reachable** (M4). Each leg reached closed forms by its own exact route and checked them against
  a second route.
  - Claude:
    - the carried drain-1 group at the outer surface, against direct quadrature at an off-axis separation, at seven
      exponent values including the special ones: residuals at most 5e-16;
    - the region momentum of the drain-1 group, at seven exponent values: residuals at most 2.4e-12;
    - its rate, against a time finite difference: residual −4.96e-11;
    - the centre transport of the drain-1 group, and the carried group of both drains, at the outer surface, each
      component: residuals at most 1.1e-14.
  - Grok:
    - the radius transport of the drain-1 group at the outer surface, with its special-exponent branches, and the
      region momentum of the drain-1 group. None of the four printed objects contains an integral.
    - The first matches a second route at one exponent value (residual `0`). The second matches the divergence theorem
      at another (`0`) and quadrature (−3.3e-24).
    - The region momentum of the drain-2 group is `0`.
- **The labels.** Grok reads them as honest for the route the engines have. Claude reads them as describing a route
  never tried on these objects. Both readings fit the text. Each label says that the object fails a syntactic test,
  and neither engine attempted an exact integration of these objects. Part B asks for the closed form, so a label does
  not discharge item 3 here.
- **The unevaluated objects are right where Claude evaluated them** (M4):
  - SymPy's printed family objects for the carried and momentum drain-1 groups and for the rate match Claude's
    reference to within 3e-10.
  - Wolfram's carried object matches to within 3e-10. Its momentum and rate match to within 4e-7, with four
    `NIntegrate::eincr` warnings in the leg's quadrature.
- **Effect.**
  - The build prints no closed-form departure and no item-4 condition for any object that depends on `d`, so the
    reading Part B asks for cannot be made from its output.
  - Both engines decline the same objects through the same kind of test, so on these objects they cannot disagree.
- **What must be true.**
  - Part B item 3 holds for every object in the gap index: each is printed in closed form, exact on the domain, with
    its dependence on `d` explicit and every special exponent value treated.
  - An object is labelled as not obtained only where that engine attempted it exactly and the attempt did not return
    it. The label says what was attempted.
  - Each engine reaches its closed forms by its own construction (Part B, "Supplied"; M1a).
  - No item-1, item-2 or item-5 object changes.

**V2. Item 4's condition is the closed object set equal to zero, and "Q under the condition" is the family's `Q`.**
(Claude F2, Grok 2.) **CONFIRMED.**
- **The contract** (Part B item 4; M1a): "the conditions under which that object vanishes for every `D` on the domain,
  exact, with every special exponent value treated. … It is printed as equations and inequalities in the named
  quantities; `Q` on the family under each condition."
- **The engines** (M5):
  - SymPy sets `condition=S.Eq(result,0,evaluate=False)` and `qclosed=closed_route(family_q)` (`.py:504`, `:509`). It
    then prints "Exact implicit equation in the named parameters; its Piecewise branches include exceptional exponent
    values. Q is printed restricted to that equation." (`:512`).
  - Wolfram sets `parameterCondition=familyResult==0` and `conditionedQ` from `familyQ` (`.wl:323`, `:327`). It then
    prints "… Q is restricted to the printed equation." (`:330`).
  - For a closed form that depends on the separation, each engine prints "NOT OBTAINED EXACTLY" without attempting the
    condition (`.py:495–499`, `.wl:318–322`). Once V1 is repaired, every such object reaches this branch.
- **The transcripts** (M5):
  - Witness, the outer centre transport, component 1, drain-2 group. The condition is the item-3 object `== 0`, and
    `Q_under_condition` is the family's `Q`, which still holds `A_2`, `B_2` and a component of `Ẋ_2`.
  - No `exponent_branches` line in either transcript contains a Piecewise.
  - The run record says "alongside Q restricted to the condition." (`S15_separation_departure_run_record.md:47`).
- **The legs.** Both computed the witness's `Q`, and `E_Q^geo` as twice it (M5). Each also states in prose the
  conditions under which the witness vanishes. The finding rests on the engines' lines, not on those statements.
- **Effect.**
  - No printed condition is reached by computation.
  - The status lines and the run record assert a Piecewise and a restriction that no computation produced (E1).
- **What must be true.**
  - Part B item 4 holds for every closed form that item 3 prints, including those that depend on the separation. The
    conditions are reached by computation from the closed form and printed as equations and inequalities in the
    named quantities, and `Q` on the family is printed under each.
  - Where an engine does not obtain the conditions exactly after attempting them, it prints a label saying so beside
    what it did print (M1a).
  - No status line or run-record sentence states a property of a printed object that no computation in that engine
    produced.

## Examined and not filed

- **SymPy's printed form varies between runs.** Both legs reran each engine under the guard (M6).
  - Wolfram reproduced its transcript byte for byte (`651c200f…`).
  - SymPy's rerun differs in 17 tags, in how terms are arranged (Claude). Each changed tag's residual against the
    committed one is `0` in both legs' checks.
  - The run record names each transcript's sha256 and does not claim that a rerun reproduces it.
  - Nothing computed or claimed changes. Carried to the reading comparator's decision list: it compares SymPy's
    objects by value, not by bytes.
- **`departure_partition/residual` is zero by construction** (Claude). It shows only that the split adds up. Carried to
  the S15 record decision list: the record cites it as a definition, not as a check.
- **Process.**
  - Claude's first derivation (`d01_general_profiles.py`) stayed in `sympy.simplify` for about 4.7 hours. The leg
    reported it and finished without it, recomputing that reference independently (`e01`).
  - Because that job held the pool's 25 GiB reservation, the leg ran its later jobs at 12 GiB. Their peaks were at most
    5.51 GiB.
  - After both reports were filed and their hashes recorded, the orchestrator stopped the agent and the job's guard
    unit (`s11c-guard-44fc98082f0b`, now inactive; M0). Neither touched the artifact.

## What the legs checked and found sound

Both legs derived first, with every field and every change live (M4, M8).
- **Item 1.**
  - All 309 object tags of the accepted S15 step are byte-identical in both engines (Claude).
  - Each transcript reproduces the accepted transcript up to its old end marker (Grok).
  - The named sums follow Part B.
- **Item 2.**
  - Drain 1's centre moves along `D`, and every changing quantity, `Ẋ_1` included, has its own live change.
  - The held quantities are kept, and `E_Q^geo` and `E_Q^resp` are as defined.
  - With general profiles, both engines' printed `E_Q^geo` and `E_Q^resp` match Claude's brute-force reference to
    within 3e-9 (M4).
  - Each existence condition is labelled sufficient.
- **Items 5 and 6.** Item 5 is printed with `Q` unchanged. The dimension residuals are computed from signatures.
- **Domain and inputs.**
  - The domain is printed.
  - No limit or expansion goes beyond Part A's labelled pressure choice.
  - No force law, coefficient, normalisation or mass is introduced, and no object is identified with a pull.
- **Emission and independence.**
  - There are no verdict tokens and no value-gated emission.
  - The Wolfram engine reads nothing from SymPy. It is a different construction, apart from sharing the kind of
    integration test that V1 concerns.
- **Ablations** (M7). All were run on copies, and every one moved the output:
  - FORM:
    - the wrong centre moved (Wolfram, Claude);
    - `Ẋ_1`'s response dropped (Wolfram, Claude);
    - the change kept inside the geometric part (SymPy, both legs; Wolfram, Grok);
    - `u_k` frozen (SymPy, Claude).
  - One-sided corruption:
    - of the supplied balance (SymPy, Claude). The groups without `G` did not move.
    - of the family (Wolfram, both legs).

## Next

1. Preserve this reviewed baseline, with the transcripts annexed to GIN.
2. Repair V1 and V2 in the builder's session. The repair brief carries the "What must be true" lines above, and none of
   the legs' closed forms or conditions.
3. Round 1 reviews the repair (fresh Claude + Grok).
