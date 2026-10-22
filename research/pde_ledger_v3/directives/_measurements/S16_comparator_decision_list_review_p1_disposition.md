# S16 comparator decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_comparator_decision_list_review_p1_lookups.md`, beside
this file (generator `_scratch/s14a/gen/s16_cmpdl_p1_lookups.sh`). "Q-n" cites a section of that file. Every computed
result quoted here is a leg's saved stdout or an engine's transcript line.

## The pass

- **Artifact.** `directives/S16_comparator_decision_list.md`, as first written, with its setup lookups
  `S16_comparator_setup_lookups.md`. Both files matched the pass baseline after the pass. The working tree was
  unchanged apart from the lookups file (Q1).
- **Prompt and legs.** `_scratch/s16/cmp/cmp_dl_review_prompt_p1.md`, sha256 `2eb8be0a…` (Q1). Both legs got the
  identical prompt, in parallel. The orchestrator wrote the list, so the legs are Codex + Grok (G1):
  - Codex (`gpt-6.1-sol`, effort xhigh, Q1) wrote `/var/projects/_review_scratch/s16-cmpdl-p1-codex/REVIEW.md`.
  - Grok (`grok-4.7`, effort high) wrote `/var/projects/_review_scratch/s16-cmpdl-p1-grok/REVIEW.md`.
  - Both exited 0 (Q1).
- **Verdicts.** Codex: two findings. Grok: CLEAR (Q2).

## Findings

**K1. Evaluation cannot see a derivative taken at another point.** (Codex 1, High.) **CONFIRMED.**
- **The list.**
  - Item 5 asks that the test functions "separate each argument's dependence: at the evaluation points, its
    derivatives in different arguments differ" (Q3).
  - Item 7's coverage controls include "forms that do not separate a function's arguments" (Q3).
  - The list never uses the words "another point" (Q3: count 0).
- **The gap.**
  - Additively separable forms meet item 5's text. Codex's three forms of degrees 1–3 have pairwise distinct partials
    at its point (Q3).
  - Under each of them, the time derivative of `F1` with its `x1` argument shifted by `1/29` minus the original is
    identically 0 (Q3). The derivative is the one in `Drate0.termRatePart.bulk.weight.localEvolution` and
    `drate1.termRatePart.bulkWeight.localEvolution`.
  - A coupled form moves it: the integrand difference is `1/29`, and its integral over the shell with radii 1 and 2 is
    `28*pi/87` (Q3).
  - So where the exact residual of a value is not resolved, and evaluation is the check, a derivative at another point
    (another radius, another centre, another field point) moves no sample. The same gap was found in the S14a part 1
    comparator's build, as H2. Its repair coupled the arguments, and the round 1 ablations then moved the samples
    (Q3). This list did not carry it.
- **What must be true.**
  - For every compared value whose exact residual is not resolved, a live function's derivative taken at another
    point, in any one argument, moves the samples in at least one form and stratum. That includes a derivative inside
    a substitution, a derivative at a composite argument, and a derivative inside a volume or surface integral.
  - A control fails when a derivative taken at another point, in any one argument, moves no sample.
  - A coverage control fails under forms in which a derivative in one argument does not depend on the others.

**K2. No stratum evaluates the unnormalized objects where `N` is zero.** (Codex 2, Medium.) **CONFIRMED.**
- **The list.** Item 5 evaluates a normalized object "where `N` is nonzero", and its strata are the generic case and
  the three motion strata (Q4). The list never names the zero-count case (Q4: no match).
- **The governing list.** The S16 decision list at `0464b80a`: "Each object normalized by `N` is stated on the domain
  where `N` is nonzero. The unnormalized balances and moments keep the full domain (D6)." (Q4)
- **The engines.** Each prints `domain.unnormalized` without a count condition and `domain.normalized` with one (Q4).
- **The gap.**
  - Codex's live witness has shell integral of `Σ` equal to 0 at `t=0`. There the unnormalized dipole component 1 is
    `124*pi/15` and `P1` is `28*pi/3`. That point is in the common unnormalized domain and outside the normalized
    one (Q4).
  - Every motion stratum of item 5 can be sampled with `N=1` (Q4).
  - So a comparator that evaluates the unnormalized objects only where `N` is nonzero meets the list. A difference
    between the engines' unnormalized objects that shows only at zero count is then never sampled. One engine
    building an unnormalized object through a normalized one would do that. The record could then claim full-domain
    agreement that the comparison did not test.
- **What must be true.**
  - The evaluation includes a stratum where `N` is zero, with the fields and the geometry live.
  - The unnormalized balances and moments are compared there. A normalized object there is printed as outside its
    domain.
  - A control: on one side, an unnormalized binding given the normalized domain. Its zero-count coverage must move.
  - A coverage control fails under evaluation of the unnormalized objects only where `N` is nonzero.

## Not findings

- **Grok's checks.** Grok checked items 1–7 against the streams and the engines and filed nothing (Q2). Facts it
  recorded, which this pass relies on:
  - The shared tag `count.surface1.rest` is SymPy's outer flux and Wolfram's inner. Joined by the shared name, its
    expanded residual is nonzero. Joined inner to inner and outer to outer, the 20 name-mapped same-surface integrand
    pairs it expanded have residual 0 (Q5). So item 4's surface rule fits these streams, and an index join does not.
  - `dipole.definition.1` is SymPy's second Cartesian component and Wolfram's first (Q5).
  - The two engines' domain lines differ in the angular and reality conjuncts, which only Wolfram prints (Q5). Item 1
    keeps content one engine prints in the residual.
- **Grok's limits.** Its symbolic check is not a residual for every object tag, and it did not produce a residual for
  the `Subs` quotient in `Nrate.termRatePart.surface0.rest.localEvolution` (Q5). These limit Grok's review. They are
  not findings on the list, which item 1 and item 7 already cover.
- **Grok on K1.** For a comparator at one point, under one form or in one stratum, Grok cited item 7's coverage
  controls, "plus forms that do not separate arguments" (Q5). K1 is a different defect: forms that meet the
  separation text still miss a derivative at another point (Q3).

## Outcome

Two findings, both confirmed. The list as reviewed is preserved, with these dispositions, before the fold. The fold
repairs item 5 and item 7 for K1 and K2. It is the one fold G2 allows. Then the comparator build.
