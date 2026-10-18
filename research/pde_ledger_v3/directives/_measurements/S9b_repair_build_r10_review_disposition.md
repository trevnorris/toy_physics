# S9b repair build, review round 10 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_repair_build_r10_review_lookups.md` beside this file
(generator `_scratch/s9b_build/gen/s9b_build_r10_lookups.sh`).

## The round
- **Artifacts.** Both ablation harnesses, the Wolfram checks script, the SymPy tag inventory and run measurements,
  the Wolfram round-10 measurement files and its rewritten round-8 measurements file, and both builders' reports.
  34 files are in the baseline `_scratch/s9b_build/s9b_repair_build_review_baseline_r10.sha256`, all OK at review.
- **Contract.** The build directive at `36ccec3e` (amendments 1–12), item 14 as amended by amendments 11 and 12, and
  the round-9 disposition's "What must be true after the repair". The spec is at `34a994d5`.
- **Builders.** The same `gpt-6.1-sol` session per engine as in rounds 8 and 9, resumed. **Item-15 deviation, recorded
  again:** item 15 names `gpt-6-astra` at high effort, in a fresh session per engine.
- **Run history.**
  - SymPy: run `full1` exited 1 in K10 (`ValueError: 2*sqrt(2) is not symbol-like`, a substitution into a
    `ConditionSet` bound symbol). The builder fixed it, passed validation run `validation4`, and ran `full2` (exit 0,
    33,340 s), checks run `checks5` (exit 0) and `evidence1` (its `finalize.py` over those two runs, exit 0). Only
    `full2`, `checks5`, `evidence1` and its scratch directory were handed to the legs.
  - Wolfram: final runs `s9b-r10-full`, `checks-final`, `engine-repo` and `engine-alone`, at 8 GiB, all exit 0.
- **Engines.** Both engines and `S9b_exports.py` are byte-identical to `86097a44`.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/s9b_repair_build_review_prompt_r10.md` (sha256 `5d43f2d8…`).
  The legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high). Their reports are
  `_scratch/s9b_build/s9b_repair_build_review_r10_claude.md` and `…_r10_grok.txt`, with evidence directories beside
  them. Grok returned **CLEAR**; Claude returned **NEEDS REVISION**. The legs ran after the machine's memory went to
  64 GB; the guard refused some launches while both legs ran, and both waited and retried.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **Each engine's K10 copy negates a different object, and the Wolfram witness is the negation itself.** The Wolfram copy negates the whole `propagationGate = rayConstruction && Inactive[ForAll][…]` (harness line 29; engine lines 262–264). That object is what the engine emits as `BRANCH_EXISTENCE` (engine line 275), so the corrupted tag is literally `Not` of the baseline tag, and the printed witness is that tag. SymPy negates `real_domain`, the along-ray condition (engine line 696), in the gate copy at line 770, after `BRANCH_EXISTENCE` is emitted at line 720; the ray root sets stay conjoined in `travel_domain` (line 771). Grok recorded the same facts and read them as each harness following its engine. | **Confirmed.** The round-9 disposition required the Wolfram K10 mutation to be "the negation of the engine's computed branch-existence condition, the same object SymPy negates, before the observables are gated on it". The spec defines that condition as the branch being real and propagating everywhere along the ray (spec line 193), with ray traversal a separate condition. Wolfram's copy also negates `rayConstruction`, and it changes the emitted `BRANCH_EXISTENCE` tag, which SymPy's copy leaves alone. So the two engines' K10 runs still do not ablate the same check. The Wolfram witness also cannot fail whenever the gate is decided, so it shows only that the mutation was applied, not that the observables depend on the gate. This is the material round 10 changed, and it is a new defect in it. |
| 2 | Claude 2 | **The SymPy K10 attempts cannot decide the gated domain parts, so "undecided" comes from the method, not from the engine's objects.** `evaluate_attempt` reduces a `ConditionSet` only when its binder is a single symbol (harness line 343), and the engine's ray root sets have tuple binders. Every attempt uses one fixed `PROFILES`/`POINT` set (lines 232–238), with no `j_n` profile and no `Phi`. The leg's counts: 1,259 of the 3,022 undecided attempts still hold tuple-binder ray sets after evaluation; `j_n` is unfixed in 1,046 attempts and `Phi` unassigned in 714. At a flat admissible profile the leg decided one gated domain part itself: the baseline copy true and the corrupted copy false, which is a witness that the harness's attempt cannot certify. Grok checked the "undecided" label against the attempt counts, not whether an attempt could succeed. | **Confirmed.** The round-9 disposition required "each K10 site part gets an attempt that can return a witness … with every attached condition decided", and item 14 requires every live profile, including `j_n` where present, and every independent input to be assigned. Round 10 fixed single-binder `ConditionSet` conditions; the tuple-binder ray sets and the missing inputs remain. This is the round-9 finding, narrowed, in the same material. |

**Sound, both legs.** Every printed non-K10 witness in both engines (K1–K9b and K11, SymPy and Wolfram) was
recomputed by each leg with its own scripts and matches. The SymPy K10 no-change decisions (1,141) compare literally
identical operands. Each leg's form ablation of the identity mark changed the printed marks in both harnesses. The
SymPy inventory matches its run. Neither harness prints a verdict or an expected value.

**Not clear.** Both findings change what the harnesses print for K10 and what may be claimed for it in both engines.
K10 has now failed two repair rounds in both engines, in the material each round changed. Under G4 the next repair
of K10, if there is one, changes the author. How to finish is the user's decision (below).

## The user's decision (2026-10-18)

The user chose through the question tool; the option texts are mine. Chosen: "Repair K10, new builders
(Recommended)", whose description reads: "Fresh Codex builders on the directive's named model (gpt-6-astra), one per
engine, K10 only. Wolfram negates the same along-ray condition SymPy does, in the gate the observables use, leaving
the emitted tags alone. SymPy's search assigns every input and solves the ray sets, so it can find the witness that
exists. Then SymPy's ~9.3 h full rerun, a 5-minute Wolfram rerun, and a fresh two-leg review: about 12–14 h of work.
It can't start until the running Wolfram #89b job frees the kernel and memory (likely several hours). Another round
is possible." The other options were to record K10 as not shown in both engines, or to fix Wolfram's K10 only.

So repair round 11 repairs K10 only, in both harnesses, with a new author per engine: a fresh Codex session on item
15's builder model (`gpt-6-astra`, high effort), each in a fresh export of the round-10 preservation commit that
holds only its own engine's S9b files. This ends the item-15 model deviation recorded since round 8.

## What must be true after the repair
- **Wolfram.** The K10 copy replaces the condition that the branch is real and propagating everywhere along the ray
  (spec line 193; in this engine, the along-ray `ForAll` inside `propagationGate`) by its logical negation in the
  gate the observables are gated on. The ray construction, the traversal condition, and the emitted
  `BRANCH_EXISTENCE`, `PATH_TRAVERSAL` and `BRANCH_TYPE` objects are unchanged in the copy. K10's site parts and
  outcome are derived from that copy.
- **SymPy.** Each K10 attempt assigns every live profile, including `j_n` where present, and every independent input
  free in either copy or its conditions, including `Phi`; an admissible profile may be chosen per part. The engine's
  own ray root sets are decided by the engine's CAS at those inputs, keeping every solution and choosing no ray
  independently, so that an attempt can return a witness at a gated domain part. Unresolved evaluations still
  furnish no witness. K10's outcome is derived from those attempts; where other knives use the changed code, they
  meet the same.
- Nothing else changes. Both harnesses rerun from final code, and the inventory, measurements and reports are
  regenerated.

The round-10 harnesses and files, as reviewed, are preserved before the repair overwrites them.
