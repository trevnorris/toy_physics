# Independent review — the REVISED S9 harness pilot plan

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/HARNESS_S9_PILOT_PLAN.md`

**Codex-written**, applied from two orchestrator decision lists across two rounds. ~145 lines changed
since the version that was last reviewed. It proposes rebuilding **how a derivation step is verified** —
not the physics — and proving it on step S9, whose values are already closed.

Background you should verify rather than assume: this project derives toy-physics results symbolically.
Each step is built twice, independently (one Wolfram engine, one SymPy), from one shared specification,
and a harness under `research/pde_ledger_v3/reduction/` compares their emitted objects. The measured
problem is that days per step go into making the two engines and the shared spec agree on **how to report**
results while the physics was right in both engines from the first build.

## ⭐⭐ The single most important thing to check

⚠⚠ **This document has now been revised twice to repair defects. In this repository, repairs have
repeatedly introduced NEW defects into the material just changed** — 3 of 9 shared-spec defects were
introduced by their own repair, and in this very document a **leak was introduced while repairing a
different leak.**

⇒ ⭐⭐ **Attack the REPAIRS, not just the document.** For each repaired area, ask: *did the fix actually
fix it, or did it relocate the problem, overshoot, or create a new one?*

⭐ **Review the DIFF as well as the file.** The pre-revision baseline blob is `173df9cf`; the
intermediate is `d75e3a73`. Use:
```
git show 173df9cf > /tmp/plan_v1.md
git diff 173df9cf -- research/pde_ledger_v3/HARNESS_S9_PILOT_PLAN.md
```
⛔ Do not assume a change is an improvement because it is labelled a retraction.

## What to check

### ① ⭐⭐ THE LEAK TEST — take the builder's position, and try to reconstruct the answers

The plan claims (§3c) that S9's result values, expected verdicts and named declarations have been moved
out to an orchestrator-only artifact, leaving the document safe for whoever builds the comparator.

⛔ **You are forbidden to read the adjudication artifact** (see *Do not read*). ⭐ **That is deliberate — it
puts you in exactly the builder's position.**

⭐ **From the plan alone, plus the repository the builder can see, how much of S9's answer can you
reconstruct?** Specifically:
- Can you determine what any S9 quantity equals, or what verdict any row reaches?
- Does the **structure** leak — a row count, a kind table, an example name, a §3c description that
  identifies a value precisely enough to target it?
- ⚠ Does the plan's own **§3c leak table** restate anything while describing it?
- ⭐ Does anything reachable in the repository (`steps/S9_light_requires_shear.md`,
  `reduction/checks_S9.yaml`, committed outputs) make the strip pointless anyway? ⭐ **If S9's answers are
  freely readable elsewhere in the tree, say so plainly** — that would mean the strip is theatre and the
  real control has to be somewhere else.

⭐ Run your own grep for residual values and paste the literal output. ⛔ Do not take the document's word
that it is clean.

### ② Are the four repaired areas actually repaired?

Each was a confirmed defect. For each: is the repair correct, complete, and free of new defects?

- **L3 (unit-covariance)** — was claimed to be a fabrication defence; demoted after counterexamples.
  ⭐ Is the demotion complete — does *any* residual text still credit it with provenance or fabrication
  detection, including in the risk table? ⭐ Are the stated limits now correct, or has it been
  **under**-claimed to the point of being pointless?
  ⭐ The counterexample block is claimed to be genuine output from a placeholder run. **Re-run it yourself
  and paste your literal stdout.** ⛔ If your values differ, that is a finding.
- **L6 (value comparison)** — was a regression from the current comparator; now inverted so symbolic
  equality is primary. ⭐ Is the inversion coherent with the rest of the document, or does some other
  section still assume value comparison replaces symbolic equality?
- **A1** — was defeatable by a comparator that never compares. Now per-row equivalence against the current
  symbolic comparator over baseline plus mutants. ⭐⭐ **Is the new A1 actually runnable and actually
  falsifiable?** Name a wrong implementation that still passes it, if one exists.
- **§3c (builder/oracle separation)** — ⭐ Does the split hold operationally, or does the plan still
  require the builder to know something only the adjudication set contains?

### ③ ⭐ The new material — none of it has ever been reviewed

- **The two-phase boundary** (§1): S9 is phase 1 and a smoke test; phase 2 adds S10/S11b object kinds.
  ⭐ Is the phase-2 list adequate, and is the prohibition on adopting from S9 alone stated strongly enough
  to actually bind?
- **The declared standing hole** (§3): the claim that **no layer L0–L7 catches a defect introduced
  identically into both engines from the shared spec.** ⭐ Verify that by walking the layers yourself.
  ⭐ Is the proposed mitigation (an oracle whose construction does not come from the shared spec) real, or
  is it a placeholder that no one can actually execute?
- **The fourth cost layer** (§2, semantic binding/domain). ⭐ Is it complete? Is there a fifth?
- **L8 cross-step consistency** — ⭐⭐ **the plan claims this is the ONLY oracle outside a single step's
  spec. Is that true?** Check `reduction/relations.yaml`, `reduction/quantities.yaml`, the
  `registry_residual` blocks in `checks_S9.yaml` and `checks_S10.yaml`, and
  `reduction/dimensional_homogeneity_gate.py`. ⭐ The plan states S9's R4 residual is nonzero and that
  S9 binds several quantities to a single tag — ⭐ **verify both from source**, and judge whether the
  proposed remedy (independent tags, blocking on nonzero) is right.
- **The hardened acceptance criteria** A2–A10. ⭐ For each, name a wrong implementation that passes.

### ④ Internal consistency — ⭐ two revision rounds is where contradictions breed

⭐ Read the whole document for statements that contradict each other. Known risk areas: the reordering
claim in §3 versus the standing hole; §3b's reduction gate versus §1's two-phase prohibition; R4's risk
row versus L3's demotion; the L0 kind table versus the fourth cost layer.

### ⑤ Is the plan worth executing at all?

⭐ Step back. The plan's own evidence is that the naming layer is load-bearing on a small number of rows,
and that cross-engine comparison sits below four cheaper checks and is blind to shared-spec defects.
⭐ **Given that, is this pilot the right next action, or is there a cheaper action that would settle more?**
⚠ Argue it out; ⛔ do not simply endorse the plan because it is careful.

## Do not read

⛔ **`/tmp/claude-1000/**` — this holds the adjudication artifact (S9's answers) and the previous review
legs' reports. Reading either destroys what this review is for.**
⛔ `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/` — raw build transcripts.
⛔ Any file matching `*_review*` or `*decision_list*` under `directives/` or `docs/` — prior review outputs
and the instructions that produced this revision. ⭐ You are reviewing the **result**, not the intent.
⛔ This prompt is not evidence. Claims restated here are restated so you know what to check.

## Required method — DOCUMENT branch

1. ⭐ Read the **repository first** and form your own view of what this harness standard should be:
   `reduction/checks_S9.yaml`, `reduction/checks_S10.yaml`, `reduction/engine_output_checks.py`
   (especially `check_cross_engine` and `symbolic_equal`), `reduction/harness_ablation.py`,
   `reduction/derived_or_declared.py`, `reduction/dimensional_homogeneity_gate.py`,
   `reduction/relations.yaml`, `reduction/measurements/declaration_load_ablation.py`, the S9 engines,
   `steps/S9_light_requires_shear.md`, `REBUILD_HANDOFF.md`, `CLAUDE.md`.
2. ⭐ Write down what you think the check stack and acceptance test should be **before** opening the plan.
3. ⭐ Then read the plan and the diff, and audit against ①–⑤.

⛔⛔ **A PROSE ASSERTION IS WORTH NOTHING.** Every finding needs a **quotation with a locus** — file and
line, or a literal command and its output. Where a claim is checkable by running something, ⭐ run it and
paste the literal stdout. Where you claim a check can be defeated, ⭐ **exhibit the defeating example.**

## Physics filter

> Report a finding only if it could change **what the project computes, what it may claim, or what method
> it adopts.** ⛔ No style, formatting, tone or structure findings.

⭐ Worth more than ten completeness notes: *"this check cannot fail"* · *"this defect class is caught by
nothing"* · *"this still leaks the answer"* · *"this repair introduced a new defect"* · *"these two
sections contradict each other."*

⚠ **A leg that returns "nothing survives the filter" is weak evidence on its own.** If that is genuinely
your conclusion, state what you checked, what you could not check, and what would have had to be true for
you to find something. ⛔ Do not manufacture findings to fill a quota.

## Operational constraints

- ⛔ **Do not modify the repository.** Read-only. Copy to `/tmp` if you need to run something — ⛔ but see
  the do-not-read rule: ⛔ do not read anything already under `/tmp/claude-1000/`.
- ⛔ **No Mathematica kernel is needed.** If you start one anyway, wrap it in `timeout 600` and ⛔ never
  run more than one at a time — the licence has two seats.
- ⭐ Save every script you write and its literal stdout to named absolute paths, and report those paths.
- ⭐ Rank findings most-severe first: the claim · the evidence (quotation + locus) · what must change.

## Quarantine rule

Not applicable — no parallel builder. Blindness comes from **reading order** and the do-not-read list,
and from your being denied the adjudication artifact so that ① is a real test.
