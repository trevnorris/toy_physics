# Decision list — propagate the quarantine cut into the four artifacts that still contradict it

## Why this exists

Measured 2026-08-12: the orchestrator moved 336 files out of `/tmp` to keep a builder from reading S11
answers, then reversed it. This is the **third** occurrence of the same shape, not the first:
`REBUILD_HANDOFF.md:857-860` records an earlier draft that said to move engines and transcripts out of the
tree, and that **both a review leg and the user rejected it**;
`memory/feedback_denylist_means_wrong_architecture.md:23` quotes rule 12 and records *"I did it anyway."*

⇒ The rule exists, is written down, and is propagated in five artifacts. It failed anyway, because four
other artifacts still instruct the cut mechanism. Evidence:
`docs/_measurements/quarantine_cut_propagation.md`.

## The authority — point at it, do not restate it

- `CLAUDE.md:60-64` — rule 12.
- `.claude/skills/build/SKILL.md:137-148` — the canonical KEEP/CUT table, and the discriminator this list
  turns on: **independent CONSTRUCTION is kept; hiding ANSWERS from the builder is cut.** It already names
  `git show` as a route quarantine has failed through.

## What must become true

1. No governing artifact instructs, as a blindness mechanism, that files be moved out of the tree or that a
   reader be told not to read a path. Both sit in the CUT column.
2. `.claude/skills/review-legs/SKILL.md` reaches the position `build` and `step-run` reached on 2026-08-04,
   by **pointing at** the canonical table rather than carrying its own account of it.
3. Every KEEP-column item survives unweakened — in particular that the two engines must be able to
   genuinely disagree, and that the disagreement is the test.
4. `STATUS.md`'s current-state block asserts no quarantine.
5. The two memories that still instruct the cut mechanism stop instructing it, without losing the
   measurement each was written to preserve.

## ⛔⛔ FOLDED 2026-08-12 — both legs said DO NOT APPLY AS-IS, and both were right

⚠ **The census above was not a census.** Its generator grepped a **path set I chose**, so it could only
find loci I had already guessed. ⇒ ⛔ "four contradictions against five clean authorities" was wrong in
**both** halves: the live set is **eight files**, and **three of the five I called canonical authority**
(`build`, `step-run`, `development_pipeline`) were themselves on the list. `build/SKILL.md` put
do-not-read lists in its CUT column at `:143` while requiring `--do-not-read` to invoke it at `:10`.

⭐ **`LAUNCH_PROMPT.md` is the locus that matters** — it is the paste-after-`/compact` file, so it seeds
every session. ⚠ Together with `STATUS.md:25`, that is the actual recurrence mechanism: **session-seeding
artifacts, written fast at the end of a session and never reviewed.**

⭐⭐ **The authority is now one sentence, and it already existed** —
`research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md:17`: the Wolfram engine imports
nothing and re-derives, **"that is the only blindness control in this design, and nothing else in it
should be built pretending to be one."** ⛔ Everything repaired points at that, rather than re-listing
banned mechanisms.

⚠ **The legs split on reading order and the split dissolved:** Grok called it KEEP, Codex called it CUT as
a control, and **both proposed the same replacement** — sources first, artifact second. ⇒ recorded as a
method instruction that is only structural for an **Agent** leg; ⛔ unenforceable in a one-shot
`codex exec`/`grok` run and ⛔ not to be counted as a control there.

⚠ **Grok's save:** the **FORM ablation** lives inside `review-legs/SKILL.md`, the file needing the most
deletion. ⭐ Verified intact after the edit (`:35`).
⚠ **Codex's save:** the load-bearing invariant this list failed to name — SymPy imports the prior ledger,
Wolfram imports nothing.

⛔ **Not repaired, by the "does anything EXECUTE this?" test** ([[feedback-quarantine-gap-governing-prose]]):
`S11b_py_directive.md:11`, `S11bB_SHARED_PHYSICS.md:43`, `S11bB_rev11:16`,
`docs/_review_prompt_method_prior_art.md:115` — records of completed work; S11b's directives are rewritten
for its rebuild.

⚠ **Acceptance is NOT a grep count.** ⛔ A count cannot tell an instruction from a declaration that the
instruction is cut — the repaired files' counts went *up* ⇒ [[feedback-grep-acceptance-dodgeable]].

## The one question this list does not decide (ANSWERED ABOVE — kept for the record)

`review-legs/SKILL.md:52-53` rests document-review blindness on **reading order** — *"read the source of
truth first … do not read them in the other order."* That is an instruction to an agent, not an absence.

Decide which column it belongs in: a method instruction that happens to be phrased as an order (like
*"derive it yourself before opening the artifact"*, which is KEEP), or a prohibition (CUT). If CUT, name
the bounded-input mechanism that replaces it. A document cannot be ablated, so this is not answerable by
the route that settles it for a script.

## Boundaries

- Do not weaken a KEEP-column item to resolve a contradiction.
- Do not add a mechanism, registry, checklist or new document.
- Three of the five items are deletions. Do not grow any file.
