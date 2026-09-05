# Task — make the orchestrator's rule list easier to consume and follow (produce a PROPOSAL, do not edit)

## Context
`/var/projects/toy_physics/CLAUDE.md` is the governing rule list an LLM orchestrator (Claude) must follow on **every**
task in this repository. It has 17 numbered rules grouped into "The method / The gates / The traps", plus a preamble
and a "Repository infrastructure" section. Each rule was added after a specific failure cost a work session, so the
war-story justifications are load-bearing *evidence*, not filler — but the document has grown organically and is now
hard for the orchestrator to hold in working memory and apply reliably under load.

A concrete, recent failure motivates this: rule 7 contains the phrase *"One pass, then fold and go — never iterated to
green,"* written about the orchestrator's **decision list**. The orchestrator mis-applied it to a physics **spec**
(`SHARED_PHYSICS`), when the correct standard for a physics-bearing artifact is **review-until-clear** (rules 8–10,
esp. rule 10). The ambiguity — *which artifact class gets which review discipline* — is exactly the kind of thing that
should be impossible to get wrong. There are likely more such latent ambiguities and overlaps.

## What to produce
Read `CLAUDE.md` in full (and, for context on the surrounding process, the skill files under
`.claude/skills/build/SKILL.md` and `.claude/skills/review-legs/SKILL.md`, which the rules cross-reference). Then
write a **streamlining proposal** to `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/claude_md_streamline_proposal.md`.

⛔ **Do NOT edit `CLAUDE.md` or any file under `.claude/` or the repo.** Your only write is the proposal file. The
human owner will review it and decide what, if anything, to apply.

The proposal must contain, in this order:

1. **A one-paragraph diagnosis** of why the current list is hard for an LLM to follow reliably (be specific and
   evidence-based, not generic).
2. **An ambiguity / conflict register** — a table of every place two rules can be read as conflicting, or where a rule
   is silent on which artifact/situation it governs. For each: the rules involved, the ambiguous case, and the exact
   disambiguation you propose. The rule-7-vs-8/10 case above must be one row; find the others.
3. **A redundancy / overlap map** — rules or clauses that say the same thing in different words (e.g. the several
   "don't state the expected value" clauses, the several "a prohibition is not a control" clauses). Propose which is
   canonical and what merges.
4. **A proposed structure** — a restructured skeleton that keeps **every** control but is easier to consume. Strongly
   consider: (a) a compact **decision table mapping artifact class → required review discipline** (decision list vs
   spec vs script vs step record → how many legs, who reviews, fold-once vs review-until-clear, commit gate); (b)
   separating the terse *rule* from its *war-story evidence* so the rule is scannable and the evidence is available
   but not in the way; (c) a short quick-reference the orchestrator can hold in one glance. Show the skeleton with
   headings and one example rule fully rewritten in the new form (rule 7 is a good candidate).
5. **A preservation checklist** — an explicit list confirming each of the 17 current rules' load-bearing content maps
   into your proposed structure, so nothing is silently dropped. ⛔ Correctness is king: the goal is easier
   consumption, **never** fewer controls or a weaker gate.

## Constraints
- ⛔ Do not weaken, remove, or narrow any control. If you think a control is redundant, propose a *merge* and show
  both sources map into the survivor — never a deletion of substance.
- ⭐ Ground every claim in the actual `CLAUDE.md` text with line numbers; quote the clauses you are merging or
  disambiguating. A recommendation without the quoted source is discarded.
- Keep the proposal itself tight and scannable — you are demonstrating the target style.
- Report the absolute path of the proposal file you wrote at the end.
