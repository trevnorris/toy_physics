# Independent review — a decision list for repairing contradictory governing prose

## Artifact

`/var/projects/toy_physics/docs/quarantine_cut_propagation_decision_list.md`

Its evidence file: `/var/projects/toy_physics/docs/_measurements/quarantine_cut_propagation.md`
(regenerate with `/var/projects/toy_physics/docs/_measurements/gen_quarantine_cut_propagation.sh`).

## What this is

`CLAUDE.md` rule 12 says blindness is enforced by absence, never by a sentence forbidding a read, and that
"quarantine is cut". On 2026-08-12 the orchestrator moved 336 files out of `/tmp` to keep a builder from
reading answers, then reversed it. The list proposes repairing the four artifacts that still instruct the
cut mechanism.

Nothing is edited until this review reports. No file named in the list has been changed.

## Required method — read in THIS order

1. **First**, form your own view from the authorities, before opening the list:
   - `CLAUDE.md` rule 12 (lines 60-64) and rules 1, 2, 6, 7, 14.
   - `.claude/skills/build/SKILL.md:137-148` — the KEEP/CUT table.
   - `.claude/skills/step-run/SKILL.md:106-113`.
   - `research/pde_ledger_v3/REBUILD_HANDOFF.md:855-862`.
2. **Then** read the decision list and the measurements file.

Reading the list first anchors you to its framing, which is the thing under test.

## What to check — these four, and anything else you find

1. **IS THE CENSUS COMPLETE?** This is the item that matters most. The list claims exactly four surviving
   contradictions and five artifacts where the cut is already propagated. **Run your own search** across
   `CLAUDE.md`, `.claude/skills/`, `docs/`, `STATUS.md`, `research/pde_ledger_v3/`, and
   `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/`. Report every locus the list
   missed, and every locus it names that on inspection does not actually contradict rule 12.
   An incomplete census is the specific defect that cost this project five review rounds in the preceding
   work; asserting a count from precedent instead of measuring it is the failure being repaired.
2. **DOES REMOVING THE MECHANISM COST ANYTHING REAL?** The architecture needs two engines that can
   genuinely disagree — the SymPy engine imports the previous ledger, the Wolfram engine imports nothing
   and re-derives. If nothing is hidden, a Wolfram builder can read the SymPy engine's committed output and
   transcribe it, and the comparator then agrees trivially. State whether item 3 of the list actually
   protects against that, what does, and whether anything in the KEEP column is load-bearing for it that
   the list fails to name.
3. **THE OPEN QUESTION.** The list declines to decide whether document-review "reading order" is a KEEP
   method instruction or a CUT prohibition. Answer it, and if CUT, name the bounded-input mechanism that
   replaces it. Say plainly if you think the list was wrong to leave it open.
4. **WOULD APPLYING ANY ITEM BREAK SOMETHING?** Take each of the five "what must become true" items,
   locate the text it would delete or rewrite, and say what is lost. Flag any item whose application would
   weaken an obligation that is doing real work.

You are not obliged to agree that the cut is correct. If the evidence says rule 12 is itself wrong, or that
the two conflated mechanisms should be separated differently, say so and show why.

## Method requirements

- **Every claim about an artifact carries the command that produced it**, with its literal output. A claim
  with no command behind it will be discarded. This is `CLAUDE.md` rule 2 and it binds reviewers.
- Quote both sides for every contradiction you report: the rule text and the contradicting text, each with
  file and line.
- Report a finding only if it catches a way a build or a review could **miss physics**, or could let a
  wrong result survive. Do not report style, wording preference, or "this could be clearer".
- Do not propose new machinery, registries or checklists. Three of the five items are deletions; a repair
  that grows the corpus is a finding against itself.

## Operational constraints — these bind both legs identically

- **READ-ONLY. Change no file in the repository.** A Codex build is writing
  `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` and
  `research/pde_ledger_v3/scripts/S11_exports.py` right now. Do not touch `research/pde_ledger_v3/scripts/`.
- Do not run any Mathematica or `wolframscript` kernel; nothing here needs one.
- Write scratch to `/tmp/qcp_leg_<yourname>/` and report the absolute paths of anything you save.

## Deliverable

- A verdict per numbered item above.
- Every additional locus you found, with file:line and the command that found it.
- Anything you think the list gets wrong, including the premise.
