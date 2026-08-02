# DIRECTIVE — three skills that encode the step recipe

Write three skill files under `.claude/skills/`. ⛔ Do not commit. ⛔ Do not edit anything else.

**Match the format of the existing PROJECT skills** — read `.claude/skills/` (e.g. `adversarial-audit`,
`redteam-audit`) and follow their frontmatter and structure exactly.

**The recipe you are encoding already exists and has run end-to-end twice.** Read
`research/pde_ledger_v3/LAUNCH_PROMPT.md` § *HOW A STEP ACTUALLY RUNS* — that is the source of truth for
ordering, and the FAILURE MODE box under it explains why the ordering matters. ⛔ Do not invent process.

---

## ⭐⭐ The design point that makes this worth building

The recipe is written down and was still skipped twice in one session — the review legs, both times.
⛔ **A skill that must be *remembered* fixes nothing.** So:

**`build` LAUNCHES THE REVIEW LEGS ITSELF.** Not "then run `/review-legs`" — the build skill, on seeing
the artifact appear, launches them. ⇒ Skipping review requires actively not using the tool, which is a
much harder failure to walk into than forgetting a step.

---

## The three skills

### 1 · `build` — launch a builder and its reviews as one operation

Steps it encodes:
- **leak-gate the directive before launch** — grep it for the answers it must not contain.
  ⚠ The prompt is snapshotted into argv, so editing after launch does nothing.
- launch Codex: `codex exec -c model_reasoning_effort=xhigh`, backgrounded, ⛔ **never** wrapped in a
  shell `timeout`, **absolute paths everywhere** (the shell cwd persists between calls, so a relative
  write plus an absolute read is a silent empty string). Use `--sandbox danger-full-access` when the
  build must run Mathematica.
- **verify the DELIVERABLE, not the session** — a run handed an empty prompt has exited 0 with a Stop
  hook and done nothing; the tell was ~3k tokens against ~37k for a real build. Check the artifact
  exists and the log is a plausible size.
- ⭐⭐ **then launch the review legs automatically** (see skill 2), before the caller reads the results.

### 2 · `review-legs` — two independent legs on one artifact

- one **fresh agent** and one **Grok**, in parallel. Grok invocation:
  `grok --prompt-file <file> --cwd /var/projects/toy_physics --model grok-4.5 --effort high
   --permission-mode bypassPermissions --output-format plain`, backgrounded, ⛔ no shell `timeout`.
- renders the review prompt from a template that **always** carries: the artifact, what to check, a
  do-not-read list (so the reviewer derives independently rather than checking against a remembered
  answer), and ⭐⭐ **the physics filter** — *"report a finding only if it catches a way the physics could
  be wrong; do not report 'the script would be wrong on a different input'."*
- ⭐ **demand ablations, not readings** — "ablate it and report the literal output" has repeatedly found
  real defects that code-reading missed.
- known failure shapes worth naming in every script review: a value verified using the predicate or
  definition that produced it (`c ≔ √(x)` then asserting `c² − x = 0`); a conclusion emitted as an
  unconditional literal; a check whose expected value lives inside the artifact it checks.
- ⚠ a **quarantined** artifact is still reviewable via its git blob (`git show <sha>:<path>`), which
  keeps a parallel builder blind — tell the reviewer never to restore it into the working tree.

### 3 · `step-run` — the whole step, in order

The ordering from `LAUNCH_PROMPT.md`, with these as hard sequence points:
- pre-register predictions → **commit** (timestamp priority) → **move out of the tree** (a "do not read"
  instruction does not survive a grep)
- ⭐ **review the build directives before any build runs** — the directive is the one artifact both
  engines share, so an error in it lands in both and dual-engine certifies wrong physics
- blind `.wl` **first**, before any `.py` exists
- arbiter re-run — ⚠ reproducing proves determinism, matching predictions proves it agrees with you;
  ⛔ neither is a review
- quarantine → build the SymPy side → restore → **verify byte-identical to the committed blob**
- ⛔ **commit before any repair or destructive edit** — a repair on an untracked file leaves no baseline
- all gates → orchestrator writes the step record → Codex writes the card → the card gets its own legs
- ⛔⛔ **filter before acting: a finding is not a mandate.** "Recorded, not acted on" is a complete
  disposition.

---

## Constraints

⛔ Keep them **short and operational**. These are checklists a future session executes, not essays.
⛔ Do not add gates, checks, or process that is not already in the recipe — this project has twice
lost days to building checks for its checks.
⭐ Where a rule exists because of a specific measured failure, say so in one clause. A rule without its
reason gets rationalised away the first time it is inconvenient.

## Report

Under 15 lines: the three file paths, and one line each on what a caller types and what happens.
⛔ Do not summarise this directive back to me.
