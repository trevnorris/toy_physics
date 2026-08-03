---
name: review-legs
description: Launch two independent PDE-ledger reviews of one artifact in parallel, choosing the two legs by WHO WROTE the artifact — Codex plus Grok for orchestrator-written plans, directives and prose; a fresh Claude agent plus Grok for Codex-written scripts and TeX. Renders a blind review prompt with the artifact, physics checks, do-not-read list, required ablations, and physics-only finding filter; enforces blindness by moving answer-bearing files out of the tree rather than by denylisting them.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Independent Review Legs

Invoke as `/review-legs <absolute-artifact-or-sha:path> --check "<physics to check>" --do-not-read <absolute-path> ...`.

## Prompt Template

Render this template to an absolute path and give the identical prompt to both legs. Do not render it
until every field, including a concrete do-not-read list, is filled.

```markdown
# Independent physics review
## Artifact
{{ABSOLUTE_PATH_OR_GIT_SHOW_SELECTOR}}
## What to check
{{PHYSICS_CLAIM_AND_ARTIFACT_ROLE}}
## Do not read
{{PREDICTIONS_SIBLING_IMPLEMENTATIONS_ANSWERS_PRIOR_REVIEWS}}
## Required method
{{SCRIPT_BRANCH_OR_DOCUMENT_BRANCH — use the one that matches the artifact}}

**If the artifact is a SCRIPT:** derive independently. Ablate every load-bearing check and report its
literal output; code-reading alone has repeatedly missed real defects. Probe for: a value verified using
the predicate or definition that produced it (`c ≔ √(x)` then asserting `c² − x = 0`); a conclusion
emitted as an unconditional literal; and a check whose expected value lives inside the artifact it checks.

**If the artifact is a DOCUMENT** (a `.tex` card, a step record, prose): read the **source of truth
first**, form your own view of what it establishes and what it does not, and **only then** read the
artifact. ⛔ Do not read them in the other order — reading the artifact first anchors you to its framing,
which is the thing under test. ⭐ Blindness for a script comes from quarantine; for a document it comes
from **reading order**, and it is just as load-bearing. Quote both sides for every finding.
⚠ Put the **build directive on the do-not-read list**: an artifact can satisfy its directive and still
misrepresent its source, and that case is exactly what this leg exists to catch.
⚠ For a `.tex` card, check `paper/macros.tex` — some fields are **suppressed in the default build**, so
reader-critical content placed in one is invisible in the PDF. This has happened.
## Physics filter
"report a finding only if it catches a way the physics could be wrong; do not report 'the script would be wrong on a different input'."
## Quarantine rule
{{GIT_BLOB_READ_ONLY_RULE_OR_NOT_APPLICABLE}}
```

For a quarantined artifact, tell both reviewers to read it only with `git show <sha>:<path>` and to
ablate a temporary copy. They must never restore it into the working tree; the blob keeps the parallel
builder blind.

## ⭐⭐ WHO REVIEWS — decided by WHO WROTE IT, ⛔ never by file type

**The rule: whatever writes it does not review it** (user, 2026-08-03). That is the whole principle; the
table is just it applied.

| the artifact was written by | the two legs are |
|---|---|
| ⭐ **the orchestrator** — plans, build directives, step records, prose | **Codex** + **Grok** |
| ⭐ **Codex** — scripts, `.tex` cards, any generated code | **a fresh Claude agent** + **Grok** |

⛔⛔ **Do NOT send an orchestrator-written directive to a fresh Claude agent as one of its two legs.** It
is the closest thing to self-review the architecture allows, and it displaces the engine that will have to
*execute* the directive — the one with the strongest reason to catch an ambiguity in it.
⚠ **Measured 2026-08-03:** the S11b directive review was run as fresh-Claude + Grok. Both legs found real
defects, so the substitution did not announce itself — ⭐ a productive review is not evidence the
composition was right.

Launch a Codex review leg read-only, at xhigh, in the background, with the same rendered prompt:

```bash
codex exec -c model_reasoning_effort=xhigh "$(</absolute/review-prompt.md)" \
  > /absolute/outside-the-repo/codex-review.txt 2>&1
```

## ⛔⛔ BLINDNESS IS ENFORCED BY ABSENCE, ⛔ NOT BY INSTRUCTION

⚠ **A do-not-read list is a denylist, and a denylist means the architecture is wrong.** If each new step
bans one more path and the next probe evades it, stop patching and move the artifact.

⭐ **Anything carrying this step's answers must be OUT OF THE TREE while a blind build or review runs** —
the pre-registration, the sibling engine's script, **the build directives**, and ⛔ **the raw build
transcripts.**

⚠⚠ **THE HOLE THIS CLOSES, measured 2026-08-03 and open for several steps:** `_scratch/` accumulates raw
Codex transcripts that contain a prior engine's **complete tag values verbatim** — `codex_s10_wl_raw.txt`
carried every `WL_S10_*` value. It is not a `.wl`, not under `mathematica/`, not named `PREREGISTERED`,
and not reachable by `git show`. ⇒ a builder could defeat quarantine **while obeying every instruction.**
⭐ **Fix: raw transcripts are written OUTSIDE the repository**, and `_scratch/<step>_*` is moved out
alongside the quarantined engine.

⛔ Keep a denylist only for large live trees that cannot be moved, and then **symmetrically in both
directives** — an entry in one and not the other silently makes one engine better-informed than the other.

## Launch in Parallel

1. Start Grok using the **Bash tool with `run_in_background: true`**, and ⛔ **no shell `timeout`**.
   ⛔ Do NOT use a shell `&` — this harness detaches the job and notifies you on exit; `&` inside a
   foreground call leaves it untracked.

   ```bash
   grok --prompt-file /absolute/review-prompt.md --cwd /var/projects/toy_physics \
     --model grok-4.5 --effort high --permission-mode bypassPermissions \
     --output-format plain > /absolute/grok-review.txt 2>&1
   ```

2. In the **same message**, launch the **second leg chosen by the authorship table above** — a Codex
   review for orchestrator-written artifacts, a fresh `general-purpose` Agent for Codex-written ones.
   ⛔ For the agent case use a **fresh** agent, never a fork — a fork inherits the caller's context,
   including the results the reviewer must not see. Give either leg no Grok output and no prior context.
3. Do not poll; the harness notifies you as each leg finishes. Preserve the fresh-agent and Grok reports
   **separately attributed**, and ⛔ do not turn either finding into an edit or rebuild — filter first.
