---
name: review-legs
description: Launch two independent PDE-ledger reviews of one artifact using one fresh agent and one Grok 4.5 process in parallel. Renders a blind review prompt with the artifact, physics checks, do-not-read list, required ablations, and physics-only finding filter; also reviews quarantined artifacts through committed git blobs without restoring them.
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
Derive independently. Ablate every load-bearing check and report its literal output; code-reading alone
has repeatedly missed real defects. In a script, probe: a value verified using the predicate or
definition that produced it (`c ≔ √(x)` then asserting `c² − x = 0`); a conclusion emitted as an
unconditional literal; and a check whose expected value lives inside the artifact it checks.
## Physics filter
"report a finding only if it catches a way the physics could be wrong; do not report 'the script would be wrong on a different input'."
## Quarantine rule
{{GIT_BLOB_READ_ONLY_RULE_OR_NOT_APPLICABLE}}
```

For a quarantined artifact, tell both reviewers to read it only with `git show <sha>:<path>` and to
ablate a temporary copy. They must never restore it into the working tree; the blob keeps the parallel
builder blind.

## Launch in Parallel

1. Start Grok using the **Bash tool with `run_in_background: true`**, and ⛔ **no shell `timeout`**.
   ⛔ Do NOT use a shell `&` — this harness detaches the job and notifies you on exit; `&` inside a
   foreground call leaves it untracked.

   ```bash
   grok --prompt-file /absolute/review-prompt.md --cwd /var/projects/toy_physics \
     --model grok-4.5 --effort high --permission-mode bypassPermissions \
     --output-format plain > /absolute/grok-review.txt 2>&1
   ```

2. In the **same message**, launch one fresh `general-purpose` Agent with the rendered prompt. ⛔ Use a
   fresh agent, never a fork — a fork inherits the caller's context, including the results the reviewer
   must not see. Give it no Grok output and no prior agent context.
3. Do not poll; the harness notifies you as each leg finishes. Preserve the fresh-agent and Grok reports
   **separately attributed**, and ⛔ do not turn either finding into an edit or rebuild — filter first.
