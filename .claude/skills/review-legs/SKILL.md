---
name: review-legs
description: Launch two independent PDE-ledger reviews of one artifact in parallel, choosing the two legs by WHO WROTE the artifact — Codex plus Grok for orchestrator-written plans, directives and prose; a fresh Claude agent plus Grok for Codex-written scripts and TeX. Renders a blind review prompt with the artifact, physics checks, do-not-read list, required ablations, and physics-only finding filter; requires a FORM ablation on every script, since a tag can be typed prose that no computation produced.
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

⛔⛔ **FOR A SCRIPT, A FORM ABLATION IS MANDATORY, ⛔ NOT OPTIONAL — AND IT IS THE ONLY THING THAT HAS EVER
CAUGHT THE WORST DEFECT.** Change the **structure** of a load-bearing object — flip a sign *and* an
off-diagonal, collapse two independent symbols into one — then re-run and report the **literal** diff.
⭐ A **COEFFICIENT** rescale tests arithmetic; only a **FORM** change tests physics, because scaling never
leaves the family.
⚠⚠ **Measured 2026-08-04:** a script `emit`ted a physics conclusion as a **typed sentence** with no
computation behind it; the ablation produced **byte-identical output**. ⭐ Confirmed by source reading at
named lines in **three independently-built steps**. ⛔ **Do not quote a FRACTION of the corpus — two review
legs rejected that as unmeasurable.** ⛔ **Eight fidelity legs missed it**, because *"does it say the right
thing"* and *"does it depend on anything"* are different questions — and a prose tag is **perfectly
faithful** to a step record that quotes it.
⇒ ⭐ **Ask of every claim: WHICH LINE COMPUTED THIS?** and give the line number or report it as uncomputed.
⇒ ⭐ Cross-check with `reduction/derived_or_declared.py`; ⚠ its limits are in
`research/pde_ledger_v3/REBUILD_HANDOFF.md` — ⛔ it is triage, not a verdict.
⚠ **An `assert` before the emit hides this** — a perturbation strong enough to flip a check kills the
process, so the leg sees only PASS-or-crash. ⭐ Report any `assert` that precedes the value it guards.

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

## ⭐⭐ REVIEWERS CARRY THE INDEPENDENT-DERIVATION BURDEN — builders do not

⭐⭐ **Give every review leg an explicit *"derive it yourself from first principles"* mandate.** The genuine
physics confirmations come from a reviewer deriving from scratch, ⛔ never from a blind builder — see
`[[feedback-supply-verified-blind-only-open]]`.

### ⛔⛔ A PROSE DERIVATION IS WORTH NOTHING — DEMAND A SCRIPT AND ITS LITERAL OUTPUT

> *"trusting grok and codex and even yourself on how they 'rederive' is not trustworthy. Unless it's in
> CAS and we can see the output from the inputs, it's not to be trusted."* — user, 2026-08-04

⇒ ⭐ A leg reporting *"I derived it independently and got X"* in prose is **a typed sentence with no
computation behind it** — ⛔ **the identical defect class this rebuild exists to remove, relocated from
the engine into the review.**

**Put this in every script-review prompt, and mean it:**
> Write your own derivation script **before** opening the artifact, and save **both the script and its
> literal stdout** to named absolute paths. ⛔ Without these, your derivation claims will be discarded.

⚠ **Two hand derivations agreeing is [[feedback-matching-number-is-not-evidence]] at full strength** —
two unverifiable claims that happen to coincide, which is how two cancelling errors survive.
⚠ **Measured value of the demand:** two legs given it produced runnable scripts and disagreed on a
load-bearing test; the disagreement was resolvable **only** because one showed the counterexample matrix
element instead of asserting a conclusion.

### ⭐⭐ ONE-SIDED CORRUPTION — the only way to test that two routes are INDEPENDENT

When an artifact claims two independent routes to the same object and emits their residual, ⛔ **a zero
residual proves nothing on its own** — it is worthless if the second route is derived from the first.

⭐ **Corrupt ONE route at a time and report which objects moved.** If breaking route A also moves route
B's object, they were never independent and the residual is decoration.
⚠ **Measured 2026-08-04:** three one-sided corruptions — a sign, a form error, and a wrong derivative
pairing — each moved **only** its own route and drove the residual nonzero. ⭐ That, and nothing weaker,
is what establishes a two-route check.

⛔ **The reviewer must NOT verify a supplied object against the artifact's own identity** — that is circular.
⭐ **Name the forbidden route explicitly in the prompt.**

⚠ **A leg that returns "nothing survives the filter" is weak evidence on its own** ⇒
`[[reference-grok-cli-review]]`.

⭐ **When the prompt asks a leg to settle a specific contested question, DEMAND COMPUTATION and say so** —
leg quality tracks that demand.

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

## ⛔⛔ SERIALIZE WHEN BOTH LEGS ABLATE MATHEMATICA — the licence has TWO seats

⚠ **Measured 2026-08-03:** two legs were launched in parallel to ablate the same `.wl`. Each builds its own
ablation harness, each spawns kernels, and the second leg **died mid-run (exit 144)** having written only
its preamble. ⛔ A dead leg is not a clean leg — the review simply did not happen.

⭐ **Rule: if the artifact is a Mathematica script and both legs will ablate it, run them ONE AT A TIME.**
Parallelism is a wall-clock optimisation; ⛔ it is not part of what makes the legs independent.
**Independence comes from freshness and from neither leg seeing the other's output** — both survive
serialization intact, so long as the second leg is given no hint of the first's findings.

⚠ Parallel launch remains correct for **document** reviews and for anything that does not spawn kernels.

## Launch in Parallel (non-Mathematica artifacts)

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
