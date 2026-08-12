---
name: build
description: Launch one Codex build and its two independent review legs as a single operation for the PDE ledger. Leak-gates the directive before launch, backgrounds an xhigh Codex run with absolute paths and no shell timeout, verifies the deliverable rather than the session, then automatically starts fresh-agent and Grok review before the caller reads the result.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Build With Automatic Review

Invoke as `/build <absolute-directive> <absolute-deliverable> --check "<physics to check>"`.
⛔ **There is no `--do-not-read` argument.** ⚠ It was one until 2026-08-12 — a denylist means the
architecture is wrong, and this file's own CUT table at §143 had already said so.

This is one operation: build, verify the deliverable, and launch both review legs. Never return a
successful build without running `.claude/skills/review-legs/SKILL.md` yourself.

## Runbook

1. Resolve the directive, deliverable, log, review prompt, and review outputs to absolute paths before
   launch. Background shells retain their cwd between calls; a relative write followed by an absolute
   read has already produced a silent empty prompt.
2. Leak-gate the directive by grepping it with fixed-string `rg` searches for every pre-registered answer it must not
   contain. If anything matches, repair and re-run the gate **before launch**. Codex snapshots the prompt
   into argv, so editing the file after launch changes nothing.
3. Launch Codex at xhigh effort using the **Bash tool with `run_in_background: true`**. ⛔ Do NOT use a
   shell `&` — this harness detaches the job itself and notifies you when it exits; a `&` inside a
   foreground call leaves the build untracked and unreported.

   ```bash
   codex exec -c model_reasoning_effort=xhigh "$(</absolute/directive.md)" \
     > /absolute/path/OUTSIDE/the/repo/codex-build.log 2>&1
   ```

   ⭐ **Write the raw transcript OUTSIDE the repository — ⚠ as TREE HYGIENE, ⛔ not as a blindness claim.**
   ⚠ A transcript carries the engine's complete tag values verbatim and is noise in the tree.
   ⛔⛔ **It is NOT a leak to be plugged by relocation:** the same measurement (2026-08-03, several such
   files sitting in `_scratch/`, reachable by none of the naming conventions) is evidence that ⭐ **hiding
   cannot work**, ⛔ not that it should be done harder ⇒ the CUT table at §137 below, and `CLAUDE.md`
   rule 12.

   Add `--sandbox danger-full-access` when the build must run Mathematica. ⛔ Never wrap the command in a
   shell `timeout` — SIGKILL has cost 300k+ tokens.
4. Do not poll for completion; the harness re-invokes you when the job exits. On that notification,
   verify the **deliverable**, not exit status: require the
   requested artifact to exist and be non-empty, and require the transcript/token usage to be plausible
   for the requested build. Exit 0 plus `hook: Stop` has accompanied an empty prompt and no work; the
   measured tell was about 3k tokens instead of 37k+.
5. Do not open the deliverable or read its results. Immediately read and execute
   `.claude/skills/review-legs/SKILL.md` with the deliverable and `--check` from this invocation.
   ⛔ **There is no `--do-not-read` argument to pass** — it was cut 2026-08-12 (rule 12).
   Launching the legs is the build skill's responsibility; do not tell the caller
   to invoke `/review-legs` later.
   ⭐ The deliverable here is **Codex-written**, so its two legs are **a fresh Claude agent + Grok**.
   ⛔ Codex does not review what Codex wrote — see that skill's authorship table.
6. Return only after both independent review legs have completed. Keep their findings separately
   attributed so the caller can filter them before acting.

## ⭐⭐⭐ EVERY SCRIPT DIRECTIVE CARRIES THESE THREE CLAUSES — non-negotiable

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload
> must be a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose
> describing a result.
> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` **is the builder writing down the
> expected output**, and it turns an informative value into a binary crash. Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

⚠⚠ **Measured 2026-08-04: named tags at named lines in three independently-built steps are typed sentences
with no CAS object behind them** (e.g. `S11bB:348-360`, and a whole transverse block whose symbol `mu_R`
appears in **no expression in the file**). ⛔ **Do not quote a fraction of the corpus** — two review legs
rejected that as unmeasurable. ⛔ Cross-engine agreement on such a tag is **vacuous** —
both engines carry the same author's sentence — and **eight fidelity review legs did not catch it**, because
*"does it say the right thing"* and *"does it depend on anything"* are different questions.
⇒ `research/pde_ledger_v3/REBUILD_HANDOFF.md` · `[[feedback-scripts-print-never-assert]]`.

⭐ **Copy the one shape that worked** (`S11bB` lines 421–443): type a candidate, compute the object
**independently**, `emit` the symbolic difference, hard-stop if nonzero.

## ⭐⭐⭐ FOUR COROLLARIES — all measured 2026-08-04, all of them defeat the three clauses as written

**1. ⛔⛔ A HAND-TYPED CAS OBJECT IS STILL HAND-TYPED.** Clause 1 bans a *prose* payload. It does **not**
stop `emit(FullSimplify[{0, alpha q2/beta}])` — a genuine CAS object, hand-authored, with **no data
dependency on the derivation**. Delete `Det` and `Solve` entirely and the output does not move. That is
the original defect wearing algebra instead of prose, and it is just as dead.

⛔⛔ **AND NOTE THE SYMBOLS ABOVE ARE PLACEHOLDERS — THEY MUST STAY THAT WAY.** ⚠ **Measured 2026-08-04:**
this example previously carried a **real step's real root**, so an author who copied the corollary verbatim
into a directive **pasted that step's answer into it** — and the two steps shared an action, so it was the
answer to the step being built. ⇒ ⭐ **Every anti-example in every skill and every directive uses
placeholder symbols.** A warning is a sentence a builder reads as attentively as any other.
⭐⭐ **THE STRUCTURAL RULE THAT CLOSES IT — put it in every script directive verbatim:**
> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
> ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
> the chain at the ACTION, ⛔ never at a result.**

**2. ⛔ THE TAG NAME IS OUTPUT TOO.** A name may name **the object**, ⛔ never its value, ratio, sign, or
the shape of the answer. A genuine CAS payload does not rescue a name that already gave the answer away.

**3. ⛔ NO TAUTOLOGICAL RESIDUAL — clause 2 does not imply "emit a triple everywhere."** Before emitting a
difference, ask: ⭐ **were these two operands produced by INDEPENDENT ROUTES?** If `q := A/B` and you emit
`A − q·B`, it is zero **by construction** and vanishes for *any* input, including a wrong one.
⚠ **Measured:** a directive written *against* tautologies contained one; and a builder obeying "emit
operands and residual" produced `EQUATION_REFERENCE = {0,0,0}` and differenced the ansatz against its own
substitution — **operand theatre**, the appearance of a check where there is none.
⇒ ⭐ Where no second route exists, **emit the objects and say so**. An honest *"there is nothing to compare
against here"* beats a zero dressed as a check.
⭐⭐ **And where a second route is worth building, build it:** derive the same object two structurally
different ways and difference them. ⚠ Verify independence by **one-sided corruption** — break route A
only; if route B moves too, they were never independent.

**4. ⛔⛔ EMISSION MUST NEVER BE CONDITIONAL ON A PAYLOAD'S VALUE.** Whether a tag appears may depend only
on **which package and which quantity** it belongs to.
⚠ **Measured, and it deleted evidence:** a builder told to remove duplicate payloads applied it *across*
control packages, so quantities that were **correctly invariant** vanished from the output. ⛔ Nothing
caught it — names stayed unique, no untagged output, exit 0.
⭐ **Why it inverts the meaning:** a value **present and identical** = **INVARIANT**, a real result;
a value **absent** = indistinguishable from *never computed*, which is the defect being hunted.
⇒ Tag **names** must be unique; **payloads may legitimately repeat**, and that repetition is the finding.

⭐ **After every script build, run `reduction/derived_or_declared.py` on the deliverable.** A script still
showing mostly CONSTANT has not been built. ⚠ It is a triage tool, ⛔ not a verdict — see its limits in the
rebuild handoff.

## ⭐⭐ SUPPLY WHAT IS ALREADY VERIFIED — withhold only the ACCEPTANCE CRITERION

⭐⭐ **Supply every already-verified object to the builders** — setup, field content, governing equations,
supplied premises. ⛔⛔ **DO NOT BLIND THE INPUTS.** ⚠ **Under-specification has cost this ledger far more
than contamination:** `∇·u` fell out of a spec **four times**, every time it was prose instead of an
equation; and a spec describing a system with **no linear coupling at all** would have had both engines
faithfully conclude *"it does not couple."* ⇒ `[[feedback-supply-verified-blind-only-open]]`.

⛔ **Withhold exactly one thing: an acceptance criterion that references an expected value.** ⭐ The risk is
**not** a faked computation — it is the **"fix until it matches" loop.** Codex iterates to exit 0, so if
*"matches the recorded value"* is the exit condition it will get there, and a genuine disagreement — the
most valuable output available — becomes silent confirmation.
⇒ **The builder's job ENDS at compute-and-print. The diff happens on our side**, where a mismatch is a
**finding**, ⛔ not a build failure.

⚠ **State IN the directive when an object is supplied** and therefore **unfalsifiable within the build**, so
a passing build does not read as if it verified it.

⚠⚠ **THE QUARANTINE APPARATUS IS CUT (2026-08-04). Two mechanisms were being conflated:**

| ⭐ **KEEP — independent CONSTRUCTION** | ⛔ **CUT — hiding ANSWERS from the builder** |
|---|---|
| `.wl` written first, barred from the registry, ⛔ never a transcription of the `.py` | moving answer-bearing files out of the tree |
| two engines that can genuinely **DISAGREE** | quarantining directives, `_scratch` transcripts, the "answer-bearing set" |
| ⭐ the disagreement **is** the test | byte-identical-restore checks, tripwires, do-not-read lists for answers |

⇒ ⭐ Clause 1 removes the slot a typed answer goes in, which is **structural**; quarantine is
**behavioural** and has failed here repeatedly, including through `git show`. ⚠ It also defended against
**anchoring**, while the measured failure was **absence of computation** — a threat quarantine never
touched. ⛔ Builds now run **in-repo**; transcripts still go outside the repo only to keep the tree clean.

⛔ **Verify a builder's compliance declaration against the DIFF, never against its report.**

## ⛔⛔ LAUNCH DISCIPLINE — two failures, both measured 2026-08-04, together ~25 minutes

**1. ⛔ NEVER CHAIN ANYTHING IN FRONT OF `codex exec`.** A multi-line `git commit -m "…"` before the
launch, in the same Bash call, makes the `"$(<file)"` substitution come back **empty**; Codex then sits on
stdin **forever** and ⛔ **never sends a completion notification**. ⚠ Absolute paths did **not** prevent
this — the second failure had one. ⇒ ⭐ **A launch goes in its own tool call, alone.**

**2. ⭐ PROVE THE PROMPT, THEN PROVE THE PROCESS — both take seconds:**
```bash
test -s /abs/directive.md && echo "PROMPT OK ($(wc -c < /abs/directive.md) bytes)"   # BEFORE launch
b=$(wc -c < /abs/log.txt); [ "$b" -lt 500 ] && echo DEAD || echo ALIVE               # AFTER launch
```
⚠ A hung Codex looks exactly like a busy one and will happily burn an hour. ⭐ **A real build is tens of
thousands of tokens; a dead one is a few hundred bytes.**

## ⛔⛔ COMMIT THE REVIEWED ARTIFACT BEFORE A REPAIR RUNS AGAINST IT

⚠ **Measured 2026-08-04, caught with seconds to spare:** a repair build was launched against an engine
that two review legs had just ablated and that was **still uncommitted**. ⛔ An uncommitted baseline is
exactly the thing that gets destroyed. ⇒ ⭐ **The artifact a review leg reviewed must be committed before
anything overwrites it** — otherwise the review has no fixed target and the finding cannot be reproduced.

## ⚠ WRITING THE DIRECTIVE — a PROHIBITION leaks the answer as surely as an assertion

⚠ **Measured three times in one session, twice inside sentences FORBIDDING the answer**, and one of those
was introduced **by the repair for the first**:
`Print["TAG: LINEAR_IN_K"] is forbidden` · `WL_S9_..._MU_OVER_RHO is forbidden` · *"if `[μ_F]` equals
`[μ_R]` the repair has failed."*
⇒ ⭐ **Write forbidden-pattern examples with PLACEHOLDER content, never the step's real content**, and
⛔ never state the pass condition for a load-test — the builder iterates toward whatever it can see.

### ⛔⛔ THE GATE MUST PROBE SYMBOLS IN PROXIMITY, ⛔ NOT ASSEMBLED EXPRESSIONS

⚠ **Measured 2026-08-04, and the gate reported CLEAN over 52 probes while the answer sat in the file.** The
probe was the fixed string `muR/rhoBr`; the text read `muR k2/rhoBr`. ⛔ **A fixed-string probe for a
composed ratio does not match that same ratio with a factor between its symbols** — and there are
unboundedly many such factors, so ⛔ enumerating spellings is the wrong architecture.

⭐ **Probe for the SYMBOLS, co-occurring, and read every hit yourself:**
```bash
rg -n "muR" $DIRECTIVES | rg "rhoBr"        # both symbols, same line — then READ the hits
```
⭐ A handful of legitimate hits you adjudicate beats a clean report you cannot trust. ⇒ ⚠ A leak gate that
returns **zero hits** on a packet that names the step's own symbols is evidence the **probes** are wrong,
⛔ not evidence the packet is clean.

## Invariants

- Use absolute paths for every file a background job reads or writes.
- Never inspect results between deliverable verification and review-leg launch.
- Never treat process success, an arbiter re-run, or agreement with predictions as independent review.
- ⭐ Run `reduction/engine_output_checks.py` on the deliverable's output; ⛔ verify the harness's own
  claims by **ablating the harness**, never by reading its self-report.
