# Measurements behind `quarantine_cut_propagation_decision_list.md`

Commands and their literal output. CLAUDE.md rule 2: a claim about an artifact carries the command that
produced it. Regenerate with the commands as written; do not transcribe.

## 1. The rule itself

```
$ grep -n "12\. \*\*A prohibition" -A 4 CLAUDE.md
60:12. **A prohibition is not a control.** Blindness is enforced by *absence* — by bounding what the builder is
61-    handed — never by a sentence forbidding a read. A do-not-read list is a denylist, and a denylist means
62-    the architecture is wrong. The measured failure is absence of computation, not anchoring; quarantine is
63-    cut and rule 2 replaced it. **This applies to these rules as well: every one is prose I drift from
64-    under load, so the ones that hold are the ones that leave an artifact whose absence you can see.**
```

## 2. Where the cut is ALREADY propagated — five artifacts

```
$ grep -rn "quarantine is\|QUARANTINE APPARATUS IS CUT\|quarantine ritual is CUT\|Quarantine is cut\|quarantine never touched\|out of the tree to manufacture" CLAUDE.md .claude/skills/build/SKILL.md .claude/skills/step-run/SKILL.md docs/development_pipeline.md research/pde_ledger_v3/REBUILD_HANDOFF.md
CLAUDE.md:62:    the architecture is wrong. The measured failure is absence of computation, not anchoring; quarantine is
.claude/skills/build/SKILL.md:141:⚠⚠ **THE QUARANTINE APPARATUS IS CUT (2026-08-04). Two mechanisms were being conflated:**
.claude/skills/build/SKILL.md:149:⇒ ⭐ Clause 1 removes the slot a typed answer goes in, which is **structural**; quarantine is
.claude/skills/step-run/SKILL.md:110:7. **Build SymPy independently of the `.wl`.** ⚠⚠ **REWRITTEN 2026-08-04 — the quarantine ritual is CUT.**
.claude/skills/step-run/SKILL.md:117:   quarantine never touched. ⇒ `research/pde_ledger_v3/REBUILD_HANDOFF.md`.
docs/development_pipeline.md:72:on a known answer — quarantine never touched it, and 2.2 kills it structurally. A do-not-read list is a
research/pde_ledger_v3/REBUILD_HANDOFF.md:857:- ⛔⛔ **DO NOT move engines, directives or transcripts out of the tree to manufacture blindness.** ⚠ An
research/pde_ledger_v3/REBUILD_HANDOFF.md:859:  `CLAUDE.md` rule 12 (*"Don't build blindness apparatus… Quarantine is cut"*), and ⛔ it **cannot** work
```

## 3. This is the THIRD occurrence — the prior two are on the record

```
$ sed -n '857,860p' research/pde_ledger_v3/REBUILD_HANDOFF.md
- ⛔⛔ **DO NOT move engines, directives or transcripts out of the tree to manufacture blindness.** ⚠ An
  earlier draft of this block said to, and **both a review leg and the user rejected it**: it contradicts
  `CLAUDE.md` rule 12 (*"Don't build blindness apparatus… Quarantine is cut"*), and ⛔ it **cannot** work
  anyway, since the PY engine is **required** to import `S10_exports.py` ⇒ [[project-export-chain-pivot]].
```

```
$ sed -n '20,26p' /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_denylist_means_wrong_architecture.md
chose to put outside the tree · stage/`git add` baselines as elaborate diff scaffolding.

⚠ **CLAUDE.md rule 12 ALREADY SAYS THIS** — *"Don't build blindness apparatus. The measured failure is
absence of computation, not anchoring. Quarantine is cut; rule 2 replaced it."* ⛔ I did it anyway across a
whole session and had to be told again. ⇒ **The control is [[feedback-scripts-print-never-assert]] and
mutation, ⛔ not custody of files.**

```

## 4. The canonical KEEP/CUT discriminator the repair points at

```
$ sed -n '137,148p' .claude/skills/build/SKILL.md

⚠ **State IN the directive when an object is supplied** and therefore **unfalsifiable within the build**, so
a passing build does not read as if it verified it.

⚠⚠ **THE QUARANTINE APPARATUS IS CUT (2026-08-04). Two mechanisms were being conflated:**

| ⭐ **KEEP — independent CONSTRUCTION** | ⛔ **CUT — hiding ANSWERS from the builder** |
|---|---|
| `.wl` written first, barred from the registry, ⛔ never a transcription of the `.py` | moving answer-bearing files out of the tree |
| two engines that can genuinely **DISAGREE** | quarantining directives, `_scratch` transcripts, the "answer-bearing set" |
| ⭐ the disagreement **is** the test | byte-identical-restore checks, tripwires, do-not-read lists for answers |

```

## 5. THE FOUR SURVIVING CONTRADICTIONS

### 5a. STATUS.md — orchestrator-written 2026-08-12, the proximate cause

```
$ sed -n '25,26p' STATUS.md
⚠ **HYGIENE, ⛔ NOT A CONTROL:** `/tmp/s11*_leg_*/`, `/tmp/f9*_leg_*/`, `/tmp/s11_fold_leg/` hold review-leg
scratch that computed real S11 physics. ⛔ Do not commit it into the tree — it is scratch, and the tree is
```

### 5b. review-legs — the one skill that never took the 2026-08-04 rewrite

```
$ grep -n "quarantin\|out of the tree\|do.not.read\|Do not read" .claude/skills/review-legs/SKILL.md
3:description: Launch two independent PDE-ledger reviews of one artifact in parallel, choosing the two legs by WHO WROTE the artifact — Codex plus Grok for orchestrator-written plans, directives and prose; a fresh Claude agent plus Grok for Codex-written scripts and TeX. Renders one identical review prompt carrying the artifact, what the leg is handed, physics checks, required ablations, and a physics-only finding filter; requires a FORM ablation on every script, since a tag can be typed prose that no computation produced. There is no do-not-read list — what a leg must not use, it is not given.
11:⛔ **There is no `--do-not-read` argument.** ⚠ It was one until 2026-08-12; a denylist means the
26:{{LIST_EVERYTHING_THE_LEG_GETS — ⛔ there is no do-not-read list; what a leg must not use, it must not be given}}
73:⛔ **There is no quarantine rule field.** ⚠ It said to hand reviewers `git show <sha>:<path>` because "the
75:is the very route quarantine has already failed through (`.claude/skills/build/SKILL.md:146`).
150:⚠ **A do-not-read list is a denylist, and a denylist means the architecture is wrong.** ⛔ If each new step
160:pre-registrations or transcripts out of the tree · the `git show <sha>:<path>` read protocol ·
```

### 5c-d. The two memories that still instruct the cut mechanism

```
$ grep -n "OUT OF THE TREE\|quarantin" /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_donotread_doesnt_survive_grep.md /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_review_agents.md
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_donotread_doesnt_survive_grep.md:33:must be quarantined too.** To make a comparison honest I wrote my predictions down **before** the scripts
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_donotread_doesnt_survive_grep.md:37:the tree** for the duration. ⛔ That is the quarantine mechanism `CLAUDE.md` rule 12 cut, and following it
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_donotread_doesnt_survive_grep.md:61:⚠ **Keep the empty output directory in place** — quarantining a whole `scripts/` dir once left the builder
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_review_agents.md:77:"keeps a parallel builder blind."** ⚠ `git show` is the very route quarantine has already failed through
```

## 6. What the orchestrator did and reversed today

```
$ ls -d /tmp/s11*_leg_* /tmp/f9*_leg_* /tmp/s11_fold_leg | wc -l ; find /tmp/s11*_leg_* /tmp/f9*_leg_* /tmp/s11_fold_leg -type f | wc -l ; ls /home/trevnorris/.s11_quarantine
20
336
ls: cannot access '/home/trevnorris/.s11_quarantine': No such file or directory
```

## 7. POST-REPAIR — the form ablation survived the deletion (the thing most at risk)

```
$ grep -n "FORM ABLATION IS MANDATORY" .claude/skills/review-legs/SKILL.md
35:⛔⛔ **FOR A SCRIPT, A FORM ABLATION IS MANDATORY, ⛔ NOT OPTIONAL — AND IT IS THE ONLY THING THAT HAS EVER
```

## 8. POST-REPAIR — every remaining match, unfiltered

⛔ THE COUNT IS NOT THE ACCEPTANCE TEST, and this section is the proof: a CUT-DECLARATION matches
the same grep as the INSTRUCTION it replaced, so the repaired files score HIGHER than before.
Each line below must be read and classified. Instruction => not repaired. Declaration that the
mechanism is cut, rule 12 quoted, or an unrelated sense of the word => repaired.

```
$ grep -nE "out of the tree|quarantin|do.not.read|tripwire|byte-identical" <the eight live files>
research/pde_ledger_v3/LAUNCH_PROMPT.md:35:   the Lagrangian/action, ⛔ do not read them off the register and agree.
research/pde_ledger_v3/LAUNCH_PROMPT.md:52:   exists.** ⇒ There is nothing to anchor to. ⛔ Do not rely on a *"do not read the `.py`"* instruction —
research/pde_ledger_v3/LAUNCH_PROMPT.md:117:7. **Have Codex build the SymPy audit + any registry insertion — ⛔ and quarantine NOTHING.** ⭐ The `.wl`
research/pde_ledger_v3/LAUNCH_PROMPT.md:121:   ⛔ **CUT: moving the `.wl` out of the tree · restore-and-verify-byte-identical · the `git show` read
research/pde_ledger_v3/LAUNCH_PROMPT.md:122:   protocol.** ⚠ All three sit in the CUT column at `.claude/skills/build/SKILL.md:143`, and quarantine
research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md:34:## 3. ⛔⛔ The tripwire — ⛔ CUT. ⭐ The measurement it produced is kept.
research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md:40:⛔ **The technique is nonetheless CUT** — `.claude/skills/build/SKILL.md:143` lists tripwires in the CUT
research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md:41:column, beside byte-identical-restore checks and do-not-read lists. ⚠ It defends against **anchoring**,
research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md:44:question a tripwire cannot ask, because a hand-typed payload passes a tripwire it never read.
research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md:100:⚠ **Worded "quarantine" until 2026-08-12.** ⛔ Same word, ⛔ unrelated mechanism — this is a **repair
research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md:102:⇒ [[feedback-quarantine-gap-governing-prose]].
.claude/skills/review-legs/SKILL.md:3:description: Launch two independent PDE-ledger reviews of one artifact in parallel, choosing the two legs by WHO WROTE the artifact — Codex plus Grok for orchestrator-written plans, directives and prose; a fresh Claude agent plus Grok for Codex-written scripts and TeX. Renders one identical review prompt carrying the artifact, what the leg is handed, physics checks, required ablations, and a physics-only finding filter; requires a FORM ablation on every script, since a tag can be typed prose that no computation produced. There is no do-not-read list — what a leg must not use, it is not given.
.claude/skills/review-legs/SKILL.md:11:⛔ **There is no `--do-not-read` argument.** ⚠ It was one until 2026-08-12; a denylist means the
.claude/skills/review-legs/SKILL.md:26:{{LIST_EVERYTHING_THE_LEG_GETS — ⛔ there is no do-not-read list; what a leg must not use, it must not be given}}
.claude/skills/review-legs/SKILL.md:41:computation behind it; the ablation produced **byte-identical output**. ⭐ Confirmed by source reading at
.claude/skills/review-legs/SKILL.md:73:⛔ **There is no quarantine rule field.** ⚠ It said to hand reviewers `git show <sha>:<path>` because "the
.claude/skills/review-legs/SKILL.md:75:is the very route quarantine has already failed through (`.claude/skills/build/SKILL.md:146`).
.claude/skills/review-legs/SKILL.md:150:⚠ **A do-not-read list is a denylist, and a denylist means the architecture is wrong.** ⛔ If each new step
.claude/skills/review-legs/SKILL.md:160:pre-registrations or transcripts out of the tree · the `git show <sha>:<path>` read protocol ·
.claude/skills/review-legs/SKILL.md:161:byte-identical-restore checks · tripwires · symmetric denylists. ⚠ Every one defends against **anchoring**,
.claude/skills/review-legs/SKILL.md:171:mechanisms. ⚠ Measured: a hand-typed payload passes every denylist, tripwire and restore check it never
.claude/skills/review-legs/SKILL.md:197:runs belongs in the rendered prompt**, which both legs receive byte-identical. A launch wrapper may carry
.claude/skills/build/SKILL.md:11:⛔ **There is no `--do-not-read` argument.** ⚠ It was one until 2026-08-12 — a denylist means the
.claude/skills/build/SKILL.md:50:   ⛔ **There is no `--do-not-read` argument to pass** — it was cut 2026-08-12 (rule 12).
.claude/skills/build/SKILL.md:145:| `.wl` written first, barred from the registry, ⛔ never a transcription of the `.py` | moving answer-bearing files out of the tree |
.claude/skills/build/SKILL.md:146:| two engines that can genuinely **DISAGREE** | quarantining directives, `_scratch` transcripts, the "answer-bearing set" |
.claude/skills/build/SKILL.md:147:| ⭐ the disagreement **is** the test | byte-identical-restore checks, tripwires, do-not-read lists for answers |
.claude/skills/build/SKILL.md:149:⇒ ⭐ Clause 1 removes the slot a typed answer goes in, which is **structural**; quarantine is
.claude/skills/build/SKILL.md:151:**anchoring**, while the measured failure was **absence of computation** — a threat quarantine never
.claude/skills/step-run/SKILL.md:82:   grows faster than the list.** ⛔ Do not grep-and-quarantine what hits.
.claude/skills/step-run/SKILL.md:110:7. **Build SymPy independently of the `.wl`.** ⚠⚠ **REWRITTEN 2026-08-04 — the quarantine ritual is CUT.**
.claude/skills/step-run/SKILL.md:114:   ⛔ **CUT: moving the `.wl` out of the tree, quarantining directives, `_scratch` transcripts and the
.claude/skills/step-run/SKILL.md:115:   answer-bearing set, and the byte-identical-restore check.** ⚠ Those defended against **anchoring**;
.claude/skills/step-run/SKILL.md:117:   quarantine never touched. ⇒ `research/pde_ledger_v3/REBUILD_HANDOFF.md`.
docs/development_pipeline.md:72:on a known answer — quarantine never touched it, and 2.2 kills it structurally. A do-not-read list is a
STATUS.md:27:the record. ⛔⛔ **It is NOT quarantine and must not be turned into any: do not move it, do not hide it from
STATUS.md:198:else**. ⇒ **Nothing remaining independently RE-DERIVES the physics outside one fresh agent.** ⛔ Do not soften this and ⛔ do not read a fix into it: it is stated here
research/pde_ledger_v3/steps/S11b_HANDOFF.md:79:- ⛔⛔ **CUT 2026-08-12 — this bullet described quarantine, sandboxing and a tripwire, ⛔ all three of which
```
