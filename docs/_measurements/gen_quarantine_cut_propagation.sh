#!/bin/bash
cd /var/projects/toy_physics || exit 1
mkdir -p docs/_measurements
OUT=docs/_measurements/quarantine_cut_propagation.md
{
cat <<'HDR'
# Measurements behind `quarantine_cut_propagation_decision_list.md`

Commands and their literal output. CLAUDE.md rule 2: a claim about an artifact carries the command that
produced it. Regenerate with the commands as written; do not transcribe.

## 1. The rule itself

```
$ grep -n "12\. \*\*A prohibition" -A 4 CLAUDE.md
HDR
grep -n "12\. \*\*A prohibition" -A 4 CLAUDE.md
cat <<'H2'
```

## 2. Where the cut is ALREADY propagated — five artifacts

```
$ grep -rn "quarantine is\|QUARANTINE APPARATUS IS CUT\|quarantine ritual is CUT\|Quarantine is cut\|quarantine never touched\|out of the tree to manufacture" CLAUDE.md .claude/skills/build/SKILL.md .claude/skills/step-run/SKILL.md docs/development_pipeline.md research/pde_ledger_v3/REBUILD_HANDOFF.md
H2
grep -rn "quarantine is\|QUARANTINE APPARATUS IS CUT\|quarantine ritual is CUT\|Quarantine is cut\|quarantine never touched\|out of the tree to manufacture" CLAUDE.md .claude/skills/build/SKILL.md .claude/skills/step-run/SKILL.md docs/development_pipeline.md research/pde_ledger_v3/REBUILD_HANDOFF.md
cat <<'H3'
```

## 3. This is the THIRD occurrence — the prior two are on the record

```
$ sed -n '857,860p' research/pde_ledger_v3/REBUILD_HANDOFF.md
H3
sed -n '857,860p' research/pde_ledger_v3/REBUILD_HANDOFF.md
cat <<'H3b'
```

```
$ sed -n '20,26p' /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_denylist_means_wrong_architecture.md
H3b
sed -n '20,26p' /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_denylist_means_wrong_architecture.md
cat <<'H4'
```

## 4. The canonical KEEP/CUT discriminator the repair points at

```
$ sed -n '137,148p' .claude/skills/build/SKILL.md
H4
sed -n '137,148p' .claude/skills/build/SKILL.md
cat <<'H5'
```

## 5. THE FOUR SURVIVING CONTRADICTIONS

### 5a. STATUS.md — orchestrator-written 2026-08-12, the proximate cause

```
$ sed -n '25,26p' STATUS.md
H5
sed -n '25,26p' STATUS.md
cat <<'H6'
```

### 5b. review-legs — the one skill that never took the 2026-08-04 rewrite

```
$ grep -n "quarantin\|out of the tree\|do.not.read\|Do not read" .claude/skills/review-legs/SKILL.md
H6
grep -n "quarantin\|out of the tree\|do.not.read\|Do not read" .claude/skills/review-legs/SKILL.md
cat <<'H7'
```

### 5c-d. The two memories that still instruct the cut mechanism

```
$ grep -n "OUT OF THE TREE\|quarantin" /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_donotread_doesnt_survive_grep.md /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_review_agents.md
H7
grep -n "OUT OF THE TREE\|quarantin" /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_donotread_doesnt_survive_grep.md /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_review_agents.md
cat <<'H8'
```

## 6. What the orchestrator did and reversed today

```
$ ls -d /tmp/s11*_leg_* /tmp/f9*_leg_* /tmp/s11_fold_leg | wc -l ; find /tmp/s11*_leg_* /tmp/f9*_leg_* /tmp/s11_fold_leg -type f | wc -l ; ls /home/trevnorris/.s11_quarantine
H8
ls -d /tmp/s11*_leg_* /tmp/f9*_leg_* /tmp/s11_fold_leg 2>/dev/null | wc -l
find /tmp/s11*_leg_* /tmp/f9*_leg_* /tmp/s11_fold_leg -type f 2>/dev/null | wc -l
ls /home/trevnorris/.s11_quarantine 2>&1
echo '```'
} > "$OUT"
wc -l "$OUT"

# appended 2026-08-12 — post-repair state, after folding both legs
{
echo
echo '## 7. POST-REPAIR — the form ablation survived the deletion (the thing most at risk)'
echo
echo '```'
echo '$ grep -n "FORM ABLATION IS MANDATORY" .claude/skills/review-legs/SKILL.md'
grep -n "FORM ABLATION IS MANDATORY" .claude/skills/review-legs/SKILL.md
echo '```'
echo
echo '## 8. POST-REPAIR — every remaining match, unfiltered'
echo
echo '⛔ THE COUNT IS NOT THE ACCEPTANCE TEST, and this section is the proof: a CUT-DECLARATION matches'
echo 'the same grep as the INSTRUCTION it replaced, so the repaired files score HIGHER than before.'
echo 'Each line below must be read and classified. Instruction => not repaired. Declaration that the'
echo 'mechanism is cut, rule 12 quoted, or an unrelated sense of the word => repaired.'
echo
echo '```'
echo '$ grep -nE "out of the tree|quarantin|do.not.read|tripwire|byte-identical" <live set> | grep -v <cut-declarations>'
grep -nE 'out of the tree|quarantin|do.not.read|tripwire|byte-identical' \
  research/pde_ledger_v3/LAUNCH_PROMPT.md research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md \
  .claude/skills/review-legs/SKILL.md .claude/skills/build/SKILL.md .claude/skills/step-run/SKILL.md \
  docs/development_pipeline.md STATUS.md research/pde_ledger_v3/steps/S11b_HANDOFF.md \

echo '```'
} >> docs/_measurements/quarantine_cut_propagation.md
