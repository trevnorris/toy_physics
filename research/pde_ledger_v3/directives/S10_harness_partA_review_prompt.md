# Review a large change to a cross-engine checking harness

You are one of two independent review legs. The other leg is a different model and you will not see its
output. ⛔ Do not try to find it.

## ⛔⛔ THE CHANGE IS NOT IN THE WORKING TREE — it was REVERTED

⚠ **This change was committed and then reverted, because it was landed without review. This review is
what decides whether it comes back.** ⇒ ⛔ **Reading `research/pde_ledger_v3/reduction/*` in the working
tree shows you the OLD code, not the change.** The commit is still reachable:

```bash
git -C /var/projects/toy_physics show 92461853 --stat
git -C /var/projects/toy_physics show 92461853 -- \
  research/pde_ledger_v3/reduction/engine_output_checks.py \
  research/pde_ledger_v3/reduction/checks_S10.yaml \
  research/pde_ledger_v3/reduction/test_engine_output_checks.py
```

⭐ **To RUN the changed harness, extract just the changed files.** ⛔ Do **not** copy the source tree, and
⛔ do not restore anything into the working tree. Every input the harness needs is an explicit path
argument, so the changed script reads the repo's committed data in place:

```bash
V3=/var/projects/toy_physics/research/pde_ledger_v3
for f in engine_output_checks.py checks_S10.yaml test_engine_output_checks.py; do
  git -C /var/projects/toy_physics show "92461853:research/pde_ledger_v3/reduction/$f" > "/tmp/partA_$f"
done
```

**Run the NEW version** (from `/tmp`, reading the repo's data). ⚠ `PYTHONPATH` is required — the harness
imports a sibling module — and ⭐ **I ran this exact command; it works**:

```bash
PYTHONPATH=$V3/reduction timeout 600 python3 /tmp/partA_engine_output_checks.py \
  --config /tmp/partA_checks_S9.yaml \
  --output wl=$V3/mathematica/out/S9_light_requires_shear_mathematica_audit.out \
  --output py=$V3/scripts/out/S9_light_requires_shear_sympy_audit.out \
  --registry-dir $V3/reduction
```

⚠ **The same command with the `S10` config and the two `S10` `.out` files exits 2 with**
`OPERATIONAL_FAILURE: cannot parse tagged output ...: rejected invalid tag line: 'Solve::svars: ...'`
⭐ **I am not telling you what to conclude from that** — reproduce it, work out where it comes from, and
say whether it is in scope for this change. `S9` is a working configuration you can use meanwhile.

**Run the OLD version** for comparison — straight from the repo, no copy, nothing modified:

```bash
timeout 600 python3 $V3/reduction/engine_output_checks.py \
  --config $V3/reduction/checks_S10.yaml \
  --output wl=$V3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out \
  --output py=$V3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out
```

⭐ Swap `S10` for `S9` in the config and both `.out` names to exercise the S9 configuration.
⭐ To ablate an **input**, extract a copy of an `.out` file to `/tmp`, corrupt the copy, and point
`--output` at it. ⛔ Never edit the committed `.out` files themselves.

**The directive it was built from:** `research/pde_ledger_v3/directives/S10_harness_repair_directive.md`

## Background

The harness compares the tagged stdout of two independently-written symbolic-algebra engines — a
Mathematica `.wl` and a SymPy `.py` — whose committed output sits in
`research/pde_ledger_v3/mathematica/out/` and `research/pde_ledger_v3/scripts/out/`.

**The defect this change fixes:** the layer that partitions main-package tags into RESPONSIVE and
INVARIANT against control packages — the harness's **ablation** test, the one thing distinguishing a
computed output from a typed one — matched **zero** tags, because it keyed on a naming convention from an
earlier step. It reported `compared=0 responsive=0 invariant=0`, ⛔ **which reads as health.**

## Operational constraints — ⛔ these bind you

- ⛔ **READ-ONLY on the repository.** Do not edit, create, or delete any file under `research/`. Do not
  `git add`, `git commit`, `git stash`, `git checkout`, or `git revert`. Scratch work lives in `/tmp/`.
- ⛔ Wrap **every** command that runs a script in `timeout 600`. A 600 s hit is a failed experiment —
  report it and move on. ⛔ **Never raise the limit.**
- ⛔ **Do not run Mathematica, `math -script`, or `wolframscript`.** ⚠ The licence has two seats and an
  orphaned kernel has previously exhausted this machine's memory. Nothing in this review needs one — the
  engine output you need is already committed as text.
- ⭐ Save every script you write **and its literal stdout** to named absolute paths under `/tmp/`, and
  report those paths. ⛔ A prose claim that you "verified" something, with nothing shown, is discarded.

## ⭐ The whole repository is open to you — with ONE withhold

⭐ Read the engines, the outputs, the configs, the tests, the step records, anything that helps. ⛔ There is
no clean-room here and building one has cost this project more than contamination ever has.

⛔ **The one withhold: `REBUILD_HANDOFF.md` and the subject line of `92461853` both state a VERDICT on this
exact change** — that the ablation layer is now alive. ⭐ **That assertion is what you are testing.** Treat
it as the claim under examination, ⛔ never as a finding you may rely on. If you have already seen it, that
is fine — just do not let it stand in for your own measurement.

## The questions, in priority order

**Q1 — ⭐ Find any remaining check in this file that can report success without having compared
anything.** This is the defect class the change exists to fix. A previous review found several instances.
⭐ **Assume the fix did not catch them all, and assume the fix introduced new ones.** For every layer the
harness reports, answer: *what does its output look like when its input is empty, mis-keyed, or filtered
to nothing, and is that distinguishable from a clean pass?* Name each with line numbers. ⛔ Do not stop at
the first.

⭐ **Establish this by ABLATION, not by reading.** Corrupt an input — empty a package, rename a tag family,
truncate a payload — re-run the harness in `/tmp/partA`, and report the **literal** before/after counters.
⚠ A layer that reports the same numbers after you have destroyed its input is the finding.

**Q2 — Is the new coverage guard honest, or can it be satisfied vacuously?** The change reports
`CONTROL_COVERAGE: main=836 compared=314 uncovered=522 uncovered_fraction=0.62 no_ablation_D=[D2,D5]`.
⭐ **Verify those numbers yourself** from the output files. Then ask: can a layer pass this guard while
still being inert over a region it claims to cover? ⭐ Construct the case if so.

**Q3 — The dimension counter was split into `compared` and `vacuous`.** ⭐ Is the split drawn in the right
place? Does `compared` count only tags where a real dimensional comparison happened? Is anything counted
as compared that was not? ⭐ Verify against actual payloads, not against the code's intent.

**Q4 — Regressions, and ⛔ verify the attributions rather than accepting them.** The S9 configuration must
still work; its recorded baseline is `CONTROL_RESPONSE: compared=170 responsive=150 invariant=20
unparsed=1 unpaired=8`. The unit suite was `2 failed, 15 passed` before this change and is now
`5 failed, 12 passed`. ⭐ **Establish the cause of each new failure yourself**, and say for each whether it
is a real regression, a test that was already wrong, or a test the change should have updated.

**Q5 — Did the change weaken anything?** Guards removed, error paths softened to warnings, a report line
that now says less than it did, a test whose assertion was loosened.

**Q6 — Is the change's own key claim true?** It asserts the ablation layer is now alive. ⭐ Decide that
from behaviour you produced, ⛔ never from the harness's self-report — a checker that mis-parses two
expressions into agreement manufactures confidence, and its counters look identical either way.

**Q7 — Anything in the change you believe is wrong.** ⭐ A leg that finds nothing is weak evidence.

## ⚠ Scope

⭐ This change is **part A of two**. A part B is planned for the Mathematica reader and the dimension
table. ⛔ Do not report part B's work as a defect in part A — but ⭐ **do** say if you think the split was
drawn in the wrong place, or if part A cannot be judged without part B.

## Output format

```
VERDICT: <one line — should this change be re-landed as written?>
BLOCKING: <numbered; each with file:line and why it blocks>
Q1 SILENT-ZERO LAYERS: <table: layer, line, ablation you ran, before/after counters, distinguishable?>
Q2 / Q3 / Q4 / Q5 / Q6 / Q7: <one short block each>
SOLID: <what you checked and found correct — with what you actually ran>
PATHS: <absolute paths of every script and stdout you saved>
```

⛔ Do not report the physics conclusions of the step; report on the instrument.
⛔ "Looks fine" without saying what you ran is not a review.
