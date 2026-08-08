# Independent review — W3/W4 round 5 · THE COMMIT GATE

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — Codex-written, sixth round

`/var/projects/toy_physics/research/pde_ledger_v3/reduction/`. Diff: `git diff HEAD -- research/pde_ledger_v3/reduction/`
**and** `git status` — some staged, some not, some untracked. Report: `reduction/W3_REGISTRY_D_LAWS_REPORT.md`.

⭐ **Small and targeted.** Rounds 1–4 were reviewed by eight legs; ⛔ their substance is settled and is not
what you are reviewing. ⭐ **Nothing else stands between this and a commit.**

## ⭐ What round 5 was told to fix

1. The pin's central claim — *"11 pairs, 0 circular"* — was **asserted, not computed**. Two one-line edits
   to its population table made it false while the pin reported `PASS` and 138 tests passed: repointing the
   `S9-wl` record at the **SymPy** transcript (one engine compared with itself), and deleting both `wl`
   records (Python-only).
2. `README.md` claimed every command here is bounded by `timeout 600`, while this change set added one
   needing ~33 minutes.
3. A two-spelling denylist over stdout reappeared in `test_dimension_laws.py` after being correctly removed
   from the ablation script — **third location for that pattern** in this file set.

## ⭐⭐ What to attack

**① Do the two weakenings fail, and does a third?** ⭐ Build both yourself and confirm. ⭐ Then build a
**third** the first two do not cover — ⛔ if you cannot, say what class you searched and why it is empty.
⚠ Harnesses for the first two exist at
`/tmp/claude-1000/-var-projects-toy-physics/36f37d88-e717-46ce-a790-6f9d1ef3d7bc/scratchpad/w3r4/`
(`mutate_circular_operand.py`, `mutate_drop_wl_engines.py`), each taking the pin file's path. ⭐ Use them or
better ones; ⛔ do not trust the artifact's own battery alone.

**② ⭐⭐ WHERE DOES THE REGRESS BOTTOM OUT — and is it honest?**
⭐ Round 4 pinned the guard clauses; nothing pinned the population table. ⭐ Round 5 pins the population
table. ⇒ ⛔ **what pins the thing that pins the population table?** ⭐ Find the level at which this stops,
and judge whether stopping there is **honest** — some property must ultimately be asserted, ⭐ and the
question is whether the artifact *says* which one, or lets a reader believe the regress closed.
⛔ If the new guard can itself be silently weakened, we are one level up with a longer script and no more
evidence.

**③ Is the new guard a measurement or a spelling check?** ⛔ A denylist over tag names or engine labels is
what round 4 was told to remove and round 5 was told not to reintroduce. ⭐ Determine what the guard
actually computes. ⭐ Could it pass on a population that *is* circular but spelled differently — or fail on
one that is legitimately independent?

**④ Does the printer whitelist bind?** ⭐ Confirm it passes on the shipped printer and fails on added prose.
⭐ Then check the harder case: prose that happens to **start with** a computed token.

**⑤ Are the README claims true?** ⚠⚠ **This file's coverage paragraph has now been wrong twice, in opposite
directions, and its timeout sentence once.** ⭐ Check every claim against measurement.

**⑥ Did anything regress?** ⭐ Re-run every engine and diff against committed outputs.
⚠ The orchestrator has independently confirmed all five byte-identical, including `S11-py` at
`sha 5ed934e5…` with empty stderr and `S9-wl`/`S10-wl` via `wolframscript`.
⚠ **The build reports `wolframscript` exiting 255 and using a different invocation for the WL engines** —
⭐ the orchestrator measured `wolframscript` working (exit 0) and believes that was transient licence-seat
contention. ⭐ **Check this yourself**: ⛔ a build substituting its own tool and headlining PASS is a known
failure mode here.
⚠ `WL_S10_RUNTIME_SECONDS` moves every run. ⚠ `wolframscript` writes 440 bytes of config warnings to
**stderr** — ⭐ separate the streams. ⛔ One Mathematica kernel at a time (two seats); ⛔ `timeout 900` each.

**⑦ Rule 2 across the round's diff.** ⛔ Any status typed rather than computed? ⛔ Any `assert` before the
value it guards? ⛔ Any expected value typed beside the value it checks? ⛔ Any denylist?

## Method

⭐ Read the sources of truth first, form your own view, ⭐ **then** the diff, then the report.
⛔ **Do not modify the working tree** — it holds substantial uncommitted work.
⛔⛔ **Never `git checkout`, `git stash`, or `git restore`.**
⛔⛔ **ABSOLUTE PATHS ONLY. ⛔ Never rely on the shell's working directory.** ⚠ That has failed **four** times
in this session, including once for the orchestrator, where it edited live files. ⭐ For a mirror, the pin
resolves paths relative to its own location and the registry validates loci into `pde_ledger_v2` and
`steps/` — ⭐ symlink those siblings or the registry will refuse to load.
⭐ Save every script and its literal stdout to named absolute paths.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"a weakening still passes"* · *"the regress does not bottom out honestly"* · *"the guard is a
spelling check"* · *"a README claim is false"* · *"the WL engines were not verified by the tool claimed."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.
⭐ **This is a commit gate — a clean result, stated plainly, is an actionable one.**

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
