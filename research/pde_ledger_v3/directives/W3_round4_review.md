# Independent review — W3/W4 round 4 · THE COMMIT GATE

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — Codex-written, fifth round on this material

`/var/projects/toy_physics/research/pde_ledger_v3/reduction/`. Diff: `git diff HEAD -- research/pde_ledger_v3/reduction/`
**and** `git status` — some staged, some not, some untracked. Report: `reduction/W3_REGISTRY_D_LAWS_REPORT.md`.

⭐ **This round is small and targeted.** Rounds 1–3 were reviewed by six legs; the substance is settled and
⛔ is not what you are reviewing. ⭐ **Nothing else stands between this and a commit**, so the question is
narrow: **did round 4 fix its three defects without breaking anything, and is every claim it now makes
true?**

## ⭐ What round 4 was told to fix

1. `engine_dimension_pin.py` printed *"is a specialisation artefact, not a physics disagreement"*
   **unconditionally** — it still printed when the residual proved a genuine disagreement.
2. The pin's completeness guards were untested: stripping them let it pass while comparing nothing, and
   the whole test population still passed.
3. `README.md` asserted the `D`-coefficient gap was not closed — **by the same round that closed it.**

## ⭐⭐ What to attack

**① Is the printer's characterisation now COMPUTED?** ⭐ Drive it with a residual that does **not** vanish
at the declared binding and confirm no exculpatory prose survives. ⛔ A status typed rather than computed is
the defect; ⭐ check the whole file, ⛔ not only the line that was reported.

**② Do the new completeness tests bind — including in COMBINATION?**
⚠ **Round 3's four single-guard weakenings each failed their test while the COMBINED removal passed
everything.** ⭐ Build the combined weakening yourself, and build one the artifact's own battery does not
contain. ⛔ A test that "covers" an invariant on one example is a demo, not a pin.

**③ Is every claim in `README.md` and the report TRUE?**
⚠⚠ **This paragraph has now been wrong TWICE, in OPPOSITE directions** — first claiming coverage that did
not exist, then denying coverage that did. ⭐ Check every coverage sentence against what you measure.
⭐ Specifically verify: `B_comp` corroborated by one engine only; `mu_R`/`B_comp` attribution carrying no
independent content because their laws are identical; and the common-mode class being genuinely open.

**④ Did anything regress?** ⭐ Re-run every engine the import fence discovers, plus the four named engines,
and diff against committed outputs. ⭐ Compare the test population against a pristine `HEAD`.
⚠ `WL_S10_RUNTIME_SECONDS` moves every run. ⚠ `wolframscript` writes config warnings to **stderr** —
⭐ separate the streams. ⛔ One Mathematica kernel at a time (two seats); ⛔ `timeout 900` each.

**⑤ Rule 2, across the round's diff.** ⛔ Any status typed rather than computed? ⛔ Any `assert` before the
value it guards? ⛔ Any expected value typed beside the value it checks? ⚠ A denylist grep was removed this
round — ⭐ confirm nothing equivalent replaced it.

## ⭐ What is already established — ⛔ do not re-litigate, ⭐ but do report if you find it false

- The pin's operands are **independent**: 11 pairs, 0 circular, established by two legs **by mutation**.
- The laws match five independently-built engines' own symbolic derivations.
- The pin **fails** on a law that agrees at the declared `D` and is wrong elsewhere.

## Method

⭐ Read the sources of truth first, form your own view, ⭐ **then** the diff, then the report.
⛔ **Do not modify the working tree** — it holds substantial uncommitted work.
⛔⛔ **Never `git checkout`, `git stash`, or `git restore`.**
⛔ **Absolute paths for everything outside the repository.** ⚠ A `cd` into a temp directory has already
failed silently in this session and edited live files. ⭐ Baseline with `git archive HEAD | tar -x -C <abs>`.
⭐ Save every script and its literal stdout to named absolute paths.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"a claim in the README is false"* · *"the combined weakening still survives"* · *"the
printer still exculpates"* · *"this regressed."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.
⭐ **If nothing survives, say so plainly — this is a commit gate and a clean result is an actionable one.**

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
