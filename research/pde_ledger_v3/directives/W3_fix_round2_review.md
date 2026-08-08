# Independent review — W3 fix round 2

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — Codex-written, third round on this material

`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` — the bound-dimension-law change and its
checks. Diff with `git diff HEAD -- research/pde_ledger_v3/reduction/` **and** `git status` (some changes
are staged, some are not, some files are untracked). Report: `reduction/W3_REGISTRY_D_LAWS_REPORT.md`.

⚠⚠ **In this repository, repairs have repeatedly introduced NEW defects into the material just changed —
including in this very file set, twice.** ⭐ **Attack the repairs.**

## What the registry change is for

`quantities.yaml` declares each quantity's dimension as an LTM exponent vector. Three brane quantities live
on a brane whose dimension `D` is itself a registry quantity, so their true exponents are **functions of
`D`**; the schema previously held only three bare integers and silently recorded the `D = 3` evaluation.
W3 lets a declaration be a **law** bound to a structural quantity.

## ⭐⭐ The history you must not repeat

⭐ **Round 1** was reviewed by two legs. Two blockers were found **after** it reported success:

1. ⛔ **It broke an engine.** `S10-python` died with `TypeError: 'BoundDimensionLaw' object is not
   iterable`. ⚠ **Round 1's regression bar was all-green and byte-identical across six checks** — gate,
   loader, acceptance, able-to-fail, `engine_output_checks.py` on both configs, and the full test
   population. ⛔ **None of them runs an engine**; `engine_output_checks.py` reads committed `.out` files.
2. ⛔ **Its new "pin" was a second copy of the answer** — the expected laws hardcoded in Python beside the
   registry's. Changing **both** copies to the same wrong law left it green.

⇒ ⭐ Round 2 was told to repair both, to re-run **all four engines** as its bar, and to **record** rather
than close the fact that the `D` coefficient is not policed inside `reduction/`.

## ⭐⭐ What to attack

**① Does the compatibility view hide anything?** ⭐ `Quantity.dimension` is now a reference-vector view
while the law lives elsewhere. ⭐ Check every consumer: which of them now silently reads a `D = 3`
specialisation where it should see the law? ⛔ The old defect was a silent specialisation; ⭐ verify this
repair did not reintroduce it one layer down under a new name.

**② Is the replacement check honest, and does it still bind?** ⭐ It now reports the `D` coefficient as not
policed and fails. ⭐ Confirm the sub-verdicts still discriminate — construct the cases and show the output.
⛔ A check that always fails is a check nobody reads: ⭐ say whether a **real** new defect in this script
would be visible against a standing red.

**③ Is the duplicated expectation genuinely gone?** ⭐ Search for any remaining typed copy of an expected
dimension anywhere under `reduction/`. ⛔ A control whose expected value is typed beside the value it checks
is the defect; ⭐ a renamed one is the same defect.

**④ Re-run all four engines yourself.** ⭐ `S9-py`, `S10-py`, `S9-wl`, `S10-wl`, diffed against the committed
outputs. ⛔ Do not take the report's word. ⚠ `wolframscript` writes configuration warnings to stderr —
⭐ separate stderr from stdout or you will read them as a diff. ⚠ `WL_S10_RUNTIME_SECONDS` moves every run.
⛔ One Mathematica kernel at a time — the licence has **two** seats; ⛔ wrap each in `timeout 900`.

**⑤ Is the owed pin recorded truthfully?** ⭐ The real pin — comparing each engine's emitted **symbolic**
dimension in `D` against the registry law — is deferred to the witness rebuild. ⭐ Check that what is
written about it is accurate, that nothing claims more coverage than exists, and that a reader meeting a
green run cannot conclude the coefficient is checked.

**⑥ Weaker fixes.** ⛔ A test that "covers" an invariant demonstrates it on one example. ⭐ Build the weaker
implementations of each object this round claims to fix and show whether the tests fail on each.
⚠ **Measured three times here:** a half-fix passed the new test, the whole suite, the full battery **and**
produced byte-identical output.

**⑦ Rule 2.** ⭐ Does every script print computed objects — both operands and the residual — then guard?
⛔ Any `assert` before the value it guards? ⛔ Any status typed rather than computed?

## Method

⭐ Read the sources of truth first — `relations.yaml`, `quantities.yaml` at HEAD
(`git show HEAD:research/pde_ledger_v3/reduction/quantities.yaml`), the engines, the committed outputs —
and form your own view of what the declarations should be and what each check can honestly claim.
⭐ **Then** read the diff, then the report.

⛔ **Do not modify the working tree.** ⚠ It holds substantial uncommitted work.
⛔⛔ **Never `git checkout`, `git stash`, or `git restore`** — an uncommitted baseline dies in the revert.
⛔ **Use absolute paths for everything outside the repository.** ⚠ A `cd` into a temp directory has already
failed silently in this session and edited live files. ⭐ Extract a pristine baseline with
`git archive HEAD | tar -x -C <absolute-dir>`.
⭐ Save every script you write and its literal stdout to named absolute paths and report them.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"a consumer still gets a silent D=3 specialisation"* · *"an expected value is still typed
beside the value it checks"* · *"an engine does not regenerate"* · *"this claims coverage it does not
have"* · *"the weaker fix passes."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
