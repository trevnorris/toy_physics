# Independent review — W3 round 3 and the W4 engine-dimension pin

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — Codex-written, fourth round on this material

`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` — chiefly `engine_dimension_pin.py`,
`registry_import_fence.py`, `dimension_law_check.py`, `test_dimension_laws.py`, `w3_acceptance_ablations.py`,
`dimension_law_able_to_fail.py`, `README.md`, and the retirement of `generate_rows.py`.
Report: `reduction/W3_REGISTRY_D_LAWS_REPORT.md`. Diff: `git diff HEAD -- research/pde_ledger_v3/reduction/`
**and** `git status` — some changes are staged, some are not, some files are untracked.

⚠⚠ **The same defect has moved three times in this file set** — a typed copy of the expected answer, first
in the production checker, then in a test table. ⭐ **Attack the repairs**; ⛔ do not assume the fourth
round removed rather than relocated it.

## What the pin is for

`quantities.yaml` declares three brane quantities' dimensions as **laws in `D`**. Nothing inside
`reduction/` could distinguish a correct law from one that agrees at the registry's own `D` and is wrong
elsewhere — the relation gate constrains only **differences**, and `D` cancels out of every difference.
The pin compares each engine's **emitted symbolic dimension** against the registry's **law**.

## ⭐⭐ ① THE QUESTION THAT DECIDES EVERYTHING — is the pin independent, or circular?

⭐ The pin's whole value is that its expected operand comes from **outside** the registry. ⛔ Verify that,
per engine, ⛔ do not assume it.

⭐ Establish for **each** transcript the pin reads: does that engine compute its emitted symbolic dimension
from its **own** derivation, or could it have obtained it from the registry it is being compared against?
⚠ Some of these engines import the registry loader and some do not — ⭐ **measure which, and trace the tag
back to what produced it.** ⛔ A pair where the engine read the registry is the registry agreeing with
itself, and must not be counted as corroboration.

⭐ Report the pin's population split into **independent** and **potentially circular** pairs, with the
evidence for each. ⭐ If some pairs are circular, say whether what remains still pins the coefficient.

## ⭐ ② Is the typed duplicate gone, or has it moved a fourth time?

⭐ Search all of `reduction/` for any typed expected dimension, exponent vector, or law value sitting beside
the thing it checks — production code **and** tests. ⭐ For each one you find, decide whether it is
derivable from an independent object or is a second copy of the answer.
⚠ Note `test_dimension_laws.py` still contains constant triples for the two speeds. ⭐ Judge whether those
are forced by structure or are duplicates.

## ⭐ ③ Does the pin bind, and can it be defeated?

⭐ Run it. Then attack it: ⛔ can it pass while comparing nothing? ⭐ Construct a manifest or transcript set
that yields zero pairs and confirm that is guarded. ⭐ Can a tag be silently missing, or matched to the
wrong quantity, and still report `PASS`? ⭐ Does it fail for the right reason, or only because a parse
broke?
⭐ **Exhibit at least one wrong law it does NOT catch**, or show that class is empty.

## ⭐ ④ The import fence

⭐ It is meant to discover, by import analysis, every engine that consumes the registry — replacing a
hand-maintained list that had already missed one. ⭐ Verify the set it discovers is right: ⛔ can an engine
consume the registry without the fence seeing it (indirect import, dynamic load, reading the YAML directly)?

## ⭐ ⑤ The recorded claims

⚠ A previous round's `README.md` asserted the gap in absolute terms and was **false**. ⭐ Read every claim
now made about coverage in `README.md`, `W3_REGISTRY_D_LAWS_REPORT.md` and `REBUILD_HANDOFF.md`, and check
each against what you measure. ⭐ Worth most: a claim of coverage that does not exist, **or** a claim of a
gap that is not a gap.
⭐ The S11 transcript carries a standing nonzero residual now recorded as a `D = 3` specialisation artefact.
⭐ Verify that characterisation — ⛔ it must not be a physics disagreement being explained away.

## ⭐ ⑥ Weaker implementations and rule 2

⛔ A test that "covers" an invariant demonstrates it on one example. ⭐ Build the weaker implementations of
each object this round claims to fix and show whether the tests fail on each — ⭐ including weaker versions
the artifact's own battery does **not** contain.
⭐ Does every script print computed objects — both operands and the residual — then guard? ⛔ Any `assert`
before the value it guards? ⛔ Any status typed rather than computed? ⚠ A previous round shipped four typed
literals formatted exactly like computed operand lines.

## ⭐ ⑦ Regression

⭐ Re-run every engine the fence discovers, plus the four named engines, and diff against committed outputs.
⚠ `WL_S10_RUNTIME_SECONDS` moves every run. ⚠ `wolframscript` writes configuration warnings to **stderr** —
⭐ separate the streams or you will read them as a diff. ⛔ One Mathematica kernel at a time — the licence
has **two** seats; ⛔ wrap each run in `timeout 900`. ⭐ Compare the test population against a pristine `HEAD`.

## Method

⭐ Read the sources of truth first — `relations.yaml`, `quantities.yaml` at HEAD, the engines, the committed
outputs — form your own view of what the laws must be and what an honest pin can claim. ⭐ **Then** read the
diff, then the report.

⛔ **Do not modify the working tree.** ⚠ It holds substantial uncommitted work.
⛔⛔ **Never `git checkout`, `git stash`, or `git restore`** — an uncommitted baseline dies in the revert.
⛔ **Use absolute paths for everything outside the repository.** ⚠ A `cd` into a temp directory has already
failed silently in this session and edited live files. ⭐ Extract a baseline with
`git archive HEAD | tar -x -C <absolute-dir>`.
⭐ Save every script and its literal stdout to named absolute paths and report them.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"pair X is circular — the engine read the registry"* · *"a typed duplicate remains"* ·
*"the pin passes while comparing nothing"* · *"this claim of coverage is false"* · *"the specialisation
artefact is a real disagreement."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
