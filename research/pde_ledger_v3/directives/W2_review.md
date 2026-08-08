# Independent review — W2 registry-dimension witness

⚠ **Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.** A warning, not a
control.

## Artifact — Codex-written

- `reduction/registry_dimension_witness.py` (679 lines)
- `reduction/registry_dimension_witness.yaml` (41)
- `reduction/registry_dimension_witness_able_to_fail.py`
- `reduction/test_registry_dimension_witness.py`

**What it claims to do:** compare each engine's **emitted** dimension vector for a registry quantity
against that quantity's **declared** exponents in `reduction/quantities.yaml` — an oracle outside any
step's shared specification. Nothing in the project did this before.

**What it reported on committed S9/S10 outputs:**
`AGREEMENT=9 · DISAGREEMENT=0 · CONVENTION_MISMATCH=0 · UNDETERMINED=0 · NOT_EMITTED=47`

## ⭐⭐ What to attack

### ① Is the convention check REAL, or an assertion about the engine?

The manifest declares `emitted_convention: LTM-exponent-vector-v1` per artifact, and its own comment says
*"Engine output carries bare triples, so its convention is an explicit input to this witness."*

⭐ **So: would this check catch an engine that actually changed its axis order?** The able-to-fail case
perturbs the **manifest's** declared convention. ⛔ That proves the comparison fires; it does **not** prove
the engine's real order is measured. ⭐ Determine what it would take to make this a measurement rather than
a declaration, and whether the current report over-states what was verified.
⚠ This project has a **measured** `M,L,T` vs `(L,T,M)` engine disagreement on record (`STATUS.md`,
stage037). ⭐ Construct the case and say whether the witness catches it.

### ② Is the row selection independent of the result?

Selection inherits `dimension_sources` from `checks_S9.yaml` / `checks_S10.yaml`, filtered to
`packages: [MAIN]` (and `dimensions: [3]` at S10), plus a hand-added `extra_dimension_sources` entry.

⭐ **Does that filter exclude rows that would disagree?** Specifically: the controls (`X1`–`X8`,
`XFORM_*`) are excluded by `packages: [MAIN]` — ⭐ is that legitimate, or does it drop comparable rows?
⛔ A selection authored by the same party that then reports 0 disagreements is a control inside the thing
it polices — ⭐ test it, don't assume it either way.

### ③ Re-derive what SHOULD be comparable

⭐ Independently of the manifest, work out from `quantities.yaml` and the committed outputs **which
registry quantities have an emitted dimension vector somewhere.** Compare your list to the witness's 9
`AGREEMENT` rows and its `NOT_EMITTED` list.
⭐ **Is `NOT_EMITTED=47` honest**, or does it absorb quantities that *are* emitted under a tag the manifest
did not look at? ⚠ Report any quantity you can compare that the witness does not.

### ④ Does the `D`-specialisation hide anything?

Engines emit dimensions symbolic in `D` (`[-D,0,1]`, `[2-D,-2,1]`); the registry declares the `D = 3`
specialisation. The witness substitutes `D ← Q.brane.D_brane`.
⭐ Could a genuinely wrong emitted vector become right under that substitution? ⭐ Construct one if so.
⛔ And check: is the substitution ever applied when the emitted vector is **not** actually `D`-symbolic?

### ⑤ Do the able-to-fail cases test the WITNESS or the CALIBRATOR?

Three cases pass with `deviation_count=0`. ⭐ Verify each perturbs the thing it claims to and that the
calibrator is not the only thing asserting the transition. ⛔ Run them yourself and paste literal stdout.
⭐ Then find a defect class none of the three would catch.

### ⑥ Rule 2 compliance

⭐ Does it **print computed objects** — both operands and the residual — or does it state conclusions?
⭐ Is there an `assert` before any value it guards? ⭐ Is any status typed rather than computed?

## Required method

⭐ Read `reduction/quantities.yaml`, `checks_S9.yaml`, `checks_S10.yaml`, the committed outputs under
`scripts/out/` and `mathematica/out/`, and `dimensional_homogeneity_gate.py` **first**. Form your own view
of what this check should compare. ⭐ **Then** read the artifact.

⛔⛔ **A prose assertion is worth nothing.** Every finding needs a quotation with a file:line locus, or a
literal command and its literal output. ⭐ Where you claim the check can be defeated, **exhibit the
defeating case and paste its output.**

⛔ **Do not modify the repository.** Copy to a temp dir outside `/tmp/claude-1000/` if you need to run
mutations.

## Physics filter

> Report a finding only if it changes **what the project computes, what it may claim, or what method it
> adopts.**

⭐ Worth most: *"this check cannot catch X"* · *"this row selection excludes the disagreements"* ·
*"AGREEMENT=9 is weaker than it reads"* · *"the convention check verifies nothing about the engine."*

⚠ A leg returning "nothing survives" is weak evidence alone — state what you checked, what you could not,
and what would have had to be true for you to find something. ⛔ Don't manufacture findings.

## Report

Ranked most-severe first: the claim · the evidence (quotation + locus, or literal output) · what must
change. ⭐ Save every script and its stdout to named absolute paths and report them.
