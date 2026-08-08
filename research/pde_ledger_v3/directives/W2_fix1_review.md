# Independent review — W2 witness, after fix round 1

⚠ **Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — Codex-written, second round

`reduction/registry_dimension_witness.py` · `.yaml` · `registry_dimension_witness_able_to_fail.py` ·
`test_registry_dimension_witness.py` · the `reduction/README.md` section.

It compares each engine's **emitted** dimension vector for a registry quantity against the **declared**
exponents in `reduction/quantities.yaml` — an oracle outside any step's shared specification.

Round 1 was reviewed by two legs and seven defects were fixed. ⚠⚠ **In this repository, repairs have
repeatedly introduced NEW defects into the material just changed.** ⭐ **Attack the repairs.** Diff:
`git diff research/pde_ledger_v3/reduction/` and `git status`.

**What it now reports on committed outputs (exit 1):**
```
AGREEMENT=3 DISAGREEMENT=0 ECHOED=2 ECHOED_MISMATCH=0 BRANCH_DIMENSION_MISMATCH=4
CONVENTION_MISMATCH=0 UNDETERMINED=4 NO_TAGS=0 NO_ROWS=0 NOT_EMITTED=57 guard_count=8
```

## ⭐⭐ ① The physics-bearing question — adjudicate it yourself

`BRANCH_DIMENSION_MISMATCH=4` fires on `Q.brane.rho_br` and `Q.brane.mu_R` at S10-wl and S10-py. Residuals
walk across S10's brane-dimension branches (`[1,0,0]` at D=2 → `[0,0,0]` at D=3 → `[-1,0,0]` → `[-2,0,0]`).

⭐⭐ **Is this an ENGINE defect, a REGISTRY defect, a WITNESS defect, or none of them?** ⛔ Do not accept
any framing — derive it. Work out from the model what `[rho_br]` and `[mu_R]` *should* be on a
D-dimensional brane, check what each engine emits, check what `quantities.yaml` declares, and say who is
wrong. ⭐ If your answer is "nobody is wrong," say what is missing and where it belongs.

⭐ Then: is `BRANCH_DIMENSION_MISMATCH` the right **status** for it — guarded, blocking, exit 1 — or is the
witness reporting a modelling gap as though it were a failure?

## ② Is the multiplier safe?

`c_gamma` and `c_L` are now compared via a **declared multiplier of 2** (their squares are what the engines
emit), giving residual `[0,0,0]`.

⛔ **A declared multiplier is exactly where a wrong declaration manufactures agreement.** ⭐ Check: is `2`
correct for each source it is declared on? Could a wrong multiplier hide a real disagreement — construct
one. ⭐ Is the multiplier printed on every row alongside declared, expected, emitted and residual, so a
reader can audit it?
⭐ And verify independently: does `[mu_R] − [rho_br]` really give `[2,-2,0]` with `D` cancelling, so that
`[c_gamma] = [1,-1,0]` holds for **all** D? Check the same for R1, R2.h0, R2.xi_h, R5.

## ③ Did `AGREEMENT` fall from 9 to 3 honestly?

⭐ Account for every one of the six rows that stopped being `AGREEMENT`. For each: is the new status right,
or has a real agreement been reclassified away? ⛔ A number that drops is not automatically more honest.

## ④ Do the new statuses bind?

- **`ECHOED`** — sources where the engine re-emits a registry declaration. ⭐ Is the DERIVED/ECHOED
  classification correct for every source? Can a genuinely derived source be misfiled as echoed, or worse,
  an echo as derived?
- **Coverage invariants** replaced pinned integers (`len(rows)==9`, `all AGREEMENT`, `47`). ⭐ Do they
  actually bind? Construct a manifest that compares nothing and confirm it is now guarded.
- **`UNDETERMINED=4`** — axis order unmeasured where the engine emits no labelled components. ⭐ Is that
  honest, or does it suppress rows that could be compared?

## ⑤ Are four calibrations enough?

Cases: `inherited-config`, `multiplier`, `echoed`, `branch-dimension`. ⭐ Run them and paste literal
stdout. ⭐ Then name a defect class none of the four would catch.

## ⑥ Rule 2

⭐ Does it print computed objects — operands, multiplier, residual — or state conclusions? Any `assert`
before a value it guards? Any status typed rather than computed?

## Method

⭐ Read `quantities.yaml`, `relations.yaml`, `checks_S9.yaml`, `checks_S10.yaml`, the committed outputs and
`dimensional_homogeneity_gate.py` **first**; form your own view of what should be compared and what each
status should mean. ⭐ **Then** read the artifact and the diff.

⛔⛔ **A prose assertion is worth nothing.** Every finding needs a quotation with a file:line locus, or a
literal command and its output. ⭐ Where you claim the check can be defeated, **exhibit the defeating case**
and paste its output. ⛔ Do not modify the repository — copy to a temp dir outside `/tmp/claude-1000/`.

## Physics filter

> Report only what changes what the project computes, what it may claim, or what method it adopts.

⭐ Worth most: *"this status is wrong about who is at fault"* · *"the multiplier hides a disagreement"* ·
*"this repair broke something round 1 had right"* · *"the coverage invariant does not bind."*

⚠ A leg returning "nothing survives" is weak evidence alone — state what you checked, what you could not,
and what would have had to be true to find something. ⛔ Don't manufacture findings.

Rank most-severe first: claim · evidence (quotation + locus, or literal output) · what must change.
