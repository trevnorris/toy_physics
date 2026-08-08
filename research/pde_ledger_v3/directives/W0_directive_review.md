# Independent review — the W0 build directive

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — orchestrator-written

`/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_emit_speed_itself.md`

It is **the single shared statement that four independent engine builds will each read** (S9-wl, S9-py,
S10-wl, S10-py). ⛔⛔ **An error in it makes all four engines wrong the same way, which defeats dual-engine
verification by construction.** That is what this review exists to catch.

## Context you must verify for yourself, ⛔ not accept

The directive claims `R4` has never been tested at S9 or S10 because the checks bind a **squared** speed to
a quantity the registry declares as a **linear** speed. ⭐ Check that claim end to end:

- `research/pde_ledger_v3/reduction/relations.yaml` — `R4` and `R5`
- `research/pde_ledger_v3/reduction/quantities.yaml` — `Q.brane.c_gamma`, `Q.brane.c_L`
- `research/pde_ledger_v3/reduction/checks_S9.yaml` and `checks_S10.yaml` — the `registry_residual` blocks
- the engines: `scripts/S9_light_requires_shear_sympy_audit.py`,
  `mathematica/S9_light_requires_shear_mathematica_audit.wl`,
  `scripts/S10_brane_mode_spectrum_sympy_audit.py`,
  `mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
- the committed outputs under `mathematica/out/` and `scripts/out/`

⭐ Is the diagnosis right? Is it **complete** — does the same defect exist anywhere the directive does not
name (other relations, other quantities, other steps, the `D≠3` cells S10's binding never reaches)?
⛔ If the diagnosis is wrong, that is the most valuable thing you can report.

## ⭐⭐ What to attack

**① Does the directive name an OBJECT, or does it specify a RECIPE?**
⚠ **Measured in this repository:** an orchestrator decision that specified a *mechanism* instead of naming
an *object* caused a builder to implement a short-circuit that left a planted defect invisible across an
entire step. ⭐ Find every place this directive tells a builder *how* rather than *what*, and say what the
object should have been.

**② Can a builder satisfy every word of it and still produce a worthless emission?**
⭐ Construct the compliant-but-useless implementation and describe it concretely. The directive names three
failure modes it intends to block — ⭐ **find a fourth**, or show one of the three is not actually blocked
by the words as written.

**③ Is the "not tautological" requirement enforceable, or only exhorted?**
The directive demands the emitted speed descend from the engine's own dispersion relation rather than from
`μ_R/ρ_br` recombined. ⭐ Could a reviewer of the resulting engine actually **tell the difference** from the
artifacts the directive asks for? If not, the requirement is decoration and the directive must ask for
something that makes it visible.

**④ Does it leak an expected result?**
⛔ A builder iterating to a green run converges on any target it can see. ⭐ Does the directive anywhere
state, imply, or let a builder infer what the residual should equal, which candidate is the transverse one,
or what the speed is? ⚠ The `R4` residual formula is quoted in it — ⭐ judge whether quoting the relation
the builder must **not** reproduce by hand is itself the leak.

**⑤ Is the regression bar sound?**
It demands the re-run diff contain only additions. ⭐ Is that achievable — do these engines emit anything
run-dependent (timing, ordering, hashes) that makes a clean diff impossible? ⭐ And is "only additions"
sufficient, or can a real defect hide inside an addition-only diff?

**⑥ Scope and consistency.**
⭐ Is the exclusion list right — does the change genuinely not force edits to things it forbids touching?
⭐ Four builders read this one file: is anything ambiguous enough that two of them would reasonably do
**different** things, manufacturing a cross-engine disagreement that is a specification artifact rather
than physics? ⚠ That exact failure has already occurred here once.

## Method

⭐ Read the sources of truth **first** and form your own view of what the defect is and what the fix must
require. ⭐ **Then** read the directive. ⛔ Do not read them in the other order — reading the artifact first
anchors you to its framing, which is the thing under test.

⛔ Do not modify the repository. ⛔ Do not write or fix the engines — this is a review of the **directive**.
⭐ Where you claim a builder could comply and still be wrong, **exhibit the complying implementation**.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"the diagnosis is wrong or incomplete"* · *"a compliant build is still tautological"* ·
*"this leaks the answer"* · *"two builders would reasonably diverge here."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
