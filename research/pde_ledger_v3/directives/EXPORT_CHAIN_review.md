# Independent review — the export-chain architecture proposal

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — orchestrator-written

`/var/projects/toy_physics/research/pde_ledger_v3/EXPORT_CHAIN_ARCHITECTURE.md`

⭐ A proposal to **replace** the registry-and-relations machinery with a dataflow chain: Python engines
import the previous step's export record, Wolfram engines stay siloed and derive everything independently,
and the disagreement between them is the check.

⛔ **Nothing has been built or deleted.** ⭐ This review decides whether it gets built.

## The context you must verify, ⛔ not accept

The proposal rests on a claim about the current state: that of five relations in `reduction/relations.yaml`,
only `R4` is wired into an executable check, and that `R4`'s binding is wrong.
⭐ **Check it end to end** — `reduction/relations.yaml`, both `reduction/checks_S*.yaml`, the four engines,
and the committed outputs under `mathematica/out/` and `scripts/out/`.
⛔ If that diagnosis is wrong, the proposal's premise fails and that is the most valuable thing you can
report.

## ⭐⭐ What to attack

**① Does the asymmetry actually work?** ⭐ The claim is that a wrong export propagates through the Python
chain but the siloed Wolfram engine disagrees. ⭐ Test that reasoning against the **actual** engines:
does each `.wl` engine independently derive the objects the `.py` chain would carry? ⛔ Where it does not,
the export is **unchecked** — ⭐ enumerate those cases from the committed outputs and say how large the gap
is. ⚠ The proposal claims the catch happens at the step of origin rather than at each step of use;
⭐ verify that is true for the objects that actually cross steps.

**② Is anything lost that the registry provided?** ⭐ Be concrete. The registry is claimed to reconcile a
derivation against a hand-copied restatement of itself, with no quantity having two independent sources.
⭐ **Check that claim** — is there any quantity in `quantities.yaml` genuinely derived in two places? ⛔ If
there is, the proposal's root-cause argument is wrong for that case.

**③ Is `import ⊄ export` sufficient to protect blind derivation?** ⭐ Construct the way around it. A step
imports what it consumes and derives what it produces — ⭐ can a builder import something that *determines*
its output without literally being it? ⛔ If so, the rule is a spelling check.

**④ Is the per-step cost estimate honest?** ⭐ The proposal claims the recurring cost is two short lists,
and that the comparison list is mechanical because S11's engines share 190 of ~225 tag suffixes.
⭐ Verify that number and judge whether it generalises — ⚠ S10's list is 3,121 lines and the proposal
attributes that to S10 predating a shared spec. ⛔ Is that the real reason, or is S10 structurally harder?

**⑤ Is deleting the homogeneity gate right?** ⭐ The proposal cites three measurements showing it passes
under a wrong `D` coefficient, under deleted laws, and under a common-mode shift. ⭐ Verify at least one.
⭐ Then judge the replacement claim: does per-engine homogeneity, computed from each engine's own
derivation, actually cover what the gate covered? ⛔ Name what is lost.

**⑥ Migration.** ⭐ Do S9, S10 and S11 already emit enough to form an export record **without engine
changes**? ⭐ Answer from the committed outputs, per step, and name what is missing. ⚠ The proposal asserts
this is an open question — ⭐ close it if you can.

**⑦ What does this cost that the proposal does not mention?** ⭐ The honest failure mode of a dataflow chain
is a silent dependency. ⭐ Name the costs it omits, and any place where it is optimistic about work that has
already proven expensive here.

## Method

⭐ Read the repository first and form your own view of what the cross-step requirement needs. ⭐ **Then**
read the proposal. ⛔ Not in the other order.

- ⛔ Do not modify the repository. ⛔ Do not build anything.
- ⭐ **Absolute paths for everything outside the repository.** ⚠ A `cd` into a temp directory has failed
  four times in this session, once editing live files.
- ⭐ Where you claim something is or is not emitted, give the literal command and its output.
- ⛔ One Mathematica kernel at a time (two seats); ⛔ `timeout 900` on each.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"the diagnosis of the current state is wrong"* · *"the WL silo does not actually check
object X, so the chain is unguarded there"* · *"a quantity IS derived twice and the root-cause argument
fails"* · *"the cost estimate is wrong"* · *"deleting the gate loses Y."*

⚠ **This proposal deletes working machinery.** ⭐ A finding that it deletes too much is as valuable as one
that it keeps too much.
⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
