# Independent review — the W0 builder and integrator statements

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifacts — Codex-written

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_builder_statement.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_integrator_statement.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_restatement_verification_report.md`

⭐ The **builder statement** will be read by **four independent engine builds** (S9-WL, S9-Python, S10-WL,
S10-Python). ⛔⛔ **An error in it makes all four engines wrong the same way**, which defeats dual-engine
verification by construction.

⚠ **History you need:** an earlier orchestrator-written directive for this same work
(`W0_emit_speed_itself.md`) was reviewed by two legs and rejected — wrong diagnosis, leaked target,
impossible regression bar, and a compliant build that could emit a frequency instead of a speed. ⭐ These
statements are the replacement. ⛔ Do not assume the replacement fixed what it claims to fix.

## ⭐⭐ What to attack

**① The leak.** ⛔ The builder statement must not let a builder infer what any check will report or what
the emitted object equals. ⭐ Read it as an adversary iterating toward a green run: what can you infer?
⚠ It names `μ_R` and `ρ_br` as quantities held fixed by a mutation. ⭐ Judge whether that is a necessary
input or a leaked answer. ⛔ A grep for forbidden tokens is **not** a leak audit — the previous acceptance
was exactly that and it is dodgeable.

**② Buildability.** ⭐ The certificate demands a selected set of **cardinality one**. ⭐ Check against the
committed outputs whether the stated predicate actually achieves that in **every cell the scope table
requires** — S9's cells and S10's `D=2,3,4,5`. ⛔ If it does not, all four builds emit failure objects and
W0 stalls. ⭐ Verify per cell, ⛔ not on one example.

**③ Can a compliant build still be worthless?** ⭐ Construct one. The statement pins the object by
dimension `[1,-1,0]` and by a defining link residual. ⭐ Find the implementation that satisfies every
sentence and still tests nothing — or show that the previous four worthless modes (registry recombination,
selection by known answer, unestablished sign, emitting a frequency) are each genuinely excluded **by the
words as written**.

**④ Are the mutations real tests?** ⭐ The post-division sentinel is supposed to separate an honest object
from one recombined from fixed inputs. ⭐ Work out whether it does — and whether a builder could satisfy it
while still recombining. ⭐ The selection mutation must make the selected root, `v_T²` and `c_T` all differ:
⛔ can it pass vacuously in a cell where only one candidate exists?

**⑤ Is the regression bar satisfiable AND sufficient?** ⭐ It excludes runtime and inventory/count metadata
from a deterministic projection, then requires an explicit manifest and a verified count delta. ⭐ Check
both halves against the actual committed outputs: can it be met, and does it catch a bad emission hiding
inside a legitimate addition?

**⑥ Four-builder divergence.** ⭐ Is anything ambiguous enough that two builders would reasonably do
different things, manufacturing a cross-engine disagreement that is a specification artifact rather than
physics? ⚠ That failure has already occurred in this repository. ⭐ Check the scope table, the certificate
field list, and whether the emitted records can later be bound by one harness row.

**⑦ The split.** ⭐ Does the integrator statement carry everything the builder statement must not, and does
the builder statement stand alone without it? ⭐ Does anything load-bearing fall in the gap between them?

**⑧ The verification report.** ⭐ It claims D1–D8 verified and no decision-list error found. ⚠ An earlier
pass on the same decision list **did** find a false claim, which was then corrected. ⭐ Spot-check its
claims against the repository. ⛔ Do not accept "verified" as evidence.

## Method

⭐ Read the sources of truth **first** — `reduction/relations.yaml`, `reduction/quantities.yaml`, both
`checks_S*.yaml`, the four engines, and the committed outputs under `mathematica/out/` and `scripts/out/` —
and form your own view of what the object must be and what the statement must require. ⭐ **Then** read the
artifacts. ⛔ Not in the other order.

⭐ `W0_decision_list.md` is the orchestrator's input to this build. ⚠ Read it **after** forming your own
view; ⛔ it is an input under test, not a source of truth.

- ⛔ Do not modify the repository. ⛔ Do not write or run any engine.
- ⭐ Where you claim a builder could comply and still be wrong, **exhibit the complying implementation**.
- ⭐ Where you check a per-cell claim, give the literal command and its output.
- ⛔ Wrap any long run in `timeout 600`; ⛔ never raise it; ⛔ never more than one Mathematica kernel at a
  time — the licence has **two** seats.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"the cardinality-one requirement fails in cell X"* · *"this leaks the target"* · *"a
compliant build is still tautological"* · *"two builders would reasonably diverge here"* · *"the regression
bar cannot be met."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
