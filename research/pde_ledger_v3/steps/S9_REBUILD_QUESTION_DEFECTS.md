# S9 REBUILD — defects found in the QUESTION, before the build ran

**Rewritten 2026-08-04, replacing a prediction table.** ⚠ The first version of this file was a table of
**expected values** for the rebuilt engine to be checked against. ⛔ **That was the wrong artifact and it
has been removed rather than annotated.**

## ⛔⛔ WHY, because the reasoning generalises to every remaining step

⭐⭐ **The entire point of writing the script is to GET the answer. The orchestrator's job is to make
sure the correct QUESTION is asked** (user, 2026-08-04).

A committed table of expected values is:

- **an acceptance criterion referencing an expected value** — the one thing the rebuild withholds,
  because a builder iterates until it exits 0, and if *"matches the recorded number"* is reachable it
  will get there, turning a genuine disagreement (the most valuable output available) into silent
  confirmation;
- **not evidence of anything.** It was derived **by hand**. ⚠ A hand derivation agreeing with another
  hand derivation is [[feedback-matching-number-is-not-evidence]] operating at full strength: two
  unverifiable claims that happen to coincide, which is exactly how two cancelling errors survive.

⛔⛔ **AND THE SAME APPLIES TO A REVIEW LEG'S "INDEPENDENT DERIVATION"** (user, 2026-08-04):
> *"trusting grok and codex and even yourself on how they 'rederive' is not trustworthy. Unless it's in
> CAS and we can see the output from the inputs, it's not to be trusted."*

⇒ ⭐ A leg reporting *"I derived it independently and got X"* in prose is **a typed sentence with no
computation behind it** — the identical defect class this whole rebuild exists to remove, relocated from
the engine into the review. ⇒ **A review leg claiming a derivation must produce a SCRIPT and its LITERAL
OUTPUT at a path the orchestrator can re-run.** Prose is not a derivation.

---

## ⭐ What survives: the question-quality defects found while working the setup

These were found by working through what the directive **asks**, ⛔ not by predicting what it would
answer. Each was a question that could not be honestly answered as posed. All four are repaired in
`directives/S9_wl_rebuild_directive.md`.

| # | the bad question | why it could not be honestly answered | found by |
|---|---|---|---|
| **1** | control **X4** asked for *"the difference between the transverse root and the non-transverse root"* | X4 collapses the root set, and the directive defined **transverse** (`M·T=0`) while never defining **longitudinal** at all ⇒ *"non-transverse"* was an **undefined term**, so the only honest output was an **invented** value ⇒ [[feedback-no-fabrication-forcing-rules]] | orchestrator (partially) + Grok leg (sharper diagnosis: `M = 0` at that root, so *every* root passes) |
| **2** | tag 14 asked for the residual `ω² − speed²·k²` after tag 13 **defined** `speed² := ω²/k²` | identically zero **by construction** ⇒ a tautology wearing the costume of a check. ⚠ Written into a directive whose own §0 bans exactly this pattern | Grok leg |
| **3** | tag **names** were unconstrained | a name can state the conclusion while the payload stays a genuine CAS object ⇒ §0 closed the payload slot and left the name slot open | Grok leg |
| **4** | `D[speedSquared, k]` differentiated with respect to a **vector** `k` | ill-defined as written; needs a scalar `kMag` | Grok leg |

⭐ **The repair for 1 that works in every case:** emit the **complete root multiset with multiplicities
and the nullity of `M` at each root**, for the main derivation and every control. It stays well posed
when roots collapse, when a root disappears, and when a polarisation subset is empty.

⭐ **The repair for 2 is to DELETE the check, not to strengthen it.** At S9 there is genuinely **no
independent second route** to the propagation speed — the script may not read the registry and there is
only one engine. ⇒ The directive now says to emit the objects and **report that no second operand
exists**. ⭐ An honest *"there is nothing to compare against here"* beats a zero that looks like a check.

## ⚠ Leak-gate result, recorded because it caught the same trap TWICE

⛔⛔ **A PROHIBITION LEAKS THE ANSWER JUST AS SURELY AS AN ASSERTION.** Both hits were inside sentences
**forbidding** the builder to state the answer:

- `Print["TAG: LINEAR_IN_K"] is forbidden` — named the shape of the dispersion while banning it;
- `WL_S9_..._MU_OVER_RHO is forbidden` — named the ratio while banning it. ⚠ This one was introduced
  **by the repair for defect 3**, i.e. while I was actively fixing a leak.

⇒ ⭐ **Write forbidden-pattern examples with placeholder content, never with this step's real content.**

## ⛔ What this build does NOT test

- **P2** — that a scalar superfluid carries no transverse mode. Cited, never computed; this build does
  not compute it either, and the step record must say so.
- **P4** — bulk shear-freeness. The bulk is absent from the action, so nothing here can test it.
- Wall width, background flow, dissipation, frequency-dependence of the moduli, amplitude — all sent to
  zero (directive §6).
- ⛔ S9 has **one** engine, so **no cross-engine disagreement is available at this step.**
