# S11b-B — WRITE REVISION 9 OF THE BUILD DIRECTIVE

⛔⛔ **DOCUMENT TASK. ⛔ Do NOT write, run, or modify any `.wl` or `.py` script.**

You authored revisions 7 and 8. ⭐ **You are the author again for rev 9.**

## Artifacts

**Edit:** `directives/S11bB_SHARED_PHYSICS.md` — and this round **also the two headers**,
`S11bB_wl_header.md` and `S11bB_py_header.md` (finding 2 is a header defect; previous rounds barred header
edits, that bar is lifted **for finding 2 only**).

Then reassemble both directives, verify byte-identity of the shared block, and report the new `sha256`.

## The review to fold

`/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11bB_dir_rev8_agent.txt`

⚠ Findings **1–6**, plus the **Q3 gap** (which carries a concrete covering check). A second reviewer
returned no blocking findings and independently reached the same Q3 conclusion. ⭐ Fix every finding that
survives your judgement; say so if one does not. ⛔ Derive remedies — do not transcribe them.

## ⛔⛔⛔ A CONTAMINATION WARNING — READ BEFORE OPENING THE REVIEW

**The review file contains a reviewer's DERIVED CANDIDATE ANSWER to a question the engines are required to
compute for themselves** — specifically an admissibility conclusion about the closure coefficients.

⛔⛔ **THAT ANSWER MUST NOT APPEAR, BE HINTED AT, BE NARROWED TOWARD, OR BE MADE EASIER TO REACH ANYWHERE IN
THE DIRECTIVE.** ⛔ Do not state which coefficients survive admissibility. ⛔ Do not state what reciprocity
forces. ⛔ Do not add an example, a "note", a suggested form, or a tag whose name encodes an outcome.
⚠ **The directive must leave that question exactly as computable — and exactly as open — as it is now.**
⭐ Use the finding to fix the *specification defect* it identifies, ⛔ never to propagate its conclusion.
⚠ **This will be checked against your diff.**

## ⭐ ONE DECISION ALREADY TAKEN — implement it, ⛔ do not re-litigate

**Finding 5 — the common `τ`.** The project owner has twice chosen, on the same class of question, to
**carry the disputed structure symbolically and let the step compute its consequence.** ⇒ apply that here:

⭐ **Carry a separate relaxation time for each coupling** — `τ_A`, `τ_V`, `τ_X` — each real and `≥ 0`, with
⛔ **`τ_A = τ_V = τ_X` recovering rev 8's specification exactly**, so nothing is asserted either way.
⚠ Propagate them wherever the single `τ` appears: the closure, the face response, the power identity, the
passivity and reciprocity tasks, the dimension list, the validity conditions, and the controls.
⭐ Require the engines to **report whether any relation among the `τ_I` is forced**, and ⛔ to report it as
an output. ⛔ Do not state or hint at what that relation is.

## The remaining findings

- **1 — the read-bar is broken.** `directives/S11b_*` does not match `S11bB_*`, so **every** `S11bB_*` file
  — including the fix directives and review prompts sitting beside the engine's own directive — is
  unbarred, and `_scratch/` and `LAUNCH_PROMPT.md` are not listed at all. ⭐ Fix §0b so the bar actually
  covers what it intends, ⚠ **stating paths in a form that cannot be defeated by a prefix mismatch.**
- **2 — the headers impose different information bars** (one says "no external file reads", the other has
  no equivalent). ⭐ Make the bar **symmetric**. ⛔ Change nothing else in either header.
- **3 — the causality diagnostic fires on the passivity computation the directive itself mandates**, since
  the required Hermitian forms necessarily contain conjugated kernels. ⚠⚠ **This one reaches the
  falsification channel**: a spurious `FAIL` would let a genuine growing root be labelled an artifact.
  ⭐ Scope the diagnostic to the objects where an advanced kernel is actually a defect.
- **4 — Control F's criterion names an artifact no engine can read** ("revision 7's mechanics"). ⭐ Restate
  it **self-containedly**. ⛔ A task whose only honest outcome is a refusal or an invented comparison is
  worse than no task.
- **6 — "an independent invariant built from `∇·u`" is singular where the basis admits several.** ⭐ Fix the
  phrasing so an engine cannot carry one and silently drop the others. ⛔ Do not enumerate the basis for
  them — §3 already requires them to construct it.

- **Q3 gap — the energy balance is checked against an independent standard only in the sub-case where the
  interface does nothing.** ⭐ The reviewer shows the restriction is **unnecessary**: isolating **by order**
  in the reciprocal coefficient, rather than setting it to zero, loses nothing. ⭐ Implement a covering
  check on that basis. ⛔ It must classify no root and must say nothing about the sign of any imaginary
  part. ⭐ State in one line what wrong derivation it catches.

## ⛔⛔ CONSTRAINTS

1. ⛔ **MINIMAL AND SURGICAL.**
2. ⛔ **NEVER SUPPLY AN EXPECTED RESULT OR ITS REASON** — see the contamination warning above.
3. ⛔⛔ **KEEP THE FALSIFICATION CHANNEL OPEN.** A **growing** root remains first-class. Any diagnostic that
   could reject one must be explicitly scope-bounded to a sub-case where growth would necessarily be a
   derivation error.
4. ⛔ **EVERY CHECK MUST BE ABLE TO FAIL**, with a one-line statement of what wrong derivation it catches.
5. ⛔ **DO NOT RE-OPEN** what independent legs cleared: the scope boundary; `B1`'s mass balance; the `A/B/C`
   split; §1b's branch prescription **including its upper-half-plane extension** (verified by computation);
   §3b's virtual-displacement rule; the **supplied affinity** `𝒜` (independently derived and confirmed by
   two separate parties this round); `B8` controls B/C/D.
6. ⛔ **TWO ENGINES MUST NOT BE ABLE TO DIVERGE.**
7. ⛔ **No new scope** beyond what a finding or the decision above requires.

## Output

The edited files, plus a report **under 60 lines**: each finding, what you changed, the `sha256`, and for
every check added/kept/modified the one-line statement of what wrong derivation it catches.
⭐ **State explicitly that you did not propagate the reviewer's derived answer**, and name where you were
tempted to.

⛔ Then stop and exit.
