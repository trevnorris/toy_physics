# Integration TODO — stop the reduction checks being one-offs

**Status: QUEUED, 2026-08-07.** ⛔ Blocked on the two review legs for the W2 witness fix round.
⚠ Deliberately short. ⛔ Do not grow this into a plan document.

## Why — measured, not asserted

- ⛔ **There is no runner.** `ls reduction/*.sh reduction/Makefile reduction/run*` → *NO RUNNER SCRIPT*.
  Every check is a separate script someone must remember to invoke.
- ⛔ **The step records cite ONE check each.** `steps/S9_light_requires_shear.md:229` and
  `steps/S10_two_transverse_photons.md:305` both cite only `engine_output_checks.py`. Neither cites
  `dimensional_homogeneity_gate.py`, `derived_or_declared.py`, `able_to_fail.py`, or the new witness.
- ⛔ **The witness ran against S9/S10 committed outputs and the result is recorded nowhere** — it exists
  only in a scratchpad transcript.

⇒ ⭐ **The testing was redone. The recording was not.** Without the four items below, this work evaporates.

---

## The four items

### 1 · A runner — ⛔ a list with exit codes, ⛔ NOT a framework

One command, given a step, that executes every `reduction/` check against it and prints each status and
exit code. ⭐ It must make **"which checks does this step actually pass"** answerable in one place.

⛔ **Scope discipline:** no config format, no plugin system, no new declaration file. If it grows past
roughly a screen, it is the wrong shape.
⚠ It must **report** each check's status, ⛔ never summarise them into a single verdict — a check that was
not run and a check that passed must not print the same way.

### 2 · Fold the witness's measurement into S9's and S10's records — ⛔ AFTER the legs clear

`REBUILD_HANDOFF.md` already names this as per-step pattern item 3: *"the step record cites what it
measured."* ⭐ Both records must gain what the witness found, including:

- the status counts, at their real width (⭐ `AGREEMENT` is **3**, not 9 — the tautology is reclassified
  `ECHOED` and axis order is `UNDETERMINED` where the engine emits no labelled components);
- ⭐ **`c_γ` and `c_L` are compared for the first time**, via a declared multiplier, residual zero;
- ⛔ **`BRANCH_DIMENSION_MISMATCH` on `rho_br` and `mu_R` at S10** — see item 4.

⛔ **Do not record any of this until both legs report.** Putting an unreviewed instrument's number into a
closed step's record is how a wrong measurement becomes a citation.

### 3 · Grow `REBUILD_HANDOFF.md`'s per-step pattern

It currently lists four items every rebuilt step must satisfy. ⭐ The witness belongs in it — otherwise
S11 and everything after runs it by memory, which is the same as not running it.

⭐ Add the one-line requirement that makes it work on **future** steps rather than only retrofits:
> **State which registry quantities this step produces and which it consumes; emit each in the form the
> registry declares it; mark each emission derived or declared.**

⚠ That is exactly what S9/S10 did not do — which is why `c_gamma` was emitted only as a *square* and R4
has been red at both steps since they closed.

### 4 · Record the registry-`D` decision

⛔ Open: `quantities.yaml` cannot express a **D-dependent** dimension. Every `Q.brane.*` entry is silently
a D=3 specialisation, while S10 varies D = 2…5 deliberately.

⚠ **On the measurement so far nobody is wrong** — the engines emit `[-D,0,1]`, correct for a
D-dimensional brane; the registry declares `[-3,0,1]`, correct at D=3. ⭐ **The gap is that nothing says
so.** A Codex consultation on the schema options is in flight; the review legs are adjudicating
independently whether this is an engine, registry, or witness fault.

⇒ ⭐ Whatever is decided lands in `quantities.yaml` **and** in the step records that relied on it.
⛔ Do not let an instrument paper over it by widening a tolerance.

---

## Sequencing

| | item | precondition |
|---|---|---|
| 1 | runner | ⭐ none — but it invokes the witness, so land it **with** item 2 |
| 2 | S9/S10 records | ⛔ **both W2 legs reported** |
| 3 | handoff pattern | none |
| 4 | registry-`D` decision | ⛔ **user decision**, informed by the legs + the Codex consultation |

⭐ Items 1 and 3 are prose/plumbing ⇒ one fresh reviewer. ⚠ Item 2 edits **physics records** ⇒ two legs.
⛔ Item 4 is not mine to decide.
