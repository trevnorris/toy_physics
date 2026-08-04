# Independent review — the S10 Mathematica engine (SCRIPT review)

⭐⭐ **THIS IS A REVIEW OF A SCRIPT, ⛔ NOT A PHYSICS INVESTIGATION.** The question is **"does this script
compute what its specification asks, and does every printed value actually depend on a computation?"**
⛔ It is **not** "what is the right answer?" — the script's output answers that, and it is read **after**
this review, ⛔ not during it.

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`

Committed at `3bf61a6b`. It is **engine 1 of two**, written **blind**: barred from the registry, from the
sibling SymPy engine (which does not exist yet), and from `steps/`.

## Its specification — ⭐ read this FIRST, in full, before opening the artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
and `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_wl_directive.md`

⭐ Form your own view of what the specification **demands** before you see how it was implemented.

## ⛔⛔ DO NOT READ

- `/var/projects/toy_physics/research/pde_ledger_v3/steps/` — **the whole directory**, every file. ⛔ Do
  not list it, grep it, or reach it via `git show` / `git log`.
- `/var/projects/toy_physics/research/pde_ledger_v3/paper/` — cards and PDF.
- ⛔ **Any earlier version of this engine**, including via `git log`/`git show` on its path. The previous
  version is the thing being replaced and it carries the defects at issue.

⭐ You **may** read the `reduction/` registry and `REBUILD_HANDOFF.md`.

---

## What to check

### ⭐⭐⭐ 1 · CONFORMANCE — is each requirement implemented, and implemented as specified?

Walk **Q1 through Q8** of the shared physics against the code. For **each** requirement report
`IMPLEMENTED` / `PARTIAL` / `ABSENT` **with the line number**, and where partial, exactly what is missing.

⭐ Pay particular attention to these, each of which exists because a weaker version of it has failed here:

- **Q4 N3** — is the transverse count the **rank of `M_r` stacked on the row `kᵀ`**? ⛔ Or is it inferred
  from the vectors `NullSpace` returned? The spec forbids the latter because `NullSpace` returns an
  arbitrary basis.
- **Q4 N7** — does it difference the **null-space routine's** basis count against the **rank routine's**
  implied count (two different algorithms)? ⛔ Or has it collapsed into `rank + (D − rank) − D`, which is
  identically zero for any rank whatsoever?
- **Q5** — is dispersion tested by **scaling `k → λk` and emitting the ratio**? ⛔ Or by substituting a
  symbol for `Σ k_m²`, which silently no-ops?
- **Q6** — are dimensions obtained by **walking the whole expression tree** and **counting field factors**?
  ⛔ Or by reading the exponents of coefficient symbols (drops explicit wavevector factors) or by summing
  `Derivative` nodes (gives no contribution for a bare undifferentiated field)?
- **Q8b** — are Q3 and Q4 **re-run on each allowed stratum**, or is the locus merely reported?
- **§7** — does **every** package re-enter at the **ACTION**, recomputing everything downstream by the same
  code path? ⛔ Or does any control re-enter at a **result**?

### ⭐⭐⭐ 2 · DEPENDENCY — ABLATION IS MANDATORY, ⛔ NOT OPTIONAL

⛔⛔ **Copy the file to `/tmp/` and ablate the copy. ⛔ Never modify the working tree.**

⭐ **A FORM ablation, ⛔ not a coefficient rescale.** Change the **structure** of a load-bearing object at
the **action**: flip a sign, change which indices are contracted, collapse two independent symbols into
one. Re-run and report the **literal diff** in the emitted tags.

⭐ Do **at least four**, each hitting a different part of the pipeline, and for each report **which tags
moved and which did not**.

⚠ **A coefficient rescale tests arithmetic only** — scaling never leaves the family. ⭐ Only a form change
tests whether the physics is being computed.

⭐⭐ **The defect this is hunting:** a tag whose payload is a **hand-typed CAS object** — genuine algebra,
authored rather than derived, with **no data dependency on the computation**. Delete `Det`/`Solve` and it
does not move. ⇒ ⭐ **For every tag that did NOT move under any ablation, say whether it is legitimately
invariant or never computed.** ⚠ Those are indistinguishable from the output alone and that is the point.

### 3 · ⭐ Two routes and one-sided corruption

Q2 emits `M_A`, `M_B` and their difference, and the spec states this tests **coding consistency only**.
⭐ Verify the script says so and does **not** present it as verifying physics. ⭐ Then **corrupt one route
only** and confirm the other does not move; if breaking route A also moves route B's object, they were
never separately constructed.

### 4 · ⛔ Prohibited constructs

Report every instance, with line numbers:
- a `Print`/`emit` payload that is **prose stating a conclusion** rather than a CAS object;
- a **verdict**, `PASS`/`FAIL`, or summary judgement;
- `Exit[1]` (or any nonzero exit) reached on a **physics outcome** rather than an operational failure;
- a **tag name** that states a value, count, sign, ratio, or the shape of an answer;
- **emission made conditional on a payload's value** — a tag that appears only when its value is
  interesting, which makes *legitimately invariant* indistinguishable from *never computed*;
- an assertion or guard placed **before** the emit of the value it guards.

### 5 · ⚠ Runtime and completeness

Run the script yourself and report wall-clock and exit code. ⛔ It must finish under **10 minutes**.
⭐ Confirm the `(package, D)` pairs it reports running are the ones §7 requires, and that nothing was
silently skipped.

---

## ⛔⛔ WHAT NOT TO DO

⛔ **Do not report what any physical quantity came out to be.** ⛔ Do not derive the spectrum, the mode
counts, or the dimensions independently and tell me the values. ⭐ If you believe the script computes
something **other than what the specification asks**, report **that discrepancy** — the method, not the
number.

⛔ Do not propose or make edits. Report only.

## Physics filter

⛔ Report a finding only if it catches a way the **physics could be computed wrongly** or a way a **printed
tag could fail to mean what it appears to mean**. ⛔ Not style, not "wrong on a different input".

## Report format — ⛔ under 40 lines

1. **Q1–Q8 conformance table**: requirement → `IMPLEMENTED`/`PARTIAL`/`ABSENT` → line → note.
2. **Ablation table**: what you changed → tags that moved → tags that did not → your read.
3. **Findings**, most severe first, each with a line number and why it matters.
4. **Tags that moved under nothing**, and for each: legitimately invariant, or never computed?
5. ⭐ Anything wrong this prompt did not ask about. ⭐ This is wanted.
