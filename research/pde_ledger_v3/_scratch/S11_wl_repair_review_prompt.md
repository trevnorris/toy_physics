# REVIEW — a repair that removes two unearned claims from the S11 Mathematica audit

**Read-only.** ⛔ Do not edit any file under `/var/projects/toy_physics`. Extract to `/tmp` if you need
to run or ablate it. Report findings as text.

Inspect the change with `git diff -- research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
(nothing here is committed). The repair spec is `_scratch/S11_wl_repair_directive.md`, which you may read.

## ⛔ Do not read

- `research/pde_ledger_v3/scripts/` — an **independent second engine** for this same physics lives there.
- `research/pde_ledger_v3/reduction/`, `V3_STEP_PLAN.md`, anything under `steps/`.

⚠ You must judge this script **on its own terms**, from what it computes. That constraint is not
bureaucratic — see P2.

---

## The two defects the repair was supposed to fix

**D1** — `WL_S11_TRANSVERSE_MATCHING` emitted a fixed conclusion token *unconditionally*: a literal
chosen when the code was written, which no computed value could change. Its accompanying residual was
also circular — a set of vectors was **selected** by a predicate, then "checked" against that same
predicate.

**D2** — the four `dimensionalClosure` checks were algebraic rearrangements of their own definitions
(a quantity defined as `a − b`, "verified" by adding `b` back). Demonstrated: setting the displacement's
assumed dimension to an absurd value produced nonsense dimension vectors and the script **still printed
`WL_S11_VERDICT: PASS`, exit 0**.

---

## What to check

**P1 · ⭐⭐ Is D2 actually able to fail now? ABLATE IT.** Extract to `/tmp`, perturb the premise the
dimensional derivation rests on to a physically absurd value, run it, and report the exact output.
**It must now FAIL with a non-zero exit.** ⚠ If it still passes, the repair did not land — blocking.
Then try a *second*, different corruption (e.g. an exponent in one derived vector) and report whether
that is caught too, or whether the guard only catches the one case it was built against.

⛔ **Also check what the guard is made of.** The directive forbade adding a table of expected dimension
values, because a hardcoded target is worse than no check. If the new guard works by comparing against
typed-in expected exponents, **that is a finding**. A legitimate guard tests an internal consistency
requirement the physics imposes — something a wrong premise would violate.

**P2 · ⭐⭐ Is D1's replacement EARNED BY THIS SCRIPT?** The repair had two honest options: derive the
conclusion from something computed, or stop concluding and instead state the **scope** of what the
calculation can settle.

⇒ Whichever it took, ask: **does every statement it now prints follow from what this script actually
computes and was actually given?** ⛔ Flag any assertion about the physics that the script has no way to
establish from its own inputs — a claim can be *true* and still be unearned here, and an unearned true
claim is exactly what D1 was.

⚠ Pay attention to whether the circular residual was fixed or merely re-worded: is anything still
"verified" using the predicate that selected it?

**P3 · ⭐ Did any COMPUTED PHYSICS VALUE change?** The repair was to touch **claims, not computations**.
Every other `WL_S11_*` tag — invariant counts, dynamical matrix, all four spectra with nullities and
orientations, cross-sector results, degeneracy condition, the four dimension tags, `KW_SQUARED`,
`BOUND_CONDITION`, both controls — must print exactly what it printed before. Diff the tag lines and
report. ⛔ **Any changed physics value is a blocking finding** — it would mean something was silently
absorbed rather than surfaced.

**P4 · Did the repair introduce a new tautology?** This project has caught the same shape twice:
a value verified using the definition or predicate that produced it. Check the new code for it.

---

## ⭐⭐ THE FILTER

Report a finding **only if it catches a way the PHYSICS or the REPAIR'S EFFECTIVENESS could be wrong.**

⛔ Not: style, naming, formatting, comments, performance, "wrong on a different input", error handling.

⚠ A short report with one real finding beats a long one. If the repair landed, say so plainly and show
the ablation output.

## Output format

```
VERDICT: REPAIR LANDED | REPAIR DID NOT LAND | BLOCKING FINDINGS
D2 ABLATION    (P1: exact perturbation, exact output, PASS or FAIL)
D2 SECOND      (P1: a different corruption — caught or missed?)
D2 GUARD KIND  (P1: internal consistency requirement, or hardcoded expected values?)
D1 EARNED      (P2: does everything it now prints follow from what it computes? circularity gone?)
UNCHANGED      (P3: confirmed, or list every tag whose value moved)
NEW TAUTOLOGY  (P4 — list, or "none")
```
