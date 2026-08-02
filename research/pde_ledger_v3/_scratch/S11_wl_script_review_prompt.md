# REVIEW — the S11 blind Mathematica audit, before anything is believed

**Read-only review.** Report findings as text.

## ⛔⛔ FIRST — a hazard specific to this review

The script under review is **deliberately absent from the working tree**. A second builder is right now
writing an independent SymPy version of the same physics, and it must not see this file.

⇒ ⛔⛔ **DO NOT restore, copy, check out, or write this file (or any part of it) anywhere inside
`/var/projects/toy_physics`.** Doing so would silently destroy the independence of a build in flight.
If you need it on disk to run it, extract it to `/tmp` and work there.

Read it with:
```
git show 911d0af8:research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
```

You may run it (`math -script`) from a `/tmp` copy. ⚠ Mathematica here has a **2-seat** licence — never
more than 2 concurrent `math -script` processes. A licence-contention run exits **40 with no output**;
⛔ do not score that as a failure of the script.

## ⛔ Do not read

- `research/pde_ledger_v3/V3_STEP_PLAN.md`, anything under `research/pde_ledger_v3/steps/`
- `research/pde_ledger_v2/scripts/ledger_stage003_*`,
  `research/pde_ledger_v2/mathematica/ledger_stage003_*`
- `research/pde_ledger_v3/_scratch/S11_directive_review_prompt.md`

These contain the step's narrative or a related prior calculation. Reading them would make you check the
script against a remembered answer instead of deriving independently — which is the point.

⭐ You **may** read `research/pde_ledger_v3/_scratch/S11_wl_directive.md` — the specification the script
was built from. Part of your job is checking the script against it.

---

## What the script is for

It is one of **two independent engines** for a derivation step. Its value is that it could **disagree**
with the other. It imports nothing and reads no registry, by design. Its own `VERDICT` line is
explicitly **not** a certificate of the physics — it only claims the script's internal checks did not
contradict each other.

---

## What to check

**R1 · Is the nullity genuinely computed?** The spec requires the **dimension of the kernel** of the
dynamical matrix at each root — ⛔ not a multiplicity read off an exponent in a factorised determinant.
Verify what the code actually does. A multiplicity that happens to equal the nullity here would still be
the wrong computation and would fail in a case with a non-trivial Jordan structure.

**R2 · ⭐ Independently verify the invariant census.** The script reports a count of independent
rotationally-invariant quadratic forms in the `D×D` matrix `∂_i u_j`, for `D = 2,3,4,5`. **Derive that
count yourself, from scratch, before looking closely at how the script gets it.** Then say whether the
script's number is right and whether its construction is a correct way to obtain it.
⚠ If your count and the script's differ for some `D`, that is exactly the kind of thing this review
exists to surface — state which `D`, and what the correct count is and why.

**R3 · Are the dimensional results derived or asserted?** Which `[L,T,M]` vectors are computed from a
stated premise, and which are typed in as literals? A premise the *directive* did not supply but the
script assumed is a finding — name it and say whether the assumed value is correct.

**R4 · Is the bulk matching computed?** Or is the classification (bound vs radiating vs threshold)
hard-coded around a pre-decided inequality?

**R5 · Tautological or vacuous checks.** Any assertion that cannot fail; any comparison of a value with
itself; any `PASS` that does not depend on a substantive computation. ⭐ Specifically: could the script
emit its `PASS` verdict with the physics wrong?

**R6 · Hardcoded expected values.** Any target, reference number, or expected result baked into the
script. It was built blind and must contain none.

**R7 · Fidelity to the spec.** Read `S11_wl_directive.md` and check the script implements what it says —
particularly the two **controls**. One is meant to change the *shape* of the answer and one only its
*arithmetic*. Verify each does what it claims: a control that doesn't actually alter the operator's
structure isn't testing the physics.

**R8 · Symbol and assumption errors.** Sign conventions, factors, an assumption (`positive`, `real`)
that quietly does work it shouldn't, a substitution valid only in a special case.

---

## ⭐⭐ THE FILTER — apply it before writing anything down

Report a finding **only if it catches a way the PHYSICS or the SCRIPT'S CORRECTNESS could be wrong.**

⛔ Do **not** report: style, naming, formatting, comments, performance, "it would break on a different
input", missing error handling, or process suggestions.

⚠ This project has repeatedly lost days to reviews producing many findings of which nearly none were
about the physics. **A short report with one real finding beats a long one.** If the script is sound,
say so plainly and stop — "no blocking findings" is a complete and welcome result.

## Output format

```
VERDICT: SOUND | BLOCKING FINDINGS
BLOCKING     (the script computes something wrong, or its PASS is unearned — list, or "none")
INVARIANT COUNT (R2: your independently derived count, and whether the script's is right)
NOT-COMPUTED (values asserted/hardcoded that should be derived — list, or "none")
TAUTOLOGICAL (checks that cannot fail — list, or "none")
SPEC DEVIATIONS (R7 — list, or "none")
WHAT IS SOUND (brief — what you verified and found correct)
```
