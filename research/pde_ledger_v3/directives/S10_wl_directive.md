# S10 — build the Mathematica engine (engine 1, BLIND)

**Write:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`

**Read first, in full, and treat as binding:**
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`

⛔ **Do not commit.** ⛔ **Do not modify any other file.**

---

## ⛔⛔ WHAT THIS ENGINE MAY NOT READ

This is **engine 1** and it is **blind by construction**. It must not read, import, grep, or otherwise
consult:

- `research/pde_ledger_v3/reduction/` — ⛔ **the registry, in any form.** This engine derives dimensions
  from the action alone. The comparison against the registry is engine 2's job, and it is only a real
  comparison if this engine never saw it.
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the sibling engine.
- `research/pde_ledger_v3/steps/` — ⛔ **all of it**, every file, without exception. ⚠ Some filenames in
  that directory state results; ⛔ **do not read the directory listing either.**
- `research/pde_ledger_v3/paper/` — the TeX cards and the built PDF.

⭐ **Two engines exist so that they can DISAGREE.** A disagreement is the single most valuable output this
build can produce. ⛔ It is destroyed the moment this engine is written against the other one's answer.

⚠ **You are overwriting a file that already exists.** ⛔ **Do not read the existing
`S10_brane_mode_spectrum_mathematica_audit.wl` and do not adapt it** — it carries the defects this
rebuild exists to remove. ⭐ Write a new file from `S10_SHARED_PHYSICS.md`.

## The mandate

Implement **Q1 through Q8** of the shared physics for **every package** in its §7 table, obeying §4 (the
structural rule), §5 (three clauses, four corollaries, no verdict), and §8 (tag grammar, `<ENGINE>` = `WL`).

## Mathematica-specific requirements

- ⭐ Emit payloads with `ToString[expr, InputForm]` so every value is re-parseable.
- ⭐ Carry the §3 assumption set — **as one joint `And[...]`, exactly as written there** — into every
  `FullSimplify` / `Simplify` / `Refine` / `Reduce` call. ⚠ A rank or a sign computed without it is a
  different quantity, and engine 2 is required to carry the identical joint set. ⭐ Emit the joint
  assumption expression itself as a tag.
- ⚠ `Solve[... , omegaSquared]` may return `ConditionalExpression`. ⛔ Do not drop the condition — emit it
  as its own tag, with the `_LOCAL_` infix per shared physics §8. If a payload would print as a bare
  `ConditionalExpression[value, condition]`, emit **both** parts separately.
- ⭐ Every locus solve uses `Reduce[… , {k1, …, kD}, Reals]` or `Solve[…, {k1,…,kD}, Reals]`. ⛔ An
  unrestricted solve returns complex wavevectors, which §3 forbids.
- ⚠ `MatrixRank` returns the **generic** rank — this is the subject of Q8, so compute the rank-drop locus
  explicitly rather than assuming genericity.
- ⭐ For Q4 **N3**, build the stacked matrix explicitly, e.g. `Join[Mr, {kvec}]`, and take its
  `MatrixRank`. ⛔ Do not attempt to infer `nu_T` from the `NullSpace` basis.
- ⭐ For Q6, walk the expression tree. `ρ_br` and `μ_R` are **unknowns to be solved for**, not looked up.
- ⛔ **No `Exit[1]` on a physics outcome.** Exit non-zero **only** if the kernel cannot complete the
  computation.
- ⛔ **No `VERDICT`, `PASS`, `FAIL`, or `checks` list anywhere in the file.**
- ⭐ Set `$HistoryLength = 0` and prefer `Together`/`Cancel` over repeated `FullSimplify` on large `D` if
  runtime becomes a problem — then record which simplifier was used, per §7's runtime rule.

## Verify before you finish

Run the script. It must complete inside **10 minutes** and exit **0**. ⭐ Paste nothing from its output into
your report.

⇒ Then answer §10 of the shared physics, **under 25 lines.**
