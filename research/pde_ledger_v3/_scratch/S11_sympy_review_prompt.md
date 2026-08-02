# REVIEW — the S11 SymPy audit and its registry insertion

**Read-only review.** ⛔ Do not edit or create any file inside `/var/projects/toy_physics`. Use `/tmp` if
you need scratch space. Report findings as text.

## What you are reviewing

- `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
- its changes to `research/pde_ledger_v3/reduction/quantities.yaml`, `relations.yaml`,
  and `acceptance_check.py` (see them with `git diff` — none of this is committed)

Its specification is `research/pde_ledger_v3/_scratch/S11_sympy_directive.md`, which you **may and
should** read.

## ⛔ Do not read

- `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` — this is the
  **other engine** for the same physics. The entire value of the pair is that they were written
  independently. ⛔ If you check one against the other you destroy exactly the property under test.
  **Derive everything yourself.**
- `research/pde_ledger_v3/V3_STEP_PLAN.md`, anything under `research/pde_ledger_v3/steps/`
- `research/pde_ledger_v2/scripts/ledger_stage003_*`,
  `research/pde_ledger_v2/mathematica/ledger_stage003_*`

---

## What to check

**S1 · ⭐⭐ The acceptance fixture — highest stakes item, do this one properly.**
`reduction/acceptance_check.py` carries `EXPECTED_MEDIUM_PAYLOAD`, a **control fixture** whose whole
purpose is to be an *independent* statement of what the registry's reduction should produce. Its own
comment says *"On any registry change, RECOMPUTE and independently re-derive; never copy forward."*

⛔⛔ **A fixture back-filled from the code's own output is the code agreeing with itself, and it would
silently destroy the only control in the harness.**

⇒ **Derive the expected payload yourself**, from the registry contents: how many independent continuous
quantities are there, how many does each admitted relation eliminate, and what does each of the three
mutation cases do to that count? Then compare with what is in the file. Report your derivation and
whether the committed numbers are right. ⚠ *"The gate passes"* is **not** an answer to this question.

**S2 · Is the physics right?** Derive independently: the dynamical matrix for a plane wave, the distinct
roots, the nullity and orientation of each root's kernel, the dimensional exponent vectors, and the
bulk-matching result. Compare with what the script computes. ⛔ Do not assume the script is right.

**S3 · ⭐ The invariant census.** The script reports counts of independent rotationally-invariant
quadratic forms in `∂_i u_j`, separately for `SO(D)` and `O(D)`, for `D = 2,3,4,5`. **Derive both counts
yourself from scratch**, then say whether the script's numbers are right and whether its method is a
correct way to obtain them. ⚠ If your counts differ from the script's for any `D`, say which and why —
that is precisely the kind of thing this review exists to surface.

**S4 · ⛔ TAUTOLOGICAL CHECKS — two specific shapes are known to occur in this project; hunt them.**
- **A closure check that cannot fail:** a quantity *defined* as `a − b` and then "verified" by adding
  `b` back. This shape has appeared before in dimensional blocks — it validates nothing about the
  formulas, only that subtraction is invertible. ⭐ Test it by **ablation**: change a premise to a
  physically absurd value and see whether the script still reports `PASS`.
- **A residual that recomputes its own selection predicate:** a mode is *selected* by testing some
  condition, and then a "check" asserts that same condition holds. It cannot fail by construction.
- Also: any conclusion emitted as an unconditional literal token rather than derived from a computed
  value.

**S5 · Are the registry dimensions DERIVED or asserted?** The directive requires the declared
`dimension.exponents` for the new quantities to be the values the script derived, and requires a
cross-check against the registry that **fails** on disagreement. Verify that cross-check is real —
⭐ **ablate it**: perturb a declared exponent and confirm the script fails rather than passes.

**S6 · Does `VERDICT: PASS` depend on anything substantive?** Could the script emit `PASS` with the
physics wrong? The directive required `PASS` to be the conjunction of an enumerated, printed list of
assertions — verify that is what it does.

**S7 · Registry insertion hygiene.** Is the new relation well-formed and complete against the shape of
the existing ones (guards, assumptions, provenance, designated output)? ⚠ **Naming:** the new
compression modulus must not be aliased to, or given the provenance loci of, `K_br` (a rejected
elastic branch's bulk modulus) or `B_eff` (`= ρ_B0²/χ_c` elsewhere in the corpus) — they are different
objects. Check it did not conflate them.

**S8 · Hardcoded values, symbol assumption errors, spec deviations.** Any target or expected result
baked in; any `positive`/`real` assumption quietly doing work it should not; any place the script
departs from its directive. ⚠ A deviation that makes the answer *more* correct is worth flagging as a
note, not as a defect — say which it is.

---

## ⭐⭐ THE FILTER — apply it before writing anything down

Report a finding **only if it catches a way the PHYSICS or the SCRIPT'S CORRECTNESS could be wrong.**

⛔ Do **not** report: style, naming, formatting, comments, performance, "it would break on a different
input", missing error handling, or process suggestions.

⚠ This project has repeatedly lost days to reviews producing many findings of which nearly none were
about the physics. **A short report with one real finding beats a long one.** If it is sound, say so
plainly and stop — "no blocking findings" is a complete and welcome result.

## Output format

```
VERDICT: SOUND | BLOCKING FINDINGS
BLOCKING          (wrong physics, or an unearned PASS — list, or "none")
ACCEPTANCE FIXTURE (S1: your independent derivation, and whether the committed numbers are right)
INVARIANT COUNT   (S3: your independently derived SO(D) and O(D) counts, and whether the script's are right)
TAUTOLOGICAL      (S4/S5/S6 — with ablation evidence where you ran one — list, or "none")
NOT-DERIVED       (values asserted that should be computed — list, or "none")
REGISTRY          (S7 — list, or "none")
WHAT IS SOUND     (brief — what you verified and found correct)
```
