# The cross-engine comparator — fix round 2. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: `scripts/S10_cross_engine_comparator.py` · its regenerated `scripts/out/*.out` · a new directive
under `directives/`. ⛔ Out of scope: both engines, every `.wl`, every step record, every `.tex`.

⭐ **Two legs reviewed the last round. ⭐⭐ No false agreement, no misclassification in either direction, and
the nullspace residual passes BOTH halves** — it fires on wrong-subspace, partial rescale, rank loss and
dropped vectors, and stays silent under whole-vector rescale, permutation and a general invertible change
of basis.

⚠ **Everything below is LATENT on this input** — none of it changes a verdict today. ⭐ It matters because
**this instrument is meant to be reused for S9 and S11**, where the inputs differ.

⛔⛔ **Do not state, anywhere in the script or your directive, what any count, tally or population size
comes out as.** ⚠ The previous decision list stated its measured counts and then forbade stating them; the
builder transcribed them forward, and a leg correctly reported that the build therefore is not evidence for
them. ⭐ **Emit them and let them be read.**

---

## F1 · A derivative must carry the field's dependence set.

The derivative normalisation reduces an applied function to its **name** and its (variable, order) pairs,
and **discards the argument list** on both sides. ⇒ two engines emitting the same field name over
**different dependence sets** would compare **equal on every derivative term**.

⚠ Dead on this input — no derivative-bearing row currently passes. ⛔ Not dead at S11.

**Decision: two derivatives compare equal only if the differentiated object is the same object**, including
what it depends on. ⚠ The two engines order their arguments differently, so ⛔ the dependence set is a
**set**, not a sequence — ⭐ do not make the argument *slots* significant, and do not make `∂x1` and `∂x2`
interchangeable.

## F2 · A payload symbol must not be captured as a CAS builtin.

The payload parse sympifies with no local namespace, so a symbol spelled `I`, `E`, `S`, `N`, `O`, `Q`,
`beta`, `gamma`, `zeta`, `oo`, `zoo`, `nan` or `Symbol` silently becomes the library's object instead of a
free symbol.

**Decision: a symbol in a transcript denotes what the engine that wrote it meant.** ⚠ One exception exists
in the current input and it is genuine — `pi` denotes π on both sides. ⭐ Keep that true.

## F3 · Delete the ordering normalisation that contradicts the stated rule.

The `dict` branch of the normaliser **sorts** its entries, against the directive's own "preserve sequence
order and shape." ⭐ Measured unreachable — the branch is never entered on this input.

**Decision: it goes.** ⛔ Do not replace it. ⚠ If some object genuinely needs unordered comparison, that is
a **decision to state**, not a default to leave lying in the path.

## F4 · The accounting property must hold unconditionally.

Duplicate names are removed from the shared set before the categories are computed, and the duplicate
counter counts **rows**, not **names** ⇒ if a shared name were ever duplicated, the categories would stop
summing to the shared-name count. ⭐ Zero duplicates today.

**Decision: every shared name carries exactly one verdict, and the arithmetic that says so cannot be made
false by an input.** ⭐ Demonstrate it by constructing the input that would break it.

## F5 · One disagreement kind is typed rather than computed.

A failing nullspace row is assigned its disagreement kind as a **literal**, unlike every other failing row,
whose kind is computed. ⚠ The residual and guard are computed; only the label is typed.

**Decision: the kind is read off the comparison, like every other row.**

---

## Record, ⛔ do not fix

⭐ The membership residual and the span residual do **different** jobs and both are load-bearing — a
mutation that leaves the kernel moves membership, a rank loss moves only the span. ⭐ **Say so**, so a
future reader does not "simplify" one of them away.

## Constraints

- ⛔ Do not modify either engine or re-run them; the transcripts are committed inputs.
- ⛔ Do not start `wolframscript`.
- ⛔ No hardcoded name→name pair table, no config file, no runner.
- ⭐ Report the complete `.out` diff and account for every changed line.

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** For each item construct the weakest change that
should be rejected and show whether it is.
⛔⛔ **For F1 and F2, the defect is invisible on the current input** — ⭐ you must **construct** the input
that exhibits it and show the repaired instrument catches it, and that the unrepaired one does not.
⭐ **And show the nullspace residual still passes both halves** — fires on a wrong subspace, silent under a
general change of basis. ⚠ A round that quietly breaks that is the failure this question exists to find.

## Deliverables

The changed files, the literal stdout, the complete diff, and every ablation script with its stdout at
named absolute paths **outside the repository**.
