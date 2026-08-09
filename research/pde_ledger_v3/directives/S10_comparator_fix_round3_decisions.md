# The cross-engine comparator — fix round 3. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: `scripts/S10_cross_engine_comparator.py` · its regenerated `scripts/out/*.out` · a new directive
under `directives/`. ⛔ Out of scope: both engines, every `.wl`, every step record, every `.tex`.

⭐ Two legs reviewed round 2. ⭐⭐ **All five repairs hold**, each shown catching an input the unrepaired
instrument passes, with the overreach controls holding — and **the round-1 nullspace property survives in
both halves**, verified on real payloads including symbolic invertible changes of basis.

⛔⛔ **One repair overreached and is a regression. ⚠ The fault is in my previous decision list, which named
an instance instead of the property.** ⭐ Read F1 with that in mind.

⛔ Do not state, anywhere in the script or your directive, what any count, tally or population size comes
out as. ⭐ Emit them and let them be read.

---

## F1 · ⛔⛔ A mathematical constant means itself. Only a NAME THAT WOULD SHADOW gets captured.

Round 2 was told to stop payload symbols being captured as library objects, and `pi` was named as the
genuine exception. ⇒ ⛔ **`pi` was protected and the other constants were captured.**

⚠ **Measured — and this is the shape of the damage:**
- Both engines write the same complex object; the repaired instrument now calls it a **content
  divergence**, and the printed residual reads `2*k1*(I - I)` — ⛔ **an auditor reads that as zero.**
- Against a genuinely real locus, the instrument proposes **the imaginary unit and a density ratio as one
  symbol under two spellings**, on exactly the object the reality filter turns on.

**Decision: a name that denotes a mathematical constant in both engines keeps denoting it; only a name that
would SHADOW a payload symbol is captured.** ⚠ The two populations are different in kind — a namespace or
a callable shadows; ⛔ a constant does not.

⛔⛔ **And the mechanism must not be a name list.** ⚠ Measured: the list covers a handful of names against
hundreds that behave the same way, and one uncovered name already reaches a shared row. ⭐ **A denylist means
the architecture is wrong** — decide it by what the object *is*, ⛔ not by whether someone remembered it.

## F2 · A proposed symbol pair must not be able to bind two objects the engines already distinguish.

The worklist pairs a py-only symbol with a wl-only symbol whenever substituting one for the other zeroes
the residual. ⚠ Measured on constructed physics-shaped input: **a repointed wavevector component, a
stiffness repointed to a density, and a replaced dependence set are all indistinguishable from a genuine
spelling difference**, and each is proposed as a rename candidate.

⭐ Nothing is silently accepted — the row still fails and the residual is printed. ⛔ But this worklist is
about to be *consumed* by a naming pass, and acting on a proposed `stiffness ↔ density` line is the
recorded failure mode where re-pointing one entry made every check in the repository pass.

**Decision: the worklist may only propose a pair the transcripts do not already refute.** ⭐ If both engines
emit both names elsewhere as distinct objects, they are not two spellings of one symbol, and the
transcripts already say so. ⛔ If you conclude no such refutation is available for some pair, **mark that
pair as unrefuted rather than dropping it, and say so** — ⛔ do not silently narrow the worklist.

## F3 · The content-divergence population is not one thing, and the step record will misread it.

⚠ Measured: the large majority are **shape or type mismatches** and a further group are **route word
tokens**; only a small remainder carry a genuine numeric or algebraic residual. ⇒ ⛔ a step record citing
the content-divergence total as "physics disagreements" would be wrong by more than an order of magnitude.

**Decision: emit the decomposition, so the number that matters can be read rather than inferred.**

---

## Record, ⛔ do not fix

⭐ **The membership residual and the span residual catch different mutations** — measured in both
directions: a rank loss moves only the span; a perturbed matrix against identical spans moves only
membership, and that case includes a repointed wavevector component **nothing else catches.** ⭐ Say so.
⭐ The unparsed rows are one CAS's parser refusing a Boolean containing a domain-membership assertion.
⭐ **State that a name that denotes a constant is a shared assumption between the two engines** — ⛔ it is
not something this instrument derives, and a future step whose engines disagree about one would be
unpoliced here.

## Constraints

- ⛔ Do not modify either engine or re-run them; the transcripts are committed inputs.
- ⛔ Do not start `wolframscript`. ⛔ No hardcoded name→name pair table, no config file, no runner.
- ⭐ Report the complete `.out` diff and account for every changed line.

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.**
⛔⛔ **F1 and F2 are invisible on the committed input** — ⭐ **construct** the inputs that exhibit them, show
the repaired instrument is right and the unrepaired one is not, and ⭐ **show a genuine constant still
compares equal across the two engines** rather than becoming a divergence.
⛔⛔ **Show the round-1 nullspace property still passes BOTH halves**, and that round 2's five repairs still
hold. ⚠ A round that quietly breaks an earlier one is the failure this exists to find, and it has now
happened once.

## Deliverables

The changed files, the literal stdout, the complete diff, and every ablation script with its stdout at
named absolute paths **outside the repository**.
