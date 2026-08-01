# REVIEW — a seven-part fold, uncommitted in the working tree. Read-only, one pass.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.
**Run `git log --oneline -3` and `git status` first.** ⛔ Do not trust a hash written in any document,
including this one. The work under review is **uncommitted** — `git diff` shows all of it.

## What was folded

1. **`CHARTER.md` §1.3** — a new statement of what would count as FAILURE. The charter previously said what
   v3 excludes and yields, but never what failure is.
2. **`V3_STEP_PLAN.md` S20a** — a **correction**. The block previously asserted `λγ = 1` is an
   experimental fact the model must reproduce and that S20a is a place the model can fail. It now says
   `λγ` is a calibration and that a calibrated value cannot falsify.
3. **`V3_STEP_PLAN.md` S11** — records that the stray longitudinal mode is load-bearing rather than a
   defect, despite its token being named `FAIL_CAUCHY_STRAY_LONGITUDINAL`.
4. **`V3_STEP_PLAN.md` S0.5** — a note that `R⟨n⟩` register edges and `R1`-the-rung are different
   namespaces.
5. **`DEFECT_REGISTER.md` `A-CAUGHT`** — a pin-shaped identification that was flagged, checked and
   refuted before being made.
6. **`DEFECT_REGISTER.md` `B2`** — a **rescope**. The row previously read as though light rested on a
   falsified foundation.
7. **`NEXT_SESSION.md`** — an eleven-item list of beliefs a fresh session gets wrong, each with a
   correction and a locus.

## What to check

⭐ These are questions, ⛔ not a list of known answers. **An empty BLOCKING list is a real verdict.**

1. ⭐⭐ **Is every claim true to its source?** Open every locus. This fold makes substantive physics claims
   — that one medium can carry both a longitudinal and a transverse mode consistently; that the stray
   longitudinal is what becomes charge; that `r_cone` and `λγ` are different objects; that `B2` rules out
   only a derivation route. ⛔ Verify each against the cited artifact, not against this summary.
2. **Is the S20a correction now right, or has it over-corrected?** It was too strong before. Check it has
   not become too weak — in particular whether anything of value was lost in the rewrite, and whether the
   observational constraint is still stated accurately.
3. **Does `B2`'s rescope soften a genuine negative result?** Closed negative results are the most valuable
   artifacts in this project and ⛔ must not be diluted. Check that the no-go itself still bites and that
   only its *scope* changed.
4. **Are the eleven items in `NEXT_SESSION.md` each accurate?** They are asserted as measured fact. Any
   that is wrong will mislead every future session, so they carry more weight than ordinary prose.
5. **Is §1.3's falsification standard consistent with the rest of the charter**, especially §1.1, §1.2 and
   §4's acceptance criteria? Does it contradict anything already stated?
6. ⭐ **Append residue.** For each of the two CORRECTIONS (S20a, `B2`), check that the pre-correction claim
   does not survive anywhere — in this file, in another file, in a heading, a table row, a summary line, or
   a quotation. This project's recurring defect is a correction added while the original assertion stays.
7. **Loci.** Citations here use two forms: name-based tags (`` `<path>.md#<tag>` ``) for files this project
   owns, and line numbers for everything else. Check both kinds resolve. `scripts/check_citation_tags.py`
   checks the tag form; ⛔ it does not check line-number citations, so check those by opening them.
8. **Anything the fold missed** that its own findings imply.

## Operating constraints

⛔⛔ **BUDGET — read this first.** This fold is large and its sources are scattered. **Two prior runs of
this exact review investigated for a long time and then died having emitted NOTHING.** ⇒ Work the
"What to check" list **in the order given** — it is ordered by value — and ⭐ **emit the report as soon as
you have covered items 1, 2, 3 and 6, even if the rest is untouched.** State plainly which items you did
not reach. ⛔ A partial report that ships is worth far more than a thorough one that never does.

- **READ ONLY.** ⛔ Do not modify, stage, or commit. ⛔ One pass, no clarifying questions.
- You may run Python/SymPy. Write scratch **outside** the repository.
- ⭐ Open the cited file rather than trusting a document's summary of it.
- ⛔ For any claim of the form "there is no X", read the WHOLE artifact first. `wc -l` before any universal
  negative.

## Output

```
# Review — <your name/model>

## VERDICT
<CLEAN / STILL BLOCKING / REGRESSED> + 2-3 sentences

## BLOCKING
<numbered; each: claim, file:line, why wrong, what to change. EMPTY IS A REAL VERDICT.>

## CLAIMS THAT DO NOT MATCH THEIR SOURCE
<table: claim | cited locus | what the source actually says>

## APPEND RESIDUE FOUND

## NON-BLOCKING

## LOCI THAT DO NOT RESOLVE
```

## Standards

A matching number is not evidence. Dimensional agreement is not physical agreement. Falsification is
welcome and first-class. Absence of a denial is not evidence. ⛔ Do not be agreeable — ⛔ but do not invent
blockers either. If it is ready, say it is ready.
