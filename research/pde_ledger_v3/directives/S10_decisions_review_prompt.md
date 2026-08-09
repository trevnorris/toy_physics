# Independent review — S10 export-chain decision list

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_export_chain_decisions.md`

It was written by the orchestrator and has **not** been built from yet. You are one of two independent
legs; the other has this identical prompt and you will not see its output. **Your job is to find what is
wrong with it before a builder spends a round discovering the same thing.**

⚠ **This is a rewrite.** An earlier version was reviewed and replaced. Its central defect was that it
**invented mechanisms for problems the S9 precedent had already solved**, and the invented mechanisms were
unsound. The rewrite's central move is *"do what S9 does."* ⇒ **That move is now the thing most worth
attacking**, see question 9.

## Context you need

- `research/pde_ledger_v3/S9_REWRITE_PLAN.md` — the closed plan the list must not exceed.
- `research/pde_ledger_v3/REBUILD_HANDOFF.md` — current state, and the D-key convention.
- `research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py` and its generated
  `scripts/S9_exports.py` — the worked precedent, closed at `83eace28`.
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the engine the list governs.
  It does not currently run.
- `CLAUDE.md` — the method the list must comply with.

## What to check

**Read the S10 engine and the S9 precedent first, and form your own view of what this change requires.
Only then judge whether the list asks for it.** Reading the list first anchors you to its framing, which
is the thing under test.

Then, specifically:

1. **Is any instruction self-contradictory?** In particular, does any scope fence forbid touching a file
   that another item in the same list requires changing? **This exact defect has shipped twice in the
   last two lists.** Check every fence against every item.

2. **Does any instruction name a *location* where something must be recorded without naming the
   *property* it must have?** *"Record it at the point of production"* is fully satisfied by typing a
   constant there. **This exact defect shipped once already** and cost a round plus two legs. Find every
   instruction with the same shape.

3. **Is anything asked for unbuildable against the actual engine?** Read the code. If an instruction
   assumes a structure S10 does not have, or assumes S9 exposes something it does not, say so with line
   numbers.

4. **Does any instruction leak an expected result?** CLAUDE.md rule 5: the spec says what to compute,
   never what anything equals, is expected, or was measured. A builder iterating to exit 0 converges on
   any target it can see. Quote anything that gives a value, a count, a sign or a dimension the run is
   supposed to produce. ⚠ Distinguish this from a *measured fact about the current artifact*, which is
   legitimate context.

5. **What decision does the list leave the builder to invent?** Anything a builder must choose that two
   reasonable people would choose differently is a missing decision, and the list is where it belongs.
   Rank these — which one, if guessed wrong, does the most damage.

6. **Does the list exceed the closed plan?** The plan's failure mode is documented as *adding* machinery,
   not omitting it. Flag anything that would create a guard for a guard, a check that polices another
   check, or a file the plan's list does not permit.

7. **Would following this list produce a check that cannot fail?** The measured pattern: two operands
   derived from one source, or a residual asserted zero that is zero by construction. The list warns
   against it — check whether its own instructions could still produce one.

8. **Is the list's governing item (D0, input-driven) actually enforceable as written, and is its
   acceptance demonstration able to fail?** D0 says every emitted value must be reachable only through
   the computation, and prescribes a perturbation matrix over the declared inputs including a form
   change. Judge it hard, against the real engine:
   - Are S10's inputs identifiable and independently perturbable, or are several separately fixed at the
     same value so that perturbing one alone crashes the run? **This exact obstacle was measured on S9**,
     where the prescribed single-structure control was not executable.
   - Could a builder satisfy the matrix while leaving a typed literal in place? If a literal and a live
     readout produce identical output at the values the run actually uses, the matrix does not separate
     them — say so, and say what would.
   - Does D0 create pressure to fabricate? It demands a matrix and then forbids manufacturing a
     perturbation or reclassifying an entry to make it come out clean. Check that the honest outcome is
     always available and clearly sanctioned — a check whose only clean outcome is invented is worse
     than no check.

9. **Is "replicate S9" actually executable against S10's code, mechanism by mechanism?** This is the
   rewrite's load-bearing claim and it is untested. Read both engines and check each mechanism the list
   points at — the output collection, the per-object production metadata, the standard-name table, the
   class assignment, the export writer, the round-trip check. For each one:
   - Does S10 have a place to put it, or does adopting it require restructuring code the list does not
     mention?
   - Does S10 have a property S9 lacked that **breaks** the mechanism rather than merely differing from
     it? The list claims S10's only genuine differences are several component counts per run, a `LOCAL_`
     convention, and more objects to name. **Test that claim** — if there is a fourth, it is the most
     valuable thing you can report.
   - Is any S9 mechanism the list tells S10 to copy one that a review already found defective on S9? Copying
     a known-bad mechanism into a second engine is worse than inventing a new one.

## What a good finding looks like

Quote the instruction, say what a builder would do with it, and say what would go wrong. A finding that
says "this could be clearer" is not a finding. A finding that says "a builder reading this would
reasonably do X, and X is wrong because Y" is.

⭐ **Rank your findings by whether they change what gets computed.** The list's job is to make a correct
artifact, not to be well written.

## Do not read

- Any file under `directives/` beginning `S9_export_dkey_` — those are the previous round's decision
  lists and review briefs. They state defects already found and would make this leg a check that two
  copies of an assumption agree.
- Any Codex or review transcript. None is in the repository.

## Operational constraints

- ⛔ **This is a document review. Do not modify the working tree.** Do not build anything.
- If you run code to establish buildability, copy to `/tmp` and run the copy, wrapped in `timeout 600`.
- ⛔ Do not start `wolframscript` — no Wolfram engine is in scope.
- Report the absolute path of anything you write.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could end up wrong.
Do not report style, tone or wording that changes nothing about what gets computed.
