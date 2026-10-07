# O2 sub-step 6: cross-engine comparator (build directive / pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-07. `AGENTS.md` and `CLAUDE.md` apply.

**Gate:** this is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2), and no
builder starts before then. The comparator it describes is an instrument written by Codex. It is reviewed by
two non-author build legs (fresh Claude + Grok) until clear, as O2 scoping §7 row 6 requires
(`directives/O2_steady_brane_balance_scoping.md:252`).

## Object

`research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py`, with synthetic tests in
`research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py`. It compares the two accepted O2 constructions
object by object and prints what it finds. It decides nothing.

## Inputs (read-only)

- The production transcripts, which are git-annex content; run `datalad get` on both first:
  - `research/pde_ledger_v3/scripts/out/O2_live_balance_sympy_audit.out`: SymPy, one line per tag, payloads in the
    engine's lossless constructor form (`LosslessReprPrinter`, audit line 60), not stock `srepr`;
  - `research/pde_ledger_v3/mathematica/out/O2_live_balance_mathematica_audit.out`: Wolfram, multi-line
    associations whose held and OPEN heads are defined in the engine source.
- The two accepted engine sources (`0f2e0af8`). These are the authority for what each emitted object, head and
  name means. Do not run either engine.
- The physics spec `directives/O2_SHARED_PHYSICS.md` (`4680e251`), §§3–9.
- Earlier comparators, as machinery to adapt: `scripts/S11c_b_cross_engine_comparator.py` and its tests
  (multi-line Wolfram reader, accounting, memory discipline). Reuse only what is consistent with items 1–8 and
  re-verify it against these streams. `scripts/S11b_cross_engine_comparator.py`'s `residual` cannot be reused as
  is: it joins by name, returns zero for equal text atoms, compares other atoms by stock `srepr`, has a per-leaf
  time alarm, and prints verdict tokens.

## What must be true

1. **Object-level join, declared and injective.** The two streams cut the same physics at different grain, so
   tags cannot be joined by name. The comparator carries an explicit join table in its source. Each row pairs a
   SymPy object (a tag, or a component or entry inside one) with a Wolfram object (a tag plus key path). Each row
   cites the construction line in each engine that produces its object. The table is injective: no object is on
   two rows. An object is unjoined only when the other stream has no component, entry or key path for the same
   spec §9 object. A different tag name, key or head makes a join row, not an unjoined object. Content that
   cannot be subtracted stays on its joined row as a coverage finding. Every emitted object is on a row or
   listed as unjoined with that reason. Nothing is skipped silently.

2. **Names are bindings, not evidence.** Any map between the two engines' symbol and function names (for example
   the live profiles, the OPEN operands and the coordinates) is an explicit, injective table. Each entry cites the
   spec object both names denote. A name match is never counted as agreement. Equal mapped names are a
   substitution, never by themselves a zero: any agreement comes from the content compared beneath them, such as
   arguments, operands and structure. A **repoint control** must show that pointing a mapped name at a different
   object moves the printed residual.

3. **Three-valued residuals, printed before any guard.** For each joined row, print the SymPy operand, the
   Wolfram operand and their residual. Use exact symbolic subtraction where feasible. Otherwise use numeric
   evaluation at several generic rational points, chosen so that no factor the objects contain vanishes there
   (`CLAUDE.md` E1; numeric-PIT practice). Either way, print the points and both operands' values. There is no
   per-leaf time budget.
   - A residual that cannot be formed (a native boolean, or structurally mismatched or unsupported content) is
     printed as that outcome. It is a coverage finding, distinct from a parse failure.
   - No applied function is collapsed to a bare symbol, and no live profile is replaced by a constant, to make
     two operands meet.

4. **OPEN content is compared, not discarded.** For OPEN actions and balance entries, print each engine's
   structure:
   - the role or head;
   - the OPEN operands it names, through the name table;
   - the live quantities it takes as arguments;
   - its orientation sign.

   Then print the differences. A difference in representation (for example, native face geometry as a general
   immersion in one engine and an OPEN `𝒥_map` action in the other) is printed as a difference. It is never
   reconciled inside the comparator.

5. **No target.** Nothing in the comparator or its report states what any residual or comparison on the
   measured streams is expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`VERDICT`/`STATUS` token. It exits 0
   whatever it finds, and nonzero only on operational failure (a missing or ungrammatical input). Interpreting
   the output belongs to the record (sub-step 7).

6. **Controls that fail when their defect goes undetected.** The tests use synthetic fixtures only and never
   load the measured streams. Each fixture is serialized in its engine's own format (a lossless SymPy line, a
   Wolfram association) and goes through the same reader, extraction, join and comparison path as the production
   run. Each control is a test that fails unless its mutation changes the printed output; the assertion names
   no value. Required mutations:
   - one-sided corruption of an operand;
   - a form change, not a rescaling;
   - a repoint of every name-table row (item 2);
   - a stripped argument of an applied function;
   - a live profile replaced by a constant;
   - a changed OPEN head, named operand, argument or orientation;
   - changed derivative or binder structure;
   - a nested sibling removed;
   - a counterpart moved under a different tag or key path (it must still join, not fall to unjoined).

   Two further tests: a duplicate join or name row is rejected, and a native boolean is rejected as an operand
   while a sibling algebraic leaf is still subtracted.

   Run against the real streams, the comparator also prints per-row accounting: joined, or unjoined with its
   reason. It also prints two leaf counts per object: one enumerated independently from each stream's parse, and
   the count actually compared. A row compared below its parsed count stays visible.

7. **Ablation transcripts are out of scope.** The comparator reads the two production transcripts only. Each
   harness's knife triples are evidence about its own engine, already reviewed. They are not cross-engine
   objects and are not joined.

8. **Running.** Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode
   (`--pool o2-comparator --memory-gib 6`). There are no time limits. Stop and report any kill or admission
   refusal. Never answer one by narrowing what is compared. Write output under your own scratch directory,
   never under `scripts/out/` or `mathematica/out/`.

## Exclusions

- **No oracle check here.** Each candidate oracle in O2 scoping §5 applies only within a declared constitutive,
  fluid or kinetic limit. O2's construction keeps that content OPEN (spec §1; scoping §6), so no oracle's domain
  matches without a premise O2 does not supply. Choosing such a limit is a premise decision for the user, not a
  comparator task.
- No new physics, closure, normalization or profile choice. No engine edit. No commit.

## Builder report (at most 40 lines)

- The join table and the name table, each row with its two construction-line citations.
- Per-row accounting and extracted-leaf counts, and every unjoined object with its reason.
- Which machinery you reused, and what you re-verified against these streams.
- The commands for, and outcomes of, the tests and the full guarded run.
- Runtime and peak memory.
- That no residual target was given.
