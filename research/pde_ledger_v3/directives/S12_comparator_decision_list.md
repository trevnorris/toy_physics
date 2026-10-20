# S12 cross-engine comparator (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-20. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). No builder starts
before then. The comparator it describes is an instrument written by Codex, reviewed by two non-author build legs
(fresh Claude + Grok) until clear.
- **Pass 1** (Codex + Grok) returned four findings, all confirmed. The list as reviewed is preserved at `63ae8c9c`,
  with `_measurements/S12_comparator_decision_list_review_p1_{disposition,lookups}.md`.
- **This revision** folds those findings into items 1, 3, 4 and 6 and the builder report. It is the one fold G2
  allows.

## Object
`scripts/S12_cross_engine_comparator.py`, with synthetic tests in `scripts/test_S12_cross_engine_comparator.py`. It
compares the two S12 engines' printed objects, object by object, and prints what it finds. It decides nothing.

## Inputs (read-only)
- The two transcripts accepted at `a3da43d9`, which are git-annex content (run `datalad get` on both if their content
  is absent):
  - `scripts/out/S12_order_source.out` (SymPy);
  - `mathematica/out/S12_order_source.out` (Wolfram).
- The two engine sources at the same commit, `scripts/S12_order_source.py` and `mathematica/S12_order_source.wl`. They
  are the authority for what each printed object, head and name means. Do not run either engine.
- The S12 decision list, `directives/S12_decision_list.md`, revision 4 (`b08673a2`): "What the build computes" and
  "How U1 enters".
- Earlier comparators, as machinery to adapt: `scripts/S9b_cross_engine_comparator.py` and
  `scripts/test_S9b_cross_engine_comparator.py`, accepted at `88bb1f41`. Reuse only what is consistent with items 1–7
  below, and re-verify it on these streams.

## What must be true

1. **Join by tag name; a name is a binding.** A tag name the two transcripts share is joined. A joined name is a
   binding, never agreement: agreement comes only from the content compared beneath it. Each joined name occurs once
   per stream, and a repeat is reported. Every tag is on a joined row, or is listed as unjoined with its engine and
   the reason. Nothing is skipped silently. No pairing is made between differently named tags.
   - A shared name can bind payloads that are not constituted the same way. Before a joined row's payloads are
     compared as values, the comparator establishes which components of the two payloads denote the same object,
     citing the line in each engine that builds them.
   - A row without such a correspondence is listed as not compared, with its reason.
   - A structural mismatch is printed as a comparison that could not be formed, for example different component
     counts or different kinds of object. No pairing is made by position alone.
   - Both operands stay printed under their original tags.

2. **Symbols and functions through an explicit table.** Any map between the engines' names (coordinates, faces,
   regional fields, the interior order, `G`, the rates, the dimension symbols) is an explicit, injective table in the
   comparator's source. Each entry cites the object both names denote and the line in each engine that introduces it.
   A mapped name is a substitution, never by itself a zero. A repoint control shows that pointing a mapped name at a
   different object moves the printed residual.

3. **Distributions compared part by part.** An object that carries steps or deltas in `w` is compared in parts: its
   regular part on each region, and each face delta's part. Each part is compared with its counterpart.
   - A delta's coefficient carries data about the distribution only on that delta's support. A delta's part is
     therefore its coefficient restricted to the delta's support (its value on the face), with the delta's
     normalization in `w`. The comparator derives it from the printed distribution itself.
   - A raw coefficient, unrestricted, may also be compared. If so, it is a separately labelled representation
     diagnostic, never the delta's part.
   - Content that falls in no such part (another singular atom, or a part the comparator cannot isolate) is printed as
     that outcome, as a coverage finding.

4. **Residuals printed before any guard.** For each compared value or part, print the SymPy operand, the Wolfram
   operand and their residual.
   - Use exact symbolic subtraction where feasible, and print the exact residual as computed. Where it is feasible,
     evaluation at points is added beside it, never in its place. Otherwise evaluate at several generic rational
     points (E1; numeric-PIT practice).
   - Under evaluation, every live function is replaced by test functions of several different forms, applied
     identically to both operands. That includes each regional field, the interior order, both faces and `G`. Print
     the test functions, the points and both operands' values.
   - For a function of several arguments, the test functions separate each argument's dependence: at the evaluation
     points, its derivatives in different arguments differ. That includes both offsets of `G`.
   - Points lie in each region and on each face, inside the domain the engines print as supplied, which includes the
     strict face order. A residual's dependence on that domain is printed: the residual as computed, and its value on
     the supplied domain, side by side.
   - The evaluation covers strata where a field's jump at a face vanishes, and where the interior order's value at a
     face vanishes, as well as the generic case. A stratum the comparator cannot evaluate is printed as unresolved
     coverage.
   - A dimension object is compared as a value. A label (a premise line, a scope line, a comparison note) is compared
     as text.
   - No applied function is collapsed to a bare symbol, and no live function is replaced by a constant, to make two
     operands meet. A difference in representation is printed as a difference, never reconciled inside the
     comparator.

5. **No target.** Nothing in the comparator or its report states what any residual on the measured streams is
   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
   and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output belongs to the
   S12 record.

6. **Controls that fail when their defect goes undetected.** The tests use synthetic fixtures only and never load the
   measured streams. Each fixture is serialized in its engine's own format and goes through the same reader,
   extraction, join and comparison path as the production run. Each control is a test that fails unless its mutation
   changes what the comparator computes for the mutated object (a residual, a part comparison, a coverage outcome or
   the accounting) while the other operand is unchanged. A change in the printed operand alone does not satisfy a
   control, and each control fails when the comparison step it targets is disabled. No assertion names a value.
   Required mutations:
   - one-sided corruption of an operand;
   - a form change, not a rescaling;
   - a coefficient or sign change;
   - a repoint of every name-table row (item 2);
   - a stripped argument of an applied function;
   - a live function replaced by a constant;
   - a face delta's coefficient moved to the other face, or removed;
   - a difference confined to one region's regular part;
   - a difference confined to a stratum where a jump or the interior order's face value vanishes;
   - a changed label;
   - a leaf that is parsed but not compared (it must show in the leaf counts);
   - changed derivative structure, including a derivative of a composite argument;
   - a nested sibling removed;
   - on one side only, a factor that equals one on the supplied domain changed or removed (the as-computed residual
     must move);
   - a repoint of a component correspondence (item 1);
   - a joined row whose payloads differ in structure (it must be printed as not formed, not compared position by
     position).

   Coverage controls: for each of these reductions there is a fixture difference that it would miss. The reductions
   are evaluation at one point only, under one test-function form only, and under forms that do not separate a
   function's arguments. Each test fails if the comparator's coverage is reduced in that way.

   One equivalence control: a delta's coefficient is changed away from its face only, so the distribution is
   unchanged. The test fails if the delta part's comparison moves.

   One further test: a duplicate join or name row is rejected.

   Run on the real streams, the comparator also prints per-row accounting (joined, or unjoined with its reason) and
   two leaf counts per joined row: one enumerated independently from each stream's parse, and the count actually
   compared. A row compared below its parsed count stays visible.

7. **Running.** Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode (`--pool s12-comparator
   --memory-gib 50`, with the log directory under `_scratch/s11c/`), one job at a time. There are no time limits.
   Stop and report any kill or admission refusal; never answer one by narrowing what is compared. Write output under
   your own scratch directory, never under `scripts/out/` or `mathematica/out/`.

## Exclusions
- No engine edit, no new physics, and no choice of form beyond the test functions of item 4. No commit.
- The engines' intermediate objects that only one engine prints are not paired with anything (item 1).

## Builder report (at most 40 lines)
- The name table, each row with its object and two engine-line citations.
- The component correspondences, each with its two engine-line citations, and every joined row that was not compared
  or could not be formed, with its reason.
- Per-row accounting and leaf counts; every unjoined tag with its reason.
- The test functions and points, and the strata covered.
- Which machinery you reused, and what you re-verified on these streams.
- The commands for, and outcomes of, the tests and the full guarded run.
- Runtime and peak memory.
- That no residual target was given.
