# S14a part 1 cross-engine comparator (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The mechanical lookups behind this list are in
`directives/_measurements/S14a_1_comparator_setup_lookups.md`; "L-n" cites a section of that file.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). No builder starts
before then. The comparator it describes is an instrument written by Codex, reviewed by two non-author build legs
(fresh Claude + Grok) until clear.
- **Pass 1** (Codex + Grok) returned one finding, K1, confirmed. The list as reviewed is preserved at `f1b41e77`, with
  `_measurements/S14a_1_comparator_decision_list_review_p1_{disposition,lookups}.md`.
- **This revision** folds K1 into item 2, item 7 and the builder report. It is the one fold G2 allows.

## Object
`scripts/S14a_1_cross_engine_comparator.py`, with synthetic tests in `scripts/test_S14a_1_cross_engine_comparator.py`
(neither exists yet, L7). It compares the two S14a part 1 engines' printed objects, object by object, and prints what
it finds. It decides nothing.

## Inputs (read-only)
- The two transcripts accepted at `94fddf2f`, which are git-annex content (run `datalad get` on both if their content
  is absent; sha256s in L1):
  - `scripts/out/S14a_1_projection.out` (SymPy);
  - `mathematica/out/S14a_1_projection.out` (Wolfram).
- The two engine sources at the same commit, `scripts/S14a_1_projection.py` and `mathematica/S14a_1_projection.wl`.
  They are the authority for what each printed object, head and name means. Do not run either engine.
- The S14a part 1 decision list, `directives/S14a_1_decision_list.md`, revision 3 (`e9eed614`): "Supplied equations"
  and "What the build computes".
- Earlier comparators, as machinery to adapt (L6): `scripts/S12_cross_engine_comparator.py` and
  `scripts/test_S12_cross_engine_comparator.py`, accepted at `fe485b73`; `scripts/S9b_cross_engine_comparator.py` and
  `scripts/test_S9b_cross_engine_comparator.py`, accepted at `88bb1f41`. Reuse only what is consistent with items 1–8
  below, and re-verify it on these streams.

## What must be true

1. **Each stream is read as its engine prints it.**
   - Every printed line is assigned to a tag's payload or to a progress marker. That includes payloads that span
     several lines (L3). A line the reader cannot assign is an operational failure.
   - Where the two engines print the same object in different notations, the reader builds the same object from both.
     Examples: a derivative at composite arguments, an unevaluated integral, a relation with its sides exchanged, a
     condition attached to an object, a list of strings (L4, L5). The exact residual is then the residual of the
     printed objects.
   - A different object is still read as different: another derivative slot, order or point; another integrand,
     variable or range; another relation or condition.
   - This translates each engine's notation. It does not reconcile two representations, and content that one engine
     prints and the other does not stays in the residual.

2. **Join through an explicit binding table; a binding is not agreement.** The two streams share no tag name (L2). A
   join is made only through a table in the comparator's source.
   - Each binding joins one end in each stream. An end is a tag, or a component of a tag's payload named by its tag
     and its place in the payload: one engine may print under several tags what the other prints as components of one
     tag. Each binding cites the object both ends denote, by the S14a part 1 decision list's item or supplied equation
     where the list names it, and the line in each engine that builds it.
   - The table is injective at its ends: no tag or component is bound twice, and a repeat is reported. No binding is
     made by similarity of names, by position, or because two payloads are identical.
   - A binding is never agreement: agreement comes only from the content compared beneath it.
   - Each tag occurs once per stream, and a repeat is reported. Every tag is bound, or is listed as unjoined with its
     engine and the reason. Every component of a bound tag is bound, or is printed unpaired with its reason. Nothing
     is skipped silently.
   - An object that both engines print and the S14a part 1 decision list names is compared. It is left unjoined only
     where the cited lines build different objects, and its reason says so.
   - An engine's own checks (its comparison notes, operand pairs and residuals) check that engine's construction. A
     binding joins two of them only where both check the same object against the same counterpart. Otherwise each is
     listed as unjoined with that reason.
   - A binding can join payloads that are not constituted the same way. Before bound payloads are compared as values,
     the comparator establishes which of their components denote the same object, citing the line in each engine that
     builds them. A component pair is made only where the cited lines build the same object; otherwise the components
     stay printed, unpaired. A binding without such a correspondence is listed as not compared, with its reason.
   - A structural mismatch is printed as a comparison that could not be formed, for example different component
     counts or different kinds of object. No pairing is made by position alone.
   - Both operands stay printed under their original tags. The comparator prints the table, each binding with its
     citations, before its comparisons.

3. **Symbols and functions through an explicit table.** Any map between the engines' names is an explicit, injective
   table in the comparator's source. That covers the coordinates of both settings and time, the faces, each
   balance's regional fields (`n`, the components of `u` and `J_χ`, `χ_B`, the rates) region by region, the shell's
   radii, and the dimension symbols.
   - Each entry cites the object both names denote and the line in each engine that introduces it.
   - Where the engines split a field into regional pieces or components differently, the table states the
     correspondence piece by piece.
   - A mapped name is a substitution, never by itself a zero. A repoint control shows that pointing a mapped name at a
     different object moves the printed residual.
   - The comparator prints this table, with its citations, before its comparisons.

4. **Distributions, integrals and conditions compared as the objects they are.**
   - An object that carries steps or deltas in `w` is compared in parts: its regular part on each region, and each
     face delta's part. Each part is compared with its counterpart.
     - A delta's coefficient carries data about the distribution only on that delta's support. A delta's part is
       therefore its coefficient restricted to the delta's support (its value on the face), with the delta's
       normalization in `w`. The comparator derives it from the printed distribution itself.
     - A raw coefficient, unrestricted, may also be compared. If so, it is a separately labelled representation
       diagnostic, never the delta's part.
     - Content that falls in no such part (another singular atom, or a part the comparator cannot isolate) is printed
       as that outcome, as a coverage finding.
   - An unevaluated integral is compared as what it denotes: its integrand, variable and range together, through the
     evaluation of item 5. Two integrals that differ by a change of integration variable denote the same object.
   - Where an integrand carries a delta at an endpoint of its integral's range, the delta is counted in full, as the
     S14a part 1 decision list's item 1 defines the closed region ("each face's concentrated part counted in full").
     That evaluation does not depend on either CAS's convention for a delta at an endpoint. Each such integral is
     counted and listed.
   - An object that an engine prints with a condition (L5) is compared as a condition and a value. The two conditions
     are compared with each other, and the values where both conditions hold. Points where one condition holds and the
     other does not are printed.

5. **Residuals printed before any guard.** For each compared value or part, print the SymPy operand, the Wolfram
   operand and their residual.
   - Use exact symbolic subtraction where feasible, and print the exact residual as computed. Where it is feasible,
     evaluation at points is added beside it, never in its place. Otherwise evaluate at several generic rational
     points (E1; numeric-PIT practice).
   - Under evaluation, every live function is replaced by test functions of several different forms, applied
     identically to both operands. That includes each regional field, the interior order, the rates and both faces.
     The test functions are not required to satisfy the supplied balances. Print the test functions, the points and
     both operands' values.
   - Distinct live functions get distinct test functions, in every form and stratum. The exception is the deliberately
     equal one-sided values of a jump stratum.
   - For a function of several arguments, the test functions separate each argument's dependence: at the evaluation
     points, its derivatives in different arguments differ. That includes each face's dependence on position and
     time.
   - An integral is evaluated after the substitution, exactly where feasible. Otherwise it is evaluated numerically, and
     the working precision and an error estimate are printed. An integral that cannot be evaluated is printed as
     unresolved coverage.
   - Points lie in each region and on each face, inside the domains the engines print as supplied. A residual's
     dependence on those domains is printed: the residual as computed, and its value on the supplied domain, side by
     side.
   - The evaluation covers strata where a field's jump at a face vanishes, where the interior order's value at a face
     vanishes, and where a condition of item 4 fails, as well as the generic case. A stratum the comparator cannot
     evaluate is printed as unresolved coverage.
   - No supplied equation (a committed balance, the identity for `Γ_B`, U1, the face factor) is applied to an operand
     to make two operands meet. Where a pair differs by a supplied equation, the residual as computed is printed. The
     residual after that equation is applied may be printed beside it, naming the equation and its source line, never
     in its place.
   - A dimension object is compared as a value. A label (a premise line, a scope line, a comparison note) is compared as
     text, after item 1's reading.
   - A binding key is neither printed nor counted as a compared value.
   - No applied function is collapsed to a bare symbol, and no live function is replaced by a constant, to make two
     operands meet. A difference in representation is printed as a difference, never reconciled inside the
     comparator.

6. **No target.** Nothing in the comparator or its report states what any residual on the measured streams is
   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
   and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output belongs to the
   S14a part 1 record.

7. **Controls that fail when their defect goes undetected.** The tests use synthetic fixtures only and never load the
   measured streams. Each fixture is serialized in its engine's own format and goes through the same reader,
   extraction, join and comparison path as the production run. Each control is a test that fails unless its mutation
   changes what the comparator computes for the mutated object (a residual, a part comparison, a condition comparison,
   a coverage outcome or the accounting) while the other operand is unchanged. A change in the printed operand alone
   does not satisfy a control, and each control fails when the comparison step it targets is disabled. No assertion
   names a value. Required mutations:
   - one-sided corruption of an operand;
   - a form change, not a rescaling;
   - a coefficient or sign change;
   - a repoint of every binding (item 2) and of every name-table row (item 3), checked by evaluation as well as by the
     exact residual;
   - two distinct live functions given the same test function;
   - a stripped argument of an applied function;
   - a live function replaced by a constant;
   - a face delta's coefficient moved to the other face, or removed;
   - a difference confined to one region's regular part;
   - a difference confined to a stratum of item 5;
   - on one side, an integral's range or variable changed with its integrand kept, and its integrand changed with its
     range kept;
   - on one side, a condition changed with its value kept;
   - a changed line within a payload that spans several lines;
   - a changed label, in each engine's own label serialization and through the production label path, for a label tag
     and for a label that is a list of strings;
   - a leaf that is parsed but not compared (it must show in the leaf counts);
   - changed derivative structure, including a derivative of a composite argument;
   - a nested sibling removed;
   - on one side only, a factor that equals one on the supplied domain changed or removed (the as-computed residual
     must move);
   - on one side only, a supplied equation applied (the as-computed residual must move);
   - a repoint of a component correspondence (item 2), through each place where the comparator pairs components;
   - split-against-grouped packaging: one side prints under several tags what the other prints as components of one
     tag, and each counterpart component is mutated on its own (that component's comparison must move);
   - a binding whose payloads differ in structure (it must be printed as not formed, not compared position by
     position).

   Coverage controls: for each of these reductions there is a fixture difference that it would miss. The reductions
   are evaluation at one point only, under one test-function form only, and under forms that do not separate a
   function's arguments. Each coverage control fails under its own reduction, run alone.

   Equivalence controls, each a test that fails if the comparison moves:
   - a delta's coefficient changed away from its face only, so the distribution is unchanged;
   - the same object in each engine's notation (item 1), for each example item 1 names;
   - an integral, against the same integral after a change of integration variable;
   - a delta at an endpoint inside an integral on one side, against its face part counted in full outside the integral
     on the other.

   One further test: a duplicate binding end or name row is rejected.

   Run on the real streams, the comparator also prints per-tag accounting (bound, or unjoined with its reason) and
   two leaf counts per binding: one enumerated independently from each stream's parse, and the count actually
   compared. A binding compared below its parsed count stays visible.

8. **Running.** Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode (`--pool s14a1-comparator
   --memory-gib 50`, with the log directory under `_scratch/s11c/`), one job at a time. There are no time limits.
   Stop and report any kill or admission refusal; never answer one by narrowing what is compared. Write output under
   your own scratch directory, never under `scripts/out/` or `mathematica/out/`.

## Exclusions
- No engine edit, no new physics, and no choice of form beyond the test functions of item 5. No commit.
- The engines' intermediate objects that only one engine prints are not paired with anything (item 2).

## Builder report (at most 40 lines)
- Where the binding table and the name table are in the source, their sizes, and where the output prints them.
- Every binding that was not compared or could not be formed, with its reason; every unjoined tag's reason, grouped;
  every object the S14a part 1 decision list names that is left unjoined, with its reason.
- Per-tag accounting and leaf-count totals, and every binding compared below its parsed count.
- The test functions and points, the strata covered, and the integrals evaluated numerically or left unresolved.
- Which machinery you reused, and what you re-verified on these streams.
- The commands for, and outcomes of, the tests and the full guarded run.
- Runtime and peak memory.
- That no residual target was given.
