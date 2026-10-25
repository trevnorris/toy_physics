# S15 cross-engine comparator (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-24. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The mechanical lookups behind this list are in
`directives/_measurements/S15_comparator_setup_lookups.md`; "L-n" cites a section of that file.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). No builder starts
before then. The comparator it describes is an instrument written by Codex, reviewed by two non-author build legs
(fresh Claude + Grok) until clear.

## Object
`scripts/S15_cross_engine_comparator.py`, with synthetic tests in `scripts/test_S15_cross_engine_comparator.py`
(neither exists yet, L10). It compares the two S15 engines' printed objects, object by object, piece by piece and
group by group, and prints what it finds. It decides nothing.

## Inputs (read-only)
- The two transcripts accepted at `3192123b` (L1). They are git-annex content; run `datalad get` on both if their
  content is absent.
  - `scripts/out/S15_two_drain_region.out` (SymPy);
  - `mathematica/out/S15_two_drain_region.out` (Wolfram).
- The two engine sources at the same commit, `scripts/S15_two_drain_region.py` and
  `mathematica/S15_two_drain_region.wl`. They are the authority for what each printed object, head, alias and name
  means. Do not run either engine.
- The S15 decision list, `directives/S15_decision_list.md`, revision 3 (`91248b67`): "Supplied equations", "The
  drains, their fields and the region" and "What the build computes" (L12).
- An earlier comparator, as machinery to adapt (L10): `scripts/S16_cross_engine_comparator.py` and
  `scripts/test_S16_cross_engine_comparator.py`, accepted at `f08db34a`. Reuse only what is consistent with items 1–8
  below, and re-verify it on these streams.

## What must be true

1. **Each stream is read as its engine prints it.**
   - Every printed line is assigned to a tag's payload or to a progress marker (`WORKING_ON`). That includes payloads
     that span several lines, and tags that contain spaces inside braces (L2, L3). A line the reader cannot assign is
     an operational failure.
   - The SymPy stream repeats some tags as records (L2, L3): `dimension_argument/name`, `/operand_a`, `/operand_b` and
     `/residuals` (208 records), and `partition/unlabelled_coefficient` (27 lines). Each record, and each repeated
     line, is read as its own object, identified by its content or its place as the source defines it, never by its
     position in the stream alone.
   - Where the two engines print the same object in different notations, the reader builds the same object from both.
     The notations differ in at least these ways (L4, L6, L7):
     - the integral heads (`Integral(…)` against `Inactive[Integrate][…]`), with each engine's chart and measure lines
       and source defining the integration variables, ranges and measure;
     - unevaluated substitutions (`Subs`) and dummy variables, and derivatives at composite arguments;
     - relations and equations (`Eq(…)` against `==`, a conjunction of inequalities against a chained inequality);
     - matrices, lists, tuples and associations;
     - labels (text against quoted strings);
     - dimensions (exponent triples over each engine's printed basis, whose names differ);
     - group names (`labels_drain1+drain2` against `labels_{drain1, drain2}`).

     The exact residual is then the residual of the printed objects.
   - A different object is still read as different. Examples: another derivative slot, order or point; another
     integrand; another surface or region; another label set; another relation or condition.
   - This translates each engine's notation. It does not reconcile two representations. Content that one engine prints
     and the other does not stays in the residual.

2. **Join through an explicit binding table; a binding is not agreement.**
   - **The table.** A join is made only through a table in the comparator's source. The streams share 401 tag names
     (L2), and sharing a name is not a binding. The engines also differ in group names, in the names of the supplied
     objects and of their own checks, and in the order they print (L2, L4, L5).
     - The table may be generated from explicit rules that map name components. Each rule cites the object both ends
       denote, by the S15 decision list's item or supplied equation where the list names it, and the line in each
       engine that builds it.
     - The generated table is printed in full before the comparisons, or written beside the output with its size and
       sha256 printed.
     - Each binding joins one end in each stream. An end is a tag, a record of item 1, or a component of a payload
       named by its tag and its place in the payload.
     - The table is injective at its ends: no end is bound twice, and a repeat is reported. No binding is made by
       similarity of names, by position, or because two payloads are identical. Whether a binding is made never
       depends on the payloads' values or on their comparison.
   - **A binding is never agreement.** Agreement comes only from the content compared beneath it.
   - **Pieces and groups one by one.** Each object the S15 decision list asks for is its own binding, per component
     and at each surface:
     - item 1's four rest pieces and two transport pieces, and each label group of each grouped piece, all seven,
       including a group an engine prints as zero;
     - the whole pressure, the expansion's groups and the whole pressure minus the expansion; `Π_o`, whole;
     - item 2's rate, each surface's parts, the integrated force density, the region's momentum and its groups, and
       the rate's groups.

     Totals are compared as well. Agreement of a total is never reported as agreement of its pieces, and agreement of
     a piece never as agreement of its groups. The engines' own partition residuals show completeness only: a term
     filed in the wrong group leaves them, and the piece, unchanged (L11). Filing is compared only group by group.
   - **Accounting.**
     - Each tag occurs once per stream, apart from progress markers and item 1's records; any other repeat is
       reported.
     - Every tag and record is bound, or is listed as unjoined with its engine and the reason. Unjoined tags may be
       listed by grouped reason, with every member countable from the output. Every component of a bound end is bound,
       or is printed unpaired with its reason. Nothing is skipped silently.
     - An object that both engines print and the S15 decision list names is compared. It is left unjoined only where
       the cited lines build objects the S15 decision list distinguishes, and its reason names both.
   - **The engines' own checks.** These are each engine's operand and residual triples, partition checks, dimension
     operands and pair residuals, comparison notes and the balance reduction label (L8). They check that engine's
     construction. A binding joins two of them only where both check the same object against the same counterpart.
     Otherwise each is listed as unjoined with that reason. The region balance operands, residual and equation are
     objects of item 2 and are compared as such.
   - **Components.**
     - A binding can join payloads that are not constituted the same way. Before bound payloads are compared as values,
       the comparator establishes which of their components denote the same object, citing the line in each engine
       that builds them. A component pair is made only where the cited lines build the same object. Otherwise the
       components stay printed, unpaired. A binding without such a correspondence is listed as not compared, with its
       reason.
     - A structural mismatch is printed as a comparison that could not be formed. Examples: different component counts,
       or different kinds of object. No pairing is made by position alone.
   - **Output.** Both operands stay printed under their original tags. The comparator prints the table, each binding
     with its citations, before its comparisons.

3. **Symbols, functions and aliases through an explicit table.** Any map between the engines' names is an explicit,
   injective table in the comparator's source (L9). It covers:
   - the coordinates and time, and the chart variables of the surfaces and of the region;
   - both centres `X_1`, `X_2`, component by component, and both radii;
   - every live field, component by component: `u_k`, `ρ_k`, `ϱ`, `W^i`, `p_br`, `η`, `ζ`, `Π_o^{ij}` with its index
     order, `f^i`, and the free `c_comp` symbol;
   - the label markers and dummy variables where they are printed;
   - the alias symbols that name printed terms, resolved through their printed definitions;
   - the dimension basis.

   - Each entry cites the object both names denote and the line in each engine that introduces it.
   - A mapped name is a substitution, never by itself a zero. A repoint control shows that pointing a mapped name at a
     different object moves the printed residual.
   - The comparator prints this table, with its citations, before its comparisons.

4. **Integrals and conditions compared as the objects they are.**
   - An integral is compared as what it denotes: its integrand, surface or region, and measure together. Each engine's
     chart and measure lines and source define them (L4). Its value is compared through the evaluation of item 5.
     - Two integrals whose integrands differ only in how they write the same surface point, normal or region point
       denote the same object.
     - An integrand-level difference may also be printed. If so, it is a separately labelled representation diagnostic,
       never the comparison of the integral.
   - Surfaces are paired by what they are, the inner sphere or the outer, as each engine's chart lines define them.
     They are never paired by index or by order.
   - **Conditions** (L6). The domain `|X_1 − X_2| > R > r_in > 0` is compared as a relation. The pressure conditions,
     the regularity and chart lines and the other prose labels are compared as text, side by side. The comparator never
     reads a condition as satisfied. A value printed under a condition is compared where the test functions of item 5
     meet that condition.

5. **Residuals printed before any guard.** For each compared value, piece or group, print the SymPy operand, the
   Wolfram operand and their residual.
   - **Exact first.** Use exact symbolic subtraction where feasible, and print the exact residual as computed. Where it
     is feasible, evaluation at points is added beside it, never in its place. Otherwise evaluate at several generic
     rational points (E1; numeric-PIT practice).
   - **Test functions.** Under evaluation, every live function is replaced by test functions of several different
     forms, applied identically to both operands. That covers each field of item 3, both centres and both radii, with
     their time derivatives.
     - The test functions are not required to satisfy the supplied balance. Print the test functions, the points and
       both operands' values.
     - Distinct live functions get distinct test functions in every form and stratum. That includes `Π_o^{ij}` and
       `Π_o^{ji}`, and the two drains' profiles.
     - In every form, no test function is constant in any of its arguments, and none vanishes identically on a surface
       or in the region.
     - For a function of several arguments, the test functions tell its arguments apart: at the evaluation points,
       its derivatives in different arguments differ.
     - A derivative taken at another point is a different object (item 1), and evaluation sees it. For every compared
       value whose exact residual is not resolved, a live function's derivative taken at another point, in any one
       argument, moves the samples in at least one form and stratum. That includes a derivative inside a
       substitution, a derivative at a composite argument, and a derivative inside a surface or volume integral.
     - The geometry stays inside the printed domain at every evaluation time: `|X_1 − X_2| > R > r_in > 0`.
   - **Integrals.** An integral is evaluated after the substitution, exactly where feasible. Otherwise it is evaluated
     numerically, and the working precision and an error estimate are printed. An integral that cannot be evaluated is
     printed as unresolved coverage.
   - **Strata.** The evaluation covers the generic case, and strata where, at the evaluation time, the velocity of
     `X_2` vanishes, where both radii's velocities vanish, where both vanish, and where the velocity of `X_1` vanishes.
     In each, the higher time derivatives stay generic. A stratum the comparator cannot evaluate is printed as
     unresolved coverage.
   - **No supplied equation applied.** No supplied equation is applied to an operand to make two operands meet. That
     covers the momentum balance, the forms of `G`, `Π` and `T`, and the split of the fields. Where a pair differs by a
     supplied equation, the residual as computed is printed. The residual after that equation is applied may be
     printed beside it, naming the equation and its source line, never in its place.
   - **Values and labels.**
     - A dimension object is compared as a value, through the basis of item 3.
     - A label (a comparison note, a reduction label, a condition, a prose line) is compared as text, after item 1's
       reading.
     - A binding key is neither printed nor counted as a compared value.
   - **No collapse.** No applied function is collapsed to a bare symbol, and no live function is replaced by a constant,
     to make two operands meet. A difference in representation is printed as a difference, never reconciled inside the
     comparator.

6. **No target.** Nothing in the comparator or its report states what any residual on the measured streams is
   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
   and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output belongs to the
   S15 record.

7. **Controls that fail when their defect goes undetected.**
   - The tests use synthetic fixtures only, and never load the measured streams.
   - Each fixture is serialized in its engine's own format, and goes through the same reader, extraction, join and
     comparison path as the production run.
   - Each control is a test that fails unless its mutation changes what the comparator computes for the mutated object,
     while the other operand is unchanged. What it computes includes a residual, a piece or group comparison, a
     condition comparison, a coverage outcome and the accounting.
   - A change in the printed operand alone does not satisfy a control. Each control fails when the comparison step it
     targets is disabled. No assertion names a value.

   Required mutations:
   - one-sided corruption of an operand;
   - a form change, not a rescaling;
   - a coefficient or sign change;
   - a repoint of every binding rule (item 2) and of every name-table row (item 3). This includes swapping the inner and
     outer surfaces, shifting a component index, transposing `Π_o`, exchanging the two drains, and pointing a group at
     another label set. Each is checked by evaluation as well as by the exact residual;
   - two distinct live functions given the same test function;
   - a stripped argument of an applied function;
   - a live function replaced by a constant;
   - on one side, a term moved from one group to another with the piece unchanged (both group comparisons must move);
   - on one side, a term moved from one piece to another with the total unchanged (both piece comparisons must move);
   - a difference confined to one group, one piece, or one component;
   - a difference confined to a stratum of item 5;
   - on one side, an integral's surface or region changed with its integrand kept, and its integrand changed with its
     surface or region kept;
   - on one side, the domain relation changed, and a condition's text changed;
   - a changed line within a payload that spans several lines;
   - a changed label, in each engine's own label serialization and through the production label path;
   - in a repeated record of item 1, one record's payload changed, and one record removed;
   - a leaf that is parsed but not compared (it must show in the leaf counts);
   - changed derivative structure, including a derivative inside a substitution and a derivative at a composite
     argument;
   - on one side, a live function's derivative taken at another point, in each of its arguments in turn: inside a
     substitution, at a composite argument, and inside a surface and a volume integral (its samples must move);
   - a nested sibling removed;
   - on one side only, a supplied equation applied (the as-computed residual must move);
   - a repoint of a component correspondence (item 2), through each place where the comparator pairs components;
   - an alias resolved to another printed term on one side;
   - split-against-grouped packaging: one side prints under several tags or records what the other prints as
     components of one tag, and each counterpart component is mutated on its own (that component's comparison must
     move);
   - a binding whose payloads differ in structure (it must be printed as not formed, not compared position by
     position).

   **Coverage controls.** Each of these reductions has a fixture difference that it would miss:
   - evaluation at one point only;
   - one test-function form only;
   - forms that do not tell a function's arguments apart;
   - forms in which a derivative in one argument does not depend on the other arguments;
   - a form in which a test function vanishes on a surface or in the region;
   - the generic stratum only.

   Each coverage control fails under its own reduction, run alone.

   **Equivalence controls.** Each is a test that fails if the comparison moves:
   - the same object in each engine's notation, for each notation item 1 names;
   - a surface integral whose integrand writes the surface point or normal another way;
   - item 1's repeated records printed in another order.

   **A further test:** a duplicate binding end, rule or name row is rejected.

   **On the real streams,** the comparator also prints per-tag accounting (bound, or unjoined with its reason). For
   each binding it prints two leaf counts: one enumerated independently from each stream's parse, and the count
   actually compared. A binding compared below its parsed count stays visible.

8. **Running.**
   - Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode (`--pool s15-comparator --memory-gib 50`),
     one job at a time. The `--log-directory` must resolve under `_scratch/s11c/`; the guard refuses any other.
   - There are no time limits. A long run prints progress.
   - Stop and report any kill or admission refusal. Never answer one by narrowing what is compared.
   - Write output under your own scratch directory, never under `scripts/out/` or `mathematica/out/`.

## Exclusions
- No engine edit, no new physics, and no choice of form beyond the test functions of item 5. No commit.
- The engines' intermediate objects that only one engine prints are not paired with anything (item 2).

## Builder report (at most 40 lines)
- Where the binding rules, the generated binding table and the name table are, their sizes, and where the output
  prints them.
- Every binding that was not compared or could not be formed, with its reason. Every unjoined tag's reason, grouped.
  Every object the S15 decision list names that is left unjoined, with its reason.
- Per-tag accounting and leaf-count totals, and every binding compared below its parsed count.
- The test functions and points, the strata covered, and the integrals evaluated numerically or left unresolved.
- Which machinery you reused, and what you re-verified on these streams.
- The commands for, and outcomes of, the tests and the full guarded run.
- Runtime and peak memory.
- That no residual target was given.
