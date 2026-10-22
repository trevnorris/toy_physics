# S16 cross-engine comparator (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-22. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The mechanical lookups behind this list are in
`directives/_measurements/S16_comparator_setup_lookups.md`; "L-n" cites a section of that file.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). No builder starts
before then. The comparator it describes is an instrument written by Codex, reviewed by two non-author build legs
(fresh Claude + Grok) until clear.
- **Pass 1** (Codex + Grok) returned two findings, K1 and K2, both confirmed. The list as reviewed is preserved at
  `958a55a9`, with `_measurements/S16_comparator_decision_list_review_p1_{disposition,lookups}.md`.
- **This revision** folds K1 and K2 into item 5 and item 7. It is the one fold G2 allows.

## Object
`scripts/S16_cross_engine_comparator.py`, with synthetic tests in `scripts/test_S16_cross_engine_comparator.py`
(neither exists yet, L10). It compares the two S16 engines' printed objects, object by object and piece by piece, and
prints what it finds. It decides nothing.

## Inputs (read-only)
- The two transcripts accepted at `1d990ee6` (L1). They are git-annex content; run `datalad get` on both if their
  content is absent.
  - `scripts/out/S16_moving_region.out` (SymPy);
  - `mathematica/out/S16_moving_region.out` (Wolfram).
- The two engine sources at the same commit, `scripts/S16_moving_region.py` and `mathematica/S16_moving_region.wl`.
  They are the authority for what each printed object, head, alias and name means. Do not run either engine.
- The S16 decision list, `directives/S16_decision_list.md`, revision 4 (`0464b80a`): "Supplied equations", "The
  region" and "What the build computes" (L12).
- The build's review dispositions, `directives/_measurements/S16_build_review_r{0,1}_disposition.md`. Their notes on
  comparison are in L11.
- Earlier comparators, as machinery to adapt (L10): `scripts/S14a_1_cross_engine_comparator.py` and
  `scripts/test_S14a_1_cross_engine_comparator.py`, accepted at `de81b9fc`. Reuse only what is consistent with items
  1–8 below, and re-verify it on these streams.

## What must be true

1. **Each stream is read as its engine prints it.**
   - Every printed line is assigned to a tag's payload or to a progress marker. That includes payloads that span
     several lines, and tags that contain spaces (L3). A line the reader cannot assign is an operational failure.
   - Where the two engines print the same object in different notations, the reader builds the same object from both.
     The notations differ in at least these ways (L4, L6, L7):
     - the integral heads (`Volume`, `Surface(k, …)` against `V`, `B[k, …]`), as each engine's `integral.*` lines and
       source define them;
     - unevaluated substitutions (`Subs`), and derivatives at composite arguments;
     - relations and definitions (an equation against a rule);
     - conditions and domains;
     - labels (text against `Missing[…]`), and lists against associations;
     - dimensions (unit symbols against exponent vectors over a printed basis).

     The exact residual is then the residual of the printed objects.
   - A different object is still read as different. Examples: another derivative slot, order or point; another
     integrand; another surface; another region; another relation or condition.
   - This translates each engine's notation. It does not reconcile two representations. Content that one engine prints
     and the other does not stays in the residual.

2. **Join through an explicit binding table; a binding is not agreement.**
   - **The table.** A join is made only through a table in the comparator's source. The streams share 786 tag names
     (L2), and sharing a name is not a binding. The engines also differ in index base, surface numbering and piece names
     (L5).
     - The table may be generated from explicit rules that map name components. Each rule cites the object both ends
       denote, by the S16 decision list's item or supplied equation where the list names it, and the line in each
       engine that builds it.
     - The generated table is printed in full before the comparisons, or written beside the output with its size and
       sha256 printed.
     - Each binding joins one end in each stream. An end is a tag, or a component of a tag's payload named by its tag
       and its place in the payload.
     - The table is injective at its ends: no tag or component is bound twice, and a repeat is reported. No binding is
       made by similarity of names, by position, or because two payloads are identical.
   - **A binding is never agreement.** Agreement comes only from the content compared beneath it.
   - **Pieces one by one.** Each piece the decision list asks for is its own binding, at every level the engines print
     it. That covers each rate's rest, centre-motion and radius-motion parts at each surface, its bulk parts, and every
     part of every printed derivative.
     - Totals are compared as well.
     - Agreement of a total is never reported as agreement of its pieces. Merging two pieces leaves the totals
       unchanged (L11).
   - **Accounting.**
     - Each tag occurs once per stream, apart from progress markers, and a repeat is reported.
     - Every tag is bound, or is listed as unjoined with its engine and the reason. Unjoined tags may be listed by
       grouped reason, with every member countable from the output. Every component of a bound tag is bound, or is
       printed unpaired with its reason. Nothing is skipped silently.
     - An object that both engines print and the S16 decision list names is compared. It is left unjoined only where
       the cited lines build different objects, and its reason says so.
   - **The engines' own checks.** These are each engine's comparison notes, dimension sums, and partition operands and
     residuals (L8). They check that engine's construction. A binding joins two of them only where both check the same
     object against the same counterpart. Otherwise each is listed as unjoined with that reason.
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
   injective table in the comparator's source. It covers:
   - the coordinates and time, and the chart and radial variables;
   - `X_c`, both radii, `N`, the dipole and `P`;
   - every live field, component by component: `Σ`, `F^i`, `G^i`, `Π^{ij}` with its index order, `f^i`, and the face
     quantities;
   - the alias symbols that name printed terms, resolved through their printed definitions;
   - the dimension basis.

   - Each entry cites the object both names denote and the line in each engine that introduces it.
   - A mapped name is a substitution, never by itself a zero. A repoint control shows that pointing a mapped name at a
     different object moves the printed residual.
   - The comparator prints this table, with its citations, before its comparisons.

4. **Integrals and conditions compared as the objects they are.**
   - An integral is compared as what it denotes: its integrand, region or surface, and measure together. Each engine's
     `integral.*` lines and source define its heads (L4). Its value is compared through the evaluation of item 5.
     - Two integrals whose integrands differ only in how they write the same surface point or normal denote the same
       object.
     - An integrand-level difference may also be printed. If so, it is a separately labelled representation diagnostic,
       never the comparison of the integral.
   - Surfaces are paired by what they are, the inner sphere or the outer, as each engine's integral and region lines
     define them. They are never paired by index.
   - An object that an engine prints with a condition or on a domain (L6) is compared as a condition and a value. The
     two conditions are compared with each other, and the values where both hold. Points where one holds and the other
     does not are printed.

5. **Residuals printed before any guard.** For each compared value or piece, print the SymPy operand, the Wolfram
   operand and their residual.
   - **Exact first.** Use exact symbolic subtraction where feasible, and print the exact residual as computed. Where it
     is feasible, evaluation at points is added beside it, never in its place. Otherwise evaluate at several generic
     rational points (E1; numeric-PIT practice).
   - **Test functions.** Under evaluation, every live function is replaced by test functions of several different
     forms, applied identically to both operands. That covers each field, `X_c`, and both radii with their time
     derivatives.
     - The test functions are not required to satisfy the supplied balances. Print the test functions, the points and
       both operands' values.
     - Distinct live functions get distinct test functions in every form and stratum. That includes `Π^{ij}` and
       `Π^{ji}`.
     - For a function of several arguments, the test functions tell its arguments apart: at the evaluation points,
       its derivatives in different arguments differ.
     - A derivative taken at another point is a different object (item 1), and evaluation sees it. For every compared
       value whose exact residual is not resolved, a live function's derivative taken at another point, in any one
       argument, moves the samples in at least one form and stratum. That includes a derivative inside a
       substitution, a derivative at a composite argument, and a derivative inside a volume or surface integral.
     - The geometry stays inside the printed domains: `0 < r_in < R`.
   - **Integrals.** An integral is evaluated after the substitution, exactly where feasible. Otherwise it is evaluated
     numerically, and the working precision and an error estimate are printed. An integral that cannot be evaluated is
     printed as unresolved coverage.
   - **Domains.** A residual's dependence on the printed domains is printed: the residual as computed, and its value on
     the supplied domain, side by side. A normalized object is evaluated where `N` is nonzero, and that is printed.
     The unnormalized balances and moments keep the full domain, including where `N` is zero (S16 decision list, the
     paragraph after item 5).
   - **Strata.** The evaluation covers the generic case, and strata where the velocity of `X_c` vanishes at the point,
     where both radii's velocities vanish, and where both vanish. In each, the higher time derivatives stay generic.
     It also covers a stratum where `N` vanishes at the point, with the fields and the geometry live. There the
     unnormalized objects are compared, and each normalized object is printed as outside its domain. A stratum the
     comparator cannot evaluate is printed as unresolved coverage.
   - **No supplied equation applied.** No supplied equation is applied to an operand to make two operands meet. That
     covers the two balances and the face factor. Where a pair differs by a supplied equation, the residual as
     computed is printed. The residual after that equation is applied may be printed beside it, naming the equation
     and its source line, never in its place.
   - **Values and labels.**
     - A dimension object is compared as a value, through the basis of item 3.
     - A label (a schema entry, a coverage line, a comparison note) is compared as text, after item 1's reading.
     - A binding key is neither printed nor counted as a compared value.
   - **No collapse.** No applied function is collapsed to a bare symbol, and no live function is replaced by a constant,
     to make two operands meet. A difference in representation is printed as a difference, never reconciled inside the
     comparator.

6. **No target.** Nothing in the comparator or its report states what any residual on the measured streams is
   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
   and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output belongs to the
   S16 record.

7. **Controls that fail when their defect goes undetected.**
   - The tests use synthetic fixtures only, and never load the measured streams.
   - Each fixture is serialized in its engine's own format, and goes through the same reader, extraction, join and
     comparison path as the production run.
   - Each control is a test that fails unless its mutation changes what the comparator computes for the mutated object,
     while the other operand is unchanged. What it computes includes a residual, a piece comparison, a condition
     comparison, a coverage outcome and the accounting.
   - A change in the printed operand alone does not satisfy a control. Each control fails when the comparison step it
     targets is disabled. No assertion names a value.

   Required mutations:
   - one-sided corruption of an operand;
   - a form change, not a rescaling;
   - a coefficient or sign change;
   - a repoint of every binding rule (item 2) and of every name-table row (item 3). This includes swapping the inner and
     outer surfaces, shifting a component index, and transposing `Π`. Each is checked by evaluation as well as by the
     exact residual;
   - two distinct live functions given the same test function;
   - a stripped argument of an applied function;
   - a live function replaced by a constant;
   - on one side, a term moved from one piece to another with the total unchanged (that piece comparison must move);
   - a difference confined to one piece of a derivative;
   - a difference confined to a stratum of item 5;
   - on one side, an integral's surface or region changed with its integrand kept, and its integrand changed with its
     surface or region kept;
   - on one side, a condition or domain changed with its value kept;
   - a changed line within a payload that spans several lines;
   - a changed label, in each engine's own label serialization and through the production label path;
   - a leaf that is parsed but not compared (it must show in the leaf counts);
   - changed derivative structure, including a derivative inside a substitution and a derivative at a composite
     argument;
   - on one side, a live function's derivative taken at another point, in each of its arguments in turn: inside a
     substitution, at a composite argument, and inside a volume and a surface integral (its samples must move);
   - an unnormalized binding evaluated under the normalized objects' domain (its zero-count coverage must move);
   - a nested sibling removed;
   - on one side only, a factor that equals one on the supplied domain changed or removed (the as-computed residual
     must move);
   - on one side only, a supplied equation applied (the as-computed residual must move);
   - a repoint of a component correspondence (item 2), through each place where the comparator pairs components;
   - an alias resolved to another printed term on one side;
   - split-against-grouped packaging: one side prints under several tags what the other prints as components of one
     tag, and each counterpart component is mutated on its own (that component's comparison must move);
   - a binding whose payloads differ in structure (it must be printed as not formed, not compared position by
     position).

   **Coverage controls.** Each of these reductions has a fixture difference that it would miss:
   - evaluation at one point only;
   - one test-function form only;
   - forms that do not tell a function's arguments apart;
   - forms in which a derivative in one argument does not depend on the other arguments;
   - the generic stratum only;
   - evaluation of the unnormalized objects only where `N` is nonzero.

   Each coverage control fails under its own reduction, run alone.

   **Equivalence controls.** Each is a test that fails if the comparison moves:
   - the same object in each engine's notation, for each notation item 1 names;
   - a surface integral whose integrand writes the surface point or normal another way.

   **A further test:** a duplicate binding end, rule or name row is rejected.

   **On the real streams,** the comparator also prints per-tag accounting (bound, or unjoined with its reason). For
   each binding it prints two leaf counts: one enumerated independently from each stream's parse, and the count
   actually compared. A binding compared below its parsed count stays visible.

8. **Running.**
   - Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode (`--pool s16-comparator --memory-gib 50`,
     with the log directory under `_scratch/s11c/`), one job at a time.
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
  Every object the S16 decision list names that is left unjoined, with its reason.
- Per-tag accounting and leaf-count totals, and every binding compared below its parsed count.
- The test functions and points, the strata covered, and the integrals evaluated numerically or left unresolved.
- Which machinery you reused, and what you re-verified on these streams.
- The commands for, and outcomes of, the tests and the full guarded run.
- Runtime and peak memory.
- That no residual target was given.
