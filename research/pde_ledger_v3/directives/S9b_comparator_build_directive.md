# S9b cross-engine comparator (build directive / pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-12. `AGENTS.md` and `CLAUDE.md` apply. **Status:** one Codex + Grok
pass (both NEEDS REVISION; findings U1–U3 accepted), folded once. Dispositions:
`directives/_measurements/S9b_comparator_directive_review_disposition.md`.

**Gate:** this is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2), and no
builder starts before then. The comparator it describes is an instrument written by Codex. It is reviewed by two
non-author build legs (fresh Claude + Grok) until clear. The engines' harness is in a separate repair round; by the
user's decision (2026-10-12) the comparator proceeds on the engines both build-review legs found sound, which that
round does not change.

## Object

`research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py`, with synthetic tests in
`research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py`. It compares the two S9b constructions object by
object and prints what it finds. It decides nothing.

## Inputs (read-only)

- The production transcripts, which are git-annex content; run `datalad get` on both first:
  - `research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out` (SymPy);
  - `research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out` (Wolfram).
- The two engine sources at the commit that adds those transcripts. They are the authority for what each emitted
  object, head and name means. Do not run either engine.
- The spec `directives/S9b_SHARED_PHYSICS.md` (`34a994d5`), Parts A–C, and the build directive
  `directives/S9b_repair_build_directive.md`: item 4 (tags and the shared vocabulary), item 5 (`NOT_ESTABLISHED` and
  booleans), item 9 (Part C domains and the labelled density parts) and item 10 (restricted and forward cases).
- Earlier comparators, as machinery to adapt: `scripts/O2_cross_engine_comparator.py` and its tests. Reuse only
  what is consistent with items 1–8 below, and re-verify it against these streams.

## What must be true

1. **Join on the shared vocabulary; the name is a binding.** Item 4 gives both engines the same `<QUANTITY>` names
   for the same spec objects. The comparator joins each shared name across the two streams. A joined name is a
   binding, never agreement: agreement comes only from the content compared beneath it. The join is injective:
   each shared name occurs once per stream, and a repeat is reported. A shared name missing from a stream is listed
   as unjoined. `_LOCAL_` tags are not joined. Each is listed with its engine. Every emitted tag is on a joined row,
   unjoined with its reason, or listed as local. Nothing is skipped silently.

2. **Symbols and functions through an explicit table.** Any map between the two engines' symbol and function names
   (the coordinates, the impact parameter and endpoints, `GM`, `c₀`, the live profiles, the OPEN quantities, the
   Part C parameters) is an explicit, injective table in the comparator's source. Each entry cites the spec object
   both names denote and the line in each engine that introduces it. A mapped name is a substitution, never by
   itself a zero. A **repoint control** must show that pointing a mapped name at a different object moves the
   printed residual.

3. **Value and domain compared separately; labelled parts compared part by part.** Where a payload attaches a
   domain to its computed object (item 14 of the build directive defines both), the comparator prints a value
   residual and a domain comparison separately. Where a Part C row carries the labelled `ρ > 0` and `ρ ≤ 0` parts
   (item 9), each part is compared with its counterpart, and the labels themselves are compared. A
   `NOT_ESTABLISHED` marker is compared as a marker together with the domain it carries.

4. **Three-valued residuals, printed before any guard.** For each compared value, print the SymPy operand, the
   Wolfram operand and their residual. Use exact symbolic subtraction where feasible. Otherwise evaluate at several
   generic rational points where both operands are defined (E1; numeric-PIT practice).
   - The evaluation covers every stratum the engines emit, for the value and for the domain, including a stratum on
     which a parameter, a profile or a density vanishes. A stratum the comparator cannot evaluate is printed as
     unresolved coverage.
   - For an object that is a functional of a live profile, evaluate under several test profiles of different
     functional forms, each admissible under the spec, applied identically to both operands. Print the profiles,
     the spec lines that make them admissible, the points and both operands' values.
   - For a domain, print both conditions' truth values at each sampled point. Stratify a condition that is a union
     of pieces before sampling, so the sample meets every piece of either condition and of its complement.
   - A boolean-valued object (a truth value, or a condition an engine returned as its computed object; build
     directive item 5) is compared as a value. Both are printed, and their comparison is part of the value
     residual.
   - A comparison that cannot be formed (structurally mismatched or unsupported content) is printed as that
     outcome. It is a coverage finding, distinct from a parse failure.
   - No applied function is collapsed to a bare symbol, and no live profile is replaced by a constant, to make two
     operands meet.
   - A difference in representation is printed as a difference. It is never reconciled inside the comparator.

5. **No target.** Nothing in the comparator or its report states what any residual or comparison on the measured
   streams is expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`VERDICT`/`STATUS` token. It exits 0 whatever it
   finds, and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output
   belongs to the S9b record.

6. **Controls that fail when their defect goes undetected.** The tests use synthetic fixtures only and never load
   the measured streams. Each fixture is serialized in its engine's own format and goes through the same reader,
   extraction, join and comparison path as the production run. Each control is a test that fails unless its
   mutation changes what the comparator computes for the mutated object (its value residual, its domain comparison,
   its comparison outcome or its accounting) while the other operand is unchanged. A change in the printed operand
   alone does not satisfy a control, and each control fails when the comparison step it targets is disabled. The
   assertion names no value. Required mutations:
   - one-sided corruption of an operand;
   - a form change, not a rescaling;
   - a coefficient or sign change;
   - a repoint of every name-table row (item 2);
   - a stripped argument of an applied function;
   - a live profile replaced by a constant;
   - an attached domain changed;
   - a domain changed on one labelled density part only;
   - the two density labels swapped;
   - a `NOT_ESTABLISHED` marker changed;
   - the domain a `NOT_ESTABLISHED` marker carries changed;
   - a boolean-valued object flipped;
   - a difference confined to a stratum on which a parameter or a density vanishes;
   - a leaf that is parsed but not compared (it must show in the leaf counts);
   - changed derivative or binder structure, including an inert integral's limits;
   - a nested sibling removed.

   One further test: a duplicate join or name row is rejected.

   Run against the real streams, the comparator also prints per-row accounting (joined, unjoined with its reason,
   or local) and two leaf counts per joined row: one enumerated independently from each stream's parse, and the
   count actually compared. A row compared below its parsed count stays visible.

7. **Out of scope.** The ablation transcripts (each harness's knife triples are evidence about its own engine) and
   the SymPy export delta `S9b_exports.py` (SymPy-only). Neither is joined.

8. **Running.** Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode
   (`--pool s9b-comparator --memory-gib 16`). There are no time limits. Stop and report any kill or admission
   refusal; never answer one by narrowing what is compared. Write output under your own scratch directory, never
   under `scripts/out/` or `mathematica/out/`.

## Exclusions

- No oracle check beyond what the engines print: Part B's residuals against the spec's references are engine
  objects and are compared as such.
- No new physics, closure, normalization or profile choice beyond the test profiles of item 4 above. No engine edit. No
  commit.

## Builder report (at most 40 lines)

- The name table, each row with its spec object and two engine-line citations.
- Per-row accounting and leaf counts; every unjoined or local tag with its reason.
- The test profiles and points, with the spec lines that make each profile admissible.
- Which machinery you reused, and what you re-verified against these streams.
- The commands for, and outcomes of, the tests and the full guarded run.
- Runtime and peak memory.
- That no residual target was given.
