# S11b SymPy engine — build directive

## Authority and boundary

Write `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py` in full. Its products are
the flushed stdout tag stream and `research/pde_ledger_v3/scripts/S11b_exports.py`. Those two products are
its only writes; every other file, including `S11_exports.py`, is read-only. The two historical S11b
engines (`scripts/S11bA_interface_response_sympy_audit.py`, `scripts/S11bB_interface_assembly_sympy_audit.py`)
and their stdout are **not** build premises — this is one unified step, not a merge of those two.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` is the sole physics authority
and the sole physics input: every equation, ansatz, premise, supplied law, derivation route and tag
obligation enters from it alone, and it wins every conflict with this directive. Implement its §§0–13
directly. Point to rather than duplicate §4; §5 and its basis construction; §6, its virtual-displacement
rule, causality diagnostics, convention cross-checks and energy accounting; §7; §8, including its
corollaries, the no-verdict rule and the locus protocol; §9; and §10. Add no expected value or acceptance
criterion (`CLAUDE.md` rule 5); the spec withholds every reference the orchestrator diffs (§11 wiring aside,
below, adds no physics).

Do not add parallel machinery for those contracts: no per-cell completion registry beyond §9's task
structure, and no directive-local exit policy. §8 corollary 4 applies as a property, without named
exceptions.

## Accumulated export

Import `LEDGER` from `research/pde_ledger_v3/scripts/S11_exports.py` as the incoming chain record, and write
`S11b_exports.py` as the accumulated model-definition and knob record the next step imports — the previous
`LEDGER`, plus this step's own entries, merged flat.

Bind the inherited objects to the imported rows per the spec's §11 mapping, ⛔ never re-declaring them under
a new identity: `c_s0` → `LEDGER['c_s0']`, `μ_R` → `LEDGER['mu_R']`, and `ρ_br⁰` → `LEDGER['rho_br']`
(the spec settles `ρ_br⁰ ≡ ρ_4D⁰ W₀` **is** the imported `rho_br`; ⛔ do not ask the engine to decide that
identity and ⛔ do not mint a second inertia knob). `ρ_m` and `v_dr` (`v_bulk_normal_0`) originate in this
step and have no upstream row; ⛔ do not manufacture one.

This step's candidate population is the computed objects the **primary** derivations of §9 emit — every §9
task **except B8** (the controls) — and the free symbols they contain (`S9_REWRITE_PLAN.md` **D1**: export
what the primary derivations emit;
**B8's control recomputations are ablation evidence, ⛔ not exports**; §10 `_LOCAL_` tags are engine-local,
⛔ not exports). Carry forward the imported `LEDGER`. A row's `class` follows from what its object **is** — a
declared symbol carries its declared class, ⛔ never a table lookup by row name; ⛔ do not author a
per-quantity class map. Use the class vocabulary at `research/pde_ledger_v3/S9_REWRITE_PLAN.md:41` and
`research/pde_ledger_v3/REBUILD_HANDOFF.md:40`. Report anything unclassifiable this way under §10 rather than
guess or add apparatus.

`research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` **F1** (keys stay FLAT; `D5` unchanged)
is the storage rule, and **F9** alone decides which key form a collision writes. **Apply both; ⛔ do not
restate them, specify a tag-to-key transform, or substitute an earlier export proposal.** Emit as computed
objects each candidate row's engine-derived key and the imported `LEDGER` keys having those spellings —
these are `CLAUDE.md` rule 2 operands, ⛔ not a claim or count; add no guard, threshold or assertion over
them. Apply **F3** (a re-derived row carries its own evidence, in the row) to the rows this step writes.
⭐ **`F9`'s collision rule is three-valued and its comparison is TOTAL over the imported row shapes**
(elementwise into containers — the import carries tuples, booleans, strings and relations, and a subtraction
residual that raises on them makes `F9b` silently never fire, `S11_export_chain_decisions_v2.md:171-186`):
**apply it whole and ⛔ do not re-render it here** — a paraphrase drops exactly that totality. Only two things
are S11b-specific and are stated because pointing cannot supply them: **(1)** on an `F9c` write the prefix is
the step tag **`s11b_`** — `F9` left the `s11_`→`<step>_` generalization open
(`S11_export_chain_decisions_v2.md:203`), and the literal `s11_` is **already taken** by S11's own `F9c`
keys in the imported `LEDGER` (`scripts/S11_exports.py:11123` ff., `s11_a1`…), so reusing it would collide;
**(2)** `F9`'s report distinction between proved-different and undecided reads through **this spec's §8 locus
protocol** — `F9` (`:157`) points at "§5's locus protocol", which is the physics spec's, and in S11b that
protocol is **§8, ⛔ not §5** (S11b's §5 is the energy basis). ⛔ Note any `F9c` write in the §13 report.

`BUILD_INPUT_DIGESTS` is a `MappingProxyType` sha256-pinning **exactly** the three inputs the spec names
(§11): `{ 'S11b_interface_coupling_law_sympy_audit.py' (this engine's own source), 'S11_exports.py',
'S11b_SHARED_PHYSICS.md' }`. Freeze the export as S11 does — outer `MappingProxyType`, per-record
`MappingProxyType`, then `del _LEDGER`.

Choose **F6**'s first branch: publish `S11b_exports.py` only if every primary derivation task of §9
(every task except B8, the controls) completed. Use §9's observed run record; add no completeness field.

## Three deviations from copying S11 verbatim (decision list `S11b_unified_decisions.md` G8)

- **(a) The cross-engine comparator is a SEPARATE downstream artifact — ⛔ NOT built here.** It is authored
  after both engines build, to the comparator contract of decision list **G8(a)**
  (`S11b_unified_decisions.md:83-87`): join by name, residual the paired payloads, three-valued undecided,
  and a repoint ablation — meeting the frozen **`T7`** contract
  (`S11_C17_C18_spec_repair_decisions_v2.md:53-60`, which states only: frozen before it sees either output,
  ⛔ reject a native boolean as a residual operand, undecided per `T5`). It is **net-new for the light
  sector** — S11 shipped without one. This engine emits no comparator and no cross-engine residual.
- **(b) Restore the `D3` in-run exec-and-residual round-trip** of the finished export module
  (`S9_REWRITE_PLAN.md:217`: reconstruct what was written and compare against the live object, carrying a
  human-readable rendering alongside). S10 had it; S11 dropped it for serialization symmetry. It returns
  here.
- **(c) Include the `_RELATIONALS` reviver block** in `S11b_exports.py` (S11 already carries it; S9 omitted
  it): S11b emits interface-law relations and inequalities (admissibility regions, any forced reciprocity
  relation), so a row's `value` may be a `sympy` relational, and `_restore` must revive it.

## Mechanical precedent

S11's exporter (`scripts/S11_stray_longitudinal_sympy_audit.py`'s export-emission and the resulting
`scripts/S11_exports.py`) is **mechanical precedent, not authority** — the row schema (`display`, `value`
as `_restore("<srepr>")`, `value_kind`, `class`, `step`, optional `dimension_key`, `f9_operands`,
`corroborated_steps`), the `_restore`/`_RELATIONALS` reviver, the self-referential digest, and the freeze
footer are the working shape. ⛔ Non-verbatim points: S11 **dropped `D3`** (restore it, above); S11 shipped
**no comparator** (it is separate and net-new, above); and any S10/S11 zero-collision-residual assertion is
⛔ not reused — `F9b` carries **no** residual (the operand pair is the object).

## Conflicts

`F9` is not an acceptance harness. Its file's **OWED TO THE BUILD REVIEW** checklist
(`S11_export_chain_decisions_v2.md:210`) belongs to the review legs run against the built script; ⛔ do not
encode that checklist as guards, registries or invariants in the engine, and ⛔ do not make a check audit the
input from which it was built.

Use the spec's §10 (and the §13 report) for any conflict, ambiguity, unavailable construction or object you
cannot emit under one-tag-per-object. ⛔ Do not fill such a gap with new physics, an expected result, or a
self-review mechanism.
