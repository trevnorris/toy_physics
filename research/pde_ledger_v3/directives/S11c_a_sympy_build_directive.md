# S11c-a SymPy engine — build directive

## Authority and boundary

Write `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` in full. Its products are the
flushed stdout tag stream and `research/pde_ledger_v3/scripts/S11c_a_exports.py`. Those two products are its
only writes; every other file, including `S11b_exports.py`, is read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (HEAD `2926c71c`) is the sole
physics authority and the sole physics input: every equation, ansatz, premise, supplied law, face map,
derivation route, control and tag obligation enters from it alone, and it wins every conflict with this
directive. Implement its §§0–8 directly. Point to rather than duplicate: §1 (the inherited S11b setup); §2
(the background ansatz, the two anchoring branches, the density representatives, admissibility premise); §3
(the face maps, the supplied laws, the shifted-trace and dynamic-window laws); §4 (the objects to compute,
T-0…T-i); §5 (the routes and controls); §6 (method, dimensions, the three script clauses, no VERDICT); §7
(names, F9, tag grammar); §8 (supplied vs computed). Add no expected value or acceptance criterion
(`CLAUDE.md` rule 5) — S11c-a withholds nothing numeric, and this §-wiring adds no physics.

Do not add parallel machinery: no per-cell completion registry beyond §4/§5's task structure, and no
directive-local exit policy. §6's script clauses and the emission-not-conditional-on-value rule apply as
stated, without named exceptions.

## ⭐ The one convention this directive pins (both engines must anchor it identically)

Two final review legs flagged that §4's T-i types the un-superscripted `μ_s=μ_θ/ρ_br,bg⁰` while §3a's
constitutive map is the branched `μ_s^α ≡ μ_θ^α/ρ_br,bg^{0,α}`. ⭐ **Use §3a's branched form.** `μ_θ` is a
**reserved S11c-b operand** (§4 T-i: keep the variable-coefficient `μ_θ` operator as S11c-b's named operand,
⛔ do not construct it), and on the MATERIAL_ADVECTED branch it is anchored at `χ` exactly as every other
`Q_bg` in §2c — i.e. carry it as `μ_θ^α` (the single branch-anchored reserved symbol). ⛔ Do not rebuild
`μ_s` from T-i's un-superscripted glyph, and ⛔ do not construct `μ_θ` from the uniform energy. This is a
naming/anchoring convention, ⛔ not new physics.

## Accumulated export

Import `LEDGER` from `research/pde_ledger_v3/scripts/S11b_exports.py` as the incoming chain record, and write
`S11c_a_exports.py` as the accumulated model-definition and knob record the next sub-step imports — the
previous `LEDGER`, plus this step's own entries, merged flat.

Bind the inherited objects to the imported rows per §1/§2a of the spec, ⛔ never re-declaring them under a
new identity: `c_s0` → `LEDGER['c_s0']`, `μ_R` → `LEDGER['mu_R']`, `ρ_br⁰`/`W̄₀`/`e_W` →
`LEDGER['rho_br']`/`LEDGER['W_0']`/`LEDGER['e_W']` (§2a: `W̄₀ ≡ W_0`, `μ̄_R ≡ mu_R` are the imported
constants), `ρ_m` → `LEDGER['rho_m']` (the bulk density is an S11b `KNOB`, `S11b_exports.py:7292`; ⛔ do
not re-originate it — a fresh declaration severs chain identity in the flux/balance/traction/closure objects
that use it), and `v_dr` → `LEDGER['v_bulk_normal_0']`. Only the **varying profiles** of §2/§7 (`W_bg`,
`mu_R_bg`, the density representatives, `e_W_bg`, …) originate in this step and have no upstream row; ⛔ do
not manufacture one, and ⛔ do not let a varying profile reuse a constant key (§7's reservations bind).

This step's candidate population is the computed objects the **primary** tasks of **§4 (T-0 … T-i)** emit and
the free symbols they contain (`S9_REWRITE_PLAN.md` **D1**: export what the primary derivations emit).
⛔ **The §5 controls (5a representation-invariance + one-sided corruption, 5b form ablation, 5c uniform-limit)
and §6's homogeneity controls are ablation evidence, ⛔ NOT exports** — they are the S11b-`B8` analog. §7's
engine-local tags are engine-local, ⛔ not exports. Carry forward the imported `LEDGER`. A row's `class`
follows from what its object **is** — a declared symbol carries its declared class, ⛔ never a table lookup by
row name. Use the class vocabulary at `research/pde_ledger_v3/S9_REWRITE_PLAN.md:41` and
`research/pde_ledger_v3/REBUILD_HANDOFF.md:40`; report anything unclassifiable under §7/§8 rather than guess.

`research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` **F1** (keys stay FLAT; `D5` unchanged)
is the storage rule, and **F9** alone decides which key form a collision writes. **Apply both whole; ⛔ do not
restate them, specify a tag-to-key transform, or substitute an earlier export proposal** — §7 already directs
the engine to apply the inherited F9 collision check. ⭐ **F9's collision rule is three-valued and its
comparison is TOTAL over the imported row shapes** (elementwise into containers; a residual that raises on a
tuple/boolean/string/relation makes `F9b` silently never fire, `S11_export_chain_decisions_v2.md:171-186`):
apply it whole. Emit as computed objects each candidate row's engine-derived key and the imported `LEDGER`
keys having those spellings (rule-2 operands, ⛔ not a claim or count; add no guard). Apply **F3** (a
re-derived row carries its own evidence, in the row). Two S11c-a-specific points pointing cannot supply:
**(1)** on an `F9c` write the prefix is the step tag **`s11c_a_`** — `s11_` and `s11b_` are already taken in
the imported `LEDGER` (`S11_export_chain_decisions_v2.md:203` left the generalization open); **(2)** S11c-a
has **no locus protocol** (it emits geometry, not spectral loci), so `F9`'s proved-different-vs-undecided
reporting uses the **inherited `T5` default**, with no spec-specific override. ⛔ Note any `F9c` write in the
§8 report.

`BUILD_INPUT_DIGESTS` is a `MappingProxyType` sha256-pinning **exactly** the three inputs the chain
requires: `{ 'S11c_a_interface_geometry_sympy_audit.py' (this engine's own source), 'S11b_exports.py',
'S11c_a_SHARED_PHYSICS.md' }`. Freeze the export as S11b does — outer `MappingProxyType`, per-record
`MappingProxyType`, then `del _LEDGER`.

Choose **F6**'s first branch: publish `S11c_a_exports.py` only if every **primary** task of §4 (T-0 … T-i)
completed. Use §4/§5's observed run record; add no completeness field.

## Three deviations from copying S11 verbatim (inherited from S11b `G8`)

- **(a) The cross-engine comparator is a SEPARATE downstream artifact — ⛔ NOT built here.** It is authored
  after both S11c-a engines build, to the frozen **`T7`** contract
  (`S11_C17_C18_spec_repair_decisions_v2.md:53-60`) and S11b's **`G8(a)`**
  (`S11b_unified_decisions.md:82-87`) — ⭐ applied **whole**, ⛔ **not restated here** (a paraphrase drops
  `G8(a)`'s **repoint ablation** and weakens `T7`'s undecided from a `T5` coverage-bounding finding to an
  ordinary third state — the same defect class as an F9 paraphrase). This engine emits no comparator and no
  cross-engine residual.
- **(b) Restore the `D3` in-run exec-and-residual round-trip** of the finished export module
  (`S9_REWRITE_PLAN.md:217`: reconstruct what was written and compare against the live object, with a
  human-readable rendering alongside).
- **(c) Include the `_RELATIONALS` reviver block** in `S11c_a_exports.py` (S11b carries it): if any emitted
  row's `value` is a `sympy` relational, `_restore` must revive it. ⚠ S11c-a's primary objects are
  shape-derivative expressions; if none is relational, the block is present and simply unused — ⛔ do not
  drop it and ⛔ do not manufacture a relational to justify it.

## Mechanical precedent

S11b's exporter (`scripts/S11b_interface_coupling_law_sympy_audit.py`'s export-emission and
`scripts/S11b_exports.py`) is **mechanical precedent, not authority** — the row schema (`display`, `value`
as `_restore("<srepr>")`, `value_kind`, `class`, `step`, optional `dimension_key`, `f9_operands`,
`corroborated_steps`), the `_restore`/`_RELATIONALS` reviver, the self-referential digest, and the freeze
footer are the working shape. ⛔ Non-verbatim: any S10/S11/S11b zero-collision-residual assertion is not
reused — `F9b` carries **no** residual (the operand pair is the object).

## Conflicts

`F9` is not an acceptance harness; ⛔ do not encode any review checklist as guards, registries or invariants
in the engine, and ⛔ do not make a check audit the input from which it was built. Use the spec's §7/§8 (and
the §8 report) for any conflict, ambiguity, unavailable construction, or object you cannot emit under
one-tag-per-object. ⛔ Do not fill such a gap with new physics, an expected result, or a self-review
mechanism. ⛔ Do not run the full task set as a completeness loop that hides a per-task failure — §6's exit
rule (emit-and-continue on a physics disagreement; nonzero only on operational failure) is exact.
