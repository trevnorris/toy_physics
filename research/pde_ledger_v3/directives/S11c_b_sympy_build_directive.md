# S11c-b SymPy engine — build directive

## Authority and boundary

Write `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` in full. Its products are the
flushed stdout tag stream and `research/pde_ledger_v3/scripts/S11c_b_exports.py`. Those two products are its
only writes; every other file, including `S11c_a_exports.py` and `S11b_exports.py`, is read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (HEAD `7023420d`) is the sole
physics authority and the sole physics input: every equation, inherited object, ansatz, premise, supplied
law, derivation route, control and tag obligation enters from it alone, and it wins every conflict with this
directive. Implement its §§0–8 directly. Point to rather than duplicate: §1 (the inherited S11b/S11c-a setup,
the DOF and the local curl/div sector labels, the stored energy and symmetry inputs, the uniform decoupling
context); §2 (the background ansatz inherited from S11c-a, which quantities vary vs stay constant, and the
activated admissibility computation); §3 (the variable-coefficient energy basis and new invariants, the
divergence-form operator, the off-diagonal coupling kernel, and the §3d background-order balance); §4 (the
objects to compute and emit); §5 (the representation-invariance / one-sided independence / form-ablation /
uniform-limit controls); §6 (method, dimensions, the three script clauses, no VERDICT); §7 (names, F9, tag
grammar); §8 (supplied vs computed). Add no expected value or acceptance criterion (`CLAUDE.md` rule 5) —
S11c-b withholds every computed value (the coupling grade/sign, the new invariants, the admissibility
residual, the basis count), and this §-wiring adds no physics.

Do not add parallel machinery: no per-cell completion registry beyond §4/§5's task structure, and no
directive-local exit policy. §6's script clauses and the emission-not-conditional-on-value rule apply as
stated, without named exceptions.

## ⭐ The conventions this directive pins (both engines must anchor them identically)

- **The sector labels are the §1a local differential operators, not a global projection.** `u_T` is the
  `∇×u` part, `u_L` the `∇·u` part; ⛔ do not build a global Helmholtz projector, plane-wave projector, or
  inverse-Laplacian (spec §1a/§3c, `N5`). The off-diagonal block is extracted by the operator's own curl/div
  structure.
- **The energy basis is a non-unique quotient; construct YOUR OWN representative — do not pin one.** The spec
  (§1c/§3a) requires constructing the variable-coefficient basis from the symmetry group with the background
  first jets as spurions; different representatives are reconciled after the run, ⛔ never by importing the
  sibling's basis or a fixed representative. ⛔ Do not obtain the energy by substituting `W_0→W_bg` into the
  uniform `U` (spec §2a/§3a, `N15`): that omits the newly admitted gradient-of-background invariants, which
  are emitted as results.
- **`mu_theta` is the reserved variable-coefficient operand (spec §1c/§3b/§4).** Keep it as a named
  variational operator `mu_theta = δU/δθ` on the variable-coefficient energy; ⛔ do not construct it by
  substitution into the uniform energy. On the MATERIAL_ADVECTED branch it is anchored at `χ` exactly as
  every other `Q_bg` in the inherited §2c, i.e. carry `mu_theta^α`.
- **Admissibility is the §3d background-order (ε⁰) balance, not the ε→0 wave operator.** Compute the operator
  operand as the background-order generalized body force + per-face traction the profile sources, in the same
  ordered pairing as `𝒮_hold⁰` (spec §2b/§3d). ⛔ Do not define it as the `ε→0` limit of the §3b first-order
  operator, and ⛔ do not insert `W_bg−W_0` into the uniform perturbation equations.

These are naming/anchoring conventions, ⛔ not new physics; every one is already in the spec.

## Accumulated export

Import `LEDGER` from `research/pde_ledger_v3/scripts/S11c_a_exports.py` as the incoming chain record (it
already carries the merged S11b `LEDGER` plus the S11c-a rows — `N1`: each sub-step imports the prior one's
exports), and write `S11c_b_exports.py` as the accumulated model-definition and knob record the next sub-step
imports — the previous `LEDGER`, plus this step's own entries, merged flat.

Bind the inherited objects to the imported rows per §1/§2 of the spec, ⛔ never re-declaring them under a new
identity: the constants `c_s0` → `LEDGER['c_s0']`, `μ_R` → `LEDGER['mu_R']`, `ρ_br⁰`/`W̄₀`/`e_W` →
`LEDGER['rho_br']`/`LEDGER['W_0']`/`LEDGER['e_W']`, `ρ_m` → `LEDGER['rho_m']` (an S11b `KNOB`; ⛔ do not
re-originate it — a fresh declaration severs chain identity in every flux/balance/traction/closure object
that uses it), `v_dr` → `LEDGER['v_bulk_normal_0']`; and the S11c-a-originated varying objects the operator
consumes — the profiles `W_bg`, `mu_R_bg`, `w1_profile`, `m1_profile`, `L_W`, `sigma_W`, `eta_bg`, the
density representatives, `e_W_bg`, and the S11c-a shape-derivative substrate rows (the `S11CA_*` T-objects) —
bind to their `S11c_a_exports` rows, ⛔ never re-originating them. Only S11c-b's **own new objects** (the
variable-coefficient operator, the coupling kernel, the new gradient-of-background invariants and any new
constant they introduce, the admissibility operands) originate in this step and have no upstream row; ⛔ do
not manufacture one, and ⛔ do not let any of them reuse a constant key (§7's reservations bind: the kernel
and new invariants take fresh names, never `W_0`/`mu_R`/`rho_br`/`e_W`/`v_0`).

This step's candidate population is the computed objects the **primary** tasks of **§4** emit (the
variable-coefficient energy basis and new invariants, the divergence-form operator and its μ_θ operand, the
off-diagonal coupling kernel, and the admissibility operands) and the free symbols they contain
(`S9_REWRITE_PLAN.md` **D1**: export what the primary derivations emit). ⛔ **The §5 controls (5a
representation-invariance + one-sided corruption, 5b per-source form ablation, 5c uniform-limit) and §6's
homogeneity controls are ablation evidence, ⛔ NOT exports** — they are the S11b-`B8` analog. §7's
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
re-derived row carries its own evidence, in the row). Two S11c-b-specific points pointing cannot supply:
**(1)** on an `F9c` write the prefix is the step tag **`s11c_b_`** — `s11_`, `s11b_`, `s11c_a_` are already
taken in the imported `LEDGER`; **(2)** S11c-b has **no locus protocol** (it emits an operator and kernel,
not spectral loci — the spectrum is S11c-d), so `F9`'s proved-different-vs-undecided reporting uses the
**inherited `T5` default**, with no spec-specific override. ⛔ Note any `F9c` write in the §8 report.

`BUILD_INPUT_DIGESTS` is a `MappingProxyType` sha256-pinning **exactly** the three inputs the chain requires:
`{ 'S11c_b_brane_operator_sympy_audit.py' (this engine's own source), 'S11c_a_exports.py',
'S11c_b_SHARED_PHYSICS.md' }`. Freeze the export as S11b/S11c-a do — outer `MappingProxyType`, per-record
`MappingProxyType`, then `del _LEDGER`.

Choose **F6**'s first branch: publish `S11c_b_exports.py` only if every **primary** task of §4 completed. Use
§4/§5's observed run record; add no completeness field.

## Three deviations from copying S11 verbatim (inherited from S11b `G8`)

- **(a) The cross-engine comparator is a SEPARATE downstream artifact — ⛔ NOT built here.** It is authored
  after both S11c-b engines build, to the frozen **`T7`** contract and S11b's **`G8(a)`** — ⭐ applied
  **whole**, ⛔ **not restated here** (a paraphrase drops `G8(a)`'s **repoint ablation** and weakens `T7`'s
  undecided from a `T5` coverage-bounding finding to an ordinary third state). This engine emits no comparator
  and no cross-engine residual.
- **(b) Restore the `D3` in-run exec-and-residual round-trip** of the finished export module
  (`S9_REWRITE_PLAN.md:217`: reconstruct what was written and compare against the live object, with a
  human-readable rendering alongside).
- **(c) Include the `_RELATIONALS` reviver block** in `S11c_b_exports.py` (S11b/S11c-a carry it): if any
  emitted row's `value` is a `sympy` relational, `_restore` must revive it. ⚠ If none of S11c-b's primary
  objects is relational, the block is present and simply unused — ⛔ do not drop it and ⛔ do not manufacture
  a relational to justify it.

## Mechanical precedent

S11c-a's exporter (`scripts/S11c_a_interface_geometry_sympy_audit.py`'s export-emission and
`scripts/S11c_a_exports.py`) is **mechanical precedent, not authority** — the row schema (`display`, `value`
as `_restore("<srepr>")`, `value_kind`, `class`, `step`, optional `dimension_key`, `f9_operands`,
`corroborated_steps`), the `_restore`/`_RELATIONALS` reviver, the self-referential digest, and the freeze
footer are the working shape. ⛔ Non-verbatim: any zero-collision-residual assertion is not reused — `F9b`
carries **no** residual (the operand pair is the object). ⛔ Do not import or transcribe any S11c-a computed
physics beyond the chain bindings named above — S11c-b's own §3–§5 objects are re-derived from its spec.

## Conflicts

`F9` is not an acceptance harness; ⛔ do not encode any review checklist as guards, registries or invariants
in the engine, and ⛔ do not make a check audit the input from which it was built. Use the spec's §7/§8 (and
the §8 report) for any conflict, ambiguity, unavailable construction, or object you cannot emit under
one-tag-per-object. ⛔ Do not fill such a gap with new physics, an expected result, or a self-review
mechanism. ⛔ Do not run the full task set as a completeness loop that hides a per-task failure — §6's exit
rule (emit-and-continue on a physics disagreement; nonzero only on operational failure) is exact.
