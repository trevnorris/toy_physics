# S11c-c2 EXPORT REPAIR — build directive (PUBLICATION-ONLY; ⛔ no physics change)

You are the builder (`gpt-6-astra`). This is a **publication-only** repair of one already-reviewed script's
**export serialization**. The physics is frozen and was independently reviewed clear
(`research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md`, commit `8f3a017f`): ⛔ **do not
change any computed/emitted object.** Your job is to make the ledger export carry **only what a later step binds**,
in a **compact, transparent** form — nothing more. **Your job ENDS at build → verify → report.**

## The problem
`scripts/S11c_c2_exports.py` is **60 MB**. Three causes, all in the export path only:
1. It exports `s11cc2SelfEnergyIncrement` (~16 MB), which is a **comparator/emit representation**, not a
   downstream binding — it must be EMIT-only (`.out`), not in the ledger.
2. Each exported row stores its value **fully `sp.expand`'d** (verbose srepr).
3. Each row also stores a `display` field = `sp.sstr(value)` — a full human-readable **duplicate** of each
   multi-MB operator that no consumer binds.

## Scope — edit ONLY these three sites
- `EXPORT_ROOTS` (line ~48).
- the `export_key` mapping inside `run()` (line ~894, the dict `{'CLOSED_SLAB_OPERATOR': …, 'CLOSED_COUPLING_KERNEL': …, 'SELF_ENERGY_INCREMENT': …}`).
- `publish(...)` (lines ~807–853).

You MAY add a new **publication helper** (a compaction/tree-walk function) *adjacent to* `publish` and call it from
`publish` — that is still publication. ⛔ **Do NOT modify** `build_case`, `build_face`, `extract`, `kernel_bridge`,
`retained_shape`, `traction_pairing`, `control`, `emit`, `grade_object`, or the physics loop (lines ~856–951)
except the one `export_key` dict. ⛔ Do NOT change any emitted object, tag, or its computation. **If the repair
appears to require touching construction code (`emit`/`build_case`/`extract`/`retained_shape`/…), STOP and report —
that means the change crossed out of publication and needs re-review, not a silent edit.**

## Requirements — what must be TRUE (⛔ not a prescription of how)

**R1 · Membership.** The written delta in `scripts/S11c_c2_exports.py` must be **exactly the bind-closure** of the
two roots `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` (each over all four `(anchoring, density)`
cases) plus their new-symbol declaration closure. `s11cc2SelfEnergyIncrement` must be **ABSENT** from the delta
(drop it from `EXPORT_ROOTS` **and** the `export_key` map) — it stays EMIT-only, still printed to `.out` at its
existing emit (line ~890). ⛔ No term_origins / parity / §3d / §5-control / increment-operand rows in the delta.
The existing `check_consumer` / `assert_delta_is_minimal` guards must still pass on the two-parent fold. ⚠ S11c-d
has **no concrete `IMPORT_KEYS` manifest yet**; both operators are the declared prospective binds (the closed slab
operator for d's closed spectral equation / resolvent, the coupling kernel for its Born/linear-mixing), so this is
D1's accepted "export the step's declared deliverable" — ⛔ neither drop one, nor let d re-close/re-extract at
consume time.

**R2 · Representation (transparent-compact, ⛔ not expanded, ⛔ not opaque, ⛔ not CSE).** ⚠ Each export root is
**one row holding all four `(anchoring, density)` cases** as a nested payload tree — `cas({case: {VALUE, MULTIGRADE,
DIMENSION_L_T_M, COMPUTED_BRANCH_BINDINGS, FOURIER_PROFILE_BINDINGS}})` — and a downstream `load_model` consumer
parses it exactly as this engine parses its parents: `fold[key]['value']` → `cases(v)` → `named(v,'VALUE')`
(script `cases`@84, `named`@76). So:
- **Compact ONLY the per-case `VALUE` operator leaves.** ⛔ Do **not** replace the row with a bare factored
  operator or any other shape — the cased `{VALUE,MULTIGRADE,DIMENSION,bindings}` tree MUST survive intact or the
  operator becomes unbindable. The non-`VALUE` metadata is small (grade tuple / dimension / Piecewise map / a few
  Integral identities); leave it unchanged.
- **Allowed transforms:** `sp.factor` / `sp.cancel` / `sp.together` / `collect` (or equivalent) that leave **one
  evaluable `sp.Basic` of the same type as the emitted leaf** (a `Matrix`/`Tuple` leaf stays that type, its
  entries compacted). ⛔ **Not `sp.expand`** (the current bloat). ⛔ **Not `sp.cse`** — its dummy replacement
  symbols are not declarations (absent from `NEW_DIMENSIONS`), so `check_consumer` treats them as free and a later
  `diff` differentiates an incomplete skeleton; and srepr does not preserve DAG sharing, so CSE only shrinks the
  file if the replacement table is stored, which is the forbidden opaque wrapper. ⛔ **No** `UnevaluatedExpr`,
  pickle, hold-wrapper, or a string that only re-parses; **no CSE/hold temporary may survive in the runtime value.**
- **Preserve the singular locus.** These operators feed d's resonances/spectrum, so the pole set matters. ⛔ Do
  **not** use a transform that cancels a denominator factor / changes the singular locus. If `cancel`/`together`
  is used, additionally verify each leaf's denominator (pole) factorization is unchanged — a leafwise value
  equality alone can accept a removable-singularity cancellation.

**Two separate guards — keep one, add one (⛔ do not merge or replace):**
1. **KEEP** the existing structural roundtrip (line ~848–851): the compact in-memory row ↔ `_restore(srepr(compact))`.
   That verifies serialization fidelity of the compact form; it stays.
2. **ADD** a separate **emitted↔compact semantic** check: compare each root's original emitted `objects` payload
   against its restored-compact form with **strict recursive container checks** — identical case-key sets, tuple
   arities, mapping keys, `Str` labels, and matrix shapes; exact equality for the non-compacted metadata; and
   **leafwise `sp.expand(decoded_leaf − emitted_leaf) == 0`** for every algebraic leaf. ⚠ Make it **Integral-aware**
   (top-level `sp.expand` leaves `Integral(0,…)` unevaluated — the exact F false-positive from the adjudication —
   so evaluate/`doit` the integrand or use `difference`'s Integral-dummy protection at the leaf). ⛔ Do **not**
   naively reuse `difference`/`tree` (script@88,100): its `zip` does **not** check tuple lengths, so a dropped
   trailing channel would silently compare as zero — add explicit length/shape checks. Print each literal result;
   **hard-stop** on any mismatch. ⚠ Run this against the **generated module** (the re-`_restore`d fold), not only
   the pre-serialization candidates.

**R3 · Hygiene — no non-binding duplicate.** The row schema must not carry a second full copy of each giant
operator. In particular the `display = sp.sstr(value)` field currently duplicates every multi-MB value; the delta
must carry the binding `value` + declarations + the equivalence evidence, ⛔ not a full human-readable dup. **Give
each operator row a SHORT bounded `display`** (e.g. the root name, or a one-line summary) rather than omitting the
field — this keeps the row schema uniform at negligible size. ⚠ This **supersedes**, for the giant operator rows
only, the retained "carry a human-readable rendering" clause (S9 D3 / bind-closure design D4); short coordinate
`display` stays as-is. ⚠ Do NOT drop a field a downstream **consumer** binds — only the non-binding review
duplicate. The `value_kind`/`class`/`step`/`route`/`dimension_key` fields and the coordinate/dimension declaration
rows stay as they are. No other field is a comparable non-binding oversized dup (`MULTIGRADE`, `DIMENSION_L_T_M`,
the branch/Fourier bindings are small generated evidence — leave them; cutting them would be an `emit` change the
fence forbids).

**R4 · Consumability + guards preserved.** After the edit: `check_consumer(combined, EXPORT_ROOTS)['closure']`
resolves; `assert_delta_is_minimal` passes; the F9 collision check passes; `BUILD_INPUT_DIGESTS` is present
(structure unchanged; it will re-hash this edited script — expected); the generated `scripts/S11c_c2_exports.py`
**imports cleanly** and `_restore` round-trips every row to a valid SymPy expression equal to the emitted object
(R2). The two closed operators must remain **bindable and differentiable** by a later `load_model` consumer.

## Build + verify (this is the whole job)
1. Make the three edits.
2. Regenerate by running the full script (it re-emits `.out` and re-runs `publish` at line ~952):
   `python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py > /tmp/S11c_c2_export_repair.out 2>&1`
   (~19 min, ~1.4 GB peak — watch RSS; ⛔ never two heavy CAS jobs at once; ⛔ no shell `timeout`/SIGKILL).
   ⚠ The run writes `_measurements/S11c_c2_sympy_guard_evidence.json` at its checkpoint calls (script@864,952) —
   that is an **expected run side-effect of the unchanged harness**, ⛔ not a builder-created review artifact.
3. Verify the DELIVERABLE, not exit status: `scripts/S11c_c2_exports.py` exists, imports, and
   - its byte-size is **dramatically** below the 60 MB baseline (report the number),
   - `s11cc2SelfEnergyIncrement` is **absent**; both closed operators present for **all four** cases,
   - the KEPT structural roundtrip and the ADDED emitted↔compact semantic check (R2) both printed clean (every
     leaf `== 0`, all container shapes matched),
   - `check_consumer` closure + `assert_delta_is_minimal` pass.
4. Confirm the physics `.out` **tag/emit-key names** are unchanged (identity follows from not touching
   construction) — ⛔ this means tag/key identity, **not** a byte-diff of the 499 MB payload.

## Report (build → verify → report — then STOP)
Write `/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md`
with: the exact diff sites; before/after export byte-size + per-object per-case value sizes; the literal
equivalence-check results; the closure/minimality/collision guard outputs; explicit confirmation of R1–R4; and a
one-line statement of what representation transform you used and how much it shrank each object.

⛔⛔ **Do NOT** launch any review, comparator, or downstream step; ⛔ do NOT read `.claude/skills/*`; ⛔ do NOT create
self-review / derived-or-declared / output-checker / finalize files. The independent review is the orchestrator's
job and runs after you report. A builder-facing directive ENDS at build → verify → report.
