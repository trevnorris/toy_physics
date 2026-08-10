# S11 PY engine rebuild — decision list

**Deliverable:** `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`, a **full rewrite**,
which on running writes `research/pde_ledger_v3/scripts/S11_exports.py`.

**The physics is `directives/S11_SHARED_PHYSICS.md` at `ab8cb50e` (1005 lines).** ⛔ This list adds no
physics, names no object the spec does not name, and states nothing about what any computation produces.

**The pattern is `S9_REWRITE_PLAN.md`**, whose `D1`–`D11` are settled. ⛔ A builder does not re-choose them.
⭐ In particular `D3` (reconstruction format, **round trip verified in the run**), `D4` (dimension written
in the shape the engine's own solve produced), `D5` (**flat** `LEDGER`, `{value, dim, class, step}`).

⛔ The existing file is replaced, ⛔ not repaired: it fails at import (`registry_read`, `:21`), it ends in a
`VERDICT` tag the spec deletes, and it predates `§8`'s one-tag-per-named-object rule.

---

## 1 · What gets exported

⭐ **Every object the `MAIN` package emits, at every `D` in its sweep.** ⛔ The seven control packages are
ablation **evidence**, ⛔ not exports.

⚠ This is `S9_REWRITE_PLAN.md#D1` carried forward, and it is the measured behaviour of the previous step:
`S10_exports.py` carries **617** entries and **zero** keyed to a control package. ⛔ Do not widen it.
⭐ `D1`'s stated reason for using package membership was that S9 had no `_LOCAL_` convention. S11 has one,
and it does not change this boundary: `_LOCAL_` is an **engine-parity** device, ⛔ not an export rule.

## 2 · Engine-local objects

⭐ A `_LOCAL_` tag's payload **is exported when it is one of `MAIN`'s derived objects**, and ⛔ is not
exported when it exists only to compare against an imported operand.
⇒ ⭐ `Q6r`'s **derived** dimension vectors export (they are `MAIN`'s own dimension solve); ⛔ `Q6r`'s
imported vectors, differences, provenance, `Q6R_RESOLVED_COEFFICIENTS`, `Q6R_UNRESOLVED_COEFFICIENTS` and
`Q6R_RESIDUAL_SCOPE` do **not** — they are a comparison against the previous step, and exporting them would
put the previous step's own values back into the chain as if this step had derived them.

## 3 · How a tag name becomes a `LEDGER` key

⛔ The `LEDGER` key is **not** the tag. `§8`'s tag is `PY_S11_<PACKAGE>_D<n>_<QUANTITY>`.

⭐ The key is `<QUANTITY>`, lowercased, with the `D<n>` scope appended as `_d<n>` — the convention the
previous step already uses (`lagrangian_d3`, `el_residual_d3`; **556** of 617 keys carry it).
⭐ Objects with no `D` scope carry no suffix. ⭐ `ROOT<r>`-scoped and `STRATUM<s>`-scoped objects keep their
index in the key, in the same position the tag has it.

⚠⚠ **This is a MECHANICAL transform of the emitted tag, ⛔ not a judgement about what an object "really
is."** ⭐ Write it as one function and route every key through it. ⛔ Do not hand-name any key, and ⛔ do not
build a rename table — S10's D12 naming pass is still open and a hand-authored map here would prejudge it.

## 4 · What is imported, and what may be overwritten

⭐ Import `S10_exports.LEDGER`, bind what this step consumes, and export the **merged** dict (`D5`).

⭐ An entry this step re-derives is **overwritten in place**, and the overwrite is annotated at the line
where it happens (`S9_REWRITE_PLAN.md#4`). ⛔ Nothing is deleted from the imported `LEDGER`, and ⛔ no entry
this step does not touch may change — that invariant is the chain's integrity check.

⛔ **Do not fence, filter or scope the import.** ⚠ The downstream step seeing the upstream one is the point:
it is what makes cross-step consistency measurable. ⛔ Do not add any control to keep this engine from
seeing S10.

## 5 · Guards the run must carry

⭐ Each is a **computed** check that prints its operands, ⛔ never an assertion whose only outcome is `0`:

- **round trip** (`D3`) — reconstruct each written object and compare against the live one;
- **chain integrity** — every imported key this step did not touch is identical to what was imported;
- **class coverage** — every declared symbol carries a class tag; an untagged symbol is a finding.

⛔ **No `VERDICT`, no `PASS`, no `FAIL`, no summary judgement** (spec `§5`). ⛔ A physics finding exits `0`;
non-zero is for operational failure only.

---

## What this list does not decide

- ⛔ Anything the spec decides. Where this list and the spec appear to conflict, **the spec wins** and the
  conflict is reported rather than resolved by the builder.
- ⛔ The Wolfram engine, the comparator, the step record, or any register.
- ⛔ `S10`'s open `D12` naming pass, and `C17`/`C18`, which stay open.
