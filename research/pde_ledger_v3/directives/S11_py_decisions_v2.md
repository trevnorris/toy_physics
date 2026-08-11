# S11 PY engine — decision list v2

**Orchestrator-written, 2026-08-11.** ⛔ v1 (`11bf8e05`) is a **record of a blocked gate**; ⛔ do not build
from it. This list folds its five findings and the ten that `S11_export_chain_decisions_v2.md` (`4d81e9de`)
settled as `F1`–`F7`.

**Deliverable:** `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`, a **full rewrite**,
which on running writes `research/pde_ledger_v3/scripts/S11_exports.py`.

**The physics is `directives/S11_SHARED_PHYSICS.md` at `cf4a21a4` (1149 lines).** ⛔ This list adds no
physics, names no object the spec does not name, and states nothing about what any computation produces.
⭐ Where this list and the spec appear to conflict, **the spec wins and the conflict is reported**, ⛔ not
resolved by the builder.

⛔ The existing file is replaced, ⛔ not repaired: `:21` is `from registry_read import …` against a
`reduction/` directory that no longer exists.

---

## ⭐⭐ The measurements this list rests on — ⛔ re-run them; ⛔ do not take them from here

Three consecutive orchestrator lists drew heavy findings for **citing a document where the artifact should
have been measured**. Each item below was measured on the committed artifact, and the command is given so a
reviewer can refute it.

| # | measured | on |
|---|---|---|
| **M1** | The row schema of the imported ledger is `{display, value, value_kind, class, step}` on **all 617** rows, `dimension_key` on **50**, `corroborated_steps` on **3**. ⛔ **ZERO rows carry `dim`** — `S9_REWRITE_PLAN.md#D5:229-231` is a **sketch that was never built**, and this list does not follow it | `scripts/S10_exports.py`, field census |
| **M2** | The ledger has **two disjoint key populations**: **46** rows whose value is a `Symbol` whose `.name` **is** the key (`mu_R`, `rho_br`, `lambdaScale`, `D` …), and the rest, which came from a tag suffix. ⭐ **All 10 uppercase keys are in the first population; every tag-derived key is already lowercase** | same, plus `S10_brane_mode_spectrum_sympy_audit.py:1830,1847,2024` |
| **M3** | The tag-derived rule S10 actually ran is `suffix.lower()` with `_d<n>` appended when the computation dimension is an integer (`:1949-1953`); **556** of 617 keys carry `_d<n>` | same |
| **M4** | Of **97** candidate `<QUANTITY>` names in the S11 spec, **3** collide with an existing ledger key — `ROOT_ORDERING` (→ `root_ordering_d2…d5`), `DIM_COEFFICIENTS`, `DIM_SOLUTION`. ⭐ **Every one is a different object under the same name**: S10's `root_ordering_d3` is `(0, mu_R*(k1**2+k2**2+k3**2)/rho_br)`, the roots of **S10's** dynamical matrix | full namespace census, spec × `S10_exports.py` |
| **M5** | `RUN_PAIRS` / `SKIPPED_PAIRS` already exist as spec-mandated **observed** objects, with `SKIPPED = declared ∖ completed`, accumulated only after a package finishes emitting | `S11_SHARED_PHYSICS.md:1037-1039` |
| **M6** | `DIM_COEFFICIENTS`, `DIM_EQUATIONS` and `DIM_SOLUTION` are emitted by **`Q6`**, which is not engine-local; ⭐ only **`Q6r`** carries the `_LOCAL_` infix, and it emits **only** the previous-step comparison | `S11_SHARED_PHYSICS.md:498-501` vs `:553-556` |

---

## The decisions

**P1 · The export set is `MAIN`'s non-`_LOCAL_` tags, at every `D` in the sweep, plus the sweep's own
`RUN_PAIRS` / `SKIPPED_PAIRS`.**
⛔ The seven control packages are ablation **evidence**, ⛔ not exports — `S10_exports.py` carries **zero**
rows keyed to a control package, and this boundary does not widen.
⭐ ⛔ **No `_LOCAL_` tag is exported**, whatever it holds. ⚠ By **M6** this loses nothing: `MAIN`'s own
dimension solve is already exported through `Q6`'s non-local tags, and every `_LOCAL_` tag is either a
comparison against an imported operand, a declaration, or one CAS's solver conditions.
⇒ ⭐ The rule is decided by the tag, ⛔ never by a judgement about what an object "really is", so ⛔ there is
no third bucket for `PREMISE_INVENTORY`, the `_LOCAL_` name-list tag or a solver-condition tag to fall into.

**P2 · A tag-derived key is the emitted tag with the engine prefix stripped, lowercased. A
symbol-identity key is the symbol's own `.name`, verbatim.**
⭐ `PY_S11_MAIN_D3_ROOT_ORDERING` → `s11_main_d3_root_ordering`. ⭐ One function; ⛔ every key routed through
it; ⛔ no key hand-named anywhere; ⛔ no rename table.
⭐ The two populations of **M2** are kept apart exactly as the imported artifact keeps them: case survives
where the key **is** a symbol name, and `Q6r`'s case-sensitive map (`μ_R → mu_R`, `B_comp → B_comp`) reads
that population, ⛔ never the tag-derived one.

⚠⚠ **This departs from `F1`'s letter — the scope is in the key — and the departure is the load-bearing
decision in this list. The reason is `M4`, and a reviewer should attack it there:**
- ⭐ `F1`'s reason is that **two steps deriving one object must be able to meet**. By `M4` **no** S11
  tag-derived key meets an imported one as the same object; all three overlaps are **different objects
  under an under-specified name**, and a flat key would hand a consumer S11's roots under the name of
  S10's.
- ⭐ The meeting `F1` protects survives intact, in the population where it is real: **symbol-identity keys
  stay flat**, so a coefficient S11 binds writes the key S9 wrote, and `F2` compares the objects. ⚠ That
  comparison is the one that can catch `Symbol("D") ≠ Symbol("D", integer=True, positive=True)`, whose
  difference prints as something that looks like zero.
- ⭐ `D11` forbids a name that varies with **which engine** computed the object. ⛔ The engine prefix is the
  one part stripped: **WL and PY produce the same key**, because §8 already fixes the tag they share.
- ⛔ This does **not** prejudge `D12`. `D12` decides what an object is **named**; ⭐ `P2` decides only that a
  name carries the scope the spec's own tag grammar already assigns it.

**P3 · The written row schema is `M1`'s, exactly** — `{display, value, value_kind, class, step}` on every
row; `dimension_key` on a row for which this engine has a dimension row; `corroborated_steps` only where
`F3` applies. ⛔ No `dim` field, ⛔ no new field invented, ⛔ no field dropped from an imported row.
⭐ `class` comes from the declaration annotation, read live; ⛔ never re-typed at the export site.
⭐ Reconstruction is verified by a **round trip in the run** (`D3`): re-read what was written and emit both
operands and the residual.

**P4 · Before writing any key, the writer compares against the imported `LEDGER` and emits what it found —
including when it found nothing.**
- ⭐ Key absent ⇒ a new row.
- ⭐ Key present, **same object** ⇒ a re-derivation: emit **both operands and the residual**, then guard.
- ⛔ Key present, **different object** ⇒ a **finding that fails loudly**, naming both. ⛔ Never a silent
  overwrite. ⛔ The builder does not resolve it by renaming.
- ⭐⭐ **Same object** means: reconstruct both to CAS objects and compare **as objects** — for a `Symbol`,
  identity **including its assumptions**; otherwise a residual that reduces to zero. ⛔ Not string equality
  of `display`, and ⛔ not equality of `srepr` alone.

⚠⚠ **Emit the number of keys this comparison examined and the number that collided, unconditionally, as
operands.** ⭐ `P2` makes a tag-derived collision impossible by construction, so ⛔ **a guard that reports
only "no collision" cannot be distinguished from a guard that never ran.** ⇒ its own vacuity must be
visible, exactly as `Q6d` requires of the homogeneity check.

**P5 · Chain integrity is measured from the two dicts, ⛔ never from the engine's record of what it
touched.**
⭐ Take a snapshot of the imported `LEDGER` at import as an **independent reconstruction**, ⛔ not a
reference to the imported mapping. Then emit, as three operands:
- the set of imported keys whose value **differs** between snapshot and output — computed by comparing the
  two dicts;
- the set of keys this step **declared** it re-derives;
- their **symmetric difference**, as the residual.

⚠ v1's guard compared against the engine's own record of what it touched, so an unintended write was
classified as touched and passed — demonstrated, `accidental_change=2->999, guard_using_actual_touched=0`.
⭐ Under `P5` an unintended write appears in the first set and not the second. ⛔ Nothing is deleted from the
imported `LEDGER`.

**P6 · An export is published only if `SKIPPED_PAIRS` contains no `MAIN` cell; and the observed
`RUN_PAIRS` / `SKIPPED_PAIRS` are written into the ledger either way.**
⭐ Both halves, ⛔ not one: the refusal stops a later step binding an object S11 never derived while every
row still reads as valid, and the rows let any consumer see the sweep that produced the file.
⭐ Both are read from the **observed** objects of `M5`; ⛔ neither is re-derived at the export site, and ⛔ a
declared sweep is not a run one.

**P7 · Import `S10_exports.LEDGER`, bind what this step consumes, export the **merged** dict.**
⛔ Do not fence, filter or scope the import. ⚠ The downstream step seeing the upstream one is the point.
⛔ Do not add any control to keep this engine from seeing S10.
⭐ A `Q6r` lookup that cannot resolve is reported as **unresolved, generically** (`F7`) — ⛔ no placeholder
vector, ⛔ no exception, ⛔ no membership of either object stated anywhere the builder can read it.

**P8 · Every declared symbol carries a class tag and an English description; an untagged symbol is a
finding the run prints.**
⛔ No `VERDICT`, no `PASS`, no `FAIL`, no summary judgement (spec §5). ⛔ A physics finding exits `0`;
non-zero is for operational failure only.

**P9 · ⭐⭐ REPORTING IS SUCCESS.** ⭐ If an item here is ambiguous, under-determined, ill-posed,
tautological, or cannot be built from what the spec supplies, ⭐ **report it and build the rest.** ⛔ Do not
invent a mechanism to cover the gap, and ⛔ do not resolve a conflict with the spec yourself.
⚠ The last builder given this instruction refused three items and was right every time.

---

## ⛔ What this list does not decide

- Anything the spec decides.
- The Wolfram engine, the comparator contract, the step record, the paper card, any register.
- `S10`'s open `D12` naming pass and its `C19` disclosure; `F4`'s regeneration of S10's export under `F3`,
  which is S10's debt and ⛔ not a precondition of this build.
- `F3`'s row shape beyond `P3`: ⚠ by `M4` and `P2` this step writes **no** re-derived tag-derived row, so
  ⛔ do not build a field for one. ⭐ A symbol-identity re-derivation is covered by `P4`'s operands.
