# Script rebuild — state and what's next

Read `CLAUDE.md` first (how we work), then this (where we are).
Process detail: `docs/development_pipeline.md`. Defects: `DEFECT_REGISTER.md`.

Last updated 2026-08-06, after S11's shared spec closed.

---

## The scope — three steps need rebuilding, and the debt is exactly these

The measured defect was in **three independently-built steps**: engines emitting physics conclusions as
typed sentences with no CAS object behind them, missed by eight review legs.

| step | state |
|---|---|
| S9 | rebuilt |
| S10 | **CLOSED** `e5a2c695` — harness, record, register pass 1, spec repair, card |
| S11 `stray_longitudinal` | ▶ **IN PROGRESS** — spec done `f49a1684`; engines not yet rebuilt |
| S11b-A `interface_response` | built under the old pattern → **rebuild** |
| S11b-B `interface_assembly` | built under the old pattern → **rebuild** |
| S11b-C non-uniform coupling | **never built** — the MacCullagh differentiator |

S12 onward does not exist yet, so every remaining step in every remaining sector gets built under the new
pattern the first time. This tax is paid once.

### The pattern each rebuilt step must follow

1. Both engines emit **CAS objects** as tagged payloads. No typed conclusions; a residual asserted zero
   always prints `0` and carries no information — emit both operands and the residual, then guard.
2. A **`reduction/checks_S<n>.yaml`** declaring cross-engine pairs, the control matrix, dimension sources
   per engine and package, and action families. This is the mechanical tie-in to the harness.
3. The harness runs and the **step record cites what it measured**, not what the author concluded.
4. **Requirements-register entries** fall out of the same read ⇒ `SUBSTRATE_REQUIREMENTS.md`.

S10 is the worked example of all four; its record is the template.

---

## S11 — where it actually is

### ⭐ Done

| | |
|---|---|
| **Shared physics spec** | `directives/S11_SHARED_PHYSICS.md`, committed **`f49a1684`**. Three review rounds, four legs. |
| **As-built baseline** | tag **`s11-as-built`** (annotated, at `e7658d3a`) — both old engines and the two directives they were built from. |

### ▶ In flight

**The comparator's derivative atoms.** `OpaqueDerivative` and `CanonicalDerivative` both carried a
process-global identity defect. Uncommitted in the working tree; four fix rounds, six review legs.
⛔ **This must land before `checks_S11.yaml` declares any rename.**

### ⛔ THE SURVEY CORRECTED THIS FILE'S OWN PREVIOUS CLAIMS — verify, do not inherit

The earlier version of this handoff was wrong on three of four points. Measured 2026-08-06:

- ⛔ **"The two engine directives are gitignored and not in the repository."** ⭐ **False.** Both are
  **tracked in HEAD** (`911d0af8`); gitignore never applies to already-tracked files. Nothing needed
  promoting.
- ⛔ **"The namespace is misaligned — WL emits `WL_S11_*`, PY emits `S11_*`."** ⭐ **Understated.** The two
  engines have different tag **vocabularies** and different **granularity** — WL bundles
  `{root → {nullity, orientation}}` into one payload where PY emits three. Stripping engine prefixes, the
  two share exactly **one** tag suffix: `VERDICT`, which the new spec deletes. ⇒ renaming prefixes yields a
  machine-comparable row count of **one**.
- ⛔ **"The `.wl` prints a boolean as `TRUE`/`FALSE` — a rule-2 slip."** ⭐ **It is specified by its own
  directive** (T4, T7, T8, T9, T10). Repairing the engine's line 9 would have fixed nothing.
- ✅ **"No committed outputs, no `checks_S11.yaml`."** Both correct.

⇒ ⭐ **S11 is a RECONSTRUCTION, not a repair.** Emission volume, for scale: S11's engines emit **23** and
**79** payloads; S10's emit **2983** and **4227** against a 690-row declaration.

### ⭐ The remaining order — do not reorder

1. Land the comparator fix (two legs, then commit).
2. **Build both engines**, independently, each from the spec plus a short engine-specific directive.
   Two legs per engine.
3. Run both into committed outputs under `mathematica/out/` and `scripts/out/`.
4. Write **`reduction/checks_S11.yaml`**. Two legs.
5. **Extend `harness_ablation.py` to cover it.** ⛔ A checks config with no battery is an unpoliced
   declaration — S10's battery is `ACCEPTANCE 1–19` and it is what let a leg find two gate holes.
6. Run the harness. Only now can it say anything about S11.
7. Step record · paper card · requirements-register pass 2 · defect-register entries.
8. **Repoint the registry provenance** for `B_comp`, `c_L` and `R5`.

### ⚠ What the next session must NOT be told, and must run itself

`W_XFORM_CURLONLY` is **literally S10's `MAIN` action** — same density, ansatz, phase average and `D`
sweep. ⇒ S11's curl-only package should reproduce S10's baseline exactly: the **same physical system
computed by four independently built engines**, which the ledger has never had.

⛔ The cross-step language was **deliberately removed from the spec** during the loosening round, because
pointing a builder at committed S10 rows is a target it can converge on. ⭐ **The comparison is the
orchestrator's to run after both engines exist.**

### ⚠ One trap for `checks_S11.yaml`

⛔ **Do not declare `N6_BASIS` as a cross-engine row.** A nullspace basis is not canonical; S10 has 11 rows
reporting `DISAGREE` on representation alone. The spec already says `N6` is display-only.

---

## ⭐ Banked for the S11 step record — findings, not fixes

- ⛔ **`V7_RESIDUAL = 0` does not validate `V1`.** A leg's own flawed route returned `V1_DIM` 2 against 4,
  and `V2`, `V6` and `V7` were all self-consistent on the wrong `V1`. `V7` tests `V2` against `V6`
  **within** `V1`; only the cross-engine `V1_BASIS` comparison tests `V1` itself.
- ⛔ **`Q.brane.B_comp` and `Q.brane.c_L` carry `source_locus` pointing at the engine being replaced**
  (`scripts/S11_stray_longitudinal_sympy_audit.py:632-664`), and `R5` came from that build. ⇒ the registry
  is **not an independent operand** for those two rows — only for `ρ_br` and `μ_R`.
- ⛔ **The old build's `N_O = 3 for every D` is a SINGLE-ENGINE claim**: the `.wl` has zero reflection /
  `O(D)` content, inside a record section headed *"two engines… agreeing on every computed value"*.
- ⚠ **The S11 pre-registration was wrong on the invariant count and the engines were right.** Not a clean
  pass; say so.
- ⛔ **A `CanonicalDerivative` collision was a false-AGREE channel in the atom every verdict runs on.**
  Found while dormant: 2072 constructions on S10 across 162 names, zero carrying two identities.

## ⛔⛔ Two process findings that cost this session, both now governing

- **A test that "covers" an invariant demonstrates it on one example.** Three times a half-fix passed the
  new test, the whole suite, the full battery **and** produced byte-identical comparator output.
  ⇒ require the builder to build the weaker implementations and show the test fails on each
  ⇒ `feedback_test_must_fail_on_weaker_fixes`.
- **An unchanged regression baseline can be worth nothing.** Measured: `unfixed`, `fixed` and a
  knowingly-wrong variant all produced identical stdout on both configs, because no separator character
  occurs in any identity field. ⛔ Never report "the baseline did not move" as evidence of correctness.

---

## ⚠ What is still open on S10, deliberately

- **`declaration_load_ablation.py` runs ~19 minutes on the S10 config.** ⚠ Under the clarified rule this is
  **acceptable if it prints incrementally** — the bar is observable progress, not elapsed time
  ⇒ `feedback_script_timeout_policy`. Check that it flushes; ⛔ do not narrow it.
- **Three sentences corrected after the last leg reported** — the eight-of-twelve verdict count and its two
  companions. They have not been through a leg.
- `DEFECT_REGISTER.md#f7`'s owed measurement. Re-scope it or discharge it.
- ⚠ **A second registry provenance dent**: `D_brane`'s `source_locus` (`reduction/quantities.yaml:177`)
  points at S9 for a quantity S10 introduces. Fix with the `B_comp`/`c_L` pass.

## Known limits of the instrument

- The canonical derivative key uses `sorted(set(arguments))`, so a repeated coordinate in a dependence list
  is not distinguished. Multiplicity is not physical; recorded rather than open.
- A stiffness form the gradient substitution cannot cover has no completeness guard.
- Exit code carries no signal on either config; both exit 2 in the healthy state. **Judge by counters.**
- `_free_symbol_names` has **no call site anywhere in the repository** — pre-existing dead code.
- ⚠ The comparator requires `--config` **and** one `--output ENGINE=PATH` per engine. It resolves the
  registry root from the config path, so a bare copy of `reduction/` fails with `RegistryValidationError`.

## Pinned

- **`s11-as-built`** (`e7658d3a`) — both S11 engines and the two directives they were built from, with the
  survey findings in the tag annotation.
- **`s10-as-built`** (`9309da70`) — the S10 spec and both engines as they ran before that rebuild.
- **`wip-2026-08-05-unreviewed`** (`92461853`) — two builds committed without review and reverted. Prior art
  only; nothing in it is known-good.
