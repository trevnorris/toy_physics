# S9 rewrite — the export chain, starting at the beginning

**Status: PLAN, 2026-08-08.** Supersedes `EXPORT_CHAIN_ARCHITECTURE.md`, which proposed a registry-shaped
design that was reviewed twice and abandoned.

⛔ **There is no YAML in this design.** ⛔ No registry, no relations document, no per-step checks file, no
declared dimensions anywhere. ⭐ Every value and every dimension is **computed by a CAS and written by the
script that computed it.**

---

## The idea

⭐ **Each step's SymPy engine writes a module that the next step imports.** Cross-step consistency is then
true by construction — S10 does not restate what S9 found, it **uses** it.

⭐ **Each step's Wolfram engine imports nothing** and derives independently from the action and its stated
premises. That is where the contention lives: py-versus-wl **within** each step.

⇒ ⭐ consistency by dataflow on one side, independent re-derivation on the other.

## Why start at S9

It is the first step of the light sector and it consumes nothing from an earlier v3 step. ⇒ ⭐ it is the
only place the chain can begin.

---

## 1 · Annotate the declarations

⭐ Every symbol declared in the engine carries a **class tag** and an **English description**.

```python
rho_br     = sp.Symbol("rho_br", positive=True)   # KNOB · brane inertia density · not derived here
mu_R       = sp.Symbol("mu_R", positive=True)     # KNOB · brane shear modulus · not derived here
D          = sp.Symbol("D", integer=True)         # STRUCTURAL · brane spatial dimension
kx, ky, kz = ...                                  # COORDINATE · wavevector components
lambda_rho = ...                                  # CONTROL · inertial ablation scale
```

⭐ Classes: `KNOB` · `STRUCTURAL` · `COORDINATE` · `CONTROL` · `PREMISE` · `DERIVED`.

### ⛔⛔ What an annotation must NEVER contain

⛔ **What the computation will produce.** ⛔ Not a value, ⛔ not a dimension, ⛔ not a count, ⛔ not a sign.

⚠ **Measured, repeatedly, in this project:** once an expected result is written where a builder can read
it, the builder iterates until the computation matches it, and the check becomes a check that two copies of
an assumption agree. ⇒ ⭐ **the script is made correct; ⛔ the answer is read afterwards.**

⭐ **An inserted equation is allowed**, because it is an input under test rather than a result:

```python
# PREMISE · stiffness density posited as curl-only
stiffness_density = mu_R * curl(u).dot(curl(u)) / 2
```

## 2 · Emit every end derivation

⭐ The engine collects **every object it ends up deriving** — ⛔ not every intermediate, ⭐ but every result
the step actually concludes — into a named collection.

⭐ **Emit more than the next step needs.** ⚠ That is deliberate: a full, visible namespace is what makes a
naming collision obvious, and later steps — including the knit phase, which assembles the brane and bulk
requirements — import from these same files.

## 3 · Generate `S9_exports.py`

⭐ Written by the run, at the end of the run. ⛔ Never hand-edited.

⭐ It contains, for each end derivation: the object in a form SymPy reconstructs exactly, and **its dimension
as computed by the engine's own dimension solve.**

⭐ It also carries the step's **knobs** — the symbols it consumed and did not derive — with their class tags,
so a consumer knows what is still undetermined.

⭐ It must be **readable by a person.** Hand-checking it is an intended use.

## 4 · Import it into S10

⭐ S10's SymPy engine imports from `S9_exports` instead of re-declaring. ⭐ It imports what it consumes;
⭐ what it derives itself, it still derives.

### ⭐⭐ The overwrite mechanism

⭐ When a later step **derives** something an earlier step exported as a knob, the import is shadowed:

```python
from S9_exports import rho_br          # arrives as a KNOB
...
rho_br = derive_from_substructure(...)  # OVERWRITTEN · this knob is retired here
```

⇒ ⭐ the knob-to-derived transition **announces itself in the source**, at the moment it happens, and is
annotated there. ⛔ No separate register to keep in sync.

## 5 · Extract the knob inventory

⭐ A script scans the engines for annotated declarations and reports every symbol by class.
⇒ ⭐ the inventory of what is still undetermined is **generated**, ⛔ never hand-authored.
⭐ A symbol with no class tag is a finding.

---

## Verification

⭐ **The `.out` tag stream is NOT required to be byte-identical.** ⚠ The output shape is changing on purpose;
a byte bar would forbid the thing being built.

⭐ What must hold instead:

1. ⭐ **The derived objects have the same values as before**, compared **by name**. That is what shows the
   restructure moved no physics.
2. ⭐ **S10 imports `S9_exports` and runs.**
3. ⭐ **Every declared symbol carries a class tag**, and the inventory extracts.
4. ⭐ **Dimensions are computed**, and py's agree with wl's where both derive the same object.
5. ⭐ **S9's Wolfram engine still runs and its DERIVATION is unchanged** — ⭐ only emitted names moved
   (**D12**) — and its results still agree with py's where both derive the same object.

## Scope

⛔ **Do not change any `.wl` engine's derivation** — its independence is the contention and is
load-bearing. ⭐ **Emitted names may change, and must**, per **D12**.
⛔ **Do not add a YAML file, a registry, or a declared dimension.**
⛔ **Do not write an expected value or dimension anywhere a builder can read it.**

## Order

1. Annotate S9's declarations.
2. Restructure S9 so end derivations land in a named collection.
3. Generate `S9_exports.py`.
4. Import into S10 and make S10 run against it.
5. Write the knob extractor.
6. Run the verification list.

---

## ⛔⛔ THE PLAN IS CLOSED

⛔⛔ **Nothing may be added to this plan by me or by a builder.** ⚠ The failure mode of this project is not
omitting machinery — ⭐ it is **adding** it. Five build rounds and ten review legs were spent on checks that
policed other checks, and none of it touched physics.

⭐ **If something appears to be missing: STOP and ask the user.** ⛔ Do not add it. ⛔ Do not "briefly" extend
scope. ⛔ Do not build a guard for a guard.

### ⛔ Files this work may create or change — ⛔ anything else is a STOP

| may be created | may be changed |
|---|---|
| `scripts/S9_exports.py` (generated) | `scripts/S9_light_requires_shear_sympy_audit.py` |
| one knob-extractor script | `scripts/S10_brane_mode_spectrum_sympy_audit.py` |
| one baseline snapshot (committed first) | `mathematica/S9_light_requires_shear_mathematica_audit.wl` — ⛔ **emitted names ONLY**, per **D12** |

⛔ **No other file.** ⛔ No new YAML. ⛔ No new directory. ⛔ No new test framework.

### ⛔ Explicitly NOT part of this work

⛔ A provenance check that an object "traces to this step's action" — **measured unbuildable**: SymPy exposes
`free_symbols`, ⛔ not a derivation path, and an expression built from symbols is indistinguishable from one
built from the determinant.
⛔ A population/coverage guard on the export. ⛔ A test that a guard exists. ⛔ A witness. ⛔ A registry.
⛔ A relations document. ⛔ Declared dimensions. ⛔ A runner. ⛔ An ablation battery for the export writer.

⭐ **Blindness comes from one thing only: the Wolfram engine re-derives the action independently.** ⛔ Nothing
else in this design is a blindness control, and ⛔ nothing should be built pretending to be one.

---

## ⛔⛔ Deletion — ⭐ ONE PASS, ⛔ committed on its own, ⛔ before any building

⭐ **Delete `reduction/` and everything that served the reconciliation design.**

⚠ Some Python engines import `registry_read`. ⭐ **That does not matter** — ⛔ those engines are being
rewritten and ⛔ nothing will run them in their current form again. ⭐ The `.wl` engines import nothing from
`reduction/` (measured: **0** registry references), so they are unaffected.

### ⭐ What must survive

| keep | why |
|---|---|
| all `.wl` engines | untouched, siloed, and they import nothing from `reduction/` |
| `mathematica/out/*.out` | the Wolfram side of every comparison |
| **S9's current `scripts/out/*.out`** | ⭐ the **D9 baseline** — ⛔ captured and committed **first** |
| the `cross_engine` pair lists from `checks_S9.yaml` / `checks_S10.yaml` | ⭐ **D11** uses them ⇒ extract the pairs, ⛔ do not keep the files for it |

### ⭐ Also clean `out/`

⭐ `scripts/out/` and `mathematica/out/` hold outputs from the abandoned design. ⭐ Keep what the table
above names; ⛔ delete the rest.

⛔ **Everything else in `reduction/` goes**, including `quantities.yaml`, `relations.yaml`,
`registry_read.py`, `registry_schema.yaml`, `dimensional_homogeneity_gate.py`, `engine_output_checks.py`,
the witness, the pin, the `w3_*`/`w4_*` runners and every `test_*.py`.

⭐ Root-level: `_review_prompt*.md` ×6 · `HARNESS_S9_PILOT_PLAN.md` · `CROSS_STEP_DIMENSION_SCOPE.md` ·
`INTEGRATION_TODO.md` · `EXPORT_CHAIN_ARCHITECTURE.md`

⛔ **Do not keep anything "in case".** ⚠ Two architectures on disk is how the next session picks the wrong
one. ⭐ Git has it all if it is ever wanted.

---

## ⭐⭐ Decisions — ⛔ a builder does not choose these

⚠ A comprehension check found **eleven** places where two reasonable implementations were possible. ⭐ Each
is settled here. ⛔ Where it says *builder's choice*, it means the choice genuinely does not matter and the
builder should not ask.

**D1 · What counts as an end derivation.** ⭐ Export **every object S9's `MAIN` package emits.** The control
packages (`X1`–`X8`) are ablation **evidence**, ⛔ not exports. ⚠ S9 has **no** `LOCAL_` convention — that is
an S11 innovation and is not available here — so package membership is the boundary. ⭐ Mechanical, ⛔ no
per-object judgement, and it errs toward exporting more.

**D2 · The named collection's form.** ⭐ Builder's choice.

**D3 · Reconstruction format.** ⭐ Use a form SymPy reconstructs **exactly**, and ⭐ **verify the round trip
in the run** — reconstruct what was written and compare against the live object. ⭐ Carry a human-readable
rendering alongside it, since hand-checking is an intended use.

**D4 · Dimension representation.** ⭐ Write it in whatever shape the engine's own dimension solve produces.
⛔ Do not reformat, normalise, or re-type it.

**D5 · Export shape — ⛔ FLAT, ⛔ not nested by step.**

```python
# S10_exports.py — GENERATED
LEDGER = {
    "rho_br":                   {"value": ..., "dim": ..., "class": "KNOB",    "step": "S9"},
    "transverse_speed_squared": {"value": ..., "dim": ..., "class": "DERIVED", "step": "S9"},
    "polarization_count":       {"value": ..., "dim": ..., "class": "DERIVED", "step": "S10"},
}
```

⭐ Each step **imports the previous step's `LEDGER`, adds its own entries, overwrites what it derives, and
exports the merged flat dict.** ⇒ the final step's file is the whole picture, and `LEDGER["x"]` is the
current value of `x` — ⛔ no walking sections to find it.

⭐ **The `step` field records who last set the entry**, so the flat structure stays self-describing without
nesting.

⭐⭐ **Each step's export is COMMITTED.** ⇒ knob history and cross-step comparison are a **diff between two
committed files**, ⛔ not a structure to navigate. An overwrite shows up as a changed entry in the next
step's file, with `class` moving `KNOB → DERIVED` and `step` moving to the step that retired it.

⇒ ⭐ **Chain integrity is then checkable for free**: every entry the next step did not touch must be
**identical** to the previous step's. A difference means a generated file was hand-edited or a downstream
run is stale against a re-derived upstream.

**D6 · Extractor scan scope.** ⭐ An explicit list of files, passed in. ⭐ For this change that list is S9's
and S10's SymPy engines. ⛔ Do not discover engines by globbing.

**D7 · Knob inventory versus overwrite.** ⭐ A knob is undetermined **as of a step**. ⇒ the per-step
inventory is what that step consumed and did not derive. ⭐ The **frontier** is the union across steps minus
anything a later step derives — and the overwrite is what marks that. ⛔ There is no global "is this a knob"
flag.

**D8 · Annotation grammar.** ⭐ `# TAG · English description`, on the declaration line, one tag per
declaration, `TAG` from the closed set. ⛔ Nothing else on that line.

**D9 · The baseline for "same values as before".** ⭐ Run S9 **as it is today** and capture its emissions
keyed by tag name **before any edit**. That snapshot is the baseline and is committed first.
⛔ Do not reconstruct a baseline afterwards from memory or from the report.

**D10 · How "same value" is decided.** ⭐ Emit **both operands and the residual**, then guard on the
residual. ⛔ Not a boolean comparison with the answer typed beside it.

**D11 · Cross-engine object matching — ⭐ BY NAME. ⛔ Superseded the pair list.**

⭐ **The Wolfram engine names a computed object exactly as the SymPy export names it.** ⇒ post-run
comparison is a lookup by name, ⛔ not a hand-written pair table.

⚠ Today the two diverge completely — `WL_S9_RHO_DIMENSION` against `PY_S9_MAIN_DIM_PRIMARY_INERTIA` — and
that divergence is the *only* reason `checks_S9.yaml`'s pair list exists. ⇒ ⭐ once the names agree, the
pair list has no content and dies with the rest of `reduction/`.

⭐ Where only one engine derives an object (WL's solver-condition tags, py-only bookkeeping), ⭐ it simply
has no counterpart. ⛔ That is not a defect and ⛔ nothing should be invented to pair it.

**D12 · What may change in a `.wl` engine.**
⭐ **Emitted names only.** ⛔ No change to the derivation, the action, the ansatz, the premises, or any
computed value.
⭐ The rename must be shown to be a rename: **same values, new names**, against the committed output.
⚠ ⭐ A name is not a value ⇒ ⛔ this does not leak anything and ⛔ does not weaken the silo. ⭐ The Wolfram
engine still imports nothing and still re-derives the action independently — ⭐ that, and only that, is what
makes it a check.
