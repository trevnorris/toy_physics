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
5. ⭐ **S9's Wolfram engine is unchanged and still runs**, and its results still agree with py's where both
   derive the same object.

## Scope

⛔ **Do not modify any `.wl` engine.** ⚠ Its independence is the contention, and it is load-bearing.
⛔ **Do not add a YAML file, a registry, or a declared dimension.**
⛔ **Do not delete the existing `reduction/` machinery in this change** — it is superseded, ⭐ but nothing is
removed until the replacement runs on S9 → S10.
⛔ **Do not write an expected value or dimension anywhere a builder can read it.**

## Order

1. Annotate S9's declarations.
2. Restructure S9 so end derivations land in a named collection.
3. Generate `S9_exports.py`.
4. Import into S10 and make S10 run against it.
5. Write the knob extractor.
6. Run the verification list.
