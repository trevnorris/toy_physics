# BUILD — the engine-output check harness

You are building a Python module that **consumes an engine script's tagged stdout** and runs automated
checks on it, so that verifying a 300-tag derivation stops depending on a human reading every tag.

**Deliverable (absolute paths):**
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py` — the module + CLI
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/test_engine_output_checks.py` — pytest tests
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S9.yaml` — the S9 configuration

⛔ Do not modify any other file. ⛔ Do not commit anything to git.

**Real input to develop against** (⛔ do NOT invent fixtures for the main path — this is a parser, and a
parser tested only on synthetic input is untested):

```
cd /var/projects/toy_physics
timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl > /tmp/s9.txt
```

⚠ Mathematica has a **two-seat licence** on this machine and seats may be busy. Run **one** kernel at a
time; on a licence error, wait and retry. You only need to produce `/tmp/s9.txt` once.

---

## ⭐⭐⭐ THE CONSTRAINT THAT GOVERNS THE WHOLE DESIGN

> ⛔⛔ **NO EXPECTED PHYSICS VALUE MAY EVER BE STORED IN A CONFIG FILE, A TEST, OR THE MODULE.**

**Why:** a stored expected value is an *acceptance criterion referencing an expected value*. Engine
scripts here are built by an agent that iterates until its build exits 0. If a file in the tree states
what a root or a speed should be, that agent can grep it and fix the computation until it matches — and
a genuine disagreement, which is the most valuable output this project can produce, becomes silent
confirmation.

⇒ ⭐ **Every comparison target must be GENERATED AT COMPARE TIME** — from the other engine's run, from
the same engine's control runs, or from the registry. **Config files map TAG NAMES, never values.**

⚠ **One narrow exception, and it must be justified in a comment where it appears:** *definitional*
dimensions of primitives — a wavenumber is an inverse length, `ω` an inverse time, a coordinate a
length. These are conventions, ⛔ not results of any step. **Anything a step derives — the dimensions of
`rho_br` and `mu_R`, any root, any speed — must come from the engine's own output, never from config.**

---

## 1. Parsing

### 1.1 Tagged output
Engine stdout is a sequence of `TAG: value` lines. Continuation lines belong to the preceding tag.
Match the grammar already used by `reduction/derived_or_declared.py` so the two agree:
tag matches `^([A-Za-z][A-Za-z0-9_]*)[ \t]*(?::|->)[ \t]*(.*)$`.

⛔ **Raise on a duplicate tag. Raise on untagged output before the first tag. Raise on empty output.**
⚠ A duplicate tag is a real defect that has already blocked this pipeline once.

### 1.2 CAS normalisation — ⛔ this is where a silent bug would be worst
Convert a Mathematica `InputForm` string into a SymPy object. The grammar you need is small:

| Mathematica | SymPy |
|---|---|
| `a^b` | `a**b` |
| `{a, b, c}` | a list |
| `{{a,b},{c,d}}` | nested list (a matrix) |
| `x -> y` | a `(lhs, rhs)` pair; `{{s -> v}}` is a solution set |
| `E^x`, `I`, `Sqrt[x]`, `Rational` | `exp(x)`, `I`, `sqrt(x)` |
| `True` / `False` | Python `bool` |
| a bare integer | Python `int` |
| `Derivative[i,j,k,l][f][t,x,y,z]` | keep as an opaque symbolic atom; ⛔ do not try to model it |

⛔⛔ **FAIL LOUDLY. NEVER SILENTLY.** Anything you cannot parse is recorded as `UNPARSED` with its raw
text and **counted in the report**. ⛔ Never fall back to string comparison, never coerce, never skip.
⚠ **A normaliser that silently mis-parses two different expressions into agreement is worse than no
harness at all**, because it manufactures confidence.

⭐ **Comparison is SYMBOLIC, never textual:** two expressions agree iff
`sympy.simplify(a - b) == 0`. Two lists agree iff same length and elementwise agreement. Order-insensitive
comparison (for root multisets) must be an explicit, separate function.

---

## 2. The four checks

Each returns structured data. ⛔ None of them returns a verdict on the physics.

### Layer 3 — CONTROL RESPONSE (build this first; it needs no configuration at all)
The engines emit control runs alongside the main run, under tag prefixes (in S9: `WL_S9_X1_` … `WL_S9_X6_`).

For each main tag, find its counterparts across control prefixes and report whether the value **differs
under any control**. Partition the tags into:

- **RESPONSIVE** — differs under at least one control ⇒ it depends on the action;
- **INVARIANT** — identical under every control ⇒ ⭐ **for adjudication**: either a genuine invariant, or
  a value **no computation produced**.

⛔⛔ **INVARIANT IS NOT A FAILURE AND MUST NEVER BE PRINTED AS ONE.** Some quantities are *supposed* to be
invariant — anything built only from definitional primitives, for instance. ⭐ This is a **triage** list
that shortens human reading from hundreds of tags to a handful; ⛔ it is not a verdict.

### Layer 4 — DIMENSIONAL HOMOGENEITY of every emitted expression
⭐⭐ **Take the dimensions of derived quantities FROM THE ENGINE'S OWN OUTPUT.** Config names which tag
carries which symbol's dimension vector — **names only**:

```yaml
derived_dimensions:
  rhoBr: WL_S9_RHO_DIMENSION
  muR:   WL_S9_MU_DIMENSION
```

Config also carries the definitional primitives (the narrow exception above), in the `[L, T, M]`
exponent convention the registry uses:

```yaml
primitive_dimensions:            # definitional conventions, NOT results of any step
  kx: [-1, 0, 0]                 # a wavenumber is an inverse length
  omega: [0, -1, 0]
  # ... and so on
dimensionless: [lambdaRho, lambdaMu, D]
```

Then, for each emitted expression, walk the SymPy tree and compute its dimension vector:
- `Add`: **every summand must have the same dimension**, else `NON_HOMOGENEOUS` — ⭐ report the summands
  and their differing dimensions, ⛔ not just a boolean;
- `Mul`: dimensions add; `Pow` with an integer/rational exponent: dimension scales;
- an unknown symbol: report `UNKNOWN_SYMBOL` and name it. ⛔ Never assume dimensionless.

⚠ Exponents may be **symbolic** (`D` appears in the brane dimensions), so dimension components are SymPy
expressions and equality is `simplify(a - b) == 0`, ⛔ not `==`.

⭐ The check asks: *are the expressions this engine computed homogeneous under the dimensions this same
engine derived?* That is genuine internal consistency and needs no external answer key.

### Layer 1 — CROSS-ENGINE AGREEMENT (build the mechanism; S9 has only one engine, so it will be unused)
Config maps a physical quantity to **one tag name per engine** — ⛔ never to a value:

```yaml
cross_engine:
  - quantity: transverse_dispersion
    wl: WL_S10_OMEGA_SQUARED
    py: S10_OMEGA_SQUARED
```

Load both engines' outputs, normalise, compare symbolically, report `AGREE` / `DISAGREE` / `UNPARSED` /
`MISSING`. ⭐ A `DISAGREE` is a **finding to be reported**, ⛔ never an error to be fixed by the harness.

### Layer 2 — REGISTRY RESIDUAL
Config names a relation id in `reduction/relations.yaml` and maps the relation's input qids to the tags
carrying the engine's computed values. Substitute and check the residual simplifies to zero. Reuse
`reduction/registry_read.py` for loading; ⛔ do not reimplement registry parsing.

⚠ `relations.yaml` uses a prefix expression language (`[Sub, [Q, ...], [Sqrt, [Div, ...]]]`). Convert it
to SymPy. ⛔ If a relation uses a head you do not support, raise — ⛔ do not skip it.

---

## 3. CLI

```
python engine_output_checks.py --config checks_S9.yaml --output /tmp/s9.txt
```

Print a **short** report a human can act on: counts, then the INVARIANT list, the NON_HOMOGENEOUS list,
the UNPARSED list, and any DISAGREE rows. ⛔ Do not print hundreds of RESPONSIVE tags.

End with one line naming what the run does **not** establish. ⭐ The harness is a **triage tool, not a
verdict**, and its report must say so in its own output — this project has previously had a tool's
output quoted as coverage it never had.

Exit non-zero only on **operational** failure (unparseable config, missing file, a check that could not
run). ⛔ **Do NOT exit non-zero on a physics disagreement or on a non-empty INVARIANT list** — those are
findings for a human, and a build agent must never be able to make them go away.

---

## 4. ⭐⭐ ABLE-TO-FAIL — required, and it is the point of the tests

⚠ A check that cannot fail is worse than no check. `test_engine_output_checks.py` must **prove each
check fires**, by corrupting real input and asserting the check catches it:

1. take the **real** S9 output, change one control's value to equal the main value ⇒ that tag becomes
   INVARIANT;
2. corrupt an expression so a sum mixes dimensions ⇒ `NON_HOMOGENEOUS`, naming the summands;
3. feed a value that is not valid CAS text ⇒ `UNPARSED`, ⛔ never a silent skip;
4. duplicate a tag ⇒ raises;
5. two engine outputs that disagree symbolically but not textually (e.g. `a/b` vs `a*b^-1`) ⇒ `AGREE`;
   and two that agree textually in part but differ symbolically ⇒ `DISAGREE`.

⭐ **Test 5 is the one that matters most** — it proves the comparison is symbolic rather than textual, in
both directions.

⚠ You may use small synthetic strings for the *unit* tests, but the **end-to-end test must run against
the real S9 output file**, not a fixture. A parser exercised only on invented input is untested.

Run `python -m pytest reduction/test_engine_output_checks.py -q` from
`/var/projects/toy_physics/research/pde_ledger_v3/` and iterate until it passes.

---

## 5. Report back — under 40 lines

1. Paths of the three deliverables and the pytest result.
2. The literal CLI report from running the harness on the **real** S9 output.
3. For each of the five able-to-fail tests, one line: what you corrupted, what fired.
4. Anything in §2 you could not implement, and what blocked it.
5. ⭐ Anything the harness reported that surprised you. This is wanted.
