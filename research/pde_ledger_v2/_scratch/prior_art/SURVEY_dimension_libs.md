# Prior-art survey: could a Python dimensional-analysis library replace `ledger_dimensions.py`?

Read-only survey. Nothing installed, nothing in the corpus modified.
Environment: `/usr/bin/python3` = **CPython 3.10.12**, `sympy 1.14.0`.

Bottom line up front: **no**. `sympy.physics.units` is the only candidate importable here and it
can hold the algebra, but it discards the one property the ledger's cross-engine control is built
on (a *declared* axis order with zero-padded, full-length exponent vectors) and it silently accepts
three classes of input the incumbent rejects. The other four libraries are not installed and cannot
be verified by execution; on upstream source evidence, one (`astropy`) *manufactures* the exact
defect this project forbids, one (`numericalunits`) has no dimension object at all, one (`unyt`)
hardcodes its base-dimension list, and one (`pint`) is workable-with-gaps but buys nothing the
corpus uses.

---

## 0. The incumbent, measured

`research/pde_ledger_v2/scripts/ledger_dimensions.py` — 222 lines, currently pinned and green:

```
$ python3 scripts/check_ledger_dimensions_pin.py
MODULE_PIN_OK|module=scripts/ledger_dimensions.py|authority=scripts/ledger_dimensions.accepted.sha256
  |accepted_sha256=7cb8b7c93d189320f3ff857a665162e3290d1628e4511b2b387dd9a4f08e5bea|consumer_artifacts=(none)
```

API surface actually provided:

| Symbol | Behaviour that matters |
|---|---|
| `DimensionBasis(*axes, render=)` | ordered, unique, **user-declared** axis labels; `render ∈ {"symbolic","tuple"}`; raises on empty/duplicate axes |
| `basis(*exponents)` | positional construction against the declared order; `basis()` = all-zero, **full length** |
| `basis.from_mapping(dict)` | **raises** on `missing`/`extra` axes — membership validation |
| `basis.render(dim)` | `"L^2 M T^-2"` (symbolic) or `"(2,1,-2)"` (tuple) |
| `Dimension.exponents` | `MappingProxyType`, **always full length including zeros** |
| `Dimension.components()` | ordered tuple in the basis's declared order |
| `Dimension.without(axis)` | zero one axis, keep the rest (used by stage004 for the M-blind check) |
| `__mul__ / __truediv__ / __pow__` | `_combine` **raises** if axis sets differ |
| `__eq__ / __hash__` | equality on the axis→exponent map; `hash(frozenset(items))` — hashable, dict-keyable |
| `_exact(value)` | `sp.Rational(value)` unconditionally — accepts `int/float/str/sp.Rational`, **never yields a Float** |
| `dim_residual(actual, expected)` | `sp.Expr` sum of squared per-axis differences, **paired by label, never by position**; raises on axis-set mismatch |
| `emit_dimension_sidecar(stage_file, {name: Dimension})` | writes `<stage>.dimensions.txt`; enforces **one** declared axis order per sidecar; header carries `source_sha256` **and** `ledger_dimensions_sha256` |

Sidecar record format (this is the contract the whole control rests on):

```
DIMENSIONS|stage=<stem>|axes=L,T,M|source_sha256=<64hex>|ledger_dimensions_sha256=<64hex>
DIM|axes=L,T,M|name=zero|exponents={0, 0, 0}
DIM|axes=L,T,M|name=L|exponents={1, 0, 0}
```

### Consumers (the required API surface)

Converted stages and their bases — **measured from the source, not the prompt**:

| Stage | Basis | `render` | dimension call sites |
|---|---|---|---:|
| 004 `gnls_action_dimensional_foundation` | `("L","T","M")` | `symbolic` | 56 |
| 011 `frozen_reduction_certificate` | `("L","M","T")` | `tuple` | 16 |
| 012 `dtn_pole_ladder_robin` | `("L","M","T")` | `tuple` | 24 |
| 013 `breathing_harmonic_mk_projection` | `("L","M","T")` | `tuple` | 27 |
| 016 `l2_so3_covariance` | `("L","M","T")` | `symbolic` | 36 |
| 018 `dtn_hankel_fingerprint` | `("L","M","T")` | `symbolic` | 21 |
| 023 `nullspace_underdetermination` | `("L","M","T")` | `tuple` | 34 |

≈ **214 call sites** across 7 stages; 7 committed `.dimensions.txt` sidecars (129 lines total).

The two hard bases are **not yet converted**. They are exercised only by
`scripts/probe_ledger_dimensions_extremes.py`:

* `probe_stage038()` → `DimensionBasis("M", "L", "T", "E", render="tuple")` — a **4-axis** basis with
  `E` genuinely independent (`mu_r = basis(2, 1, -4, -1)`, `a_e = basis(0, 1, 0, 1)`).
* `probe_stage042()` → `DimensionBasis("stiffness", "L", "T", render="tuple")` — a **non-physical**
  base, with `charge = basis(Fraction(1,2), Fraction(3,2), Fraction(-1))` and an explicit assertion
  `assert not any(exponent.has(sp.Float) for exponent in charge.exponents.values())`.

Per-stage term-by-term walker (each converted stage carries its own `dim_of`, structurally identical):

```python
if expr.is_Symbol:
    if expr not in symbol_dims:
        raise DimError(f"missing dimension for symbol {expr}")   # <-- the load-bearing line
    return symbol_dims[expr]
if expr.is_Mul:   ...product of dim_of(arg)...
if expr.is_Pow:   ...dim_of(base) ** sp.Rational(exponent), raises on non-numeric power...
if expr.is_Add:   ...all arg dims equal else DimError...
raise DimError(f"unsupported expression in dimension checker: {expr}")
```

stage016 additionally has `assert_no_float(name, expr)` which **recurses into
`Dimension.exponents`** and raises `AuditFailure` on any `sp.Float` atom.

### Tooling that consumes the format

* `compare_dimension_artifacts.py` — parses the Python sidecar **and** the Mathematica `.out` into
  `LabelledDimension(axes: tuple[str,...], exponents: Mapping[str, Fraction])`. It **rejects**
  `len(axes) != len(exponent_values)` and any `DIM` line whose `axes` differ from the header's.
  It enforces `source_sha256` freshness, `ledger_dimensions_sha256` freshness, and the module pin,
  and fails if `compared == 0`.
* `generate_canonical_dimension_table.py` — renders every quantity in a **canonical** `(M,L,T)`
  label order regardless of the stage's declared order; `dimension_key` = `tuple(sorted(exponents.items()))`
  so comparison is position-independent. Emits `CANONICAL_DIMENSIONS.md` (207 lines, 122 rows,
  7 of 30 dimension-bearing stages represented).
* `check_ledger_dimensions_pin.py` — whole-file SHA-256 of `ledger_dimensions.py` against
  `ledger_dimensions.accepted.sha256`; `--accept` re-baselines. Wired at `scripts/run_all_audits.sh:20`,
  i.e. a stale pin fails the **entire** audit runner, not one stage.
* **The Mathematica half is hand-rolled and independent.** 7 `.wl` scripts `Print` the same records:
  `Print["DIM|axes=L,M,T|name=a|exponents=", ToString[InputForm[dimRules[a]]]]`. The library question
  is therefore a *one-engine* question; the WL side would have to be kept byte-compatible either way.

---

## 1. Availability

| Library | `import` result |
|---|---|
| `sympy` | ✅ **1.14.0** |
| `pint` | ❌ `ModuleNotFoundError: No module named 'pint'` |
| `unyt` | ❌ `ModuleNotFoundError: No module named 'unyt'` |
| `astropy` | ❌ `ModuleNotFoundError: No module named 'astropy'` |
| `numericalunits` | ❌ `ModuleNotFoundError: No module named 'numericalunits'` |

Bounded filesystem search found no copies under `/usr/lib/python3*/dist-packages`,
`/usr/local/lib/python3.10/dist-packages`, or `/home/trevnorris` (depth 6).
Installation was forbidden, so **four of five verdicts below cannot rest on execution.** They rest on
upstream source and documentation, cited line by line, and are flagged `⚠ NOT EXECUTED`.

---

## 2. Capability match

### 2a. `sympy.physics.units` — MEASURED (probes run against sympy 1.14.0)

**What works:**

| Requirement | Result |
|---|---|
| Non-physical base `(stiffness, L, T)` | ✅ `DimensionSystem((Dimension("stiffness"), L, T))` builds. `charge = S**Rational(1,2)*L**Rational(3,2)*T**-1` → `{L: 3/2, T: -1, stiffness: 1/2}`, value types `Rational` / `NegativeOne` / `Half` |
| 4-axis `(M, L, T, E)` with `E` independent | ✅ `mu_r = M**2*L*T**-4*E**-1` → `{L: 1, E: -1, M: 2, T: -4}`; `is_dimensionless(E/E)` → `True` |
| Exact rationals when *given* rationals | ✅ zero `Float` atoms; `L**sp.Rational('1/2')` → `{L: 1/2}` |
| Hashability / dict key | ✅ `L*T**-1 == L/T` → `True`, hashes equal, `{L*T**-1: "speed"}[L/T]` → `"speed"`, `len({L*T**-1, L/T}) == 1` |
| Add-homogeneity checking | ✅ (private API) `SI._collect_factor_and_dimension(meter + second)` → `ValueError: Dimension of "second" is Dimension(time, T), but it should be Dimension(length, L)`; also works over a *custom* non-physical system (`ValueError: Dimension of "tau" is Dimension(tau), but it should be Dimension(ell*k_stiff)`) |
| Extras the incumbent lacks | `equivalent_dims`, `dim_vector` (a `Matrix`), `can_transf_matrix`, `inv_can_transf_matrix`, `is_consistent`, `print_dim_base`, `dimsys_SI`'s derived-dimension table |

**Blocking gaps — every one measured, none inferred:**

1. **Declared axis order is discarded.** `DimensionSystem` sorts its bases.
   ```
   declared ('M','L','T') -> base_dims (L, M, T)
   declared ('L','T','M') -> base_dims (L, M, T)
   declared ('L','M','T') -> base_dims (L, M, T)
   ```
   All three collapse to sorted `(L, M, T)`; `dim_vector` is likewise in sorted order.
   stage004's committed sidecar declares `axes=L,T,M`. The library **cannot** produce that label
   or that positional vector; the order must be re-imposed by a wrapper the library does not supply.

2. **Zero-exponent axes are dropped.** `DimensionSystem((M,L,T)).get_dimensional_dependencies(L)`
   → `{Dimension(L): 1}`, **not** `{L:1, M:0, T:0}`. Every `DIM|` record needs a full-length vector,
   and `compare_dimension_artifacts.build_dimension` raises when
   `len(axes) != len(exponent_values)`. Reconstruction via `.get(axis, 0)` is mandatory.

3. **No membership validation — typos become new base dimensions.**
   ```
   dsA = DimensionSystem((L,M,T))
   dsA.get_dimensional_dependencies(L * Dimension("stiffness")) -> {L: 1, stiffness: 1}   # no error
   dsA.get_dimensional_dependencies(Dimension("typo_axis"))     -> {typo_axis: 1}         # no error
   dimsys_SI.get_dimensional_dependencies(Dimension("stiffness")) -> {stiffness: 1}       # no error
   DimensionSystem((stiffness,L,T)).get_dimensional_dependencies(M) -> {M: 1}             # no error
   ```
   The incumbent's `from_mapping` raises `missing=[...], extra=[...]`. This is a **control regression**,
   not a convenience gap: a misspelled axis would silently pass and silently change the answer.

4. **Float exponents are silently accepted.** `get_dimensional_dependencies(L**0.5)`
   → `{Dimension(L): 0.500000000000000}`, `sympify(v).atoms(sp.Float)` = `{0.5}`.
   The incumbent's `_exact` coerces to `sp.Rational` at construction, so a Float can never enter.
   Under the library, the earliest possible catch is stage016's `assert_no_float` at *use* time —
   or nowhere at all in the six stages that don't have that helper.

5. **Cross-basis contamination.** `Dimension("L")` is a global singleton keyed by name, so two stages
   with different bases share it and can be combined without complaint. Worse,
   `Dimension("mass") == sympy's SI mass dimension` → `True`, so a stage that names an axis `mass`
   silently aliases the SI dimension (measured; note `Dimension("M") == mass` is `False`, so the
   aliasing is name-dependent and therefore easy to trip over by accident).
   The incumbent's `Dimension.__eq__`/`_combine` compare axis-keyed maps and raise on differing sets.

6. **The term-by-term walker keys off `Quantity`, not `Symbol`.**
   ```
   SI._collect_factor_and_dimension(sp.Symbol('x')) -> (x, Dimension(1))
   ```
   **An unmapped symbol is silently DIMENSIONLESS.** The corpus's `dim_of` raises
   `DimError("missing dimension for symbol …")`. To use sympy's walker every physics symbol would
   have to be wrapped in a `Quantity` with a `set_quantity_scale_factor` (a scale factor no stage has
   any use for), and the walker itself is private: `dir(SI)` exposes only `get_dimensional_expr`,
   `get_dimension_system`, `get_quantity_dimension`, `set_quantity_dimension` publicly —
   `_collect_factor_and_dimension` carries the underscore.

7. **No serialisation, no ordered rendering.** `str(L**2*M/T**2)` = `"Dimension(L**2*M/T**2)"`;
   `hasattr(d, "components")` and `hasattr(d, "exponents")` are both `False`; there is no tuple
   render, no `DIM|` emitter, no digest, no per-basis label set.

### 2b. `pint` — ⚠ NOT EXECUTED (upstream source evidence)

* **Non-physical base: YES.** `docs/advanced/currencies.rst` defines `USD = [currency]` and
  `EUR = [currency_EUR]`; `docs/advanced/defining.html` states *"primary dimensions don't need to be
  declared; they can be defined for the first time as part of a unit definition."*
  So `[stiffness]` and a fourth `[E]` axis are expressible — the **hard cases pass on paper**.
* **Exact rationals: NOT BY DEFAULT.** `pint/registry.py` — `UnitRegistry.__init__(..., non_int_type=float, ...)`.
  `pint/util.py UnitsContainer`:
  ```python
  if not isinstance(value, Number):
      raise TypeError(f"value must be a number, not {type(value)}")
  if not isinstance(value, int) and not isinstance(value, self._non_int_type):
      d[key] = self._non_int_type(value)
  ```
  With the default registry, an exponent of `1/2` becomes **float `0.5`**.
  Fixable via `UnitRegistry(non_int_type=Fraction)` — but `non_int_type` is **registry-wide** and
  applies to quantity magnitudes too, so every numeric that touches the registry must become a
  `Fraction`. That is a real coupling, not a flag flip.
* **Declared axis order: NO.** `.dimensionality` is a `UnitsContainer` (an immutable dict). Ordering
  and zero-padding would still be the project's job.
* **Hashability:** `UnitsContainer` is an immutable mapping and defines `__hash__` — expected to work,
  **not verified**.
* **SymPy attachment: NONE.** No expression walker. The stages' `dim_of` would be kept verbatim with a
  pint dimensionality object as the dict value. Exponents would be `Fraction`, not `sp.Rational`;
  `dim_residual` returns `sp.Expr` and `assert_no_float` sympifies, so a `Fraction` converts cleanly —
  but it is a new conversion boundary in the one place the project cares most about exactness.
* **Cost of entry:** each stage's basis becomes a `.txt` definitions file, and the library is a new
  **unpinned** third-party dependency.

### 2c. `unyt` — ⚠ NOT EXECUTED (upstream source evidence)

`unyt/dimensions.py` defines base dimensions as module-level SymPy Symbols:
```python
mass = Symbol("(mass)", positive=True)
...
base_dimensions = [mass, length, time, temperature, angle,
                   current_mks, dimensionless, luminous_intensity, logarithmic]
```
`base_dimensions` is a **hardcoded module-level list**. `define_unit`/`UnitRegistry.add` register
*units* whose dimensions must be expressed in the existing base dimensions; no public API to register
a **new base dimension** was found. Minting `Symbol("(stiffness)", positive=True)` yourself produces a
symbol that unyt's registry machinery does not know about — at which point you are using raw SymPy,
not unyt, i.e. option 2a. Monkey-patching a third-party global list is strictly worse than 222 pinned
lines for a control that is supposed to be byte-reproducible.

unyt's actual value is fast unit-aware arithmetic on NumPy arrays. This corpus has no arrays and
performs no unit conversion.

### 2d. `astropy.units` — ⚠ NOT EXECUTED (upstream source evidence) — ⭐ decisive

`astropy/units/utils.py`, `sanitize_power` — *"Convert the power to a float, an integer, or a Fraction."*
```python
if denom == 1:
    return int(p.numerator)
elif (denom & (denom - 1)) == 0:     # denominator is a power of two
    p = float(p)
```
**Denominator 2 is a power of two, so `1/2` and `3/2` become float `0.5` and `1.5` by design.**

That is precisely the corpus's hardest case: `probe_stage042` builds
`charge = basis(Fraction(1,2), Fraction(3,2), Fraction(-1))` and then asserts *no* `sp.Float` is
present; stage004 raises dimensions to `sp.Rational(1, 2)` for the sound-speed and healing-length
checks. astropy's exponent normalisation would **manufacture** the exact defect this project treats
as fatal, and would do so silently, in the library, below any wrapper.

Custom irreducible bases via `def_unit` (no `represents=`) appear possible, so the
`(stiffness, L, T)` and `(M, L, T, E)` *label* structure is probably expressible — but the exponents
inside them would not survive.

### 2e. `numericalunits` — ⚠ NOT EXECUTED (PyPI documentation)

> *"A complete set of independent base units (meters, kilograms, seconds, coulombs, kelvins) are
> defined as randomly-chosen positive floating-point numbers."*
> *"In a dimensionally-correct calculation, the units all cancel out, so the final answer is
> deterministic, not random. In a dimensionally-incorrect calculation, there will be random factors
> causing a randomly-varying final answer."*
> *"If you have a quantity with units, you cannot directly see what the units are. You are supposed
> to already know what the units are."*

There is no dimension object, therefore no exponents (exact or otherwise), no axis labels, no
hashable dimension key, no serialisable record, and nothing to attach to a SymPy expression.
Detection requires running the calculation **twice in separate Python sessions** and comparing —
i.e. it is stochastic and session-dependent, which is flatly incompatible with a SHA-pinned,
byte-reproducible sidecar and a deterministic cross-engine comparison.

---

## 3. What each library gives that the incumbent does NOT

| Library | Genuine additions |
|---|---|
| `sympy.physics.units` | `equivalent_dims`, `dim_vector`, `can_transf_matrix`/`inv_can_transf_matrix` (change-of-basis), `is_consistent`, an Add-homogeneity-checking walker (private), `Quantity` scale factors + real unit conversion, `dimsys_SI`'s derived-dimension table. **Already a dependency → zero new pin surface.** |
| `pint` | A definitions-file DSL (a stage's basis becomes *data*, not code), unit parsing/formatting, contexts/equivalency machinery, a large physical-unit database, `@ureg.check(...)` decorators, NumPy support |
| `unyt` | `unyt_array` — unit-aware NumPy arrays with negligible overhead; astro/physics unit sets |
| `astropy.units` | `physical_type` naming, equivalencies (`u.spectral()` etc.), `decompose()`, format round-tripping (FITS/CDS/OGIP/LaTeX), `Quantity` on ndarray |
| `numericalunits` | Zero runtime overhead; works inside any numeric library without wrapping |

**None of the five provides:** a per-stage *declared* axis ORDER, zero-padded full-length exponent
vectors, a `DIM|axes=…|name=…|exponents={…}` record, or a source + module SHA-256 self-attestation.
Those four are exactly what `compare_dimension_artifacts.py`,
`generate_canonical_dimension_table.py`, and the seven `.wl` scripts consume. Of everything in the
table above, the only item the corpus has no substitute for is
`dim_vector`/`can_transf_matrix` — and no consumer currently needs a change of basis.

---

## 4. Verdicts

| Library | Verdict |
|---|---|
| `sympy.physics.units` | **COULD REPLACE WITH GAPS** — as a *storage/algebra core only*, under a wrapper that re-adds: (1) declared axis order, (2) zero-padded full-length vectors, (3) basis-membership validation, (4) exactness coercion, (5) per-basis isolation (global `Dimension("L")` singleton; `Dimension("mass")` aliases SI), (6) ordered/tuple rendering + the sidecar emitter, (7) `Symbol`-keyed walking (its walker treats an unmapped `Symbol` as dimensionless). Hard cases both **PASS**. |
| `pint` | **COULD REPLACE WITH GAPS** ⚠ not executed — gaps: default `non_int_type=float` silently floats `1/2`; the fix (`non_int_type=Fraction`) is registry-wide and drags magnitudes with it; no declared axis order; no zero-padded vector; no SymPy interop; a `.txt` definitions file per basis; a new unpinned dependency. Non-physical `[stiffness]` and a 4th `[E]` axis are documented as supported. |
| `unyt` | **CANNOT** (without monkey-patching) ⚠ not executed — `unyt/dimensions.py` `base_dimensions` is a hardcoded module-level list of 9 SymPy Symbols; no public API registers a *new base dimension* (`define_unit` registers units over existing bases). Adding `stiffness` means patching a third-party global, which defeats the point of a byte-pinned control. |
| `astropy.units` | **CANNOT** ⚠ not executed — `astropy/units/utils.py sanitize_power` deliberately converts any power whose denominator is a power of two to a **float** (`elif (denom & (denom-1)) == 0: p = float(p)`). Both the `(stiffness,L,T)` charge `(1/2, 3/2, -1)` and stage004's `sp.Rational(1,2)` powers hit this branch. In a project where a float exponent is a defect, the library's normalisation *is* the defect. |
| `numericalunits` | **CANNOT** — no dimension object exists ("*you cannot directly see what the units are*"); units are random floats; error detection is stochastic across two separate Python sessions. Fails every one of the five capability requirements, and its non-determinism is incompatible with SHA-pinned artifacts. |

---

## 5. Migration cost

### Mechanical (a Codex build session, largely shim-able)

* **7 basis declarations** — `DIMENSION_BASIS = DimensionBasis(...)` + `Dim = DIMENSION_BASIS` in
  stages 004, 011, 012, 013, 016, 018, 023.
* **≈214 dimension call sites** across those 7 files (56/16/24/27/36/21/34). Almost all are
  `Dim(a,b,c)` constructions and `dim_residual(...)` calls; a shim preserving the names
  `DimensionBasis` / `Dimension` / `dim_residual` / `emit_dimension_sidecar` makes most of them no-ops.
* **7 sidecar regenerations.** Unavoidable: each header embeds `ledger_dimensions_sha256`, so
  changing the module invalidates all 7 at once; `compare_dimension_artifacts.require_fresh_python_sidecar`
  fails loudly until they are regenerated.
* **Pin re-baseline** — one command: `python3 scripts/check_ledger_dimensions_pin.py --accept`.
* **`CANONICAL_DIMENSIONS.md` regeneration** — 207 lines / 122 quantity rows / 6 candidate groups /
  3 `NEEDS_ADJUDICATION`, via `generate_canonical_dimension_table.py`.
* **Acceptance re-run** — `python3 scripts/compare_dimension_artifacts.py NNN` must exit 0 for each of
  004/011/012/013/016/018/023, plus `bash scripts/run_all_audits.sh` (the pin gate is at line 20).

### Needs judgement (the expensive part)

1. **The pin's meaning changes — this is a trust-root decision, not a code edit.**
   Today `ledger_dimensions.accepted.sha256` pins *all* dimension semantics into 222 hashed lines, and
   the sidecar self-attests that digest. Under a library, the SHA covers only the wrapper; the
   library's version becomes an **unpinned input the sidecar does not attest**. Either add a second
   pin (library version + its own hash, and a way to detect a silently upgraded site-package), or
   accept that the control no longer covers the code computing the exponents. `sympy` is the least-bad
   case here because it is *already* an unpinned input to every stage.
2. **Re-probe the hard bases BEFORE converting anything.** `(stiffness,L,T)` and `(M,L,T,E)` exist
   only in `probe_ledger_dimensions_extremes.py`; 23 of 30 dimension-bearing stages are unconverted.
   Choosing a library without re-running that probe risks blocking stages 038 and 042 mid-rewrite.
3. **The Mathematica half doubles the per-stage cost.** 7 `.wl` scripts hand-`Print` the same
   `DIM|axes=…|name=…|exponents={…}` records, and the comparator requires the two engines' axis lists
   to match exactly. Any format change (reordered axes, dropped zero exponents, a different exponent
   spelling) requires matching WL edits. The failure mode is at least **loud** —
   `load_dimensions` raises on both `axes != declared_axes` and a length mismatch — but it is 7 more
   files, in a second language, subject to the dual-engine rule.
4. **Float policy has to be relocated.** The incumbent forbids floats at *construction* (`_exact`).
   With pint (default) or astropy, floats arrive from the library and the earliest catch is
   stage016's `assert_no_float` at use time — six stages have no such helper. Coercing in the wrapper
   restores the guarantee, but then the coercion is the thing to trust — which is `_exact` again.
5. **`dim_residual`'s label pairing must be rebuilt.** It pairs by axis label and raises on axis-set
   mismatch. Every candidate returns a *sparse* mapping, so the pairing must be re-derived from the
   declared axis list before subtraction — reintroducing exactly the code being replaced.

### Sizing

The incumbent is **222 lines, pinned, green, and dual-engine-matched across 7 stages**. A
`sympy.physics.units`-backed replacement meeting the same contract removes roughly the 35 lines of
`_combine`/`__mul__`/`__truediv__`/`__pow__` exponent arithmetic and must add **seven** guard layers
(axis order, zero-padding, membership validation, exactness, per-basis isolation, rendering, symbol
walking) plus the unchanged sidecar emitter. It is not smaller, it is not simpler, and it imports a
namespace where `Dimension("mass")` silently aliases SI's mass.

**Recommendation: keep the incumbent.** The only borrowing worth revisiting is
`DimensionSystem.dim_vector` / `can_transf_matrix` if a future stage ever needs a change of basis;
that can be added on top later without touching the storage type or the pin.

---

## 6. What could not be verified

* **`pint`, `unyt`, `astropy`, `numericalunits` were never executed.** They are not installed and
  installation was forbidden. Their verdicts rest on upstream source (`pint/util.py`,
  `pint/registry.py`, `unyt/dimensions.py`, `astropy/units/utils.py`) and documentation (PyPI,
  readthedocs), cited inline above. Every such claim is marked `⚠ NOT EXECUTED`.
* Specifically **not** verified by running code:
  `UnitsContainer.__hash__`; the exact shape of `Quantity.dimensionality`; whether
  `UnitRegistry(non_int_type=Fraction)` actually round-trips `[stiffness]**Fraction(1,2)` without
  a float anywhere; whether astropy `def_unit` can build a genuinely irreducible non-physical base;
  whether unyt tolerates a foreign `Symbol` in a dimension expression.
* `pint`'s user-facing documentation page for `non_int_type` was not located (two candidate URLs
  404'd / did not mention it). The evidence is the `UnitRegistry.__init__` signature
  (`non_int_type=float`) and the `UnitsContainer` coercion code.
* The `numericalunits` GitHub README 404'd; the quotations are from the PyPI project page.
* This was a read-only survey: the 7 stage scripts, `compare_dimension_artifacts.py`, and
  `generate_canonical_dimension_table.py` were **read**, not run.
  `check_ledger_dimensions_pin.py` was run **without** `--accept` (read-only) and reported
  `MODULE_PIN_OK`.
* The 23 unconverted dimension-bearing stages were not audited for further exotic bases beyond the
  two recorded in `probe_ledger_dimensions_extremes.py`; a third hard case could still be lurking.
