# S11 — build the SymPy engine (engine 2)

**Write:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
**And:** a sibling `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.premises`
**Run it into:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`

**Read first, in full, and treat as binding — it is the only artifact the two engines share:**
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`

⛔ **Do not commit.** ⛔ **Do not modify any file** other than the three named above.

---

## ⛔⛔ WHAT THIS ENGINE MAY NOT READ

- ⛔⛔ `research/pde_ledger_v3/mathematica/` — **the sibling engine**, and every other Wolfram engine in the
  ledger. ⛔ Do not open it, do not grep it, do not read it via `git show`.
- ⛔ `research/pde_ledger_v3/steps/` — **all of it, every file, without exception.** ⚠ Some filenames in that
  directory state results; ⛔ **do not read the directory listing either.**
- ⛔ `research/pde_ledger_v3/paper/` — the TeX cards and the built PDF.
- ⛔⛔ **`research/pde_ledger_v2/` — ALL OF IT, every file, without exception.** ⚠ That tree is this
  ledger's **predecessor**: it holds earlier engines *and their committed output files* for the same
  physics, under different names. ⛔ Do not read it, do not list it, do not `git show` it.
- ⛔⛔ **Every orientation, status and planning document in the repository**, at the repository root and at
  the ledger root — including `STATUS.md`, `REBUILD_HANDOFF.md`, `DEFECT_REGISTER.md`, `V3_STEP_PLAN.md`
  and any `SESSION_*.md`. ⚠ **These state results outright.** They are written for the orchestrator, and
  one of them exists specifically to tell the orchestrator things a builder must not be told.
- ⛔⛔ **`research/pde_ledger_v3/_asbuilt/` — all of it.** ⚠ It holds a **byte-identical frozen copy of the
  engine you are replacing**, kept only so the registry's line-range loci stay resolvable across the
  rebuild. ⛔ Reading it is reading the file you were told not to read, at a different path.
- ⛔ **Any other step's shared-physics file, engine, committed output, or checks config** — including the
  earlier-step engines and `out/` files sitting **in the very directory you are writing to**, and their
  `__pycache__` artifacts. ⚠ Several steps in this ledger share machinery with this one, and their
  committed outputs carry values for objects you are being asked to compute. ⭐ Compute yours; ⛔ do not
  look at theirs.
- ⛔ Any build transcript, log, or review report from any previous build.

⚠ **This list is not a spelling test.** ⛔ It is not satisfied by avoiding the exact paths named while
reaching the same content another way — `git show`, `git cat-file`, `git archive`, `git log -p`, a
`__pycache__` artifact, a symlink, or a file that merely quotes one of the above. ⭐ **The rule is the
intent: you compute every number in this build yourself.**

⭐⭐ **THE ENTIRE VALUE OF THIS ENGINE IS THAT IT CAN DISAGREE WITH ENGINE 1.** A transcription agrees
**vacuously** and is worth nothing. ⭐ Derive everything from `S11_SHARED_PHYSICS.md`.

⚠ **You are overwriting a file that already exists.** ⛔ **Do not read the existing
`S11_stray_longitudinal_sympy_audit.py`, do not `git show` it, and do not adapt it.** It carries the defects
this rebuild exists to remove — it ends in a `VERDICT` tag and its tag decomposition is not the one §8 now
requires. ⭐ Write a new file from the shared physics.

## ⭐ WHAT THIS ENGINE **MAY** READ, and must

- ✅ **The registry**, and only the registry, for §Q6r — `reduction/quantities.yaml` and its schema, reached
  through that directory's own reader, `registry_read.load_registry`. ⛔ Never parse the YAML by hand.

⚠ **`reduction/` is NOT wholesale readable.** It also holds the harness checks configs
(`reduction/checks_S*.yaml`), the comparator, and `reduction/measurements/`. ⛔ Those are **other steps'
machinery and other steps' measured output** and they fall under the ban above. ⭐ The registry is the
allowance; the directory is not.

This is the one asymmetry between the engines, and it is deliberate.

### ⛔⛔ AN ALLOWLIST, ⛔ NOT A BAN — because the registry ROUTES AROUND the ban above

⚠ **`reduction/quantities.yaml` carries `source_locus` and provenance entries that are file paths WITH LINE
NUMBERS pointing into `research/pde_ledger_v3/steps/`** — including into files that state this step's
results. ⇒ ⛔ **A builder that reads the registry and follows its provenance lands on the answer while
obeying every instruction above.** A "do not read `steps/`" rule cannot fix that, because the registry hands
you the path.

⭐ **So this is an allowlist. From every registry quantity you may use EXACTLY these, and nothing else:**

| ✅ may use | where it actually lives | how |
|---|---|---|
| `symbol_name` | `Quantity.symbol_name` | to match a registry row against a coefficient in `COEFFICIENT_ORDERING` |
| `dimension` | `Quantity.dimension` — a 3-tuple | as the **declared** vector in §Q6r |
| the declared row's locus | ⚠ **`Quantity.raw["dimension"]["provenance"]["source_locus"]`** — the **singular** field, and only it | ⭐ **as an opaque VALUE, to emit as part of the §Q6r payload — and for nothing else** |

⚠ **There is no flat `Quantity.source_locus` attribute** — it is nested under `raw`, as above. ⛔ Reaching
it through `raw` does **not** put the rest of `raw` in scope: `Quantity.value`, the aliases, the full
provenance and every relation field remain **out of scope**.

⛔⛔ **The row also carries a PLURAL `Quantity.raw["source_loci"]` list. It is OUT OF SCOPE and must not be
emitted.** ⚠ It is a **different object** — for some rows the singular field is merely one member of it.
§Q6r names the **declared dimension's singular `source_locus`**, and §8 forbids bundling two named objects
into one payload. ⇒ ⛔ Emitting the list, or both, emits the wrong object.

⛔ **Do not open, follow, `cat`, `grep`, `git show`, or in any way read the file that a locus names.** ⭐ It
is emitted **because §Q6r requires the provenance to travel with the residual**, so the step record can tell
an independent operand from a circular one — ⛔ it is not a pointer for you to follow.

### ⚠ The reader opens those files itself, and that is NOT you reading them

⚠ **Measured:** `load_registry` also loads `reduction/relations.yaml`, and its `_validate_loci` **opens and
line-counts every locus file it names** — including files under `steps/` and the engine you are replacing.

⭐ **That is the reader's own integrity check. It returns no file content to you, and it does not put any of
those files in scope.** ⇒ ⛔ Calling `load_registry` is **not** a violation of the bans above; **reading
what it touched would be.** ⛔ Do not monkeypatch, disable or work around that validation, and ⛔ do not
treat `relations.yaml` being loaded as permission to inspect it.

⛔ If you believe a computation needs a field or a file outside this allowlist, ⭐ **stop and report that in
your §10 answer instead.** That is a finding about the specification, and it is wanted.

⚠ **If `load_registry` raises `RegistryValidationError`, ⛔ do not "fix" the registry and ⛔ do not patch
around it.** ⭐ Emit what you can, and **report the exact error in your §10 answer** — that is an
orchestrator problem, not yours.

## The mandate

Implement **`Q1` through `Q11`** of the shared physics, for **every package** in its §7 table, at **every
`D`** in that package's sweep, obeying:

- §4, the structural rule;
- §5, the three clauses, the five corollaries, **no verdict**, and ⭐⭐ **the locus protocol**;
- §6's ordering requirement — ⭐ **`Q9` is computed first at each `D`**, because one package's action is
  built from its output;
- §7's `P_D` rule;
- §8's tag grammar, with `<ENGINE>` = **`PY`**.

⚠ **§Q6r is engine-local to this engine.** Emit it under the `_LOCAL_` infix of §8, and emit the tag listing
every `_LOCAL_` name this engine produced.

⚠ **If a derived and a declared dimension differ, that is a FINDING, ⛔ not a build failure.** Emit the
difference and **exit 0**. ⛔ Do not adjust the derivation to match the registry, and ⛔ do not filter,
judge or omit any row on the basis of what its `source_locus` says.

## SymPy-specific requirements

### Emission

- ⭐ Emit payloads via `str(expr)` (or `sp.srepr` where `str` would be ambiguous) so they re-parse. One line
  per tag, `TAG: <payload>`, per §8.
- ⭐⭐ **Emit each tag the moment it is computed, and flush** — `print(..., flush=True)`. ⛔ Do not
  accumulate tags and print them at the end of a package, of a dimension, or of the run.
- ⛔ **No `assert` before the value it guards.** Compute → emit → *then* guard. ⚠ An `assert` that precedes
  an `emit` converts an informative value into a crash and hides it from every ablation.
- ⛔ **No `VERDICT` tag, no `PASS`/`FAIL` tag, no summary.** `sys.exit(1)` **only** on an exception or a
  genuine operational failure — ⛔ never on a physics outcome.
- ⛔ **Every emitted tag name must be unique.** ⚠ A previous engine in this ledger emitted the same tag name
  twice, which broke automated consumption of its output entirely.

### Assumptions

- ⛔⛔ **PER-SYMBOL ASSUMPTIONS ARE NOT ENOUGH, AND THIS WILL CAUSE A FALSE CROSS-ENGINE DISAGREEMENT.**
  `sp.Symbol("k1", real=True)` does **not** give SymPy `Σ k_m² > 0`, and a sign test that Wolfram answers
  definitely will return `None` here. ⇒ ⭐ **Build §3's assumption set as an explicit JOINT predicate**
  (`sp.Q.positive(...) & …`), with the package's own domain premises joined in when that package is
  selected, and pass it to `sp.refine` / `sp.ask` for **every** sign and simplification test. ⭐ Emit the
  joint predicate itself as a tag.
- ⛔ **Add no premise relating the stiffness coefficients to one another** — not genericity, not
  distinctness, not an ordering. §Q8 asks for exactly the locus such a premise would delete.
- ⛔ Add no assumption about the sign of `ω²` anywhere in the construction of `M`.

### Solving, and the locus protocol

- ⚠⚠ **`_SOLUTION` is emitted with the domain LEFT UNRESTRICTED**, per §5's protocol table — the reality
  question is answered by `_REAL_ADMISSIBLE`, which is a separate object with its own operands. ⛔ Do not
  fold `domain=sp.S.Reals` into the solve and call the protocol satisfied.
  ⚠ §5 exists because it was **measured** on this CAS that a multivariate solve ignores real declarations,
  the single-variable real solver rejects a tuple of variables, and the ordinary solver returns the **same
  empty token** for a system that is identically satisfied and one that is inconsistent.
- ⭐ `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` are computed **from `_EQUATIONS`**, ⛔ never read off an
  empty solver token. Emit all five objects for **every** locus the spec requests, unconditionally.
### ⚠⚠ One trap that would hit BOTH engines identically

§Q3 asks for `ROOT_SOLUTION_SET` **"as returned, retaining multiplicity"** and `ROOT_DISTINCT` **"after
de-duplication"** as **two different objects**, with `ROOT_COUNT_ALL` and `ROOT_COUNT_DISTINCT` each counted
from its own list.

⚠ **Measured on this CAS: the obvious call de-duplicates silently.** For a cubic with a repeated root,
`sp.solve` returns two entries — with or without `dict=True`. ⇒ using it for `ROOT_SOLUTION_SET` would make
those two objects identical **by construction** and their two counts equal for reasons that have nothing to
do with the physics. ⛔ That is not a result; it is the solver's default.

⇒ ⭐ **Produce `ROOT_SOLUTION_SET` by whatever route actually retains multiplicity.** ⛔ This directive
names no API for it — ⭐ if you cannot obtain a multiplicity-retaining solution set, ⛔ do not fake one:
emit what you have and **report it under §10 item 3**.

### The computations

- ⭐ Use `sympy.Matrix.rank()` for `Q4`; for **N3** build the stacked matrix explicitly, e.g.
  `Mr.col_join(sp.Matrix([[*kvec]]))`, and take **its** rank. ⛔ Do not infer `nu_T` from `nullspace()` —
  §Q4 says why, and §7 says which package that warning is aimed at.
  ⭐ For **N7**, the basis count is `len(Mr.nullspace())` — a **different algorithm** from `rank()`, which is
  the whole point of that check.
- ⚠ `Matrix.rank()` assumes generic symbols — that is exactly `Q8`. Compute the rank-drop loci explicitly,
  ⭐ **including the solve-variable sets over the stiffness coefficients**, and re-run `Q3`/`Q4` on each
  stratum per `Q8b`.
- ⭐ `M.det(method="berkowitz")` is usually the workable route here; ⚠ the default path may not finish on
  the larger `D`. ⭐ Either is the same determinant, so this is a route choice, ⛔ not an object choice.
  ⚠⚠ **§Q3 requires `DET_M` FACTORED.** ⛔ Emitting an unfactored determinant instead is emitting a
  different object than the spec asks for, and the other engine will not be doing it — ⛔ do not substitute
  one. ⭐ If factoring genuinely cannot complete for a cell, that cell goes in `SKIPPED_PAIRS` and into your
  §10 report.
- ⭐ For `Q6`, walk the expression tree. `COEFFICIENT_ORDERING` is built from the action's own terms per
  §Q6, ⛔ never from a fixed list you type — §Q6d's `DIM_UNKNOWN_COUNT` is counted from it.
- ⭐ For `Q7`, compute `c_i` from the Levi-Civita definition (`sp.LeviCivita`) and compute `c·c` from
  **that result**. ⛔ Do not write out its expanded polynomial.
- ⭐ For `Q9`, emit `V1_BASIS`, `V2_BASIS` and `V6_BASIS` as the **reduced row echelon form** of the
  coefficient matrix in the pinned `MONOMIAL_ORDERING` coordinates — `Matrix.rref()` over exact rationals.
  ⚠ `Matrix.rref()` returns the pair `(rref_matrix, pivot_columns)`, ⛔ not a bare matrix — emit the matrix.
  ⛔ Do not substitute a numerical or floating-point solve for the exact one.
- ⭐ `V6_OPERATOR` is the map `p ↦ p(R₀ G R₀ᵀ)` **as a linear operator on the `V1` space** — emit the
  operator, then its `(−1)`-eigenspace. ⛔ Do not obtain `V6` by subtracting `V2` from `V1`, and ⛔ do not
  obtain `V2` as the `(+1)`-eigenspace of `V6_OPERATOR`; §Q9 says what either shortcut destroys.

## ⭐ If an object cannot be built the way the spec asks

⚠ Some objects in `Q9` and `Q11` may have only one obvious route in a given CAS, and §Q9 and §5 corollary 3
**forbid** particular routes because they make a residual identically zero. ⛔ **If the only route your CAS
offers for an object is one the spec forbids, do not take it.** ⭐ Emit what you can build honestly and
**report it under §10 item 3** — that is a finding about the specification and it is worth more than a
complete tag set obtained by taking the banned route.

⚠ The same applies to any `<QUANTITY>` name in §6/§7 you cannot emit under §8's one-tag-per-object rule.

## The `.premises` file

One line per supplied premise from shared physics §3 and §9, each stating the premise and that it is
**unfalsifiable within this build**.

## Runtime — ⛔ there is no wall-clock deadline

⚠⚠ **The rule is OBSERVABLE PROGRESS, ⛔ not total elapsed time.** A run that visibly completes
`(package, D)` cells for hours is **acceptable and expected**; §Q9's unknowns grow as `D²(D²+1)/2`.

- ⛔ **Do not narrow the declared sweep**, drop a `(package, D)` cell, or reformulate a requested object
  into a cheaper one to save time.
- ⭐ `RUN_PAIRS` and `SKIPPED_PAIRS` are **observed, not declared** — accumulate a pair only **after** that
  package has finished emitting, and emit both **after** the sweep, with `SKIPPED = declared ∖ completed`.
- ⭐ If a cell genuinely cannot complete, record it in `SKIPPED_PAIRS` **and report it in §10**. ⛔ Never
  drop one silently.

## Verify before you finish

Run the script over the **complete declared sweep**, with stdout redirected to the `out/` path above. It
must exit **0**. ⭐ **Paste nothing from its output into your report.**

## Operational

- ⛔ **Never `pkill` or `killall`** — shared machine. Kill by explicit PID only.
- ⛔ **Do not launch Mathematica or `wolframscript`.**
- ⛔ No `git add`, no `git commit`, no other git write.

⇒ Then answer **§10 of the shared physics, under 25 lines.**
