# S10 — build the SymPy engine (engine 2)

**Write:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`
**And:** a sibling `…_sympy_audit.premises` file listing every supplied premise, one per line.

**Read first, in full, and treat as binding:**
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`

⛔ **Do not commit.** ⛔ **Do not modify any other file** except the two named above.

---

## ⛔⛔ WHAT THIS ENGINE MAY NOT READ

- ⛔⛔ `research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl` — **the sibling
  engine.** ⛔ Do not open it, do not grep it, do not read it via `git show`.
- ⛔ `research/pde_ledger_v3/steps/` — all of it, every file, without exception. ⚠ Some filenames in that
  directory state results; ⛔ **do not read the directory listing either.**
- ⛔ `research/pde_ledger_v3/paper/` — the TeX cards and the built PDF.
- ⛔ Any build transcript or log from the engine-1 build.

⭐⭐ **THE ENTIRE VALUE OF THIS ENGINE IS THAT IT CAN DISAGREE WITH ENGINE 1.** A transcription agrees
**vacuously** and is worth nothing. ⭐ Derive everything from `S10_SHARED_PHYSICS.md`, which is the only
artifact the two engines share.

⚠ **You are overwriting a file that already exists.** ⛔ **Do not read the existing
`S10_brane_mode_spectrum_sympy_audit.py` and do not adapt it** — it carries the defects this rebuild
exists to remove, including a duplicated emitted tag name. ⭐ Write a new file from the shared physics.

## ⭐ WHAT THIS ENGINE **MAY** READ, and must

- ✅ `research/pde_ledger_v3/reduction/` — **the registry.** Use `registry_read.load_registry`.

This is the one asymmetry between the engines, and it is deliberate.

## The mandate

Implement **Q1 through Q8** of the shared physics for **every package** in its §7 table, obeying §4 (the
structural rule), §5 (three clauses, four corollaries, no verdict), and §8 (tag grammar, `<ENGINE>` = `PY`).

### Additionally — the registry comparison, which engine 1 cannot make

Having **derived** `[ρ_br]`, `[μ_R]` and the brane dimension's role from the action (Q6):

- Load the declared dimensions of `Q.brane.rho_br` and `Q.brane.mu_R` and the declared value of
  `Q.brane.D_brane` from the registry.
- Emit, as **three separate objects** per quantity (clause 2): the **derived** dimension specialised to the
  registry's declared `D`, the **declared** dimension, and their **difference**.
- ⛔ Do **not** substitute the registry's `D` before deriving. Derive symbolically in `D` first, emit that,
  **then** specialise.
- Emit the registry's declared `D_brane` value and its dimension as their own tags.

⚠ **If the derived and declared dimensions differ, that is a FINDING, ⛔ not a build failure.** Emit the
difference and **exit 0**.

## SymPy-specific requirements

- ⭐ Use `sympy.Matrix.rank()` for Q4; for **N3** build the stacked matrix explicitly, e.g.
  `Mr.col_join(sp.Matrix([[*kvec]]))`, and take **its** rank. ⛔ Do not infer `nu_T` from `nullspace()`.
- ⚠ `Matrix.rank()` assumes generic symbols — this is exactly Q8. Compute the rank-drop locus explicitly.
- ⭐ Declare symbols with their assumptions at creation (`sp.Symbol("rho_br", positive=True)` etc.) so
  `simplify` and sign tests can use them. Emit the assumption set as a tag.
- ⭐ Use `sp.solve(..., omega_squared, dict=True)` and emit the raw solution list **before** de-duplication,
  then the de-duplicated list, then the count.
- ⛔ **No `assert` before the value it guards.** Compute → emit → *then* guard. ⚠ An `assert` that precedes
  an `emit` converts an informative value into a crash and hides it from every ablation.
- ⛔ **No `VERDICT` tag, no `PASS`/`FAIL` tag, no summary.** ⛔ `sys.exit(1)` **only** on an exception or a
  genuine operational failure.
- ⛔ **Every emitted tag name must be unique.** ⚠ The previous version of this file emitted the same tag
  name twice, which broke automated consumption of its output entirely.
- ⭐ Emit payloads via `str(expr)` (or `sp.srepr` where `str` would be ambiguous) so they re-parse.

## The `.premises` file

One line per supplied premise from shared physics §3 and §9, each stating the premise and that it is
**unfalsifiable within this build**.

## Verify before you finish

Run the script. It must complete inside **10 minutes** and exit **0**. ⭐ Paste nothing from its output into
your report.

⇒ Then answer §10 of the shared physics, **under 25 lines.**
