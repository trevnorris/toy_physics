# S11 SymPy engine — fix round 2. ⭐ Eight items. ⛔ THREE BLOCK; ⭐ one is spec compliance.

## Authority and boundary

Edit `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` **only**. ⛔ Change no other
file. `CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics
authority. `research/pde_ledger_v3/directives/S11_sympy_build_directive.md` is **closed**.

⭐ **Baseline is `19d4a44a`.** ⭐ Round 1 succeeded at what it targeted: `run_cell('MAIN',5)` completes, the
rank construction was verified by two independent legs to be the exact algebraic rank of the live `M_r`
(34/34 agreements), and a **FORM ablation** confirmed nothing is typed. ⛔ Do not revisit any of that.

Every claim below was measured by a review leg; probes at `/home/trevnorris/.s11_build/leg_probes/`.
Evidence: `directives/_measurements/S11_engine_fix_round2_brief.md`.

## ⛔⛔ 1 · BLOCKING — `write_exports` CANNOT RUN

`:1693` builds `REPO_ROOT / 'pde_ledger_v3' / 'directives' / 'S11_SHARED_PHYSICS.md'`. `REPO_ROOT` is
`SCRIPT_DIR.parents[2]`, so that path is missing a component and the file does not exist. ⚠ Pre-existing
and **masked for the project's whole history** — `MAIN` D5 never completed, so the guard was never
satisfied and this line never ran. ⇒ measured end-to-end: `EXIT=1`, `FileNotFoundError`, no export.
⭐ **Make the digest read the file that exists.**

## ⛔⛔ 2 · BLOCKING — AN EXPORT FAILURE IS RECORDED AS A PHYSICS-CELL FAILURE

The mid-loop publish sits **inside the cell's `try`** (`:1742`/`:1756`/`:1759`), so a failure in
`write_exports` emits `CELL_EXCEPTION` for a physics cell while that cell **stays in `completed_pairs`** and
`PY_S11_RUN_PAIRS` lists it as run. ⇒ ⛔ **a run record that contradicts itself**, which is what round 1's
brief forbade. ⚠ Independent of item 1: *any* future export failure is mislabelled this way.
⭐ **A publish failure must be attributable to the publish**, ⛔ never to a cell.
⚠ The end-of-run `write_exports` (`:1772`) raises **uncaught** and escapes `run()`, so §10 is never emitted
and the stale-export deletion never runs. ⇒ ⭐ **a run must still explain its own failure.**

## ⛔⛔ 3 · BLOCKING REGRESSION — A COMPUTED NULL VECTOR BECAME EMPTY

`:998` now requires *proved* non-zero (`algebraic_zero_test(minor) is False`). `algebraic_zero_test`
(`:878-903`) returns `None` whenever `sp.Poly` raises at `:897` into the bare `except` at `:899-900` — which
happens for a radical of a compound expression. ⇒ at `XKIN_ANISO` D3 roots 2 and 3 **every** candidate
minor is `None`, no source minor is selected, and the payload is now `Tuple()` where the pre-fix engine
emitted a real null vector.

⛔⛔ **This is the anisotropic-inertia control — the package that tests the isotropic-inertia premise — and
two of its three modes now have no polarisation vector.** ⭐ Restore the computed object. ⛔ Do not restore
it by reverting to a structural test; ⭐ state the property: **an undecided zero test must not silently
discard a candidate.** ⛔ Do not specify the routine (rule 3).

## ⛔ 4 · THE MID-LOOP PUBLISH WRITES 18 FALSE ROWS — ⭐ carried over UNFIXED from round 1

`:1713` computes `skipped = declared_pairs − completed_pairs`, so at MAIN-completion the export asserts
**18 control cells were skipped when not one had been attempted**, with `value_kind: 'COMPUTED_OBJECT'`.
⚠ Round 1's brief required *no row claims a run record it does not have*, and said to **report the conflict
under §10 and leave the placement alone**. ⛔ It published anyway.
⇒ ⭐ **Both properties, or neither.** ⛔ If they cannot both hold, ⭐ report under §10 and revert to
end-of-run publication only.

## ⛔⛔ 5 · SPEC COMPLIANCE — A COMPUTATION IS BEING REFUSED BY A HARD-CODED COUNT

`:770` (`if len(equations) > 4`) and `:748` (`if len(residuals) > 4`) **short-circuit before any CAS call**
and emit `SOLVE_UNAVAILABLE_EXACT_SYSTEM_SIZE` / `UNAVAILABLE_EXACT_SYSTEM_SIZE`.

⚠ **A leg ran one of the refused systems and it returned in well under a second.** ⇒ the refusal is not a
capability limit.
⛔ Against the physics authority: `S11_SHARED_PHYSICS.md:245` requires `_SOLUTION` be *the solution set
exactly as your CAS returns it*; `:285` permits `NOT_APPLICABLE` **only** where a residual is
non-polynomial — ⚠ these are polynomial; and neither token is defined by the spec.
⚠ Scale: **9 / 27 / 36 / 36** such tokens at `MAIN` D2 / D3 / D4 / D5 — ⛔ this is the majority of the Q8
locus machinery at **every** dimension, not a D5 effect. ⚠ `_REAL_STATUS` reads `PROVED_NONEMPTY` beside
them.

⇒ ⭐ **Emit what the CAS returns.** ⛔ Do not substitute an undefined token for a computation that
completes. ⭐ If some system genuinely does not return, ⭐ that case reports under §10 — ⛔ a size count is
not evidence that it will not.

## ⛔ 6 · THE §10 CHANNEL IS OVER QUOTA BEFORE THE SECOND PACKAGE RUNS

`:1786` renders `ISSUES[:20]`. Measured: `MAIN` alone produces **72** entries (D2=6, D3=18, D4=24, D5=24),
and `MAIN` runs first. ⇒ ⛔ **every issue round 1 added (`:916 :932 :1016 :1512`, and `CELL_EXCEPTION` at
`:1764`) can never reach the report.** ⚠ `:749` writes its entry with no package/`D` prefix, so those are
unattributable, and the mid-loop `merged_export` duplicates every export-path entry.
⇒ ⭐ §10 must be able to carry what the run actually recorded, ⛔ and each entry must say which cell it came
from.

## ⭐ 7 · A THIRD WALL — ⛔ MEASURE IT, ⛔ do not assume it

Both the pre-fix and post-fix engines emit `..._KW_ZERO_LOCUS_EQUATIONS` and then stall in the next
statement, `sp.solve` at `:776`, inside `XKIN_ANISO` D3/D4. ⚠ Never reached by any earlier run because they
all died at `MAIN` D5 first. ⚠ **600 s was a review cap, ⛔ not proof of non-termination.**
⇒ ⭐ Determine whether it terminates. ⭐ Then either make it complete, or report it under §10 as an
unavailable construction — ⛔ do not leave it unmeasured, and ⛔ do not gate it by a count (see item 5).

## ⭐ 8 · DEAD CODE THE LAST ROUND EDITED

`nullspace_basis` (`:978-980`) has **no caller**; the live path is `generic_nullspace_vectors` (`:1245`).
⚠ Round 1's brief flagged this and it was edited rather than removed. ⭐ Remove it or give it a caller.
⚠ Minor: `exact_determinant:926` and `exact_domain_matrix:907` each call `reduced_matrix`, so
`factor(cancel(·))` runs twice per entry.
⚠ Latent, ⛔ not on the live `MAIN` path: `DomainMatrix` over the `EX` domain uses a `cancel`-based zero
test with no simplify fallback and can mis-rank an algebraically-zero unreduced entry; `XKIN_ANISO`'s
radical roots land in `EX`. ⭐ At minimum the stream should let a consumer tell the branches apart.

## Boundaries

- ⛔ No memory cap, no timeout, no handler that swallows a failure to make a run finish.
- ⛔ Do not change what is computed, `PACKAGE_DIMS`, the `D` range, or any package.
- ⛔ **No expected value and no acceptance criterion referencing one** (rule 5). ⛔ Do not treat any prior
  output as a reference — ⭐ a changed value is a **finding to report under §10**, ⛔ never something to tune
  away.
- ⛔ Do not add a registry, completeness map or exit policy.

## Acceptance

1. `/home/trevnorris/.s11_build/repro_d5.py` runs `run_cell('MAIN', 5)` to completion; report literal stdout.
2. A run whose declared cells complete **publishes `S11_exports.py`**, and one that fails **still emits its
   §10 report**. ⭐ Demonstrate both on a reduced `PACKAGE_DIMS` copy under `/tmp` — ⛔ never by editing the
   engine's own declarations.
3. `XKIN_ANISO` D3 roots 2 and 3 emit a non-empty basis, or §10 says why not.
4. ⛔ Do not run the full package loop; it has OOM-killed this machine.

## Deliverable

The edited script, and a note saying per item what changed, plus every value that moved anywhere and where
you reported it.
