---
unit_id: 164
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 164

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**

Two edits to `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl`, matching the directive exactly (git diff: 1 file, +51/-1):

- **Change A (banner):** `wl:32` changed from `banner["STAGE 147 — …"]` to `banner["STAGE 164 — …"]`. The transliteration-tell wrong stage number is corrected.
- **Change B (independent series route):** A new block inserted at `wl:118-166`, immediately after the healing-product asserts (`wl:115-116`) and before the `delta_perp on the healing-locked branch` section (`wl:168`), exactly where the directive placed it. It defines a multiplicative perturbation `pertRule` (`p -> p (1 + eps dlnP)` for `zq, csw, rhoW, tm, vw0, a, lw, cs`), forms `firstRatio`/`secondRatio` from the already-asserted `firstHealExpected`/`secondHealExpected` monomials, extracts the O(eps) coefficient via `Coefficient[Normal[Series[…,{eps,0,1}]],eps]` into `firstHealSeries`/`secondHealSeries`, and reconciles each against separately-typed `firstHealHand`/`secondHealHand` vectors (`wl:147-148`). It then builds `deltaPerpSeries = gstar*firstHealSeries + bWeight*secondHealSeries` and asserts it equals the compressed `A_*/B_*/C_*` target `deltaPerpSeriesExpected` (`wl:163`). A trailing `Clear[…]` (`wl:165-166`) restores a clean scope so the later `b`/`firstHeal`/`secondHeal`/`deltaPerp` section is untouched. The directive's instruction to use the separately-named `bWeight` (not clobbering the later `b`) and to leave the later section intact was followed.

No collateral edits: the git diff shows only the banner line and the inserted block; the SymPy script (`…_sympy_audit.py`) is unchanged.

**Assessment:**

Correct and it addresses the finding. The new route is a genuinely independent derivation: it obtains the log-channel coefficient vectors by `Series`/`Coefficient` perturbation of the explicit monomials — machinery with no counterpart in the SymPy script, which instead hand-writes `first_heal`/`second_heal` directly. So cross-engine agreement on the load-bearing coefficient vectors is no longer structurally guaranteed.

The new asserts are non-tautological. `firstHealSeries`/`secondHealSeries` are computed by Mathematica's series machinery from the monomial exponents (e.g. `csw^3`, `lw^2`, `a^2`), and compared against the independently-typed `firstHealHand`/`secondHealHand`. If any monomial exponent were wrong, the extracted coefficient (e.g. the `+3 dlncsw`, `-2 dlnLw`, `-2 dlna` terms) would differ from the hand vector and the residual would be nonzero. The `delta_perp via series route` assert compares the series-derived vectors (weighted by `gstar` and `bWeight`) against the hand-compressed `A_*/B_*/C_*` form, so it cross-checks both the series extraction and the compression algebra. The exec log confirms all three residuals are `0` with `PASS` lines (`first channel via series route`, `second channel via series route`, `delta_perp via series route`).

All pre-existing asserts remain in place and still PASS: parent-ratio identities, uniform products, healing products, `delta_perp compressed form`, and the five `expectApprox` numeric checks.

## Exec log assessment

**SymPy:** exit=n/a. The SymPy script was not modified by this directive (target was Mathematica-only), so no sympy exit is relevant to F1. The SymPy output `.txt` was nonetheless refreshed (mtime 2026-05-28 16:10) and continues to print the expected `= 0` residuals; the SymPy file uses `assert`-based `expect_zero` (no `PASS` text), and the run shows no AssertionError.

**Mathematica:** exit=0. Notable lines from the captured `.txt` (regenerated 2026-05-28 16:11):
- `STAGE 164 — MICROSCOPIC LOG-IMBALANCE CHANNELS` (banner now correct)
- `first channel via series route = 0` / `PASS: first channel via series route`
- `second channel via series route = 0` / `PASS: second channel via series route`
- `delta_perp via series route = 0` / `PASS: delta_perp via series route`
- Pre-existing PASS lines all present: `PASS: uniform first/second product`, `PASS: healing first/second product exact formula`, `PASS: delta_perp compressed form`, and `PASS:` for all five numeric checks (A_*, B_*, c_sw, v_w0, L_W) with diffs ≤ 4.4e-16.

**Output freshness:** confirmed. `.wl` mtime 2026-05-28 15:54:48; its `.txt` output mtime 2026-05-28 16:11:15 (newer than the script). SymPy `.txt` mtime 2026-05-28 16:10:12 (newer than the unchanged `.py`, mtime 2026-05-11). Both outputs are post-fix.

## Material-change assessment

`material_change`: false.

The edit added an independent verification route and corrected a banner label. It introduced no new constant and changed no derived monomial, coefficient vector, or numeric value (all asserts compare against the same already-asserted `firstHealExpected`/`secondHealExpected` and the same compressed target). The reconciled forms reproduce the pre-existing results bit-for-bit. No downstream unit's inputs change.

## Side observations (non-blocking)

- The SymPy script still carries the original `STAGE 147` banner at `…_sympy_audit.py:30` (visible in its `.txt` output as `STAGE 147`). The directive deliberately scoped the banner fix to the `.wl` only ("Edit only the Mathematica script"), and the original report's required change named only `wl:L32`. So the SymPy banner mislabel is out of scope for this directive and was correctly left untouched. It is a cosmetic codebase-wide renumbering artifact (offset ~17) noted in prior reports, not a math error. Flagging only so the user is aware the SymPy banner still reads 147 if banner consistency is later desired.

## Verdict justification

The single finding (F1) is fully resolved: both directive changes (banner 147→164 and the independent `Series`/`Coefficient` route) are in place exactly as specified, with no collateral edits. The new route is a real second-source derivation (machinery absent from the SymPy script) and its three new asserts are non-tautological, comparing series-extracted coefficient vectors against separately-typed hand forms and the compressed target. The Mathematica run exits 0 with all three new PASS lines plus every pre-existing PASS line, and outputs were regenerated post-fix. No regressions in the diff or log. Verdict: verified.
