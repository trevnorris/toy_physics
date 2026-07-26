# Rewrite feasibility measurement — shared dimension module for `scripts/ledger_stage*.py`

Measured 2026-07-26 on `ledger-v2-rebuild` @ `35dc6aa0`, tree clean. **Nothing modified or executed** beyond read-only greps and one classifier under `_scratch/pass1_dim_survey/measure/`. Paths relative to `research/pde_ledger_v2/`; `.wl` = `mathematica/ledger_stageNNN_<name>_mathematica_audit.wl`.

## 1. HEADLINE

**Dimension machinery = 3,752 / 39,511 lines = 9.5 %** (3,595 code-only, 9.1 %), in **30 of 43** scripts. **Coverage: COVERED 30 · PARTIAL 5 · UNCOVERED 8.** Two independent methods agree exactly: the 30 COVERED stages are precisely the 30 with non-zero `.py` dimension machinery; every PARTIAL/UNCOVERED stage has 0 dimension lines, so there is no attack surface there *today*.

**The `.wl` is NOT a sufficient gate for this rewrite.** It controls only *one-sided* edits — a declaration changes but its assertion target does not. The refactor's premise is that declaration and target both come from one module: the exact case the control cannot see. Only 9 of 30 `.out` files print a computed exponent triple; the other 21 reduce the whole dimensional firewall to a `PASS <TOOTH>` line. And the `.wl` never reads the `.py` — `scripts/run_all_audits.sh` and `mathematica/run_all_audits.sh` write to separate trees and **never diff**, so "reveal" means a human eyeballs two outputs.

## 2. REWRITE SIZE

**Criterion** (`measure/classify_dim_lines.py`, AST). A line is dimension machinery if it falls in a *dimension region*: (R1) the full span of any `class`/`def` named in the file's dimension-identifier set **D**; (R2) any module-level binding whose target/annotation/value references **D**; (R3) any other line referencing a **D** token at a word boundary. **D** is AST-seeded from `/dim/i` class/def names, module-level `_DIM`/`DIMENSIONLESS`/`Dim`/`Dimension` annotations, dimension-named locals (catches the stage032 `dim_A`/`dim_E` idiom) and the un-named aliases `dadd`/`dscale`/`homogeneity_residual`, then closed under one propagation round. **Excluded as non-physical** (verified by reading): `groebner_dimension`, `dimension_record` and the bare `dimension`/`dimension_status` locals in 040/043 — Krull dimension of an ideal; nullspace/`MatrixRank`/`RegionDimension` likewise.

**43 scripts, 39,511 total lines; dimension machinery 3,752 (9.5 %)**, code-only 3,595 (9.1 %). Per script: absolute min **0** · median **81** · max **351**; fraction min **0.0 %** · median **8.6 %** · max **45.5 %**. **13 scripts have zero** (001, 014, 015, 017, 019, 020, 022, 024, 025, 026, 028, 033, 043); across the 30 non-zero, median **117.5** lines / **12.4 %**.

**Heaviest five** (lines / % of file): `021` 351/45.5 · `007` 256/23.1 · `016` 247/27.3 · `005` 218/25.3 · `012` 202/19.3. By fraction: 021, 016 27.3, 004 26.4, 005 25.3, 002 23.7. **Lightest five non-zero**: `032` 9/0.7 · `039` 18/2.0 · `010` 22/3.1 · `041` 28/1.4 · `009` 53/8.6 (13 scripts sit lower still, at exactly 0).

**⚠ The line count understates the work — there is no single idiom to consolidate.** Base-dimension order and arity differ per script: **3-slot `(L,T,M)`** = 002, 004-007, 009, 030-032, 034-037, 044 (14); **3-slot `(L,M,T)`** = 011-013, 016, 018, 021, 023, 027 (8); **3-slot `(M,L,T)`** = 003, 010, 039-042 (6); **2-slot `(L,T)`** = 008 (`…stage008_…py:485-486  LENGTH = Dim(1, 0)`); **4-slot `(M,L,T,E-charge)`** = 038 (`…stage038_…py:692`). 021 is `(L,M,T)` internally, transposed to `(L,T,M)` for display (`…stage021_…py:225`). Carriers differ (`int`, `sp.Rational`, `fractions.Fraction`, `sp.Matrix`) and exponents are not integral — 042 `charge_dim = (1/2, 3/2, -1)` (`…stage042_…py:826-830`), 027 `I25: (5/2,0,0)`, 023 `gU: (-1/2,1/2,-2)`. **22 of 30 stages need their literals transposed**, and a transposition error is symmetric (§6.2).

## 3. `.wl` COVERAGE — **COVERED 30 · PARTIAL 5 · UNCOVERED 8**

Six independent read-only agents, one per stage band, judging from `.wl` content (never the title) plus `mathematica/out/*.out`. **Sanity check passes**: both purpose-built dimensional stages, 004 and 021, come out COVERED. `py` = dimension-machinery lines.

| stage | py | `.wl` justification (cited) |
|---|---|---|
|002|155|wl:42-44 `residualForUnits`; :54-73 l/t/mass triples; :76-131 9 residual keys; :390-395 `Scan[expectZero]`|
|003|81|wl:75 `dimResidual`; :613-635 8 base triples `(M,L,T)`; :637-660 18 checks; :661-672 3 ablations|
|004|175|wl:88-96 `unitVector` via **`UnitDimensions`**; :98-100 `expectDim`; :128-143 primitives; :209-228 8 asserts; :296-325 14 rows|
|005|218|wl:119-131 `unitVector`/`dimResidual`; :155-207 dictionary (`rho4 = -4 ell`); :209-305 32 rows; :327-332 count guard|
|006|162|wl:94-122 `dimAdd/dimSub/dimScale/dimPow`; :155-161 derived densities; :195-217 cited absolute triples|
|007|256|wl:98-127 dim algebra; :290-370 firewall (:310 `dK = M L^18 T^-2`, :362 `[rho_br]=M L^-3`); :874-875 drift table|
|008|77|wl:447-473 `expressionDim` walker + `kernelDimensionResidual`; :479-490 `{a}={1,0}` + per-ℓ asserts; :500 corrupt-`[a]` `expectFail`|
|009|53|*weak* — wl:426 `dimResidualVec`; :436-449 constant triples + z/tau/Z; :456-459 4 `expectZero`. Constants only, booleans only|
|010|22|*weak* — wl:490-510 `dimResidualVec`, `dimRho={0,-3,0}`, `dimGreenZeroStatic`, 5 `expectZero`. Constants only|
|011|151|wl:116 `dimensionalOk`→`FAILDIMENSIONAL`; :207-230 `dimOf` walker; :239-251 rules; :253 dim of **live** `cSSquaredBulk`; :406-419|
|012|202|wl:140 dim-gated verdict; :215-238 `dimOf` + dimensionless-fn guard; :250-285 `[K]→[c_S²]→[c_S]→[k]→[Z00]` chain; :562-584|
|013|162|wl:95-121 gated `computeVerdict`+`dimOf`; :244 `dimRules`; :248-260 dims of **live** `MEntries`/`KEntries` + corrupt-`[Tw]`; :442-463|
|016|247|wl:23 `zeroDim`; :90 `dimResidualVec`; :128-152 `dimOf`; :228-241 `makeDimRules`; :243-287 `evalDimensional`; :481-497 9 asserts; :293-300 corruption probes|
|018|101|wl:84-89 `dimL/dimSpeed/dimT2/T4/T5`; :102-121 `dimOf`; :228-239 rules+gate; :385-395 asserts + `expectFail`|
|021|351|wl:83 "KEEP-NATIVE"; :84-107 `dimOf`; :237-248 `rawDims` (`G→{3,-1,-2}`); :265-291 gate + 4 corruption maps. `.out:18-19` asserts `[N0]=(-1,1,0)`, `[D0]=(-1,1,-2)` **by value**|
|023|140|wl:209-235 `dimOf`; :237-245 21 base triples (`gU→{-1/2,1/2,-2}`); :247-252 7 targets; :254-287 audit + corruption; :607-624 gate|
|027|81|**strongest** — wl:183-190 `baseDims` (`gravG→{3,-1,-2}`, `I25→{5/2,0,0}`); :192-208 `unitDimension` by a *different algorithm* (monomial rescaling); :439 assert; **`.out:115` prints the computed triple** `[N0_den]={-1,1,0}`|
|030|98|wl:135 `dimResidual`; :168-171 12 triples; :218-370 gated asserts; :503-528 8-term firewall; :97-121 both mutation teeth recompute|
|031|124|wl:119 `dimResidual`; :445-468 21 `unitRows`; :470-485 loop + 6 live chains. **20 of 21 rows are declaration-vs-itself** (§6.3)|
|032|9|wl:836-842 `dimL={1,0,0}; dimE={2,-2,1}; dimA=dimE+dimL`, then derived `[A]/[U]/[F]` vs literals|
|034|109|wl:294-322 `unitScaleRules` + `dimensionUnderScaling` (log-derivative exponent extraction — independent method); :350-352, :683-687, :721-728|
|035|115|wl:406-443 `scaledRatio`/`unitScalingObject` (9 keys); :794-814 `expectedScaling` + assert; :490-491 mirrored inventory gate|
|036|113|wl:233-297 `unitRules`/`scalingResidual`/`dimensionResiduals` over the **live** Darwin kernel; :653-658 assert|
|037|133|**strongest** — wl:610-633 `unitRules`; :636-641 `scalingResidual`; :643-706 21-entry `unitContract` + `unitViolations`; :1268-1275; applied to `.wl`-derived expressions|
|038|69|wl:670-673 4-vector algebra; :675-719 `unitState` (`muDim={2,1,-4,-1}`); :721-730 expected; :1092-1103 assert|
|039|18|*narrow* — wl:408-424 `lengthDimension`/`multiplyDimensions`/`unitState`; :776-787 asserts only an **equality** `[curl u_T]==[b_T]`|
|040|59|wl:804-861 `dimensionState` (11 symbols + derived `K`, lockA/lockB/stability homogeneity); :1621-1628. (`.out:50 DIMENSION=10→8` is **Krull codimension**, not units)|
|041|28|wl:1197-1216 `dimensionState`; :2075-2088 asserts **absolute** triples `{{1,-3,0},{1,-3,0},{1,-1,-2}}`|
|042|123|wl:819-908 `dimensionGuard` (`chargeDim={1/2,3/2,-1}`, `physicalPower={1,2,-3}`); :1845-1853; :262-276 `speedFalloff` prints a dimension-derived **number**|
|044|120|wl:482-587 `dimensionalTooth` — 27 symbols, 24 integrand dims + measures, `bad = Select[actionDims, # =!= action &]`; :585 assert|

**PARTIAL (5) — all 0 `.py` dimension lines. The `.wl` independently restates dimension-*dependent* expressions, so monomial-exponent drift shows but no declaration error can.** `020` wl:93-96 `aPower` (single symbol `a`) + :237-246 `scalingOk`, :467-483 — slips: any `[G]`/`[c_s]`/`[c]` declaration (none exists either side); wl:812 "dim-like name present" is a substring check on identifiers, satisfiable by renaming. `024` wl:141 restates `gBase = Sqrt[rhoEff] cS^2 I25 xiQ/a^(7/2)` via a different route (:219-228). `025` wl:119-121 restates `citedN0Den`; the rest is a provenance-tag graph. `026` wl:92-94 `citedN0Den`, :203-210 `aPower === -7/2` gate (never printed to `.out`). `028` wl:99-102, :232-263 restate `Kbar_n` monomials + INV5 invariant residuals; no declaration table either side.

**UNCOVERED (8) — all 0 `.py` dimension lines.** `001` no dim content either side. `014` fully nondimensionalised; `.wl` "rank" is `MatrixRank`. `015` wl:196,:201 `MatrixRank`. `017` wl:34 `cited016DimensionalOk = True`, asserted wl:560 — a literal. `019` units-free abstract algebra; ω powers are polynomial degrees. `022` wl:588 asserts only `liveNames === {z}`. `033` wl:184-240 `NullSpace`/`MatrixRank` (Dirac constraint counting). `043` both sides print `DIMENSIONAL_HOMOGENEITY=N/A_INTEGER_COUNT_STAGE` (py:1193 / wl:1224); all "dimension" there is Krull/CAD.

## 4. WHAT THE `PASS` GATE IS WORTH

- **Verified: 16 of 43 emit no `PASS tally:` line** — exactly 001, 002, 003, 016–028. Matches `baseline/BASELINE_NOTE.md:26`; recomputed independently from `baseline/logs/`.
- **The ordered marker list is NOT stable.** Run 1 vs run 2 of the recorded baseline differ: `ledger_stage007_…` swaps two adjacent lines — `PASS  post-D16 function subcount computed from survivor enumeration` moves from position 117 to 118. All 42 other logs are byte-identical. An *ordered* gate would fire a spurious failure on 007; a **sorted/multiset** gate is stable. Counts were stable in both runs.
- **Markers are not unique**: stage013 emits 5 duplicate `PASS` strings (79 lines, 74 distinct) and stage042 emits 1 — a pure *set* comparison silently loses those 6.
- **Would a script that stopped checking early still produce a plausible count? Yes.** A failing check raises and exits non-zero, so the gate catches a check that *fires*; it does not catch a check that becomes **vacuous** — a comparison of a value against itself still prints `PASS` and keeps the count. Only **2 of 30** dimension-bearing scripts guard their own dimensional check arity: `…stage004_…py:502` (`… - 17`) and `…stage005_…py:843-846` (21/5/11/7). The other 28 have none.

## 5. COUPLING RISK

- **Filesystem writes: exactly one** — `…stage044_…py:1343-1346` writes `_scratch/stage044/verdict_py.json` and creates the directory. Corpus grep for `open(|write_text|write_bytes|json.dump|mkdir|makedirs` found no other writer.
- **Filesystem *reads*: one latent path** — `…stage007_…py:41-45` builds `SCRIPT_PATH`/`REPO_ROOT`/`REPORT_ROOT` → `software/stage1_solver/reports/pathA_24_T0_freeze.md`, `pathA_35_G0_freeze.md`. A shared module that relocates the script shifts `parents[3]`; same hazard at 044:1343.
- **No cross-script imports.** All 43 are standalone (`sympy`, stdlib, `numpy`/`scipy` in one). "Consumed from ledger_stageNNN" is prose inside `print()`, not an import — values are hand-transcribed, so a shared module cannot propagate a change automatically.
- **`ACTIVE_MUTATION` at module scope: 15 scripts** (030–044). Each sets `MUTATION_ENV = "LEDGER_STAGENNN_MUTATION"` then `ACTIVE_MUTATION = os.environ.get(...)` at import time (e.g. `…stage040_…py:28-29`). A shared module imported before those lines cannot see the mutation; moving the read later, or caching across it, silently disarms the C7 mutation harness. Also `BUILD_LOG: list[str] = []` (`…stage037_…py:242`, appended :251-252, read into the verdict :1094 — order-sensitive), and 040/043 wrap `sat_status`/`dimension_record`/`provenance` in `@lru_cache` keyed on arguments but **not** on `ACTIVE_MUTATION`.
- **Suite runtime ≈ 3–3.5 min sequential — ESTIMATE, not a measurement.** From mtime spans of the recorded baseline logs (suite not run): run 1 spans 181.7 s across 43 logs, run 2 spans 193.7 s; dominated by stage033 ~82 s, then 022 ~21 s, 018 ~11 s, 016 ~10 s. The `.wl` control costs far more — its runner sleeps 10 s between scripts and needs a Mathematica seat.

## 6. WHAT WOULD **NOT** BE CAUGHT

1. **Symmetric edits — the central hole, and exactly what this refactor produces.** In 21 of 30 COVERED stages the `.out` prints only `PASS <TOOTH>`; computed exponents are never rendered, so those stages diverge *only* when the `.py`'s declared dimension and its assertion target disagree. A shared module supplies **both** from one place, so a module error moves them together: `.py` stays green, `.out` stays byte-identical, and the `.wl` — holding its own hand-typed copy — silently disagrees with no observable. Concretely: 023 (`SOURCED_DIMS` + `EXPECTED_DIMS` rescaled together), 002/003/005/006, 009/010, 030 (a primitive plus `BULK_DENSITY`/`BRANE_DENSITY`), 038 (`unit_state` + `EXPECTED_UNIT_STATE`), 044 (a carrier plus its measure).
2. **A transposition error is symmetric by construction.** 22 of 30 stages need literals transposed to a common order (§2); a wrong permutation rewrites declaration and target identically. 008 (2-slot) and 038 (4-slot) additionally need arity changes.
3. **Declaration-vs-itself rows survive anything.** `…stage031_…py:553-580` / `wl:448-468`: 20 of 21 `PASS_UNITS_*` rows compare an expression to a second copy of itself — under the refactor both copies read one shared constant, so the row cannot fail for *any* value. Stage 017 is already there on both sides: `…stage017_…py:62 CITED_016_DIMENSIONAL_OK = True`, asserted py:814, mirrored wl:34 / wl:560 — `assert True` in two languages. Same at `…stage034_…wl:327-339` vs `:673-682`, and 044's `Z_chi` master row (`py:1038-1039` / `wl:1102-1110`).
4. **Quantities that cancel out of a homogeneity net.** 040: `[rho]` enters `dims["K"]` at `-4` and lock_A at `+4`, so `[rho] = M L^-3 → M L^-4` is invisible in *both* files; same for `[m]`, `[mu_R]`. 042: `[omega]`/`[d]` cancel out of `ratio_dim`. 039 asserts only `[curl u_T] == [b_T]`, so scaling `LENGTH` and `INVERSE_LENGTH` together passes. 044 tolerates compensating carrier shifts within a term.
5. **The 13 zero-machinery stages have no control at all** — not at risk today, but if the refactor wires them in their `.wl` files contain nothing to diverge from (022 = symbol-name sets, 025/026 = provenance-tag graphs, 033 = `NullSpace`/`MatrixRank`, 043 = `N/A_INTEGER_COUNT_STAGE`). **020 is the sharpest**: units-bearing (`a`, `c_s`, `c`, `G` all live) with a `.wl` that power-counts in `a` alone.
6. **A check made vacuous keeps its `PASS`.** 28 of 30 dimension-bearing scripts have no arity guard, so dropping a dimension-dictionary entry (shrinking a loop) or turning a comparison into an identity keeps exit 0 and, in the identity case, an unchanged count.
7. **The `PASS`-count gate cannot see any of 1-6.** All preserve exit 0, the tally, `PASSlines`, and the marker strings.

### Commands run
`git status --short`; `git log --oneline -1`; `ls`/`wc -l` over `scripts/` and `mathematica/`; `grep` surveys for `class Dim|Dimension = tuple|dim_of|dim_add|dim_scale|_DIM|UNITS|MUTATION|inspect|hashlib|open(|write_text|mkdir|lru_cache|PASS tally|count is`; `diff` of `baseline/logs` vs `baseline/logs2`; `python3 measure/classify_dim_lines.py`; `python3` one-liners over the baseline logs for tally/uniqueness/mtime statistics. Six read-only sub-agents did the `.wl` pass (bands 001-007, 008-014, 015-021, 022-028, 030-036, 037-044). **The 43-script suite was not run.**
