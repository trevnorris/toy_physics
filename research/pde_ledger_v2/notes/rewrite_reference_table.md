# Rewrite reference table — 43 SymPy stage scripts onto a shared dimension module

Measured **2026-07-26** on `ledger-v2-rebuild` @ **`adcfbdfd`** ("ledger v2: shared dimension module v0,
proven on stage004 + stage011"), working tree clean. Paths relative to `research/pde_ledger_v2/`.
**Corpus read-only**: nothing under `scripts/`, `mathematica/` or `manifests/` was modified. This file is
the only artefact written.

> ⚠ **State change during measurement.** The session opened at `89bcccdb` with stage004/011 mid-rewrite;
> the build committed as `adcfbdfd` partway through. **stages 004 and 011 are DONE** — already on
> `scripts/ledger_dimensions.py`. 28 stages remain, not 30. Their rows below are the *post*-rewrite state
> and double as the worked example for every other stage.

---

## 0. Headline numbers (all re-measured, not inherited)

| measure | value |
|---|---|
| scripts | 43 (`scripts/ledger_stage*.py`; there is no stage 029) |
| total lines | **39,444** |
| dimension machinery | **3,675 lines = 9.3 %** (code-only 3,534 = 9.0 %) |
| scripts with dimension machinery | **30 of 43** — confirms the prior count |
| already rewritten | **2** (004, 011) → **28 to go** |
| distinct bases | `(L,T,M)`×14 · `(L,M,T)`×8 · `(M,L,T)`×5 · `(L,T)` 2-axis ×1 (008) · `(M,L,T,E-charge)` ×1 (038) · `(stiffness,L,T)` ×1 (042) |
| scripts with fractional exponents | **8** — confirms the prior count (004, 005, 006, 012, 021, 023, 027, 042) |
| scripts with no `PASS tally:` line | **16** — exactly 001–003 and 016–028, confirms the prior count |
| `.out` files that render a computed dimension value | **9** — 004, 011, 012, 013, 016, 018, 021, 023, 027 |
| `.py` transcripts that render a computed dimension value | **8** — 011, 012, 013, 016, 018, 021, 023, 027 |

⭐ **The two sets differ by exactly one stage, and it is 004** — the completed pilot. Its `.wl` gained 20
`DIM|axes=L,T,M|name=…|exponents={…}` rows (`out:20-39`); its `.py` emits **zero** exponent triples. See §5.1.
*(Verified by a tight triple/monomial regex over fresh transcripts; hardcoded prose like `PASS target
[rho_br]=M L^-3` in 005/007/034/038/039 is a marker label, not a rendered value, and is excluded.)*

**Criterion for "has dimension machinery"** (column 2): the committed AST classifier
`_scratch/pass1_dim_survey/measure/classify_dim_lines.py`, re-run unmodified. A line counts if it is in a
*dimension region*: the full span of any `class`/`def` in the file's dimension-identifier set **D**; any
module-level binding referencing **D**; or any other line with a **D** token at a word boundary. **D** is
AST-seeded from `/dim/i` names plus the un-named aliases `dadd`/`dscale`/`homogeneity_residual`, and
excludes Krull/Gröbner/nullspace "dimension" (`groebner_dimension`, `dimension_record` in 040/043,
`MatrixRank`). **I agree with the prior 30-of-43.** The 13 zeros are 001, 014, 015, 017, 019, 020, 022,
024, 025, 026, 028, 033, 043 — the same 13 the `.wl` survey independently classified PARTIAL(5)/UNCOVERED(8).

**Line/percent deltas vs the prior measurement** are confined to 004 (175→128) and 011 (151→121): the
rewrite shrank them. Total moved 3,752→3,675 for the same reason. No other stage changed.

---

## 1. The per-stage table

Columns: **dim?** = has dimension machinery · **lines / %** = dimension machinery · **idiom** = how a
dimension is carried · **frac** = fractional exponents · **PASS** = real marker count (lines beginning
`PASS ` other than `PASS tally:`) · **tally** = the value on the `PASS tally:` line, `—` = no such line ·
**.out** = `mathematica/out/ledger_stageNNN_*.out` line count + sha256 first 12.

⭐ **The raw-grep caveat is confirmed and is exactly +1, never more.** `grep -c '^PASS'` also matches
`PASS tally:` (which begins `PASS ` too, so `grep -c '^PASS '` does not help). For all **27** tally-emitting
stages, `real == tally` exactly; raw = real + 1. For the 16 no-tally stages, raw = real. **No stage
disagrees by anything other than 1.**

| # | dim? | axis order — locus | lines / % | idiom | frac | PASS | tally | .out lines / sha256 | hazards |
|---|---|---|---|---|---|---|---|---|---|
|001|no|n/a — no dimension machinery|0 / 0.0 %|n/a|no|11|—|28 / `c1efe6d7edf4`|no tally line|
|002|yes|`(L,T,M)` — `…stage002_…py:199-201` (`l_unit=dim(1,0,0)`, `t_unit=dim(0,1,0)`, `m_unit=dim(0,0,1)`)|155 / 23.7 %|type alias `Dim = tuple[Rational×3]` (`:22`) + lowercase `dim()` factory + free `dim_add/dim_sub/dim_scale/dim_neg`|no|25|—|54 / `bd86cc022a2f`|no tally; order **never stated in prose** — only inferable from the three unit constants|
|003|yes|`(M,L,T)` — `…stage003_…py:866-880` (`d_rho_br=dim(1,-3,0)`=M L⁻³, `d_grad=dim(0,-1,0)`, `d_dt=dim(0,0,-1)`)|81 / 8.5 %|same as 002 (alias `:18` + `dim()` factory)|no|84|—|119 / `b4b3907c98f7`|no tally; order never stated; `.wl` is also `(M,L,T)` (see §5.3)|
|004|**yes — DONE**|`(L,T,M)` — `…stage004_…py:103` `DimensionBasis("L","T","M", render="symbolic")`|128 / 20.2 %|**shared module** `ledger_dimensions.DimensionBasis`; `Dim = DIMENSION_BASIS` (`:104`)|**yes** (¹⁄₂ powers, integer results)|49|49|133 / `c29a167b0472`|2 preserved tautologies (`:302` compares `d["hbar"]` to `d["action"]` — same object; `:308` RHS reduces to LHS). `.py` emits **no** value lines|
|005|yes|`(L,T,M)` — `…stage005_…py:183-184` (`class Dim` docstring "in {L, T, M} order")|218 / 25.3 %|frozen dataclass `Dim(l,t,m)` + `__mul__`/`__truediv__`/`__pow__`|**yes** (¹⁄₂ powers, integer results)|80|80|173 / `070f50cb1565`|one of only **2** stages with a dimension-check **arity guard** (`:843-846`)|
|006|yes|`(L,T,M)` — `…stage006_…py:115-116`|162 / 18.1 %|frozen dataclass `Dim(l,t,m)`|**yes** (`:468` ¹⁄₂ power, integer result)|121|121|193 / `3d565aa5e616`|—|
|007|yes|`(L,T,M)` — `…stage007_…py:146-147`|256 / 23.1 %|frozen dataclass `Dim(l,t,m)`|no|142|142|272 / `af0dfca6c63d`|⚠ **reads the repo filesystem** `:41-45` (`SCRIPT_PATH.parents[3]` → `software/stage1_solver/reports/pathA_24_T0_freeze.md`, `pathA_35_G0_freeze.md`) and SHA-256s them (`:273`) — relocating the script breaks it. ⚠ **marker ORDER is not run-stable** (two adjacent lines swap between runs); only a *multiset* gate holds|
|008|yes|**`(L,T)` 2-axis** — `…stage008_…py:461-462` (`class Dim` with only `l`,`t`), constants `:485-486` (`LENGTH=Dim(1,0)`, `TIME=Dim(0,1)`)|77 / 11.9 %|frozen dataclass `Dim(l,t)`|no|53|53|173 / `e83d7bac4713`|arity change, not transposition. Omitting M is **not defined to mean M=0** anywhere — needs an explicit policy|
|009|yes|`(L,T,M)` — `…stage009_…py:442-443`|53 / 8.6 %|frozen dataclass `Dim(l,t,m)`; note the named constant `MOMENT0 = Dim(0,-1,0)` = T⁻¹ (`:467`)|no|48|48|149 / `e0a30d5c8d13`|`M0` conflict source (§3.9)|
|010|yes|`(M,L,T)` — `…stage010_…py:556-571` (`dim_rho=[0,-3,0]`, `dim_cS=[0,1,-1]`; prose `:571` "dimensions ordered M,L,T")|22 / 3.1 %|**bare `sp.Matrix` column vectors** — no `Dim` type at all|no|71|71|184 / `ca14e6eaeb96`|lightest real machinery; `M0` conflict source (§3.9)|
|011|**yes — DONE**|`(L,M,T)` — `…stage011_…py:145` `DimensionBasis("L","M","T", render="tuple")`|121 / 17.2 %|**shared module**; stage-local `dim_of` SymPy walker retained|no|60|60|159 / `543c805b64ae`|walker deliberately **not** generalised into the module|
|012|yes|`(L,M,T)` — `…stage012_…py:192-193` (docstring "in (L, M, T) order"), printed `:807`|202 / 19.3 %|frozen dataclass `Dim(l,m,t)` + `dim_of` + `arg_dims` + `build_dimensional_block`|**yes** (`:580` ¹⁄₂ scale, integer result)|82|82|198 / `9f4eb95f1849`|closest sibling of the completed 011|
|013|yes|`(L,M,T)` — `…stage013_…py:175-176`, printed `:608`|162 / 20.2 %|frozen dataclass `Dim(l,m,t)` + `build_dimensional_block`|no|78|78|191 / `86df62922c91`|⚠ **5 duplicate PASS markers** (78 lines, 73 distinct) — a *set* comparison silently loses them; `K_eta`/`μ_η`/`T_w` conflict source (§3.3, §3.8)|
|014|no|n/a|0 / 0.0 %|n/a|no|93|93|223 / `1a19ba4a3e1a`|only stage importing `numpy`/`scipy` (`quad`, `eigh`) — float path|
|015|no|n/a|0 / 0.0 %|n/a|no|95|95|196 / `4ae6b0d6d80c`|register claims `Vp0/ℓ_c` and `{κ,χ,σ_a,σ_L}` are verified here; **neither engine computes them** (§3.6)|
|016|yes|`(L,M,T)` — `…stage016_…py:157-158`, printed `:713`|247 / 27.3 %|frozen dataclass `Dim(l,m,t)` + `DimError` + `make_dim_rules` + `dimension_eval`|no|82|—|176 / `8812f1db619d`|no tally; hashes expression strings (`:288`); **3-way `K_eta` and line-vs-volume `μ_η`/`T_w` conflict source** (§3.3, §3.8)|
|017|no|n/a|0 / 0.0 %|n/a|no|118|—|226 / `2294c0e5507a`|⚠ `CITED_016_DIMENSIONAL_OK = True` (`:62`) asserted at `:814`, mirrored `wl:34`/`:560` — `assert True` in two languages (§3.2)|
|018|yes|`(L,M,T)` — `…stage018_…py:162-163`, printed `:584`|101 / 14.1 %|frozen dataclass `Dim(l,m,t)` + `DimError`|no|59|—|145 / `74b1a69a6d6d`|no tally; `@lru_cache` `:274`|
|019|no|n/a|0 / 0.0 %|n/a|no|18|—|101 / `a08f907a21ce`|no tally|
|020|no|n/a|0 / 0.0 %|n/a|no|74|—|214 / `59c8600b4a26`|no tally; `inspect.signature` `:1159`. **Sharpest zero**: `a`, `c_s`, `c`, `G` are all live but nothing declares a dimension|
|021|yes|⭐ **TWO ORDERS.** internal `(L,M,T)` — `…stage021_…py:162`, `:557-558`, `:742`; **display `(L,T,M)`** — `:227` and `:238`|**351 / 45.5 %** (heaviest)|type alias `Dim = tuple[…]` (`:143`) + free `dim_add/dim_sub/dim_scale` + recursive evaluator `dim_of` + `TrackingDimensions`|**yes — genuine fractional VALUE** (`[mu_hat0] = (-1,-1/2,-1)`)|42|—|134 / `c48713f5cb94`|see §2 — three renderings, two conventions, in one record|
|022|no|n/a|0 / 0.0 %|n/a|no|55|—|134 / `32ce234f77a5`|no tally; **4× `@lru_cache`** (`:201`, `:244`, `:251`, `:300`); ~21 s runtime|
|023|yes|`(L,M,T)` — `…stage023_…py:171`, banner `:868`|140 / 13.1 %|type alias `Dim = tuple[Rational×3]` + `SOURCED_DIMS`/`EXPECTED_DIMS` dicts|**yes — literal** `gU`,`gW` = `(-1/2, 1/2, -2)` (`:473-474`)|111|—|214 / `f6e9c49207f5`|no tally; `@lru_cache` `:650`; `SOURCED_DIMS`+`EXPECTED_DIMS` are the textbook symmetric-edit pair; conflict source for `K_eta`, `T_Omega`, `M0`, `Z0_ret/Z1_ret` (§3.3, §3.5, §3.9)|
|024|no|n/a|0 / 0.0 %|n/a|no|7|—|64 / `2efd19431a39`|no tally; fewest markers in the corpus|
|025|no|n/a|0 / 0.0 %|n/a|no|18|—|74 / `13a97596e62a`|no tally; `inspect.signature` `:288`|
|026|no|n/a|0 / 0.0 %|n/a|no|27|—|107 / `7a4cbec22f2c`|no tally; `inspect.signature` `:354`|
|027|yes|`(L,M,T)` — `…stage027_…py:199`, printed `:793` ("dimension tuple convention: (L,M,T)")|81 / 9.4 %|type alias `Dim = tuple[Rational×3]` + `BASE_DIMS` dict + a second `unit_dimension` by monomial rescaling|**yes — literal** `I25 = (5/2, 0, 0)` (`:208`)|55|—|138 / `4c826275e42c`|no tally; `@lru_cache` `:328`; `inspect.signature` `:484`. **Strongest `.wl`** — independent algorithm, and `out:115` prints the computed triple|
|028|no|n/a|0 / 0.0 %|n/a|no|26|—|132 / `8be02c6136b4`|no tally; `inspect.signature` `:298`|
|030|yes|`(L,T,M)` — `…stage030_…py:172-173` (docstring "in [L,T,M] order"), `__str__` `:197-198` renders `[l,t,m]_[L,T,M]`|98 / 14.8 %|frozen dataclass `Dim(l,t,m)`|no|16|16|53 / `87fc8be638ce`|`ACTIVE_MUTATION` at `:107` (import time)|
|031|yes|`(L,T,M)` — `…stage031_…py:140-141`|124 / 19.5 %|frozen dataclass `Dim(l,t,m)`|no|50|50|79 / `1828ea5fc62f`|`ACTIVE_MUTATION` `:94`. ⚠ **20 of 21 `PASS_UNITS_*` rows compare an expression to a copy of itself** (`:553-580`) — cannot fail (§3.1)|
|032|yes|`(L,T,M)` — `…stage032_…py:1046-1047` (`dim_L=(1,0,0)`, `dim_E=(2,-2,1)`=L²T⁻²M)|**9 / 0.7 %** (lightest non-zero)|**bare inline tuples**, no type, no helpers — 9 lines inside one function|no|57|57|91 / `4902fb2de895`|`ACTIVE_MUTATION` `:28`; three SHA-256 manifest digests (`:694`, `:807`, `:1095`) — check whether any dimension row feeds them before changing a rendering|
|033|no|n/a|0 / 0.0 %|n/a|no|33|33|100 / `96d3046a9182`|`ACTIVE_MUTATION` `:25`; `random.Random(seed)` `:656` (**seeded**, deterministic); **slowest script ≈ 82 s**|
|034|yes|`(L,T,M)` — `…stage034_…py:329-336`|109 / 12.3 %|type alias `Dimension = tuple[Rational×3]` + `dim_add`/`dim_scale`|no|12|12|70 / `1e3f59dc9263`|`ACTIVE_MUTATION` `:27`; manifest SHA `:599`; `wl:327-339` vs `:673-682` is a declaration-vs-itself pair|
|035|yes|`(L,T,M)` — `…stage035_…py:339` (comment) + `:342-345`|115 / 11.5 %|type alias `Dim` + `add_dims`/`scale_dim`|no|12|12|76 / `de6a2f4421ed`|`ACTIVE_MUTATION` `:26`; manifest SHA `:687`|
|036|yes|`(L,T,M)` — `…stage036_…py:282` (comment) + `:285-290`|113 / 12.4 %|type alias `Dim` + `add_dims`/`scale_dim`|no|8|8|64 / `61afbc7ccbb4`|`ACTIVE_MUTATION` `:25`; fewest markers of any dim-bearing stage; `A_E`/`A_V` conflict source (§3.7)|
|037|yes|`(L,T,M)` — `…stage037_…py:604` (comment) + `:607-618`|133 / 9.8 %|type alias `Dim` + 21-entry `UNIT_EXPECTATIONS` contract|no|16|16|95 / `94cc46061a05`|`ACTIVE_MUTATION` `:28`. ⚠ **order-sensitive global** `BUILD_LOG: list[str] = []` (`:242`), ordinal read at `:251`, appended `:252`, frozen into the verdict at `:1094` — any reordering of construction changes the verdict payload. `A_E`/`r_BA` conflict source (§3.4, §3.7)|
|038|yes|**4-axis `(M,L,T,E-charge)`** — `…stage038_…py:692` (`Dimension = tuple[int,int,int,int]  # M, L, T, E-charge`)|69 / 5.9 %|type alias, 4-tuple of `int`|no|10|10|66 / `0609d45625d4`|`ACTIVE_MUTATION` `:30`; `EXPECTED_UNIT_STATE` (`:741`) is a **positional, unlabelled 8-tuple of 4-tuples** compared against `unit_state()` (`:708`) in the same file — moves together under any transposition. Projecting to `{L,M,T}` by dropping `E` would turn inequivalent units into a false match. `μ_R`/`A_E`/`r_BA` conflict source (§3.4, §3.7)|
|039|yes|`(M,L,T)` — `…stage039_…py:433` (`Dimension = tuple[int,int,int]  # M, L, T`)|18 / 2.0 %|type alias, 3-tuple of `int`|no|10|10|61 / `e58bf0595b97`|`ACTIVE_MUTATION` `:26`; asserts only the **equality** `[curl u_T] == [b_T]` — scaling both together passes|
|040|yes|`(M,L,T)` — `…stage040_…py:799` + the `dims` map `:810-820` (`rho=(1,-3,0)`, `c_gamma=(0,1,-1)`)|59 / 3.5 %|type alias `Dim = tuple[int×3]` + `dadd`/`dscale` + a `dims` dict|no|28|28|61 / `3ed0aa076a53`|`ACTIVE_MUTATION` `:29`. **4× `@lru_cache`** — `sat_status` `:502`, `dimension_record` `:548`, `provenance` `:593`, `run_case` `:672`. ✅ **VERIFIED argument-pure**: none reads `ACTIVE_MUTATION`; the mutation is injected as the `mutate_c_e=` kwarg at `:1497` into the **uncached** `dimension_state()`. `dimension_record` is **Krull** dimension, not units. `[rho]`, `[m]`, `[mu_R]` cancel out of the homogeneity net|
|041|yes|`(M,L,T)` — `…stage041_…py:1084-1088` (`dim_mass=Matrix([1,0,0])`, `dim_bulk_density=[0,-4,0]`, `dim_speed=[0,1,-1]`)|28 / 1.4 %|**bare `sp.Matrix` column vectors**, converted to `int` tuples on return|no|35|35|72 / `f0516b9ee719`|`ACTIVE_MUTATION` `:25`|
|042|yes|⭐ **stiffness basis `(stiffness, L, T)`** — `…stage042_…py:808` + `:820-830` (`stiffness_dim=(1,0,0)`, `speed_dim=(0,1,-1)`, `length_dim=(0,1,0)`)|123 / 6.0 %|type alias `Dim = tuple[Fraction×3]` (**`fractions.Fraction`, not `sp.Rational`**) + `dadd`/`dscale`|**yes — literal** `charge_dim = (1/2, 3/2, -1)` (`:826-827`)|50|50|100 / `abbdef310e2f`|`ACTIVE_MUTATION` `:28`; **1 duplicate PASS marker**. Axis 1 is **not mass** — `B`,`C`,`K`,`M_h`,`d`,`r`,`Q_E` all carry different dimensions here than elsewhere. **No normalisation convention is stated anywhere**, including `manifests/DIM_ORDER_DECISION.md`. ⛔ **the `.wl` comment `:816` says "MLT" — a mislabel** for a basis carrying `{1/2,3/2,-1}`. `[omega]`/`[d]` cancel out of `ratio_dim`; guard runs twice with a mutation flag (§5.3)|
|043|no|n/a — prints `DIMENSIONAL_HOMOGENEITY=N/A_INTEGER_COUNT_STAGE` (`py:1193`, `wl:1224`)|0 / 0.0 %|n/a|no|20|20|55 / `18cf107d2a8c`|`ACTIVE_MUTATION` `:30`. ✅ **VERIFIED: zero `@lru_cache`** — the prior claim that 043 caches `sat_status`/`dimension_record`/`provenance` is **wrong**; those live only in 040|
|044|yes|`(L,T,M)` — `…stage044_…py:416` (`# Dimension order is [L,T,M].`), also as dict keys `"[L,T,M]"` at `:1012`, `:1039`, `:1090`|120 / 8.8 %|**bare `int` 3-tuples** as locals inside `dimensional_tooth()` + `add_dim`/`scale_dim` (`:163`, `:167`) + `"[L,T,M]"`-keyed dict rows|no|20|20|47 / `2930e8ca06db`|⛔ **the only filesystem WRITER**: `:1343-1346` `Path(__file__).resolve().parents[3]` → `research/pde_ledger_v2/_scratch/stage044/verdict_py.json`, with `mkdir(parents=True)`. Relocating the script moves the write target. `ACTIVE_MUTATION` `:31`. `Z_chi` master row (`:1038-1039` / `wl:1102-1110`) is declaration-vs-itself. Also frozen/paused pending the **044-v2 un-freeze**|

### Verified counts, cross-cutting

- **`ACTIVE_MUTATION` at import time: exactly 15** — 030–044 inclusive. Confirms the prior count. Each is
  `MUTATION_ENV = "LEDGER_STAGENNN_MUTATION"` then `ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()`
  at module scope. **13 read it at lines 25–31; two read it late** — 030 at `:107` and 031 at `:94`.
  ⛔ *A module import placed above those lines is fine; caching a value across them, or moving the read
  later, silently disarms the C7 mutation harness.*
- **`@lru_cache`: 5 stages, 11 sites** — 018 (1), **022 (4)**, 023 (1), 027 (1), **040 (4)**. The prior
  report named only 040 and 043; 043 has **none**, and 018/022/023/027 were missed. Of these, 018/023/027
  are dimension-bearing.
- **Filesystem writes: exactly one** — 044. **Filesystem reads: exactly one** — 007.
- **Cross-script imports: none.** No stage imports another stage. The only intra-corpus import is
  `from ledger_dimensions import …` in 004 (`:16`) and 011 (`:18`).
- **`inspect.signature`: 5 stages** — 020 `:1159`, 025 `:288`, 026 `:354`, 027 `:484`, 028 `:298`. Only
  **027** is dimension-bearing; the other four are zero-machinery. ⚠ *If a dimension helper ever becomes a
  parameter of an inspected function, moving it to a module changes the signature these stages assert on.*
- **Duplicate PASS markers: 2 stages, 6 lines** — 013 (5) and 042 (1). A *set* comparison loses all six;
  use a **multiset**.
- **Marker order instability: 1 stage** — 007. An *ordered* gate false-alarms; a multiset gate holds.
- **Arity guards on the dimension check: 2 of 30** — 004 `:460` (`… - 17`; the pre-rewrite locus `:502` has moved) and 005 `:843-846` (four counts: 21/5/11/7). The other 28 keep their
  `PASS` if a check is made vacuous.
- **SHA-256 of canonical manifest text: 6 stages** — 032, 033, 034, 035, 036, 037, 038. None hashes its own
  source, but any of them could carry a dimension rendering into a frozen digest; check before changing
  `__str__`/`repr` output.

### Fractional exponents — a distinction the prior count flattens

Eight stages, confirmed. But they split into two very different jobs:

| kind | stages | what the module must do |
|---|---|---|
| fractional **power operation**, integer result | 004, 005, 006, 012 | support `** Rational(1,2)` / `dim_scale(…, 1/2)`; stored exponents stay integral |
| fractional **stored value** | **021, 023, 027, 042** | store and render half-integers: `[mu_hat0]=(-1,-1/2,-1)` (021), `gU/gW=(-1/2,1/2,-2)` (023 `:473-474`), `I25=(5/2,0,0)` (027 `:208`), `charge_dim=(1/2,3/2,-1)` (042 `:826-827`) |

⚠ **On "four of them in `{L,T,M}`":** that phrasing does not survive checking. Three of the eight are
`(L,T,M)`-ordered (004, 005, 006). The **four** are the `(L,M,T)`-ordered ones — 012, 021, 023, 027 —
which is what `_scratch/dim_harness/HARNESS_DESIGN_REVIEW.md:69-72` actually says. 042 is the only one
outside the L/M/T axis set entirely.

### `scripts/output/*.txt` is a partial, not a baseline

`scripts/run_all_audits.sh` writes tracked transcripts to `scripts/output/`, but **only 001–028 exist**
(28 stage files + `midway_knob_audit_codimension_sympy.txt` + a one-line `_summary.txt` that records a
single stage028 run). **Stages 030–044 have no committed `.txt` at all.** I verified the 28 that exist are
current: stripping the `#` header and `EXIT_CODE:` footer, all 28 bodies are byte-identical to a fresh run
today. ⚠ *So the Python-side re-baseline step has an artefact for 001–028 and nothing for 030–044.*

---

## 2. ⭐ stage021 — the two-order case, and the only one

**Checked for others and found none.** Two independent scans:

1. Textual — only `…stage021_…py` mentions more than one axis order (`(L,M,T)`×6, `(L,T,M)`×1).
2. Structural — a scan for index-permuting renderers (`("L", dim[0]), ("T", dim[2]), ("M", dim[1])`) hits
   **only** stage021, at `:227` and `:238`. Every other stage renders through **named attributes**
   (`self.l`, `self.m`, `self.t`) or a `zip(("L","M","T"), components())`, so no positional permutation is
   possible.

Inside stage021 the situation is worse than "two orders" — one `dim_record` (`:246-255`) emits **three
renderings under two conventions**:

| field | function | convention |
|---|---|---|
| `dimension` | `dim_to_text` `:224-232` | **transposed to L,T,M** |
| `dimension_monomial` | `dim_to_monomial` `:213-214` (`L**d[0] * M**d[1] * T**d[2]`) | **native (L,M,T)** |
| `dimension_vector_LMT` | raw tuple | **native (L,M,T)** |

Plus `dim_to_pretty_text` `:235-243` (transposed) and `dim_vector_text` `:257` (native). This is live in the
transcript: the same quantity prints as `[mu_hat0]=L⁻¹ T⁻¹ M⁻¹ᐟ²` (`py:656`, L,T,M) and as
`"back-solved mu_hat0 has dimension (-1,-1/2,-1)"` (`py:670`, L,M,T) **eleven lines apart**. A rewrite that
binds one basis and renders through it will change one of those two strings unless both renderings are kept
explicitly.

---

## 3. The nine cross-stage conflicts, mapped to stages

From `manifests/PIVOT_TO_REWRITE.md` §3, with loci recovered from
`notes/measure_register_sufficiency.md:75-118`. **"First surfaces at"** is keyed to the rewrite order in §4.

| # | conflict | stages spanned | first surfaces at |
|---|---|---|---|
|1|**stage031: 20 of 21 `PASS_UNITS_*` rows compare an expression to a copy of itself** (`py:553-580`; `wl:447-470`). They cannot fail for any value.|031 only|**step B6** — late. Self-contained; adjudicate in place, do not repair inside the refactor|
|2|**stage017 asserts a hardcoded `True`** in both engines (`py:62`, `:814`; `wl:34`, `:560`). The values it claims live at 016.|017 (claim) ← 016 (values)|**step A3 (016)** — 016 is rewritten early; 017 has no machinery so it is never rewritten. Decide then whether 017's citation gets a real target|
|3|**`K_eta` carries three different dimensions**: `M L⁻¹T⁻²` (`013:412`) · `M L⁻³T⁻²` (`016:362`) · `M T⁻²` (`023:466`). The register documents 013-vs-023 only.|013, 016, 023|⭐ **step A2 (013)** — the earliest adjudication in the whole plan. All three stages are in the first group, so it is fully in view by **step A5 (023)**|
|4|**Two stages compute the same ratio `r_BA = q_T²c_γ²/(μ_R A_E)` in incompatible unit systems.** 037 uses `μ_R=M L⁻¹T⁻²`, `A_E=M L³T⁻²`; 038 uses `μ_R=M²L T⁻⁴E⁻¹`, `A_E=L·E`. **Both land dimensionless**, so nothing catches it.|037, 038 (036 shares 037's convention)|**step B12 (037)**, resolved at **step D3 (038)**|
|5|**`ε0/ε1` vs `Z0_ret/Z1_ret`** — the register (`:169`) calls these the *same two dofs*, but `ε0/ε1` are dimensionless `1` at 009 and `Z0_ret/Z1_ret` are `M T⁻²` at 023 (`py:162`).|009, 023|**step A5 (023)** if you take 023 in group A; otherwise **step B3 (009)**. Either way both halves are visible before group C|
|6|**Register locus mis-attribution.** Rows `:182-185` claim "stage 017, dual-engine verified" for `T_Ω`, `β₂`, `M̃/K̃/T̃_Ω`; 017 verifies nothing — the values are at `016:355-366`. Same shape at `:174` (`Vp0/ℓ_c`, attributed to 015, which no engine computes) and `:181` (`{κ,χ,σ_a,σ_L}`, "stage 015", no triple given).|017→016, 015 (nothing), 023 (`T_Omega`)|**step A3 (016)**. Not a physics dispute — a documentation repair, already queued as PIVOT §4 item 2|
|7|**`A_E`**: `M L³T⁻²` at `036:289` and `037:614` vs `L·E` at `038:715`. **No register row exists at all.** `A_V` collides the same way (register `:227` for 032 = `M L³T⁻²` vs `036:395`/`037:680` = `L²T⁻²`).|036, 037, 038 (+032 for `A_V`)|**step B11 (036)**, forced at **step D3 (038)**|
|8|**`μ_η` / `T_w` line-vs-volume convention clash**: `M L⁻¹`, `M L T⁻²` at `013:411-412` vs `M L⁻³`, `M L⁻¹T⁻²` at `016:360-361`. Register records only 013. Likely *legitimate* (line- vs volume-density) — needs a stated convention, not a fix.|013, 016|**step A2 (013)**, decided by **step A3 (016)**|
|9|**`M0`**: `T⁻¹` at `009:467` and `010:556` vs `M T⁻¹` at `023:460`. No register row arbitrates.|009, 010, 023|**step A5 (023)** — but 009 is B3 and 010 is C2, so this one is *half-visible* when 023 lands. ⚠ **Read 009/010 before adjudicating 023.**|

⚠ **There are more than nine.** `notes/measure_register_sufficiency.md:78-85` lists four further same-name /
different-dimension pairs that PIVOT §3 folds away: **`μ_R`** (`M L⁻¹T⁻²` at 003/007/037/040/044 vs
`M²L T⁻⁴E⁻¹` at 038), **`A_V`** (032 vs 036/037), **`T_Omega`** (`M L⁻³T⁻²` at `016:363` vs `M T⁻²` at
`023:467` — the register *does* carry both), and **`Q_E`** (register `:143` asserts *no* dimension;
`042:845` assigns `charge_dim`). Counting those, it is **13**. Budget for 13 adjudications, not 9.

⭐ **The single biggest one is not on either list**: stage042's whole basis. `B`, `C`, `K`, `M_h`, `d`, `r`,
`Q_E` are all `{L,T,M}` quantities elsewhere and stiffness-basis quantities here, and **no normalisation
convention is written down anywhere in the repo** — not in the register, not in
`manifests/DIM_ORDER_DECISION.md`. That is a design decision the module needs before 042 can be touched.

---

## 4. Suggested rewrite order

The requested order was `(L,M,T)` → `(L,T,M)` → `(M,L,T)` → awkward bases → 021 last. **The data
supports the grouping but wants two changes**, both for the same reason: the completed pilots already
de-risked *idiom families*, and the idiom families cut across the axis groups differently than expected.

**Finding 1 — the two pilots between them cover the entire frozen-dataclass family.** Pre-rewrite,
stage004 was `@dataclass class Dim` with fields `(l, t, m)` (`git show d9544a62:…:106-111`) and stage011
was the same class with `(l, m, t)` (`git show 89bcccdb:…:148-153`). Eleven remaining stages use that
identical class — 005, 006, 007, 009, 030, 031 in `(L,T,M)`; 012, 013, 016, 018 in `(L,M,T)`; 008 as a
2-field variant. **Those eleven are near-mechanical repeats of work already done.** 012 and 013 are the
closest of all: they share 011's `build_dimensional_block` / `arg_dims` / `dim_of` scaffolding
(`grep -l "def build_dimensional_block"` → 011, 012, 013 and nothing else).

**Finding 2 — the `(L,M,T)` group is *exactly* the group whose `.out` already renders values.** The nine
stages whose `.out` renders a computed dimension are 004, 011, 012, 013, 016, 018, 021, 023, 027 — i.e.
**all eight `(L,M,T)` stages plus 004**. Every `(L,T,M)` and `(M,L,T)` stage other than 004 needs a `.wl`
edit, an `.out` re-baseline, and therefore a **Mathematica seat**. Doing `(L,M,T)` first buys seven stages
with little or no `.wl` work. This is a stronger reason than convention-tidiness for the same ordering, so
**keep `(L,M,T)` first**.

**Change 1:** put **021 at the end of group A, not the end of everything.** It is `(L,M,T)`, its `.out`
already renders values, and its two-order problem is *internal to one file* — it does not interact with any
other stage. Deferring it to the very end means carrying the hardest single file past 26 unrelated rewrites
with no compounding benefit. Do it while the `(L,M,T)` conventions are still loaded.
*(If you would rather bank easy wins first, the fallback is to leave it last as originally planned — the
cost is only sequencing, not correctness.)*

**Change 2:** split the `(L,T,M)` group by **idiom**, not by stage number — dataclass siblings of 004
first, then the type-alias-tuple family, then the ad-hoc ones. The `(M,L,T)` group is small and mostly
ad-hoc, so it merges naturally into that tail.

### The order

**Group A — `(L,M,T)`, 7 stages. `.out` already renders values; 011 is the worked example.**

| step | stage | one-line reason |
|---|---|---|
|A1|**018**|101 lines, same dataclass as 011, `.out` prints a symbolic dimension (`out:77`) — smallest true repeat, proves the loop|
|A2|**013**|011's direct sibling (`build_dimensional_block`); `.out:86-87` already prints the sourced tuples; **forces the `K_eta` / `μ_η` / `T_w` adjudication first, while it is cheap**|
|A3|**016**|the other side of that adjudication plus the register's false-provenance rows; heaviest of the four clones at 247 lines|
|A4|**012**|202 lines, same scaffolding, adds the first fractional *power* op (`:580`); `.out:110` prints triples|
|A5|**023**|leaves the dataclass family for the tuple alias; first **fractional literal** (`gU/gW`); brings `M0` and `Z0_ret/Z1_ret` into view — ⚠ read 009/010 first|
|A6|**027**|tuple alias + `BASE_DIMS` dict + `I25 = 5/2`; the **strongest `.wl` in the corpus** (independent algorithm, `out:115` prints the computed triple) — the best available cross-engine check, so spend it on a fractional stage|
|A7|**021**|351 lines / 45.5 %, two internal orders, three renderings, a genuine half-integer value. Hardest file in the corpus — do it while `(L,M,T)` is fresh|

**Group B — `(L,T,M)`, 12 stages. Every one needs a `.wl` print + `.out` re-baseline.**

| step | stage | one-line reason |
|---|---|---|
|B1|**005**|direct dataclass sibling of 004, and one of only two stages with an arity guard (four counts, `:843-846`) — the loop's own regression test|
|B2|**006**|same dataclass; `.wl` has **13 top-level dimension globals**, the easiest print insertion in the corpus|
|B3|**009**|same dataclass, only 53 lines; establishes `M0 = T⁻¹` before 010 and before the 023 adjudication lands|
|B4|**030**|same dataclass; already carries an explicit `[L,T,M]` label in `__str__` (`:198`), so the rendering contract is pre-stated|
|B5|**002**|switches to the tuple-alias + lowercase `dim()` factory family; 155 lines; ⚠ its order is *never written down* — record it|
|B6|**031**|dataclass, but its 20 self-comparing rows mean the rewrite has **no working check** — do it after 030 so the idiom is settled, and adjudicate §3.1 rather than repair|
|B7|**007**|256 lines, dataclass, **plus** the filesystem read and the unstable marker order — the awkward one in this family|
|B8|**034**|`Dimension` alias variant; `.wl` `dimensionObject[]` returns the whole map — **the only clean stage004-shaped cross-check left in group B**|
|B9|**032**|9 lines of bare inline tuples, `.wl` has no `Module` at all (all globals) — trivially small and trivially checkable; sits inside a 1,221-line file with three manifest digests|
|B10|**035**|type alias + `add_dims`/`scale_dim`; ⚠ `.wl` yields **monomials, not vectors** — a print needs an `Exponent[…]` conversion|
|B11|**036**|same alias family; opens the `A_E` / `A_V` conflict. ⛔ `.wl` stores **no computed exponent vector** — cross-check would re-transcribe literals|
|B12|**037**|same family, 21-entry contract; **`BUILD_LOG` order-sensitivity** and the `r_BA` conflict. ⛔ `.wl` has **zero integer tuples** — no cross-check available without adding a derivation|

⚠ **Group B tail re-ordered by `.wl` checkability, not by size.** The original 035→036→037→034→032 run put
the three stages with **no usable `.wl` value surface** (035, 036, 037) before the two with a clean one
(034, 032). Reversed above: bank the checkable ones first, then take 035/036/037 knowing they are being
rewritten **without a cross-engine check** — and record that fact rather than discovering it mid-run.

**Group C — `(M,L,T)`, 5 stages.**

| step | stage | one-line reason |
|---|---|---|
|C1|**039**|18 lines, plain `int` tuples; asserts only an equality. `.wl` constants are top-level and `unitState` returns its vectors (`wl:423`) — the cheapest place to establish the `(M,L,T)` binding, and the only fully-checkable one in the group|
|C2|**010**|22 lines, `sp.Matrix` carrier; `.wl` already prints at `:511` so one inserted line renders all 11. ⚠ its `.wl` **declares `M,L,T` but no vector evidences it** (M slot is 0 throughout) — record the order, do not infer it. Completes the `M0` picture with 009|
|C3|**041**|28 lines; first `sp.Matrix` carrier — decide the Matrix→module conversion here. ⚠ `.wl` returns only 3 of 7 vectors and **the mutation-carrying `dimSpeed` is one of the missing four**|
|C4|**040**|59 lines + **4 `lru_cache`d functions** (verified argument-pure) + quantities that cancel out of the net; keep the Krull/`RegionDimension` `dimension_record` well away from the units path. ⚠ `.wl`'s 12 per-quantity vectors never leave `dimensionState`; only 7 lock/stability combinations escape|
|C5|**003**|81 lines, tuple-alias family, and the largest `(M,L,T)` job; `.wl` is also `(M,L,T)` (`wl:621`, `:627`, `:628`) so no cross-engine transposition — but **neither file states it**, and a print emitting `axes=L,T,M` here would corrupt every triple. Do it last in the group, when the convention is written down|

**Group D — the bases that are not three L/M/T axes, 4 stages. Each needs a decision before code.**

| step | stage | one-line reason |
|---|---|---|
|D1|**044**|`(L,T,M)` but with bare `int` tuples and `"[L,T,M]"`-keyed dict rows, **and the only filesystem write** — and it is separately frozen pending 044-v2. Its `.wl` has the corpus's **largest** computed dimension set (~95 vectors) and **exactly 1 is reachable** (`wl:1104`). ⚠ *Consider deferring D1 until 044-v2 lands, or the rewrite will be redone.*|
|D2|**008**|2-axis `(L,T)`. **Policy needed**: does omitting M mean M=0? The script does not say so, and asserting it silently is a physics claim|
|D3|**038**|4-axis `(M,L,T,E-charge)`. Module supports ≥4 axes structurally but the pilot **explicitly declined to claim four-axis semantics**. Projecting to `{L,M,T}` is forbidden without an explicit `E` conversion. ✅ **Good news**: its `.wl` `unitState` returns all 8 vectors as a live top-level `List` (`wl:709-718`), so the cross-check is a drop-in once the semantics are decided|
|D4|**042**|stiffness basis + `fractions.Fraction` carrier + fractional literals + **no normalisation convention written anywhere** + the `.wl`'s own basis comment is **wrong** (`wl:816` "MLT"). Its 14 base vectors are unreachable and the guard runs twice under mutation, so there is **no safe print-only cross-check**. The single largest open design question; also the stage the register misses entirely|

> **The pilot's own advice was to design against 042 and 038 *first*** (`PIVOT_TO_REWRITE.md:139`), because
> they are the extremes and they broke the register. The module was then built without them
> (`adcfbdfd`: "no claim yet for four-axis or stiffness-basis semantics"). **Resolve that tension before
> group A, not at D3.** Deciding the ≥4-axis and stiffness semantics on paper up front costs one session;
> discovering at D3 that the v0 `DimensionBasis` contract cannot express them costs a re-run of A–C.

---

## 5. ⛔ Where the per-stage loop will NOT work

The loop is: **(a)** add print-only output to the `.wl` → **(b)** re-baseline `.out` → **(c)** rewrite the
`.py` onto the module → **(d)** compare `.py` values against `.wl` printed values, axis-labelled.

### 5.1 Step (d) has no `.py` side — including for the completed pilot

⭐ **stage004's `.py` prints no dimension value at all.** Its `.wl` emits 20 `DIM|axes=L,T,M|name=…|exponents={…}`
rows (`out:20-39`), but a fresh `python3` run of the `.py` emits **zero** lines matching `DIM|` or any
exponent triple. The 29-value comparison in `_scratch/dim_harness/PILOT_004_011_REPORT.md:34` was produced
by an **external comparison script that imported the module and read `Dimension.exponents` in-process** —
not by diffing two transcripts. stage011 is the opposite: its `.py` prints `(2,1,-2)` (`py-stdout:62`) and
its `.wl` prints `{2, 1, -2}` (`out:78`) — same values, **different bracket and spacing conventions**, so
even there the "comparison" needs a parser per format.

**Consequence:** step (d) is not a diff. It requires, for every stage, (i) a `.py`-side value emitter or an
importable accessor, (ii) a `.wl`-side parser, and (iii) a **semantic name map** (the pilot needed
`zero`↔`1`, `mGNLS`↔`m_GNLS`). None of that generalises for free. **Budget it as real work per stage, and
do not assume 004 is a template for it — 004 is a template for (a)–(c) only.**

### 5.2 Step (a) is blocked or vacuous for 13 stages

The 13 zero-machinery stages have nothing to rewrite in (c) and nothing to compare in (d):

- **UNCOVERED (8) — the `.wl` computes nothing dimensional either.** 001 (no dim content on either side) ·
  014 (fully nondimensionalised; `.wl` "rank" is `MatrixRank`) · 015 (`wl:196`, `:201` `MatrixRank`) ·
  017 (`wl:34` binds `cited016DimensionalOk = True`, asserted `:560` — a literal) · 019 (units-free
  abstract algebra; ω powers are polynomial degrees) · 022 (`wl:588` asserts only `liveNames === {z}`) ·
  033 (`NullSpace`/`MatrixRank`, Dirac constraint counting) · 043 (both engines print
  `DIMENSIONAL_HOMOGENEITY=N/A_INTEGER_COUNT_STAGE`; all "dimension" there is Krull/CAD).
  ⛔ **Steps (a), (c) and (d) are all no-ops. Do not wire these in.**
- **PARTIAL (5) — the `.wl` restates dimension-*dependent* expressions but no declarations.** 020
  (`wl:93-96` `aPower` over the single symbol `a`; `:812`'s "dim-like name present" is a **substring check
  on identifiers**, satisfiable by renaming) · 024 (`wl:141` restates `gBase` by a different route) ·
  025 (`wl:119-121` restates `citedN0Den`) · 026 (`wl:92-94`, `:203-210` gates `aPower === -7/2`, **never
  printed to `.out`**) · 028 (`wl:99-102`, `:232-263` restate `Kbar_n` monomials).
  ⛔ Monomial-exponent drift would show; **no declaration error can**. ⚠ **020 is the sharpest**: `a`,
  `c_s`, `c`, `G` are all live and units-bearing, and the `.wl` power-counts in `a` alone.

### 5.3 ⭐ Step (a) — `.wl` value reachability, all 21 stages surveyed

Three independent read-only `.wl` surveys (bands 002–009, 010–036, 037–044). **`REACHABLE`** = an existing
named map or return value is live at a call site, so a print appended at end-of-file or at the call site
works — the stage004 pattern. **`LOCAL_ONLY`** = values die at a `Module` boundary; a print must go
**inside** the body. **`NOTHING_TO_PRINT`** = the file computes no exponent vector at all.

| stage | verdict | locus | # values | `.wl` axis order | note |
|---|---|---|---|---|---|
|005|**REACHABLE**|`wl:155` `dimensionDictionary[]`, map `:184-205`; bound as `dims` `:608`|21|declared `L,T,M` (`:4`)|exact stage004 pattern — drop-in|
|007|**REACHABLE**|`wl:290` `runDimensionalFirewall[]`, returned map `:650-662`; bound `:941`|11 exported|declared `L,T,M` (`:4`)|⚠ the map is a deliberate **subset**; ~49 intermediates stay `Module`-local|
|006|**REACHABLE (partial)**|13 plain top-level globals `wl:150-162`|13 + ~25 local|declared `L,T,M` (`:92-93`)|easiest insertion in the corpus for the 13; the ~25 derived row-dims are local to `runLegADimensions[]` `:221`|
|032|**REACHABLE**|`wl:836-842`|3 real (A, U, F)|implied `L,T,M`|⭐ **no `Module` anywhere** — driver is a bare `ok = Catch[` `:602`, so `aDimension`/`uDimension`/`fDimension` are globals surviving past `:955`. Cheapest of all|
|034|**REACHABLE**|`dimensionObject[]` `wl:324-362` returns the whole map; call sites `:670`, `:721` assign to globals|10|implied `L,T,M` (`:321`)|closest analogue to stage004. ⚠ `FieldIdentity` `:332-339` is a *declared* table; only `DensityDimension`/`ActionDimension` are genuinely computed (log-derivative, `:316-322`)|
|038|**REACHABLE**|`unitState[corrupt_]` `wl:675`; **returns all 8 as a bare `List` `:709-718`**; bound `:1092`|8 (+ literal twin `:721-730`)|4-axis `(M,L,T,E)`, decodable only from a hardcoded string `:1104`|⭐ **upgrade** — the prior "computes but cannot print" read was wrong. `muDim={2,1,-4,-1}` *is* live at top level; `out:31` omits `mu_R` only because `:1103-1105` prints a **hardcoded string** instead of the computed list|
|039|**REACHABLE**|constants `wl:408-410` top level; `unitState` `:415` returns `{curlUDimension, bDimension}` `:423`; call site `:776`|~5 named / 3 distinct|implied `(M,L,T)` (`:409`)|easiest of band C; same hardcoded-string issue in `out:34`|
|002|**LOCAL_ONLY**|`wl:46` `dimensionalResiduals[]` returns **scalar residuals**; the 20 vectors are locals `:54-73`|20 base + ~19 term|**not stated** — inferable `L,T,M` from `:54-56`|insert at `wl:74`. ⚠ `residualForUnits` `:42` is also called from the mutation probes `:355` — a print there **double-emits**|
|003|**LOCAL_ONLY**|`wl:606` `runDimensions[]`; vectors `:613-635`|23 + 18 rows|**not stated — and it is `(M,L,T)`**|insert at `wl:636`. ⭐ **LANDMINE**: `dRhoBr={1,-3,0}` `:621` cross-checks against stage006's asserted `[rho_br]=M L⁻³` `:210`. Emitting `axes=L,T,M` here silently corrupts every triple|
|008|**LOCAL_ONLY / near-empty**|`wl:475` `runDimensionalBlock[dtn_]`; vectors `:478-482`|**3**|**not stated**, only 2 axes|`axes=L,T,M\|exponents={1,0}` is malformed. Needs a distinct `axes=L,T` header, or skip|
|009|**LOCAL_ONLY**|`wl:428` `runDimensionalBlock[]`; vectors `:435-449`|~15|**not stated** — inferable `L,T,M`; **M slot is identically 0 in every vector**|no map, no consumer; insert at `wl:450`|
|010|**LOCAL_ONLY** (trivial insert)|`runDimensionalBlock[]` body `wl:492-511`; vectors `:496-505`|11 exact|**declared `M,L,T`** (`:511`)|the body *already* prints at `:511`, so one line at `:510` renders all 11. ⚠ slot 1 is 0 in every vector and `dimM0={0,0,-1}` = T⁻¹ — the declared order is **not evidenced by any vector**|
|030|**LOCAL_ONLY**|`wl:169-185` (basis + derived) and `:503-512` (`dimensionTerms`)|~29|declared `[L,T,M]` (`:499`)|locals of the terminal `Module[` `:140`, body ends in `Exit[]`. Two in-body prints (after `:185`, after `:512`) reach everything|
|031|**LOCAL_ONLY**|`unitRows` assigned `wl:447-469`, declared local `:145` of the terminal `Module[` `:121`|26|declared `[L,T,M]` (`:444`)|consumed by the `Do` `:470-476`, discarded. One print at `:470` covers all 21 rows|
|040|**LOCAL_ONLY**|`dimensionState[mutateCE_]` `wl:807`; 11-entry map `:816-828`; returns **aggregates only** `:854-859`|~20|**not stated** — inferable `(M,L,T)` (`:817`, `:821`)|only the 7 lock/stability *combinations* reach `dimState` `:1621`. ⚠ `regionDimensionFor` `:531` / `dimensionRecord` `:546` are `RegionDimension`, **not** units|
|041|**LOCAL_ONLY (partial)**|`dimensionState[corrupt_]` `wl:1197`; bases `:1207-1211`; **returns 3 of 7** `:1214`|7 (+2 top-level literals)|**not stated** — inferable `(M,L,T)` (`:1207`, `:1211`)|⚠ **the mutation-carrying value `dimSpeed` `:1211` is precisely one of the unreachable four**|
|042|**LOCAL_ONLY**|`dimensionGuard[mutate_]` `wl:819`; bases `:827-832`; 14-entry map `:833-848`; returns aggregates `:898-907`|~34|⛔ comment `:816` says **"MLT" — a MISLABEL**; the real basis is `(stiffness,L,T)` (`stiffnessDim={1,0,0}` `:828`, `chargeDim={1/2,3/2,-1}` `:832`)|~13 *derived* vectors do escape (`scalar_power`, `em_power`, `ratio`, the `*_terms` lists); the **14 base per-quantity vectors** are lost|
|044|**LOCAL_ONLY**|`dimensionalTooth[]` `wl:482`; 27-entry map `:493-506`; 24 `integrandDims` `:526+`; **returns `Null` `:587`**|**~95** — largest in the corpus|⭐ **explicitly declared** `{L, T, M}` (`:576`), corroborated `:1104`, `:1183`|⛔ **exactly 1 of ~96 vectors is protocol-reachable** — the `Z_chi` master row `"[L,T,M]" -> {-2,0,1}` `:1104`, which reaches the global `evidence` via `:1192`. Everything else is discarded|
|035|**REACHABLE, but no exponent vectors exist**|`unitScalingObject[]` `wl:411-441` returns the whole map; global `restored` `:792`|9 keys / 20 entries|**none** — values are **scale monomials** in `lengthScale/timeScale/massScale`|⚠ emitting `exponents={a,b,c}` needs an inline `Exponent[#, …]` extraction — a **rendering conversion**, not a dump|
|036|**NOTHING_TO_PRINT**|`unitRules[]` `wl:232-245`; `dimensionResiduals[]` `:252-296` returns **residuals only** (all zero)|12 inputs carry a dim; ~29 checks use inline literals|**none** — monomials only|⛔ derived-quantity dimensions (`I_ij`, `D_V`, `A_V`, `U_A`, `F_A2`) exist **only** as literal argument monomials `:253-296` and hardcoded prose `:660-662`. Printing them means **re-transcribing constants, not rendering a computed value**|
|037|**NOTHING_TO_PRINT**|`unitContract` `wl:643-699` (top-level assoc)|21 named, **0 exponent vectors**|**none** — `massScale/lengthScale/timeScale` monomials (`:585-634`)|⛔ the file contains **zero integer-tuple objects**. `unitViolations` `:701` returns only `Keys`. A `DIM|` line requires an `Exponent[#, {massScale,lengthScale,timeScale}]&` extraction the file never performs — **that is a derivation, not a print**|

**Two cross-cutting corrections to the working assumption:**

1. ⭐ **Only 030, 031, 040, 041, 042, 044, and the 002/003/008/009/010 band actually have the
   "terminal `Module` swallows the map" problem.** 032, 034, 035, 036 contain **no top-level `Module`** —
   their drivers are bare `ok = Catch[…]` (`:602`, `:546`, `:641`, `:510`), so every assignment inside is a
   global that outlives the run and is printable from a line appended at end-of-file.
2. ⛔ **In-body prints are NOT safe for 040, 041, 042.** Mathematica `Module` bodies are
   `CompoundExpression` sequences, so inserting `Print[…];` between statements is semantically inert — but
   these guards are **invoked more than once, with mutation flags set**: `042:1845` passes
   `activeMutation === "DIMENSION_HOMOGENEITY"`, `040:1621` likewise, and `041` calls `dimensionState`
   twice (`:2076` mutated, `:2101` clean). An in-body print would emit **duplicate and deliberately
   corrupted `DIM|` lines** into the `.out` — worse than emitting none. Getting one clean top-level
   emission for these three requires the restructuring that is out of scope.

⛔ **Policy question the next session must settle before step (a):** does "print-only" permit in-`Module`
insertion? If yes, 002/003/008/009/010/030/031/044 become mechanical. If no, only 005, 006 (top tier), 007,
032, 034, 038, 039 qualify — and 035/036/037 never do.

### 5.4 Step (b) needs a Mathematica seat for ~21 stages

Only 9 `.out` files render a computed dimension value today (004, 011, 012, 013, 016, 018, 021, 023, 027).
The other 21 COVERED stages reduce the whole dimensional firewall to `PASS <TOOTH>`, so each needs a `.wl`
edit **and** a re-run. ⚠ The licence cap is **2 seats**, `mathematica/run_all_audits.sh` sleeps 10 s between
scripts, and the `.wl` runner is far slower than the Python suite (~3–3.5 min sequential for all 43 `.py`,
dominated by stage033 ≈ 82 s). **Batch the `.wl` work; do not interleave it one stage at a time.**

### 5.5 Step (d) is structurally incapable of failing for some stages

Even with both sides emitting values, these comparisons cannot detect a module bug:

- **031** — 20 of 21 `PASS_UNITS_*` rows compare an expression to a copy of itself (`py:553-580`,
  `wl:447-470`). Under the module both copies read one shared constant, so the row cannot fail for *any*
  value. **The rewrite of 031 has no working check.**
- **017** — `assert True` in two languages. Nothing to compare.
- **034** — `wl:327-339` vs `:673-682` is the same declaration-vs-itself shape.
- **044** — the `Z_chi` master row (`py:1038-1039` / `wl:1102-1110`), same shape.
- **040** — `[rho]` enters `dims["K"]` at −4 and `lock_A` at +4, so `M L⁻³ → M L⁻⁴` is invisible **in both
  files**; same for `[m]`, `[mu_R]`.
- **042** — `[omega]` and `[d]` cancel out of `ratio_dim`.
- **039** — asserts only `[curl u_T] == [b_T]`, so scaling `LENGTH` and `INVERSE_LENGTH` together passes.
- **038** — `unit_state()` (`py:708`) vs `EXPECTED_UNIT_STATE` (`py:741`) are both in the same file and move
  together; the `.wl` never prints `mu_R` to contradict them.
- **⚠ Sum-of-squares residuals are permutation-invariant.** `dim_residual` in
  `scripts/ledger_dimensions.py:169-179` pairs by axis label and is safe, but every *stage-local* residual
  that zips positionally is not. A transposition is symmetric by construction.

### 5.6 The five stages where step (a) is genuinely impossible, ranked

Everything else is a cost question. These five are the ones to plan around:

| stage | why | consequence |
|---|---|---|
|**037**|`.wl` contains **zero integer-tuple objects**. Its firewall is purely symbolic monomial substitution.|A `DIM|` line would be a **new derivation added to the `.wl`** — not print-only. Either accept that (and re-review the `.wl`), or accept that 037 has no cross-engine dimension check at all|
|**036**|Derived dimensions exist only as **literal arguments** and hardcoded prose.|Printing them re-transcribes constants. A transcription bug would be invisible **by construction**|
|**035**|Values are scale monomials, not vectors.|Needs an `Exponent[…]` conversion inside the print — a small derivation, reviewable, but not a dump|
|**044**|**1 of ~96** vectors is reachable; `dimensionalTooth[]` returns `Null`.|The largest computed dimension set in the corpus is 99 % unobservable. Compounded by 044 being frozen pending **044-v2**|
|**042**|14 base vectors local **and** the guard runs twice with a mutation flag.|In-body print emits corrupted duplicates. No safe print-only route exists|

---

## 6. Commands run

```
git status --short ; git log --oneline ; git diff --stat
git show HEAD:<path> ; git show d9544a62:<path> ; git show 89bcccdb:<path>   # pre-rewrite 004/011 idioms
python3 _scratch/pass1_dim_survey/measure/classify_dim_lines.py              # unmodified, re-run
ls scripts/ledger_stage*.py | grep -vE 'stage004|stage011|stage044' \
  | xargs -P 6 -I{} python3 {}                                              # 40 scripts, all exit 0
sha256sum mathematica/out/ledger_stage*.out ; wc -l mathematica/out/ledger_stage*.out
grep -c '^PASS' <stdout> ; grep -oE '^PASS tally: [0-9]+' <stdout>          # per stage
grep -ln '^ACTIVE_MUTATION' scripts/ledger_stage*.py
grep -n 'lru_cache|open\(|write_text|mkdir|Path\(__file__\)|inspect\.|hashlib\.' scripts/ledger_stage*.py
grep -n '^import |^from ' scripts/ledger_stage*.py | grep -iE 'ledger_stage|ledger_dimensions'
diff <(strip-header scripts/output/*.txt) <fresh stdout>                    # 26/28 byte-identical bodies
```

Plus **three parallel read-only `.wl` surveys** (bands 002–009, 010–036, 037–044) for §5.3, each answering
"could a print-only statement render every per-quantity dimension vector this file computes, without
restructuring?" — no Mathematica was executed by any of them.

**Deviations, stated for the record.** Stages **004**, **011** and **044** were not executed in place.
Instead `git show HEAD:<path>` copies were written into the session scratchpad and run there — the corpus
was never touched. The stage044 copy was nested so that its `Path(__file__).resolve().parents[3]` write
landed inside the scratchpad, not the project. 004/011 needed
`PYTHONPATH=research/pde_ledger_v2/scripts` to import `ledger_dimensions`; that is a read-only import.
**No `math -script` was run.** Working tree verified clean at `adcfbdfd` after all measurement.
