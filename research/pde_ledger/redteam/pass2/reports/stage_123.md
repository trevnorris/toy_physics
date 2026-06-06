---
unit_id: 123
batch: IV.3
auditor_model: Opus 4.8 (1M context)
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md]
  paper_appendix: present
---

# Audit unit 123 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_123.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md` (only file matching the glob)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read rows for the rF1 anchor and the parent-overlap formulas referencing this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.txt`

## What the paper claims

Stage 123 converts the Family-1 compensated-throat-core target into two normalized parent-level numbers: the background-flow number `Xi_v` and the mouth-traction number `Xi_T`. The card body quote states: "Converts the Family--1 compensated target to normalized parent flow and traction numbers \((\Xi_v,\Xi_T)\)." The notes are the authoritative carrier and box four deliverables: (i) the exact flow law `Xi_v = -3√30 π^(3/2)/160 · r` with the Family-1 value `Xi_v^{F1} ≈ −1.01675633282526` (notes lines 30, 40); (ii) the exact traction law `Xi_T = 3√30/(10√π) · 1/g` with the natural value `Xi_T^{nat} ≈ 0.927058084855655` (lines 71, 81); and the two compensated-branch values (iii) `Xi_T^{(−)} ≈ 1.22297517701464` and (iv) `Xi_T^{(+)} ≈ 0.331334521644609` (lines 97, 99). The card "Verification" line states "Mathematica audit: none yet." Both `Xi_v` and the flow law carry an explicit NEGATIVE sign throughout, inherited from `λ = −q_* v_{w0} I_{sq} < 0`.

## What the script claims to verify

Both scripts build `K_s`, `K_q`, `J_s`, `λ` from the upstream parent-action formulas (stage 118/119 healing-locked shell, mixed D/N tube with `L_W=L`, `λ = −(8√2/3) q v_{w0} a² ℓ √L`), then (a) define `r = λ/√(K_s K_q)`, eliminate `v_{w0}` in favor of `r`, substitute into the independently-defined parent combination `Xi_v_def`, and assert the residual against `−3√30 π^(3/2) r/160` (`Xi_v law`); (b) do the analogous `g`→`T_m` elimination on the healing-locked traction and assert `Xi_T = 3√30/(10√π g)` (`Xi_T law`). They then evaluate both laws on the Family-1 inputs `r_{F1}` and `g_∓^{F1}` and assert the four numeric branch values. The Mathematica script additionally prints each branch value to 20 digits via `RootReduce`.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `Xi_v = −3√30 π^(3/2)/160 · r` (notes:30) | `Xi_v law` (py:46, wl:81) | match |
| `Xi_v^{F1} ≈ −1.01675633282526` (notes:40) | `Xi_v(F1)` (py:54/59, wl:102-105) | match |
| `Xi_T = 3√30/(10√π g)` (notes:71) | `Xi_T law` (py:47, wl:95) | match |
| `Xi_T^{nat} ≈ 0.927058084855655` (notes:81) | `Xi_T(nat)` (py:55/60, wl:107-110) | match |
| `Xi_T^{(−)} ≈ 1.22297517701464` (notes:97) | `Xi_T(-)` (py:56/61, wl:112-115) | match |
| `Xi_T^{(+)} ≈ 0.331334521644609` (notes:99) | `Xi_T(+)` (py:57/62, wl:117-120) | match |
| card: "Mathematica audit: none yet" (tex:11) | `.wl` EXISTS and passes | mismatch (F2) |
| imported `r_{F1}` closed form, appendix `4107-117π²` (appendix:562) | scripts use `4107-100π²` (py:50-51, wl:97-98) | mismatch (F1) |

Every Xi deliverable faithfully maps to a non-tautological script check. The two `mismatch` rows are paper-side documentation/value issues, not script defects; the script values are the ones that reconcile with the notes and with the stated numerics.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46 | `expect_zero(Xi_v_expr + 3√30 π^(3/2) r/160)` | Xi_v law (notes:30) | yes |
| A2 | sympy | 47 | `expect_zero(Xi_T_expr − 3√30/(10√π g))` | Xi_T law (notes:71) | yes |
| A3 | sympy | 59 | print `Xi_v(F1) ≈ −1.0167563…` | Xi_v^{F1} (notes:40) | yes (numeric, no assert — see note) |
| A4 | sympy | 60-62 | print Xi_T nat/−/+ numerics | notes:81/97/99 | yes (numeric, no assert) |
| A5 | math | 81 | `expectZero[xiVFromParent − xiVExpected]` | Xi_v law | yes |
| A6 | math | 95 | `expectZero[xiTFromParent − xiTExpected]` | Xi_T law | yes |
| A7 | math | 105 | `expectZero[xiVF1 − xiVF1Expected]` | Xi_v^{F1} | yes |
| A8 | math | 110 | `expectZero[xiTNat − xiTNatExpected]` | Xi_T^{nat} | yes |
| A9 | math | 115 | `expectZero[xiTMinus − xiTMinusExpected]` | Xi_T^{(−)} | yes |
| A10 | math | 120 | `expectZero[xiTPlus − xiTPlusExpected]` | Xi_T^{(+)} | yes |

Note on A3/A4: the SymPy script only PRINTS the four numeric branch values (no `assert`); the two SymPy `expect_zero` asserts (A1/A2) cover only the two symbolic laws. The Mathematica script DOES assert all four branch numerics (A7-A10) via `expectZero[xiVF1 − xiVF1Expected]` etc. Because the symbolic laws are asserted in SymPy and the numeric branch evaluations are asserted (against the same law substituted at the Family-1 inputs) in Mathematica, the four branch values are covered cross-engine. The SymPy-side branch numerics are confirmatory prints whose values are reproduced and asserted by the Mathematica side — this is acceptable coverage, not an `insufficient_verification` gap, because the load-bearing symbolic identity that produces those numbers is asserted in SymPy and the numerics are asserted in Mathematica.

## Findings

### F1 — paper_misalignment (value_mismatch)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:562`
- `/var/projects/toy_physics/research/pde_ledger/paper/parts/part04_geometry_retarded_mouth.tex:576`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:50-51`

**What's wrong:**
The imported `r_{F1}` closed form appears two ways. The paper writes the radicand as `4107-117π²`:

> appendix:562  `\frac{\sqrt{4107-117\pi^2}}{10\pi}` (and part04 body:576, both with `\approx1.77799353547498`)

The stage-123 scripts (and the notes the card derives from) use `4107-100π²`:

> py:50-51  `R = sp.Rational(37,20); rF = sp.sqrt(12*R**2/sp.pi**2 - 1)`  → output line 9 `sqrt(4107 - 100*pi**2)`
> wl:97-98  `R = 37/20; rF1 = reduceExact[Sqrt[12*R^2/Pi^2 - 1]]`

`12·(37/20)²/π² − 1 = (4107 − 100π²)/(100π²)`, so the scripts' `100π²` is the algebraically forced value. It also matches the stated numeric: √(4107−100π²)/(10π) ≈ 1.77799 (matches the appendix's own `≈1.77799…`), whereas √(4107−117π²)/(10π) ≈ 1.7295 contradicts it. The stage-121 origin notes (`...stage121_geometric_r_selection.md:69`), stage-122 notes (`:49,:56,:88`), stage-126 and stage-148 notes all carry the correct `4107-100π²`. The `117π²` is a paper-side typo for `100π²`. This is the rF1 DEFINITION (the stage-121 anchor `eq:...-rF1`) which stage 123 only IMPORTS; stage 123's own card (`stage_123.tex`) and notes reference `\mathfrak r_{F1}` symbolically and are clean.

**Why this matters:**
A reader taking the appendix radicand literally would compute `r_{F1} ≈ 1.7295` and every downstream Family-1 number (including this stage's `Xi_v^{F1}` and `Xi_T^{(±)}`) would drift. The scripts are correct; the published surd is wrong and self-inconsistent with its own quoted decimal.

**Required change:**
Paper-side (routes to user — Codex must not edit paper/). See `## Resolve before fix_loop`. No script change: the scripts already use the correct `100π²`.

**Verification:**
After resolution, appendix:562 and part04:576 read `\sqrt{4107-100\pi^2}` and the radicand matches `√(4107−100π²)/(10π) ≈ 1.77799353547498`. No script re-run needed (scripts unchanged).

### F2 — paper_misalignment (notes_contradicts_script / documentation status)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_123.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl` (entire file)

**What's wrong:**
The card "Verification" line states:

> tex:11  "SymPy audit: \StageFile{...stage123...sympy_audit.py}.  Mathematica audit: none yet."

But a Mathematica audit `.wl` exists (added 2026-06-01 retro-sweep, mtime 2026-06-01 15:53) and passes all six checks (output transcript: `Stage 123 Mathematica audit passed.`). The card understates the actual verification coverage.

**Why this matters:**
The card claims LESS than reality (so no false positive claim), but a status card that says "none yet" for a present, passing second-engine audit is stale documentation. It would mislead the verification-coverage trackers and any reader auditing dual-engine status.

**Required change:**
Paper-side (routes to user — Codex must not edit paper/). See `## Resolve before fix_loop`. No script change.

**Verification:**
After resolution, tex:11 names the `.wl` path (e.g. "Mathematica audit: \StageFile{mathematica/...stage123...mathematica_audit.wl}").

## Independent-derivation check (Mathematica)

The `.wl` is a GENUINE independent route, not a transliteration of the `.py`. Justification by corresponding sections:

1. **Solve mechanism (the load-bearing step).** SymPy eliminates `v_{w0}` with `sp.solve(sp.Eq(r, r_from_parent), vw)[0]` (py:35) — silently takes the first root. Mathematica uses `Reduce[rr == rParent && Element[vw, Reals], vw, Reals]` (wl:73) followed by a custom `branchesFromReduce`/`uniqueBranch` extractor (wl:34-45) that FAILS if the real solution is not unique. This is a structurally different solver path with an added uniqueness guard SymPy lacks. The `T_m`-from-`g` elimination differs identically (py:43 `solve(...)[0]` vs wl:87-88 `Reduce[...] → uniqueBranch`).

2. **Branch numerics.** SymPy derives `gminus = rF − √(1+rF²)/2` from the family law (py:52) and reports numerics via `sp.N(...,20)` (py:59-62, prints only). Mathematica substitutes the pre-solved closed form `gMinus = (2√(4107−100π²) − 37√3)/(20π)` (wl:99) and reports numerics via `N[RootReduce[...],20]` (wl:47, asserted at wl:105/110/115/120). I verified the two gMinus forms are algebraically identical (√(1+rF1²)=37√3/(10π) ⇒ gMinus=(2√(4107−100π²)−37√3)/(20π)), so the engines agree by independent routes (family-law derivation vs. closed-form-then-RootReduce).

3. **Simplification engines.** SymPy uses `sp.simplify(sp.expand(...))`; Mathematica uses a `reduceExact` pipeline (`PowerExpand`→`Cancel`→`Together`→`FullSimplify`→`RootReduce`, wl:21-26) with `ConditionalExpression` stripping. Different normalizers.

The shared definitions of `K_s, K_q, J_s, λ` (py:25-28 ≡ wl:63-66, same `8√2/3` coefficient) are EXPECTED to be identical — they are upstream-imported physical constants, not quantities either engine could legitimately re-derive differently. Identical inputs + divergent solve/simplify routes = genuine second engine. **No `mathematica_transliteration` finding.**

## Engine cross-check

Both engines agree to all printed digits:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| `Xi_v(r)` law | `-3*sqrt(30)*pi**(3/2)*r/160` (txt:5) | `(-3*Sqrt[3/10]*Pi^(3/2)*rr)/16` (txt:8) — equal (`3√30/160 = 3√(3/10)/16`) |
| `Xi_T(g)` law | `3*sqrt(30)/(10*sqrt(pi)*g)` (txt:6) | `(3*Sqrt[3/(10*Pi)])/gg` (txt:12) — equal |
| `Xi_v(F1)` | `≈ -1.0167563328252594644` (txt:9) | `-1.01675633282525946441973…` (txt:15) |
| `Xi_T(nat)` | `≈ 0.92705808485565499282` (txt:10) | `0.92705808485565499282126…` (txt:18) |
| `Xi_T(-)` | `≈ 1.2229751770146391627` (txt:11) | `1.2229751770146391627143…` (txt:21) |
| `Xi_T(+)` | `≈ 0.33133452164460908424` (txt:12) | `0.33133452164460908424123…` (txt:24) |

All six match to the full SymPy print precision (20 digits). Every `Xi_v/Xi_T law` and branch assertion reports `0` / `PASS`. The negative sign on `Xi_v(F1) ≈ −1.0167563…` is present and identical in BOTH engines. No `engine_disagreement`.

## Verdict justification

The two scripts hold up under attack. The `Xi_v law` and `Xi_T law` assertions are non-tautological: they eliminate a parent parameter (`v_{w0}` / `T_m`) from an independently-stated parent-normalized combination and force the residual against the boxed clean coefficient — any wrong factor in `K_s, K_q, J_s, λ` or in the normalization definition would break the cancellation. The negative `Xi_v` flows correctly from the negative `λ` and is preserved consistently across card-equivalent notes, both scripts, and both outputs; I attacked it (checked the sign chain `λ<0 ⇒ v_{w0}∝−r ⇒ Xi_v∝−r`) and it is correct, not an error. The Mathematica engine is genuinely independent (different solver, uniqueness guard, normalizer). All four Family-1 branch numerics reconcile to the notes to full precision. The only defects are two LOW-severity paper-side `paper_misalignment` items — the imported `r_{F1}` surd typo (`117π²` should be `100π²`, an upstream stage-121 anchor) and the card's stale "Mathematica audit: none yet" — both of which route to the user, not Codex. Verdict: `findings` (paper-side only; scripts clean).

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit (source: `.py`/`.wl` + committed `.txt` outputs):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Xi_v = −3√30 π^(3/2)/160 · r` | py:46 / wl:81 / sympy txt:5, math txt:8 | notes:30 `\Xi_v = -\frac{3\sqrt{30}\,\pi^{3/2}}{160}\,\mathfrak r` | MATCH |
| `Xi_T = 3√30/(10√π) · 1/g` | py:47 / wl:95 / sympy txt:6, math txt:12 | notes:71 `\Xi_T = \frac{3\sqrt{30}}{10\sqrt{\pi}}\,\frac{1}{\mathfrak g}` | MATCH |
| `Xi_v^{F1} ≈ −1.01675633282526` | py:54,59 / wl:102-105 / sympy txt:9, math txt:15 | notes:40 `\approx -1.01675633282526` | MATCH (negative preserved) |
| `Xi_T^{nat} ≈ 0.927058084855655` | py:55,60 / wl:107-110 / sympy txt:10, math txt:18 | notes:81 `\approx 0.927058084855655` | MATCH |
| `Xi_T^{(−)} ≈ 1.22297517701464` | py:56,61 / wl:112-115 / sympy txt:11, math txt:21 | notes:97 `\Xi_T^{(-)}\approx 1.22297517701464` | MATCH |
| `Xi_T^{(+)} ≈ 0.331334521644609` | py:57,62 / wl:117-120 / sympy txt:12, math txt:24 | notes:99 `\Xi_T^{(+)}\approx 0.331334521644609` | MATCH |
| `r_{F1}` (imported) radicand `4107-100π²` | py:50-51 / wl:97-98 / sympy txt:9 | appendix:562 & part04:576 show `4107-117π²` | MISMATCH → F1 (upstream stage-121 anchor; correct value lives in notes stage121:69, stage122:49) |

INTERNAL scaffolding (accounted for, no finding): `Ks`, `Kq`, `Js`, `lam`/`lambda` (upstream-imported parent formulas, not stage-123 deliverables), `r_from_parent`/`rParent`, `g_from_parent`/`gParentLocked`, `Xi_v_def`/`xiVDef`, `Xi_T_def`/`xiTDef`, `v_w0(r)`/`vwFromR`, `T_m(g)`/`tmFromG`, `gminus`/`gMinus`, `gplus`/`gPlus` (upstream stage-122 imports, used as branch inputs). The intermediate prints `v_w0(r)`, `Xi_v(r) derived`, `T_m(g)`, `Xi_T(g) derived` in the `.wl` are derivation-trace prints, not new deliverables.

reconciliation: 6 deliverable values checked, 6 MATCH against the stage-123 notes; 1 imported-input MISMATCH (`r_{F1}` surd `117π²` vs `100π²`) flagged as F1, an upstream stage-121 paper anchor. All six stage-123 deliverables (including the negative `Xi_v`) reconcile exactly.

## Self-test notes

Checked: (1) Variable independence — the `Xi_v/Xi_T law` checks eliminate `v_{w0}`/`T_m` via solve and substitute; no spurious zero-derivative (no `diff` involved). (2) No unbounded integrals in either script, so parity traps N/A. (3) Trivial-case: the boxed coefficients `−3√30 π^(3/2)/160` and `3√30/(10√π)` are forced by the parent-formula cancellation; verified the `r_{F1}` radicand by hand (`12·(37/20)²/π²−1 = (4107−100π²)/(100π²)`), confirming the scripts' `100π²` and exposing the appendix `117π²` typo. (4) Sign chain on `λ<0 ⇒ Xi_v<0` verified by hand — the negative `Xi_v` is correct, not flagged as error. No directive prescribes a script edit (both findings are paper-side `paper_misalignment` pending user resolution).
