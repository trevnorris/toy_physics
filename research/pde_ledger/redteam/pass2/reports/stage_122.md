---
unit_id: 122
batch: IV.3
auditor_model: Claude Opus 4.8 (1M context)
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md]
  paper_appendix: present
---

# Audit unit 122 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_122.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: `\input{stages/stage_122}` at L1278; the MTDC-T8 derivation prose at L538–573 carries the stage's r_F1 / g_± results)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.txt`

## What the paper claims

The stage tests whether the simplest concrete "equal-normalized" local mouth source (no channel favoritism, giving normalized mouth-coupling ratio `𝔤_nat = 1`) already lands on the Stage-119 compensation family `1+𝔯² = 4(𝔤−𝔯)²`, equivalently `𝔤_comp^± = 𝔯 ± ½√(1+𝔯²)`. Inserting the Stage-121 geometric value `𝔯_F1 = √(4107−100π²)/(10π) ≈ 1.778` gives the two exact compensated couplings `𝔤_-^F1 ≈ 0.758035078944663` and `𝔤_+^F1 ≈ 2.79795199200529`. The verdict (card `quote`): "Equal-normalized mouth source misses the lower compensated branch but only by a modest traction renormalization." Distinct deliverables: (D1) `𝔤_±^F1` closed forms and numerics; (D2) the compensation defect `𝒞_nat = 1+𝔯_F1²−4(1−𝔯_F1)² = (−12321+80π√(4107−100π²))/(100π²) ≈ 1.74017` is nonzero (natural branch off-family); (D3) the nearest-branch gaps `Δ𝔤_- ≈ 0.24196`, `Δ𝔤_+ ≈ 1.79795`; (D4) the traction renormalizations `𝒯_m^(±)/𝒯_m^nat = 1/𝔤_±^F1 ≈ 1.31920` and `0.357404`, derived from the Stage-119 law `𝔤 ∝ 1/𝒯_m`. The card body is terse; the notes file is the authoritative carrier and enumerates all four deliverables with boxed values.

## What the script claims to verify

Both scripts: (i) build `𝔯_F1` from the geometric radius `L/a = 37/20`; (ii) obtain the two compensated branches `g_-`, `g_+` (sympy via the quadratic-formula closed form `𝔯∓½√(1+𝔯²)`; mathematica via `Solve` of the compensation quadratic with sign-based branch selection); (iii) certify each against the boxed exact surd forms `(2√(4107−100π²)∓37√3)/(20π)`; (iv) confirm both branches satisfy the compensation quadratic at `𝔯_F1`; (v) compute the natural-branch defect, certify its closed form, and assert it is nonzero; (vi) derive the traction ratios `𝒯_m^(±)/𝒯_m^nat` with the constant `C` carried as a free positive symbol, then assert each equals `1/g_±` only after substituting the natural-branch ansatz `g_nat=1` last. Print-only emissions cover the numerics for the defect, the Δ𝔤 gaps, and the traction ratios.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1: `𝔤_±^F1` closed forms `(2√(4107−100π²)∓37√3)/(20π)` & numerics | sympy L62–63 `expect_zero gminus/gplus exact form`; wl L85–86; numerics printed | match |
| D1': both branches solve compensation quadratic at `𝔯_F1` | sympy L65–66; wl L88–89 | match |
| D2: defect `𝒞_nat` closed form & nonzero (off-family) | sympy L69 `defect closed form`, L71 `natural off compensation`; wl L98–99 | match |
| D3: gaps `Δ𝔤_- ≈ 0.24196`, `Δ𝔤_+ ≈ 1.79795` | sympy L40–41 (print-only); wl: not printed | partial (print-only, no assertion — acceptable; these are derived numerics of already-asserted g_±) |
| D4: traction ratios `1/𝔤_±^F1 ≈ 1.31920 / 0.357404` | sympy L73–74; wl L121–128 (de-tautologized via free C + last-substituted g_nat) | match |
| `𝔯_F1 = √(4107−100π²)/(10π)` geometric form | sympy L24–25; wl L52–59 (`r_F1 geometric relation` asserted) | match (scripts/notes); **mismatch vs appendix L562 which states `√(4107−117π²)`** |

`paper_alignment: partial` — the card body and notes are fully aligned with the scripts; the single misalignment is the appendix prose carrying a wrong symbolic `𝔯_F1` (117 vs 100). Two additional paper-side staleness items (card calls this "Stage~139"; card says "Mathematica audit: none yet" though the `.wl` now exists) are noted under Verdict but are known numbering/staleness prose, not value defects.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `expect_zero(gminus − exact)` | D1 | yes |
| A2 | sympy | 63 | `expect_zero(gplus − exact)` | D1 | yes |
| A3 | sympy | 65 | `expect_zero(1+rF²−4(gminus−rF)²)` | D1' | yes |
| A4 | sympy | 66 | `expect_zero(1+rF²−4(gplus−rF)²)` | D1' | yes |
| A5 | sympy | 69 | `expect_zero(comp_def − defect_exact)` | D2 | yes |
| A6 | sympy | 71 | `expect_nonzero(comp_def)` | D2 | yes |
| A7 | sympy | 73 | `expect_zero(T_ratio_minus − 1/gminus)` | D4 | yes (g_nat=1 baked at L34, C free; see F-note) |
| A8 | sympy | 74 | `expect_zero(T_ratio_plus − 1/gplus)` | D4 | yes |
| B1 | wl | 59 | `expectZero(rF1² − (12·rad²/π²−1))` | r_F1 geom | yes |
| B2 | wl | 62–63 | `Solve` quadratic, assert exactly 2 roots | D1' | yes |
| B3 | wl | 85 | `expectZero(gMinus − exact)` | D1 | yes |
| B4 | wl | 86 | `expectZero(gPlus − exact)` | D1 | yes |
| B5 | wl | 88 | `expectZero(compResidual[rF1,gMinus])` | D1' | yes |
| B6 | wl | 89 | `expectZero(compResidual[rF1,gPlus])` | D1' | yes |
| B7 | wl | 98 | `expectZero(naturalDefect − defectExact)` | D2 | yes |
| B8 | wl | 99 | `expectNonzero(naturalDefect)` | D2 | yes |
| B9 | wl | 121–124 | `expectZero(residual_minus /. gNat→1)` | D4 | yes (residual = (20(gNat−1)π)/… stays symbolic before substitution) |
| B10 | wl | 125–128 | `expectZero(residual_plus /. gNat→1)` | D4 | yes |

## Findings

### F1 — paper_misalignment (value_mismatch)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:560-563`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:24-25` (and `notes/...stage122...md:49`)

**What's wrong:**
The part-IV appendix states the Family-1 branch value as

> `\mathfrak r_{F1} = \frac{\sqrt{4107-117\pi^2}}{10\pi} \approx 1.77799353547498`  (appendix L560–563)

but the scripts and the stage notes use `4107 − 100π²`:

> sympy L24–25: `R = sp.Rational(37,20); rF = sp.sqrt(12*R**2/sp.pi**2 - 1)` → `√(4107−100π²)/(10π)` (confirmed in sympy output L5–7, wl output L5).
> notes L49: `\mathfrak r_{F1}=\frac{\sqrt{4107-100\pi^2}}{10\pi}`.

The script form is the correct one: `12·(37/20)² = 41.07`, so `r_F1² = 41.07/π² − 1 = (4107 − 100π²)/(100π²)`. The appendix's `117π²` is arithmetically wrong as a symbolic form. Notably the appendix's own quoted numeric `≈1.77799353547498` corresponds to `4107 − 100π²` (gives 1.77799…), NOT to `4107 − 117π²` (which gives ≈1.7294…), so the appendix is internally inconsistent between its symbol and its numeric. The downstream `g_-^F1 ≈ 0.758035…` and `g_+^F1 ≈ 2.79795…` quoted at appendix L571/L573 also match the `100π²` form, confirming `117` is a stray typo, not an alternate convention.

**Why this matters:**
A reader pulling the symbolic `r_F1` from the appendix (the part-level result carrier for anchor MTDC-T8) would get a different surd than the one the verified scripts certify, and one that does not even reproduce the appendix's own numerics. This is a paper-side correctness defect on a load-bearing geometric input quoted forward to Stages 123/125–139.

**Required change:**
Paper-side; route to user. Do NOT have Codex edit the script — the script is correct. See `## Resolve before fix_loop` in the directive. Most likely direction: fix appendix L562 `4107-117\pi^2` → `4107-100\pi^2` to match scripts, notes, and the appendix's own quoted numerics.

**Verification:**
After resolution, appendix L562 reads `\frac{\sqrt{4107-100\pi^2}}{10\pi}`, matching script output L5 (`sqrt(-1 + 4107/(100*pi**2))` ⇔ `√(4107−100π²)/(10π)`) and notes L49. No script change; scripts already exit 0.

## Independent-derivation check (Mathematica)

The retro-sweep `.wl` (commit 251639c, 2026-06-01) is a **genuinely independent** route, not a transliteration of the `.py`. The decisive difference is how the central quantities `g_±` are obtained:

- sympy `.py` L26–27 writes the branches directly from the quadratic-formula closed form:
  ```python
  gminus = sp.simplify(rF - sp.sqrt(1+rF**2)/2)
  gplus  = sp.simplify(rF + sp.sqrt(1+rF**2)/2)
  ```
- mathematica `.wl` L61–73 instead **solves** the compensation quadratic and selects branches by sign:
  ```wolfram
  compensationEquation = compensationResidual[rF1, g] == 0;
  rootRules = Solve[compensationEquation, g, Reals];
  If[Length[rootRules] =!= 2, fail[...]];
  roots = reduceExact /@ (g /. rootRules);
  gMinus = SelectFirst[roots, TrueQ[FullSimplify[# - rF1 < 0, ...]] &];
  gPlus  = SelectFirst[roots, TrueQ[FullSimplify[# - rF1 > 0, ...]] &];
  ```
  The `.wl` does not hard-type the `𝔯∓½√(1+𝔯²)` form; it lets `Solve` produce the roots and disambiguates the ± by the sign of `g−r_F1`. The closed surd `(2√(4107−100π²)∓37√3)/(20π)` is then a separately-derived *certificate* it checks the Solve output against (L85–86), the mirror image of where the `.py` starts.

- The traction block is also independently structured: the `.wl` keeps `gNat` a free symbol through L101–113 and prints the *symbolic* residual `(20(gNat−1)π)/(…)` (output L29–30) before substituting `gNat→1` at L115/L123, whereas the `.py` bakes `g_nat=1` at L34 and works numerically. Same conclusion, different choreography.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every shared quantity:

| Quantity | sympy output | mathematica output |
|---|---|---|
| `r_F1` (symbolic) | `sqrt(-1 + 4107/(100*pi**2))` (L5) | `Sqrt[4107 - 100*Pi^2]/(10*Pi)` (L5) — same value |
| `g_minus` | `(-37*sqrt(3) + 2*sqrt(4107 - 100*pi**2))/(20*pi)` (L6) | identical (L8) |
| `g_plus` | `(37*sqrt(3) + 2*sqrt(4107 - 100*pi**2))/(20*pi)` (L7) | identical (L9) |
| numeric g_- | 0.758035078944663 (notes/appendix) | 0.75803507894466282692 (L10) |
| numeric g_+ | 2.79795199200529 | 2.79795199200529341011 (L11) |
| defect | `(-12321 + 80*pi*sqrt(4107 - 100*pi**2))/(100*pi**2)`, numeric 1.7401652472273881 (L8–9) | identical, 1.74016524722738812852 (L20–21) |
| T ratio (-) | 1.3192001633911203307 (L14) | 1.31920016339112033068 (L31) |
| T ratio (+) | 0.35740427386078899977 (L15) | 0.35740427386078899977 (L32) |

No `engine_disagreement`.

## Verdict justification

The scripts hold up under attack. Attacks tried and failed: (1) **Traction-ratio tautology** — I checked whether `T_m^(±)/T_m^nat = 1/g_±` is an `x−x` self-check. It is not: both engines carry the Stage-119 background constant `C`/`cStage` as a *free positive symbol* (sympy L48, wl L47–48) and define `T_m = C/g`; the cancellation of `C` is computed, and the `.wl` exposes the pre-substitution residual `(20(gNat−1)π)/(…)` (output L29–30) which is nonzero for general `gNat` and collapses to 0 only when the equal-normalized ansatz `g_nat=1` is substituted last (wl L123/L127). The check genuinely exercises the ansatz value, consistent with the IV.3 de-tautologization on record. (2) **`g_± exact form` tautology** — not tautological because the `.wl` *derives* `g_±` by `Solve` and checks against the independent surd certificate; the `.py` derives them from the quadratic formula and checks against the same certificate; in neither case is the asserted form the same object as the derived one by construction. (3) **`expect_nonzero(comp_def)`** genuinely fails if the natural branch were on-family; the residual is a concrete nonzero surd. (4) **Symbol domains** — wl `$Assumptions` (Reals, `cStage>0`, `gNat>0`, `4107−100π²>0`) are all justified by the setup; branch sign-selection (`g−r_F1 ≶ 0`) is correct (g_-−r_F1 ≈ −1.02, g_+−r_F1 ≈ +1.02). The one finding is paper-side only: the part-IV appendix's symbolic `r_F1` reads `4107−117π²` where every script, the notes, and the appendix's own quoted numerics require `4107−100π²`. Because this is a paper↔script value mismatch needing direction from the user (and Codex must not silently edit paper prose), the verdict is `findings` with `needs_user_resolution: true`; no script-side change is warranted. I confirm I read the card, notes, and the part04 appendix rows, and that the scripts' verified claim matches the notes (the authoritative carrier) exactly.

## Value Reconciliation (pass-2 augmentation)

Saved outputs are FRESH (sympy .txt mtime 2026-05-29 13:14 > .py 12:54; mathematica .txt 2026-06-01 16:05 > .wl 15:49), so the committed outputs are the authoritative record of emitted values.

Natural carrier = the notes `.md` (the card is deliberately terse and carries no numerics). All deliverable values match the notes. The only mismatch is in the part-IV appendix's symbolic `r_F1`.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_F1 = √(4107−100π²)/(10π)` (symbolic) | py out L5 / wl out L5 | notes L49 `4107-100\pi^2` MATCH; **appendix L562 `4107-117\pi^2`** | MISMATCH (F1, value_mismatch) |
| `g_-^F1 = (2√(4107−100π²)−37√3)/(20π)` | py out L6 / wl out L8 | notes L56 | MATCH |
| `g_+^F1 = (2√(4107−100π²)+37√3)/(20π)` | py out L7 / wl out L9 | notes L56 | MATCH |
| `g_-^F1 ≈ 0.758035078944663` | py L6/wl out L10 | notes L63; appendix L571 | MATCH |
| `g_+^F1 ≈ 2.79795199200529` | py L7/wl out L11 | notes L65; appendix L573 | MATCH |
| defect `𝒞_nat = (−12321+80π√(4107−100π²))/(100π²)` | py out L8 / wl out L20 | notes L88 | MATCH |
| defect numeric `≈ 1.74016524722739` | py out L9 / wl out L21 | notes L89 | MATCH |
| `Δg_- ≈ 0.241964921055337` | py out L10 | notes L104 | MATCH |
| `Δg_+ ≈ 1.79795199200529` | py out L11 | notes L110 | MATCH |
| `T_m(-)/T_m(nat) ≈ 1.31920016339112` | py out L14 / wl out L31 | notes L134 | MATCH |
| `T_m(+)/T_m(nat) ≈ 0.357404273860789` | py out L15 / wl out L32 | notes L138 | MATCH |
| `g_nat = 1` (natural ansatz) | py L34 / wl L91 | notes L27 `𝔤_nat=1` | MATCH |

INTERNAL (scaffolding, no prose expected, no finding): banner strings; `C`/`cStage` free symbol; the symbolic pre-substitution residuals `(20(gNat−1)π)/(…)`; `expect_zero/expect_nonzero/expectZero/expectNonzero` flags; the `Solve` rootRules / branch-selection booleans; `mouthRadius = 37/20` (this is the upstream `L/a` input, also stated in notes L48/appendix L557).

reconciliation: 12 deliverable values checked, 1 misaligned (the appendix symbolic `r_F1`; numerics and notes all reconcile).

## Self-test notes

Checked: (1) Variable independence — no `sp.diff`/`D[...]` derivatives in either script, so the zero-derivative trap does not apply. (2) Branch/sign selection — verified the wl `SelectFirst` sign predicates pick the right roots (g_-−r_F1<0, g_+−r_F1>0 both hold numerically). (3) Trivial-case pre-check — the de-tautologized traction residual `(20(gNat−1)π)/(…)` is manifestly nonzero for `gNat≠1` and exactly 0 at `gNat=1`, confirming the assertion is non-trivial; `expect_nonzero(comp_def)` evaluates to a concrete nonzero surd (≈1.7402), so it genuinely could fail. (4) The single finding is paper-side (appendix value_mismatch) requiring user resolution, so no script-side `## Required change` self-test for Codex is needed; the script is correct and must not be altered.
