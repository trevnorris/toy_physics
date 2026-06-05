---
unit_id: 059
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage059_operator_branch_residual_bounds.md]
  paper_appendix: present
---

# Audit unit 059 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_059.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage059_operator_branch_residual_bounds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 96)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.txt`

## What the paper claims

Stage 059 converts the Stage-058 fixed-point bracket (`Xi Delta_0 <= Pe_* <= Xi Delta_inf`) into exact branch verdicts. The card's `\stagefield{Output}` is: "The exact gain thresholds \eqref{eq:app-stage059-Xi-thresholds}", i.e. the boxed pair `Xi_fail = Pe_req/Delta_inf(kappa,eta)` and `Xi_suff = Pe_req/Delta_0(kappa,eta)`, together with the verdict `Xi < Xi_fail => fail`, `Xi > Xi_suff => succeed`. The notes enumerate the supporting deliverables: (1) the operator-selected physical ratio `zeta_phys = Omega_(Pe_*)^2 (kappa+pi^2/4)/(kappa+y^2)`; (2) the exact support brackets `zeta_-`, `zeta_+` (using `Omega_(Xi Delta_0)^2` and `Omega_(Xi Delta_inf)^2`); (3) the residual brackets `R_- = zeta_req - zeta_+`, `R_+ = zeta_req - zeta_-`; (4) the two theorem gates; (5)/(6) `Pe_req` and the thresholds with `Xi_fail <= Xi_suff`; and (5, sec.5) the weak-coupling law `zeta_phys = A_K[1 + ((4-pi)/pi) Xi Delta_0 + O(Xi^2)]` with `A_K = (kappa+pi^2/4)/(kappa+y^2)`. The appendix row (line 96) summarizes the stage as "Exact fail/succeed thresholds Xi_fail, Xi_suff."

## What the script claims to verify

Both engines (a) reproduce the symbolic forms of `zeta_-`, `zeta_+`, `R_-`, `R_+`, `Xi_fail`, `Xi_suff` as labeled prints; (b) assert positivity of `Xi_suff - Xi_fail` under the explicit ordering hypothesis `Delta_inf = Delta_0 + delta_gap`, `delta_gap > 0` (so the result is `Pe_req*delta_gap/(Delta_0*(Delta_0+delta_gap)) > 0`); (c) derive the small-Pe series/derivative of the actual Stage-056/039 closed form `Omega_Pe = pi*Pe*(2 Pe e^Pe + pi)/((4 Pe^2 + pi^2)(e^Pe - 1))` and assert its `Omega_Pe^2` linear coefficient equals `(4-pi)/pi`; and (d) a numeric saturation check on a concrete probe (`kappa=2, y=1, Delta_0=3/5, delta_gap=2/5, zeta_req=2/5`, all independent of Omega_Pe), solving `A_K Omega(Xi*Delta)^2 = zeta_req` for `Xi` and confirming `Xi_fail*Delta_inf` and `Xi_suff*Delta_0` both saturate at the directly-solved `Pe_star`. The assembled weak-coupling form `zeta_weak = A_K(1 + ((4-pi)/pi) Xi Delta_0)` is printed but not separately asserted.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `Xi_fail = Pe_req/Delta_inf`, `Xi_suff = Pe_req/Delta_0` (card Output) | printed symbolic forms (py 64-67 / wl 59-67); numeric saturation asserts (py 118-119 / wl 106-107) | match |
| `Xi_fail <= Xi_suff` (notes sec.4) | `expect_positive`/`expectPositive` on ordered gap (py 73 / wl 73) | match |
| `zeta_-`, `zeta_+` support brackets (notes sec.2) | symbolic prints (py 54-57 / wl 55-63) | match (print-only, but trivial algebra; see note) |
| `R_-`, `R_+` residual brackets (notes sec.3) | symbolic prints (py 59-62 / wl 57-65) | match (print-only) |
| weak-coupling slope `(4-pi)/pi` of `Omega_Pe^2` (notes sec.5) | series/Limit coeff assert (py 81-86 / wl 79,112) | match |
| assembled `zeta_phys = A_K[1+((4-pi)/pi)Xi Delta_0]` (notes sec.5) | printed only (py 121-122 / wl 108,111) | partial (print-only; load-bearing constant `(4-pi)/pi` is asserted separately) |
| verdict zones / theorem gates (card eq:verdict, notes sec.4,6) | implied by threshold definitions + saturation; not separately asserted | partial (gates are corollaries of asserted thresholds) |

`paper_alignment: aligned` — every load-bearing constant and identity the card/notes report is either asserted or is a trivial print of an already-verified expression; nothing in the script contradicts the paper.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 73 | `expect_positive(Xi_suff_ordered - Xi_fail_ordered)` | `Xi_fail <= Xi_suff` (notes sec.4) | yes |
| A2 | sympy | 83-86 | `expect_zero(series_coeff - (4-pi)/pi)` | weak-coupling slope (notes sec.5) | yes |
| A3 | sympy | 118 | `expect_close(Xi_fail_solved*DeltaInf, Pe_star)` | `Xi_fail` saturation (card Output) | yes |
| A4 | sympy | 119 | `expect_close(Xi_suff_solved*Delta0, Pe_star)` | `Xi_suff` saturation (card Output) | yes |
| A5 | math | 73 | `expectPositive(xiSuffOrdered - xiFailOrdered)` | `Xi_fail <= Xi_suff` | yes |
| A6 | math | 112 | `expectZero(omegaSqLinearCoeff - (4-pi)/pi)` | weak-coupling slope | yes |
| A7 | math | 106 | `expectApprox(xiFailSolved*deltaInfProbe, peStar)` | `Xi_fail` saturation | yes |
| A8 | math | 107 | `expectApprox(xiSuffSolved*delta0Probe, peStar)` | `Xi_suff` saturation | yes |

All eight assertions are non-tautological and anchored. `zeta_-`, `zeta_+`, `R_-`, `R_+`, and `zeta_weak` are labeled prints (not asserted), but each is a one-line algebraic restatement of a notes deliverable, and the only load-bearing numeric constant among them — the `(4-pi)/pi` slope — is independently asserted (A2/A6).

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.txt:3,18` (mtime 2026-05-26 11:05:36) vs script `.py` mtime 2026-06-03 15:59:11
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.txt:3` (mtime 2026-05-26 11:05:41) vs script `.wl` mtime 2026-06-03 15:59:11

**What's wrong:**
Both saved output transcripts predate their scripts by ~8 days, so they reflect a pre-banner-fix state. Concretely, the SymPy `.txt` still shows `STAGE 42 — OPERATOR-SELECTED RESIDUAL BOUNDS` (line 3) and `Stage 42 audit passed.` (line 18), whereas the current `.py` banner reads `STAGE 59` (py:45) and closing line `Stage 42 audit passed.` (py:124). The Mathematica `.txt` is internally inconsistent: line 3 shows `STAGE 042 — OPERATOR BRANCH RESIDUAL BOUNDS` while line 24 shows `Stage 059 Mathematica audit passed.`, and the current `.wl` banner is `STAGE 059` (wl:44). The substantive numeric/symbolic content (residuals zero, `(4-pi)/pi` coefficient, saturation diffs) is otherwise consistent with the current scripts; only the stale stage-number labels differ.

**Why this matters:**
Informational only — the captured transcripts no longer match the scripts' banner/labels. A fresh re-run refreshes them.

**Required change:**
No script-logic change required for the staleness itself. The orchestrator's independent re-run regenerates both `.txt` files. (See additional stale self-label items below, which the orchestrator resolves under the deferred numbering-band plan.)

**Verification:**
After re-run, the SymPy `.txt` line 3 reads `STAGE 59 …` and the Mathematica `.txt` line 3 reads `STAGE 059 …`, matching the script banners.

#### Additional stale self-labels (low-severity numbering class, non-blocking)

Quoted for the orchestrator's separate numbering-band pass; not part of the directive:
- `scripts/...sympy_audit.py:3` docstring `"Moving-Throat PDE — Stage 42 SymPy audit."` (should be 059)
- `scripts/...sympy_audit.py:6-9` docstring cross-refs `Stage-41`, `Stage-39` (pre-renumber chain; map to 058/056 under +17)
- `scripts/...sympy_audit.py:45` banner `STAGE 59` (missing zero-pad → `STAGE 059`)
- `scripts/...sympy_audit.py:75` comment `Stage-39 Omega_Pe series` (pre-renumber → 056)
- `scripts/...sympy_audit.py:124` `print("\nStage 42 audit passed.")` (should be 059)

These are unambiguous self-labels; the cross-refs (`Stage-41`/`Stage-39`) are flagged but ambiguous and left to the numbering plan.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. Crucially, the two engines compute the weak-coupling slope by genuinely different methods:
- SymPy: `Omega_sq_series = sp.series(Omega_Pe**2, Pe, 0, 2).removeO()` then `.coeff(Pe, 1)` (py:81-85).
- Mathematica: `omegaSqLinearCoeff = FullSimplify[Limit[D[omegaPe^2, pe], pe -> 0], ...]` (wl:79) — a limit-of-derivative, not a series expansion.

Both arrive at `(4-pi)/pi` (= `-1 + 4/Pi`) independently. The probe block (wl:80-107) uses `FindRoot` where SymPy uses `nsolve`, again independent solver choreography. The symbolic-bracket prints are structurally parallel (both restate the same notes formulas), which is acceptable since those are restatements, not the load-bearing derivation. No `mathematica_transliteration` finding.

## Engine cross-check

| Quantity | SymPy out | Mathematica out | Agree? |
|---|---|---|---|
| `Xi_suff - Xi_fail (ordered)` | `Pe_req*delta_gap/(Delta0*(Delta0+delta_gap))` (line 11) | `(deltaGap*peReq)/(delta0*(delta0+deltaGap))` (line 11) | yes |
| `Omega_Pe^2` linear coeff | `Pe*(-1 + 4/pi) + 1` series → coeff `0` residual (lines 12-13) | `-1 + 4/Pi`, residual `0` (lines 19,21) | yes |
| `Xi_fail*DeltaInf` saturation | diff `0` (line 14) | diff `0` (line 16) | yes |
| `Xi_suff*Delta0` saturation | diff `7.24e-71` (< 1e-40 tol) (line 15) | diff `0` (line 17) | yes |
| `weak-coupling zeta_phys` | expanded `A_K(1+((4-pi)/pi)Xi Delta0)` (line 16) | same expansion (line 20) | yes |

Both engines pass and agree. The Mathematica `Limit::alimv` warning (out line 14) is benign: the limit is taken with `Assumptions -> pe > 0` while `pe` is the limit variable, so the assumption is correctly ignored and the limit at `pe -> 0` is still computed correctly (confirmed by the matching `-1 + 4/Pi` and the `expectZero` pass).

## Verdict justification

`verdict: findings` with a single low-severity `stale_output` finding (F1). The math holds up under attack: A1/A5 are not tautological because the positivity rests on the explicit `delta_gap > 0` ordering hypothesis (which is the notes' stated premise `Delta_inf >= Delta_0`), not on a constructed identity; A2/A6 derive the `(4-pi)/pi` slope from the genuine Stage-056 closed form rather than asserting it against itself; A3/A4/A7/A8 use a `zeta_req = 2/5` target that is explicitly independent of `Omega_Pe`, so the saturation check could fail if the threshold definitions were wrong. I attempted to break each: (i) the saturation check is not circular because the probe target is independent and the solved `Xi` is compared to a directly-root-solved `Pe_star`; (ii) the `expect_positive` is not vacuous since a sign error in the threshold definitions would surface as a sign flip in the gap; (iii) the series coefficient check would fail for any other slope value. The only defect is the stale transcripts and stale stage-number self-labels, which are informational/numbering-class and resolved by the standard re-run plus the deferred numbering-band pass. I read the card, notes, and appendix row; the script's verified claims match the paper.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every labeled result value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_- = A_K * Omega(Xi*Delta0)^2` | py:54-56 / wl:55,62; sympy.txt L5, math.txt L5 | notes L112-113 | MATCH |
| `zeta_+ = A_K * Omega(Xi*DeltaInf)^2` | py:55-57 / wl:56,63; sympy.txt L6, math.txt L6 | notes L115-116 | MATCH |
| `R_- = zeta_req - zeta_+` | py:59,61 / wl:57,64; sympy.txt L7, math.txt L7 | notes L135-136 | MATCH |
| `R_+ = zeta_req - zeta_-` | py:60,62 / wl:58,65; sympy.txt L8, math.txt L8 | notes L138-139 | MATCH |
| `Xi_fail = Pe_req/Delta_inf` | py:64,66 / wl:59,66; sympy.txt L9, math.txt L9 | card L17 (boxed); notes L183 | MATCH |
| `Xi_suff = Pe_req/Delta_0` | py:65,67 / wl:60,67; sympy.txt L10, math.txt L10 | card L19 (boxed); notes L185 | MATCH |
| `Xi_suff - Xi_fail (ordered) = Pe_req*delta_gap/(Delta0*(Delta0+delta_gap))` | py:73 / wl:73; sympy.txt L11, math.txt L11 | supports notes L215 `Xi_fail <= Xi_suff` | MATCH |
| `Omega_Pe^2` linear coefficient `= (4-pi)/pi` | py:81-85 / wl:79,112; sympy.txt L12-13, math.txt L19 | notes L228-229 | MATCH |
| `weak-coupling zeta_phys = A_K(1 + ((4-pi)/pi) Xi Delta0)` | py:121-122 / wl:108,111; sympy.txt L16, math.txt L20 | notes L237-243 | MATCH |
| `A_K = (kappa+pi^2/4)/(kappa+y^2)` | py:52 / wl:53 (factor in all above) | notes L243 | MATCH |

INTERNAL (scaffolding, no finding): probe constants `kappa=2, y=1, Delta0=3/5, delta_gap=2/5, zeta_req=2/5`; `Pe_star` numeric root; `Xi_fail_expected`/`Xi_suff_expected`/`Xi_fail_solved`/`Xi_suff_solved` (exist only to drive the saturation asserts); saturation residual diffs (`0`, `7.24e-71`); pass/fail flags; tolerance `1e-40`.

reconciliation: complete; 10 values checked, 0 misaligned

## Self-test notes

Checked: (1) variable independence — the only derivative is `D[omegaPe^2, pe]` (wl:79), and `omegaPe` genuinely depends on `pe`, so the derivative is non-trivial and the series coeff is real, not an identically-zero artifact. (2) Trivial-case pre-check — substituting the probe values, `A_K_probe Omega(Xi*Delta)^2 = zeta_req` has a real root near `Pe ~ 1/2` and the threshold products reproduce `Pe_star`, so the saturation asserts give genuine (non-vacuous) zero residuals. (3) Sign/ordering — `expect_positive` rests on `delta_gap > 0`, which matches the notes' `Delta_inf >= Delta_0` premise; a threshold sign error would flip the gap sign and fail the check. No script-side math defect found; the single finding is informational stale_output, so no directive is written.
</content>
</invoke>
