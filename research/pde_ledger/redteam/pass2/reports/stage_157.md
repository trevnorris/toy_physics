---
unit_id: 157
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md]
  paper_appendix: present
---

# Audit unit 157 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_157.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_157}` at line 1348; no separate prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.txt`

## What the paper claims

Stage 157 is a `\StatusNumerical / \StatusOpen` co-evolving core–mouth transport
ledger entry. Per the card body and the quote block, it "Identifies the reduced
co-evolving target and leaves the deviation-to-normalization map as the next
task." The notes flesh this out: on the analyzed positive branch, the Family-1
closure has (i) an exact compensation identity `R(g_*) = 1/4` for the
lower-compensated moment `g_* = r − √(1+r²)/2`; (ii) a carried Stage-155 frozen
canonical point (`Sigma0_star≈1.80594`, `T_hat_star≈0.901484`, `Pi_star`) that
is no longer exactly compensated; and (iii) a carried Stage-156 renormalized
canonical tuple (`Sigma0_can≈4.651034`, `T_hat_can≈1.446708`, `Pi_can≈3.871564`)
that sits strictly above the original canonical point and restores exact
compensation. The `\stagefield{Checks}` list requires: deviations taken about the
renormalized canonical point; even-preservation constraints imposed before
reading the residual odd defect; and tangent motion on the parent compensation
family giving `delta_perp = 0`. The card is explicitly OPEN: it does not assert
the full PDE realizes the target, and (by design) DEFERS the
deviation-to-normalization map to Stage 158.

## What the script claims to verify

The docstring enumerates seven checks: (1) `g=g_* ⟺ R=1/4`; (2) the self-matched
traction law `Sigma0 = 20·T_hat²/9`; (3) the carried Stage-156 renormalized tuple
satisfies the exact branch identities (`rF1` radical, `g_*` lower branch,
`R_can=1/4`, `Pi_can` identity); (4) the renormalized point is strictly above the
original canonical point; (5) tangent motion on the lower compensated family
keeps `delta R = 0`; (6) canonical-even preservation pins `deltaC = delta_kappa_W
= 0` via a non-degenerate kernel, with the tangent/family deviation-to-
normalization map DEFERRED to Stage 158; (7) the Stage-158 expansion point is the
renormalized canonical tuple (linearized `Pi` packet). Constants are loaded from
the upstream Stage-155/156 numerical harness JSON, not hardcoded in-stage.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `R(g_*) = 1/4` compensation identity | sympy:62 / wl:48 `R(g_*) − 1/4 == 0` | match |
| `g_* = r − √(1+r²)/2` lower branch | sympy:61,66,78 / wl:47,52,64 | match |
| `rF1` exact radical = carried value | sympy:65,77 / wl:51,63 (`√(4107−100π²)/(10π)`) | match |
| Self-matched traction `Sigma0 = 20 T_hat²/9` | sympy:79 / wl:65 | match |
| `R_can = 1/4`, `Pi_can = Sigma0_can(1 − S_can/4)` | sympy:81-83 / wl:67-69 | match |
| Renormalized point strictly above canonical | sympy:93-94 / wl:79-81 | match |
| Tangent on family keeps `delta R = 0` (`delta_perp=0` check) | sympy:101-102 / wl:89-90 | match |
| Even-preservation pins `deltaC = delta_kappa_W = 0` | sympy:108-111 / wl:96-100 | match |
| Even-preservation non-degeneracy (de-tautologized) | sympy:127-131 / wl:111-114 | match |
| Stage-158 expansion point = linearized `Pi` packet | sympy:135-139 / wl:119-123 | match |
| Deviation-to-normalization map | DEFERRED to Stage 158 (by design) | match (deferral honored) |

`paper_alignment: aligned`. Every card deliverable maps to a non-tautological
script check; the Stage-158 deferral is explicit on both sides.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `expect_zero(R(g_*) − 1/4)` | compensation identity | yes |
| A2 | sympy | 77 | `expect_close(rF1_exact, rF1)` | radical vs carried value | yes |
| A3 | sympy | 78 | `expect_close(g_star_exact, g_star_num)` | lower-branch moment | yes |
| A4 | sympy | 79 | `expect_close(Sigma0_can, 20/9·T_hat_can²)` | self-matched traction law | yes |
| A5 | sympy | 82 | `expect_close(R_can, 1/4)` | renormalized compensation | yes |
| A6 | sympy | 83 | `expect_close(Pi_can, Sigma0_can(1−S_can/4))` | Pi identity | yes |
| A7 | sympy | 93-94 | `assert ordering Can > Star` | renormalized above canonical | yes |
| A8 | sympy | 102 | `expect_zero(dR on family)` | tangent `delta_perp=0` | yes |
| A9 | sympy | 110-111 | `assert kernel == {0,0}` | even-preservation | yes |
| A10 | sympy | 128-131 | `expect_zero(even_det + 27σ)` | non-degeneracy (load-bearing) | yes |
| A11 | sympy | 139 | `expect_zero(Pi_lin − Pi0 − dPi_expected)` | Stage-158 packet | yes |
| B1-B11 | mathematica | 48,63,64,65,68,69,79-81,90,98-100,113-114,123 | parallel `expectZero/expectApprox/If-fail` | same claims, independent FullSimplify | yes |

All rows "yes": no tautological or orphaned assertions.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is judged a genuine second engine, not a transliteration. Justification:

1. **Different simplification machinery on shared identities.** SymPy uses
   `sp.simplify(sp.expand(...))`; the `.wl` uses
   `FullSimplify[Together[Expand[...]], Assumptions -> $Assumptions]`
   (wl:20-24). For the load-bearing algebraic identities (`R(g_*)−1/4`, tangent
   `dR`, the Stage-158 packet, the determinant), both engines must legitimately
   land on the same symbolic zero — that *is* the intended cross-check, not echoed
   choreography.

2. **Genuinely different symbol domains.** SymPy declares `sigma_star` as bare
   `real=True` (py:104), while the `.wl` constrains `0 < sigmaStar < 1`
   (wl:93) and solves `Solve[..., Reals]`. The even-preservation `Solve` route
   (wl:96) and the determinant route (wl:112) are the natural independent ways to
   establish triviality of the kernel, not a line-by-line port of py intermediate
   variables.

3. **Constants are loaded from a shared JSON, not re-derived in either engine.**
   Both read `stage155_156_fixedpoint_samples.json`; that is a carry-forward
   (the upstream Stage-155/156 harness `FixedPointHarness` is the derivation
   source), so cross-engine agreement on the *carried* numbers is expected and
   correct rather than a transliteration tell.

The structural parallelism (same 4 banners, same check set) is dictated by the
fixed list of identities the stage must verify, which is the right shape for an
algebraic-identity capstone. No `mathematica_transliteration` finding.

## Engine cross-check

Both saved outputs agree on every check:

| Check | SymPy out | Mathematica out |
|---|---|---|
| `R(g_*) − 1/4` | `= 0` (line 9) | `= 0`, PASS (line 9-10) |
| `rF1` radical diff | 1.88e-15 (l.14) | 1.95e-15, PASS (l.15-16) |
| `g_*` diff | 1.73e-16 (l.15) | 2.16e-16, PASS (l.17-18) |
| self-matched traction diff | 1.32e-15 (l.16) | 1.20e-15, PASS (l.19-20) |
| `R_can=1/4` diff | 4.36e-16 (l.17) | 4.33e-16, PASS (l.21-22) |
| `Pi_can` identity diff | 1.69e-15 (l.18) | 1.51e-15, PASS (l.23-24) |
| tangent `dR=0` | `= 0` (l.30) | `= 0`, PASS (l.36-37) |
| even-preservation kernel | `[{deltaC:0, delta_kappa:0}]` (l.31) | `{{deltaC->0, dKappa->0}}`, PASS (l.38) |
| non-degeneracy `det+27σ` | `= 0` (l.32) | `= 0`, PASS (l.39-40) |
| Stage-158 packet | `= 0` (l.37) | `= 0`, PASS (l.45-46) |
| `dPi_tan` coefficients | `−1.16276·δS + 0.83241·δΣ0` (l.38) | identical to 30 digits (l.47) |

The reported tuple (`Sigma0_star`, `T_hat_star`, `Pi_star`, `Sigma0_can`,
`T_hat_can`, `S_can`, `Pi_can`) matches to full carried precision. Engines agree.

## Verdict justification

`clean`. I read the card, the notes, and the appendix `\input` row, then attacked
the scripts. The compensation identity `R(g_*)=1/4` is non-tautological (it
genuinely substitutes `g_*` and collapses `(1+r²)/4 / (1+r²)`); the tangent
`dR=0` is the directional derivative of `R` along the family curve and is
substantive; the self-matched traction, `R_can`, and `Pi_can` identities check
the carried tuple against exact branch laws to ~1e-15; the strict-ordering assert
exercises "renormalized above canonical." The previously-flagged
de-tautologization is intact: both engines assert the *load-bearing* reason the
even-preservation kernel is trivial — `det([[1,−9σ],[5,−72σ]]) = −27σ ≠ 0` for
σ≠0 (py:127-131, wl:111-114) — rather than re-solving the homogeneous pair, and
the SymPy docstring item 6 (py:12-14) correctly states the deviation-to-
normalization map is "deferred to Stage 158," matching the card's `\StatusOpen`
deferral. The watch item is clean: the radical is `√(4107 − 100π²)/(10π)` (the
correct, non-stale `100π²`; numerically `1.778`, matching `rF1`); no `168π²` /
`168%` anywhere. The Mathematica script is an independent second engine, not a
transliteration. Outputs are fresh (txt mtime 23:02 > script mtime 22:53). No
attack succeeded.

## Self-test notes

I checked: (1) variable-independence — every `diff(R, g)` / `diff(R, r)` and the
Stage-158 `Pi` linearization depend on the differentiated variable, so no
identically-zero-derivative trap. (2) Trivial-case substitution — `R(g_*)`
reduces to `(1+r²)/4 / (1+r²) = 1/4` symbolically, and the self-matched/`Pi`
numeric identities reproduce the carried tuple to 1e-15. (3) The non-degeneracy
`det+27σ` is a genuine 2×2 determinant check that fails if the carried
coefficients lose rank — load-bearing, not the homogeneous kernel restated. No
trap fired.

## Value Reconciliation (pass-2 augmentation)

I enumerated every RESULT/deliverable value the two engines emit (constants
loaded from `stage155_156_fixedpoint_samples.json` and the derived/printed
symbolic results) and located each in the `.tex` card and `.md` notes. The notes
file is the natural carrier; the terse card legitimately omits the numbers.

Main deliverable table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `R(g_*) = 1/4` (compensation identity) | py:62 / wl:48; sympy out l.9 | notes l.16-18 (`R[Σ]` def; `R=1/4` compensation) | MATCH |
| `g_* = r − √(1+r²)/2` | py:61,66 / wl:47,52 | notes l.16-21 (lower-branch moment) | MATCH |
| `rF1 ≈ 1.77799353547498` = `√(4107−100π²)/(10π)` | py:65,67 / wl:51,53; out l.14 | notes l.20 (`r_F1≈1.77799353547498`) | MATCH |
| `Sigma0_star ≈ 1.80594111095636` | py:69 / wl:55; out l.19 | notes l.43 (`Σ_0^*≈1.80594111095636`) | MATCH |
| `T_hat_star ≈ 0.901484054174204` | py:70 / wl:56; out l.20 | notes l.46 (`T_{m,*}≈0.901484054174204`) | MATCH |
| `Sigma0_can ≈ 4.651033550168867` | py:72 / wl:58; out l.22 | notes l.61 (`Σ_0^{can}≈4.651033550168867`) | MATCH |
| `T_hat_can ≈ 1.4467083664567613` | py:75 / wl:61; out l.23 | notes l.63 (`T_{m,can}≈1.4467083664567613`) | MATCH |
| `Pi_can ≈ 3.871564377479002` | py:74 / wl:60; out l.25 | notes l.65 (`Π_{can}≈3.871564377479002`) | MATCH |
| self-matched traction `Sigma0 = 20 T_hat²/9` | py:79 / wl:65 | notes (Family-1 closure block, implied by branch); card Purpose | MATCH (symbolic law, lives in notes narrative) |
| `Pi_can = Sigma0_can(1 − S_can/4)` identity | py:83 / wl:69 | notes l.22-29 (Φ potential with `−R·T_q`, R=1/4) | MATCH |
| `dPi_tan = −1.16276 δS + 0.83241 δΣ0` (Stage-158 packet) | py:140 / wl:124; out l.38 | Stage-158 deferral; printed status output only | INTERNAL (handoff print, deferred to 158 by design) |

Carried-input / internal items (accounted for, NO finding):

- `Pi_star ≈ 1.50882951349316` (py:71/wl:57; out l.21) — Stage-155 frozen
  *canonical traction* bias, a carried INPUT used only for the ordering test
  `Pi_can > Pi_star`. It is upstream-owned (Stage 155); the notes carry the
  distinct *solved* frozen fixed point `Π_fp≈1.4885734438` (l.53). Not a
  stage-157 deliverable → reconciles at its home stage, no finding.
- `S_can ≈ 0.6703621156734617` (py:73/wl:59; out l.24) — renormalized source
  moment, an intermediate that closes the `Pi_can` identity and is derivable from
  the reported tuple (`S_can = 4(1 − Pi_can/Sigma0_can)`); the three reported
  tuple members it depends on are all in the notes. Not an independent stated
  deliverable → no finding.
- Verification scaffolding excluded per the augmentation: all `diff` residuals
  (`1e-15`-scale), pass/fail flags, kernel solution `{deltaC:0, dKappa:0}`,
  `det+27σ=0`, tolerances.

`reconciliation: complete; 11 deliverable values checked, 0 misaligned` (plus 2
carried-input constants accounted for, 0 misaligned).
