---
unit_id: 086
batch: III.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage086_family1_loading_ratio_window.md]
  paper_appendix: present
---

# Audit unit 086 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_086.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage086_family1_loading_ratio_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 150 and 290)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.txt`

## What the paper claims

Stage 086 converts the Stage-085 Family-1 cancellation theorem into a *pure window for `rho_alpha`*. The `\stagefield{Output}` is the Family-1 loading-ratio window at `lambda_mu=1`, `eps_blk=0`:

- `rho_alpha <= 3.46622291347846` → guaranteed success (eq. app-stage086-success-window);
- `rho_alpha >= 3.46752913273870` → guaranteed failure;
- `rho_alpha < 3.46752922945601` → absolute constructive ceiling (eq. app-stage086-fail-window).

The notes additionally state (a) the structural identity `rho = Q(zeta;eps_blk) = (1+(1-2 eps) zeta)/(1 - eps zeta)`, (b) the unblocked reduction `Q(zeta;0) = 1+zeta`, (c) the blocking cap `eps_blk < 1/zeta_max^(F1) ≈ 0.405263689711371`, and (d) monotonicity `d rho_max / d eps_blk = zeta_max (zeta_max - 1)/(1 - eps_blk zeta_max)^2 > 0`. The four numerical values are obtained from carry-forward zetas via `rho = 1 + zeta` (eps=0).

## What the script claims to verify

Both scripts treat the four `zeta` values as upstream constants from Stages 63–64 and verify the Q-map machinery used to convert them into rho-values: the unblocked reduction `Q(zeta;0) = 1+zeta`, the closed-form derivatives `dQ/dzeta = (1-eps)/(1-eps*zeta)^2` and `d rho_max / d eps = zeta_max(zeta_max-1)/(1-eps*zeta_max)^2`, and the numeric values `rho = 1+zeta` for each of the four carry-forward zetas plus the `eps_blk` cap. Mathematica asserts the two derivative formulas via symbolic equality and the four rhos via `expectApprox` to 10^-14; SymPy prints the derivatives without an assertion and uses redundant `expect_zero(rho_X - (1+zeta_X))` checks.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `Q(zeta;0) = 1+zeta` unblocked reduction | sympy L28, mma L46 | match |
| `dQ/dzeta = (1-eps)/(1-eps*zeta)^2` (notes implicit) | mma L45 (formula equality) | match (sympy only prints, no assertion) |
| `rho_suff^(chi)(1;0) ≈ 3.46622291347846` | sympy L46, mma L63 | match (both tautological — see F1, F2) |
| `rho_fail^(chi)(1;0) ≈ 3.46752913273870` | sympy L47, mma L64 | match (both tautological) |
| `rho_suff^(J)(1;0) ≈ 3.44257571477179` | sympy L48, mma L65 | match (both tautological) |
| `rho_max^(F1)(0) ≈ 3.46752922945601` | sympy L49, mma L66 | match (both tautological) |
| `eps_blk cap = 1/zeta_max ≈ 0.405263689711371` | sympy L52 (print only), mma L72 | partial (sympy lacks assertion) |
| `d rho_max / d eps_blk > 0` exact formula | sympy L57 (print only), mma L71 | partial (sympy lacks assertion) |

`paper_alignment: aligned` — every paper-side deliverable has at least one script-side check, the numerical constants in the script match the paper's 14-digit literals exactly, and the Mathematica targets at extended precision are just floating-point extensions of those same literals. No `paper_misalignment` triggered.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 28 | `expect_zero(Q(zeta;0)-(1+zeta))` | unblocked reduction | yes |
| A2 | sympy | 46 | `expect_zero(rho_suff_chi - (1+zeta_suff_chi))` | rho_suff^(chi) value | no (tautological — F1) |
| A3 | sympy | 47 | `expect_zero(rho_fail_chi - (1+zeta_fail_chi))` | rho_fail^(chi) value | no (tautological — F1) |
| A4 | sympy | 48 | `expect_zero(rho_suff_J - (1+zeta_suff_J))` | rho_suff^(J) value | no (tautological — F1) |
| A5 | sympy | 49 | `expect_zero(rho_max - (1+zeta_max))` | rho_max value | no (tautological — F1) |
| A6 | mma | 45 | `expectZero[dQ - (1-eps)/(1-eps*zeta)^2]` | dQ/dzeta formula | yes |
| A7 | mma | 46 | `expectZero[Q(zeta;0) - (1+zeta)]` | unblocked reduction | yes |
| A8 | mma | 63 | `expectApprox[rhoSuffChi, 3.466222913...]` | rho_suff^(chi) value | no (tautological — F2) |
| A9 | mma | 64 | `expectApprox[rhoFailChi, 3.467529132...]` | rho_fail^(chi) value | no (tautological — F2) |
| A10 | mma | 65 | `expectApprox[rhoSuffJ, 3.442575714...]` | rho_suff^(J) value | no (tautological — F2) |
| A11 | mma | 66 | `expectApprox[rhoMaxNum, 3.467529229...]` | rho_max value | no (tautological — F2) |
| A12 | mma | 71 | `expectZero[d rho_max/d eps - closed-form]` | monotonicity formula | yes |
| A13 | mma | 72 | `expectApprox[1/zetaMaxNum, 0.40526368971...]` | eps_blk cap | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py:36-49`

**What's wrong:**
The SymPy script defines

```
rho_suff_chi = sp.N(Q.subs({zeta: zeta_suff_chi, eps: 0}), 30)
...
rho_max     = sp.N(Q.subs({zeta: zeta_max,       eps: 0}), 30)
```

at lines 36–39 and then asserts (lines 46–49)

```
expect_zero("rho_suff^(chi) - (1+zeta_suff)", sp.N(rho_suff_chi - (1 + zeta_suff_chi), 20))
expect_zero("rho_fail^(chi) - (1+zeta_fail)", sp.N(rho_fail_chi - (1 + zeta_fail_chi), 20))
expect_zero("rho_suff^(J) - (1+zeta_suff_J)", sp.N(rho_suff_J  - (1 + zeta_suff_J),   20))
expect_zero("rho_max - (1+zeta_max)",         sp.N(rho_max     - (1 + zeta_max),      20))
```

These cannot fail. The script already verified at line 28 that `Q(zeta;eps=0) - (1+zeta) == 0` symbolically. Since `rho_X` is defined as `Q.subs({zeta: zeta_X, eps: 0})`, by direct substitution `rho_X = 1 + zeta_X` and therefore `rho_X - (1+zeta_X)` is identically zero. The assertion exercises no new physics — it re-derives the same reduction with concrete numeric `zeta_X` substituted in, after the symbolic reduction has already been certified. No upstream value, no value of `Q`, no transcription of `zeta_X` can make these four assertions fail (a transcription error in `zeta_X` would silently propagate to both sides and still give zero).

**Why this matters:**
These four "checks" do not anchor the SymPy script to the paper's numerical window. If a future edit broke the `zeta` literals or swapped them between variables, these assertions would still pass; the script would silently report `rho_max = 3.466…` while paper says `3.4675…`. The SymPy side currently does **not** anchor any of the four numeric rho-values to the paper-stated constants.

**Required change:**
Replace the four tautological assertions at lines 46–49 with `expect_zero` checks against the paper-stated numeric targets, mirroring Mathematica's `expectApprox`. Concretely, add a helper `expect_close(name, value, target, tol)` (or reuse `expect_zero` on `sp.N(value - target, 30)` against a small tolerance) and assert each `rho_X` differs from the corresponding 14-digit paper constant by at most `1e-13`.

**Verification:**
After the patch, the SymPy transcript must contain four new lines anchoring `rho_suff_chi`, `rho_fail_chi`, `rho_suff_J`, `rho_max` to the paper-stated numerical values `3.46622291347846`, `3.46752913273870`, `3.44257571477179`, `3.46752922945601` respectively. Manually flipping `zeta_suff_chi` and `zeta_fail_chi` in the script must now make the script exit nonzero.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.wl:48-66`

**What's wrong:**
The Mathematica `expectApprox` checks at lines 63–66 compare `rhoSuffChi`, `rhoFailChi`, `rhoSuffJ`, `rhoMaxNum` against extended-precision targets such as `3.46622291347846012143918414949`. The script computes each via `N[qMap /. {zeta -> zetaSuffChi, epsBlk -> 0}, 30]`. Because `qMap /. epsBlk -> 0` symbolically reduces to `1 + zeta` (verified at line 46), `rhoSuffChi` is exactly `1 + zetaSuffChi`. The targets at lines 63–66 are the floating-point representations of `1 + 2.46622291347846`, `1 + 2.46752913273870`, `1 + 2.44257571477179`, `1 + 2.46752922945601` carried to 30 digits — i.e., the same arithmetic the script just performed. The "diff" is therefore at machine-precision noise (the captured output shows ~1e-16 to 1e-17), guaranteed to be below the 1e-14 tolerance regardless of whether the upstream `zeta` literals are correct.

**Why this matters:**
Like F1, these checks do not actually pin `zetaSuffChi` etc. to the paper's stated 14-digit zeta values from Stages 63/64. If the literal at line 48 were mistyped (say `2.46622291347946`), then `rhoSuffChi` and the target `3.466…` would both shift consistently and the check would still pass. The Mathematica audit lacks an anchor that the four `zeta` literals match the paper-stated 14-digit constants.

**Required change:**
Add four `expectApprox` (or `expectZero` against numeric difference) checks that compare the four `zeta` literals directly to the paper-stated 14-digit zeta values from Stage 63/64 (i.e., assert `Abs[zetaSuffChi - 2.46622291347846] <= 10^-14`, and likewise for the other three). Place these immediately after lines 48–51, before the rho computations. Keep the existing rho `expectApprox` lines — they are now non-tautological given the upstream zeta anchor.

**Verification:**
After the patch, the Mathematica transcript must contain four new `zeta` numeric-anchor PASS lines preceding the `rho_*` checks. Manually perturbing any of `zetaSuffChi`, `zetaFailChi`, `zetaSuffJ`, `zetaMaxNum` by 1e-13 must now make the script exit 1.

## Independent-derivation check (Mathematica)

The Mathematica script is not a transliteration of the SymPy script. It introduces two extra symbolic-identity assertions that SymPy lacks:

- `expectZero["dQ exact formula", dQ - (1 - epsBlk)/(1 - epsBlk*zeta)^2]` (line 45)
- `expectZero["d rho_max / d eps exact formula", dQMax - zetaMax*(zetaMax - 1)/(1 - epsBlk*zetaMax)^2]` (line 71)

It also uses `expectApprox` against extended-precision targets rather than re-evaluating the same `(1+zeta)` reduction with numeric substitution. The structural choreography differs (Mathematica binds derivative closed-forms by hand-written formula and checks via FullSimplify; SymPy uses `simplify(diff(...))` and only prints). This satisfies the independent-engine requirement.

## Engine cross-check

Both engines compute the same four rho values; both transcripts agree to all displayed digits (`3.46622291347846012143918414949`, etc.). Both engines confirm `Q(zeta;0) - (1+zeta) = 0`. The Mathematica-only derivative-formula checks are not contradicted by SymPy — SymPy prints `dQ/dzeta = (1 - eps)/(eps**2*zeta**2 - 2*eps*zeta + 1)` which is algebraically identical to `(1-eps)/(1-eps*zeta)^2`. Engines agree.

## Verdict justification

The Mathematica script substantively verifies the structural paper claims (Q-formula reduction, derivative closed forms, eps_blk cap) and presents the four rho-window values, but the four `expectApprox` numeric "checks" on rho-values are tautological because `rho = 1+zeta` symbolic reduction already certified by the script makes the comparison trivially true. The SymPy script suffers a worse version of the same defect: its four numeric `expect_zero` checks are pure tautologies, and it has *no* assertion-level anchor at all between the script's `zeta` literals and the paper's stated 14-digit values. Both findings are correctable by adding direct numeric anchors on the upstream `zeta` literals. Paper alignment is otherwise correct: constants match, formula matches, deliverable set covered. No `stop_cold` warranted.

Attacks I tried that failed: (1) sign error in `Q` — symbolic check at L28 / L46 catches it; (2) wrong derivative formula — Mathematica L45 catches it; (3) wrong monotonicity sign — Mathematica L71 catches it; (4) paper says `rho < 3.4675292...` ceiling but script only verifies equality, not inequality — confirmed this is acceptable since the ceiling value is a definitional output, not an inequality to test; (5) paper says `eps_blk < 1/zeta_max`, script uses `eps_blk*zeta < 1` assumption — consistent.

## Self-test notes

(1) Variable independence: I checked that `qMap` depends on both `zeta` and `epsBlk`, so `D[qMap, zeta]` and `D[qMax, epsBlk]` are nontrivial (qMax substitutes `zeta -> zetaMax` so it only depends on `epsBlk` and the symbol `zetaMax`, and `D[qMax, epsBlk]` is the meaningful derivative — correct). (2) The proposed F1/F2 fixes use `expect_close` / `expectApprox` against fixed numeric literals, not derivatives, so the derivative-trivialization trap doesn't apply. (3) Trivial-case mental check: for `F1`, plugging `zeta_suff_chi = 2.46622291347846` into the new check gives `|3.46622291347846 - (1+2.46622291347846)| = 0 < 1e-13`, pass; for F2, plugging `zetaSuffChi = 2.46622291347846` gives `|2.46622291347846 - 2.46622291347846| = 0 < 1e-14`, pass. (4) Paths in directive use `scripts/` for `.py` and `mathematica/` for `.wl`. (5) The new checks do not introduce constants beyond those already stated in the paper card (lines 24–28 of stage_086.tex) and notes section 2 — no new paper-misalignment risk.
