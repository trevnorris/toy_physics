---
unit_id: 055
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage055_explicit_lowest_lane_reachability.md]
  paper_appendix: present
---

# Audit unit 055 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_055.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage055_explicit_lowest_lane_reachability.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 88, plus `\input` line 228)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt`

## What the paper claims

Stage 055 combines the Stage-053 overlap boost (ceiling `A_I <= pi^2/4`, via `Omega_exp(alpha)^2`) with the Stage-054 Robin-compliance softening (ceiling `A_K <= 4/(4-x)`) into a single lowest-lane reachability window. The card factors the physical support ratio as `zeta_0^phys = A_K(eta) Omega_0^2` (eq. app-stage055-zeta-factor) and states the boxed `\stagefield{Output}`: the combined reachability window `1 <= zeta_0^phys <= (pi^2/4)(4/(4-x))` (eq. app-stage055-reachability-window), i.e. algebraically `pi^2/(4-x)`. The notes add three further deliverables consistent with the card: (i) the explicit combined family `zeta_0^(exp+R) = Omega_exp(alpha)^2 / [1 - x/4 + x y(eta)^2/pi^2]` with `y tan y = eta`; (ii) the exact reachability floor `x >= 4 - pi^2/zeta_req`; and (iii) the equivalent stiffness-ratio form `K_X/K_W^eff <= pi^2/(4 zeta_req)`, plus the A/B/C regime split. The appendix row (line 88) summarizes the stage as "Combined overlap/compliance reachability window for the lowest support lane," status \StatusExactClosure. The lower endpoint is the symmetric twin (`alpha=0, eta=+infinity -> zeta=1`); the upper endpoint is the closure (`alpha->+infinity, eta->0+ -> (pi/2)^2 * 4/(4-x) = pi^2/(4-x)`).

## What the script claims to verify

Both engines build `Omega_exp(alpha)`, `A_K(y,x)`, and `zeta_family = Omega_exp^2 * A_K` symbolically, then assert: the symmetric twin value equals 1 (`alpha->0`, `y->pi/2`, i.e. `eta->+inf`); the closure maximum equals `pi^2/(4-x)` (`alpha->infinity`, `y->0`); the reachability floor `x_floor = 4 - pi^2/zeta_req` solved from `zeta_max = zeta_req`; and the stiffness-ratio identity `(1 - x/4)|_{x=x_floor} = pi^2/(4 zeta_req)`. The Mathematica script additionally asserts the round-trip `zeta_max(x_floor) - zeta_req == 0`. The regime split (A/B/C) and the final ledger summary are printed, not asserted (they are corollaries of the endpoint values). This is what the verdict applies to.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| factorization `zeta_0^phys = A_K Omega_0^2` | `zeta_family = Omega_exp**2 * AK` (py:39 / wl:43) | match |
| lower endpoint `zeta = 1` (twin) | `expect_zero("twin value", zeta_twin - 1)` (py:49 / wl:56) | match |
| upper endpoint `pi^2/(4-x)` | `expect_zero("closure maximum", zeta_max - pi**2/(4-x))` (py:50 / wl:57) | match |
| window `1 <= zeta <= pi^2/(4-x)` (boxed Output) | endpoints verified; monotone interior argued upstream (053/054), not re-derived here | partial |
| reachability floor `x >= 4 - pi^2/zeta_req` | `expect_zero("x floor = 4 - pi^2/zeta_req", ...)` (py:56 / wl:58) | match |
| stiffness-ratio form `K_X/K_W^eff <= pi^2/(4 zeta_req)` | `expect_zero("KX/KW equivalence", (1/AK)|_{y=0,x=x_floor} - pi**2/(4 zeta_req))` (py:60 / wl:60) | match |
| A/B/C regime split | printed corollary (py:67–69) | match (informational) |

The `\stagefield{Output}` box `(pi^2/4)(4/(4-x))` and the notes' `pi^2/(4-x)` are the same number; the script's `pi**2/(4-x)` agrees with both. The "window interior" row is `partial` only in the standard sense that the script verifies the two endpoints and the monotone interpolation is carried from Stages 053/054 (where `Omega_exp` and `A_K` ranges are established); this is the consistent convention of the ledger and not a defect against this card.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 49 | `expect_zero(zeta_twin - 1)` | lower endpoint = 1 | yes |
| A2 | sympy | 50 | `expect_zero(zeta_max - pi**2/(4-x))` | upper endpoint | yes |
| A3 | sympy | 56 | `expect_zero(x_floor - (4 - pi**2/zeta_req))` | reachability floor | yes |
| A4 | sympy | 60 | `expect_zero((1/AK)|_{y=0,x=x_floor} - pi**2/(4 zeta_req))` | stiffness-ratio form | yes |
| A5 | math | 56 | `expectZero[zetaTwin - 1]` | lower endpoint | yes |
| A6 | math | 57 | `expectZero[zetaMax - Pi^2/(4-x)]` | upper endpoint | yes |
| A7 | math | 58 | `expectZero[xFloor - (4 - Pi^2/zetaReq)]` | reachability floor | yes |
| A8 | math | 59 | `expectZero[(zetaMax /. x->xFloor) - zetaReq]` | floor round-trip | yes |
| A9 | math | 60 | `expectZero[((1/aK)/.y->0/.x->xFloor) - Pi^2/(4 zetaReq)]` | stiffness-ratio form | yes |

All assertions are non-tautological: each target value (`1`, `pi^2/(4-x)`, `4 - pi^2/zeta_req`, `pi^2/(4 zeta_req)`) is an independently-stated closed form, and the LHS is constructed from `Omega_exp`/`A_K`/`zeta_family` (built from the upstream functional forms, not from the target). The `x_floor` in A3 is obtained by `solve`, then A4/A9 substitute it back into an *independent* expression `1 - x/4`, giving a genuine chained check rather than `x == x`.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt:1` (mtime May 26 03:02) vs script mtime Jun 3 15:59
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt:1` (mtime May 26 03:02) vs script mtime Jun 3 15:59

**What's wrong:**
Both committed transcripts predate their scripts and their content diverges from what the current scripts emit. The SymPy transcript header reads `STAGE 38 — EXPLICIT LOWEST-LANE REACHABILITY` (line 3) while the current script banner is `STAGE 55` (py:25). The Mathematica transcript header reads `STAGE 038 — ...` (line 3) while the current script banner is `STAGE 055` (wl:32). The numerical/symbolic results in the transcripts are otherwise consistent with the current math (twin = 1, closure max = `pi^2/(4-x)`, x_floor = `4 - pi^2/zeta_req`, all PASS), so the divergence is the stale banner label, not a result change.

**Why this matters:**
The committed output no longer reflects the current script banner; the verifier's fresh re-run will regenerate it. Informational only — no result is wrong.

**Required change:**
Re-run both scripts to refresh the saved transcripts (the orchestrator's independent re-run does this).

**Verification:**
After re-run, both transcript headers read `STAGE 55`/`STAGE 055` and mtimes postdate the scripts.

### F2 — stale_output (numbering self-labels)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:3` — docstring `"Stage 38 SymPy audit: ..."`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:73` — `"... reaches the Stage-35 threshold iff ..."`

**What's wrong:**
Two stale pre-renumber self-labels survive in the SymPy script. Line 3 says `Stage 38` (canonical is 055; `38 + 17 = 55`). Line 73 says `Stage-35 threshold`, but the threshold is the Stage-052 demand `zeta_req` per the notes ("reaches a Stage-052 demand `zeta_req`", "Stage 052 demands"); `35 + 17 = 52`. The banner (py:25) and the Mathematica script (wl:32) already carry the correct `STAGE 55`/`055` label, so this is an incomplete renumber of the two text strings only, with no effect on the math.

**Why this matters:**
Self-label drift only; the known low-severity `stale_output` numbering class. Per the in-loop Reading-2 policy, a stage that already has a (verdict:findings) directive gets its unambiguous self-labels fixed.

**Required change:**
py:3 `Stage 38` -> `Stage 055`; py:73 `Stage-35 threshold` -> `Stage-052 threshold`. (Both are unambiguous content-keyed renumbers via the established +17 chain confirmed by the notes; not an offset sweep.) Refresh the SymPy transcript afterward.

**Verification:**
`grep -n "Stage 38\|Stage-35" <py>` returns nothing; line 3 reads `Stage 055`, line 73 reads `Stage-052 threshold`.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses Mathematica-native `FullSimplify`/`Together`/`Limit` with `$Assumptions`, a distinct `expectZero` that strips `ConditionalExpression` wrappers (wl:20–30), and a different limit choreography: the closure max is taken as a *nested* limit `Limit[Limit[zetaFamily, alpha->Infinity], y->0, Direction->"FromAbove"]` (wl:50), whereas SymPy uses `sp.limit(zeta_family, alpha, sp.oo).subs(y, 0)` (py:46). The `.wl` also adds an assertion absent from the `.py`: the floor round-trip `expectZero["zeta_max(x_floor) - zeta_req", ...]` (wl:59). The shared functional definitions of `Omega_exp`/`A_K` are dictated by the upstream physics (Stages 053/054), so their commonality is expected, not echoing. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at the level they claim:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| twin value | `1` (PASS, residual 0) | `1` (PASS) |
| closure maximum | `-pi**2/(x-4)` = `pi^2/(4-x)` (residual 0) | `Pi^2/(4-x)` (PASS) |
| x_floor | `4 - pi**2/zeta_req` (residual 0) | `ConditionalExpression[4 - Pi^2/zetaReq, 4 zetaReq > Pi^2]` -> stripped -> 0 (PASS) |
| KX/KW equivalence | residual 0 | `0` (PASS) |
| floor round-trip | (not checked in py) | `0` (PASS) |

`A_K` is printed in two algebraically-equal forms (`4*pi**2/(4*x*y**2 + pi**2*(4 - x))` in SymPy vs `(1 + x*(-1/4 + y^2/Pi^2))^(-1)` in Mathematica); both equal `1/(1 - x/4 + x y^2/pi^2)`. The `Limit::alimv` warnings in the Mathematica transcript (lines 9–15) are benign — Mathematica ignores the limit-variable assumptions, which does not affect these one-sided endpoint limits, and all five checks PASS.

## Verdict justification

The scripts hold up against the paper. Every boxed/notes deliverable — the factorization, both window endpoints (`1` and `pi^2/(4-x)`), the reachability floor `4 - pi^2/zeta_req`, and the stiffness-ratio form `pi^2/(4 zeta_req)` — is exercised by a non-tautological, well-anchored assertion in both engines, and the two engines agree. Attacks I tried that failed: (a) checking whether the closure-max assertion is trivially satisfied — it is not, the endpoint is a genuine nested limit of a non-trivial expression; (b) checking whether `x_floor` substitution is circular — it is not, A4/A9 substitute the solved floor into an independent expression `1 - x/4`; (c) checking the twin-point `y -> pi/2` substitution against the notes' `eta -> +infinity` — it is the correct principal-branch image; (d) checking the box `(pi^2/4)(4/(4-x))` against the script's `pi^2/(4-x)` — algebraically identical. The verdict is `findings` only because of two low-severity `stale_output` items: the committed transcripts predate the scripts with stale `STAGE 38/038` banners (F1), and two stale pre-renumber self-labels remain in the SymPy source (`Stage 38` at py:3, `Stage-35` at py:73; F2). No result is wrong; `material_change` is false. I read the paper card, the notes, and the appendix row, and the script's verified claim matches the paper's stated claim.

## Self-test notes

I checked: (1) variable independence — no spurious zero-derivatives are introduced (no new `diff` checks proposed; F2 is a string-label edit only). (2) Symmetry/parity — n/a, no symmetric-domain integrals in scope. (3) Trivial-case pre-check — the existing endpoint limits reduce correctly (`Omega_exp(0)=1`, `Omega_exp(inf)=pi/2`, `A_K(y=pi/2)=1`, `A_K(y=0)=4/(4-x)`), confirming the asserted endpoint values. (4) Path specifications — directive targets the `.py` under `scripts/` only; no new file. (5) Paper round-trip — the F2 label fix (`Stage 38`->`055`, `Stage-35`->`Stage-052`) introduces no new value and is content-keyed to the notes' explicit "Stage-052 demand" wording, so it creates no new `paper_misalignment`.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Omega_exp(alpha) = pi alpha(2 alpha e^a + pi)/((4 alpha^2+pi^2)(e^a-1))` | py:34–37 / wl:38–41; sympy out L5, math out L5 | notes L47–48 (overlap branch) | MATCH |
| `A_K(y,x) = 1/(1 - x/4 + x y^2/pi^2)` | py:38 / wl:42; sympy out L6, math out L6 | notes L60 (softening branch), card eq L13–17 | MATCH |
| `zeta_family = Omega_exp^2 A_K` | py:39 / wl:43; sympy out L7, math out L7 | notes L70–72, card eq:zeta-factor L13–17 | MATCH |
| symmetric twin value `= 1` (lower endpoint) | py:45,49 / wl:49,56; sympy out L8, math out L16 | card box L21–26 (lower `1`), notes L86, L102 | MATCH |
| closure maximum `= pi^2/(4-x)` (upper endpoint) | py:46,50 / wl:50,57; sympy out L9, math out L17 | card box `(pi^2/4)(4/(4-x))` L21–26 (= pi^2/(4-x)), notes L98,L102 | MATCH |
| reachability floor `x_floor = 4 - pi^2/zeta_req` | py:54,56 / wl:51,58; sympy out L16, math out L18 | notes L122 | MATCH |
| stiffness-ratio `K_X/K_W^eff = pi^2/(4 zeta_req)` at floor | py:60 / wl:60; sympy out L18, math out L27 | notes L166 | MATCH |
| pure-overlap ceiling `pi^2/4` | py:64; sympy out L23 | notes L8, L66 (A_K(0)=4/(4-x)), L132 (Regime A) | MATCH |
| regime split A/B/C thresholds (`pi^2/4`, `pi^2/(4-x)`) | py:67–69; sympy out L25–27 | notes §3 L128–150 | MATCH |

Internal scaffolding (no prose expected): pass/fail flags (`PASS`/`FAIL`), `expect_zero`/`expectZero` residual prints (all `= 0`), `banner`/`fmt`/`pass`/`fail` helpers, the unused `KX_over_KW` symbol (py:59), `$Assumptions` domain declarations.

reconciliation: complete; 9 values checked, 0 misaligned
