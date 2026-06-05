---
unit_id: 075
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T21:20:26Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage075_family1_threshold_window.md]
  paper_appendix: present
---

# Audit unit 075 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_075.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 128)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.txt`

## What the paper claims

Stage 075 evaluates the exact operator-threshold functions on the explicit Family-1 / healing-locked branch with inputs `kappa = 12321/5`, `eta = 37`, `Lambda_ell = 37`, and `Upsilon_w = 100 Theta_w` (i.e. `alpha_r = 10`, `Upsilon_w = alpha_r^2 Theta_w`). The card states the endpoint kernel scales `Delta_0(12321/5,37) ≈ 1.73302079021525e-4` and `Delta_inf(12321/5,37) ≈ 2.01447565540522e-2`. The `\stagefield{Output}` is the boxed explicit wall-depth window: `Theta_w <= 3.62605617972939e-4 Pe_req ⇒ fail` and `Theta_w >= 4.21495341569977e-2 Pe_req ⇒ succeed`. The notes add the intermediate `Upsilon` thresholds (`Upsilon_fail ≈ 0.0362605617972939 Pe_req`, `Upsilon_suff ≈ 4.21495341569977 Pe_req`), the equivalent `Xi` thresholds (`Xi_fail ≈ 49.6407091004953 Pe_req`, `Xi_suff ≈ 5770.27122609299 Pe_req`), the large-`alpha` interpretation (`alpha = sqrt(kappa) = 111/sqrt(5) ≈ 49.6407091`, `Delta_inf ~ 1/alpha`, `Xi_fail ≈ alpha`), and the wall-amplitude reduction `Upsilon_w = alpha_r^2 Theta_w` with `alpha_r = 10`. The appendix row (128) summarizes it as "Explicit `Theta_w` fail/succeed window."

## What the script claims to verify

Both engines (a) pin `alpha_r^2 = 100` as a lock on the paper Inputs line; (b) build the closed forms for `Delta_0` and `Delta_inf` from `(kappa, eta)`, confirm they satisfy their defining algebraic identities (explicitly flagged in-script as tautological), and verify two non-trivial asymptotic limits (`alpha*Delta_inf -> 1` as `alpha -> oo`, `Delta_0 -> 1/2` as `alpha -> 0`); (c) form `Upsilon_fail/suff`, `Xi_fail/suff`, `Theta_fail/suff` and print their closed forms and 20-digit numerics; and (d) assert a `Upsilon = alpha_r^2 * Theta` round-trip. The Mathematica script additionally adds eight independent `expectApprox` numeric anchors that compare the computed `Delta_0`, `Delta_inf`, `Upsilon`, `Xi`, and `Theta` numerics against hardcoded literals matching the paper/notes values (the SymPy script only *prints* these numerics, it does not assert them).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Delta_0 ≈ 1.73302079021525e-4` | closed form + `expectApprox` (math:119); sympy prints `0.00017330207902152514906` | match |
| `Delta_inf ≈ 2.01447565540522e-2` | closed form + `expectApprox` (math:120); sympy prints `0.020144756554052159427` | match |
| `Upsilon_w = 100 Theta_w` (`alpha_r=10`) | `assert alpha_r**2 == 100` (sympy:28), `expectZero["alpha_r^2-100..."]` (math:41) | match |
| `Theta_fail = 3.62605617972939e-4 Pe_req` (fail box) | `Theta_fail` formed + `expectApprox` (math:125); sympy only prints (no assert) | match (sympy: partial — printed not asserted) |
| `Theta_suff = 4.21495341569977e-2 Pe_req` (succeed box) | `Theta_suff` formed + `expectApprox` (math:126); sympy only prints (no assert) | match (sympy: partial — printed not asserted) |
| `Upsilon_fail/suff` (notes) | formed + `expectApprox` (math:121-122) | match |
| `Xi_fail/suff` (notes) | formed + `expectApprox` (math:123-124) | match |
| large-`alpha`: `Xi_fail ≈ alpha`, `Delta_inf ~ 1/alpha` (notes) | `alpha*Delta_inf -> 1` limit (sympy:74-76, math:104-108) | match |

Every paper/notes value reconciles with the scripts (see Value Reconciliation section). `paper_alignment: aligned`. The only weakness is on the SymPy side, where the two boxed `\stagefield{Output}` endpoints (`Theta_fail`, `Theta_suff`) are printed but not asserted, and the one assertion the SymPy script *does* aim at the reduction is tautological (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 28 | `assert alpha_r**2 == 100` | `Upsilon_w = 100 Theta_w` (alpha_r=10) | yes (locks the only real content of the reduction) |
| A2 | sympy | 65-66 | `assert delta0/deltainf_identity == 0` | Delta closed forms | no (self-flagged tautological by construction) |
| A3 | sympy | 76 | `assert large_alpha_check == 1` | `Delta_inf ~ 1/alpha`, `Xi_fail ≈ alpha` | yes (non-trivial asymptotic) |
| A4 | sympy | 81 | `assert small_alpha_check_delta0 == 1/2` | Delta_0 closed-form shape | yes (non-trivial asymptotic) |
| A5 | sympy | 124-125 | `assert simplify(Upsilon_fail - alpha_r^2*Theta_fail) == 0` | `Upsilon_w = alpha_r^2 Theta_w` | **no — tautological (Theta_fail ≡ Upsilon_fail/100 by line 99)** |
| A6 | math | 41 | `expectZero[alphaR^2 - 100]` | reduction constant | yes |
| A7 | math | 88-93 | `expectZero[Delta_0/inf algebraic identity]` | Delta closed forms | no (self-flagged tautological) |
| A8 | math | 104-113 | `Limit` checks (`alpha*Delta_inf->1`, `Delta_0->1/2`) | asymptotics | yes |
| A9 | math | 116-117 | `expectZero[Upsilon - alphaR^2*Theta]` | reduction | **no — tautological (thetaFail ≡ upsilonFail/100 by line 69)** |
| A10 | math | 119-126 | 8× `expectApprox` vs paper literals | all Delta/Upsilon/Xi/Theta numerics incl. both boxed Theta endpoints | yes (independent numeric pins) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:117-125`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:79-117` (lines 116-117 are the assertion; 79-83 is the misleading comment context)

**What's wrong:**
In SymPy, `Theta_fail` and `Theta_suff` are *defined* as `Upsilon_fail/alpha_r**2` and `Upsilon_suff/alpha_r**2`:

```python
Theta_fail = sp.simplify(Upsilon_fail / alpha_r**2)   # line 99
Theta_suff = sp.simplify(Upsilon_suff / alpha_r**2)   # line 100
```

The "round-trip" assertions then multiply that definition straight back by `alpha_r**2`:

```python
Upsilon_fail_from_Theta = sp.simplify(alpha_r**2 * Theta_fail)            # line 120 == Upsilon_fail
assert sp.simplify(Upsilon_fail - Upsilon_fail_from_Theta) == 0          # line 124  (Upsilon_fail - Upsilon_fail)
```

This is `100 * (Upsilon_fail/100) - Upsilon_fail == 0`, i.e. `Upsilon_fail - Upsilon_fail == 0`, which is algebraically guaranteed for *any* `Upsilon_fail` regardless of the physics. The inline comment at lines 117-119 explicitly asserts the opposite — *"Test the round-trip on the actually-constructed fail and suff branches, **not the trivial identity 100\*Theta == 100\*Theta**"* — but it IS exactly that trivial identity, merely routed through the intermediate `Theta_fail = Upsilon_fail/100`. The Mathematica script reproduces the same pattern: `thetaFail = upsilonFail/alphaR^2` (line 69) then `expectZero["Upsilon_fail - alphaR^2 * Theta_fail", upsilonFail - alphaR^2*thetaFail]` (line 116) = `upsilonFail - upsilonFail`.

The genuine content of the reduction `Upsilon_w = alpha_r^2 Theta_w` is the single constant `alpha_r^2 = 100`, and that is already locked non-tautologically by the separate `assert alpha_r**2 == 100` (sympy:28) / `expectZero[alphaR^2 - 100]` (math:41). The round-trip therefore verifies nothing beyond those locks while carrying a comment that falsely advertises it as a non-trivial branch check.

**Why this matters:**
A reviewer reading the comment would believe the boxed `Theta_fail`/`Theta_suff` window endpoints are independently cross-checked against the `Upsilon` thresholds, when in fact they cannot fail. On the SymPy side this is the *only* assertion aimed at the reduction/`Theta` endpoints — the two boxed `\stagefield{Output}` values (`Theta_fail = 3.626e-4 Pe_req`, `Theta_suff = 4.215e-2 Pe_req`) are otherwise only `print`-ed (sympy:114-115, 132-133), never asserted. So the paper's bottom-line deliverable has no genuine SymPy-side anchor; the Mathematica `expectApprox` block (lines 125-126) is currently the sole assertion pinning those endpoints to the paper values.

**Required change:**
Replace the tautological round-trip in the SymPy script with a genuine numeric anchor of the boxed window endpoints to the paper-stated literals (mirroring the Mathematica `expectApprox` block, which already does this independently), and fix the misleading comment. Concretely (SymPy lines 117-125):
- Remove the false "not the trivial identity" comment and the two round-trip asserts.
- Add asserts that `Theta_fail/Pe_req` and `Theta_suff/Pe_req` (and, since they are cheap and already computed, `Upsilon_fail/suff`, `Delta_0`, `Delta_inf`) numerically match the paper/notes literals to a tight tolerance, e.g. `assert abs(sp.N(Theta_fail/Pe_req, 30) - sp.Rational(...)) < 1e-18` using the paper literals `3.62605617972939e-4` and `4.21495341569977e-2`. This is non-tautological: the literal is independent of the computed closed form, so a wrong `Delta`/`Upsilon`/factor would make it fail.
- Keep the existing `assert alpha_r**2 == 100` (sympy:28) as the reduction-constant lock — it already carries the genuine content of `Upsilon_w = alpha_r^2 Theta_w`.
- The Mathematica round-trip `expectZero` at lines 116-117 is likewise tautological and should be dropped (or replaced with the same independent literal anchor it effectively already has at lines 125-126); at minimum the misleading "genuine independent check" framing in the comment block (lines 79-83 / 96-99) should not be applied to the round-trip.

**Verification:**
After the fix, the SymPy output must show new PASS lines for `Theta_fail`/`Theta_suff` (and the other `Delta`/`Upsilon` numerics) checked against the paper literals, and the round-trip `Upsilon - alpha_r^2*Theta == 0` lines (sympy output lines 35-36) should no longer be presented as the reduction's verification. The script must still exit 0. The verifier confirms the new numeric asserts reference `3.62605617972939e-4` / `4.21495341569977e-2` (the paper's boxed values) and that removing the tautological round-trip did not drop coverage of the reduction (the `alpha_r**2 == 100` lock remains).

## Independent-derivation check (Mathematica)

The `.wl` is structurally parallel to the `.py` — same `Delta_0`/`Delta_inf` closed forms, same algebraic-identity check, same two asymptotic limits, same `Upsilon`/`Xi`/`Theta` construction. This parallelism is **acceptable and expected here**, not a `mathematica_transliteration` finding, for two reasons: (1) `Delta_0` and `Delta_inf` are *upstream-defined results* (Stage 058/059/072 machinery), so both engines legitimately must start from the same physical closed forms — neither is meant to re-derive them from scratch in this stage; the genuine cross-check is whether each engine *independently* simplifies and evaluates those forms to the same numerics, which they do (SymPy via `sp.simplify`/`sp.limit`/`sp.N`, Mathematica via `FullSimplify`/`Limit`/`N`). (2) The Mathematica script is not a mere echo: it adds an eight-line independent `expectApprox` numeric-anchor block (lines 119-126) tying every deliverable to a hardcoded literal, which the SymPy script lacks as assertions. The two engines' simplified closed forms even land in different surface forms (compare sympy output lines 15-16 vs math output lines 11-12), confirming independent simplification rather than line-by-line porting. No transliteration finding.

## Engine cross-check

Both engines agree to full precision on every deliverable:

| value | SymPy output | Mathematica output |
|---|---|---|
| Delta_0 | 0.00017330207902152514906 | 0.00017330207902152514905715619654992403 |
| Delta_inf | 0.020144756554052159427 | 0.02014475655405215942710329560991777563 |
| Upsilon_fail/Pe | 0.036260561797293886969 | matches (diff 0, math:37) |
| Upsilon_suff/Pe | 4.2149534156997728721 | matches (diff 0, math:39) |
| Xi_fail/Pe | 49.640709100495331260 | matches (diff 0, math:41) |
| Xi_suff/Pe | 5770.2712260929890619 | matches (diff 0, math:43) |
| Theta_fail/Pe | 0.00036260561797293886969 | matches (diff 0, math:45) |
| Theta_suff/Pe | 0.042149534156997728721 | matches (diff 0, math:47) |

All `expectApprox` diffs in the Mathematica output are `0` to displayed precision; the asymptotic limits agree (`1` and `1/2`). Engines agree.

## Verdict justification

The stage holds up against the paper: every constant in the card and notes (`Delta_0`, `Delta_inf`, the `Upsilon`/`Xi`/`Theta` thresholds, the `alpha_r^2 = 100` reduction) is reproduced by both engines to full precision, the two non-trivial asymptotic limits pass, and the engines agree. The single defect is F1: the SymPy and Mathematica "round-trip" assertions `Upsilon - alpha_r^2*Theta == 0` are tautological (Theta is *defined* as Upsilon/100), and the SymPy comment falsely advertises them as non-trivial branch checks. Because the genuine content (`alpha_r^2 = 100`) is locked elsewhere and the Mathematica `expectApprox` block independently pins the boxed endpoints to the paper values, the defect is low-severity script-hygiene, not a math error — hence `verdict: findings`, low, no `stop_cold`. Attacks that failed: I checked the small-/large-alpha limits against the closed forms by hand (`Delta_0 -> 1/2` and `alpha*Delta_inf -> 1` are both correct and genuinely sensitive to a wrong factor), confirmed `alpha = sqrt(12321/5) = 111/sqrt(5)`, confirmed every paper literal matches the script literals digit-for-digit, and confirmed the two engines' simplified closed forms differ in surface form (ruling out transliteration). I read the paper card, the single notes file, and appendix row 128; the script's verified claim matches the paper's claim.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit against the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Delta_0 ≈ 1.73302079021525e-4` | py out L17, wl out L13 | tex:17 `1.73302079021525e-4`; md:39 `1.73302079021525e-4` | MATCH |
| `Delta_inf ≈ 2.01447565540522e-2` | py out L18, wl out L13 | tex:22 `2.01447565540522e-2`; md:41 `2.01447565540522e-2` | MATCH |
| `Upsilon_fail ≈ 0.0362605617972939 Pe_req` | py out L28 | md:45 `0.0362605617972939 * Pe_req` (not in tex — intermediate) | MATCH |
| `Upsilon_suff ≈ 4.21495341569977 Pe_req` | py out L29 | md:47 `4.21495341569977 * Pe_req` | MATCH |
| `Xi_fail ≈ 49.6407091004953 Pe_req` | py out L30 | md:51 `49.6407091004953 * Pe_req` | MATCH |
| `Xi_suff ≈ 5770.27122609299 Pe_req` | py out L31 | md:53 `5770.27122609299 * Pe_req` | MATCH |
| `Theta_fail ≈ 3.62605617972939e-4 Pe_req` (boxed fail) | py out L32, wl out L45 | tex:28 `3.62605617972939e-4 Pe_req`; md:126 `3.62605617972939e-4 * Pe_req` | MATCH |
| `Theta_suff ≈ 4.21495341569977e-2 Pe_req` (boxed succeed) | py out L33, wl out L47 | tex:34 `4.21495341569977e-2 Pe_req`; md:128 `4.21495341569977e-2 * Pe_req` | MATCH |
| `alpha_r^2 = 100` (reduction `Upsilon_w = 100 Theta_w`) | py L28, wl L41 | tex:7 `Upsilon_w=100 Theta_w`, tex:24 `alpha_r=10`; md:108 `Upsilon_w = 100 Theta_w` | MATCH |
| `alpha = sqrt(kappa) = 111/sqrt(5)` | py out L14, wl out L10 | md:63 `alpha := sqrt(kappa) = 111/sqrt(5)` | MATCH |
| `kappa = 12321/5`, `eta = 37`, `Lambda_ell = 37` | py L31-33, wl L42-44 | tex:7 (kappa,eta); md:11/31-35 (all three) | MATCH |

INTERNAL (verification scaffolding, no prose expected): `delta0_identity`/`deltainf_identity` residuals (=0), `large_alpha_check`/`small_alpha_check_delta0` limit literals (`1`, `1/2`), `Upsilon_fail_from_Theta`/`Upsilon_suff_from_Theta` round-trip residuals (the tautology in F1), all `expectApprox` `diff` values, all PASS flags, `Pe_req` symbol.

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

## Self-test notes

I checked: (1) Variable independence — no `sp.diff`/`D[]` derivatives in this unit, so the zero-derivative trap does not apply. (2) The asymptotic limits by hand: `Delta_0 -> 1/2` (since `cosh-1 ~ alpha^2/2`, denominator `~ alpha^2*eta`) and `alpha*Delta_inf -> 1` (since `Delta_inf ~ cosh/(alpha*sinh) ~ 1/alpha`) — both correct and factor-sensitive, so A3/A4/A8 are genuine. (3) Trivial-case: I confirmed the F1 round-trip reduces to `Upsilon - Upsilon == 0` identically (the tautology) and that the proposed replacement (numeric anchor of `Theta_fail/suff` against the paper literals `3.626e-4`/`4.215e-2`) is non-tautological because the literal is independent of the computed closed form. (4) Path: F1 edits existing `scripts/*.py` and `mathematica/*.wl`, no new file path. (5) Paper round-trip: the proposed fix anchors to the paper's own boxed literals and keeps the `alpha_r**2 == 100` lock, so it introduces no new `paper_misalignment`.
