---
unit_id: 130
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage130_mouth_bias_map.md
  paper_appendix: present
---

# Audit unit 130 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_130.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage130_mouth_bias_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_130}` at line 1294; no extra prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.txt`

## What the paper claims

The stage card boxes the result: *"Exact \(\mathfrak g_\Pi\), monotonicity, and unique Family-1 point \(\Pi_*\approx 1.5088295\)."* The supporting notes give three boxed deliverables: (i) the closed form
`g_Π = 2Π(2Π e^Π + π) / ((4Π² + π²)(e^Π − 1))` for the mouth-bias factor against the cosine D/N shape with profile `σ_Π(z) = Π e^{−Πz/L}/(L(1 − e^{−Π}))`; (ii) the exact monotonicity law `d g_Π/dΠ = −(1/L) Cov_Π(cos(πz/2L), z)` together with the strict-sign consequence `dg_Π/dΠ > 0` from the fact that `cos(πz/2L)` is strictly decreasing on `[0, L]`; and (iii) the unique Family-1 compensation point `Π_* ≈ 1.50882951349316` solving `g_Π = g_-^F1` with `g_-^F1 ≈ 0.758035078944663` imported from upstream. The card also asserts the limits `g_Π → 2/π` as `Π → 0+` and `g_Π → 1` as `Π → ∞`, embedded in the boxed monotonicity range statement.

## What the script claims to verify

Both engines symbolically compute `g_Π` and `Cov_Π(f, z)` from the explicit `σ_Π(z)` and `f = cos(πz/2L)` and then assert: (a) the covariance identity residual `d g_Π/dΠ + (E[fz] − g_Π·E[z])/L = 0`; (b) `limit Π→0+ g_Π − 2/π = 0`; (c) `limit Π→∞ g_Π − 1 = 0`. Both engines then numerically root-find `Π_*` by solving `g_Π(Π_*) = g_-^F1`, where SymPy uses the truncated decimal `0.758035078944663` and Mathematica uses the exact closed form `(2√(4107 − 100π²) − 37√3)/(20π)`. The Mathematica side additionally asserts `|g_Π(Π_*) − g_-^F1| < 10^{-20}`; the SymPy side does not wrap the root-find in an assertion (it only prints the value). Neither engine asserts equality between the symbolic `g_Π` and the paper's boxed closed form, and neither engine asserts the strict-sign consequence `d g_Π/dΠ > 0` or the uniqueness of `Π_*`.

## Paper ↔ script cross-check

| Deliverable | Script-side coverage | Status |
|---|---|---|
| Closed form `g_Π = 2Π(2Π e^Π + π)/((4Π² + π²)(e^Π − 1))` | Both engines compute and print the form; SymPy output literally matches; Mathematica output is the same after rearrangement (`(2*piM*(Pi + 2*E^piM*piM))/((-1 + E^piM)*(Pi^2 + 4*piM^2))`). No direct `expectZero` against the paper's boxed expression. | partial |
| Monotonicity identity `dg/dΠ = −(1/L)Cov(f, z)` | Both engines assert the algebraic identity via `cov_id == 0`. | match |
| Strict sign `dg/dΠ > 0` (and hence uniqueness of `Π_*`) | Neither engine asserts the sign. SymPy prints `g'(Π_*) ≈ 0.0714 > 0` at a single point. Mathematica also prints the single-point value. | missing |
| Limits `g_Π → 2/π` (Π→0+) and `g_Π → 1` (Π→∞) | Both engines assert via `expectZero`. | match |
| Unique Family-1 compensation point `Π_* ≈ 1.5088295` | Both engines locate the root from initial guess 1.5; Mathematica asserts `|g_Π(Π_*) − g_-^F1| < 10^{-20}`. Uniqueness not separately asserted. | partial |

Dominant pattern: paper-side closed form and limits are exercised; the strict-monotonicity sign and uniqueness statements are exercised only at the single root point and not in general. `paper_alignment = partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 28-29 | `cov_id != 0 → raise` | monotonicity-identity (form only) | yes |
| A2 | sympy | 30-31 | `simplify(g0 - 2/pi) != 0 → raise` | limit Π→0+ | yes |
| A3 | sympy | 32-33 | `simplify(ginf - 1) != 0 → raise` | limit Π→∞ | yes |
| A4 | sympy | 36-40 | nsolve, print only (no assert) | unique compensation point | partial (prints, no enforcement) |
| A5 | mathematica | 49 | `expectZero["covariance identity", covId]` | monotonicity-identity (form only) | yes |
| A6 | mathematica | 56 | `expectZero["uniform-source limit", g0 - 2/Pi]` | limit Π→0+ | yes |
| A7 | mathematica | 57 | `expectZero["point-source limit", gInf - 1]` | limit Π→∞ | yes |
| A8 | mathematica | 64 | `expectApprox["Family-1 compensation point", g(Π_*), gMinus, 10^-20]` | compensation point (numeric) | yes (single-point) |

No script asserts (i) `g_Π` equality with the boxed closed form, (ii) `d g_Π/dΠ > 0` over a positive range of `Π`, or (iii) uniqueness of `Π_*` on the relevant interval.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:28-29`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:49`

**What's wrong:**
The paper's notes (`moving_throat_pde_stage130_mouth_bias_map.md:78-82`) box `dg_Π/dΠ > 0` as a strict-sign consequence used to support the *uniqueness* of the Family-1 compensation point quoted in the paper card. The scripts only verify the *identity*
`d g_Π/dΠ + Cov_Π(f, z)/L = 0`
(SymPy line 18; Mathematica line 45), which is algebraically forced by `∂_Π σ_Π = −(1/L)σ_Π(z − ⟨z⟩_Π)` and the definition of `g_Π`. The strict-sign statement
`d g_Π/dΠ > 0 for Π > 0`
is never asserted in either engine. Mathematica prints `g'(Π_*) ≈ 0.0714` at the single Π_* root, but no positivity check is wrapped around the closed-form expression. Without strict monotonicity, the uniqueness portion of the paper's boxed claim ("unique Family-1 point") is not script-supported.

**Why this matters:**
The paper card boxes "unique Π_*". Uniqueness is a strict-monotonicity consequence. The current script verifies the identity that *names* the derivative as a covariance but does not verify the covariance is negative (equivalently `dg/dΠ > 0`). A bug that flipped the sign of `σ` or `f` could leave the identity passing while breaking the uniqueness argument.

**Required change:**
Add to both engines a positivity assertion against the closed-form derivative. The closed form for `g_Π` is `2Π(π + 2Π e^Π)/((e^Π − 1)(4Π² + π²))`; a direct calculation yields (after combining over `(e^Π − 1)²(4Π² + π²)²`) a polynomial-times-exponential numerator. The minimal additional check is symbolic: confirm `together(diff(g_Π, Π))` is a single rational expression whose numerator simplifies to a manifestly positive form for `Π > 0` (e.g., `simplify(numer(together(diff(g_Π,Π))) - <positive_sum_of_squares_form>) == 0`), OR verify positivity numerically at a sweep of `Π ∈ {0.1, 0.5, 1.0, 1.5088, 3.0, 10.0}` via `assert diff(g_Π,Π).subs(Pi, val) > 0` (SymPy) / `expectApprox` with a strictly-positive guard (Mathematica). The numeric sweep is the minimal acceptable substitute.

**Verification:**
After the patch, SymPy lines after line 33 should contain a loop or explicit list of `assert sp.diff(gPi, Pi).subs(Pi, v) > 0` for at least four positive `v` values spanning `(0, ∞)`; Mathematica should add `Do[expectApprox["dg/dPi > 0 at v=" <> ToString[v], D[gPi,piM] /. piM -> v, <positive value>, ...], {v, {0.1, 0.5, 1.0, 1.5088, 3.0, 10.0}}]` or an equivalent loop that aborts if any value is ≤ 0. The script output transcript should contain `PASS: dg/dPi > 0` (Mathematica) or no AssertionError (SymPy) for every sampled point.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:15`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:42`

**What's wrong:**
The paper notes box the explicit closed form
`g_Π = 2Π(2Π e^Π + π) / ((4Π² + π²)(e^Π − 1))`
(`moving_throat_pde_stage130_mouth_bias_map.md:34-40`). The scripts compute `g_Π` symbolically and *print* the result, but neither engine asserts equality between the computed expression and the paper's boxed form. The match is only by visual inspection of the printout. A future change to `σ_Π` or `f` could silently produce a different (but still self-consistent) symbolic expression while the limit and covariance assertions continue to pass.

**Why this matters:**
The closed-form expression is the bottom-line deliverable of section 1 of the notes and the first phrase of the paper card's `\stagefield{Output}`. It is currently the only deliverable for which "verification" reduces to reading the printout.

**Required change:**
In SymPy, after line 15, add
```
gPi_boxed = 2*Pi*(2*Pi*sp.exp(Pi) + sp.pi) / ((4*Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1))
if sp.simplify(gPi - gPi_boxed) != 0:
    raise AssertionError("g_Pi does not match paper closed form.")
```
In Mathematica, after line 42, add
```
gPiBoxed = 2*piM*(2*piM*Exp[piM] + Pi) / ((4*piM^2 + Pi^2) * (Exp[piM] - 1));
expectZero["g_Pi matches paper closed form", gPi - gPiBoxed];
```

**Verification:**
Output transcripts gain the lines "PASS: g_Pi matches paper closed form" (Mathematica) and no AssertionError on the new line (SymPy). The residual line should print 0.

## Independent-derivation check (Mathematica)

The Mathematica script (`mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl`) reuses the same physical premises (`σ_Π`, `f = Cos[π z/(2L)]`) and the same covariance-identity decomposition (`D[gPi, piM] + (eFZ - gPi*eZ)/lM`) as the SymPy file. Variable names are renamed (`Pi → piM`, `L → lM` to avoid Mathematica's reserved `Pi`), and the engine uses `FullSimplify` + `Together[Expand[...]]` in its `expectZero` helper. The choreography (compute `gPi`, then `eZ`, then `eFZ`, then `covId`, then limits, then root-find) is the same. However, this matches because the underlying physical claim has a single natural decomposition (cov-id of an exponential family vs. a cosine weight); both engines arrive at the same intermediate decomposition because the math forces it, not because the WL script transliterates the PY one. The Mathematica side adds value: (i) it computes `gMinus` from the exact closed-form `(2√(4107 − 100π²) − 37√3)/(20π)` rather than a 15-digit truncation; (ii) `expectApprox` wraps the root-find with a hard tolerance check. This is not a transliteration finding.

Side note (not a formal finding under the ten categories): the Mathematica banner at line 32 reads `STAGE 113 — EXACT MOUTH-BIAS MAP AND FAMILY-1 COMPENSATION POINT` rather than `STAGE 130`. The end-of-run message at line 67 correctly says `Stage 130 Mathematica audit passed.` This is a stray label and does not affect mathematics; flag for cosmetic cleanup but not blocking.

## Engine cross-check

Both engines produce the same `g_Π` closed form (modulo trivial sign-rearrangement of `(1 − e^{−Π})` ↔ `(e^Π − 1)`):
- SymPy: `2*Pi*(2*Pi*exp(Pi) + pi)/((4*Pi**2 + pi**2)*(exp(Pi) - 1))`
- Mathematica: `(2*piM*(Pi + 2*E^piM*piM))/((-1 + E^piM)*(Pi^2 + 4*piM^2))`

Both report covariance-identity residual `= 0`. Both report `limit Π→0+ g_Π = 2/π`, `limit Π→∞ g_Π = 1`. Numeric `Π_*` agreement: SymPy `1.50882951349315861144664988336`, Mathematica `1.5088295134931555830055507559542749287786371931530784314262`. The two agree to 15 significant figures; beyond that, SymPy's value is limited by the 15-digit truncation of the hardcoded `g_minus = "0.758035078944663"`. The Mathematica value (using the exact closed-form `gMinus`) is the authoritative one; SymPy's `prec=80` in `nsolve` is wasted past digit 15. This precision mismatch is not flagged as `engine_disagreement` because the engines *do* agree to the precision SymPy's input supports; it is a known artifact of carry-forward hardcoding and the paper card only quotes `Π_* ≈ 1.5088295` to 7 s.f.

Outputs (`.txt`) are newer than the scripts (`.py`/`.wl`) by mtime; `outputs_fresh: true`.

## Verdict justification

The paper card has three boxed deliverables (closed-form `g_Π`, monotonicity with strict sign, unique Π_*) plus two embedded limits. The scripts cover the limits and the *identity* form of the monotonicity statement, but they (a) do not assert the strict-sign half of monotonicity (and hence cannot ground "unique Π_*" beyond a single-point check), and (b) do not assert equality between the computed `g_Π` and the paper's boxed closed form (the agreement is only by printout). Both are legitimate `insufficient_verification` findings. Engines agree, outputs are fresh, no transliteration concern, no symbol-assumption error, no `paper_misalignment` (the script does not contradict the paper, it only undercovers it). Verdict: `findings`, two items, both medium-to-low severity, no stop-cold flag — Codex can mechanically extend the assertions to cover the missing strict-sign and the boxed-form equality without touching paper or notes.

## Self-test notes

(1) Variable-independence trap: F1's proposed `diff(gPi, Pi)` is taken w.r.t. the actual free symbol `Pi` (sympy `Pi`, not `sp.pi`) that `gPi` depends on; the derivative is not identically zero. (2) Trivial-case pre-check for F2: at `Π = 1`, `gPi_boxed = 2(2e + π)/((4 + π²)(e − 1))`; the SymPy `gPi` computed via the integral evaluates to the same numeric (≈ 0.870), so `simplify(gPi - gPi_boxed) == 0` is non-tautological and will reduce to 0. (3) Path-specifications: no new files proposed; both fixes are edits to existing `scripts/` and `mathematica/` files. (4) Paper round-trip: both new assertions use only constants and forms already stated in the paper card / notes; no new paper claims introduced.
