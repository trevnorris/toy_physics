---
unit_id: 050
batch: III.2
auditor_model: claude-opus-4-7
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 050 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.txt`

## What the script claims to verify

The scripts (per the SymPy docstring) claim five things:
(1) equivalence `zeta_req <= 1 <=> S_req <= 2` for the lowest twin lane,
(2) the doubling theorem `S(1; eps) = 2`,
(3) same-operator twin threshold formula `x_max(n; zeta_req) = (1/((2n+1)^2 zeta_req) - 1)/(n(n+1))` such that `zeta_n^(twin) >= zeta_req iff x <= x_max`,
(4) the impossibility bound `zeta_req > 1/(2n+1)^2` makes the twin family insufficient for n>=1,
(5) the higher-harmonic enhancement ceiling `S_n^(max) = 1 + (1-eps)/((2n+1)^2 - eps)`.
Both engines run on the closed-form expression `zeta_n^(twin) = 1/((2n+1)^2 (1 + x n(n+1)))`. In SymPy this is imported from stage 049's `twin_support_ratio`; in Mathematica it is redeclared directly.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 47 | `S(1;eps) - 2 == 0` | yes |
| A2 | sympy | 51-54 | `zeta_req - 1 - (1-eps)(S_req-2)/(1+eps(S_req-2)) == 0` | yes |
| A3 | sympy | 64-67 | `x_eq - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)] == 0` (from `sp.solve`) | partial |
| A4 | sympy | 73 | `(2n+1)^2 zeta_n - 1/(1+x n(n+1)) == 0` | no (tautological) |
| A5 | sympy | 82 | `S_n^(twin)(x=0) - S_n^(max) == 0` | partial |
| A6 | mathematica | 38 | `S(1;eps) - 2 == 0` | yes |
| A7 | mathematica | 41-44 | `zeta_req - 1 - (1-eps)(S_req-2)/(1+eps(S_req-2)) == 0` | yes |
| A8 | mathematica | 51 | `(zetaN /. x -> xEq) - zetaReq == 0` | yes |
| A9 | mathematica | 52-55 | `xEq - (1/(((2n+1)^2) zetaReq) - 1)/(n(n+1)) == 0` | no (tautological by definition) |
| A10 | mathematica | 59 | `(2n+1)^2 zetaN - 1/(1+x n(n+1)) == 0` | no (tautological) |
| A11 | mathematica | 65 | `S_n^(twin)(x=0) - S_n^(max) == 0` | partial |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:71-73`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:57-59`

**What's wrong:**
The "suppression factor" check (`supp - 1/(1 + x n(n+1)) == 0`) is an algebraic tautology by construction. `zeta_n` is set to the closed form `1/((2n+1)^2 (1 + x n(n+1)))` (imported from stage049's `twin_support_ratio` in SymPy, written out directly on line 46 in Mathematica). Then `supp = (2n+1)^2 * zeta_n` simplifies to `1/(1 + x n(n+1))` by direct cancellation — no algebra to verify. The assertion cannot fail no matter what physics the underlying formula was supposed to encode.

**Why this matters:**
The docstring labels this as proof of the "impossibility bound from higher-harmonic suppression" (claim 4), but the only thing the assertion proves is that one can cancel `(2n+1)^2` against itself. The impossibility bound itself (the statement that `zeta_req > 1/(2n+1)^2` admits no x >= 0) is never asserted — only printed as text on SymPy line 74. A reviewer skimming the PASS line would believe the impossibility bound was verified when it was not.

**Required change:**
Replace the trivial cancellation check with a non-tautological assertion of the impossibility bound itself. The substantive statement is: the equation `zeta_n^(twin)(x) = zeta_req` has a solution with `x >= 0` if and only if `zeta_req <= 1/(2n+1)^2`. Concretely, assert that the numerator of `x_max` (which is `1/((2n+1)^2 zeta_req) - 1`) flips sign exactly at `zeta_req = 1/(2n+1)^2`. Add an assertion of the algebraic identity `((2n+1)^2 * zeta_req - 1) + (2n+1)^2 * zeta_req * (x_max * n * (n+1)) == 0` after expanding `x_max`, so that admissibility (`x_max >= 0`) reduces to `(2n+1)^2 * zeta_req <= 1` for any `zeta_req > 0`.

**Verification:**
A new assertion appears in both scripts that does not reduce to a `(2n+1)^2 / (2n+1)^2` cancellation. The saved output should show a non-trivial residual being simplified to 0, anchored to the inequality `zeta_req <= 1/(2n+1)^2`.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:47,52-55`

**What's wrong:**
On line 47 the Mathematica script defines
`xEq = FullSimplify[(1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1)), Assumptions -> $Assumptions];`
On lines 52-55 it then asserts
`xEq - (1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1)) == 0`.
This compares xEq to its own definition. It is a textbook `x = expr; assert x == expr` tautology and cannot fail under any value of zetaReq, n, or eps.

**Why this matters:**
The SymPy counterpart at lines 62-67 derives `x_eq` via `sp.solve(sp.Eq(zeta_n, zeta_req), x)[0]` and then compares to the closed form — which, while not deep, at least exercises that `sp.solve` returns the claimed parameterization. The Mathematica script skips that derivation entirely and ships a definitional identity instead. The Mathematica script needs an analog of the SymPy `solve` step, or a verification step that `zetaN[x_max] == zetaReq` (which DOES exist on line 51 — that one is fine). The duplicate tautological assertion on lines 52-55 should either be removed or replaced with an independent derivation.

**Required change:**
Replace the assertion on lines 52-55 with a derivation: use `Solve[zetaN == zetaReq, x]` in Mathematica, extract the solution, and assert it equals the claimed closed form `(1/(((2n+1)^2) zetaReq) - 1)/(n (n+1))`. This makes the Mathematica check an independent algebraic re-derivation rather than a definitional restatement.

**Verification:**
The Mathematica script now contains a `Solve[...]` step for `x_max`. The saved output shows the solve returning the expected formula, with a residual of 0 against the closed form on a separate line. The tautological line 52-55 form is gone.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:76-85`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:61-67`

**What's wrong:**
The "higher-harmonic enhancement ceiling" check (claim 5) only verifies `S_n^(twin)(x=0) == S_n^(max)`. It never verifies that this value is actually a *ceiling* — i.e., that `S_n^(twin)(x) <= S_n^(max)` for all admissible `x >= 0`. The check at the single point x=0 is consistent with either monotonic decrease, monotonic increase, or non-monotonic behavior in x. The name "S_n^(max)" and the docstring claim "enhancement ceiling" therefore go unverified.

**Why this matters:**
A reader sees `S_n^(twin)(x=0) - S_n^(max) = 0 ... PASS` and concludes that S_n_max is established as an upper bound. The assertion does not establish that — it establishes only that the two expressions coincide at one point of the domain. Without a monotonicity (or sign-of-derivative) check, the "ceiling" claim is unverified.

**Required change:**
Add an assertion that `S_n^(twin)(x) - S_n^(max)` has a definite sign for `x > 0` under the script's stated assumptions `0 < eps < 1, n >= 1`. The cleanest form: compute `diff(S_n^(twin), x)` (sympy) / `D[sN, x]` (mathematica), `FullSimplify`/`simplify` it, and assert the simplified result, when written as `numerator/denominator^2`, has a numerator whose sign is fixed (negative) under the assumptions — i.e., assert `numerator(diff(S_n,x)) + |something_nonneg| == 0` so SymPy/Mathematica reduce it to 0. Alternatively, factor `S_n^(max) - S_n^(twin)` and assert it equals `(positive_factor) * x * n * (n+1)` so the sign is manifest for x>=0.

**Verification:**
A new assertion appears that exercises the sign of `S_n^(max) - S_n^(twin)(x)` for `x > 0`, not just the value at x=0. The saved output shows the derivative or factored form being reduced to a manifestly-signed expression.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:56-67`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:46-55`

**What's wrong:**
Claim 3 in the SymPy docstring is "`zeta_n^(twin) >= zeta_req gives x <= x_max`" — a directional inequality. The actual assertion only solves for *equality* (`zeta_n = zeta_req`) and pattern-matches the resulting x against the closed form. The direction of the inequality (that `zeta_n^(twin)` is monotonically decreasing in `x`, so that the threshold has the claimed sign) is never asserted. The Mathematica analog checks `(zetaN /. x -> xEq) - zetaReq == 0` (line 51), which is again equality only.

**Why this matters:**
Without a monotonicity check, swapping `x <= x_max` for `x >= x_max` in the docstring would not be caught by the assertions. The script's stated claim is stronger than what is verified.

**Required change:**
Add a one-line assertion that `d zeta_n^(twin)/dx` has fixed (negative) sign for `x >= 0, n >= 1` — for instance, factor the derivative and assert it equals `-n(n+1)/((2n+1)^2 (1 + x n(n+1))^2)`, which is manifestly negative for the stated assumption domain. This pins down the direction of the threshold.

**Verification:**
A new assertion appears that exercises the derivative of `zeta_n^(twin)` with respect to `x` (or the sign of `zeta_n^(twin)(x1) - zeta_n^(twin)(x2)` for an ordered pair). The saved output shows the residual being reduced to 0 against a manifestly-negative reference expression.

## Independent-derivation check (Mathematica)

The Mathematica script is *not* a line-by-line transliteration of the SymPy script. It diverges in several substantive places: (a) it redeclares `zetaN` directly rather than importing it from stage049, (b) it defines `xEq` as the closed form directly and adds the substantive `zetaN[xEq] - zetaReq == 0` check (line 51), which the SymPy script lacks, (c) it skips the SymPy `sp.solve` step. Variable choreography is similar (zetaReq, sEnhance, zetaN, sN, sNMax) but that similarity is forced by the underlying claim's algebraic structure. I do not file `mathematica_transliteration` here — the Mathematica derivation has independent structural choices, even though the algebraic skeleton is shared. The duplicate tautological check (lines 52-55) is captured under F2 as a tautology, not a transliteration.

## Engine cross-check

Both engines produce the same final closed forms:
- `zeta_req = (S_req - 1)/(1 + eps(S_req - 2))` (SymPy line 17 of output; Mathematica line 13)
- `S(zeta;eps) = 1 + zeta(1-eps)/(1 - eps zeta)` (algebraic equivalents)
- `zeta_n^(twin) = 1/((2n+1)^2 (1 + x n(n+1)))`
- `(2n+1)^2 zeta_n^(twin) = 1/(1 + x n(n+1))`
- `S_n^(twin)(x=0) = 1 + (1-eps)/((2n+1)^2 - eps)` = `S_n^(max)`
- `S_1^(max) = 2(eps-5)/(eps-9)` (SymPy) = `1 + (-1+eps)/(-9+eps)` (Mathematica) — algebraically identical
- `S_2^(max) = 2(eps-13)/(eps-25)` (SymPy) = `1 + (-1+eps)/(-25+eps)` (Mathematica) — algebraically identical

Engines agree. No `engine_disagreement` finding.

## Verdict justification

The scripts establish the algebraic identities they assert. The doubling theorem `S(1;eps) = 2` and the `zeta_req - 1` factorization are substantive non-tautological checks that both engines exercise. However, two of the five docstring claims (the impossibility bound, claim 4; and the directional threshold + ceiling claims, parts of claims 3 and 5) are only printed or only checked at single points — they are not exercised as asserted. In addition, both scripts contain a `(2n+1)^2 / (2n+1)^2` cancellation that is labeled as the proof of claim 4, and the Mathematica script contains a definitional tautology (xEq vs its own definition). Verdict: `findings`. Output files are fresh (script mtime 11:56; outputs 12:43-12:51), engines agree, both scripts present.

## Self-test notes

Checked: (1) variable independence — F3 and F4 propose derivatives `d/dx` of expressions that explicitly contain `x`, so the derivatives are non-zero; (2) parity/symmetry — not applicable, the domain is `x >= 0` not symmetric; (3) trivial-case pre-check — at n=1, eps=1/2, zeta_req=1/9: x_max = 0, and the impossibility bound `(2n+1)^2 zeta_req = 1` is saturated as required, consistent with F1's proposed identity; at the same point S_n^(max) = 1 + (1/2)/(9-1/2) = 1 + 1/17, matching the script's formula. (4) Path specifications — both scripts already exist at the expected paths; no missing-script directive required.
