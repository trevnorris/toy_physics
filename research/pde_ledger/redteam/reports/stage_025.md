---
unit_id: 025
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 8
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 025 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.txt`

## What the script claims to verify

Per the SymPy docstring, the scripts purport to verify, for the minimal isotropic single-support/single-port closure, (i) exact zero-frequency coefficients `B0 = C^2/varpi^2`, `Z0 = Q/Delta`, `N0 = P^2/Delta^2` with `Delta = OmegaU^2 OmegaW^2 - R^2`, `Q = GU^2 OmegaW^2 + 2 GU GW R + GW^2 OmegaU^2`, `P = OmegaU^2 GW + R GU`; (ii) the exact closed formula `P0 = N0/D0` with `D0 = K - B0 - Z0`, expressed in compact form `P0 = P^2/(Delta*(K*Delta - Delta*C^2/varpi^2 - Q))`; (iii) the "exact target equation" `mhat^2 * P0 = 54 G c_s^5 / (5 a^5 c^5)`; (iv) the "stability-domain positivity structure"; and (v) "exact monotonic derivatives" of `P0` with respect to `K` and to `X = C^2/varpi^2`. The assertions actually executed are: in II.1 that two algebraic rearrangements of the same `N0/D0` are equal; in III that `Delta*D0 = K*Delta - Delta*C^2/varpi^2 - Q` and `N0 = P^2/Delta^2`; in IV that `dP0/dK = -N0/(K - X - Q/Delta)^2`, `dP0/dX = +N0/(K - X - Q/Delta)^2`, and that the two sum to zero. Section II.2 prints the residual `mhat^2*P0_compact - target` but executes no assertion on it. Sections III and IV print prose about positivity and monotonicity but never test them.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 104 | `expect_zero("P0 - P0_compact", P0 - P0_compact)` | **no — tautological (both forms are simplify of the same N0/D0)** |
| A2 | sympy | 111 | `sp.pprint(equation_residual)` (no assertion) | **no — no assertion at all on the "exact target equation"** |
| A3 | sympy | 123 | `expect_zero("Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", ...)` | **no — tautological by definition of D0** |
| A4 | sympy | 124 | `expect_zero("N0 - P^2/Delta^2", ...)` | **no — tautological by definition of N0 at line 76** |
| A5 | sympy | 146 | `expect_zero("dP0/dK + N0/(K - X - Q/Delta)^2", ...)` | **no — tautological because N0 has no K-dependence and P0 is defined explicitly as N0/(K-X-Q/Delta)** |
| A6 | sympy | 147 | `expect_zero("dP0/dX - N0/(K - X - Q/Delta)^2", ...)` | **no — tautological by same construction (X dependence linear in denominator)** |
| A7 | sympy | 148 | `expect_zero("dP0/dX + dP0/dK", ...)` | **no — follows from A5+A6 by linearity** |
| M1 | mathematica | 65 | `expectZero["P0 - P0_compact", p0 - p0Compact]` | no — same content as A1 |
| M2 | mathematica | 71 | `Print["Target residual = ", ...]` (no assertion) | no — same gap as A2 |
| M3 | mathematica | 78 | `expectZero["Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", ...]` | no — same as A3 |
| M4 | mathematica | 79 | `expectZero["N0 - P^2/Delta^2", ...]` | no — same as A4 |
| M5-M7 | mathematica | 94-96 | three `expectZero` for derivative identities | no — same as A5-A7 |

Every "yes" row is absent. Every executed assertion in both engines is either an algebraic identity by construction or a derivative that is forced by how `P0` was defined two lines above. The script's single physically meaningful relation (the Section II.2 "exact target equation") is never asserted in either engine.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:97-104`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:58-65`

**What's wrong:**
Section II.1 defines `P0 = sp.simplify(N0/D0)` (sympy line 98) and `P0_compact = sp.simplify(P**2/(Delta*(K*Delta - Delta*C**2/varpi**2 - Q)))` (line 99), then asserts `P0 - P0_compact == 0` (line 104). Since `D0 = K - B0 - Z0 = K - C^2/varpi^2 - Q/Delta = (K*Delta - Delta*C^2/varpi^2 - Q)/Delta`, the rearrangement `N0/D0 = (P^2/Delta^2)*Delta/(K*Delta - Delta*C^2/varpi^2 - Q) = P^2/(Delta*(K*Delta - Delta*C^2/varpi^2 - Q))` is an algebraic identity guaranteed by construction. `simplify` is built to find exactly this kind of identity. The Mathematica file at lines 59-65 performs the identical move. No physics is checked; only that SymPy/Mathematica can each cancel a Delta.

**Why this matters:**
The script's docstring says it verifies an "exact closed formula for P0 = N0/D0". The check tests neither the *value* of P0 nor that the compact form is the right rewrite for any downstream consumer — it only tests the trivial algebraic identity between two forms the script itself just wrote down.

**Required change:**
Replace the tautological II.1 check with a substantive test: pick at least one concrete, internally consistent numerical sample point (`K=2`, `varpi=1`, `C=1`, `OmegaU=2`, `OmegaW=2`, `R=1`, `GU=1`, `GW=1` — chosen so `Delta = 15 > 0` and `D0 > 0`) and verify (a) that `P0.subs(sample)` evaluates to a single positive rational and (b) that the same numerical value is reproduced both by the `N0/D0` form and the compact `P^2/(Delta*(K*Delta - Delta*X - Q))` form. The check still passes only if the algebraic rearrangement is correct AND the formula produces a sensible positive value on the stability branch.

**Verification:**
After the fix, the saved output must show a new line `P0_numerical = <rational>` matching for both forms, and the existing `expect_zero("P0 - P0_compact", ...)` may remain as a redundant check but a new `expect_positive("P0 on sample point", ...)` line must appear.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:121-124`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:74-79`

**What's wrong:**
Section III's two assertions are
```
expect_zero("Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", compact_denom - (K*Delta - Delta*C**2/varpi**2 - Q))
expect_zero("N0 - P^2/Delta^2", N0 - P**2/Delta**2)
```
`D0` is defined at line 77 as `K - B0 - Z0 = K - C^2/varpi^2 - Q/Delta`, so `Delta*D0 = K*Delta - Delta*C^2/varpi^2 - Q` is an arithmetic restatement. `N0` is defined at line 76 as `P**2/Delta**2`, so `N0 - P**2/Delta**2 == 0` is a definition tautology. Neither tests the "stability-domain positivity structure" the section banner claims to address. The Mathematica file at lines 77-79 repeats the same two tautologies.

**Why this matters:**
The section is titled "STABILITY AND POSITIVITY STRUCTURE" and the print at sympy line 126 says "If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0" — but this implication is asserted only in prose, never tested. If the formula for P0 had a sign error (e.g., `N0 = -P^2/Delta^2`) the current Section III checks would still pass, while the actual positivity statement would be falsified.

**Required change:**
Replace (or augment) the two tautological identities with a real positivity test: with the same sample point as F1, substitute and verify that
- `Delta.subs(sample) > 0`
- `D0.subs(sample) > 0`
- `P0.subs(sample) > 0`
Use an `assert simplified.subs(sample) > 0` or equivalent (sympy `> 0` returns `BooleanTrue`/`BooleanFalse` on numeric inputs). In Mathematica, use `If[NumericQ[val] && val > 0, pass, fail]`. This converts the prose claim into an executed check.

**Verification:**
After the fix, saved output must show three new lines `Delta on sample = <positive rational>`, `D0 on sample = <positive rational>`, `P0 on sample = <positive rational>`, each followed by a PASS marker.

### F3 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:133-148`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:83-96`

**What's wrong:**
In Section IV `P0` is redefined as `sp.simplify(N0 / (K - X - Q/Delta))` (line 137). Both `N0` and `Q/Delta` depend on `GU, GW, R, OmegaU, OmegaW` but not on `K` or on the new symbol `X`. By the elementary quotient rule, `dP0/dK = -N0/(K-X-Q/Delta)^2` and `dP0/dX = +N0/(K-X-Q/Delta)^2` are forced by the way `P0` was just written. The three assertions at lines 146-148 are guaranteed by `sp.diff` correctly differentiating a literal expression of the form `A/(K-X-B)`; they tell us nothing about the physics. The Mathematica file at lines 86-96 performs the same syntactic move.

Worse: the section banner says "EXACT MONOTONIC DERIVATIVES" but the checks never test the *sign* of those derivatives (which is the actual monotonicity content). The sign of `dP0/dK = -N0/(K-X-Q/Delta)^2` is `-sign(N0)` (denominator squared is non-negative). Since the script's assumptions allow `Delta < 0` (R is unrestricted real and OmegaU^2*OmegaW^2 - R^2 can be either sign), `N0 = P^2/Delta^2` is non-negative but `Q/Delta` can flip sign with Delta. None of this is exercised.

**Why this matters:**
A reader of the docstring would believe "P0 decreases in K, increases in X" is a verified claim. It is not. If a downstream unit relies on `dP0/dK < 0` for an existence argument or a stability condition, that downstream argument is unsupported by this script.

**Required change:**
Replace (or augment) the three tautological derivative identities with three numerical-sign checks at the sample point from F1:
- `assert dP0_dK.subs(sample) < 0`
- `assert dP0_dX.subs(sample) > 0`
- `assert (dP0_dK + dP0_dX).subs(sample) == 0` (this is the only structural one worth keeping; it tests that the only K/X dependence enters via the combination K - X)
In Mathematica, mirror via `If[val < 0, pass, fail]` style.

**Verification:**
After the fix, saved output must show `dP0/dK on sample = <negative rational>`, `dP0/dX on sample = <positive rational>`, and the sum-to-zero check.

### F4 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:106`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:67`

**What's wrong:**
The "exact target equation" right-hand side `target = 54*G*c_s**5/(5*a**5*c**5)` appears as a literal numerical fraction `54/5` with no derivation in this script and no citation to an upstream unit's symbolic computation. The Mathematica file uses the identical literal at line 67: `target = FullSimplify[54*gConst*cs^5/(5*a^5*cSpeed^5), ...]`. Per category definition this is a pre-baked symbolic form used as "the answer" without in-script derivation.

**Why this matters:**
If the upstream-derived coefficient is actually e.g. 27/5 or 54/25 in the upstream source, both engines would silently agree on the wrong target because each independently hard-codes the same literal. The two engines cannot catch each other on this.

**Required change:**
Either (a) add a comment block immediately above sympy line 106 (and mathematica line 67) citing the specific upstream stage script and the line at which `54/5` is *derived* (not just stated), e.g. `# 54/5 = ... (see scripts/moving_throat_pde_stageNNN_*.py line MM where Gamma_5 = a^5/(27 c_s^5) and the missing factor of (5/2)^... arises)`, OR (b) reconstruct the coefficient from the underlying upstream constants symbolically within this script (e.g. `target = sp.Rational(2,5) * mhat_coefficient_from_radiation * stuff` where each factor is named and tied to a stated origin). Option (a) is sufficient; option (b) is preferred. Use the same fix in both engines.

**Verification:**
The new comment or symbolic derivation must reference a specific file path and line; the verifier confirms the cited file/line exists and the coefficient claimed there matches.

### F5 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:106-111`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:67-71`

**What's wrong:**
Section II.2 builds `equation_residual = sp.simplify(mhat**2 * P0_compact - target)` (sympy line 108) and then merely `sp.pprint(equation_residual)` (line 111). There is no assertion of any kind. The Mathematica file at line 71 likewise only `Print[...]`s the residual. The docstring labels this as the "exact target equation" — the most physically meaningful claim of the entire script — yet the engines verify nothing about it.

**Why this matters:**
The only non-trivial content in the file (the target equation) is treated as decorative output. If the residual were generically nonzero but a typo flipped a sign in `target` or `P0_compact`, neither engine would alert. A reader of the PASS line in the saved output would believe the equation has been verified; it has not.

**Required change:**
Convert Section II.2 into a real assertion. Two acceptable forms:
1. **Solvability assertion (preferred):** solve `mhat**2 = target/P0_compact` for `mhat**2`, then assert at the F1 sample point that `mhat**2 > 0`. In sympy: `mhat_sq = sp.simplify(target / P0_compact); assert mhat_sq.subs(sample) > 0`. This catches sign errors and zero-divisor errors and is a real test (would fail if either side had a sign flip).
2. **Parameter-elimination assertion:** if the upstream context defines specific overlap amplitudes that DO null the residual, substitute them and `expect_zero("II.2 target residual on closure", equation_residual.subs(closure_subs))`. Only acceptable if the closure substitution is named in an accompanying comment with a file:line citation.

Pick form (1) unless a closure substitution is already documented in this script's comments (it is not).

**Verification:**
After the fix, saved output must show a new line `mhat^2 on sample = <positive rational>` and a PASS marker, OR (if form 2 is used) a `expect_zero("II.2 target residual on closure", 0)` PASS line.

### F6 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:126`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:80`

**What's wrong:**
The print at sympy line 126 / Mathematica line 80 states "If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0" but no executed line in either engine tests this implication. It is prose advertising a verification that does not happen.

**Why this matters:**
Same as F2 from a different angle: a sign error in `P^2/Delta^2` would not be caught. The implication itself is trivially true given the explicit formula `P0 = P^2/(Delta * D0)` if Delta > 0 and D0 > 0 — but only if the formula is the one currently written. The check should confirm the formula in force satisfies it.

**Required change:**
Already absorbed by F2 (add a numerical positivity check on the sample point). This finding documents the prose-vs-assertion gap; no additional Codex edits are required beyond F2's. If F2 lands, this finding is closed by the same patch — Codex should mark this directive block "Applied: F6" pointing to the F2 patch.

**Verification:**
The same sample-point positivity check that satisfies F2.

### F7 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:142-148`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:90-96`

**What's wrong:**
Section IV is titled "EXACT MONOTONIC DERIVATIVES" and the section banner is reproduced verbatim in Mathematica. The word "monotonic" implies a sign claim, but neither engine asserts a sign anywhere. The three executed identities are about *form* (the derivative equals a specific expression) not about *sign* (the derivative is negative / positive).

**Why this matters:**
Same as F3 from a different angle. If a downstream unit cites "from Stage 025: P0 is monotonically decreasing in K" as a fact, that fact is unsupported by the current script.

**Required change:**
Already absorbed by F3 (sign checks on the sample point). No additional Codex edits beyond F3.

**Verification:**
The same sample-point sign checks that satisfy F3.

### F8 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl` (whole file)
- compared to `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py` (whole file)

**What's wrong:**
The `.wl` file is structurally a line-by-line port of the `.py` file. Evidence:

1. Same module decomposition: `zero_frequency_coefficients` ↔ `zeroFrequencyCoefficients` (py:67-87 / wl:37-54), `normalization_formula` ↔ `normalizationFormula` (py:94-111 / wl:56-72), `stability_and_positivity` ↔ `stabilityAndPositivity` (py:118-126 / wl:74-81), `monotonic_derivatives` ↔ `monotonicDerivatives` (py:133-148 / wl:83-97). Same four functions in same order.

2. Identical literal banner strings, e.g. py:68 `banner("SECTION I — MINIMAL ISOTROPIC ZERO-FREQUENCY COEFFICIENTS")` ↔ wl:38 `banner["SECTION I — MINIMAL ISOTROPIC ZERO-FREQUENCY COEFFICIENTS"]`, including the same em-dash. Same in II.1, II.2, III, IV.

3. Identical symbol choreography in zero-frequency block:
   - py:70 `Delta = sp.simplify(OmegaU**2 * OmegaW**2 - R**2)` ↔ wl:39 `delta = FullSimplify[omegaU^2*omegaW^2 - r^2, ...]`
   - py:71 `Q = sp.simplify(GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2)` ↔ wl:40 `q = FullSimplify[gU^2*omegaW^2 + 2*gU*gW*r + gW^2*omegaU^2, ...]`
   - py:72 `P = sp.simplify(OmegaU**2 * GW + R * GU)` ↔ wl:41 `p = FullSimplify[omegaU^2*gW + r*gU, ...]`
   The term orderings (`omegaU^2 * omegaW^2` first, `R^2` subtracted; `gU^2*omegaW^2 + 2*gU*gW*r + gW^2*omegaU^2` in exactly that order) match exactly.

4. Same assertion list in same order in each section, with the same three formulas in Section IV (py:146-148 ↔ wl:94-96).

5. Section II.2 in both engines prints the residual without asserting it (py:111 / wl:71). An independent re-derivation would likely reach for `Resolve` / `Reduce` or `Solve` in Mathematica to actually test the equation; instead the wl file echoes the print.

This is precisely the mathematica_transliteration pattern: both engines walk the same algebraic moves rather than independently arriving at the same conclusion through different machinery. Both engines are vulnerable to the same systematic errors (the same hardcoded `54/5`, the same missing II.2 assertion, the same tautology pattern in every section).

**Why this matters:**
The two-engine policy exists so that an error in one engine's algebra (or in the author's translation of the physics) shows up as engine disagreement. When the two engines are line-by-line ports, the policy is defeated: they agree because they walk the same path, not because the answer is right. The current "engines_agree: true" status reflects translation fidelity, not independent confirmation.

**Required change:**
Re-derive the Mathematica side independently from the script's *stated premises* (the closed forms of `Delta, Q, P, B0, Z0, N0, D0` as scalar functions of the named parameters). Concretely:
1. Rewrite `zeroFrequencyCoefficients` to *not* mirror the SymPy term ordering. Use Mathematica idiom: `delta = Factor[omegaU^2 omegaW^2 - r^2]` (note: `Factor` not `FullSimplify`; sympy doesn't factor by default).
2. Compute `P0` via `Together[n0/d0]` and `Apart[Together[n0/d0], k]` rather than mirroring the python compact-form choice. Cross-check by `FullSimplify[%1 - %2] === 0`.
3. For Section IV, compute `dP0dK` via `Limit[(p0[k + h] - p0[k])/h, h -> 0]` (definition-of-derivative form) instead of `D[p0, k]`. Then test `FullSimplify[dP0dK - D[p0, k]] === 0`.
4. For Section II.2, add an independent solvability test via `Reduce[mhat^2 == target/p0Compact && mhat > 0, mhat, Reals]` and assert the resulting condition is not `False`.

If a re-derivation in this style produces the same final symbolic result as the SymPy script, engine agreement becomes meaningful. If not, an `engine_disagreement` finding gets raised in a future audit pass.

**Verification:**
The new Mathematica script must not contain the literal banner strings copied from the Python file verbatim (different em-dashes or restated section labels acceptable), must use at least one Mathematica primitive (`Factor`, `Apart`, `Limit`, `Reduce`) absent from the Python side, and must reach `Print["Stage 8 Mathematica audit passed."]` (or equivalent) only after independent assertions land.

## Independent-derivation check (Mathematica)

No. The Mathematica script is a transliteration of the SymPy script. Item-by-item evidence:

- Both compute `Delta` as `OmegaU^2*OmegaW^2 - R^2` with the same term ordering. `FullSimplify` in Mathematica naturally factors this as `(omegaU*omegaW - r)*(omegaU*omegaW + r)` (visible in the wl output line 17), which confirms Mathematica *could* have offered an independent factored form — but the script chooses to keep the unfactored definition matching SymPy's, then never uses the factorization downstream.

- Both define `P0_compact` with the same nesting `P^2/(Delta*(K*Delta - Delta*C^2/varpi^2 - Q))` (py:99, wl:60). A Mathematica-native re-derivation would more naturally `Together[n0/d0]` and inspect the resulting canonical form; instead the wl file hard-codes the same nested rewrite SymPy uses.

- Both engines have Section II.2 as `Print` only (no assertion). An independent author writing Mathematica from scratch would reach for `Solve` or `Reduce` rather than echo the SymPy choice to just pretty-print a residual.

- Same em-dash characters in banners and same exact section titles. Same four-Module / four-function decomposition. Same triple of derivative checks in identical order.

This is mathematica_transliteration (F8 above).

## Engine cross-check

Both engines reach `EXIT_CODE: 0` and report PASS on every assertion that exists. The expressions printed in Section II.2 are algebraically equivalent (sympy:35-48 vs mathematica:36) — the Mathematica form is the same rational function expressed via `(omegaU*omegaW - r)*(omegaU*omegaW + r)` rather than `(omegaU^2*omegaW^2 - r^2)`, but `FullSimplify[sympy_form - mathematica_form] === 0`. Engine agreement is real for what is computed; the problem is that what is computed in both engines is largely tautological. No engine_disagreement finding.

## Verdict justification

The scripts execute, exit 0, and print PASS for every check they make. They do not, however, verify the physics they claim to verify. Six of the seven SymPy assertions are tautological-by-construction (algebraic identities forced by the immediately preceding definition); the seventh (Section II.2's "exact target equation") is not asserted at all, only printed. The Mathematica script repeats the same content with the same gaps and uses identical banner strings and term orderings, qualifying it as a transliteration rather than an independent second engine. Attacks I tried that would have failed but did not: (a) injecting a sign flip in `N0` and re-reading — Section II.1, II.2, III, IV.1 would all still pass (because each is anchored to the same flipped expression); (b) replacing `54/5` with `27/5` — Section II.2 prints a residual either way, no assertion catches it; (c) replacing `D0 = K - B0 - Z0` with `D0 = K + B0 - Z0` — Section II.1's tautology still passes because P0_compact is defined to match `N0/D0` after the change. The verdict is `findings`, not `clean`. `stop_cold` does not apply: nothing here is unfixable, and the proposed fixes do not propagate to break downstream units (they add positivity / sign / solvability checks; they do not alter the symbolic form of P0 or the target equation coefficient).

## Self-test notes

- **Variable independence (directive step 1):** For the proposed F3 sign checks on `dP0/dK` and `dP0/dX`, I confirmed `N0 = P^2/Delta^2` depends only on GU, GW, OmegaU, OmegaW, R (none of K, X), and the denominator `K - X - Q/Delta` depends linearly on K and X. So `dP0/dK` and `dP0/dX` are nontrivially nonzero at the sample point (denominator squared is positive; numerator is `P^2 > 0` whenever P != 0). No identically-zero-derivative trap.
- **Symmetry/parity:** No unbounded integrals in this unit; the parity trap does not apply.
- **Trivial-case pre-check:** For the F1-F3 sample point `K=2, varpi=1, C=1, OmegaU=OmegaW=2, R=1, GU=GW=1`: `Delta = 4*4 - 1 = 15`, `Q = 1*4 + 2*1*1*1 + 1*4 = 10`, `P = 4*1 + 1*1 = 5`, `B0 = 1`, `Z0 = 10/15 = 2/3`, `D0 = 2 - 1 - 2/3 = 1/3`, `N0 = 25/225 = 1/9`, `P0 = (1/9)/(1/3) = 1/3 > 0`. Positive, sensible. `dP0/dK = -N0/(K-X-Q/Delta)^2 = -(1/9)/((2-0-2/3)^2)`; with X = C^2/varpi^2 = 1 the denominator (K - X - Q/Delta) is `2 - 1 - 2/3 = 1/3` and `dP0/dK = -(1/9)/(1/9) = -1 < 0`. `dP0/dX = +1 > 0`. Sum is zero. All proposed sign-checks evaluate to definite, nonzero, correctly-signed values at this sample point — no silent-pass risk.
- **Path specifications:** No missing-script findings; both engines exist; F4 and F8 cite paths to existing files.
