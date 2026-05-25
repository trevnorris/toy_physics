---
unit_id: 073
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: false
  outputs_fresh: true
---

# Audit unit 073 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`
- sympy output: (missing)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.txt`

## What the script claims to verify

Both scripts assert an explicit Family-1 geometry map: with `epsilon_r = ell/a = 1/20` and a "carried reference branch" `L/a = 37/20`, the ratio `Lambda_ell = (L/a)/(ell/a)` equals 37. They also assert that, under the local Robin mouth closure `K_m = T_X/ell`, the Robin variable `eta = K_m L / T_X` reduces algebraically to `L/ell` and (under the substitution `L/ell -> 37`) numerically equals 37. The "verification" is therefore a 3-step arithmetic check on two hand-typed rational numbers plus one cancellation of `T_X`. No derivation of `epsilon_r`, `L/a`, or the closure `K_m = T_X/ell` is performed in-script.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 47 | `expect_zero("Lambda_ell - 37", Lambda_ell - 37)` where `Lambda_ell = (37/20)/(1/20)` | no |
| A2 | sympy | 53 | `expect_zero("eta - L/ell", eta - L/ell)` where `eta = simplify((T_X/ell)*L/T_X)` | no |
| A3 | sympy | 54 | `expect_zero("eta(reference) - 37", eta.subs({L/ell: 37}) - 37)` | partial |
| A4 | mathematica | 37 | `expectZero["Lambda_ell - 37", lambdaEll - 37]` where `lambdaEll = (37/20)/(1/20)` | no |
| A5 | mathematica | 46 | `expectZero["eta - L/ell", eta - len/ell]` where `eta = FullSimplify[(tx/ell)*len/tx]` | no |
| A6 | mathematica | 47 | `expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37]` | **broken: trivially 0 by precedence bug** |

A1/A4 are arithmetic identities between literals (`(37/20)/(1/20)` vs `37`), with no physics input. A2/A5 are algebraic identities by construction: `eta` is defined as `(T_X/ell)*L/T_X` and immediately simplified, so `eta - L/ell` is zero by cancellation, independent of the geometry. A6 is broken (see F1). Only A3 actually tests a substitution.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:47`

**What's wrong:**
The line reads
```
expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37];
```
Mathematica's `Rule` (`->`) has lower precedence than `Plus`, so the RHS of the rule is `37 - 37`, which evaluates to `0`. The expression handed to `expectZero` is therefore `eta /. (len/ell) -> 0`, which simplifies to `0` (since `eta == len/ell`). The check trivially passes, and the captured output line `eta(reference) - 37 = 0` reflects the substituted-to-zero residual, **not** the intended `(eta at L/ell = 37) - 37`. The mathematica script never actually exercises the numerical claim `eta(reference) = 37`.

The corresponding sympy line 54
```
expect_zero("eta(reference) - 37", eta.subs({L / ell: 37}) - 37)
```
substitutes `L/ell -> 37` and then subtracts 37, which is the intended check. So the two engines produce `0` for *different* expressions, masking a real bug behind a coincidence (both happen to be zero).

**Why this matters:**
The mouth-closure numerical fingerprint `eta = 37` is the one nontrivial output of this unit. The Mathematica audit's assertion that this number is correct is vacuous: it would pass for any `eta` whatsoever, including `eta = 0`, `eta = -100`, or `eta = anything/anything`. The unit's `is_checkpoint: False` status still requires both engines to verify substantive claims; this engine does not.

**Required change:**
At `mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:47`, replace
```
expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37];
```
with
```
expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];
```
The added parentheses bind the `- 37` outside the rule's RHS so the substitution `L/ell -> 37` happens first, after which `37 - 37 == 0` is the residual the engine actually checks.

**Verification:**
After the fix, re-running the .wl produces the output `eta(reference) - 37 = 0` (same surface text) but the residual is now `(len/ell -> 37) - 37 = 37 - 37 = 0`, i.e. a real check. A sanity perturbation: replacing the `37` inside the rule with `36` should fail (residual `-1`); under the current buggy line, that change has *no effect* on the residual because the RHS is still `36 - 37 = -1`, which substitutes `len/ell -> -1`, giving `(-1) - 37` — the bug would surface as residual `-38` rather than `-1`, demonstrating the precedence trap.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:52-53`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:44-46`

**What's wrong:**
Both scripts define
```
eta = sp.simplify((T_X / ell) * L / T_X)        # sympy line 52
eta = FullSimplify[(tx/ell)*len/tx, ...]         # mathematica line 44
```
and then immediately assert
```
expect_zero("eta - L/ell", eta - L / ell)        # sympy line 53
expectZero["eta - L/ell", eta - len/ell]         # mathematica line 46
```
But `(T_X/ell) * L / T_X` is algebraically identical to `L/ell` for any nonzero `T_X` — the cancellation is automatic in both engines. The assertion can never fail; it confirms that the simplifier cancelled `T_X`, not that the Robin closure `K_m = T_X/ell` is correctly applied to the definition of `eta`.

**Why this matters:**
The physical content of the line is supposed to be "under the closure `K_m = T_X/ell`, the Robin variable `eta = K_m L / T_X` reduces to `L/ell`". A non-tautological check would carry `eta` symbolically as a function of an unspecified `K_m` (e.g., `eta_sym = K_m * L / T_X`) and only then substitute `K_m -> T_X/ell` to verify the cancellation. The current form bakes the closure into the definition before the assertion, so the assertion only re-checks elementary simplification.

**Required change:**
At sympy line 52, separate the definition from the closure substitution:
```python
K_m, T_X, L, ell = sp.symbols('K_m T_X L ell', positive=True, real=True)
eta_sym = K_m * L / T_X
eta = sp.simplify(eta_sym.subs(K_m, T_X / ell))
```
At mathematica line 44, do the same:
```
Clear[km, tx, len, ell];
$Assumptions = Element[{km, tx, len, ell}, Reals] && km > 0 && tx > 0 && len > 0 && ell > 0;
etaSym = km*len/tx;
eta = FullSimplify[etaSym /. km -> tx/ell, Assumptions -> $Assumptions];
```
Leave the subsequent `expect_zero("eta - L/ell", eta - L/ell)` / `expectZero["eta - L/ell", eta - len/ell]` line in place. Now the residual is `(K_m L/T_X with K_m = T_X/ell) - L/ell`, which exercises the substitution rather than a self-cancellation.

**Verification:**
The line number of the assertion doesn't move, but the LHS of the residual is now built via `subs(K_m, T_X/ell)` instead of via literal substitution at definition time. A sanity perturbation: changing the closure to `K_m -> 2 T_X/ell` should fail with residual `L/ell` (in sympy) / `len/ell` (in Mathematica), confirming the assertion is now load-bearing.

### F3 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:36-47`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:28-37`

**What's wrong:**
Both scripts hard-code
```
epsilon_r = 1/20            # ell/a, given without derivation
Lambda_star = 37/20         # L/a, given without derivation
```
and then assert `Lambda_ell = (L/a)/(ell/a) - 37 == 0`. This is the arithmetic identity `(37/20)/(1/20) = 37`. Neither script cites where `1/20` and `37/20` come from (e.g., a previous stage, a paper equation, a chosen reference geometry), nor does either expose the cancellation algebraically: with a symbolic `Lambda_star` and a symbolic `epsilon_r`, the identity `Lambda_ell = Lambda_star / epsilon_r` would be the substantive content.

**Why this matters:**
The "ledger" entry `Lambda_ell = 37` is the one number this unit injects into downstream stages. As written, the check confirms only that python/Mathematica can compute `(37/20)/(1/20)`. If a later edit changed `Lambda_star` to, say, `38/20` to match a doc change, the assertion would silently follow along (the assert is `Lambda_ell - 37 == 0`, but the constant `37` is itself derived from the hard-coded inputs — so editing both `Lambda_star` and the literal `37` together would still pass while flipping the downstream value). At minimum the script should record the symbolic identity `Lambda_ell == Lambda_star * a / ell` before specializing.

**Required change:**
At sympy lines 36-47, add a symbolic check above the numerical one:
```python
# Symbolic identity: Lambda_ell = (L/a) / (ell/a) = L/ell, independent of the
# specific reference branch chosen.
L_sym, a_sym, ell_sym = sp.symbols('L a ell', positive=True)
Lambda_star_sym = L_sym / a_sym
ell_over_a_sym = ell_sym / a_sym
Lambda_ell_sym = sp.simplify(Lambda_star_sym / ell_over_a_sym)
expect_zero("Lambda_ell - L/ell (symbolic)", Lambda_ell_sym - L_sym / ell_sym)

# Numerical specialization at the reference branch (epsilon_r = 1/20, L/a = 37/20).
epsilon_r = sp.Rational(1, 20)
Lambda_star = sp.Rational(37, 20)
...
```
At mathematica lines 28-37, mirror this with `Clear[lSym, aSym, ellSym]; lambdaStarSym = lSym/aSym; ...; expectZero["Lambda_ell - L/ell (symbolic)", lambdaEllSym - lSym/ellSym];` *before* the numerical block. Keep the existing numerical assertion `Lambda_ell - 37 == 0` afterward (it remains arithmetic but is no longer the only line of defense).

**Verification:**
After the patch, both scripts contain a new symbolic-check assertion above the numerical one. The verifier looks for a line `Lambda_ell - L/ell (symbolic) = 0` in the saved output. A sanity perturbation: replacing the symbolic numerator with `2 L_sym/a_sym` should fail with residual `L/ell` (sympy) or `lSym/ellSym` (Mathematica), confirming the symbolic check is load-bearing.

### F4 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:26-47`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:34-54`

**What's wrong:**
The `.wl` is a line-by-line transliteration of the `.py`, not an independent re-derivation. Compare:

Sympy lines 36-47:
```python
epsilon_r = sp.Rational(1, 20)
Lambda_star = sp.Rational(37, 20)
ell_over_a = epsilon_r
Lambda_ell = sp.simplify(Lambda_star / ell_over_a)
...
expect_zero("Lambda_ell - 37", Lambda_ell - 37)
```
Mathematica lines 28-37:
```
epsilonR = 1/20;
lambdaStar = 37/20;
ellOverA = epsilonR;
lambdaEll = FullSimplify[lambdaStar/ellOverA];
...
expectZero["Lambda_ell - 37", lambdaEll - 37];
```
Sympy lines 50-54:
```python
K_m, T_X, L, ell = sp.symbols('K_m T_X L ell', positive=True, real=True)
eta = sp.simplify((T_X / ell) * L / T_X)
expect_zero("eta - L/ell", eta - L / ell)
expect_zero("eta(reference) - 37", eta.subs({L / ell: 37}) - 37)
```
Mathematica lines 39-47:
```
Clear[km, tx, len, ell];
$Assumptions = ... km > 0 && tx > 0 && len > 0 && ell > 0;
eta = FullSimplify[(tx/ell)*len/tx, Assumptions -> $Assumptions];
expectZero["eta - L/ell", eta - len/ell];
expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37];
```
Same variable choreography, same intermediate names mapped 1:1 (`epsilon_r ↔ epsilonR`, `Lambda_star ↔ lambdaStar`, `K_m, T_X, L, ell ↔ km, tx, len, ell`), same order of assertions, same banner text. The second engine adds no independent path.

**Why this matters:**
The second-engine policy exists so that a bug in either engine's algebra cannot pass undetected. Transliteration with identical structure defeats that: the precedence bug in F1 above is a perfect example — it survived only because the Mathematica script wasn't re-derived from the physics, only re-typed from the SymPy listing. (Whether F1 alone would have been caught by an independent derivation is uncertain, but the transliteration made it more likely to slip through.)

**Required change:**
This finding is informational and is correctable only by rewriting the Mathematica script as an independent re-derivation, which goes beyond the auditor's mandate ("Do not propose new features, refactors, or scope extensions"). Codex should **not** rewrite the file; the user handles structural transliteration concerns manually. This finding is recorded so that the post-batch tracker (`MATHEMATICA_MIRROR_POLICY.md`) reflects the issue. Apply only F1, F2, F3 mechanically.

**Verification:**
No mechanical verification — informational only.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py` (see F4). Both scripts hard-code `1/20` and `37/20`, both define `eta = (T_X/ell) * L / T_X` and simplify, both check `eta - L/ell` and a numerical specialization. The only difference is syntax (`sp.Rational(1,20)` vs `1/20`, `sp.simplify` vs `FullSimplify`, `expect_zero` vs `expectZero`). No independent path is taken; the Mathematica precedence bug in line 47 (F1) would have been visible immediately if the Mathematica author had derived the check from scratch.

## Engine cross-check

Sympy output is missing, so direct side-by-side comparison is not possible from saved files. From reading the scripts, the expected outputs are:

| Check | Sympy expected residual | Mathematica saved residual |
|---|---|---|
| `Lambda_ell - 37` | 0 | 0 (line 17 of `.txt`) |
| `eta - L/ell` | 0 | 0 (line 21) |
| `eta(reference) - 37` | 0 (correctly: `37 - 37`) | 0 (incorrectly: substituted `len/ell -> 0`, see F1) |

Both engines produce zero for `eta(reference) - 37`, but for **different expressions** — the sympy check is real, the Mathematica check is vacuous. This is recorded as `engines_agree: false` in the front-matter despite the matching residuals, because the assertion in the Mathematica script does not evaluate the same residual the assertion in the SymPy script evaluates.

## Verdict justification

The unit's claims are very narrow (three lines of arithmetic and one cancellation), but the assertions are weak: A1/A4 are arithmetic between literals, A2/A5 are algebraic identities by construction, and A6 contains a Mathematica precedence bug that makes the only nontrivial-looking check trivially pass. A3 is the one real assertion. The unit is `is_checkpoint: False` and `is_status_only_candidate: False`, so both engines need substantive checks. Findings F1 (bug, high) and F2/F3 (tautology/hardcoding, medium) are mechanically applicable; F4 (transliteration) is informational.

I attacked the scripts by: (a) tracing the simplification of `eta` to confirm `(T_X/ell)*L/T_X → L/ell` is automatic in both engines (it is); (b) checking Mathematica operator precedence for `Rule` vs `Plus` to confirm `(len/ell) -> 37 - 37` parses as `(len/ell) -> 0` (it does — `Plus` precedence 310 is higher than `Rule` precedence 120); (c) checking that the sympy `.subs({L/ell: 37})` succeeds against the simplified `eta = L/ell` (it does — sympy can pattern-match the `L*ell**-1` factor); (d) checking whether the assumption set `positive=True, real=True` in sympy and `Reals && >0` in Mathematica are compatible with the geometric setup (they are — lengths and stresses positive). No `UNFIXABLE` or `CRITICAL_DOWNSTREAM` flags warranted: `Lambda_ell = 37` and `eta = 37` are the cited numbers downstream, and the fixes preserve those numbers (only the *checks* become substantive).

## Self-test notes

I checked (1) variable independence: F2 and F3 require introducing symbols `K_m, T_X, L, ell` (already declared) and new symbols `L_sym, a_sym, ell_sym` (for symbolic Lambda_ell check); each `subs` operates on a symbol the target expression actually depends on. (2) Parity / domain: no integrals or unbounded domains in this unit, so parity arguments don't apply. (3) Trivial-case pre-check: with the F2 patch, substituting `K_m -> T_X/ell` into `K_m*L/T_X` gives `(T_X/ell)*L/T_X = L/ell`, so `eta - L/ell = 0` (correct PASS); substituting `K_m -> 2*T_X/ell` would give `2L/ell - L/ell = L/ell ≠ 0` (correct FAIL under perturbation). With the F1 patch, `(eta /. (len/ell) -> 37) - 37 = 37 - 37 = 0` (correct PASS); replacing `37` inside the rule with `36` gives `36 - 37 = -1 ≠ 0` (correct FAIL under perturbation). (4) Path specifications: no `missing_verification_script` findings, so target paths are existing-file edits at named line numbers.
