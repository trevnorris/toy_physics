---
unit_id: 125
batch: IV.3
created_at: 2026-05-29T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T20:30:50Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 125

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:83-85` (the `g_a_large` definition, its print line, and the weak `abs(...) < 1/20` assertion)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:93` — the peaked-source branch (the `Limit[gA, aSym -> Infinity]` check); new range assertions are inserted immediately after this line.

**Issue:**
The SymPy peaked-source assertion at line 85 is
```python
expect_true("g[peaked@L proxy a=100] < 0.05", bool(abs(g_a_large) < sp.Rational(1, 20)))
```
This is a *smallness* test: it only checks proximity to zero. Because it wraps `g_a_large` in `abs(...)`, it accepts small NEGATIVE moments — e.g. a value of `-0.01` passes — even though such a value violates the paper's lower-bound half of the theorem `0 <= g[sigma] <= 1` (paper/stages/stage_125.tex:16; notes §"Exact positivity theorem" :50-55). The peaked source `sigma_a` at `a=100` is exactly the family member that drives `g` toward its lower endpoint, so it is the natural place to exercise the lower bound `g >= 0`; the current check does not. The uniform-source checks at lines 87-88 already assert a genuine lower bound (`>= 0`) and upper bound (`<= 1`); the peaked source needs the analogous genuine range assertion.

Note on non-tautology: `g_a` is the closed-form result of `sp.integrate(sigma_profile * cos(pi z / 2L), ...)` then `.subs(sigma_a, 100)` and `sp.N(...)` (lines 76, 83). Although the *abstract* integral of a nonnegative profile against a nonnegative kernel is `>= 0` by construction, the script does NOT integrate numerically at assertion time — it evaluates a SymPy-derived closed form. A sign error, a wrong kernel, a wrong normalization, or a mis-substituted `a` in that closed form would yield a negative numeric `g_a_large`, which the genuine range check below would FAIL. (Line 79's `g_a.subs(0) - 2/pi == 0` only pins the uniform endpoint; the peaked endpoint is currently unguarded against a sign error.) The check is therefore non-trivial as an implementation guard, which is exactly the finding's concern.

**Required change (SymPy):**

Replace lines 83-85.

Before (lines 83-85):
```python
g_a_large = sp.N(g_a.subs(sigma_a, 100))
print("g_a at sigma_param = 100 (peaked-at-L proxy) =", g_a_large)
expect_true("g[peaked@L proxy a=100] < 0.05", bool(abs(g_a_large) < sp.Rational(1, 20)))
```

After (lines 83-85):
```python
g_a_large = sp.N(g_a.subs(sigma_a, 100))
print("g_a at sigma_param = 100 (peaked-at-L proxy) =", g_a_large)
# Genuine range check of the paper bound 0 <= g[sigma] <= 1 at the peaked endpoint.
# NOT wrapped in abs(): a small NEGATIVE moment (a sign error in the closed form)
# now FAILS the lower bound, instead of passing a mere smallness test.
expect_true("g[peaked@L proxy a=100] >= 0", bool(g_a_large >= 0))
expect_true("g[peaked@L proxy a=100] <= 1", bool(g_a_large <= 1))
# Retain the smallness fact (peaked source biases g toward the lower endpoint):
# a strictly positive value below 1/20 confirms the trend without admitting negatives.
expect_true("g[peaked@L proxy a=100] < 1/20", bool(g_a_large < sp.Rational(1, 20)))
```

Rationale for the tolerances: `0` and `1` are the literal paper bounds (notes :52-54, paper :16), not fabricated numbers. The `< 1/20` line preserves the original `0.05` smallness fact from the prior assertion (the saved transcript shows `g_a_large = 0.0153964...`, well inside `[0, 1/20)`), but it is now an ADDITION to — not a substitute for — the genuine `>= 0` / `<= 1` range. The load-bearing new assertion is `g[peaked@L proxy a=100] >= 0`: it is the half the old check could not detect, and it fails for any negative `g_a_large`.

**Required change (Mathematica):**

Target: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:93`.

The Mathematica peaked-source check at line 93 is NOT the same weak `abs(...) < tol` smallness test the SymPy script had — it is a genuine endpoint identity:
```wolfram
expectZero["moment g[peaked@L] limit", Limit[gA, aSym -> Infinity]];
```
This asserts the exact limit `g -> 0` as `a -> Infinity`, which already exercises the lower endpoint. Keep this line unchanged. However, to MIRROR the SymPy strengthening — a genuine `0 <= g <= 1` *range* bound at the peaked source (parallel to the uniform-source range bounds at lines 95-96) — add two range assertions for the same finite `a = 100` proxy used on the SymPy side, immediately AFTER line 93.

Insert after line 93 (before the blank `Print[""]`-less gap and the uniform-source range checks at lines 95-96):
```wolfram
gApeaked = N[(gA /. aSym -> 100), 30];
Print["g_a at aSym -> 100 (peaked-at-L proxy) = ", fmt[gApeaked]];
expectTrue["g[peaked@L proxy a=100] >= 0", gApeaked >= 0];
expectTrue["g[peaked@L proxy a=100] <= 1", gApeaked <= 1];
```

Notes on the Mathematica idioms used here:
- `(gA /. aSym -> 100)` is parenthesized BEFORE applying `N[..., 30]`, because `Rule` (`->`) binds looser than function application/arithmetic; without the parentheses the substitution would not be grouped as intended.
- `N[..., 30]` forces a 30-digit numeric value so the `>= 0` / `<= 1` comparisons resolve to explicit `True`/`False` under `expectTrue`'s `FullSimplify`/`TrueQ` gate, rather than staying symbolic.
- `gApeaked >= 0` is NOT wrapped in `Abs[...]`: a negative numeric value (sign error in `gA`) FAILS the lower bound, matching the SymPy guard.
- No comment body contains the substring `*)`.

**Claim manifest:**
- **M1** — `expect_true("g[peaked@L proxy a=100] >= 0", ...)` (SymPy) and `expectTrue["g[peaked@L proxy a=100] >= 0", ...]` (Mathematica) exercise the LOWER half of the paper's positivity theorem `0 <= g[sigma]` for the peaked positive normalized source `sigma_a(z) = (a+1)(z/L)^a/L` at `a=100`. Paper anchor: `paper/stages/stage_125.tex:16` (boxed `0 <= g[sigma] <= 1`). Notes anchor: `notes/stages/moving_throat_pde_stage125_positive_source_theorem.md:50-55` (boxed `0 <= g[sigma] <= 1`), with the kernel-in-`[0,1]` justification at `:45-49`.
- **M2** — `expect_true("g[peaked@L proxy a=100] <= 1", ...)` (both engines) exercises the UPPER half `g[sigma] <= 1` for the same peaked source. Same anchors as M1 (`stage_125.tex:16`; `..._positive_source_theorem.md:50-55`).
- **M3** — `expect_true("g[peaked@L proxy a=100] < 1/20", ...)` (SymPy only) retains the original smallness fact that the peaked source biases `g` toward the lower endpoint; this is a trend witness, not a paper-quoted constant, and is subordinate to M1/M2. The `1/20` value is carried over verbatim from the pre-edit assertion (no new number introduced).

**Verification command:**
The verifier will run `redteam exec-sympy 125` and `redteam exec-mathematica 125`. Expected:
- SymPy transcript shows the unchanged print line `g_a at sigma_param = 100 (peaked-at-L proxy) = 0.0153964...` followed by THREE new check lines:
  - `g[peaked@L proxy a=100] >= 0 = True`
  - `g[peaked@L proxy a=100] <= 1 = True`
  - `g[peaked@L proxy a=100] < 1/20 = True`
  and the old single line `g[peaked@L proxy a=100] < 0.05 = True` is GONE.
- Mathematica transcript shows the unchanged `moment g[peaked@L] limit = 0` / `PASS` line, then the new print `g_a at aSym -> 100 (peaked-at-L proxy) = ...` and two new `PASS:` lines:
  - `PASS: g[peaked@L proxy a=100] >= 0`
  - `PASS: g[peaked@L proxy a=100] <= 1`
- All other check lines unchanged; both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl`
- summary: Replaced the peaked-source SymPy smallness-only assertion with genuine range checks plus retained smallness, and added matching finite-peaked-source range checks to the Mathematica audit.
- deviation: The Mathematica finite peaked-source value uses extra precision and extracts the real part to avoid cancellation artifacts in the closed-form incomplete-gamma expression.
