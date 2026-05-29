---
unit_id: 073
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T23:08:19-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 073

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check (Mathematica precedence bug)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:47`

**Issue:**
The line
```
expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37];
```
parses as `eta /. (len/ell) -> (37 - 37)` because `Rule` has lower precedence than `Plus` in Mathematica. That makes the substitution `len/ell -> 0`, so the residual handed to `expectZero` is `0` regardless of the value 37. The check is vacuous and the saved output line `eta(reference) - 37 = 0` is misleading.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`, replace the single line at line 47:

Before:
```
expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37];
```

After:
```
expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];
```

Add parentheses so the rule's RHS is just `37`, and the `- 37` is applied to the result of the substitution. No other lines change.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 073` and confirm:
- Line 47 of the .wl now reads `expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];`
- The saved output still shows `eta(reference) - 37 = 0` (residual genuinely 0, not vacuous-0)
- The script still exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`
- summary: Corrected the eta reference substitution so the final `- 37` is outside the Mathematica rule.
- deviation: none

## F2 — tautological_check (eta - L/ell is automatic cancellation)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:50-54`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:39-47`

**Issue:**
Both scripts define `eta = (T_X/ell)*L/T_X` and immediately assert `eta - L/ell == 0`. The `T_X` cancels in any CAS regardless of what physical content the closure `K_m = T_X/ell` is supposed to encode. To make the assertion load-bearing, define `eta_sym = K_m * L / T_X` symbolically first, then substitute `K_m -> T_X/ell`, then assert.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`, replace the block at lines 49-54.

Before (lines 49-54):
```python
# Robin closure K_m = T_X/ell.
K_m, T_X, L, ell = sp.symbols('K_m T_X L ell', positive=True, real=True)
eta = sp.simplify((T_X / ell) * L / T_X)
print("eta under K_m = T_X/ell ->", eta)
expect_zero("eta - L/ell", eta - L / ell)
expect_zero("eta(reference) - 37", eta.subs({L / ell: 37}) - 37)
```

After:
```python
# Robin closure K_m = T_X/ell. Build eta symbolically in K_m first so the
# assertion eta - L/ell == 0 actually exercises the closure substitution and
# not a trivial T_X cancellation.
K_m, T_X, L, ell = sp.symbols('K_m T_X L ell', positive=True, real=True)
eta_sym = K_m * L / T_X
eta = sp.simplify(eta_sym.subs(K_m, T_X / ell))
print("eta under K_m = T_X/ell ->", eta)
expect_zero("eta - L/ell", eta - L / ell)
expect_zero("eta(reference) - 37", eta.subs({L / ell: 37}) - 37)
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`, replace the block at lines 39-47.

Before (lines 39-47):
```
Clear[km, tx, len, ell];
$Assumptions =
  Element[{km, tx, len, ell}, Reals] &&
  km > 0 && tx > 0 && len > 0 && ell > 0;

eta = FullSimplify[(tx/ell)*len/tx, Assumptions -> $Assumptions];
Print["eta under K_m = T_X/ell -> ", fmt[eta]];
expectZero["eta - L/ell", eta - len/ell];
expectZero["eta(reference) - 37", eta /. (len/ell) -> 37 - 37];
```

After (note: the F1 fix is folded into this block; if F1 was applied before F2, line 47 already has the parentheses — keep them):
```
Clear[km, tx, len, ell];
$Assumptions =
  Element[{km, tx, len, ell}, Reals] &&
  km > 0 && tx > 0 && len > 0 && ell > 0;

etaSym = km*len/tx;
eta = FullSimplify[etaSym /. km -> tx/ell, Assumptions -> $Assumptions];
Print["eta under K_m = T_X/ell -> ", fmt[eta]];
expectZero["eta - L/ell", eta - len/ell];
expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 073` and `redteam exec-mathematica 073` and confirm:
- The sympy script defines `eta_sym = K_m * L / T_X` and computes `eta` via `.subs(K_m, T_X/ell)` before the `expect_zero("eta - L/ell", ...)` line.
- The mathematica script defines `etaSym = km*len/tx` and computes `eta` via `/. km -> tx/ell` before `expectZero["eta - L/ell", ...]`.
- Both saved outputs still show `eta - L/ell = 0` and `eta(reference) - 37 = 0`.
- Both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`
  - `mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`
- summary: Built eta from symbolic `K_m`/`km` first, then substituted the Robin closure before checking `eta - L/ell`.
- deviation: none

## F3 — hardcoded_result (geometry-map literals lack a symbolic identity)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:34-47`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:26-37`

**Issue:**
Both scripts hard-code `epsilon_r = 1/20` and `Lambda_star = L/a = 37/20`, then assert `Lambda_ell = (L/a)/(ell/a) - 37 == 0`. The check is pure arithmetic between hand-typed rationals. The substantive content — that `Lambda_ell = (L/a)/(ell/a)` regardless of the specific reference branch — is never asserted. Add a symbolic identity check above the numerical specialization.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`, insert a new block between line 34 (the `banner` call) and line 36 (the `epsilon_r = sp.Rational(1, 20)` line). Specifically, after the existing line 34 (`banner("STAGE 56 — FAMILY-1 GEOMETRY MAP")`), add a blank line and then:

```python

# Symbolic identity: Lambda_ell = (L/a) / (ell/a) = L/ell, independent of the
# specific reference-branch values chosen below.
L_sym, a_sym, ell_sym = sp.symbols('L a ell', positive=True)
Lambda_star_sym = L_sym / a_sym
ell_over_a_sym = ell_sym / a_sym
Lambda_ell_sym = sp.simplify(Lambda_star_sym / ell_over_a_sym)
expect_zero("Lambda_ell - L/ell (symbolic)", Lambda_ell_sym - L_sym / ell_sym)

```

Leave the existing block from line 36 onward (`epsilon_r = sp.Rational(1, 20)` etc.) unchanged. The new block adds the symbolic check above; the existing numerical check stays below.

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`, insert a new block between line 26 (the `banner["STAGE 056 — FAMILY-1 GEOMETRY MAP"];` call) and line 28 (the `epsilonR = 1/20;` line). After line 26, add a blank line and then:

```
Clear[lSym, aSym, ellSym];
Module[{lambdaStarSym, ellOverASym, lambdaEllSym},
  lambdaStarSym = lSym/aSym;
  ellOverASym = ellSym/aSym;
  lambdaEllSym = FullSimplify[lambdaStarSym/ellOverASym,
    Assumptions -> lSym > 0 && aSym > 0 && ellSym > 0];
  expectZero["Lambda_ell - L/ell (symbolic)", lambdaEllSym - lSym/ellSym];
];

```

Leave the existing block from line 28 onward unchanged.

**Verification command:**
After Codex applies, the verifier will run both engines and confirm:
- A new line `Lambda_ell - L/ell (symbolic) = 0` appears in each saved output, followed by `PASS: Lambda_ell - L/ell (symbolic)`.
- The existing numerical line `Lambda_ell - 37 = 0` and `PASS: Lambda_ell - 37` still appear after the new symbolic line.
- Both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`
  - `mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`
- summary: Added the symbolic `Lambda_ell = (L/a)/(ell/a) = L/ell` identity check before the numerical specialization in both scripts.
- deviation: none

## F4 — mathematica_transliteration (INFORMATIONAL — DO NOT APPLY)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl` (entire file)

**Issue:**
The `.wl` is a line-by-line port of the `.py`: same variable choreography (`epsilon_r ↔ epsilonR`, `Lambda_star ↔ lambdaStar`, `K_m, T_X, L, ell ↔ km, tx, len, ell`), same banner text, same assertion order. The Mathematica engine does not derive the result independently from the physical premises.

**Required change:**
**Codex: do nothing for F4.** Mechanically refactoring the Mathematica script into an independent derivation would constitute a "scope extension" (forbidden by the directive preamble) and requires the user's judgment about what an independent derivation should look like for this particular unit. Skip F4 and append an `## Applied: F4` block stating `files_changed: []`, `summary: informational only — no changes per directive`, `deviation: none`.

**Verification command:**
None — informational.

## Applied: F4

- files_changed: []
- summary: informational only — no changes per directive
- deviation: none
