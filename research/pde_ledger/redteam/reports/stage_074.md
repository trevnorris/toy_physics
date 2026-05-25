---
unit_id: 074
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 074 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt`

## What the script claims to verify

Per docstring, the unit verifies the healing/compliance width lock: starting from the GNLS healing length `ell = hbar/(2 m c_s)`, the dimensionless support scale `chi_s` reduces to `Lambda_ell/2` (where `Lambda_ell = L/ell`), and the Family-1 branch coefficient evaluates `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2 = (9/5) Lambda_ell^2`. The reference branch (`Lambda_ell = 37`) is then checked to give `chi_s = 37/2`, `kappa = 12321/5`, `alpha = sqrt(kappa) = 111*sqrt(5)/5`. The assertions are: `chi_s - Lambda_ell/2 == 0`, `kappa - (9/5) Lambda_ell^2 == 0`, `chi_ref - 37/2 == 0`, `kappa_ref - 12321/5 == 0`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 39 | `simplify(chi_lock - Lambda_ell/2) == 0` where `chi_lock := Lambda_ell/2` | no (tautological) |
| A2 | sympy | 40 | `simplify(kappa_lock - (9/5) Lambda_ell^2) == 0` with `kappa_lock = 4*(Lambda_ell/2)^2 + (4/5)*Lambda_ell^2` | partial (real arithmetic identity, but no physics derivation) |
| A3 | sympy | 54 | `simplify(chi_ref - 37/2) == 0` where `chi_ref = (Lambda_ell/2).subs(Lambda_ell, 37)` | no (tautological substitution) |
| A4 | sympy | 55 | `simplify(kappa_ref - 12321/5) == 0` where `kappa_ref = ((9/5)*Lambda_ell^2).subs(Lambda_ell, 37)` after the A2 identity collapses it | partial (numeric specialization of A2) |
| A5 | mathematica | 42 | `expectZero[chi_s - Lambda_ell/2]` where `chiLock` is derived from `m c_s L / hbar` with `c_s -> hbar/(2 m ell)` and `L/ell -> Lambda_ell` | yes |
| A6 | mathematica | 43 | `expectZero[kappa - (9/5) Lambda_ell^2]` with `kappaLock = 4*chiLock^2 + (4/5)*lambdaEll^2` and derived `chiLock` | yes (since chiLock came from a real substitution chain) |
| A7 | mathematica | 57 | `expectZero[chi_ref - 37/2]` | partial |
| A8 | mathematica | 58 | `expectZero[kappa_ref - 12321/5]` | partial |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:33,39`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:43,54`

**What's wrong:**
Line 33 defines

    chi_lock = sp.simplify(Lambda_ell / 2)

and then line 39 asserts

    expect_zero("chi_s - Lambda_ell/2", chi_lock - Lambda_ell / 2)

This expands to `Lambda_ell/2 - Lambda_ell/2 == 0`, which is algebraically guaranteed by the definition on line 33. The assertion cannot fail regardless of the physics. The output transcript line `chi_s - Lambda_ell/2 = 0` simply records that subtraction.

The same pattern repeats for the reference branch: line 43 sets `chi_ref = chi_lock.subs(Lambda_ell, 37) = 37/2` and line 54 asserts `chi_ref - 37/2 == 0` — also tautological.

In contrast, the Mathematica script (lines 33-37 of the `.wl`) builds `chi_s` from the physical definition `m*c_s*L/hbar` and applies the healing-length substitution `c_s -> hbar/(2*m*ell)` followed by `L/ell -> Lambda_ell`, so its identical assertion at line 42 is non-tautological — it tests that the chain of substitutions actually lands on `Lambda_ell/2`. The SymPy script skips that chain entirely.

**Why this matters:**
The docstring's first non-trivial claim is "use the exact GNLS healing/compliance width `ell = hbar/(2 m c_s)`, verify `chi_s = Lambda_ell/2`". The SymPy script never introduces `hbar`, `m`, or `c_s` symbols and never performs the substitution; it asserts the conclusion against itself. If, hypothetically, the healing-length relation produced `chi_s = Lambda_ell/3` or `chi_s = 2*Lambda_ell`, the SymPy script would still print `PASS` because it has hard-coded `chi_lock = Lambda_ell/2` from the start. The reference-branch assertions inherit the same vacuity.

**Required change:**
Mirror the Mathematica derivation in the SymPy script. Introduce `hbar, m_psi, c_s, ell, L` as positive symbols, build `chi_def = m_psi*c_s*L/hbar`, apply `subs(c_s, hbar/(2*m_psi*ell))` to get an intermediate `chi_in_ell = L/(2*ell)`, then apply `subs(L/ell, Lambda_ell)` (or equivalently introduce `Lambda_ell = L/ell` and reduce the expression). The resulting `chi_lock` should equal `Lambda_ell/2` only after that substitution chain. Then keep the assertions on lines 39, 54 unchanged — they will now be non-tautological because `chi_lock` arrives from a real algebraic path rather than from a literal definition.

**Verification:**
After the fix, the SymPy script's intermediate prints should show `chi (from healing) = L/(2*ell)` and `chi_lock (after Lambda_ell substitution) = Lambda_ell/2`. The assertion at line 39 (`chi_s - Lambda_ell/2`) should still evaluate to 0 but now would fail if any step in the substitution were wrong. The output file's PASS line for `chi_s - Lambda_ell/2 = 0` remains, but is now preceded by the derived intermediate form.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:1-60`

**What's wrong:**
The docstring at lines 2-10 explicitly advertises four checks: (a) use the GNLS healing/compliance width `ell = hbar/(2 m c_s)`, (b) verify `chi_s = Lambda_ell/2`, (c) verify `kappa = (9/5) Lambda_ell^2` on the Family-1 branch, (d) evaluate the reference branch. The script never introduces `hbar`, `m`, or `c_s` as symbols anywhere — grep on the source shows zero occurrences. Claim (a) is therefore unverified entirely; the script jumps straight to the conclusion `Lambda_ell/2` on line 33 without performing any healing-length substitution.

Furthermore, the kappa relation `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2` is hard-coded as the literal expression on line 34 with no derivation or comment on where the `4` and `4/5` come from (Family-1 branch coefficients). The check at line 40 then reduces to the schoolbook identity `(1/2)^2 * 4 + 4/5 = 9/5`, which is correct arithmetic but tests no physics — it would PASS regardless of whether the underlying Family-1 EL equation actually produces those coefficients.

**Why this matters:**
The audit is a verification artifact. Verifying only an arithmetic identity, while skipping the physical substitution the docstring promises, gives false confidence that the GNLS healing-length lock was checked. A reader inspecting the saved output sees PASS lines for "chi_s - Lambda_ell/2" and assumes the GNLS chain was validated, when in fact only the trivial algebra `Lambda_ell/2 - Lambda_ell/2` was performed.

**Required change:**
In addition to the substitution chain mandated by F1, add an explicit symbolic statement of the Family-1 branch formula `kappa(chi_s, Lambda_ell) = 4 chi_s^2 + (4/5) Lambda_ell^2` as a named expression `kappa_family1`, then evaluate it at the derived `chi_lock = Lambda_ell/2` and assert it equals `(9/5)*Lambda_ell^2`. This preserves the current check at line 40 but documents the origin of the coefficients. A short inline comment above line 34 stating the source unit / branch for `4` and `4/5` is sufficient anchoring.

**Verification:**
After the fix, the script source contains the symbols `hbar`, `m_psi` (or `m`), `c_s`, and `ell`; the derivation `chi_def -> subs(c_s) -> subs(L/ell)` is visible; and a comment above the kappa definition cites the Family-1 EL coefficient source. The assertions on lines 39, 40, 54, 55 continue to PASS. The output file gains intermediate `chi (healing chain)` lines.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration; it is in fact the *more substantive* of the two. Compare:

- SymPy line 33: `chi_lock = sp.simplify(Lambda_ell / 2)` — pure definition.
- Mathematica lines 33-37:
  ```
  chiFromHealing = FullSimplify[(mpsi*cSw*len/hbar) /. cSw -> hbar/(2*mpsi*ell), ...]
  chiLock = FullSimplify[chiFromHealing /. (len/ell) -> lambdaEll, ...]
  ```
  builds `chi_s` from the physical definition `m c_s L / hbar`, substitutes the healing-length identity for `c_s`, then re-expresses `L/ell` as `Lambda_ell`.

The Mathematica script independently derives the result the SymPy script merely assumes. Both then compute `kappa = 4 chi^2 + (4/5) Lambda_ell^2` and check it equals `(9/5) Lambda_ell^2`, and both confirm the reference branch. No transliteration finding; the asymmetry runs the other way (SymPy under-implements).

## Engine cross-check

Both engines agree on all four numeric/symbolic outputs:

| Quantity | SymPy | Mathematica |
|---|---|---|
| chi_s (locked) | `Lambda_ell/2` | `lambdaEll/2` |
| kappa | `9*Lambda_ell**2/5` | `(9*lambdaEll^2)/5` |
| chi_ref (Lambda=37) | `37/2` | `37/2` |
| kappa_ref | `12321/5` | `12321/5` |
| alpha | `111*sqrt(5)/5` | `111/Sqrt[5]` (= `111*sqrt(5)/5`) |
| alpha numeric | `49.640709100495331260` | `49.64070910049533126028...` |

No `engine_disagreement` finding. Note that the alpha forms `111*sqrt(5)/5` and `111/Sqrt[5]` are algebraically identical (`111/sqrt(5) = 111*sqrt(5)/5`); the numeric values match to 20 digits.

## Verdict justification

The math itself holds up: `chi_s = Lambda_ell/2` and `kappa = (9/5) Lambda_ell^2` are correct given the GNLS healing-length lock and the Family-1 branch formula `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`. The Mathematica script verifies these non-tautologically. The verdict is `findings` (not `clean`) because the SymPy script is structurally hollow on the central claim — it defines `chi_lock = Lambda_ell/2` then asserts `chi_lock - Lambda_ell/2 == 0`, which is an identity by construction, and never introduces the `hbar / m / c_s / ell` symbols its own docstring promises to use. Engine-cross-check confirms no numerical disagreement, and outputs are fresh (SymPy: 2026-05-11, Mathematica: 2026-05-11, both newer than the corresponding script mtimes Apr 1 and May 11). No `stop_cold`: the fix is local to the SymPy script and does not alter any downstream-quoted constant.

## Self-test notes

I checked: (1) variable independence — proposed F1 fix introduces `hbar, m_psi, c_s, ell, L` as positive symbols and the substitution chain `c_s -> hbar/(2*m_psi*ell)` followed by `L/ell -> Lambda_ell` reduces `m_psi*c_s*L/hbar` to `Lambda_ell/2` correctly (manual algebra: `m * (hbar/(2 m ell)) * L / hbar = L/(2 ell) = Lambda_ell/2`); the resulting `chi_lock - Lambda_ell/2` is then a real zero rather than a definitional zero. (2) Symmetry/parity is not applicable — no integrals over symmetric domains here. (3) Trivial-case substitution: with `Lambda_ell = 37`, the derived `chi_lock` yields `37/2`, matching the existing reference assertion; `kappa = 4*(37/2)^2 + (4/5)*37^2 = 1369 + 5476/5 = 6845/5 + 5476/5 = 12321/5`, matching the existing reference assertion. (4) Path specifications n/a — no missing-script findings.
