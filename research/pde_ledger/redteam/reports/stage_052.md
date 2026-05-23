---
unit_id: 052
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 052 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.txt`

## What the script claims to verify

The scripts assert exact symbolic identities for the non-twin asymmetry threshold in stage 35/052: (a) the closed form `zeta_req = (Pi_tr - C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]`, including its boundary values at `Pi_tr = C_mix` (zero) and `Pi_tr = 2 C_mix` (unity); (b) the closed form for `d zeta_req/d Pi_tr` equal to `C_mix(1-eps)/[C_mix - eps(2 C_mix - Pi_tr)]^2`; (c) the excess `Delta_zeta = zeta_req - 1` matching `(1 - eps)(Pi_tr - 2 C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]`; (d) that `zeta_0^(phys) = K_W Omega_0^2 / K_phi0` reduces to `zeta_req` when either `Omega_0^2` is set to `zeta_req K_phi0/K_W` or `K_phi0` is set to `K_W Omega_0^2/zeta_req`; and (e) the softening fraction `1 - K_phi0^req/K_W` at equal overlap matches `(1-eps)(Pi_tr - 2 C_mix)/(Pi_tr - C_mix)`. The symmetric twin diagnostic `zeta_0^(twin) = K_W*1/K_W = 1` is also asserted.

## Assertion inventory

| #  | Script      | Line  | Form                                                                            | Anchored to claim? |
|----|-------------|-------|----------------------------------------------------------------------------------|--------------------|
| A1 | sympy       | 52    | `expect_zero(zeta_req at Pi_tr = C_mix, ...)`                                    | yes                |
| A2 | sympy       | 53    | `expect_zero(zeta_req at Pi_tr = 2 C_mix minus 1, ...)`                          | yes                |
| A3 | sympy       | 61    | `expect_zero(dzeta_req/dPi - expected, ...)`                                     | yes                |
| A4 | sympy       | 69    | `expect_zero(Delta_zeta - expected, ...)`                                        | yes                |
| A5 | sympy       | 85-88 | `expect_zero(threshold equality at fixed stiffness, ...)`                        | no (tautological)  |
| A6 | sympy       | 89-92 | `expect_zero(threshold equality at fixed overlap, ...)`                          | no (tautological)  |
| A7 | sympy       | 96    | `expect_zero(zeta_0^(twin) - 1, ...)`                                            | partial            |
| A8 | sympy       | 112   | `expect_zero(softening fraction - expected, ...)`                                | yes                |
| B1 | mathematica | 38    | `expectZero[zeta_req at Pi_tr = C_mix, ...]`                                     | yes                |
| B2 | mathematica | 39    | `expectZero[zeta_req at Pi_tr = 2 C_mix minus 1, ...]`                           | yes                |
| B3 | mathematica | 44    | `expectZero[dzeta_req/dPi - expected, ...]`                                      | yes                |
| B4 | mathematica | 49    | `expectZero[Delta_zeta - expected, ...]`                                         | yes                |
| B5 | mathematica | 58    | `expectZero[threshold equality at fixed stiffness, ...]`                         | no (tautological)  |
| B6 | mathematica | 59    | `expectZero[threshold equality at fixed overlap, ...]`                           | no (tautological)  |
| B7 | mathematica | 67    | `expectZero[zeta_0^(twin) - 1, ...]`                                             | partial            |
| B8 | mathematica | 71    | `expectZero[softening fraction - expected, ...]`                                 | yes                |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py:77-92`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl:52-59`

**What's wrong:**
The two "threshold equality" assertions are algebraic identities by construction. In SymPy:

```
Omega0_req_sq = sp.simplify(zeta_req * Kphi0 / KW)             # line 77
Kphi0_req    = sp.simplify(KW * Omega0**2 / zeta_req)           # line 78
...
expect_zero("threshold equality at fixed stiffness",
            zeta0_phys.subs(Omega0**2, Omega0_req_sq) - zeta_req)   # line 85-88
expect_zero("threshold equality at fixed overlap",
            zeta0_phys.subs(Kphi0, Kphi0_req) - zeta_req)            # line 89-92
```

`Omega0_req_sq` is *defined* as `zeta_req * Kphi0 / KW`, so substituting it into `zeta0_phys = KW*Omega0**2/Kphi0` necessarily yields `KW*(zeta_req*Kphi0/KW)/Kphi0 = zeta_req`. The residual is identically zero regardless of any physics; the assertion cannot fail. The Mathematica `.wl` lines 58-59 use the same construction with the same defect.

These two checks together with their associated lines `Omega0_req_sq = ...` and `Kphi0_req = ...` constitute the only verification of the "lowest-lane rescue thresholds" theorem in the ledger; replacing them with non-circular checks is the only way to substantiate those thresholds.

**Why this matters:**
The script's STAGE-35-THEOREM-LEDGER block (lines 114-125) advertises the rescue thresholds `Omega_0^2 >= zeta_req K_phi0 / K_W` and `K_phi0 <= K_W Omega_0^2 / zeta_req`. Currently nothing in the script proves these are equivalent to the support-ratio condition `zeta_0^(phys) = zeta_req`; the existing assertions just confirm the algebraic re-arrangement closes on itself. A genuine check is to solve `zeta_0^(phys) - zeta_req = 0` for the rescue variable and confirm the closed-form result matches.

**Required change:**
Replace lines 85-88 and 89-92 in the SymPy script and lines 58-59 in the Mathematica script with non-circular checks. Use `sp.solve` / `Solve` on `zeta0_phys - zeta_req` to derive `Omega0**2` (and `Kphi0`) and verify the resulting expression equals `zeta_req * Kphi0 / KW` (and `KW * Omega0**2 / zeta_req` respectively). The pre-existing `Omega0_req_sq`/`Kphi0_req` definitions stay; only the assertions that compare via substitution should change to comparisons of solved-for expressions against those definitions.

**Verification:**
After Codex applies, the .py and .wl output transcripts should contain two new lines reading
`solve(zeta_phys = zeta_req) for Omega0^2 - expected = 0` and `solve(zeta_phys = zeta_req) for Kphi0 - expected = 0`,
each printed as `PASS` and each with the prior tautological lines removed. The pipeline must still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl:33-71`

**What's wrong:**
The `.wl` script is a near-line-by-line transliteration of the `.py`. Compare:

SymPy lines 44-45, 60, 68, 111:
```
Sreq = sp.simplify(Pi_tr / Cmix)
zeta_req = sp.simplify((Sreq - 1) / (1 + eps * (Sreq - 2)))
dz_expected = sp.simplify(Cmix * (1 - eps) / (Cmix - eps * (2 * Cmix - Pi_tr)) ** 2)
Delta_expected = sp.simplify((1 - eps) * (Pi_tr - 2 * Cmix) / (Cmix - eps * (2 * Cmix - Pi_tr)))
soft_expected = sp.simplify((1 - eps) * (Pi_tr - 2 * Cmix) / (Pi_tr - Cmix))
```

Mathematica lines 33-34, 42, 47, 65:
```
sReq = FullSimplify[piTr/cMix, ...];
zetaReq = FullSimplify[(sReq - 1)/(1 + eps (sReq - 2)), ...];
dZExpected = FullSimplify[cMix (1 - eps)/(cMix - eps (2 cMix - piTr))^2, ...];
deltaExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(cMix - eps (2 cMix - piTr)), ...];
softExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(piTr - cMix), ...];
```

Variable names are simply camelCased, the construction order is identical, and the "expected" closed forms are hardcoded character-for-character. The Mathematica engine never re-derives any of the three closed forms (`dZExpected`, `deltaExpected`, `softExpected`) from the physical premise — it just types in the same answer the SymPy file does and lets `FullSimplify` confirm the difference is zero. Under the second-engine policy this is a paper-thin cross-check: both engines simplify the same identity to zero.

**Why this matters:**
The whole point of a second engine is to catch errors a single CAS can hide — wrong factor, wrong sign, wrong simplification under aggressive assumptions. If both engines start from the same hand-typed answer, any error in that answer survives. The current `.wl` does not exercise the engine's independent algebra.

**Required change:**
Re-derive the three "expected" closed forms inside the Mathematica script using independent Mathematica operations rather than transcribed strings. Specifically:

1. Build `zetaReq` via `Solve` on the support-ratio equation `(s - 1) - zeta (1 + eps (s - 2)) == 0`, with `s -> piTr/cMix`, returning `zetaReq` as the unique closed-form solution rather than typing `(sReq - 1)/(1 + eps (sReq - 2))`.
2. Construct `dZExpected` by `Together[D[zetaReq, piTr]]` and compare against the explicitly factored form via `FullSimplify[dZExpected - cMix (1 - eps)/(cMix - eps (2 cMix - piTr))^2] == 0`, but keep `dZdPi` derived from a *different* path — e.g., logarithmic differentiation: `dZdPiAlt = zetaReq (D[Log[Numerator[Together[zetaReq]]], piTr] - D[Log[Denominator[Together[zetaReq]]], piTr]) // FullSimplify` — and compare those two engine-derived results against each other before comparing to the hardcoded form. The cross-check is `dZdPi - dZdPiAlt == 0`.
3. Construct `deltaExpected` by `Together[zetaReq - 1]` directly and compare against the hardcoded form via a second `expectZero`. The independently-derived `Together` result is the engine's own algebra; the hardcoded form is what we are auditing against.
4. Construct `softExpected` by `Together[1 - 1/zetaReq]` rather than typing the right-hand side; then verify it equals the hardcoded `(1 - eps) (piTr - 2 cMix)/(piTr - cMix)`.

Each of these gives Mathematica an independent algebraic path so the cross-engine check actually carries information. The existing assertions (`expectZero[..., softFrac - softExpected]` etc.) can remain unchanged; what must change is how the `dZExpected`, `deltaExpected`, and `softExpected` (and `zetaReq`) are *produced* inside the `.wl`.

**Verification:**
After Codex applies, the `.wl` script should contain a `Solve` call deriving `zetaReq` from `(s-1) - zeta (1 + eps (s-2)) == 0`, a derivation of `dZdPiAlt` via logarithmic differentiation distinct from the direct `D[zetaReq, piTr]`, and `Together`-based constructions for `deltaExpected` and `softExpected`. The output transcript should print these derived forms before the existing `expectZero` calls and include a new `PASS: dZdPi vs dZdPiAlt` (or analogous) line. Exit 0 required.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py`. Both scripts hard-code the same three closed-form "expected" expressions (`dZExpected`, `deltaExpected`, `softExpected`) verbatim. Variable names are merely re-cased; the construction order is identical; no independent algebraic path is used. See F2 for line-by-line correspondence.

## Engine cross-check

Both engines simplify each residual to `0` and exit `0`. The outputs agree on:

- `zeta_req = -(Cmix - Pi_tr)/(Cmix - 2*Cmix*eps + eps*Pi_tr)` (SymPy out line 19-21; Mathematica out line 14).
- `d zeta_req / d Pi_tr = -Cmix*(eps-1)/(2*Cmix*eps - Cmix - Pi_tr*eps)^2` (SymPy out 26-29; Mathematica out 19).
- Delta_zeta numerator agrees; softening fraction agrees; threshold-equality residuals both `0` (but tautologically — see F1).

No `engine_disagreement` finding. Both transcripts are fresh (sympy out 2026-05-11T12:43:30; mathematica out 2026-05-11T12:52:07; script mtimes both earlier than these).

## Verdict justification

The substantive identities are real: the closed-form `zeta_req`, its derivative, the excess beyond the symmetric twin, and the softening fraction all check out under both engines, and I verified each by hand from the definition `zeta_req = (Sreq - 1)/(1 + eps(Sreq - 2))` with `Sreq = Pi_tr/Cmix`. What does not hold up is (a) the pair of "threshold equality" checks at lines 85-92 (.py) and 58-59 (.wl), which are algebraic identities by construction and cannot fail, and (b) the Mathematica script's transliteration of the hardcoded closed-form expecteds, which negates the value of a second engine. Attacks I tried that failed: substituting `eps=1` in the denominator (produces a defensible singular limit, not a wrong answer); checking parity/sign of `Delta_zeta` for `Pi_tr < 2 Cmix` (the script never claims sign, only the closed form); searching for missing branches in `Sreq <= 1` (the script's docstring is honest about which branch each formula applies to — no false generality claim).

## Self-test notes

- **Variable independence**: For the proposed F2 `D[zetaReq, piTr]` and `D[Log[Numerator[Together[zetaReq]]], piTr]`, `zetaReq` is `(piTr - cMix)/(cMix(1 - 2 eps) + eps piTr)`, which has explicit `piTr` dependence in both numerator and denominator, so both derivatives are nonzero. Confirmed.
- **Symmetry/parity**: No unbounded integrals in this unit; trap not applicable.
- **Trivial-case pre-check**: With `Pi_tr=3, Cmix=1, eps=1/2`, `zeta_req = (3-1)/(1 - (1/2)(2-3)) = 2/(3/2) = 4/3`. SymPy `sp.solve(KW*Omega0**2/Kphi0 - 4/3, Omega0**2)` returns `4*Kphi0/(3*KW) = zeta_req*Kphi0/KW`, matching the proposed F1 reformulation. The hardcoded `softExpected = (1-eps)(Pi_tr-2 Cmix)/(Pi_tr-Cmix) = (1/2)(3-2)/(3-1) = 1/4`, which matches `1 - 1/(4/3) = 1/4`. Confirmed.
- **Path specifications**: No new scripts proposed; existing paths are valid.
