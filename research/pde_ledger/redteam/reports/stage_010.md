---
unit_id: 010
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 010 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.txt`
- mathematica output: `(missing)`

Manifest entry confirms `is_status_only_candidate: false` and `is_checkpoint: false` (MANIFEST.yaml lines 270-271), so both engines are required.

## What the script claims to verify

The script targets the "step_08_projected_maxwell_push_bundle_master_notes" content. Concretely the assertions exercise four claim clusters: (1) first-order eps-expansion of the bundle P-slot ratios (P0, P2, P4) for a perturbed denominator/numerator pair, with an explicit closed form only for `dP0` and structural symbol-presence checks for `dP2`/`dP4`; (2) two ways of constructing the "one-pole compatibility" surface (eliminating K vs. competing K surfaces) and verifying their first variation in eps agrees, with z0 cancellation in the compatibility relation and z0 retention in the normalization K; (3) the analogous "transported-target" compatibility surface and its first variation, again showing z0 cancellation; (4) Gaunt-coefficient ratios for real Y20-square lanes (m=0,1,2) giving (1, 1/2, -1) and the resulting trace/anomaly decomposition (xbar, ax, bx) satisfying b = 3a along the weak-axisymmetric line; and (5) a primitive static "Xi" formula relating two ratio-perturbations under a specific N0sym substitution. Two `assert_nonzero` "mutation" guards intentionally break a sign and check the residual is nonzero, confirming the prior sign-correct assertions are not vacuous.

## Assertion inventory

| #  | Script | Line  | Form                                                                                                    | Anchored to claim? |
|----|--------|-------|---------------------------------------------------------------------------------------------------------|--------------------|
| A1 | sympy  | 58    | `assert_zero("delta P0", dP0 - (n0/D0 + N0*z0/D0**2))`                                                  | yes                |
| A2 | sympy  | 59-60 | `if not {z0,z2,n0,n2}.issubset(dP2.free_symbols): raise`                                                | partial            |
| A3 | sympy  | 61-62 | `if not {z0,z2,z4,n0,n2,n4}.issubset(dP4.free_symbols): raise`                                          | partial            |
| A4 | sympy  | 93    | `assert_zero("compatibility surface after eliminating K", (K_norm_p - K_one_pole_p) - compat_direct_p)` | yes                |
| A5 | sympy  | 94    | `assert_zero("one-pole K shift", dK_one_pole - (z0 + 6*S*z2/T - 3*S**2*z4/T**2))`                       | yes                |
| A6 | sympy  | 95    | `assert_zero("normalization K shift", dK_norm - (z0 + n0/Ptarget))`                                     | yes                |
| A7 | sympy  | 96    | `assert_zero("compatibility first variation from competing K surfaces", dcompat - (dK_norm - dK_one_pole))` | yes            |
| A8 | sympy  | 97    | `assert_zero("compatibility first variation from eliminated surface", dcompat - dcompat_direct)`        | yes                |
| A9 | sympy  | 98    | `assert_zero("compatibility first variation", dcompat_direct - (n0/Ptarget - 6*S*z2/T + 3*S**2*z4/T**2))` | yes              |
| A10| sympy  | 99    | `assert_zero("transported-target normalization K surface", K_norm_transport_p - (B0+Z0slot+eps*z0+D0target))` | yes           |
| A11| sympy  | 100-103| `assert_zero("transported-target compatibility surface", compat_transport_p - (D0target - 3*(S+eps*z2)**2/(T+eps*z4)))` | yes    |
| A12| sympy  | 104-107| `assert_zero("transported-target compatibility first variation", dcompat_transport - (-6*S*z2/T + 3*S**2*z4/T**2))` | yes        |
| A13| sympy  | 108   | `assert_zero("fixed-target compatibility z0 cancellation", sp.diff(dcompat_direct, z0))`                | yes                |
| A14| sympy  | 109   | `assert_zero("transported-target compatibility z0 cancellation", sp.diff(dcompat_transport, z0))`       | yes                |
| A15| sympy  | 110   | `assert_nonzero("normalization K surface still carries z0", sp.diff(dK_norm, z0) - 0)`                  | yes (weak)         |
| A16| sympy  | 111-114| `assert_nonzero("mutated compatibility transport should fail", dcompat_direct - (n0/Ptarget - 6*S*z2/T - 3*S**2*z4/T**2))` | yes  |
| A17| sympy  | 115-118| `assert_nonzero("mutated transported-target compatibility should fail", dcompat_transport - (-6*S*z2/T - 3*S**2*z4/T**2))` | yes  |
| A18| sympy  | 124   | `assert_zero("Y20 overlap lane 20", lam20 - 1)`                                                         | yes                |
| A19| sympy  | 125   | `assert_zero("Y20 overlap lane 21", lam21 - 1/2)`                                                       | yes                |
| A20| sympy  | 126   | `assert_zero("Y20 overlap lane 22", lam22 + 1)`                                                         | yes                |
| A21| sympy  | 131   | `assert_zero("weak-axisymmetric trace", xbar - x0)`                                                     | yes                |
| A22| sympy  | 132   | `assert_zero("weak-axisymmetric a anomaly", ax - eps*x1/4)`                                             | yes                |
| A23| sympy  | 133   | `assert_zero("weak-axisymmetric b anomaly", bx - 3*eps*x1/4)`                                           | yes                |
| A24| sympy  | 134   | `assert_zero("weak-axisymmetric line b=3a", bx - 3*ax)`                                                 | yes                |
| A25| sympy  | 144-148| `assert_zero("primitive static Xi formula", Xi_static.subs(...) - (2*p1/P - 2*d1/Delta + ...))`        | yes                |

The "partial" rows (A2, A3) flag a `insufficient_verification` candidate: the symbol-set check is far weaker than the explicit closed-form check used for dP0.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/` (no `*_stage010_*_mathematica_audit.wl` present)

**What's wrong:**
The manifest at `/var/projects/toy_physics/research/pde_ledger/redteam/MANIFEST.yaml` lines 263-289 records this unit with `is_checkpoint: false` and `is_status_only_candidate: false`, and the `mathematica` entry has `path: null, exists: false`. The SymPy script ships every nontrivial algebraic identity (compatibility first variations, z0 cancellation, transported-target structure, Gaunt-coefficient overlap ratios, trace/anomaly grouping, primitive static Xi formula) without a second-engine cross-check. Subtype: `missing_mathematica`.

**Why this matters:**
The unit's main asserts are non-trivial closed-form identities involving rational functions of D0, N0, S, T, Ptarget, the eps perturbations z_i / n_i, and Gaunt-coefficient overlaps. A solo SymPy run cannot catch (a) a sign error consistent between premise and asserted right-hand side (the author wrote both), (b) a `sp.simplify` bug specific to SymPy's normalization of nested rationals, or (c) a mis-encoding of the Y20-Y2m-Y2(-m) overlap convention. The second-engine policy exists precisely to provide an independent algebraic witness for these. Without a `.wl`, the audit cannot certify the assertions independently — the manifest's required-engines policy is violated.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl`. The script must independently re-derive (not transliterate) the same physical claims and exit nonzero on any failure. See claim manifest in the directive.

**Verification:**
Verifier runs `redteam exec-mathematica 010`. The new `.wl` must exit code 0 and its assertions must mirror the claim manifest (M1-M14 in the directive). Verifier also checks the script is not a line-by-line transliteration of the SymPy script (different intermediate variable names and a different derivation order are expected).

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:59-62`

**What's wrong:**
Lines 54-56 compute `dP0`, `dP2`, `dP4` — the first-order eps coefficients of the three projected push slots P0, P2, P4. Only `dP0` is then anchored to a closed form (line 58: `assert_zero("delta P0", dP0 - (n0/D0 + N0*z0/D0**2))`). For `dP2` and `dP4` the script only checks symbol-presence:

```python
    if not {z0, z2, n0, n2}.issubset(dP2.free_symbols):
        raise AssertionError("delta P2 is missing one of the advertised first-order slots")
    if not {z0, z2, z4, n0, n2, n4}.issubset(dP4.free_symbols):
        raise AssertionError("delta P4 is missing one of the advertised first-order slots")
```

A free-symbol test cannot detect a sign error, a factor-of-2 error, a wrong numerator power, or a miscoupling of the slot perturbations. Any expression with the right symbols (even physically wrong) passes. The script's banner ("Checked bundle perturbation slots") and the master-note title both advertise verification of the bundle perturbation slots P0/P2/P4 together; the actual verification only covers P0.

**Why this matters:**
The compatibility-surface and K-shift assertions further downstream (lines 93-109) depend on the P2- and P4-pole structure of the bundle via the symbols `D2, D4, N2, N4` and via the eliminated K-equation `(K - B0 - Z0slot - eps*z0)*(T + eps*z4) = 3*(S + eps*z2)**2` (line 66-67), where the `T`, `S`, and `z2`, `z4` slots correspond to the P2/P4 coefficients. A wrong `dP2` or `dP4` would not be caught here, leaving the upstream "bundle perturbation slot" claim only structurally verified.

**Required change:**
At `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:59-62`, replace the two `if not {...}.issubset(...)` blocks with explicit closed-form `assert_zero` calls of the same style as the existing `dP0` check (line 58). Derive the closed forms by hand from the definitions at lines 51-52:

- `P2p = (D0p*N2p - 2*D2p*N0p)/D0p**2`, so `dP2 = d/deps[P2p]|_{eps=0}` equals
  `n2/D0 + 2*N2*z0/D0**2 + 2*z2*N0/D0**2 - 4*D2*N0*z0/D0**3 - 2*D2*n0/D0**2`
  (verify by hand and via the script before committing).
- `P4p = (D0p**2*N4p - 2*D0p*(D2p*N2p + D4p*N0p) + 3*D2p**2*N0p)/D0p**3`, so `dP4` should be derived analogously.

Add:
```
assert_zero("delta P2", dP2 - <closed_form_dP2>)
assert_zero("delta P4", dP4 - <closed_form_dP4>)
```
and delete the two `if not {...}.issubset(...)` blocks. If the closed form differs from what you derive on paper, the existing assertion was wrong and `dP2`/`dP4` itself needs re-derivation in the script.

**Verification:**
After the fix, `assert_zero("delta P2", ...)` and `assert_zero("delta P4", ...)` should appear at the previously-occupied lines 59-62, and `redteam exec-sympy 010` must exit 0. The script's printed banner stays the same.

## Independent-derivation check (Mathematica)

Not applicable — no Mathematica script exists for this unit. See F1.

## Engine cross-check

Not applicable — single engine. The SymPy output transcript at `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.txt` (mtime 2026-05-11, exit code 0, banner "STATUS: PASS") is fresh relative to the script (script mtime 2026-05-04), so no `stale_output` finding.

## Verdict justification

Verdict is `findings`, not `clean`, because of F1 (missing Mathematica engine, required by the manifest's non-status, non-checkpoint flags) and F2 (the script's stated "bundle perturbation slots" coverage is honored only for the P0 slot — the P2 and P4 slots only get symbol-presence stubs at lines 59-62). The math of every `assert_zero` that I checked by hand (dP0, K shifts, compatibility surfaces in both fixed and transported targets, z0 cancellation, Gaunt-overlap ratios, weak-axisymmetric grouped trace, primitive static Xi) holds up: I tried to break sign and factor-of-2 attacks on the K-shift, the compatibility first-variation (note the `+3*S**2*z4/T**2` sign in line 98 is intentionally opposite to the `-3*S**2*z4/T**2` mutated version at line 113), and the trace-decomposition coefficients (xbar coefficient sum 1+2+2=5 across the lanes; ax, bx weights 2,-1,-1 and 0,1,-1), and all of them resolve correctly under the lane multipliers (1, 1/2, -1). The two `assert_nonzero` "mutation" guards (lines 111-114 and 115-118) protect against trivial passes of the sign-flipped variants. Neither finding triggers a stop-cold verdict: F1 is the standard second-engine gap, and F2 only tightens an under-verified portion of the same unit (no downstream sign or constant depends on whether the symbol-presence stub is replaced by a closed-form assertion).
