---
unit_id: 010
batch: I.1
created_at: 2026-05-20T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-21T11:29:33-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 010

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl` (new file)

**Issue:**
Manifest entry `'010'` in `/var/projects/toy_physics/research/pde_ledger/redteam/MANIFEST.yaml` (lines 263-289) has `is_checkpoint: false` and `is_status_only_candidate: false`. The `mathematica` slot is `path: null, exists: false`. The current SymPy script verifies non-trivial algebraic identities (compatibility first variations, z0 cancellation, transported-target structure, Gaunt overlap ratios, trace anomaly grouping, primitive static Xi formula) without a second-engine cross-check. The second-engine policy requires an independent re-derivation, not a transliteration. Subtype: `missing_mathematica`.

**Required change:**
Create the file at the absolute path above. Structure it as standard Mathematica/Wolfram Language (`.wl`), with `Print[...]` for the banner and `If[FullSimplify[lhs - rhs] =!= 0, Print["FAIL: ", label]; Exit[1]]` for each assertion. The script must independently verify each manifest item below; intermediate variable naming and derivation order MUST differ from the SymPy script (do not, e.g., name the perturbed variables `D0p, D2p, ...` and do not solve the K equation in the same form). Specifically:
- Use Mathematica-natural constructs: `Series[..., {eps, 0, 1}]` and `Coefficient[..., eps, 1]` instead of `sp.diff(...).subs(eps, 0)`.
- Use `ThreeJSymbol` or `ClebschGordan` (composed into the standard Gaunt formula) rather than calling a Gaunt routine directly.
- Express the bundle slots P0, P2, P4 in their closed form first, then take the eps-coefficient, rather than building `P0p, P2p, P4p` as intermediate symbols.
- Build the K-elimination compatibility from the relation `(K - B0 - Z0slot - eps z0) (T + eps z4) == 3 (S + eps z2)^2` solved for K with `Solve` (Mathematica), then independently re-solve the normalization equation `(N0 + eps n0)/(K - B0 - Z0slot - eps z0) == Ptarget` for K. After both are obtained, take their difference and its eps-coefficient as `dcompat`.

**Claim manifest:**
The new script must contain assertions (each in the form `If[FullSimplify[expr] =!= 0, Print["FAIL: M<n>"]; Exit[1]]`) verifying each of the following identities independently. The symbols are the same as in the SymPy script; the algebra need not look the same.

- M1: `dP0 = n0/D0 + N0 z0/D0^2` where `P0(eps) = (N0 + eps n0)/(D0 - eps z0)` and `dP0 = D[P0(eps), eps] /. eps -> 0`.
- M2: `dP2 = n2/D0 + 2 N2 z0/D0^2 + 2 N0 z2/D0^2 - 2 D2 n0/D0^2 - 4 D2 N0 z0/D0^3` where `P2(eps) = ((D0 - eps z0)(N2 + eps n2) - 2(D2 - eps z2)(N0 + eps n0))/(D0 - eps z0)^2`. (Codex: verify this hand-derivation by independently expanding via `Series[P2(eps), {eps, 0, 1}]`. If it differs, M2 is whatever Mathematica's expansion gives — but Codex must NOT mark the assertion satisfied unless both engines agree on the closed form, and must escalate via a `## Blocked: F1` note if there's a discrepancy with the parallel SymPy fix from F2.)
- M3: `dP4 = ` first-order eps coefficient of `P4(eps) = ((D0 - eps z0)^2 (N4 + eps n4) - 2(D0 - eps z0)((D2 - eps z2)(N2 + eps n2) + (D4 - eps z4)(N0 + eps n0)) + 3(D2 - eps z2)^2 (N0 + eps n0))/(D0 - eps z0)^3`. Codex must derive the closed form via `Series[P4(eps), {eps, 0, 1}]` and use that as the M3 right-hand side; the same cross-check rule as M2 applies.
- M4: `K_one_pole(eps) = B0 + Z0slot + eps z0 + 3 (S + eps z2)^2/(T + eps z4)` is the unique solution of `(K - B0 - Z0slot - eps z0)(T + eps z4) = 3 (S + eps z2)^2`.
- M5: `dK_one_pole = z0 + 6 S z2/T - 3 S^2 z4/T^2`.
- M6: `K_norm(eps) = B0 + Z0slot + eps z0 + (N0 + eps n0)/Ptarget`, the unique solution of `(N0 + eps n0)/(K - B0 - Z0slot - eps z0) = Ptarget`.
- M7: `dK_norm = z0 + n0/Ptarget`.
- M8: With `compat(eps) = K_norm(eps) - K_one_pole(eps)` and `compat_direct(eps) = (N0 + eps n0)/Ptarget - 3 (S + eps z2)^2/(T + eps z4)`, the identity `compat(eps) - compat_direct(eps) = 0` holds at all orders in eps.
- M9: `dcompat = n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2` (note the `+` on the z4 term).
- M10: With `Ptarget_transport(eps) = (N0 + eps n0)/D0target`, the K-equation solved for K gives `K_norm_transport(eps) = B0 + Z0slot + eps z0 + D0target`.
- M11: `compat_transport(eps) := K_norm_transport(eps) - K_one_pole(eps) = D0target - 3 (S + eps z2)^2/(T + eps z4)`. In particular the z0 dependence cancels.
- M12: `dcompat_transport = -6 S z2/T + 3 S^2 z4/T^2`.
- M13: With `lam_{2m} = (-1)^m * Gaunt(2,2,2,0,m,-m)/Gaunt(2,2,2,0,0,0)` for the appropriate Gaunt convention, `lam_{20} = 1`, `lam_{21} = 1/2`, `lam_{22} = -1`. Also `Gaunt(2,2,2,0,m,m) = 0` for `m = 1, 2`.
- M14: With `x_{2m} = x0 + eps * lam_{2m} * x1` and `(xbar, ax, bx) = ((x_{20} + 2 x_{21} + 2 x_{22})/5, (2 x_{20} - x_{21} - x_{22})/10, (x_{21} - x_{22})/2)`: `xbar = x0`, `ax = eps x1/4`, `bx = 3 eps x1/4`, `bx = 3 ax`.
- M15: With `z0_prim = (Delta q1 - Q d1)/Delta^2`, `n0_prim = 2 P (Delta p1 - P d1)/Delta^3`, and the substitution `N0sym -> P^2/Delta^2`, the combination `n0_prim/N0sym + z0_prim/D0sym` equals `2 p1/P - 2 d1/Delta + (Delta q1 - Q d1)/(D0sym Delta^2)`.

In addition, M16/M17 are negative-control mutations matching the existing SymPy `assert_nonzero` guards: assert that `dcompat_direct - (n0/Ptarget - 6 S z2/T - 3 S^2 z4/T^2)` is NOT zero (it should equal `6 S^2 z4/T^2`), and `dcompat_transport - (-6 S z2/T - 3 S^2 z4/T^2)` is NOT zero. In Mathematica, encode these as `If[FullSimplify[expr] === 0, Print["FAIL: M16 mutation passed unexpectedly"]; Exit[1]]`.

Finish with `Print["STATUS: PASS"]` and exit 0.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 010` and confirm (a) the new file exists at the named path, (b) it exits 0, (c) its assertions cover M1-M17, and (d) variable naming and ordering are sufficiently different from the SymPy script to not constitute a transliteration.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl`
- summary: Created an independent Mathematica audit covering M1-M17 with native series coefficients, Solve-based K surfaces, ThreeJSymbol Gaunt overlaps, and negative-control mutations.
- deviation: Created the `.wl` under `mathematica/` per `.redteam-config.yaml` runner layout instead of the finding's `scripts/` target.

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:59-62`

**Issue:**
Lines 59-62 verify only that `dP2` and `dP4` contain the expected free symbols, not that they equal any specific closed form. This is much weaker than the line-58 check on `dP0`, and lets sign errors, factor-of-2 errors, or wrong-power errors slip through. The script's banner and master-note title both advertise verification of the bundle perturbation slots P0/P2/P4 together.

**Required change:**
Replace lines 59-62 with explicit closed-form `assert_zero` calls in the same style as line 58. Specifically:

Before (lines 58-62):
```python
    assert_zero("delta P0", dP0 - (n0 / D0 + N0 * z0 / D0**2))
    if not {z0, z2, n0, n2}.issubset(dP2.free_symbols):
        raise AssertionError("delta P2 is missing one of the advertised first-order slots")
    if not {z0, z2, z4, n0, n2, n4}.issubset(dP4.free_symbols):
        raise AssertionError("delta P4 is missing one of the advertised first-order slots")
```

After:
```python
    assert_zero("delta P0", dP0 - (n0 / D0 + N0 * z0 / D0**2))
    assert_zero(
        "delta P2",
        dP2 - (
            n2 / D0
            + 2 * N2 * z0 / D0**2
            + 2 * N0 * z2 / D0**2
            - 2 * D2 * n0 / D0**2
            - 4 * D2 * N0 * z0 / D0**3
        ),
    )
    assert_zero(
        "delta P4",
        dP4 - (
            n4 / D0
            + 3 * N4 * z0 / D0**2
            - 2 * (D2 * n2 + N2 * z2 + D4 * n0 + N0 * z4) / D0**2
            - 4 * (D2 * N2 + D4 * N0) * z0 / D0**3
            + 6 * D2 * N0 * z2 / D0**3
            + 3 * D2**2 * n0 / D0**3
            + 9 * D2**2 * N0 * z0 / D0**4
        ),
    )
```

Important: before committing, run the closed forms through a hand-derivation. The dP2 closed form follows from
`P2p = (D0p N2p - 2 D2p N0p)/D0p**2` with `D0p = D0 - eps z0`, `N2p = N2 + eps n2`, `D2p = D2 - eps z2`, `N0p = N0 + eps n0`:
- Numerator d/deps at eps=0: `(-z0) N2 + D0 n2 - 2 (-z2) N0 - 2 D2 n0 = -z0 N2 + D0 n2 + 2 z2 N0 - 2 D2 n0`
- Denominator at eps=0: `D0^2`
- Quotient-rule second term: `(D0 N2 - 2 D2 N0) * (-2 D0 z0)/(D0^2)^2 = -2(D0 N2 - 2 D2 N0) z0/D0^3`
- Sum: `n2/D0 - z0 N2/D0^2 + 2 z2 N0/D0^2 - 2 D2 n0/D0^2 - 2(D0 N2 - 2 D2 N0) z0/D0^3 = n2/D0 + 2 N2 z0/D0^2 + 2 N0 z2/D0^2 - 2 D2 n0/D0^2 - 4 D2 N0 z0/D0^3`. (The `-z0 N2/D0^2` and `-2 D0 N2 z0/D0^3 = -2 N2 z0/D0^2` combine to give `-3 N2 z0/D0^2`. Wait — let me redo this.)

Codex action: if the hand-derivation above yields a different result, USE the result Codex derives by `dP2_target = sp.simplify(sp.diff(P2p, eps).subs(eps, 0))` mentally — that is what `dP2` already equals. The point is the right-hand side of the `assert_zero` must equal `dP2`, with `D0, D2, N0, N2, z0, z2, n0, n2` as the only free symbols. **Codex: if you cannot determine the exact closed form by hand without running Python, append `## Blocked: F2` with the question "Closed form of dP2 and dP4 — please supply" instead of guessing.**

Same procedure for `dP4`. The proposed closed form above is a derivation hint, not authoritative — if Codex's hand-derivation disagrees with it, trust the hand-derivation. If unsure, BLOCK rather than guess.

**Verification:**
After Codex applies, lines 58-62 should contain three `assert_zero` calls of the same form, no `if not {...}.issubset(...)` blocks. `redteam exec-sympy 010` must exit 0; if `assert_zero` raises on the new dP2 or dP4 form, the closed-form right-hand side is wrong and Codex must correct it before declaring the finding applied.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py`
- summary: Replaced the dP2 and dP4 symbol-presence checks with explicit closed-form `assert_zero` residual checks.
- deviation: Used expansion-derived dP2/dP4 targets because the illustrative forms in the directive disagreed with the perturbation definitions.
