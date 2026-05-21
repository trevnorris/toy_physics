---
unit_id: 012
batch: I.1
created_at: 2026-05-20T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T11:40:52-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 012

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl` (new file)

**Issue:** Unit 012 is not status-only and not a checkpoint carve-out. Only the SymPy script exists. A second-engine Mathematica script is required so the engine cross-check is non-trivially satisfied. The new `.wl` must independently re-derive the claims, not transliterate the `.py`.

**Required change:**
Create the new `.wl` file with the structure below. The script must terminate with `Print["STATUS: PASS"]` only after all claim checks pass; each individual check must call `Exit[1]` on failure. Use Mathematica primitives `Series`, `Coefficient`, `D`, `Solve`, `Simplify`, `FullSimplify` rather than mirroring the SymPy choreography. In particular, do NOT define a `dlin` helper that mirrors `sp.series(...).removeO()` — instead compute first-order corrections directly via `Coefficient[Series[expr /. subs, {ell, 0, 1}], ell, 1]` or via `D[expr /. subs, ell] /. ell -> 0`, and additionally derive `z0..n4` from `D` of `Z0..N4` with respect to each primitive and contraction with its slippage (so the cross-check between the series route and the partial-derivative route is meaningful when compared across engines). Use distinct local variable names to avoid line-by-line correspondence.

Claim manifest (each claim must have at least one explicit `If[..., Exit[1]]` check):

M1. Primitive one-port formulas (typed as inputs, same closed forms as SymPy):
- `Z0 = Q/Delta`
- `Z2 = (Q S2 - H Delta)/Delta^2`
- `Z4 = (Q (S2^2 - Delta) - S2 H Delta)/Delta^3`
- `N0 = P^2/Delta^2`
- `N2 = 2 P (P S2 - Delta Gw)/Delta^3`
- `N4 = (Delta^2 Gw^2 - 2 Delta P^2 - 4 Delta P S2 Gw + 3 P^2 S2^2)/Delta^4`

M2. First-order corrections (must equal these closed forms after the substitutions `Q -> Q + ell q1`, etc., to first order in `ell`):
- `z0 = (Delta q1 - Q d1)/Delta^2`
- `z2 = (-Delta^2 h1 + Delta (H d1 + Q s1 + S2 q1) - 2 Q S2 d1)/Delta^3`
- `z4 = (-Delta^2 H s1 - Delta^2 S2 h1 - Delta^2 q1 + 2 Delta H S2 d1 + 2 Delta Q S2 s1 + 2 Delta Q d1 + Delta S2^2 q1 - 3 Q S2^2 d1)/Delta^4`
- `n0 = 2 P (Delta p1 - P d1)/Delta^3`
- `n2 = -(2 Delta^2 (Gw p1 + P g1) - 2 Delta P (2 Gw d1 + P s1 + 2 S2 p1) + 6 P^2 S2 d1)/Delta^4`
- `n4 = 2 (Delta^3 Gw g1 - Delta^2 Gw^2 d1 - 2 Delta^2 Gw P s1 - 2 Delta^2 Gw S2 p1 - 2 Delta^2 P S2 g1 - 2 Delta^2 P p1 + 6 Delta Gw P S2 d1 + 3 Delta P^2 S2 s1 + 3 Delta P^2 d1 + 3 Delta P S2^2 p1 - 6 P^2 S2^2 d1)/Delta^5`

M3. Static Xi1 closed form: `n0/N0 + z0/D0 == 2 p1/P - 2 d1/Delta + (Delta q1 - Q d1)/(D0 Delta^2)`.

M4. Compatibility surface from K-elimination (fixed Ptarget): after solving
`(K - B0 - (Z0slot + ell z0))(T + ell z4) == 3 (S + ell z2)^2` and
`(N0base + ell n0)/(K - B0 - (Z0slot + ell z0)) == Ptarget` for `K`,
the difference of solutions equals `(N0base + ell n0)/Ptarget - 3 (S + ell z2)^2/(T + ell z4)`. Verify by *substituting each solution back into its defining equation* (round-trip) and asserting the residual is zero — not by subtracting the two K's and renaming.

M5. Linear compatibility shift (fixed Ptarget): the coefficient of `ell^1` in M4's compatibility surface equals `n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2`.

M6. Transported-target normalization K surface: with `Ptarget_transport = (N0base + ell n0)/D0target`, solving the normalization equation for K gives `K = B0 + Z0slot + ell z0 + D0target`. Verify by round-trip into the original normalization equation.

M7. Transported-target compatibility surface and linear shift: `compat_transport = D0target - 3 (S + ell z2)^2/(T + ell z4)`, and the `ell^1` coefficient equals `-6 S z2/T + 3 S^2 z4/T^2`.

M8. z0 cancellation in compatibility: after substituting the derived `z0 = (Delta q1 - Q d1)/Delta^2` into both `K_norm` and `K_one` (fixed Ptarget and transported variants), the `q1` and `d1` partial derivatives of the difference receive contributions ONLY through `n0, z2, z4` and the explicit denominator `Delta` slippage in those quantities — there is no surviving `z0` channel. Equivalently, `D[K_norm - K_one, q1]` after the substitution must equal `D[compat_direct, q1]` (and similarly for `d1`), where `compat_direct` is the compatibility surface from M4. Conversely, `D[K_norm, q1]` and `D[K_norm, d1]` (without subtraction) must be nonzero — i.e. the normalization surface alone still carries the `z0` channel.

M9. Mutation sanity: flipping the sign on the `z4` term in either compatibility-shift closed form (M5 or M7) must make the residual nonzero.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 012` and confirms the new `.wl` exists and exits 0, and that each of M1-M9 produces at least one PASS-style print line.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl`
- summary: Created the Mathematica audit with independent primitive, correction, compatibility, z0-channel, and mutation checks for M1-M9.
- deviation: Placed the `.wl` under `mathematica/` per `.redteam-config.yaml` and the missing-script layout rule rather than the directive's `scripts/` target path.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:124-125`

**Issue:** The assertion at line 124, `assert_zero("primitive compatibility surface after eliminating K", (K_norm_p - K_one_p) - compat_direct_p)`, is algebraically guaranteed because `compat_direct_p` at line 111 is defined as exactly the difference of the closed forms that `K_norm_p` (line 107-110) and `K_one_p` (line 103-106) take after sympy's `solve` returns. Line 125's `dCompat - dCompat_direct == 0` is a direct corollary (series expansion is linear). Neither assertion can fail under any physically-relevant perturbation.

**Required change:**

1. Replace line 124 with a pair of *round-trip* residual checks that exercise sympy's `solve` against the original defining equations. Concretely, after line 119 (the `compat_transport = K_norm_transport_p.subs(ell, 0) - K_one_p.subs(ell, 0)` line — wait, re-read; actually after line 117, where `compat_transport_p` is defined), and before the previous line 124 assertion, add:

```
assert_zero(
    "K_one solve round-trip",
    (K_one_p - B0 - (Z0slot + ell * z0)) * (T + ell * z4) - 3 * (S + ell * z2) ** 2,
)
assert_zero(
    "K_norm solve round-trip",
    (N0base + ell * n0) / (K_norm_p - B0 - (Z0slot + ell * z0)) - Ptarget,
)
```

Then DELETE the line 124 assertion (`assert_zero("primitive compatibility surface after eliminating K", ...)`) entirely. The substantive check (A10) at line 126-129 remains and now stands on top of round-trip-verified K solutions.

2. DELETE the line 125 assertion (`assert_zero("primitive compatibility shift from eliminated surface", dCompat - dCompat_direct)`). It is a corollary of the deleted line 124.

**Verification:**
After patching, the output transcript no longer prints `primitive compatibility surface after eliminating K` or `primitive compatibility shift from eliminated surface` labels. Instead it prints `K_one solve round-trip` and `K_norm solve round-trip`. Re-running the script still exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`
- summary: Replaced the eliminated-surface tautology with round-trip checks for the solved one-pole and normalization K equations.
- deviation: none

## F3 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:139-148`

**Issue:** The block declares a fresh symbol `z0_probe` and types `+ ell*z0_probe` identically into `K_norm_probe`, `K_one_probe`, and `K_norm_transport_probe`. The assertions at lines 143 and 144-147 that `z0_probe` cancels from differences are therefore guaranteed by the typed-in form alone. The `assert_nonzero` at line 148 reduces to `ell != 0`. None of these tests the derived bundle correction `z0` from line 76 or the q1/d1 channel the prose ascribes to it.

**Required change:**

Replace lines 139-148 with checks that use the *derived* `z0 = (Delta*q1 - Q*d1)/Delta**2` substituted into the actual `K_norm_p` and `K_one_p` from lines 103-110 (and their transport counterpart from line 113-116). Concretely, replace those ten lines with:

```
    # The derived z0 = (Delta*q1 - Q*d1)/Delta**2 enters K_norm_p and K_one_p
    # through identical +ell*z0 terms, so it must cancel from their difference.
    # The substantive check is that the q1 and d1 partials of the compatibility
    # difference receive NO contribution from the z0 channel: they must equal
    # the q1 and d1 partials of compat_direct_p, which depends only on n0, z2, z4.
    for slip in (q1, d1):
        assert_zero(
            f"primitive fixed-target compatibility has no z0 channel in {slip}",
            sp.diff(K_norm_p - K_one_p, slip) - sp.diff(compat_direct_p, slip),
        )
        assert_zero(
            f"primitive transported-target compatibility has no z0 channel in {slip}",
            sp.diff(K_norm_transport_p - K_one_p, slip) - sp.diff(compat_transport_p, slip),
        )
    # Conversely, the normalization K surface alone DOES retain a q1/d1 channel
    # through the derived z0.
    assert_nonzero(
        "primitive normalization K surface retains q1 channel from z0",
        sp.diff(K_norm_p, q1),
    )
    assert_nonzero(
        "primitive normalization K surface retains d1 channel from z0",
        sp.diff(K_norm_p, d1),
    )
```

This uses the *derived* z0 (already substituted into `K_norm_p` and `K_one_p` via `ell * z0` at lines 103-110 / 113-116), not a fresh probe symbol. The first four checks verify the cancellation claim against the actual compatibility surface; the last two verify the "but normalization still carries z0" claim.

**Verification:**
After patching, the output transcript prints labels containing `primitive fixed-target compatibility has no z0 channel in q1`, `... in d1`, the two transported-target counterparts, and the two `retains ... channel from z0` lines, all passing. The `z0_probe` symbol no longer appears in the script.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`
- summary: Replaced the probe-symbol cancellation block with q1/d1 derivative checks against the derived z0-bearing K surfaces and direct compatibility surfaces.
- deviation: none

## F4 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:81-93`

**Issue:** Assertions A1-A6 (lines 88-93) compare `dlin(expr)` against `frechet(expr)` where both compute the first-order Taylor coefficient of the same `expr` by two sympy-internal routes. They cannot detect a typo in the typed-in closed forms `Z0..N4` at lines 56-62 — both routes would silently linearize whichever (possibly wrong) primitive was typed. The closed forms for `z0..n4` printed in section 2 of the output are the actual claim and are unverified against any external anchor.

**Required change:**

After line 81 (immediately following `n4 = dlin(N4)`), and before line 83 (`slips = {...}`), insert the following six closed-form anchors. The RHSs match the printed forms in `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.txt` lines 36-47:

```
    assert_zero("z0 closed form", z0 - (Delta * q1 - Q * d1) / Delta**2)
    assert_zero(
        "z2 closed form",
        z2 - (-Delta**2 * h1 + Delta * (H * d1 + Q * s1 + S2 * q1) - 2 * Q * S2 * d1) / Delta**3,
    )
    assert_zero(
        "z4 closed form",
        z4
        - (
            -Delta**2 * H * s1
            - Delta**2 * S2 * h1
            - Delta**2 * q1
            + 2 * Delta * H * S2 * d1
            + 2 * Delta * Q * S2 * s1
            + 2 * Delta * Q * d1
            + Delta * S2**2 * q1
            - 3 * Q * S2**2 * d1
        )
        / Delta**4,
    )
    assert_zero("n0 closed form", n0 - 2 * P * (Delta * p1 - P * d1) / Delta**3)
    assert_zero(
        "n2 closed form",
        n2
        - (
            -(
                2 * Delta**2 * (Gw * p1 + P * g1)
                - 2 * Delta * P * (2 * Gw * d1 + P * s1 + 2 * S2 * p1)
                + 6 * P**2 * S2 * d1
            )
            / Delta**4
        ),
    )
    assert_zero(
        "n4 closed form",
        n4
        - 2
        * (
            Delta**3 * Gw * g1
            - Delta**2 * Gw**2 * d1
            - 2 * Delta**2 * Gw * P * s1
            - 2 * Delta**2 * Gw * S2 * p1
            - 2 * Delta**2 * P * S2 * g1
            - 2 * Delta**2 * P * p1
            + 6 * Delta * Gw * P * S2 * d1
            + 3 * Delta * P**2 * S2 * s1
            + 3 * Delta * P**2 * d1
            + 3 * Delta * P * S2**2 * p1
            - 6 * P**2 * S2**2 * d1
        )
        / Delta**5,
    )
```

Leave the existing A1-A6 Frechet-vs-series checks in place; they remain useful as method-consistency tests on top of the closed-form anchors.

**Verification:**
After patching, the script runs six new `assert_zero` calls labelled `z0 closed form` ... `n4 closed form`, and the saved output transcript will (on re-run) contain those labels. Re-running the script still exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`
- summary: Added closed-form anchors for z0, z2, z4, n0, n2, and n4 before the Frechet-vs-series consistency checks.
- deviation: none
