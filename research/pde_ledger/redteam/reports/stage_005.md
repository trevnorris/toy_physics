---
unit_id: 005
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 005 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.txt`
- mathematica output: `(missing)`

(SymPy script mtime 2026-05-04 12:00, output mtime 2026-05-11 12:38 — output is newer than the script; no staleness. Unit manifest entry has `is_checkpoint: false` and `is_status_only_candidate: false`, so both engines are required.)

## What the script claims to verify

The script's docstring states it derives "the exact projected inhomogeneous Maxwell law from the localized 4+1 bulk equation." Concretely the assertions verify, for a Gaussian transverse kernel `W(w) = exp(-w^2)/sqrt(pi)`: (i) projection commutes with the four brane-direction derivatives `partial_{t,x,y,z}`; (ii) for a polynomial profile `Q = w^3 + w` the boundary term `[W Q]_{-infty}^{infty}` vanishes and the projection-IBP identity `Proj_W[partial_w Q] = -Proj_{W'}[Q]` holds (with the sign-flipped mutant failing as a guard); (iii) for hand-chosen polynomial test fields `G^{0,x,y,z,w}_nu`, `Gamma_nu`, applying the bulk equation operator and projecting term-by-term agrees with the projected form `partial_a Proj_W[G^a_nu] - Proj_{W'}[G^w_nu] + Proj_W[Gamma_nu]`, with explicit Gaussian boundary discharge; (iv) the analytic transverse-leakage value `-Proj_{W'}[w] = 1`; and (v) the projected charge-continuity equation `partial_t Proj_W[J^0] + partial_a Proj_W[J^a] = Proj_{W'}[J^w]` for an analogous polynomial source choice. Sections 1-4 print the *generic* projected law in abstract `sp.Function(...)` notation but contain **no assertions**; only the concrete Section 5 Gaussian-kernel block makes assertions.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 226 | `simplify(Pg(diff(Phi, t)) - diff(Pg(Phi), t)) == 0` (commutation of P_g with brane derivative) | yes |
| A2 | sympy | 230 | `simplify(boundary_value(Wg * (w^3+w))) == 0` (Gaussian boundary term vanishes) | yes |
| A3 | sympy | 231 | `Pg(diff(Q, w)) - (boundary_Q - Pgp(Q)) == 0` (IBP with explicit boundary; reduces to A4 after A2) | partial (redundant with A4 given A2) |
| A4 | sympy | 232 | `Pg(diff(Q, w)) + Pgp(Q) == 0` (boundary-decay IBP identity) | yes |
| A5 | sympy | 233 | `assert_nonzero` mutated IBP sign | yes (mutant-style guard) |
| A6 | sympy | 260 | `boundary_value(Wg * Gwc) == 0` (Gaussian boundary discharge for G^w_nu = w) | yes |
| A7 | sympy | 261 | `project_bulk_lhs - projected_lhs_with_boundary == 0` (projected inhomogeneous law with boundary term, concrete fields) | yes |
| A8 | sympy | 262 | `project_bulk_lhs - projected_lhs == 0` (projected inhomogeneous law in boundary-decay form, concrete fields) | yes |
| A9 | sympy | 263-275 | `assert_nonzero` mutated leakage sign | yes (mutant-style guard) |
| A10 | sympy | 278 | `simplify((-Pgp(w)) - 1) == 0` (analytic transverse-leakage value) | yes |
| A11 | sympy | 279 | `assert_nonzero(-Pgp(w))` (sanity: leakage is not vacuously zero) | yes |
| A12 | sympy | 301 | `boundary_value(Wg * Jwc) == 0` (Gaussian boundary discharge for J^w = w) | yes |
| A13 | sympy | 302 | `project_bulk_cont - (projected_cont + boundary_Jw) == 0` (continuity-with-boundary, concrete sources) | yes |
| A14 | sympy | 303 | `project_bulk_cont - projected_cont == 0` (continuity boundary-decay form, concrete sources) | yes |

There is **no Mathematica script**, so rows B1..Bn do not exist.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- expected: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl` — **does not exist**
- expected: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.txt` — **does not exist**

**What's wrong:**
Unit 005 has `is_checkpoint: false` and `is_status_only_candidate: false` in `redteam/MANIFEST.yaml`, so both a SymPy and a Mathematica audit script are required by the two-engine policy. No `.wl` file matching the naming convention `moving_throat_pde_stage005_*_mathematica_audit.wl` exists anywhere under `/var/projects/toy_physics/research/pde_ledger/mathematica/` (confirmed by directory listing; the convention is observed for stages 001, 002, 003 and many later stages). Consequently the SymPy assertions A1-A14 in `moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py` are not cross-checked by a second engine. Subtype: `missing_mathematica`.

This is *not* a `status_only` carry-forward situation: the SymPy script of stage 005 derives new identities (the projected inhomogeneous Maxwell law with explicit `Proj_{W'}[Z F^{w nu}]` leakage term, and the open-system charge continuity law `partial_t Proj_W[J^0] + partial_a Proj_W[J^a] = Proj_{W'}[J^w]`) and verifies them numerically on a Gaussian kernel. These are first-time claims for this unit, not consolidations of upstream results, so they require independent re-derivation by Mathematica.

**Why this matters:**
The unit's SymPy script asserts non-trivial integration-by-parts identities, an explicit `-Proj_{W'}[w] = 1` Gaussian-moment value, and the sign on the leakage term `-Proj_{W'}[Z F^{w nu}]`. A sign or factor mistake propagated through the SymPy choice of Gaussian profile and polynomial test fields would not be caught: the entire assertion battery is internally consistent with itself by construction (the test fields are simple polynomials, the kernel is a single Gaussian, and `sp.simplify` is applied to both sides). A genuinely independent Mathematica derivation — using e.g. `Integrate[]` directly with `Assumptions -> {w \[Element] Reals}` or different test profiles — would expose any such error. Without it, the PASS in the SymPy output is internal consistency, not verification.

**Required change:**
Create a new Mathematica audit script at:

`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl`

that independently verifies the same physical claims (M1..M5 listed in the directive) by deriving them from first principles in Mathematica, **not** by transliterating the SymPy choreography. Use `Integrate[..., {w, -Infinity, Infinity}, Assumptions -> ...]` for the projections, choose at least one test profile that differs from the SymPy polynomial choice (e.g. mix a Gaussian-decaying profile such as `w Exp[-w^2/2]` for `Q`, and an odd-times-even product such as `Sin[w] Exp[-w^2]` for `G^w_nu`), and assert each identity with `If[Simplify[lhs - rhs] =!= 0, Print["FAIL: ..."]; Exit[1]]`. The Codex directive (file `redteam/directives/stage_005.md`) gives the claim manifest M1..M5.

**Verification:**
The verifier will run `redteam exec-mathematica 005`, confirm the new `.wl` file exists at the path above, that its output file appears at `mathematica/output/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.txt`, and that the script exits 0 after asserting all five claim items. The auditor should then re-run unit 005 with both engines present; engine cross-check (residuals, leakage value `1`, sign on `-Proj_{W'}[J^w]`) should agree.

## Independent-derivation check (Mathematica)

Not applicable: no Mathematica script exists. This finding is filed as `missing_verification_script` rather than `mathematica_transliteration`; once a Mathematica script is added, the next audit pass should verify that it is independent (not a `.py` -> `.wl` port).

## Engine cross-check

Not applicable: only one engine present. After F1 is resolved, the auditor should check that the new Mathematica script's transverse-leakage value `-Proj_{W'}[w]` for a Gaussian kernel equals `+1` (matching SymPy line 278), that the projected charge-continuity residual vanishes for an analogous polynomial source, and that the sign on the leakage term in the inhomogeneous-law check is `-Proj_{W'}[Z F^{w nu}]` (not `+Proj_{W'}[...]`, which the SymPy mutant on lines 264-275 explicitly rules out).

## Verdict justification

`findings` (one finding, `missing_mathematica`). The SymPy script itself holds up under attack on its own terms: I checked that the integration-by-parts identity is exercised non-tautologically (the `assert_nonzero` mutants on lines 233 and 263-275 would catch a sign flip on the `Proj_{W'}` term and on the leakage sign respectively, and the analytic Gaussian moment `-Proj_{W'}[w] = 1` on line 278 is anchored to a real `sp.integrate` call rather than to a hardcoded substitution); the symbol domains (`t, x, y, z, w` real; `mu0, xi` nonzero) are appropriate for the bulk-PDE setting and `xi` is unused in any assertion so no hidden positivity issue exists; the Gaussian kernel decays at `|w| -> infinity` and lines 230 / 260 / 301 verify the boundary-discharge premise rather than assuming it; output mtime (2026-05-11) is newer than script mtime (2026-05-04), so no `stale_output`. Sections 1-4 of the script print generic-form identities without asserting them, which would be `insufficient_verification` if Section 5 did not exercise the same identities concretely — but Section 5 does, so it is not a separate finding. The single defect is structural: the two-engine policy is violated by the absence of a Mathematica script, and that's where this unit fails verification.
