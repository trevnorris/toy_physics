---
unit_id: 027
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 027 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.txt`

## What the script claims to verify

The scripts verify the first nonconstant finite-throat brane family built on exact
N/N (`u0`, `u1`) and D/N (`f0`) axial bases. Concretely they assert: (i) orthonormality
of the basis functions on `s in [0, L]`; (ii) the closed-form overlap constants
`kappa0 = 2*sqrt(2)/pi` and `kappa1 = -4/(3*pi)`, and the linear overlap law
`kappa(theta) = kappa0*cos(theta) + kappa1*sin(theta)` together with its norm
`rho = 2*sqrt(22)/(3*pi)`; (iii) the blind-angle value `kappa(blind) = 0` and the
max-coupling value `kappa(max) = rho` with `sin^2(theta_max) = 2/11`; (iv) the
wall-stiffness expectation `K_geo(theta) = K_eta + 6*T_Omega + T_w*pi^2*sin^2(theta)/L^2`
obtained from the wall operator `G_eta = -T_w*chi_ss + (K_eta + 6*T_Omega)*chi`; and
(v) the profile-dressed Stage-8/9 branch quantities `C, G_U, G_W, R, Delta, Q, P,
B0, Z0, N0, D0` and the blind-angle no-go that forces the outgoing-quadrupole
normalization residual to zero (against a strictly positive target).

## Assertion inventory

| #   | Script      | Line     | Form                                                                                  | Anchored to claim? |
|-----|-------------|----------|---------------------------------------------------------------------------------------|--------------------|
| A1  | sympy       | 87-90    | `int u_i u_j - delta_ij == 0`, `int f0^2 - 1 == 0`                                    | yes                |
| A2  | sympy       | 114-120  | `kappa0 - 2*sqrt(2)/pi`, `kappa1 + 4/(3*pi)`, `kappa - (k0 c + k1 s)`, `rho - ...`    | yes                |
| A3  | sympy       | 141-143  | `kappa(blind) == 0`, `kappa(max) - rho == 0`, `sin^2(theta_max) - 2/11 == 0`          | yes                |
| A4  | sympy       | 161, 170-173 | `K_geo - (K_eta + 6 T_Omega + T_w pi^2 sin^2/L^2) == 0`, `K_geo(theta_max) - ...` | yes                |
| A5  | sympy       | 221-224  | `Delta/Q/P/B0 - expected == 0`                                                        | partial (substitution chain) |
| A6  | sympy       | 230-238  | `kappa(0) - kappa0`, `K_geo(0) - (K_eta+6T_Omega)`, `Delta(0) - ...`                  | yes                |
| A7  | sympy       | 271-273  | `P(blind) == 0`, `N0(blind) == 0`, `lhs(blind) == 0`                                  | yes                |
| B1  | mathematica | 49-52    | `int u_i u_j` and `int f0^2 - 1` checks                                               | yes                |
| B2  | mathematica | 70-73    | overlap constants and law                                                             | yes                |
| B3  | mathematica | 85-87    | blind/max trigonometric values                                                        | yes                |
| B4  | mathematica | 98, 103  | `kGeo - kGeoExpected == 0`, `K_geo(theta_max) - ...`                                  | **no (tautological — see F1)** |
| B5  | mathematica | 142-145  | `delta/q/p/b0 - expected == 0`                                                        | partial (substitution chain) |
| B6  | mathematica | 150-155  | recovery at theta=0                                                                   | partial (relies on hard-coded kGeo) |
| B7  | mathematica | 171-172  | `P(blind) == 0`, `N0(blind) == 0`                                                     | yes                |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl:91-105`

**What's wrong:**

The Mathematica `wallStiffness[]` module does NOT perform the wall-stiffness
integral. It simply assigns the closed-form answer to both sides of the assertion:

```
kGeo = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;
kGeoExpected = kGeo;
...
expectZero["K_geo - expected", kGeo - kGeoExpected];
```

`kGeo - kGeoExpected` is identically zero by construction — no integral is taken,
no derivative of `chi` is computed, no wall operator `G_eta` is built. The
assertion cannot fail no matter what the correct integral evaluates to. This
violates the docstring claim that the script verifies "exact wall stiffness
expectation `K_geo(theta)`" and contradicts what the SymPy script actually does
at lines 156–158:

```
G_eta = -T_w * sp.diff(chi, s, 2) + (K_eta + 6 * T_Omega) * chi
K_geo = sp.simplify(sp.integrate(chi * G_eta, (s, 0, L)))
K_geo_expected = sp.simplify(K_eta + 6 * T_Omega + T_w * sp.pi**2 * sp.sin(theta) ** 2 / L**2)
```

The `K_geo(theta_max)` check at line 103 of the .wl is also tautological for
the same reason: it substitutes into the hard-coded expression, never into a
computed integral.

The `(TrigExpand[kGeo] /. theta0) - (kEta + 6*tOmega)` check at line 151
inherits the same defect.

**Why this matters:**

The Mathematica engine is supposed to provide an independent re-derivation of
the unit's main physical claims. Section III's wall-stiffness expectation
`K_geo(theta) = K_eta + 6*T_Omega + T_w*pi^2*sin^2(theta)/L^2` is one of those
claims (the dependence on `theta` through `sin^2(theta)` is the central new
result of this finite-throat family — it is what makes the throat "nonconstant
axial"). Right now Mathematica certifies nothing here. The unit's two-engine
guarantee is broken for this claim; if the SymPy derivation contained a sign,
prefactor, or eigenvalue error in `chi_ss`, Mathematica would silently agree.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl`,
replace the body of `wallStiffness[]` (currently lines 91–105) so that `kGeo`
is computed from the explicit wall operator and basis. Specifically:

1. Build `gEta = -tW*D[chi, {s, 2}] + (kEta + 6*tOmega)*chi`.
2. Compute `kGeo = FullSimplify[Integrate[chi*gEta, {s, 0, l}], Assumptions -> $Assumptions]`.
3. Set `kGeoExpected = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2` (literal closed form, separate symbol from `kGeo`).
4. Keep the `expectZero["K_geo - expected", kGeo - kGeoExpected]` assertion exactly as written; it will now exercise the derivation.
5. Keep the `K_geo(theta_max)` assertion at line 103 as is; with `kGeo` now derived, this becomes a real check of the substitution.

Do not touch any other module. Do not rename existing variables.

**Verification:**

After Codex applies, the verifier runs `redteam exec-mathematica 027`. The
saved output must (a) still exit 0, (b) still print `PASS: K_geo - expected`,
and (c) the new `kGeo` value printed by `Print["K_geo(theta) = ", fmt[kGeo]]`
must show the explicit symbolic form
`kEta + 6*tOmega + (Pi^2*tW*Sin[theta]^2)/l^2` as the *result of FullSimplify on
the integral*, not as a direct assignment.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a section-for-section, function-for-function
port of the SymPy script:

- `finiteThroatBasis[]` mirrors `finite_throat_basis()` with identical `u0, u1, f0, chi` definitions and identical orthonormality assertions.
- `overlapLaw[]` mirrors `overlap_law()` with the same intermediate names (`kappa0, kappa1, kappa, rho`), identical `blindSubs`/`maxSubs` dictionaries, and identical assertion order.
- `branchSubstitution[]` mirrors `branch_substitution()` with the same `C, G_U, G_W, R, Delta, Q, P, B0, Z0, N0, D0` variable choreography.
- `normalizationAndNoGo[]` mirrors `normalization_and_no_go()` step by step.

Most of this mirroring is unavoidable because the physical premises (orthonormal
basis, the kappa overlap integral, the closed-form `Delta = Omega_U^2*Omega_W^2 - R^2`)
admit only one symbolic decomposition. The Section IV "expected" checks
(`Delta - deltaExpected`, etc.) are substitution-chain sanity checks in both
engines and share that character by construction.

The one place where genuinely independent computation can disambiguate the two
engines — Section III's wall-stiffness integral — is precisely where the
Mathematica script abandons computation in favour of a hard-coded answer. That
is the substantive transliteration concern and is addressed by F1. After F1
lands, the Mathematica script will contain at least one independently derived
result.

No separate `mathematica_transliteration` finding is filed because (a) F1 already
forces the only physically meaningful re-derivation in this unit, and (b)
inventing additional "differently structured" computations for the overlap
integrals would amount to a stylistic refactor outside this auditor's scope.

## Engine cross-check

Both engines run to completion (exit code 0). The closed-form expressions they
print agree:

- `kappa0 = 2*sqrt(2)/pi` (both)
- `kappa1 = -4/(3*pi)` (both)
- `kappa(theta) = (6*sqrt(2)*cos(theta) - 4*sin(theta))/(3*pi)` (both, equivalent forms)
- `rho = 2*sqrt(22)/(3*pi)` (both)
- `sin^2(theta_max) = 2/11` (both)
- `K_geo(theta) = K_eta + 6*T_Omega + pi^2*T_w*sin^2(theta)/L^2` (SymPy derives it; Mathematica asserts it — engines agree on the form but only one verifies it)
- `K_geo(theta_max) = K_eta + 6*T_Omega + 2*pi^2*T_w/(11*L^2)` (both)
- Section IV Delta/Q/P/B0/Z0/N0/D0: identical symbolic forms in both engines after rewriting `(2 sin(theta) - 3 sqrt(2) cos(theta))^2 = (6 sqrt(2) cos(theta) - 4 sin(theta))^2 / 4` (which holds: `(6√2 c - 4 s)^2 = 4 (3√2 c - 2 s)^2 = 4 (2 s - 3√2 c)^2`).

No engine disagreement.

## Verdict justification

The SymPy script holds up to attack. The orthonormality integrals are real
integrals on a finite interval, the overlap-law identity is non-tautological
(LHS is computed by integration, RHS is the proposed closed form), the
blind/max-angle algebraic substitutions correctly use `cos^2 + sin^2 = 1` for
both branches, the K_geo derivation is genuine (the wall operator is built and
integrated against `chi`), and the blind-angle no-go is driven by `kappa(blind) = 0`
which propagates to `P, N0, lhs`.

The Mathematica script is sound on every step *except* Section III, where the
wall-stiffness "verification" is a self-comparison of a hard-coded expression.
This is a `tautological_check` and the second engine's certification of the
unit's main new claim (the explicit `theta`-dependence of `K_geo`) is therefore
vacuous. F1 captures this and the directive lands a real Integrate-based
derivation. No `stop_cold` flag: fixing F1 does not change any downstream
result — both the symbolic form and the `theta=0` and `theta_max` reductions
remain identical post-fix.

Attacks tried that the scripts survived:

- Sign of `chi_ss`: `u1_ss = -(pi/L)^2 * u1`, so `-T_w * <chi, chi_ss> = +T_w * sin^2(theta) * pi^2/L^2`. Matches the asserted sign.
- Orthogonality of `u0` and `u1`: `int_0^L cos(pi s/L) ds = 0`. Holds.
- Blind angle consistency: `cos^2 + sin^2 = 2/11 + 9/11 = 1`. Holds.
- Max-coupling consistency: `9/11 + 2/11 = 1`. Holds.
- `kappa(blind)`: `(2*sqrt(2)/pi)*(sqrt(2)/sqrt(11)) + (-4/(3*pi))*(3/sqrt(11)) = 4/(pi*sqrt(11)) - 4/(pi*sqrt(11)) = 0`. Holds.
- `kappa(max) - rho`: `(2*sqrt(2)/pi)*(3/sqrt(11)) + (-4/(3*pi))*(-sqrt(2)/sqrt(11)) = 6*sqrt(2)/(pi*sqrt(11)) + 4*sqrt(2)/(3*pi*sqrt(11)) = (18 + 4)*sqrt(2)/(3*pi*sqrt(11)) = 22*sqrt(2)/(3*pi*sqrt(11)) = 2*sqrt(22)/(3*pi) = rho`. Holds.
- `P(blind) = kappa(blind) * (Omega_U^2*lambda_W + lambda_R*lambda_U) = 0 * (...) = 0`. Holds as an identity in symbols.
- `lhs(blind)` denominator nonzero generically: `K_geo(blind) - B0(blind) - Q(blind)/Delta(blind)` depends on free parameters `K_eta, T_Omega, T_w, lambda_*, varpi, Omega_*`, so it is generically nonzero; `simplify` correctly returns 0 for `mhat^2 * 0 / D` as an identity. Holds.

## Self-test notes

I checked the Required change in F1 against three traps:

1. **Variable independence:** `chi` depends on `s` (through `Cos[Pi*s/l]`) and on
   `theta`. `D[chi, {s, 2}]` is not identically zero — it equals
   `Sin[theta]*Sqrt[2/l]*(-(Pi/l)^2)*Cos[Pi*s/l]`. The integrand `chi*gEta` is
   non-trivial in `s` and the integral over `[0, l]` returns a nonzero
   `Sin[theta]^2*Pi^2*tW/l^2` term plus the constant `kEta + 6*tOmega`.
2. **Parity/symmetry of the integral:** On `[0, l]` (not a symmetric interval),
   `int Cos[Pi s/l] ds = 0`, so cross terms vanish; `int Cos[Pi s/l]^2 ds = l/2`,
   so the `Sin[theta]^2` term picks up `(2/l)*(l/2) = 1` and a factor of `Pi^2/l^2`
   from the second derivative — exactly the asserted `tW*Pi^2*Sin[theta]^2/l^2`.
3. **Trivial-case pre-check:** At `theta=0`: `kGeo = kEta + 6*tOmega` (matches IV.2
   assertion). At `theta=theta_max` with `Sin[theta]^2 = 2/11`:
   `kGeo = kEta + 6*tOmega + 2*tW*Pi^2/(11*l^2)` (matches line 103). The
   Integrate-based derivation reproduces both reductions consistently.

Path specification: F1 targets only the `.wl` file at the absolute path under
`mathematica/`; no missing-script subcase applies.
