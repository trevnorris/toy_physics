---
unit_id: 054
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 054 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage054_robin_softening_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage054_robin_softening_mathematica_audit.txt`

## What the script claims to verify

The two scripts purport to verify a Robin-compliance "softening" analysis for a 1D wave-equation eigenproblem on `s in [0, L]`. Starting from `psi = A cos(ks) + B sin(ks)` with a Neumann condition at `s = L`, they claim to derive (a) the Robin characteristic equation `k tan(kL) = h` and its dimensionless form `y tan y = eta`; (b) the closed-form softening factor `A_K = 1 / [1 - x/4 + x y^2 / pi^2]` from an effective stiffness ratio `K_W^eff / K_phi,0^eff`; (c) the boundary cases `A_K(y=pi/2)=1` (D/N limit) and `A_K(y->0+) = 4/(4-x)` (soft-mouth maximum); and (d) the saturation floor `x_floor = 4 - 4/zeta_req` that follows from solving `A_K,max = zeta_req`. The assertions test residuals against these closed forms.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero("Robin equation -> k tan(kL) - h", char_eq / A + h - k tan(kL))` | yes (B derived via `sp.solve`) |
| A2 | sympy | 42 | `expect_zero("dimensionless form", ...)` | yes |
| A3 | sympy | 59 | `expect_zero("K_W identity", KW - (KX_from_x + ...))` | yes |
| A4 | sympy | 62 | `expect_zero("A_K x-form", AK_x - 1/(1 - x/4 + x y^2/pi^2))` | yes |
| A5 | sympy | 69 | `expect_zero("DN limit", AK_DN - 1)` | yes |
| A6 | sympy | 70 | `expect_zero("soft-mouth limit", AK_soft - 4/(4-x))` | yes |
| A7 | sympy | 85 | `expect_zero("x floor = 4 - 4/zeta_req", x_floor - (4 - 4/zeta_req))` | yes (x_floor derived via `sp.solve`) |
| B1 | mathematica | 40 | `expectZero["Robin equation -> k tan(kL) - h", charEq/a + h - k Tan[k ell]]` | **no — bExpr hardcoded at line 34** |
| B2 | mathematica | 44 | `expectZero["dimensionless form", ...]` | yes |
| B3 | mathematica | 61 | `expectZero["K_W identity", ...]` | yes |
| B4 | mathematica | 71 | `expectZero["A_K x-form", aKX - aKSym]` | yes |
| B5 | mathematica | 72 | `expectZero["DN limit", aKDN - 1]` | yes |
| B6 | mathematica | 73 | `expectZero["soft-mouth limit", aKSoft - 4/(4 - x)]` | yes |
| B7 | mathematica | 83 | `expectZero["x floor = 4 - 4/zeta_req", xFloor - (4 - 4/zetaReq)]` | **no — xFloor hardcoded at line 78** |
| B8 | mathematica | 84 | `expectZero["A_K,max(x_floor) - zeta_req", (aKMax /. x -> xFloor) - zetaReq]` | partial (depends on B7) |

## Findings

### F1 — hardcoded_result

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:34`

**What's wrong:**
Line 34 writes the Neumann-derived coefficient as a hardcoded literal:

```
bExpr = FullSimplify[a Tan[k ell], Assumptions -> $Assumptions];
```

There is no call to `Solve` or any other derivation step that obtains `b = a Tan[k ell]` from the bottom Neumann condition `D[psi, s] /. s -> ell == 0`. The very next assertion at line 40

```
expectZero["Robin equation -> k tan(kL) - h", charEq/a + h - k Tan[k ell]];
```

is therefore tautological by construction: with `bExpr = a Tan[k ell]` substituted in, `psiBN = a Cos[k s] + a Tan[k ell] Sin[k s]`, so `D[psiBN, s] /. s -> 0 = a k Tan[k ell]` and `charEq/a = k Tan[k ell] - h`, which trivially cancels with the offset `+h - k Tan[k ell]`. The check cannot fail no matter what physics it purports to encode. The corresponding SymPy script at line 33 does `sp.solve(sp.Eq(sp.diff(psi, s).subs(s, L), 0), B)[0]`, an actual derivation from the Neumann condition; the Mathematica script skips this step and asserts the answer.

**Why this matters:**
The Robin characteristic equation `k tan(kL) = h` is the central physical result of the first half of the unit. With `b` written in by hand, the Mathematica engine no longer independently verifies the Neumann derivation; it only restates it. If the Neumann derivation were ever wrong (sign error, missing factor of k, etc.), the Mathematica script would still pass.

**Required change:**
Replace line 34 so that `bExpr` is derived by solving the Neumann condition explicitly. Use Mathematica's `Solve` on `D[psi, s] /. s -> ell == 0` for `b`.

**Verification:**
After the fix, the `Robin equation -> k tan(kL) - h` assertion at line 40 still passes (output line 16 of the .txt remains `PASS: Robin equation -> k tan(kL) - h`), and the printed `B from Neumann bottom = a*Tan[ell*k]` (output line 13) is now the result of a Solve call rather than a hand-stated value. The verifier confirms by reading the updated source line 34: it must contain `Solve[...]` operating on the Neumann condition.

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:78,83`

**What's wrong:**
Line 78 defines:

```
xFloor = FullSimplify[4 - 4/zetaReq, Assumptions -> zetaReq > 0];
```

This is the answer the script claims to derive from `A_K,max = zeta_req`. The SymPy counterpart at line 83 actually derives it via `sp.solve(sp.Eq(AK_max, zeta_req), x)[0]`. The Mathematica version writes the closed form directly, then at line 83 asserts

```
expectZero["x floor = 4 - 4/zeta_req", xFloor - (4 - 4/zetaReq)];
```

i.e. `(4 - 4/zetaReq) - (4 - 4/zetaReq) == 0`, which is a tautology and cannot fail. The downstream assertion at line 84 also rides on this hardcoded `xFloor`.

**Why this matters:**
The saturation threshold `x_floor = 4 - 4/zeta_req` is the final ledger result of the unit. If the inversion of `A_K,max = zeta_req` were ever miscarried (e.g., the wrong root, a sign error), the Mathematica check could still pass because nothing is being solved.

**Required change:**
Replace line 78 so that `xFloor` is obtained by `Solve` on `aKMax == zetaReq` for `x`, instead of writing the literal closed form. Line 83's assertion then becomes a real test that the solve returned the claimed form.

**Verification:**
After the fix, output line 37 should still read `x floor at saturation = 4 - 4/zetaReq` and the PASS at line 39 remains. The verifier confirms the source line 78 now contains a `Solve[aKMax == zetaReq, x]` call (or equivalent) rather than a hand-stated `4 - 4/zetaReq`.

## Independent-derivation check (Mathematica)

The `.wl` script mirrors the `.py` script in structure: identical variable choreography (`psi`, `bExpr`, `psiBN`, `charEq`; then `kPhi`, `kWEff`, `aK`; then `kXFromX`, `tXFromX`, `aKX`, `aKSym`, `aKDN`, `aKSoft`; then `ineqRhs`, `yReqSq`, `aKMax`, `xFloor`). Most assertions are independently meaningful (B3–B6 substitute one closed form and check it equals another, which is genuine algebra), but at two key derivation steps — solving the Neumann condition for `b` (line 34) and inverting `A_K,max = zeta_req` for `x` (line 78) — the Mathematica script hardcodes the answer instead of invoking `Solve`, while the SymPy version uses `sp.solve` in both places. The corresponding assertions (B1, B7) are then tautological. This is not a wholesale transliteration finding because the rest of the algebra (the change-of-variables identity, the closed-form match, the boundary limits) is substantive; the two hardcoded answers are the localized failure mode and are captured by F1 and F2.

## Engine cross-check

Both engines pass all their assertions and arrive at the same algebraic content:

- SymPy: `A_K in x,y form = 4*pi**2/(4*x*y**2 - pi**2*(x - 4))`, `A_K,max = -4/(x - 4)`, `x floor at saturation = 4 - 4/zeta_req`.
- Mathematica: `A_K in x,y form = kWBar/(kWBar + kWBar*x*(-1/4 + y^2/Pi^2))`, `A_K,max = -4/(-4 + x)`, `x floor at saturation = 4 - 4/zetaReq`.

These agree on the same closed forms up to algebraic rearrangement. No `engine_disagreement` finding. Note also output freshness: SymPy script mtime Apr 1, output mtime May 11 12:43 (fresh); Mathematica script mtime May 11 11:56, output mtime May 11 12:52 (fresh). No `stale_output` finding.

## Verdict justification

The math content of unit 054 is real and both engines agree on the closed forms. However, the Mathematica engine fails to independently derive two key intermediate results — the Neumann-derived coefficient `b = a Tan[k ell]` and the saturation root `x_floor = 4 - 4/zetaReq` — by hardcoding both answers at lines 34 and 78. The two assertions that depend on these (B1, B7) are therefore tautological, even though they print PASS. The SymPy script does derive both via `sp.solve` and is clean. Two findings; verdict `findings`; no stop-cold (the fixes are local Solve substitutions that will yield the same algebraic answers and not propagate to downstream units' inputs).

## Self-test notes

- **Variable independence**: neither finding introduces a new `D[...]` derivative, so no zero-by-independence trap. The existing `D[psi, s]` at line 36 already takes the derivative w.r.t. `s` of an expression that does depend on `s` — fine.
- **Parity / symmetric-domain integrals**: no integrals over symmetric unbounded domains appear; nothing to check.
- **Trivial-case pre-check for F1**: with `a=1, k=1, ell=Pi/4`, the proposed `Solve[D[a Cos[k s] + b Sin[k s], s] /. s -> ell == 0, b]` yields `b = Tan[Pi/4] = 1`, matching the hardcoded value but now derived; the downstream `charEq/a + h - k Tan[k ell]` reduces to 0 as before. For F2, with `zetaReq=2`, `Solve[4/(4-x) == 2, x]` yields `x = 2`, equal to `4 - 4/2 = 2`. Both substitutions preserve PASS while making the assertions non-tautological.
- **Path specifications**: F1 and F2 target an existing `.wl` file at `mathematica/`; F1/F2 are not `missing_verification_script`, so no new file paths are introduced.
