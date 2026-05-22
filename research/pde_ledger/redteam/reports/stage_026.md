---
unit_id: 026
batch: II.1
auditor_model: claude-opus-4-7
audit_date: 2026-05-21T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 026 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.txt`

## What the script claims to verify

The docstring says the script verifies (a) the N/N constant zero mode `u0 = 1/sqrt(L)` and the D/N half-wave family `f_n = sqrt(2/L) sin((n+1/2) pi s/L)` are unit-normalised, (b) the general overlap law `kappa_n = sqrt(2)/((n+1/2) pi)` and the lowest-branch `kappa = 2 sqrt(2)/pi`, (c) that the four "axial overlaps" `I_(eta,phi)`, `I_(eta,u)`, `I_(eta,w)`, `I_(u,w)` evaluate on the chosen branch and feed compact closed forms for `C, G_U, G_W, R, Delta, Q, P, B0, Z0, N0, D0`, and (d) an exact branch-level normalisation condition `mhat^2 N0/D0 = target` is solvable for the required wall stiffness `K_req`. (a)-(b) are genuine, but (c)-(d) collapse into tautological reverse-substitution checks (see findings) and the Mathematica `.wl` is a line-by-line transliteration of the SymPy script rather than an independent re-derivation.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 83 | `expect_zero("int u0^2 - 1", ...)` | yes |
| A2 | sympy | 84 | `expect_zero("int f0^2 - 1", ...)` | yes |
| A3 | sympy | 105 | `expect_zero("kappa_n - expected", ...)` | yes |
| A4 | sympy | 109 | `expect_zero("kappa - 2*sqrt(2)/pi", ...)` | yes |
| A5 | sympy | 136 | `expect_zero("I_(eta,phi) - kappa", ...)` | partial (kappa already by construction) |
| A6 | sympy | 137 | `expect_zero("I_(eta,u) - 1", ...)` | yes |
| A7 | sympy | 138 | `expect_zero("I_(eta,w) - kappa", ...)` | no (tautological) |
| A8 | sympy | 139 | `expect_zero("I_(u,w) - kappa", ...)` | no (tautological) |
| A9 | sympy | 181 | `expect_zero("Delta - Delta_expected", ...)` | no (reverse substitution) |
| A10 | sympy | 182 | `expect_zero("Q - Q_expected", ...)` | no (reverse substitution) |
| A11 | sympy | 183 | `expect_zero("P - P_expected", ...)` | no (reverse substitution) |
| A12 | sympy | 184 | `expect_zero("B0 - B0_expected", ...)` | no (reverse substitution) |
| A13 | sympy | 185 | `expect_zero("P0 - P0_expected", ...)` | no (reverse substitution) |
| A14 | sympy | 214 | `expect_zero("K_req - expected", ...)` | no (algebraic rearrangement) |
| A15 | mathematica | 48 | `expectZero["int u0^2 - 1", ...]` | yes |
| A16 | mathematica | 49 | `expectZero["int f0^2 - 1", ...]` | yes |
| A17 | mathematica | 63 | `expectZero["kappa_n - expected", ...]` | yes |
| A18 | mathematica | 67 | `expectZero["kappa - 2*sqrt(2)/pi", ...]` | yes |
| A19 | mathematica | 87 | `expectZero["I_(eta,phi) - kappa", ...]` | partial |
| A20 | mathematica | 88 | `expectZero["I_(eta,u) - 1", ...]` | yes |
| A21 | mathematica | 89 | `expectZero["I_(eta,w) - kappa", ...]` | no (tautological — same integrand as A19) |
| A22 | mathematica | 90 | `expectZero["I_(u,w) - kappa", ...]` | no (tautological — same integrand) |
| A23 | mathematica | 127-131 | `expectZero["Delta/Q/P/B0/P0 - expected", ...]` | no (reverse substitution) |
| A24 | mathematica | 151 | `expectZero["K_req - expected", ...]` | no (algebraic rearrangement) |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:125-139`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:76-90`

**What's wrong:**
The script claims to compute four distinct axial overlaps `I_(eta,phi)`, `I_(eta,u)`, `I_(eta,w)`, `I_(u,w)` corresponding to four different mode pairs (eta-phi, eta-u, eta-w, u-w). In the SymPy script:

```python
I_eta_phi = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
I_eta_u   = sp.simplify(sp.integrate(u0 * u0, (s, 0, L)))
I_eta_w   = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
I_u_w     = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
```

Three of the four integrals are byte-identical: `integrate(u0 * f0, ...)`. The subsequent assertions

```python
expect_zero("I_(eta,phi) - kappa", I_eta_phi - kappa)
expect_zero("I_(eta,w) - kappa",   I_eta_w   - kappa)
expect_zero("I_(u,w) - kappa",     I_u_w     - kappa)
```

cannot fail because all three left-hand sides were computed by the same call as `kappa`. The Mathematica script repeats this exactly (lines 76-90: `iEtaPhi`, `iEtaW`, `iUW` all equal `Integrate[u0*f0, {s, 0, l}]`). The check `I_(eta,phi) - kappa == 0` is the only one of the three that has even residual content (it verifies that `integrate(u0*f0)` equals the symbolic form `2 sqrt(2)/pi`); the other two are pure copies and add no information.

**Why this matters:**
The script's docstring lists "exact substitution of the concrete overlaps into the Stage-8 quantities (C, G_U, G_W, R, Delta, Q, P, B0, Z0, N0, D0)" as a primary deliverable. With three of four overlaps collapsing into the same integral, the downstream coefficients `G_W` and `R` are forced to be proportional to `kappa lambda_W` and `kappa lambda_R` by construction, not because the physics demands it. A reviewer who reads only the assertions would conclude that four independent overlap pairings have been verified when in fact only two distinct integrals (`u0*u0` and `u0*f0`) have been computed.

**Required change:**
Either (a) collapse the three duplicate definitions into a single named variable and document explicitly that on the chosen branch all three overlaps reduce to `kappa` because the script identifies eta with u0 and (phi, w, and the f-side of (u,w)) with f0; or (b) if the four pairs are physically distinct, define the four distinct mode functions (eta, phi, u, w) and integrate the correct product for each. The minimal scoped fix is (a): replace lines 125-128 with a single overlap definition and rewrite assertions to make the identification explicit, so the check `I - kappa == 0` is exercised exactly once (and clearly labelled as the one non-trivial check). Apply the same change to the `.wl` script at lines 76-79 and 87-90.

**Verification:**
After the fix, the SymPy script should contain exactly one `integrate(u0 * f0, ...)` call followed by exactly one `expect_zero(... - 2*sqrt(2)/pi, ...)` (or `... - kappa`) assertion. The reduced-coefficient computations `G_W`, `R`, `C` should reference the single overlap variable. Output file should show the single overlap line, not three identical ones.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:141-185`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:92-131`

**What's wrong:**
Section III.2 in the SymPy script defines

```python
C  = lambda_B * I_eta_phi
GU = lambda_U * I_eta_u
GW = lambda_W * I_eta_w
R  = lambda_R * I_u_w
Delta = Omega_U**2 * Omega_W**2 - R**2
Q     = GU**2 * Omega_W**2 + 2 * GU * GW * R + GW**2 * Omega_U**2
P     = Omega_U**2 * GW + R * GU
B0    = C**2 / varpi**2
```

then defines

```python
Delta_expected = Omega_U**2*Omega_W**2 - kappa**2*lambda_R**2
Q_expected     = lambda_U**2*Omega_W**2 + 2*kappa**2*lambda_U*lambda_W*lambda_R + kappa**2*lambda_W**2*Omega_U**2
P_expected     = kappa*(Omega_U**2*lambda_W + lambda_R*lambda_U)
B0_expected    = kappa**2*lambda_B**2/varpi**2
```

and asserts each Delta/Q/P/B0 equals its `_expected`. But because `I_eta_phi = I_eta_w = I_u_w = kappa` and `I_eta_u = 1` were already established (or, per F1, by construction), the `_expected` forms are *purely* what you get by substituting `kappa` for the overlap symbols in the prior definitions. The simplify call cannot fail unless the substitution itself is wrong. The check verifies that SymPy can substitute, not any physical claim.

The same applies to `P0 - P0_expected`: `P0 = N0/D0 = P^2/(Delta^2 (K - B0 - Q/Delta)) = P^2/(Delta (K Delta - Delta B0 - Q))`, which is exactly `P0_expected` after algebraic rearrangement.

The Mathematica script duplicates this pattern verbatim (lines 121-131).

**Why this matters:**
These assertions appear in the output as `PASS: Delta - Delta_expected`, etc., giving the impression that nontrivial closed-form formulas have been independently rederived. They have not. Anyone relying on this output as evidence that the Stage-8 compact forms are correct is misled — the script only checks that substituting `I_(...) -> kappa` into one expression and into another expression yields the same result.

**Required change:**
Remove the `_expected` constructions and the corresponding `expect_zero` lines (sympy 171-185, wl 121-131). The closed forms for Delta, Q, P, B0, P0 in terms of `kappa, lambda_*, Omega_*, varpi, K` are produced by the prior `print` statements (sympy 161-168, wl 108-119) and are sufficient as displayed output. If the intent was to verify the Stage-8 compact forms, that verification must come from an independent source (e.g., a hand-derived literal substituted at the assertion level — but per scope rules, do not import new physics). The honest fix is deletion of the tautological assertions; what remains is the algebraic-display block plus the non-tautological checks A1-A4, A6.

**Verification:**
After the fix, the SymPy script should not contain the symbols `Delta_expected`, `Q_expected`, `P_expected`, `B0_expected`, `P0_expected`. The output file should no longer show `PASS: Delta - Delta_expected` etc. The `.wl` file should similarly drop `deltaExpected`, `qExpected`, `pExpected`, `b0Expected`, `p0Expected`.

### F3 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:199-214`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:140-151`

**What's wrong:**
The script defines

```python
residual = sp.simplify(mhat**2 * N0 / D0 - target)
K_req = sp.solve(sp.Eq(residual, 0), K)[0]
K_req_expected = sp.simplify(B0 + Q/Delta + mhat**2 * P**2 / (target * Delta**2))
expect_zero("K_req - expected", K_req - K_req_expected)
```

With `N0 = P^2/Delta^2` and `D0 = K - B0 - Z0 = K - B0 - Q/Delta`, the equation `mhat^2 N0/D0 = target` is linear in K and yields exactly

```
K = B0 + Q/Delta + mhat^2 P^2 / (target * Delta^2)
```

i.e., `K_req_expected`. The assertion `K_req - K_req_expected == 0` is the statement "SymPy's `solve` performs algebra correctly on a one-line linear equation". No physics is exercised. The Mathematica script reproduces the same construction (wl:145-151) and the same tautology.

**Why this matters:**
The docstring's "exact solution of the normalization equation for the required wall stiffness K_req on the branch" is presented as a deliverable. The current check only validates the symbolic solver, not the claim that this `K_req` actually achieves the target on the branch. The non-tautological version would substitute `K -> K_req` back into the original residual and check that the residual vanishes — that is a different (and substantively weaker) statement, but at least it is not a guaranteed-pass.

**Required change:**
Replace the `K_req_expected` rebuild with a back-substitution check. Specifically:

- Compute `K_req = sp.solve(sp.Eq(residual, 0), K)[0]` as before.
- Drop the `K_req_expected` construction at lines 207-209.
- Replace `expect_zero("K_req - expected", K_req - K_req_expected)` with `expect_zero("residual @ K_req", residual.subs(K, K_req))`, which asks: does the residual actually vanish when K is replaced by the solver's output? This still does not test any physics beyond `solve`, but it is no longer comparing the solver's output to a hand-rearranged form of the same equation.

Apply the parallel edit to the `.wl` script (replace lines 147-151 with a `residual /. k -> kReq` check).

**Verification:**
After the fix, the SymPy script should contain `residual.subs(K, K_req)` (or equivalent) inside an `expect_zero`, and should not define `K_req_expected`. The output file should show `residual @ K_req = 0` instead of `K_req - expected = 0`.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:38-155` (whole file)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:72-218` (whole file)

**What's wrong:**
The Mathematica script is a one-to-one transliteration of the SymPy script, not an independent re-derivation. Concretely:

- Same section partition (I/II/III/IV) with identical banners.
- Same function decomposition: `concrete_modes` / `concreteModes`, `overlap_law` / `overlapLaw`, `branch_substitution` / `branchSubstitution`, `normalization_test` / `normalizationTest`, each returning the same tuple of objects.
- Same variable choreography. SymPy has `I_eta_phi, I_eta_u, I_eta_w, I_u_w` (script lines 125-128); Mathematica has `iEtaPhi, iEtaU, iEtaW, iUW` (wl lines 76-79), with the same three duplicates.
- Same `_expected` constructions in the same order. SymPy lines 171-179 define `Delta_expected, Q_expected, P_expected, B0_expected, P0_expected`; wl lines 121-125 define `deltaExpected, qExpected, pExpected, b0Expected, p0Expected` with the same algebraic content.
- Same solve-then-rebuild structure for `K_req`. SymPy lines 206-209 do `solve` then construct a hand expected form; wl lines 145-147 do the same.

Sample side-by-side:

SymPy (lines 207-209):
```python
K_req_expected = sp.simplify(
    B0 + Q / Delta + mhat**2 * P**2 / (target * Delta**2)
)
```

Mathematica (line 147):
```
kReqExpected = FullSimplify[b0 + q/delta + mhat^2*p^2/(target*delta^2), Assumptions -> $Assumptions];
```

This is symbol-for-symbol translation. An independent Mathematica derivation would, at minimum, perform the integrals using a different ansatz (e.g., orthogonality via Fourier series rather than direct `Integrate`), simplify under a different normal form, or attack the K_req solve from a different angle (e.g., construct `mhat^2 N0` and `D0 * target` separately and equate). The current file does none of this.

**Why this matters:**
The second-engine policy exists to catch SymPy-specific transcription or simplification errors. A transliteration cannot detect such errors because both scripts run the same arithmetic. The user's notes/MATHEMATICA_MIRROR_POLICY explicitly forbids this pattern (per project memory references; the auditor has not read that file but the framework cited it).

**Required change:**
Re-author the Mathematica script so that it derives the key results independently. Minimum acceptable independence:

- Replace `Integrate[u0*fN, {s, 0, l}]` with an explicit Fourier-sine-series computation: write `1 = Sum[c_k * Sin[(k+1/2) Pi s/l], {k, 0, Infinity}]` weights via inner-product projection, and identify `kappa_n` as the coefficient `c_n` (multiplied by the normalisation `Sqrt[2/l]`). This forces an independent path to `kappa_n = Sqrt[2]/((n+1/2) Pi)`.
- For the reduced coefficients, do not define `_expected` variables (these are dropped under F2 anyway). Independence here is automatic once F2 is applied.
- For `K_req`, instead of `Solve[residual == 0, k]` followed by a hand-rebuilt expected form, construct `lhs = mhat^2 N0` and `rhs = target * D0`, then `kReq = First[k /. Solve[lhs == rhs * (k - b0 - z0)/((k - b0 - z0)), k]]` (or equivalent) — the point is to not literally copy the SymPy `solve` call.

Apply this rewrite confined to the `.wl` file. Keep all assertions that survive F1-F3.

**Verification:**
After the fix, the `.wl` script should compute `kappaN` by a method whose source code does not match `Integrate[u0*fN, ...]` line-for-line with the SymPy `sp.integrate(u0 * f_n, ...)`. The mtime of the `.wl` should be newer than the `.py`. The auditor will re-read both files and check that section structure differs (different intermediate variable names *and* different computational paths).

## Independent-derivation check (Mathematica)

The Mathematica script is **not** an independent derivation. See F4 above for line-by-line correspondence. Sample triple:

| SymPy (lines) | Mathematica (lines) |
|---|---|
| `f_n = sp.sqrt(2 / L) * sp.sin((n + sp.Rational(1, 2)) * sp.pi * s / L)` (76) | `fN = Sqrt[2/l]*Sin[(n + 1/2)*Pi*s/l];` (41) |
| `kappa_n = sp.simplify(sp.integrate(u0 * f_n, (s, 0, L)))` (97) | `kappaN = FullSimplify[Integrate[u0*fN, {s, 0, l}], ...];` (56) |
| `K_req_expected = sp.simplify(B0 + Q/Delta + mhat**2 * P**2 / (target * Delta**2))` (207-209) | `kReqExpected = FullSimplify[b0 + q/delta + mhat^2*p^2/(target*delta^2), ...]` (147) |

Same symbols, same operations, same order. This is `mathematica_transliteration`.

## Engine cross-check

Both engines produce the same final symbolic expressions:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `kappa_n` | `2*sqrt(2)/(pi*(2*n + 1))` | `(2*Sqrt[2])/(Pi + 2*n*Pi)` |
| `kappa` | `2*sqrt(2)/pi` | `(2*Sqrt[2])/Pi` |
| `Delta` | `Omega_U**2*Omega_W**2 - 8*lambda_R**2/pi**2` | `omegaU^2*omegaW^2 - (8*lambdaR^2)/Pi^2` |
| `B0` | `8*lambda_B**2/(pi**2*varpi**2)` | `(8*lambdaB^2)/(Pi^2*varpi^2)` |
| `K_req` (compact form) | identical structure | identical structure |

The engines agree. (Disagreement would have been a positive sign of independence; the agreement here is consistent with the transliteration finding in F4.)

## Verdict justification

Verdict: `findings`. The script verifies the unit-normalisation integrals (A1, A2), the general overlap law `kappa_n = sqrt(2)/((n+1/2) pi)` (A3), and the lowest-branch value `kappa = 2 sqrt(2)/pi` (A4) — these are real checks against the closed forms of the trigonometric integrals. The remainder of the script (Section III.2 reduced coefficients, Section IV K_req) is structured so that every assertion compares an expression to a hand-rearranged version of itself; none of those checks can fail without breaking SymPy/Mathematica's algebra engine. Three of four "axial overlaps" are byte-identical integrals labelled differently (F1). The `.wl` script is a line-by-line port (F4). No engine disagreement, no symbol-assumption errors, no stale outputs, no missing scripts. No `UNFIXABLE` or `CRITICAL_DOWNSTREAM` flag — the corrections are local script rewrites with no downstream impact (the corrected K_req formula and Stage-8 reduced coefficients in the script's printed output remain unchanged; only the assertions and `_expected` rebuilds change).

Attacks attempted and what they showed:
- Tried to make A5/A7/A8 fail by hunting for a typo in the integrand → all three are the literal same integral, so failure impossible (F1).
- Tried to find a non-trivial residual in A9-A13 by varying which symbol enters `_expected` differently → every `_expected` is built by symbol substitution from the same defining expressions; the difference is identically zero (F2).
- Tried to find a hidden non-tautology in A14 by checking whether `sp.solve` could produce a different branch → the equation is linear in K, only one branch exists, and the "expected" is the literal rearrangement of the same equation (F3).
- Checked symbol assumptions for hidden positivity/real conflicts: `s, L > 0`, `n integer >= 0`, all couplings real, Omegas/varpi/G/c/cs/a/mhat positive. The trigonometric integrals are well-posed under these, and `Delta` is allowed to be sign-indefinite (as it should be). No assumption error.
- Checked output freshness: both `.txt` files are newer than their scripts. Not stale.

## Self-test notes

Walked through each proposed directive change as if executing it: (1) F1's single-overlap collapse leaves `kappa = 2 sqrt(2)/pi` as the unique residual check, which we already know holds; substitution path for `G_W`, `R` is preserved (`G_W = lambda_W * kappa`, `R = lambda_R * kappa`). (2) F2's deletion of `_expected` rebuilds removes assertions; what remains is non-tautological (A1-A4, A6). (3) F3's back-substitution `residual.subs(K, K_req)` mentally reduces: `K_req - B0 - Q/Delta = mhat^2 P^2/(target Delta^2)`, so `D0 @ K_req = mhat^2 P^2/(target Delta^2)`; then `mhat^2 N0 / D0 @ K_req = mhat^2 (P^2/Delta^2) / (mhat^2 P^2/(target Delta^2)) = target`, so residual is 0. The back-substitution check passes by construction but is at least not a literal algebraic restatement. (4) F4's independent integration via Fourier-projection would produce `kappa_n = Sqrt[2]/((n+1/2) Pi)` by a different computational route; the result matches but the path differs. No `assert_nonzero`/`assert_zero` proposals are in the directive (this is a corrective audit, not new-script), so the variable-independence and parity traps from the prompt do not apply. Path specifications for the directive use the absolute paths read above; no new files are created (no `missing_verification_script` finding).
