---
unit_id: 026
batch: II.1
created_at: 2026-05-21T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T16:59:28-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 026

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:125-139`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:76-90`

**Issue:** The four "axial overlaps" `I_(eta,phi)`, `I_(eta,u)`, `I_(eta,w)`, `I_(u,w)` are computed using only two distinct integrals: `integrate(u0*u0)` and `integrate(u0*f0)`. Three of the four are byte-identical instances of `integrate(u0*f0)`, so the assertions that each equals `kappa` cannot fail.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`, replace lines 124-139 (from `# Concrete overlaps on the chosen branch:` through `expect_zero("I_(u,w) - kappa", I_u_w - kappa)`) with:

```python
    # Concrete overlaps on the chosen branch. On this branch, eta is identified
    # with the constant zero mode u0, and (phi, w, and the f-leg of the (u,w)
    # pair) are identified with the lowest D/N mode f0. Three of the four
    # nominal axial overlaps therefore reduce to the single integral
    # int_0^L u0(s) f0(s) ds; the fourth is the eta self-overlap.
    overlap_u0_f0 = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
    overlap_u0_u0 = sp.simplify(sp.integrate(u0 * u0, (s, 0, L)))

    I_eta_phi = overlap_u0_f0
    I_eta_u   = overlap_u0_u0
    I_eta_w   = overlap_u0_f0
    I_u_w     = overlap_u0_f0

    subbanner("III.1 — Explicit overlap integrals on the branch")
    print("overlap_u0_f0 =", overlap_u0_f0)
    print("overlap_u0_u0 =", overlap_u0_u0)
    print("I_(eta,phi)   =", I_eta_phi)
    print("I_(eta,u)     =", I_eta_u)
    print("I_(eta,w)     =", I_eta_w)
    print("I_(u,w)       =", I_u_w)

    expect_zero("overlap_u0_f0 - kappa", overlap_u0_f0 - kappa)
    expect_zero("overlap_u0_u0 - 1", overlap_u0_u0 - 1)
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`, replace lines 76-90 (from `iEtaPhi = ...` through `expectZero["I_(u,w) - kappa", iUW - kappa];`) with:

```
  overlapU0F0 = FullSimplify[Integrate[u0*f0, {s, 0, l}], Assumptions -> $Assumptions];
  overlapU0U0 = FullSimplify[Integrate[u0*u0, {s, 0, l}], Assumptions -> $Assumptions];

  iEtaPhi = overlapU0F0;
  iEtaU = overlapU0U0;
  iEtaW = overlapU0F0;
  iUW = overlapU0F0;

  subbanner["III.1 — Explicit overlap integrals on the branch"];
  Print["overlap_u0_f0 = ", fmt[overlapU0F0]];
  Print["overlap_u0_u0 = ", fmt[overlapU0U0]];
  Print["I_(eta,phi)   = ", fmt[iEtaPhi]];
  Print["I_(eta,u)     = ", fmt[iEtaU]];
  Print["I_(eta,w)     = ", fmt[iEtaW]];
  Print["I_(u,w)       = ", fmt[iUW]];

  expectZero["overlap_u0_f0 - kappa", overlapU0F0 - kappa];
  expectZero["overlap_u0_u0 - 1", overlapU0U0 - 1];
```

The downstream definitions of `C, GU, GW, R` (sympy 141-144; wl 92-95) continue to reference `I_eta_phi`, `I_eta_u`, `I_eta_w`, `I_u_w` unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 026` and `redteam exec-mathematica 026`, confirm both scripts exit 0, and confirm that the output files now show two overlap assertions instead of four (`overlap_u0_f0 - kappa` and `overlap_u0_u0 - 1`), with the I_(...) lines printed but no longer asserted.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- summary: Replaced per-label overlap assertions with two explicit branch overlap checks while retaining the I_(...) printed values.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:170-185`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:121-131`

**Issue:** The `_expected` forms for Delta, Q, P, B0, P0 are constructed by substituting `kappa, lambda_*, Omega_*, varpi` into the same defining expressions used moments earlier. The assertion `Delta - Delta_expected == 0` (etc.) only verifies that SymPy/Mathematica can perform symbolic substitution. No physics is exercised.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`, delete lines 170-185 (from `# Check against the hand-derived compact formulas.` through `expect_zero("P0 - P0_expected", P0 - P0_expected)`). The block ends just before `return kappa, Delta, Q, P, B0, Z0, N0, D0` at line 187, which stays. After the deletion, line 169 ends with `print("P0    =", P0)` and the next non-blank line is `return kappa, Delta, Q, P, B0, Z0, N0, D0`.

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`, delete lines 121-131 (from `deltaExpected = ...` through `expectZero["P0 - P0_expected", p0 - p0Expected];`). The function continues at the existing `{kappa, delta, q, p, b0, z0, n0, d0}` return on line 133.

Additionally, in the `Module[{...}]` symbol list at wl line 72, remove the unused symbols `deltaExpected, qExpected, pExpected, b0Expected, p0Expected`. The `Module` declaration becomes:

```
branchSubstitution[] := Module[{u0, f0, kappa, iEtaPhi, iEtaU, iEtaW, iUW, cCoupling, gUeff, gWeff, rEff, delta, q, p, b0, z0, n0, d0, p0},
```

**Verification command:**
After Codex applies, the verifier will confirm:
- SymPy output no longer contains `Delta - Delta_expected = 0`, `Q - Q_expected = 0`, `P - P_expected = 0`, `B0 - B0_expected = 0`, or `P0 - P0_expected = 0`.
- Mathematica output no longer contains `PASS: Delta - Delta_expected` etc.
- Both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- summary: Deleted the substituted expected-form coefficient checks and removed their unused Mathematica module symbols.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:205-214`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:145-151`

**Issue:** `K_req` is obtained by `sp.solve(residual == 0, K)`, then `K_req_expected` is built as the literal algebraic rearrangement of the same equation. The check `K_req - K_req_expected == 0` only confirms the solver. A weaker but non-tautological alternative is a back-substitution check.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`, replace lines 205-214 (from `# Solve exactly for the required wall stiffness.` through `expect_zero("K_req - expected", K_req - K_req_expected)`) with:

```python
    # Solve exactly for the required wall stiffness and verify by back-substitution
    # that the residual vanishes when K is set to the solver's output.
    K_req = sp.solve(sp.Eq(residual, 0), K)[0]

    subbanner("IV.2 — Exact required wall stiffness")
    print("K_req =")
    sp.pprint(sp.simplify(K_req))
    expect_zero("residual @ K_req", residual.subs(K, K_req))
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`, replace lines 145-151 (from `kReq = First[Solve[residual == 0, k]];` through `expectZero["K_req - expected", kReq - kReqExpected];`) with:

```
  kReq = k /. First[Solve[residual == 0, k]];
  kReq = FullSimplify[kReq, Assumptions -> $Assumptions];

  subbanner["IV.2 — Exact required wall stiffness"];
  Print["K_req = ", fmt[kReq]];
  expectZero["residual @ K_req", residual /. k -> kReq];
```

Also remove `kReqExpected` from the `Module[{...}]` symbol list at wl line 136. The declaration becomes:

```
normalizationTest[] := Module[{kappa, delta, q, p, b0, z0, n0, d0, target, residual, kReq, kGeom},
```

**Verification command:**
After Codex applies, the verifier will confirm:
- SymPy output shows `residual @ K_req = 0` (or equivalent) instead of `K_req - expected = 0`.
- Mathematica output shows `PASS: residual @ K_req` instead of `PASS: K_req - expected`.
- Both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- summary: Replaced the expected stiffness rearrangement checks with residual back-substitution checks.
- deviation: none

## F4 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:53-70`

**Issue:** The `.wl` script computes `kappaN` by `Integrate[u0*fN, {s, 0, l}]`, identical in form and choreography to the SymPy `sp.integrate(u0 * f_n, (s, 0, L))`. This is the load-bearing analytical step of the unit; performing it the same way in both engines provides no independent confirmation. After F1-F3 land, the remaining transliterated step is this integral; once it is computed by an independent route, the engine-cross-check has real content.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`, replace the body of `overlapLaw[]` between lines 53 and 70. Specifically, replace:

```
overlapLaw[] := Module[{u0, fN, f0, kappaN, kappa, kappaExpected},
  banner["SECTION II — EXACT OVERLAP LAW"];
  {u0, fN, f0} = concreteModes[];
  kappaN = FullSimplify[Integrate[u0*fN, {s, 0, l}], Assumptions -> $Assumptions];
  kappa = FullSimplify[kappaN /. n -> 0, Assumptions -> $Assumptions];
  kappaExpected = FullSimplify[2*Sqrt[2]/Pi, Assumptions -> $Assumptions];

  subbanner["II.1 — General D/N overlap with the constant zero mode"];
  Print["kappa_n = ", fmt[kappaN]];
  Print["kappa_n_expected = ", fmt[FullSimplify[Sqrt[2]/((n + 1/2)*Pi), Assumptions -> $Assumptions]]];
  expectZero["kappa_n - expected", kappaN - Sqrt[2]/((n + 1/2)*Pi)];

  subbanner["II.2 — Lowest-branch overlap"];
  Print["kappa = ", fmt[kappa]];
  expectZero["kappa - 2*sqrt(2)/pi", kappa - kappaExpected];
  Print["kappa numeric = ", fmt[N[kappa, 15]]];
  {u0, f0, kappa}
];
```

with:

```
overlapLaw[] := Module[{u0, fN, f0, indef, kappaN, kappaNViaFundamentalThm, kappa, kappaExpected},
  banner["SECTION II — EXACT OVERLAP LAW"];
  {u0, fN, f0} = concreteModes[];

  (* Independent path: compute the indefinite integral of u0 * fN with respect
     to s, then evaluate at the boundary s = l and s = 0 by the fundamental
     theorem of calculus. This avoids the single Integrate[..., {s,0,l}] call
     used by the SymPy script and exercises a different code path in Mathematica
     (Integrate without limits + boundary evaluation, vs Integrate with limits). *)
  indef = FullSimplify[Integrate[u0*fN, s], Assumptions -> $Assumptions];
  kappaNViaFundamentalThm = FullSimplify[
    (indef /. s -> l) - (indef /. s -> 0),
    Assumptions -> $Assumptions
  ];

  (* Independent path 2: derive kappa_n algebraically from the trigonometric
     antiderivative of Sin[(n+1/2) Pi s/l], which is -Cos[(n+1/2) Pi s/l] *
     l/((n+1/2) Pi). At s = l the cosine is Cos[(n+1/2) Pi] = 0 for integer n;
     at s = 0 it is 1. So the boundary contribution is l/((n+1/2) Pi), and
     multiplying by the mode normalisation Sqrt[2/l] * (1/Sqrt[l]) = Sqrt[2]/l
     yields Sqrt[2]/((n+1/2) Pi). This is the analytic short form. *)
  kappaN = FullSimplify[
    Sqrt[2]/((n + 1/2)*Pi),
    Assumptions -> $Assumptions
  ];

  subbanner["II.1 — General D/N overlap with the constant zero mode"];
  Print["kappa_n (analytic short form) = ", fmt[kappaN]];
  Print["kappa_n (via fundamental thm) = ", fmt[kappaNViaFundamentalThm]];
  expectZero["kappa_n (analytic) - (fundamental thm)", kappaN - kappaNViaFundamentalThm];

  subbanner["II.2 — Lowest-branch overlap"];
  kappa = FullSimplify[kappaN /. n -> 0, Assumptions -> $Assumptions];
  kappaExpected = FullSimplify[2*Sqrt[2]/Pi, Assumptions -> $Assumptions];
  Print["kappa = ", fmt[kappa]];
  expectZero["kappa - 2*sqrt(2)/pi", kappa - kappaExpected];
  Print["kappa numeric = ", fmt[N[kappa, 15]]];
  {u0, f0, kappa}
];
```

This replaces the single `Integrate[..., {s,0,l}]` call with two independent computations: (i) an indefinite-integral + boundary-evaluation path and (ii) the analytic short form derived from the cosine boundary values. The assertion now compares them, which only passes if both paths agree.

**Verification command:**
After Codex applies, the verifier will confirm:
- The `.wl` file no longer contains the substring `Integrate[u0*fN, {s, 0, l}]` (with explicit limits). It should contain `Integrate[u0*fN, s]` (indefinite) and a boundary-evaluation step.
- Output shows `PASS: kappa_n (analytic) - (fundamental thm)` instead of `PASS: kappa_n - expected`.
- Script exits 0.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- summary: Reworked the Mathematica overlap law to compare an analytic short form against an indefinite-integral boundary evaluation.
- deviation: none
