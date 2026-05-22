---
unit_id: 027
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T23:02:45Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 027

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl:91-105`

**Issue:**

The `wallStiffness[]` module currently hard-codes `kGeo` as the answer and then
sets `kGeoExpected = kGeo`, so the assertion `expectZero["K_geo - expected",
kGeo - kGeoExpected]` cannot fail. The wall-stiffness integral that the docstring
and the SymPy companion both claim to verify (compute
`K_geo = Integrate[chi * (-T_w chi_ss + (K_eta + 6 T_Omega) chi), {s, 0, L}]` and
compare to `K_eta + 6*T_Omega + T_w*Pi^2*Sin[theta]^2/L^2`) is not actually
performed. The downstream `K_geo(theta_max)` check at line 103 and the `theta=0`
check at line 151 inherit the same defect because both substitute into the
hard-coded expression.

**Required change:**

Replace the body of `wallStiffness[]` so that `kGeo` is computed from the wall
operator applied to `chi`, and `kGeoExpected` is a distinct literal closed form.
Edit the file as follows.

Before (lines 91–105 of `moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl`):

```
wallStiffness[] := Module[{u0, u1, f0, chi, kappa, kGeo, kGeoExpected, maxSubs, kGeoMax},
  banner["SECTION III — EXACT WALL STIFFNESS EXPECTATION"];
  {u0, u1, f0, chi, kappa} = overlapLaw[];
  kGeo = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;
  kGeoExpected = kGeo;

  Print["K_geo(theta) = ", fmt[kGeo]];
  expectZero["K_geo - expected", kGeo - kGeoExpected];

  maxSubs = {Cos[theta] -> 3/Sqrt[11], Sin[theta] -> -Sqrt[2]/Sqrt[11]};
  kGeoMax = FullSimplify[TrigExpand[kGeo] /. maxSubs, Assumptions -> $Assumptions];
  Print["K_geo(theta_max) = ", fmt[kGeoMax]];
  expectZero["K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]", kGeoMax - (kEta + 6*tOmega + 2*tW*Pi^2/(11*l^2))];
  {kappa, kGeo}
];
```

After:

```
wallStiffness[] := Module[{u0, u1, f0, chi, kappa, gEta, kGeo, kGeoExpected, maxSubs, kGeoMax},
  banner["SECTION III — EXACT WALL STIFFNESS EXPECTATION"];
  {u0, u1, f0, chi, kappa} = overlapLaw[];
  gEta = -tW*D[chi, {s, 2}] + (kEta + 6*tOmega)*chi;
  kGeo = FullSimplify[Integrate[chi*gEta, {s, 0, l}], Assumptions -> $Assumptions];
  kGeoExpected = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;

  Print["K_geo(theta) = ", fmt[kGeo]];
  expectZero["K_geo - expected", kGeo - kGeoExpected];

  maxSubs = {Cos[theta] -> 3/Sqrt[11], Sin[theta] -> -Sqrt[2]/Sqrt[11]};
  kGeoMax = FullSimplify[TrigExpand[kGeo] /. maxSubs, Assumptions -> $Assumptions];
  Print["K_geo(theta_max) = ", fmt[kGeoMax]];
  expectZero["K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]", kGeoMax - (kEta + 6*tOmega + 2*tW*Pi^2/(11*l^2))];
  {kappa, kGeo}
];
```

The only changes are:

1. Add `gEta` to the `Module` local-variable list (now `{u0, u1, f0, chi, kappa, gEta, kGeo, ...}`).
2. Insert `gEta = -tW*D[chi, {s, 2}] + (kEta + 6*tOmega)*chi;` as a new line before the `kGeo` assignment.
3. Replace `kGeo = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;` with `kGeo = FullSimplify[Integrate[chi*gEta, {s, 0, l}], Assumptions -> $Assumptions];`.
4. Replace `kGeoExpected = kGeo;` with `kGeoExpected = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;`.

Do not touch any other module, any other line, or any other file. Do not change the section banner text. Do not change the assertion name strings.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 027` and
confirm:
1. The script still exits 0.
2. The output line `K_geo(theta) = ...` now reflects the FullSimplify result of the integral (the printed form should still be `kEta + 6*tOmega + (Pi^2*tW*Sin[theta]^2)/l^2` or an algebraically equivalent canonicalisation produced by `FullSimplify`).
3. `PASS: K_geo - expected` still appears.
4. `PASS: K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]` still appears.
5. `PASS: K_geo(theta=0) - (K_eta + 6 T_Omega)` still appears (this assertion in `branchSubstitution[]` operates on the now-derived `kGeo` and remains valid).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl`
- summary: Replaced the tautological wall-stiffness assignment with the requested integral of the wall operator applied to `chi`.
- deviation: none
