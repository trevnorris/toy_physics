---
unit_id: 070
batch: III.3
created_at: 2026-05-22T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-22T20:06:25-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 070

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:26-77`

**Issue:**
The `.wl` script reproduces the SymPy script's variable choreography line-for-line (same primitive definitions, same intermediate symbols `Hw, Nphiphi, Gphiphi, Tx, Kx, kappa, J1, Wwall, gphi, I1, Xi`, same hardcoded "expected" closed forms, same evaluation order). This means both engines are evaluating the same algebra and any hidden simplification error would propagate identically through both. The Mathematica script must derive the three claims through an independent algebraic route, not echo the SymPy script.

**Required change:**

Replace lines 26–77 of the `.wl` script with an independent derivation that lands the same three assertions (`kappa` closed form, `W_wall` closed form, `Xi == W_wall`) via a different algebraic path. Concretely:

1. Keep the `banner`, `pass`, `fail`, `expectZero`, `fmt` helper definitions (lines 1–24) and the final `Exit[0]` unchanged.

2. After the `banner` call at the top, declare symbols and assumptions as currently done, but reverse the verification direction so the script is not a re-evaluation of the SymPy assembly:

   - Define the *closed forms* `kappaClosed = 4*(m*cSw*L/hbar)^2 + (Ig/IfMoment)*(L/ell)^2` and `WwallClosed = 4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2)` directly from primitives.
   - Define `H_w`, `N_(phi phi)`, `G_(phi phi)`, `T_X`, `K_X` from their physical definitions, then compute `kappaAssembled = Kx*L^2/Tx` and `WwallAssembled = 4*Pi*a^2*L^2*(IfMoment/Hw)*V0^2/(Tx*ell)`.
   - Use `expectZero["kappa: assembled - closed", kappaAssembled - kappaClosed]` (i.e., compare assembled to closed, just as now).
   - Use `expectZero["W_wall: assembled - closed", WwallAssembled - WwallClosed]`.
   - For the `Xi == W_wall` check, derive `Xi` independently as follows: instead of using the SymPy grouping `gphi^2 * I1 * L^2 / Tx`, expand `Xi` directly in primitives as `Xi = (V0/ell)^2 * (4*Pi*a^2*ell*IfMoment*rhoW/(m*cSw^2)) * L^2 / Tx`, then verify `expectZero["Xi - W_wall", Xi - WwallAssembled]`. This is algebraically equivalent but does not introduce the `J1 = IfMoment/Hw` intermediate symbol shared with the SymPy script — the Mathematica script substitutes the `1/Hw = rhoW/(m*cSw^2)` factor inline rather than re-using SymPy's variable name `J_1`.

3. To make the independence more concrete, drop the intermediate symbols `J1`, `gphi`, `I1` from the Mathematica script (they are SymPy-mirror names) and inline their definitions where used. Keep only `Hw`, `Nphiphi`, `Gphiphi`, `Tx`, `Kx` as named intermediates because those are physical primitives whose names are conventional rather than mirrored.

4. Retain the three `expectZero` calls (one each for `kappa`, `W_wall`, and `Xi - W_wall`) so the assertion count and physical content match. The banner / theorem-ledger text after line 73 may remain identical to the SymPy ledger since it is descriptive only.

Concrete before/after for the affected `Xi`/`Wwall` block (current Mathematica lines 54–71):

Before (lines 54–71):
```
J1 = FullSimplify[IfMoment/Hw, Assumptions -> $Assumptions];
Wwall = FullSimplify[4*Pi*a^2*L^2*J1*V0^2/(Tx*ell), Assumptions -> $Assumptions];
WwallExpected = FullSimplify[
  4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2),
  Assumptions -> $Assumptions
];

Print["J_1 = ", fmt[J1]];
Print["W_wall = ", fmt[Wwall]];
expectZero["W_wall - expected", Wwall - WwallExpected];

gphi = FullSimplify[V0/ell, Assumptions -> $Assumptions];
I1 = FullSimplify[Nphiphi/Hw, Assumptions -> $Assumptions];
Xi = FullSimplify[gphi^2*I1*L^2/Tx, Assumptions -> $Assumptions];

Print["g_phi = ", fmt[gphi]];
Print["Xi = ", fmt[Xi]];
expectZero["Xi - W_wall", Xi - Wwall];
```

After:
```
Wwall = FullSimplify[
  4*Pi*a^2*L^2*(IfMoment*rhoW/(m*cSw^2))*V0^2/(Tx*ell),
  Assumptions -> $Assumptions
];
WwallExpected = FullSimplify[
  4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2),
  Assumptions -> $Assumptions
];

Print["W_wall = ", fmt[Wwall]];
expectZero["W_wall - expected", Wwall - WwallExpected];

Xi = FullSimplify[
  (V0/ell)^2 * (4*Pi*a^2*ell*IfMoment*rhoW/(m*cSw^2)) * L^2 / Tx,
  Assumptions -> $Assumptions
];

Print["Xi = ", fmt[Xi]];
expectZero["Xi - W_wall", Xi - Wwall];
```

The key difference: `Wwall` and `Xi` are now built by inlining `1/Hw -> rhoW/(m*cSw^2)` rather than by introducing `J1` and `I1` as named intermediates. This breaks the structural mirror with the SymPy script while preserving the physics and all three assertions.

**Claim manifest** (for the rewritten Mathematica script — three claims, each must be verified via `expectZero`):
- M1: `K_X L^2 / T_X == 4 (m c_sw L / hbar)^2 + (I_g / I_f) (L / ell)^2`, where `T_X = hbar^2 N_(phi phi)/(4 m rho_w)`, `K_X = H_w N_(phi phi) + hbar^2 G_(phi phi)/(4 m rho_w)`, `N_(phi phi) = 4 pi a^2 ell I_f`, `G_(phi phi) = 4 pi a^2 I_g / ell`, `H_w = m c_sw^2/rho_w`.
- M2: `4 pi a^2 L^2 (I_f / H_w) V_0^2 / (T_X ell) == 4 rho_w^2 V_0^2 L^2 / (hbar^2 c_sw^2 ell^2)`.
- M3: `(V_0/ell)^2 (N_(phi phi)/H_w) L^2 / T_X == 4 pi a^2 L^2 (I_f/H_w) V_0^2/(T_X ell)` (i.e., `Xi == W_wall`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 070` and confirm that (a) the three `expectZero` checks (`kappa - expected`, `W_wall - expected`, `Xi - W_wall`) still appear and all PASS, and (b) the script exits 0. The reviewer additionally inspects the diff to confirm that the `.wl` script no longer mirrors the SymPy variable choreography line-for-line (in particular, that `J1`, `gphi`, `I1` are no longer named intermediates and that `Wwall`/`Xi` are constructed via inlined `1/Hw` substitution).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl`
- summary: Rewrote the Mathematica wall-shell checks to compare assembled quantities against primitive closed forms and inline the wall response factors for `W_wall` and `Xi`.
- deviation: none
