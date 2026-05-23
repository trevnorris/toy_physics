---
unit_id: 070
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 070 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.txt`

## What the script claims to verify

The script reduces a GNLS wall-shell setup to closed-form expressions for the kinetic moment `T_X`, support moment `K_X`, dimensionless ratio `kappa = K_X L^2 / T_X`, wall potential weight `W_wall`, and an alternative form `Xi` built from a different grouping (gradient amplitude `g_phi`, moment `I_1`). Three assertions are made: (i) the assembled `kappa` matches the closed form `4 (m c_sw L / hbar)^2 + (I_g / I_f)(L/ell)^2`; (ii) the assembled `W_wall` matches `4 rho_w^2 V0^2 L^2 / (hbar^2 c_sw^2 ell^2)`; (iii) `Xi == W_wall`. Each of these is a non-trivial algebraic consequence of the symbolic definitions for `H_w`, `N_(phi phi)`, `G_(phi phi)`, `T_X`, `K_X`, `J_1`, and `I_1` introduced earlier in the script.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 48 | `simplify(kappa - kappa_expected) == 0` | yes |
| A2 | sympy | 55 | `simplify(Wwall - Wwall_expected) == 0` | yes |
| A3 | sympy | 62 | `simplify(Xi - Wwall) == 0` | yes |
| A4 | mathematica | 52 | `FullSimplify[kappa - kappaExpected] === 0` | yes |
| A5 | mathematica | 63 | `FullSimplify[Wwall - WwallExpected] === 0` | yes |
| A6 | mathematica | 71 | `FullSimplify[Xi - Wwall] === 0` | yes |

All six assertions can in principle fail: each compares an assembled rational expression against either a separately-written closed form (A1, A2, A4, A5) or against a value derived from a *different* grouping of the same primitives (A3, A6 — `Xi` is built from `gphi^2 * I1 * L^2 / Tx`, while `W_wall` is built from `J1 * V0^2 * L^2 / (Tx * ell)`; their equality is not algebraically forced by the definitions until one cancels the `Hw` factor through `I_f`). None is tautological.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:26-77`

**What's wrong:**
The Mathematica script is a line-by-line port of the SymPy script's algebraic choreography rather than an independent re-derivation. The two scripts share an identical sequence of intermediate symbols (`Hw`, `Nphiphi`, `Gphiphi`, `Tx`, `Kx`, `kappa`, `J1`, `Wwall`, `gphi`, `I1`, `Xi`), constructed in the same order, with identical right-hand sides, identical hardcoded "expected" closed forms, and even matching print labels and theorem-ledger banner text. There is no independent derivation pathway in the `.wl` script — it simply re-evaluates the same algebra in a different CAS.

Concrete corresponding sections:

SymPy (lines 32–40):
```
Hw = sp.simplify(m * c_sw**2 / rho_w)
Nphiphi = sp.simplify(4 * pi * a**2 * ell * If)
Gphiphi = sp.simplify(4 * pi * a**2 * Ig / ell)
Tx = sp.simplify(hbar**2 * Nphiphi / (4 * m * rho_w))
Kx = sp.simplify(Hw * Nphiphi + hbar**2 * Gphiphi / (4 * m * rho_w))
kappa = sp.simplify(Kx * L**2 / Tx)
```

Mathematica (lines 33–40):
```
Hw = FullSimplify[m*cSw^2/rhoW, ...];
Nphiphi = FullSimplify[4*Pi*a^2*ell*IfMoment, ...];
Gphiphi = FullSimplify[4*Pi*a^2*Ig/ell, ...];
Tx = FullSimplify[hbar^2*Nphiphi/(4*m*rhoW), ...];
Kx = FullSimplify[Hw*Nphiphi + hbar^2*Gphiphi/(4*m*rhoW), ...];
kappa = FullSimplify[Kx*L^2/Tx, ...];
```

And later, SymPy (lines 47–48) vs Mathematica (lines 48–52), the same hardcoded `kappaExpected` form is reproduced; similarly for `Wwall_expected` (SymPy 52 vs Mathematica 56–59), and for the `gphi/I1/Xi` block (SymPy 57–62 vs Mathematica 65–71).

**Why this matters:**
The two-engine policy exists so that one CAS's algebraic quirk (e.g., a hidden simplification that masks a sign or branch error) cannot pass through to the other. When the second script merely re-evaluates the first script's algebra in a different syntax, both engines necessarily reach the same answer — including the same hidden errors. The agreement reported in the outputs is not independent confirmation.

**Required change:**
Replace the `.wl` script's body with an independent derivation that arrives at the same three claims (`kappa` closed form, `W_wall` closed form, `Xi == W_wall`) through a different algebraic route. Specifically: derive `kappa` and `W_wall` by writing them directly from the primitive parameters (`m, c_sw, ell, L, hbar, I_f, I_g, V0, rho_w`) and verifying against the assembled `K_X L^2 / T_X` and `4 pi a^2 L^2 J_1 V_0^2 / (T_X ell)` forms — i.e., reverse the verification direction. For `Xi == W_wall`, derive each side from primitives independently and confirm the residual is zero, without using the SymPy variable choreography. Do not import the SymPy script's variable names or grouping. The targeted assertions remain:
1. `kappa == 4 (m c_sw L / hbar)^2 + (I_g/I_f) (L/ell)^2`
2. `W_wall == 4 rho_w^2 V_0^2 L^2 / (hbar^2 c_sw^2 ell^2)`
3. `Xi == W_wall` (where `Xi = (V_0/ell)^2 * (N_phiphi / H_w) * L^2 / T_X` and `W_wall = 4 pi a^2 L^2 (I_f/H_w) V_0^2 / (T_X ell)`).

**Verification:**
After Codex edits, the verifier confirms the `.wl` script still has three `expectZero` calls covering the same three claims and exits 0. The reviewer additionally checks that the `.wl` script no longer mirrors the `.py` script's variable choreography line-for-line.

## Independent-derivation check (Mathematica)

Not independent. As documented in F1, the `.wl` script reproduces the SymPy script's variable choreography step-for-step: identical primitive definitions, identical intermediate combinations, identical hardcoded "expected" forms, identical print statements. The Mathematica output's agreement with SymPy's is therefore guaranteed by construction rather than by independent verification.

## Engine cross-check

Both outputs agree at the level claimed:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `N_(phi phi)` | `4*pi*I_f*a**2*ell` | `4*a^2*ell*IfMoment*Pi` |
| `T_X` | `pi*I_f*a**2*ell*hbar**2/(m*rho_w)` | `(a^2*ell*hbar^2*IfMoment*Pi)/(m*rhoW)` |
| `kappa` | `4*L**2*c_sw**2*m**2/hbar**2 + I_g*L**2/(I_f*ell**2)` | `L^2*(Ig/(ell^2*IfMoment) + (4*cSw^2*m^2)/hbar^2)` |
| `W_wall` | `4*L**2*V0**2*rho_w**2/(c_sw**2*ell**2*hbar**2)` | `(4*L^2*rhoW^2*V0^2)/(cSw^2*ell^2*hbar^2)` |
| `Xi - W_wall` | `0` | `0` |

Agreement is bit-for-bit (modulo CAS syntax). However, given the transliteration in F1, this is not independent evidence.

## Verdict justification

The three assertions in the SymPy script (A1–A3) are mathematically sound and non-tautological — each compares an algebraic assembly against a separately-written closed form (A1, A2) or against a different grouping of the same primitives (A3). Symbol assumptions are physically appropriate (all parameters positive real); no branch errors are hiding under `simplify`. The outputs are fresh relative to script mtimes (sympy: script Apr 1, output May 11; mathematica: script Apr 21, output May 11). The only real defect is that the `.wl` script is a transliteration of the `.py` script rather than an independent second-engine derivation, which voids the policy purpose of having two engines. Hence verdict: `findings`, one finding, severity medium. No downstream propagation risk: fixing the Mathematica derivation will not change the values of `kappa`, `W_wall`, or `Xi`; it will only re-establish that the agreement is genuine. Stop-cold is null.

## Self-test notes

I checked: (1) algebraic re-derivation of `kappa` by substituting `Hw`, `Nphiphi`, `Gphiphi`, `Tx`, `Kx` into `Kx L^2 / Tx`, confirming the two terms `4 (m c_sw L/hbar)^2` and `(I_g/I_f)(L/ell)^2` emerge from the `Hw*Nphiphi` and `hbar^2 * Gphiphi/(4 m rho_w)` summands respectively; (2) algebraic re-derivation of `Xi == W_wall` from primitives — `Xi = (V0/ell)^2 (4 pi a^2 ell I_f / Hw) L^2 / Tx`, `W_wall = 4 pi a^2 L^2 (I_f/Hw) V0^2/(Tx ell)`, and the ratio reduces to 1 algebraically — so A3/A6 are genuine, non-tautological coincidences; (3) symbol positivity assumptions are consistent with the physical setup (lengths, densities, sound speeds, moments all positive). No symmetry/parity traps apply (no integrals in this unit). No variable-independence traps (no derivative-zero hazards).
