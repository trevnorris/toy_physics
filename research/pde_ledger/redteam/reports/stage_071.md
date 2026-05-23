---
unit_id: 071
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

# Audit unit 071 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.txt`

## What the script claims to verify

The scripts audit the canonical tanh-wall branch built from `f(xi) = (1 + tanh(xi))/2`. They evaluate two shape integrals `I_f = ∫(f')^2 dxi = 1/3` and `I_g = ∫(f'')^2 dxi = 4/15` (and the ratio `I_g/I_f = 4/5`), both by direct integration on the real line and via the substitution `t = tanh(xi)` mapping to `[-1,1]`. They then assemble the wall energetics `T_X`, `K_X`, `J_1`, `W_wall`, the local closure `K_m = T_X/ell` with `eta = K_m L / T_X`, the reduced stiffness ratio `kappa = K_X L^2 / T_X = 4 chi_s^2 + (4/5) Lambda_ell^2`, and the reduced wall energy `W_wall = Upsilon_w Lambda_ell^2`. The Mathematica side additionally pins exact closed forms for `T_X`, `K_X`, and `J_1`, which the SymPy side does not check.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero("I_f - 1/3", If - 1/3)` | yes |
| A2 | sympy | 49 | `expect_zero("I_f direct - substitution", If - If_sub)` | partial (both integrals are evaluated, but they encode the same change of variable; primarily a self-consistency probe) |
| A3 | sympy | 50 | `expect_zero("I_g - 4/15", Ig - 4/15)` | yes |
| A4 | sympy | 51 | `expect_zero("I_g direct - substitution", Ig - Ig_sub)` | partial (see A2) |
| A5 | sympy | 52 | `expect_zero("I_g/I_f - 4/5", Ig/If - 4/5)` | yes |
| A6 | sympy | 69 | `expect_zero("eta - L/ell", eta - L/ell)` | NO — tautological (see F1) |
| A7 | sympy | 79-80 | `expect_zero("kappa_reduced - [4 chi_s^2 + 4/5 Lambda_ell^2]", kappa_red - (4*chi_s**2 + 4/5*Lambda_ell**2))` | yes (depends on `subs` pattern-matching succeeding, but output confirms it does) |
| A8 | sympy | 87 | `expect_zero("W_wall_reduced - Upsilon_w Lambda_ell^2", W_red - Upsilon_w*Lambda_ell**2)` | yes (same caveat) |
| B1 | mathematica | 50 | `expectZero["I_f - 1/3", ifDirect - 1/3]` | yes |
| B2 | mathematica | 51 | `expectZero["I_f direct - substitution", ifDirect - ifSub]` | partial |
| B3 | mathematica | 52 | `expectZero["I_g - 4/15", igDirect - 4/15]` | yes |
| B4 | mathematica | 53 | `expectZero["I_g direct - substitution", igDirect - igSub]` | partial |
| B5 | mathematica | 54 | `expectZero["I_g/I_f - 4/5", igDirect/ifDirect - 4/5]` | yes |
| B6 | mathematica | 70 | `expectZero["T_X exact formula", tx - Pi*a^2*ell*hbar^2/(3*m*rhoW)]` | yes (independent closed-form pin) |
| B7 | mathematica | 71-74 | `expectZero["K_X exact formula", kx - 4*Pi*a^2*(5*m^2*cSw^2*ell^2 + hbar^2)/(15*ell*m*rhoW)]` | yes |
| B8 | mathematica | 75 | `expectZero["J_1 exact formula", j1 - rhoW/(3*m*cSw^2)]` | yes |
| B9 | mathematica | 81 | `expectZero["eta - L/ell", eta - L/ell]` | NO — tautological (see F1) |
| B10 | mathematica | 89 | `expectZero["kappa reduced law", kappa - kappaExpected]` | yes (kappaExpected is built from `(m cSw L/hbar)^2` and `(L/ell)^2` as an independent symbolic target) |
| B11 | mathematica | 92 | `expectZero["W_wall reduced law", wwall - wExpected]` | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py:65-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl:77-81`

**What's wrong:**

The SymPy script defines

```python
Km = sp.simplify(Tx / ell)
eta = sp.simplify(Km * L / Tx)
...
expect_zero("eta - L/ell", eta - L / ell)
```

and the Mathematica script does the same:

```mathematica
km = FullSimplify[tx/ell, Assumptions -> $Assumptions];
eta = FullSimplify[km*L/tx, Assumptions -> $Assumptions];
...
expectZero["eta - L/ell", eta - L/ell];
```

By construction, `Km` is set equal to `Tx/ell`, so `eta = Km*L/Tx = (Tx/ell)*L/Tx = L/ell` is an algebraic identity that holds regardless of the value or physical content of `Tx`. The assertion `eta - L/ell = 0` is guaranteed to pass even if `Tx` were replaced with any nonzero expression — e.g. `Tx = 42*xi**3 + cos(L)` — so it does not test the physics of the canonical tanh-wall branch. The script prints in the theorem ledger "eta = L/ell under K_m = T_X/ell", but the `under` clause is the very definition that makes the equality automatic; the check therefore has no failure mode.

**Why this matters:**

The unit's theorem ledger lists `eta = L/ell under K_m = T_X/ell` alongside non-trivial results (`I_f = 1/3`, `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`). A reader scanning the asserted checks would reasonably believe `eta = L/ell` is an independently verified consequence of the wall branch's wall energetics. It is not — it is a definitional rewrite. Leaving the assertion in place gives the audit a false sense of coverage on this row of the ledger.

**Required change:**

Replace the tautological identity check with a substantive verification that actually exercises the closed-form structure of `T_X`. Concretely, verify (a) that `K_m * ell` equals the closed form of `T_X` (so the definition `K_m = T_X/ell` is the correct one given the explicit `T_X = pi a^2 ell I_f hbar^2 / (m rho_w)`), and (b) that `eta * ell` equals `L` after substitution — i.e. the check should pin `K_m` to the actual closed form rather than to `Tx/ell` by fiat.

In SymPy at line 65-69 of `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`, replace the existing block with:

```python
Km = sp.simplify(Tx / ell)
eta = sp.simplify(Km * L / Tx)
print("K_m =", Km)
print("eta =", eta)
# Independent pin: K_m equals the closed form pi a^2 hbar^2 / (3 m rho_w)
Km_expected = sp.pi * a**2 * hbar**2 / (3 * m * rho_w)
expect_zero("K_m - pi a^2 hbar^2 / (3 m rho_w)", Km - Km_expected)
# eta equality reduces to L/ell only when K_m is built from T_X/ell;
# pin eta to L/ell against the closed-form K_m, not against the tautological Tx/ell.
expect_zero("eta - L/ell (from closed-form K_m)", (Km_expected * L / Tx) - L / ell)
```

In Mathematica at lines 77-81 of `mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl`, replace the existing block with:

```mathematica
km = FullSimplify[tx/ell, Assumptions -> $Assumptions];
eta = FullSimplify[km*L/tx, Assumptions -> $Assumptions];
Print["K_m = ", fmt[km]];
Print["eta = ", fmt[eta]];
kmExpected = Pi*a^2*hbar^2/(3*m*rhoW);
expectZero["K_m - pi a^2 hbar^2 / (3 m rhoW)", km - kmExpected];
expectZero["eta - L/ell (from closed-form K_m)", (kmExpected*L/tx) - L/ell];
```

**Verification:**

After Codex applies the fix, the verifier confirms that:
1. In the SymPy script, the line containing `expect_zero("eta - L/ell", eta - L / ell)` no longer exists; in its place the two new `expect_zero` calls (`K_m - pi a^2 hbar^2 / (3 m rho_w)` and `eta - L/ell (from closed-form K_m)`) appear.
2. In the Mathematica script, the line `expectZero["eta - L/ell", eta - L/ell];` no longer exists; in its place the two new `expectZero` calls appear.
3. Both scripts re-run via `redteam exec-sympy 071` and `redteam exec-mathematica 071` and exit 0, with new output lines `K_m - pi a^2 hbar^2 / (3 m rho_w) = 0` (or Mathematica equivalent) and `eta - L/ell (from closed-form K_m) = 0` appearing in the saved transcripts.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally parallel to SymPy (same `f`, `fp`, `fpp`, same `ifDirect`/`ifSub`/`igDirect`/`igSub`, same `T_X`, `K_X`, `J_1`, `W_wall` definitions). However it is not a pure transliteration: it adds three independent closed-form pins not present in SymPy — `T_X exact formula`, `K_X exact formula`, `J_1 exact formula` (Mathematica lines 70-75) — and it constructs `kappaExpected` and `wExpected` from the *symbolic* expressions `4*(m*cSw*L/hbar)^2 + (4/5)*(L/ell)^2` and `(4*rhoW^2*V0^2/(hbar^2*cSw^2))*(L/ell)^2` and checks them by direct difference, whereas SymPy relies on `subs`-based substitution and verifies the substituted form. These are genuinely different verification paths for `kappa` and `W_wall`. I judge this not a `mathematica_transliteration` finding.

Quoted corresponding sections:
- SymPy lines 73-80 (`kappa_red = sp.simplify(kappa.subs({m * c_sw * L / hbar: chi_s, L / ell: Lambda_ell})); expect_zero("kappa_reduced - [4 chi_s^2 + 4/5 Lambda_ell^2]", kappa_red - (4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2))`)
- Mathematica lines 83-89 (`kappa = FullSimplify[kx*L^2/tx, ...]; kappaExpected = FullSimplify[4*(m*cSw*L/hbar)^2 + (4/5)*(L/ell)^2, ...]; expectZero["kappa reduced law", kappa - kappaExpected]`)

The SymPy path probes whether the `subs` pattern matches; the Mathematica path probes whether the closed-form combination matches. They are not the same algebraic manipulation.

## Engine cross-check

Both engines produce identical numeric/symbolic forms:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `f'(xi)` | `1/(2*cosh(xi)**2)` | `Sech[xi]^2/2` |
| `f''(xi)` | `-sinh(xi)/cosh(xi)**3` | `-(Sech[xi]^2*Tanh[xi])` |
| `I_f` | `1/3` | `1/3` |
| `I_g` | `4/15` | `4/15` |
| `T_X` | `pi*a**2*ell*hbar**2/(3*m*rho_w)` | `(a^2*ell*hbar^2*Pi)/(3*m*rhoW)` |
| `K_X` | `4*pi*a**2*(5*c_sw**2*ell**2*m**2 + hbar**2)/(15*ell*m*rho_w)` | `(4*a^2*(hbar^2 + 5*cSw^2*ell^2*m^2)*Pi)/(15*ell*m*rhoW)` |
| `J_1` | `rho_w/(3*c_sw**2*m)` | `rhoW/(3*cSw^2*m)` |
| `W_wall` | `4*L**2*V0**2*rho_w**2/(c_sw**2*ell**2*hbar**2)` | `(4*L^2*rhoW^2*V0^2)/(cSw^2*ell^2*hbar^2)` |
| `kappa` (raw) | `4*L**2*c_sw**2*m**2/hbar**2 + 4*L**2/(5*ell**2)` | `(4*L^2*(ell^(-2) + (5*cSw^2*m^2)/hbar^2))/5` |
| `eta` | `L/ell` | `L/ell` |

Engines agree at the level claimed. No `engine_disagreement` finding.

## Verdict justification

The non-trivial content of the unit — the two shape integrals (`I_f = 1/3`, `I_g = 4/15`), the ratio (`4/5`), the closed forms of `T_X`/`K_X`/`J_1`/`W_wall`, and the reductions `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`, `W_wall = Upsilon_w Lambda_ell^2` — holds up under attack. I tried (a) recomputing `I_f` and `I_g` by hand via `u = tanh(xi)` substitution and got 1/3 and 4/15 respectively; (b) re-deriving `kappa` by direct algebra `K_X L^2 / T_X` and got `4 chi_s^2 + (4/5) Lambda_ell^2` with `chi_s = m c_sw L / hbar`, `Lambda_ell = L/ell`; (c) checking that the SymPy `subs` pattern in line 73-76 is in fact matching the squared form in the printed output and is therefore not silently leaving the expression unreduced (the printed `kappa reduced = 4*Lambda_ell**2/5 + 4*chi_s**2` confirms the substitution succeeded). The Mathematica script independently pins `T_X`, `K_X`, `J_1` to their closed forms, which strengthens the audit.

The single defect is the `eta - L/ell` assertion in both engines: it is an algebraic consequence of the definition `K_m := T_X/ell` and therefore cannot fail. This is a clear `tautological_check` and warrants a fix that pins `K_m` to its closed form before forming `eta`. It is not a `stop_cold` issue — replacing the tautology with a substantive check does not propagate to downstream units' results; only this row of the local ledger gains a real verification.

Verdict: `findings`, one finding, severity medium, fixable in-place.

## Self-test notes

Checked: (i) the substitution-form integrals `I_f_sub` and `I_g_sub` evaluate correctly under `t = tanh(xi)` and reproduce 1/3 and 4/15 — not a derivative-of-constant trap; (ii) the proposed replacement check `Km - pi*a^2*hbar^2/(3*m*rho_w)` exercises the closed form of `Tx` non-trivially (substituting the explicit `Tx = pi*a^2*ell*hbar^2/(3*m*rho_w)` gives `Km = pi*a^2*hbar^2/(3*m*rho_w)`, matching `Km_expected`, so the new `assert_zero` is satisfied for the correct physics and would fail if `I_f != 1/3` or if a factor of `ell` were dropped); (iii) the second proposed check `(Km_expected * L / Tx) - L/ell` evaluates to `(pi*a^2*hbar^2/(3*m*rho_w))*L / (pi*a^2*ell*hbar^2/(3*m*rho_w)) - L/ell = L/ell - L/ell = 0`, so it is true for the correct closed form and would fail under any sign or factor error in `Tx`. The `.py` target lives in `scripts/`, the `.wl` target lives in `mathematica/` — paths checked.
