---
unit_id: 252
batch: VIII.1
created_at: 2026-06-03T15:32:13Z
findings_count: 4
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
applied_at: 2026-06-03T19:13:59Z
findings_applied: 4
findings_blocked: 0
---

# Codex directive — unit 252

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## F1 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl`
(NEW file. `.wl` files live in `mathematica/`. The committed output will land at `mathematica/output/<same stem>.txt`.)

**Issue:** Stage 252 has no Mathematica engine, but it is not status-only and not a carry-forward; every claim is independently derivable in Mathematica. Author a NEW independent-route `.wl` — NOT a transliteration of the `.py`. Use native Mathematica primitives and a *different decomposition* than the SymPy script wherever practical. Specifically:
- Derive the partition fractions by simplifying `E_lat/E_exp` and `E_vac/E_exp` directly with `Together`/`Cancel`/`FullSimplify` (do NOT mirror the `.py`'s `.subs(I3, rV*I2)` choreography — instead substitute `I3 -> rV I2` once and `Cancel`).
- Obtain the drift law by `D[flat, rV]` followed by `Together`, comparing to the paper's closed form via `expectZero`.
- Take endpoint limits with `Limit[flat, rV -> 0]` and `Limit[flat, rV -> Infinity]`.
- Compute the exponential shape integrals with `Integrate[D[Vin Exp[s t], {t, k}]^2, {t, 0, T}]` (k=1,2,3), with `Assumptions -> {s > 0, T > 0, Vin > 0}`.
- Verify the safe-edge energy theorem AND the safe-edge rate identity by substituting the safe equality `Gamma3 sc^3 + Gamma5 sc^5 -> mueta (s0^2 - sc^2)` (use a replacement rule on the unexpanded combo, e.g. multiply/divide so the literal `Gamma3 sc^3 + Gamma5 sc^5` is present before applying the rule, then `Simplify`).
- Verify the 3:1 iff via the route in F3 (relate `flat[sc] - 3/4` numerator to `Gamma3lat + Gamma5lat sc^2 - 3(Gamma3vac + Gamma5vac sc^2)`), NOT by re-expanding a polynomial against its own expansion.
- Verify the φ-family reduction `flat -> phi` under `{Gamma3lat -> phi G3T, Gamma3vac -> (1-phi) G3T, Gamma5lat -> phi G5T, Gamma5vac -> (1-phi) G5T}`.
- Reproduce the Session-IV numeric benchmark with `N[...]` and `expectApprox`, deriving the 1/4, 3/4 weights from the φ=3/4 family (see F4), not as bare literals.

Use the project's standard Mathematica check helpers. Preemptively strip `ConditionalExpression[0, ...]` from any `Solve`/`Reduce`/`Limit` result before comparison (project idiom). Each helper call must be a real check that aborts with `Exit[1]` on failure.

**Symbol assumptions for the .wl:** declare `Gamma3vac, Gamma3lat, Gamma5vac, Gamma5lat >= 0`; `I2, rV, s, T, Vin, sc, s0, mueta > 0`. Do not over-constrain (e.g., do not assume `flat < 1` — that must emerge, not be assumed).

**Claim manifest** (each must be an independent check; ≥9 to match the .py's 7-item docstring with the safe-edge split counted):
- **M1** (energy ledger + fractions): with `Evac = Gamma3vac I2 + Gamma5vac I3`, `Elat = Gamma3lat I2 + Gamma5lat I3`, `Eexp = Evac + Elat`, and `I3 -> rV I2`:
  `Cancel[Evac/Eexp] == (Gamma3vac + Gamma5vac rV)/((Gamma3vac+Gamma3lat) + (Gamma5vac+Gamma5lat) rV)` and the lattice analogue, and `fvac + flat == 1`.
- **M2** (drift law): `Together[D[flat, rV]] == (Gamma5lat (Gamma3vac+Gamma3lat... )...)` — i.e. `(Gamma5lat Gamma3vac - Gamma3lat Gamma5vac)/((Gamma3vac+Gamma3lat) + (Gamma5vac+Gamma5lat) rV)^2`.
- **M3** (endpoint limits): `Limit[flat, rV->0] == Gamma3lat/(Gamma3vac+Gamma3lat)`, `Limit[flat, rV->Infinity] == Gamma5lat/(Gamma5vac+Gamma5lat)`.
- **M4** (exponential shape integrals): `I2exp == Vin^2 s^3 (Exp[2 s T]-1)/2`, `I3exp == s^2 I2exp`, `I3exp/I2exp == s^2`.
- **M5** (event-equivalent rates): `Evac /. {I2->I2exp, I3->I3exp}` equals `(Gamma3vac s^2 + Gamma5vac s^4) I1exp`, lattice analogue, where `I1exp == Vin^2 s (Exp[2 s T]-1)/2`.
- **M6** (safe-edge energy theorem): `Eexp` at `s->sc, T->1/sc` (so `Exp[2]` appears) equals `Vin^2 (E^2-1)(Gamma3 sc^3+Gamma5 sc^5)/2`, and under the safe rule → `Vin^2 (E^2-1) mueta (s0^2-sc^2)/2`.
- **M7** (safe-edge rate identity): `(Gamma3 sc^3+Gamma5 sc^5)/sc == Gamma3 sc^2+Gamma5 sc^4`, AND under the safe rule equals `mueta (s0^2-sc^2)/sc`. (This is the substance F2 fixes on the SymPy side; the .wl must carry the safe rule, not assert a self-identity.)
- **M8** (3:1 coefficient surface iff): numerator of `Together[flat /. rV->sc^2] - 3/4` is a positive multiple of `Gamma3lat + Gamma5lat sc^2 - 3(Gamma3vac + Gamma5vac sc^2)` (pin the exact factor in-script); and the φ-family check `flat -> phi`.
- **M9** (Session-IV benchmark): `N` reproduction of `sc≈0.5489386551`, `gamma_eff^eq,safe≈87.26925235`, prefactor `≈153.0353549`, `Vin_match≈8.21771260e-3`, `Evac≈0.00258365`, `Elat≈0.00775095`, with the 1/4, 3/4 weights taken from the φ=3/4 family.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 252` and confirm the new `.wl` appears, contains ≥9 substantive checks (M1–M9), and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M9, including partition fractions, drift and endpoint limits, exponential integrals, safe-edge identities, phi-family reduction, and the Session-IV benchmark.
- deviation: none

---

## F2 — tautological_check

**Target:** `scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py:141`

**Issue:** L122 defines `gamma_safe_eq = sp.simplify(mu_eta * (s0**2 - sc**2) / sc)`, and L141 asserts `sp.simplify(gamma_safe_eq - mu_eta * (s0**2 - sc**2) / sc) == 0` — a vacuous `X − X == 0`. The paper deliverable is the bridge `Gamma3 sc^2 + Gamma5 sc^4 = mu_eta (s0^2 - sc^2)/sc` *under the safe equality* `Gamma3 sc^3 + Gamma5 sc^5 = mu_eta (s0^2 - sc^2)` (card lines 359-365). The safe equality (`safe_eq`, L114) is defined but never used in any assertion.

**Required change:**
Replace L141. Carry the safe equality through `gamma_safe_eq_expected` (which is `safe_combo/sc = Gamma3 sc^2 + Gamma5 sc^4`, L121). Robust substitution form (do NOT rely on `simplify` having preserved the `safe_combo` subexpression — substitute on the explicit combo built fresh):

Before (L141):
```python
    assert sp.simplify(gamma_safe_eq - mu_eta * (s0**2 - sc**2) / sc) == 0
```
After:
```python
    # carry the Stage-251 safe equality through the rate, so the bridge identity
    # gamma_eff,safe^eq = mu_eta (s0^2 - s_c^2)/s_c is actually exercised (not X-X).
    gamma_safe_eq_bridged = sp.simplify(
        ((G3 * sc**3 + G5 * sc**5) / sc).subs(
            G3 * sc**3 + G5 * sc**5, mu_eta * (s0**2 - sc**2)
        )
    )
    assert sp.simplify(gamma_safe_eq_bridged - mu_eta * (s0**2 - sc**2) / sc) == 0
```
If the `.subs(G3*sc**3 + G5*sc**5, ...)` pattern does not fire (SymPy may auto-distribute the `/sc`), use the algebraically equivalent route: build `lhs = G3*sc**2 + G5*sc**4`, substitute the safe equality solved for one coefficient, or substitute `mu_eta*(s0**2-sc**2) -> G3*sc**3+G5*sc**5` in reverse and compare `lhs*sc` to `G3*sc**3+G5*sc**5`. The assertion MUST reference `mu_eta`, `s0`, AND the Γ-coefficients so it can fail if the bridge is wrong. Keep L140 unchanged.

**Self-test:** `gamma_safe_eq_bridged` depends on `mu_eta, s0, sc` after the safe-rule fires; with the safe rule NOT applied the residual would be `(G3 sc^3+G5 sc^5)/sc - mu_eta(s0^2-sc^2)/sc ≠ 0`, so the check is non-vacuous. ✓

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py`
- summary: Replaced the tautological safe-rate assertion with a check that carries the Stage-251 safe equality through the Gamma-coefficient rate expression.
- deviation: none

---

## F3 — insufficient_verification

**Target:** `scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py:148,164`

**Issue:** L164 asserts `sp.simplify(split_surface) == G3l + G5l*sc**2 - 3*G3v - 3*G5v*sc**2`, where `split_surface` (L148) is the `expand` of that same expression — `expand(X) == expand(X)`, verifying nothing. The paper deliverable (6) is the iff `f_lat(s_c)=3/4 ⟺ Gamma3lat + Gamma5lat sc^2 = 3(Gamma3vac + Gamma5vac sc^2)` (card lines 388-395). The check must tie the partition fraction `fl_sc` (already defined at L124 as `fl.subs(rV, sc**2)`) to the surface.

**Required change:**
Replace the L164 assertion (keep `split_surface` defined for the print at L161, or redefine `surface` cleanly). Add an iff check linking `fl_sc - 3/4` to the surface:

Before (L164):
```python
    assert sp.simplify(split_surface) == G3l + G5l * sc**2 - 3 * G3v - 3 * G5v * sc**2
```
After:
```python
    # the 3:1 split is the iff f_lat(s_c)=3/4  <=>  surface = 0; tie the fraction to it.
    surface = sp.expand((G3l + G5l * sc**2) - 3 * (G3v + G5v * sc**2))
    resid = sp.together(fl_sc - sp.Rational(3, 4))
    num = sp.numer(sp.cancel(resid))
    # numerator of f_lat(s_c) - 3/4 over common denom 4*(Gamma3 + Gamma5 sc^2):
    #   4*(Gamma3lat + Gamma5lat sc^2) - 3*(Gamma3 + Gamma5 sc^2)
    #   = (Gamma3lat + Gamma5lat sc^2) - 3*(Gamma3vac + Gamma5vac sc^2) = surface
    assert sp.simplify(num - surface) == 0
```
**Self-test of the factor:** `fl_sc = (G3l + G5l sc^2)/(G3 + G5 sc^2)` with `G3=G3v+G3l`, `G5=G5v+G5l`. Then `fl_sc - 3/4 = [4(G3l+G5l sc^2) - 3(G3+G5 sc^2)] / [4(G3+G5 sc^2)]`. Numerator `= 4G3l+4G5l sc^2 - 3(G3v+G3l) - 3(G5v+G5l)sc^2 = (G3l+G5l sc^2) - 3(G3v+G5v sc^2) = surface`. So after `cancel`, `numer(resid)` equals `surface` exactly (factor 1, denom common factor 4 cancels into the rational). If SymPy keeps the denominator's 4 inside the numerator instead (giving `4*surface`), the assert will fail — in that case change the line to `assert sp.simplify(num - 4*surface) == 0` OR normalize via `assert sp.simplify(sp.factor(num)/sp.factor(surface)) == sp.simplify(sp.factor(num)/sp.factor(surface))`-style ratio-is-constant check. Codex: run it, observe the actual `num`, and pick the matching factor; the load-bearing requirement is that `num` is a *nonzero constant multiple of `surface`* (not that it is a self-expansion). Do NOT let the assertion degenerate back into comparing `surface` to its own `expand`.

**Trivial-case pre-check:** set `G3l=3*G3v, G5l=3*G5v` → `fl_sc = (3G3v+3G5v sc^2)/(4G3v+4G5v sc^2) = 3/4`, so `resid=0` and `num=0=surface`. Confirms the iff and that the check fires correctly. ✓

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py`
- summary: Replaced the self-expansion assertion with a check tying the numerator of `f_lat(s_c) - 3/4` to the microscopic 3:1 coefficient surface.
- deviation: none

---

## F4 — hardcoded_result

**Target:** `scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py:174-175`

**Issue:** `frac_v_num = 0.25` / `frac_l_num = 0.75` are typed-in literals (L174-175) and drive the rate/energy products (L180-181, L184-185). The paper obtains these specifically from the φ=3/4 speed-independent family (`f_vac=1−φ`, `f_lat=φ`; card §5, lines 458-465). The script already proves `fl_phi == phi` at L165 but never uses it to source the numeric weights, so the calibration step is decoupled arithmetic.

**Required change:**
Derive `f_vac` on the φ-family and source the weights from φ=3/4. After L165 (which proves `fl_phi - phi == 0`), add the vacuum analogue, then set the numeric fractions from φ=3/4:

After L165, add:
```python
    fv_phi = sp.simplify(
        fv.subs(
            {
                G3l: phi * G3T,
                G3v: (1 - phi) * G3T,
                G5l: phi * G5T,
                G5v: (1 - phi) * G5T,
            }
        )
    )
    assert sp.simplify(fv_phi - (1 - phi)) == 0  # vacuum fraction on the phi-family
    phi_val = sp.Rational(3, 4)  # the speed-independent 3:1 microscopic family
```
Then change L174-175:
```python
    frac_v_num = 0.25
    frac_l_num = 0.75
```
to:
```python
    frac_v_num = float(1 - phi_val)  # = f_vac on the phi=3/4 family
    frac_l_num = float(phi_val)      # = f_lat on the phi=3/4 family
```
Leave `t_cross_num`, `s0_num`, `E_diss_num` and ALL the literal assert targets at L197-206 unchanged — those are paper-stated benchmark values (card §5/§5.1).

**Self-test:** `fv_phi` numerator `= (1-phi)G3T + (1-phi)G5T rV = (1-phi)(G3T+G5T rV)`, denom `= G3T+G5T rV`, so `fv_phi = 1-phi`; with `phi=3/4`, `frac_v_num=0.25`, `frac_l_num=0.75`, so all existing numeric asserts (L197-206) still hold. ✓

**Verification command:**
After Codex applies F2/F3/F4, the verifier runs `redteam exec-sympy 252` and confirms: the script exits 0; L141 now references `mu_eta`/`s0`/Γ (F2); the 3:1 check references `fl_sc`/`surface` (F3); `frac_v_num`/`frac_l_num` derive from `phi_val` (F4); and the assertion count has not dropped.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py`
- summary: Added the vacuum phi-family fraction check and sourced the benchmark vacuum/lattice weights from `phi_val = 3/4`.
- deviation: none
