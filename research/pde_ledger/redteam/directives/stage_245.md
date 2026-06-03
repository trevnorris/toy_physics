---
unit_id: 245
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 5
stop_cold: null
applied: true
applied_at: 2026-06-03T10:01:53-06:00
findings_applied: 5
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 245

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

There are NO paper_misalignment findings here — every paper/script value agrees. All fixes are script-side.

## F1 — insufficient_verification (vacuous support-blindness self-test)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:156-178` (Section 6)

**Issue:** Section 6 differentiates `U_sol, V_sol, eps_eta, R_ratio, D_UV, Delta_Keta, Delta_mu` w.r.t. fresh symbols `Lam, varrho`. None of those expressions contain `Lam` or `varrho`, so `sp.diff(obj, Lam)` is identically 0 by construction and the asserts cannot fail. The paper claims support-blindness *conditional* on `f_U` and `chi_UV` being orbit-side data (card line 107). The check must (a) carry `f_U`, `chi_UV` as orbit-side functions independent of `Lam, varrho`, and (b) include a positive control that a support-contaminated forcing genuinely produces a nonzero derivative, proving the check discriminates.

**Required change:**
Keep the existing zero-asserts on the seven objects (lines 172-178) — those remain correct under the orbit-side reading. ADD, after the existing loop (after line 178), a positive-control block. Concretely:

```python
    # Positive control: a support-contaminated forcing must NOT be support-blind.
    # If f_U secretly carried a support coordinate, U would depend on Lam, so the
    # support-blindness check above must be capable of detecting it.
    f_U_bad = f_U + Lam            # orbit forcing contaminated by support coordinate Lam
    U_bad = sp.simplify(a_V * f_U_bad / Delta_UV)
    dU_bad_dLam = sp.simplify(sp.diff(U_bad, Lam))
    print("control d/dLam U_bad  =", dU_bad_dLam)
    assert dU_bad_dLam != 0
    assert sp.simplify(dU_bad_dLam - a_V / Delta_UV) == 0
```

This shows the differentiation route registers a real support dependence when one exists (so the zero results above are meaningful, not vacuous). `Lam` is already declared at line 160.

**Self-test (already performed by auditor):** `diff(f_U + Lam, Lam) = 1`, so `diff(a_V*(f_U+Lam)/Delta_UV, Lam) = a_V/Delta_UV != 0`. The control assert fires correctly; the equality assert pins the exact value. The pre-existing zero checks are unaffected because the real objects still contain no `Lam`/`varrho`.

**Verification command:**
The verifier runs `redteam exec-sympy 245`; the new `control d/dLam U_bad` print and the two control asserts must be present and the script must exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py`
- summary: Added a support-contaminated forcing positive control after the existing support-blindness derivative checks.
- deviation: none

## F2 — insufficient_verification (R_target not derived from selected-branch identity; T^2/eps_eta checks tautological)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:90-109` (Section 3)

**Issue:** `R_ratio` is written directly as the answer (line 99) and `R_exact_check` (line 100) round-trips it to `1-1=0`, so A12 cannot fail. The paper (eq:app-part08-stage245-selected-identity → eq:app-part08-stage245-rtarget-finite) requires `R_target/R_ref` to FOLLOW from `R_target*T^2 = Lambda_0(1-eps_eta)`. Derive it from that identity instead.

**Required change:**
Within Section 3, after the existing `T2`, `eps_eta` definitions (lines 97-98), ADD a derivation of `R_target/R_ref` from the selected-branch identity and assert it equals the paper's boxed form. Introduce `Lambda_0` as a positive symbol:

```python
    Lambda_0 = sp.symbols("Lambda_0", positive=True, real=True)
    # Selected-branch identity: R_target * T^2 = Lambda_0 * (1 - eps_eta).
    R_target_from_id     = Lambda_0 * (1 - eps_eta) / T2
    R_target_ref_from_id = Lambda_0 * (1 - eps_eta_ref) / T2_ref
    R_ratio_derived = sp.simplify(R_target_from_id / R_target_ref_from_id)
    R_ratio_paper = sp.simplify(((1 - eps_eta_ref * sp.exp(V_sol)) / (1 - eps_eta_ref)) * sp.exp(-U_sol))
    print("R_target/R_ref (from identity) =", R_ratio_derived)
    assert sp.simplify(R_ratio_derived - R_ratio_paper) == 0
```

Keep the existing A10/A11 asserts (lines 107-108) — they are weak but harmless. The new assert is the load-bearing one: it proves the implication from the identity. `Lambda_0` cancels in the ratio (auditor-verified), so no new free constant enters the result.

**Self-test (already performed by auditor):** `R_target_from_id / R_target_ref_from_id = [(1-eps_ref*exp(V_sol))/(1-eps_ref)] * (T2_ref/T2) = [(1-eps_ref*exp(V_sol))/(1-eps_ref)] * exp(-U_sol)`, which equals `R_ratio_paper`. So `R_ratio_derived - R_ratio_paper` simplifies to 0. Non-tautological because `R_ratio_derived` is built from the identity, not from `R_ratio`.

**Verification command:**
`redteam exec-sympy 245`; the `R_target/R_ref (from identity)` print and the `assert sp.simplify(R_ratio_derived - R_ratio_paper) == 0` must appear, script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py`
- summary: Added the selected-branch identity derivation of `R_target/R_ref` and asserted it matches the paper form.
- deviation: none

## F3 — tautological_check (Session-I rebuild is an inverse round-trip; session anchors print-only)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:200-218` (Section 7)

**Issue:** `U_rebuilt`/`V_rebuilt` (lines 200-201) invert the very formulas just used to infer `chi_UV_obs`, `fU_obs`, so the asserts at lines 217-218 are tautological (`x = g(g^{-1}(x))`). Meanwhile the genuine session numbers `R_target/R_ref ~ 0.87984149` and `R_1^nr ~ -0.12762119` are printed (lines 205, 214) but never asserted.

**Required change:**
1. Demote the two tautological asserts at lines 217-218 to comments or remove them (keep the `U_rebuilt`/`V_rebuilt` computations and prints — they document the implied parameters). Replace with assertions against the session-recorded values:

   Before (lines 217-218):
   ```python
       assert abs(float(U_rebuilt - U_obs)) < 1e-12
       assert abs(float(V_rebuilt - V_obs)) < 1e-12
   ```
   After:
   ```python
       # Round-trip through the inverses is identically exact and is not a physics check;
       # assert instead against the independently-recorded Session-I numbers.
       assert abs(float(R_ratio_obs) - 0.87984149) < 5e-9
       assert abs(float(R1_obs) - (-0.12762119)) < 5e-9
   ```

Keep the existing `eps_eta_obs` assert (line 216) unchanged — it is a correct numeric anchor.

**Self-test (already performed by auditor):** Script already prints `R_target/R_ref(obs) = 0.8798414919352429` (within 5e-9 of 0.87984149) and `R1_obs = -0.1276211900000000` (within 5e-9 of -0.12762119), so both new asserts pass. `R_ratio_obs` and `R1_obs` are already defined at lines 192 and 202, so they are in scope.

**Verification command:**
`redteam exec-sympy 245`; the two new numeric asserts against `0.87984149` and `-0.12762119` must appear, the `< 1e-12` round-trip asserts removed/demoted, script exits 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py`
- summary: Replaced the inverse round-trip assertions with direct Session-I numeric assertions for `R_target/R_ref` and `R1_obs`.
- deviation: none

## F4 — insufficient_verification (drain nonnegativity not asserted)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:78-88` (Section 2)

**Issue:** The paper's deliverable is that the drain is nonnegative "even when U and V have opposite signs" (card line 73, notes §5). The script checks only the closed form (A9, line 88), never the sign. Add a concrete opposite-sign (chi_UV < 0) admissible-point check that the drain is strictly positive.

**Required change:**
After the existing assert at line 88, ADD:

```python
    # Nonnegativity: the drain stays > 0 even on an opposite-sign branch (chi_UV < 0),
    # which is the Session-I branch (U > 0, V < 0).  This is the stated physical claim.
    drain_neg_chi = D_expected.subs({a_U: sp.Float("2.5"), a_V: sp.Float("3.0"),
                                     chi_UV: sp.Float("-0.76"), f_U: sp.Float("0.33")})
    print("D_UV at opposite-sign point =", sp.N(drain_neg_chi, 16))
    assert float(drain_neg_chi) > 0
```

**Self-test (already performed by auditor):** With those values `Delta_UV = 2.5*3.0 - 0.76**2 = 6.9224`, numerator `0.76**2 * 3.0 * 0.33**2 = 0.18871...`, so `D = 0.18871/6.9224**2 ≈ 0.00394 > 0`. `chi_UV = -0.76 < 0` and `Delta_UV = 6.9224 > 0` (admissible). Assert fires True.

**Verification command:**
`redteam exec-sympy 245`; the `D_UV at opposite-sign point` print and `assert float(drain_neg_chi) > 0` must appear, script exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py`
- summary: Added the concrete opposite-sign admissible-point drain positivity check.
- deviation: none

## F5 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.wl` (CREATE — new file in `mathematica/`)

**Issue:** Stage 245 is single-engine (SymPy only) but every claim is independently verifiable in Mathematica. Create an INDEPENDENT second engine — NOT a transliteration of the `.py`. Use native primitives (`Solve`, `Det`, `Series`/`Normal`/`Coefficient`, `D`, `N`, `Simplify`, `FullSimplify`) and a DIFFERENT decomposition where it matters: derive `R_target` from the selected-branch identity (M4), build the support-blindness check WITH a positive control (M7), and assert the drain's nonnegativity (M3). Use the project's standard helper idioms (`expectZero`, `expectTrue`, `expectApprox`) — strip any `ConditionalExpression[0, ...]` from `Solve`/`Reduce` outputs, and prefer `1/expr == 0` style pole checks if any arise (none expected here).

**Claim manifest** (each must be an independent `expectZero`/`expectTrue`/`expectApprox`, NOT echoing the `.py` algebra):

- **M1 — stationary packet.** From `F_nr = (1/2) a_U U^2 + (1/2) a_V V^2 - chi_UV U V - f_U U`, solve `D[F_nr,U]==0, D[F_nr,V]==0` via `Solve` for `{U,V}`. `expectZero` that `U_sol - a_V f_U/(a_U a_V - chi_UV^2) == 0` and `V_sol - chi_UV f_U/(a_U a_V - chi_UV^2) == 0`. Also `expectZero[V_sol/U_sol - chi_UV/a_V]` and `expectZero[Det[Hessian] - (a_U a_V - chi_UV^2)]` where Hessian `= {{a_U,-chi_UV},{-chi_UV,a_V}}`. Limit cases: `expectZero` of `U_sol /. f_U->0`, `V_sol /. f_U->0`, `V_sol /. chi_UV->0`, and `(U_sol /. chi_UV->0) - f_U/a_U`.
- **M2 — drain form.** `expectZero[chi_UV U_sol V_sol - chi_UV^2 a_V f_U^2/(a_U a_V - chi_UV^2)^2]`.
- **M3 — drain nonnegativity (opposite-sign).** Substitute `{a_U->2.5, a_V->3.0, chi_UV->-0.76, f_U->0.33}` into the drain closed form; `expectTrue` the numeric value `> 0`.
- **M4 — finite physical compiler from the selected-branch identity.** Define `T2 = T2ref Exp[U_sol]`, `epseta = epsref Exp[V_sol]`. Then define `Rtarget = Lambda0 (1 - epseta)/T2`, `Rtargetref = Lambda0 (1 - epsref)/T2ref` (Lambda0 a positive symbol), and `expectZero[FullSimplify[Rtarget/Rtargetref - ((1 - epsref Exp[V_sol])/(1 - epsref)) Exp[-U_sol]]]`. Also `expectZero[T2/T2ref - Exp[U_sol]]` and `expectZero[epseta/epsref - Exp[V_sol]]`. (Deriving R_target from the identity is the route the .py omits — this is what makes the .wl independent rather than a port.)
- **M5 — dependent microscopic correction.** `expectZero` of the vector `{0, -V_sol, U_sol - V_sol} - {0, -chi_UV f_U/(a_U a_V - chi_UV^2), (a_V - chi_UV) f_U/(a_U a_V - chi_UV^2)}` componentwise.
- **M6 — first-order packet.** Set `U->eps u1`, `V->eps v1`. Using `Series[..., {eps,0,1}]` // `Normal` // `Coefficient[#, eps]&` (native, different from the .py's custom `coeff_linear`): `expectZero[Coefficient[Normal[Series[Log[Exp[eps u1]], {eps,0,1}]], eps] - u1]`; `expectZero` for `Log[(1 - epsref Exp[eps v1])/(1 - epsref)]` coefficient `- (-epsref v1/(1-epsref))`; `expectZero` for `Log[((1 - epsref Exp[eps v1])/(1 - epsref)) Exp[-eps u1]]` coefficient `- (-u1 - epsref v1/(1-epsref))`; and the consistency `expectZero[R1 + Xi1 - (R1plusXi1)]`.
- **M7 — support-blindness WITH positive control.** With `Lam, varrho` symbols: `expectZero[D[U_sol, Lam]]`, `expectZero[D[U_sol, varrho]]`, and likewise for `V_sol, epseta, Rtarget/Rtargetref, drain, Delta_Keta, Delta_mu`. THEN positive control: `expectTrue[D[a_V (f_U + Lam)/(a_U a_V - chi_UV^2), Lam] =!= 0]` and `expectZero[D[a_V (f_U + Lam)/(a_U a_V - chi_UV^2), Lam] - a_V/(a_U a_V - chi_UV^2)]`.
- **M8 — Session-I readback.** With `U_obs=0.14313458, V_obs=-0.03619791, epsref=0.3`: `expectApprox[epsref Exp[V_obs], 0.28933482, 5*10^-9]`, `expectApprox[((1 - epsref Exp[V_obs])/(1 - epsref)) Exp[-U_obs], 0.87984149, 5*10^-9]`, `expectApprox[-U_obs - epsref V_obs/(1 - epsref), -0.12762119, 5*10^-9]`. (Assert the recorded session numbers directly — do NOT reconstruct parameters and round-trip them, which is the .py defect F3.)

Declare symbols with assumptions matching the physics: `a_U, a_V, T2ref, epsref, Lambda0` positive; `chi_UV, f_U, U, V, u1, v1, eps, Lam, varrho` real. `0 < epsref < 1` where needed for the `R_target` ratio (so `1 - epsref > 0`).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 245` (or `math -script` on the new `.wl`) and confirms M1–M8 checks appear AND the script exits 0.

## Applied: F5

- files_changed:
  - `mathematica/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.wl`
- summary: Created an independent Mathematica audit script covering manifest checks M1 through M8.
- deviation: none
