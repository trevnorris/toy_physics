---
unit_id: 251
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-03T13:13:55-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 251

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the targets named. Do NOT change any symbolic form or numeric constant that matches the paper card/notes — only change the *checks*.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl`

**Issue:** Stage 251 has no Mathematica engine. Every claim is reachable with native Mathematica primitives, so the dual-engine rule requires an independent `.wl`. It must NOT transliterate the SymPy script: use a *different decomposition* for each shared claim (e.g., `Series` where SymPy used `Limit`; `Reduce`/`Resolve` for sign/uniqueness data that SymPy only checked as an algebraic form or at sample points).

**Required change:**
Create the file at the target path (note: `.wl` files live in `mathematica/`). Use a small helper-assertion style consistent with the repo's other `.wl` audits (`expectZero[expr_] := ...`, `expectTrue[...]`, `expectApprox[a_,b_,tol_]`), strip `ConditionalExpression` from any `Solve`/`Reduce` output, and `Exit[1]` on any failed check. Independently derive and assert the claim manifest M1-M7 below. End by printing a PASS banner.

Implementation notes to keep it a genuine second route (not a port):
- M1: build `N0[ω] = (A0[ω] gW0[ω])^2 / (A0[ω] W0[ω] - R0^2)^2` and extract the cubic-order content via `Series[N0[ω],{ω,0,4}]` then `Coefficient[..., ω, 2]` — a Series route, not the SymPy `Limit[N0/ω^2,ω->0]` route.
- M5: get monotonicity/uniqueness via `Reduce[ForAll[s, s>0, F'[s]>0], {Γ3,Γ5,μ,...}>0]`-style positivity (or `Resolve`), plus `expectTrue[F[0] < 0]` under positive `κV` and `expectTrue[Limit[F[s],s->Infinity]==Infinity]`. SymPy only checked the `F'` *form* and sampled two points; do the general positivity here.
- M6: derive the slowdown by `Series[F[s0+ε] /. κV->μ s0^2 /. {Γ3->g Γ3, Γ5->g Γ5}, {g,0,1}]`, solve the `O(g)` coefficient for `ε`, and `expectZero` against `-(Γ3 s0^2 + Γ5 s0^4)/(2 μ)`.
- M7: `expectApprox` the benchmark numbers against the paper card values (289.61004918, 961.09429528, weight 0.3013336471) using the same input data `t_cross=1.82169718`, `γ_crit=6.94311167`.

**Claim manifest** (the new `.wl` must independently verify each):
- **M1** (cubic): with `A0(ω)=Ω_{U0}^2-ω^2`, `W0(ω)=Ω_{W0}^2-ω^2`, `g_{W,0}(ω)=η0 ω`, `Δ0=Ω_{U0}^2 Ω_{W0}^2 - R0^2`, and `N0(ω)=(A0 g_{W,0})^2/(A0 W0 - R0^2)^2`, the `ω^2` coefficient of `N0` is `η0^2 Ω_{U0}^4/Δ0^2`, and `Γ_{3,0}=γ1·(that)`, `Γ3=Π_{V0}^2 Γ_{3,0}`.
- **M2** (quintic structure): `Γ_{5,-}=a^5/(27 c_s^5)·P_{0,-}`, `P_{0,-}=β0 s_-/λ_-`, `Γ5=Π_{V-}^2 Γ_{5,-}`. Verify the projection homogeneity `Γ5(k·Π_{V-}) = k^2 Γ5(Π_{V-})` and that the `s^5` coefficient of `K_exp = Γ3 s^3 + Γ5 s^5` is exactly `Γ5` — do NOT verify by re-expanding `Γ5` against itself.
- **M3** (kernel): `K_exp(s)=Γ3 s^3 + Γ5 s^5`; `Σ_exp(ω) = -I Γ3 ω^3 - I Γ5 ω^5` (the `s^3`,`s^5` coefficients of `K_exp` are `Γ3`,`Γ5`; even powers absent).
- **M4** (Schott): with `F_odd = Γ3 V'''(t) - Γ5 V^{(5)}(t)` and `S_odd = Γ3 V' V'' - Γ5 (V' V^{(4)} - V'' V''')`, the identity `V'(t)·F_odd - D[S_odd,t] = -Γ3 V''(t)^2 - Γ5 V'''(t)^2` holds (so `P_exp = Γ3 V''^2 + Γ5 V'''^2 ≥ 0`).
- **M5** (characteristic polynomial): `F(s)=Γ5 s^5 + Γ3 s^3 + μ s^2 - κV`; for `Γ3,Γ5,μ,κV>0`: `F(0)=-κV<0`, `Limit[F,s->∞]=+∞`, and `F'(s)=5Γ5 s^4+3Γ3 s^2+2μ s>0` for all `s>0` (use `Reduce`/`Resolve`, establishing exactly one positive root).
- **M6** (slowdown): with `κV=μ s0^2`, the unique positive root expands as `s_+ = s0 - (Γ3 s0^2 + Γ5 s0^4)/(2μ) + O(Γ^2)` (derive via Series in kernel strength `g`, then `expectZero` the `O(g)` coefficient against the boxed expression).
- **M7** (safe surface + benchmark): solving `F(s_c)=0` for `Γ3` with `κV=μ s0^2` gives `Γ̂3 + s_c^2 Γ̂5 = (s0^2 - s_c^2)/s_c^3` (`Γ̂_n=Γ_n/μ`); and with `t_cross=1.82169718`, `s_c=1/t_cross`, `s0=γ_crit=6.94311167`: `s_c^2≈0.3013336471`, rhs `≈289.61004918`, `Γ̂_{5,safe}=(s0^2-s_c^2)/s_c^5≈961.09429528`. Also confirm via `NSolve` that the normalized polynomials `s^2 + Γ̂3_safe s^3 - s0^2` and `s^2 + Γ̂5_safe s^5 - s0^2` each have exactly one positive real root, equal to `s_c`.

**Verification command:**
The verifier runs `redteam exec-mathematica 251`, confirms the `.wl` exists at the target path, the M1-M7 checks (`expectZero`/`expectTrue`/`expectApprox`) appear, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl`
- summary: Created the independent Mathematica audit covering M1-M7 with Series, Resolve, bookkeeping-strength, benchmark, and root checks.
- deviation: none

---

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py:81`

**Issue:** The assertion at line 81 subtracts the assembled `Gamma5` (built at lines 69-73 from exactly `PiVm**2 * a**5/(27*cs**5) * beta0*sminus/lamminus`) from a literal copy of those same factors. It is identically zero by construction and exercises no physics. The cubic check at line 59 is genuine (a `limit` extraction); the quintic check is not.

The paper card (`eq:app-part08-stage251-gamma5`) and notes (Section 9) carry `Γ_{5,-}=a^5/(27 c_s^5)P_{0,-}`, `P_{0,-}=β0 s_-/λ_-` as imported placeholders — so do NOT change the symbolic value. Replace the *tautological assertion* with checks that can actually fail if the prefactor power or projection power is wrong.

**Required change:**
Replace the single assertion at line 81 with non-tautological structural checks. Suggested (adjust names to existing symbols):
1. Build the kernel coefficient independently and extract the quintic piece, then compare to `Gamma5`:
   ```python
   s_kernel = sp.symbols("s_kernel", positive=True)
   K_exp = sp.expand(Gamma3 * s_kernel**3 + Gamma5 * s_kernel**5)
   coeff_s5 = sp.Poly(K_exp, s_kernel).coeff_monomial(s_kernel**5)
   assert sp.simplify(coeff_s5 - Gamma5) == 0          # quintic sits at s^5, no even/cubic leakage
   assert sp.Poly(K_exp, s_kernel).coeff_monomial(s_kernel**3) == Gamma3  # and cubic at s^3
   ```
2. Projection homogeneity (degree-2 in the projection amplitude), which fails if the `Π_{V-}^2` power is mistyped:
   ```python
   k = sp.symbols("k", positive=True)
   assert sp.simplify(Gamma5.subs(PiVm, k*PiVm) - k**2 * Gamma5) == 0
   ```
3. First-order dependence of `P0_minus` on each microscopic factor (fails if a power slips):
   ```python
   assert sp.simplify(P0_minus.subs(beta0, k*beta0) - k * P0_minus) == 0
   assert sp.simplify(P0_minus.subs(sminus, k*sminus) - k * P0_minus) == 0
   assert sp.simplify((P0_minus.subs(lamminus, k*lamminus)) - P0_minus / k) == 0
   ```
Keep `Gamma5` itself (and its symbolic value) unchanged; only the assertions change. Do not reuse `Gamma5`'s own assembled literal on both sides of any single check.

**Verification:**
Line 81's single self-equality is replaced by the checks above; `coeff_s5` is extracted from `K_exp` rather than copied from `Gamma5`'s factors; script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py`
- summary: Replaced the tautological quintic equality with kernel coefficient extraction, projection homogeneity, and microscopic-factor scaling checks.
- deviation: none

---

## F3 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py:126`

**Issue:** The small-kernel slowdown (boxed paper deliverable, card eq `app-part08-stage251-small-kernel`) is only hand-written and printed (`root_shift`, lines 126/135) — no assertion, and not derived from `F(s)`. It restates the paper's answer rather than verifying it.

**Required change:**
Derive the slowdown coefficient from the characteristic polynomial and assert it equals the printed `root_shift`. Add near line 126 (after `F` is defined at line 119; reuse `s0`, `mu_eta`, `G3`, `G5`):
```python
g = sp.symbols("g", positive=True)            # bookkeeping kernel strength
eps = sp.symbols("eps")
F_kappa = F.subs({kappaV: mu_eta * s0**2})    # impose kappa_V = mu_eta s0^2
F_weak = F_kappa.subs({G3: g*G3, G5: g*G5, s: s0 + eps})
# leading balance: O(g^0) gives eps=0; solve O(g^1) for the first-order shift
ser = sp.series(F_weak, g, 0, 2).removeO()
c0 = ser.coeff(g, 0)                          # = mu_eta((s0+eps)^2 - s0^2)
c1 = ser.coeff(g, 1)                          # = G3 (s0+eps)^3 + G5 (s0+eps)^5 (at eps->0 to leading order)
# first-order root: eps = eps1 g + O(g^2); substitute and collect O(g)
eps1 = sp.symbols("eps1")
balance = sp.series((c0 + g*c1).subs(eps, g*eps1), g, 0, 2).removeO().coeff(g, 1)
eps1_sol = sp.solve(balance, eps1)[0]
print("delta s (derived)     =", eps1_sol)
assert sp.simplify(eps1_sol - root_shift) == 0
```
If the exact series choreography above is awkward in-tree, an equivalent and simpler route is acceptable: substitute `s = s0 + g*eps1` into `F.subs(kappaV, mu_eta*s0**2).subs({G3: g*G3, G5: g*G5})`, expand, take the `O(g^1)` coefficient, solve for `eps1`, and assert it equals `root_shift`. The key requirement: the asserted shift must be *derived from F*, not the hand-written `root_shift`, and the two must be shown equal. Keep the existing `print` of `root_shift`.

**Self-test (already done by auditor):** `F(s0)|_{κ=μ s0^2} = Γ3 s0^3 + Γ5 s0^5`, `F'(s0)|_{g=0} = 2 μ s0`, so `eps1 = -F(s0)/F'(s0) = -(Γ3 s0^2 + Γ5 s0^4)/(2 μ)` — exactly `root_shift`. The derivative is in variables `F` genuinely depends on (`s`, and `g` via the substituted coefficients); no identically-zero/independence trap.

**Verification:**
A new assertion `assert sp.simplify(eps1_sol - root_shift) == 0` (or the equivalent) appears near line 126; the derived coefficient prints and matches; script exits 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py`
- summary: Derived the weak-export root shift from the characteristic polynomial with a bookkeeping strength and asserted it matches the printed slowdown.
- deviation: none
