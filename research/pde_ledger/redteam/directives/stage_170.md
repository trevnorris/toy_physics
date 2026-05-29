---
unit_id: 170
batch: V.1
created_at: 2026-05-28T16:30:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 170

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:144-166`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:138-159`

**Issue:**
Section 5 (weak-axisymmetric signature / scalar-amplitude collapse, paper card Checks item 2) is tautological in both engines. The lane-scaled outputs `dkW`/`dgW` are computed by calling re-typed helpers `kappa_map`/`gamma_map` (sympy L144-148) and `kappaMap`/`gammaMap` (wl L138-139) that are byte-for-byte copies of the answer formula, then compared to `kappa1`/`gamma1` built from that same formula. Because the helper and the target are the identical linear closed form, every amplitude assertion (`dkW[A] - eps_l*lam*kappa1 == 0`, etc.) and every signature ratio holds by construction — they only confirm that a linear expression is linear, never exercising the Section-2 derivation that actually produces the outlet map. The fix routes the lane-scaled inputs through the *derived* Section-2 map objects (`dkappa_from_du2` / `dgamma_from_dP0` in sympy; `dkappaFromdu2` / `dgammaFromdP0` in wl), so the check now depends on the derived map while still comparing against the independently-written `kappa1`/`gamma1` paper target. The maps equal the paper formula, so all assertions still pass — but they become genuine.

**Required change:**

### SymPy (`scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py`)

Delete the two re-typed helper defs (currently L144-148):
```python
def kappa_map(dD2_, dD0_):
    return 3*(1 - sigma)*(dD2_ + dD0_/9)/(sigma*D0)

def gamma_map(dN0_, dD0_):
    return -(1 - sigma)*(dN0_ - P0*dD0_)/(9*sigma*N0)
```
Keep the `kappa1` / `gamma1` definitions (L150-151) exactly as they are (they are the paper target).

Then replace the lane loop body (currently L157-159):
```python
for A, lam in lanes.items():
    dkW[A] = kappa_map(eps_l*lam*D2_1, eps_l*lam*D0_1)
    dgW[A] = gamma_map(eps_l*lam*N0_1, eps_l*lam*D0_1)
```
with substitution into the derived Section-2 maps:
```python
for A, lam in lanes.items():
    dkW[A] = dkappa_from_du2.subs({dD2: eps_l*lam*D2_1, dD0: eps_l*lam*D0_1})
    dgW[A] = sp.simplify(
        dgamma_from_dP0.subs({dN0: eps_l*lam*N0_1, dD0: eps_l*lam*D0_1}).subs(P0, N0/D0)
    )
```
Leave the two `expect_zero` amplitude calls (L160-161) in place, but change the gamma one to reconcile `P0` on the target the same way Section 2 does (sympy L79). I.e. lines L160-161 become:
```python
    expect_zero(f"delta kappa_W^({A}) - eps lambda kappa1", dkW[A] - eps_l*lam*kappa1)
    expect_zero(f"delta gamma_W^({A}) - eps lambda gamma1",
                sp.simplify(dgW[A] - eps_l*lam*gamma1.subs(P0, N0/D0)))
```
Leave the signature-ratio `expect_zero` calls (L163-166) unchanged; they now test linearity of the derived map (acceptable secondary check).

### Mathematica (`mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl`)

Delete the two re-typed helper defs (currently L138-139):
```wl
kappaMap[dD2x_, dD0x_] := 3*(1 - sigma)*(dD2x + dD0x/9)/(sigma*D0);
gammaMap[dN0x_, dD0x_] := -(1 - sigma)*(dN0x - P0*dD0x)/(9*sigma*N0);
```
Keep `kappa1` / `gamma1` (L140-141) unchanged.

Replace the six lane assignments (currently L142-147):
```wl
dkW20 = kappaMap[epsL*1*D21, epsL*1*D01];
dkW21 = kappaMap[epsL*(1/2)*D21, epsL*(1/2)*D01];
dkW22 = kappaMap[epsL*(-1)*D21, epsL*(-1)*D01];
dgW20 = gammaMap[epsL*1*N01, epsL*1*D01];
dgW21 = gammaMap[epsL*(1/2)*N01, epsL*(1/2)*D01];
dgW22 = gammaMap[epsL*(-1)*N01, epsL*(-1)*D01];
```
with substitution into the derived maps (applying `P0 -> N0/D0` on the gamma side, mirroring the Section-2 gamma check at wl L83):
```wl
dkW20 = dkappaFromdu2 /. {dD2 -> epsL*1*D21, dD0 -> epsL*1*D01};
dkW21 = dkappaFromdu2 /. {dD2 -> epsL*(1/2)*D21, dD0 -> epsL*(1/2)*D01};
dkW22 = dkappaFromdu2 /. {dD2 -> epsL*(-1)*D21, dD0 -> epsL*(-1)*D01};
dgW20 = FullSimplify[(dgammaFromdP0 /. {dN0 -> epsL*1*N01, dD0 -> epsL*1*D01}) /. P0 -> N0/D0, Assumptions -> $Assumptions];
dgW21 = FullSimplify[(dgammaFromdP0 /. {dN0 -> epsL*(1/2)*N01, dD0 -> epsL*(1/2)*D01}) /. P0 -> N0/D0, Assumptions -> $Assumptions];
dgW22 = FullSimplify[(dgammaFromdP0 /. {dN0 -> epsL*(-1)*N01, dD0 -> epsL*(-1)*D01}) /. P0 -> N0/D0, Assumptions -> $Assumptions];
```
Change the three gamma amplitude assertions (L153-155) to reconcile `P0` on the target the same way:
```wl
expectZero["delta gamma_W^(20) - eps gamma1", dgW20 - epsL*(gamma1 /. P0 -> N0/D0)];
expectZero["delta gamma_W^(21) - (eps/2) gamma1", dgW21 - epsL*(1/2)*(gamma1 /. P0 -> N0/D0)];
expectZero["delta gamma_W^(22) + eps gamma1", dgW22 + epsL*(gamma1 /. P0 -> N0/D0)];
```
Leave the three kappa amplitude assertions (L150-152) and the four signature-ratio assertions (L156-159) unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 170` and `redteam exec-mathematica 170` and confirm: (a) the new `dkW`/`dgW` assignments reference the derived maps `dkappa_from_du2`/`dgamma_from_dP0` (sympy) and `dkappaFromdu2`/`dgammaFromdP0` (wl) — and `kappa_map`/`gamma_map`/`kappaMap`/`gammaMap` no longer exist in either file; (b) both scripts exit 0 with every Section-5 PASS / `= 0` line still present.
