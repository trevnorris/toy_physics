---
unit_id: 003
batch: I.1
created_at: 2026-05-20T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T07:04:35Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 003

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:130-151`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:189-210`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:78-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:115-122`

**Issue:**

The effective wall kernel `D_eff(omega) = K - omega^2 M - C (Omega_m^2 - omega^2 I)^{-1} C^T` is asserted in both engines but never derived by actually eliminating the matter amplitudes (xa, xb) from the coupled Euler-Lagrange equations. The current `D_eff - manual = 0` check only verifies that SymPy/Mathematica matrix algebra agrees with a hand-expanded 2x2 form of the same formula. A sign error in the claimed elimination result would pass. The same hazard applies to the one-mode dispersion `(K - M w2)(varpi2 - w2) - g^2 = 0` in Section II — also asserted, not derived.

**Required change:**

In `moving_throat_pde_stage003_bdg_sympy_audit.py`, inside `axisymmetric_matrix_kernel_audit()`, after the existing I.1 block (which already declares `EL_qa, EL_qL, EL_xa, EL_xb`) and *before* the `subbanner("I.2 ...")` at line 130, add an explicit elimination subsection. Step-by-step:

1. Introduce frequency-space amplitude symbols `Qa, QL, Xa, Xb = sp.symbols("Qa QL Xa Xb")`.
2. Substitute the ansatz `qa -> Qa*sp.exp(-sp.I*omega*t)`, `qL -> QL*sp.exp(-sp.I*omega*t)`, `xa -> Xa*sp.exp(-sp.I*omega*t)`, `xb -> Xb*sp.exp(-sp.I*omega*t)` into each `EL_*.lhs`. Multiply through by `sp.exp(+sp.I*omega*t)` and `sp.simplify` to obtain linear equations in the amplitudes.
3. Solve the two `EL_xa`, `EL_xb` equations for `Xa, Xb` in terms of `Qa, QL` using `sp.solve`.
4. Substitute back into the `EL_qa`, `EL_qL` equations to get two equations linear in `Qa, QL` only.
5. Read off the resulting 2x2 matrix `D_derived` from the coefficients of `Qa, QL`.
6. Assert `expect_zero("derived D0_eff vs Deff", sp.simplify(D_derived - Deff))` (Deff is computed at line 135).

Insert this elimination block between current line 129 (blank after I.1) and line 130 (`subbanner("I.2 ...")`).

Inside `one_mode_pole_shift_audit()` at line 189, *before* the existing `dispersion = ...` line at line 197, add an explicit Lagrangian derivation:

1. Define `q = sp.Function("q")(t)`, `x = sp.Function("x")(t)`.
2. Build `L_one = M*sp.diff(q,t)**2/2 - K*q**2/2 + sp.diff(x,t)**2/2 - varpi2*x**2/2 + g*q*x` (note: here `varpi2` plays the role of `varpi^2` already in scope; if the parametrization uses `varpi` not `varpi^2`, follow the existing convention).
3. Derive `EL_q = sp.euler_equations(L_one, q, [t])[0]` and `EL_x = sp.euler_equations(L_one, x, [t])[0]`.
4. Substitute frequency-space ansatz `q -> Q*exp(-I*omega*t)`, `x -> X*exp(-I*omega*t)`, divide by `exp(-I*omega*t)`, simplify.
5. Solve EL_x for X in terms of Q.
6. Substitute back into EL_q; the prefactor of Q is the derived dispersion (up to a `(varpi2 - omega^2)` denominator factor).
7. Assert `expect_zero("derived dispersion vs (K - M w2)(varpi2 - w2) - g^2", ...)` after clearing the denominator (replace `omega**2 -> w2` to compare to the existing `dispersion` variable).

Now mirror both changes in `moving_throat_pde_stage003_bdg_mathematica_audit.wl`:

For Section I, insert between line 77 (after kinetic coefficient checks; note that lines 64-72 will themselves be replaced per F2 to provide EL equations) and line 78 (`mMat = ...`):

```
(* derive D_eff by eliminating Xa, Xb from the EL equations *)
elQa = D[D[lRed, D[qa, t]], t] - D[lRed, qa];
elQL = D[D[lRed, D[qL, t]], t] - D[lRed, qL];
elXa = D[D[lRed, D[xa, t]], t] - D[lRed, xa];
elXb = D[D[lRed, D[xb, t]], t] - D[lRed, xb];
ansatz = {qaFun[t] -> Qa Exp[-I omega t], qLFun[t] -> QL Exp[-I omega t],
          xaFun[t] -> Xa Exp[-I omega t], xbFun[t] -> Xb Exp[-I omega t]};
elQaF = FullSimplify[elQa /. ansatz, Assumptions -> $Assumptions]/Exp[-I omega t] // FullSimplify;
elQLF = FullSimplify[elQL /. ansatz, Assumptions -> $Assumptions]/Exp[-I omega t] // FullSimplify;
elXaF = FullSimplify[elXa /. ansatz, Assumptions -> $Assumptions]/Exp[-I omega t] // FullSimplify;
elXbF = FullSimplify[elXb /. ansatz, Assumptions -> $Assumptions]/Exp[-I omega t] // FullSimplify;
xsol = Solve[{elXaF == 0, elXbF == 0}, {Xa, Xb}][[1]];
elQaRed = FullSimplify[elQaF /. xsol, Assumptions -> $Assumptions];
elQLRed = FullSimplify[elQLF /. xsol, Assumptions -> $Assumptions];
dDerived = {{Coefficient[elQaRed, Qa], Coefficient[elQaRed, QL]},
            {Coefficient[elQLRed, Qa], Coefficient[elQLRed, QL]}};
expectMatrixZero["derived D0_eff vs Deff", dDerived - dEff];
```

For Section II, insert between line 114 and line 115:

```
(* derive dispersion from a single-mode Lagrangian *)
Clear[qFun, xFun, Q, X];
lOne = m/2 D[qFun[t], t]^2 - k/2 qFun[t]^2 + 1/2 D[xFun[t], t]^2 - varpi2/2 xFun[t]^2 + g qFun[t] xFun[t];
elQ = D[D[lOne, D[qFun[t], t]], t] - D[lOne, qFun[t]];
elX = D[D[lOne, D[xFun[t], t]], t] - D[lOne, xFun[t]];
ansatz1 = {qFun[t] -> Q Exp[-I omega t], xFun[t] -> X Exp[-I omega t]};
elQF = FullSimplify[elQ /. ansatz1, Assumptions -> $Assumptions]/Exp[-I omega t] // FullSimplify;
elXF = FullSimplify[elX /. ansatz1, Assumptions -> $Assumptions]/Exp[-I omega t] // FullSimplify;
xSol1 = Solve[elXF == 0, X][[1]];
elQRed = FullSimplify[elQF /. xSol1, Assumptions -> $Assumptions];
dispersionDerived = FullSimplify[Together[elQRed (varpi2 - omega^2) / Q], Assumptions -> $Assumptions];
expectZero["derived dispersion vs (k - m w2)(varpi2 - w2) - g^2",
           (dispersionDerived /. omega^2 -> w2) - dispersion];
```

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 003` and `redteam exec-mathematica 003` and confirm:
- SymPy output contains a new line `derived D0_eff vs Deff = 0` and a new line `derived dispersion vs ... = 0`, both before any `EXIT_CODE: 0`.
- Mathematica output contains new lines `derived D0_eff vs Deff = ...` and `derived dispersion vs (k - m w2)(varpi2 - w2) - g^2 = 0` with `PASS:` prefixes.
- Both scripts still exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- summary: Added explicit frequency-space elimination checks for the two-mode kernel and one-mode dispersion in both audit engines.
- deviation: Used the SymPy Euler-operator sign and defined the target dispersion before comparing so the requested assertions execute.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:64-72`

**Issue:**

Lines 70-72 extract literal kinetic coefficients of `vqa^2`, `vqL^2`, `vqa vqL` from a Lagrangian that was just written with those exact coefficients (`1/2 maa D[qa, t]^2 + ...`). These three assertions cannot fail; they convey no verification information.

**Required change:**

Replace lines 61-72 (currently `staticL`, `staticTmp`, `staticBack`, the velocity-renaming `lVel`, and the three `expectZero` kinetic-coefficient calls) with explicit Euler-Lagrange equation checks that match the SymPy script's I.1 block (lines 107-128 of the .py).

Specifically, replace lines 61-72 of the .wl with:

```
elQa = D[D[lRed, D[qaFun[t], t]], t] - D[lRed, qaFun[t]];
elQL = D[D[lRed, D[qLFun[t], t]], t] - D[lRed, qLFun[t]];
elXa = D[D[lRed, D[xaFun[t], t]], t] - D[lRed, xaFun[t]];
elXb = D[D[lRed, D[xbFun[t], t]], t] - D[lRed, xbFun[t]];

expectZero["qa equation",
  elQa - (maa qaFun''[t] + maL qLFun''[t] + kaa qaFun[t] + kaL qLFun[t]
          - c1a xaFun[t] - c1b xbFun[t])];
expectZero["qL equation",
  elQL - (maL qaFun''[t] + mLL qLFun''[t] + kaL qaFun[t] + kLL qLFun[t]
          - c2a xaFun[t] - c2b xbFun[t])];
expectZero["xa equation",
  elXa - (xaFun''[t] + wa^2 xaFun[t] - c1a qaFun[t] - c2a qLFun[t])];
expectZero["xb equation",
  elXb - (xbFun''[t] + wb^2 xbFun[t] - c1b qaFun[t] - c2b qLFun[t])];
```

Note the sign convention: Mathematica's `D[D[L,q'],t] - D[L,q]` gives the EL operator with the convention `EL = d/dt(dL/dq') - dL/dq`. The expected expressions above are written so that adding them to `elQa` (etc.) gives zero — i.e., the SymPy convention is `dL/dq - d/dt(dL/dq')` and the .py adds `+ Maa qa'' + ...`. Here, Mathematica's convention is the opposite sign, so we subtract `(Maa qa'' + ...)`. Verify by hand that for the explicit `lRed`, `elQa = maa qa'' + maL qL'' + kaa qa + kaL qL - c1a xa - c1b xb` and `elQa - (that same expression) = 0`.

Then renumber: the subsequent `mMat = ...` block (currently starting at line 74) follows immediately. The F1 block (elimination of Xa, Xb in freq space) is inserted right after these EL checks.

**Verification command:**

After Codex applies, verify that the .wl output transcript for Section I includes four PASS lines: `qa equation`, `qL equation`, `xa equation`, `xb equation` — matching the SymPy output lines 17-20. The three tautological kinetic-coefficient PASS lines should no longer appear.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- summary: Replaced the tautological kinetic coefficient checks with four Euler-Lagrange residual checks.
- deviation: Completed `lRed` inside the replacement block first because the existing Mathematica multiline assignment evaluated only the kinetic terms.

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:322-331`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:182-190`

**Issue:**

The selection-rule claim ("the angular overlap is diagonal in (l,m)") requires more orthonormality checks than the script provides. Missing: `<Y00, Y21c>`, `<Y00, Y22c>`, `<Y21c, Y22c>`, and the norms of Y00, Y21c, Y22c.

**Required change:**

In `.py` at line 322 (after `I00_20 = sp.simplify(...)`) and before line 326, add:

```python
I00_21c = sp.simplify(sp.integrate(sp.integrate(Y00 * Y21c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
I00_22c = sp.simplify(sp.integrate(sp.integrate(Y00 * Y22c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
I21c_22c = sp.simplify(sp.integrate(sp.integrate(Y21c * Y22c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
N00 = sp.simplify(sp.integrate(sp.integrate(Y00 * Y00 * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
N21c = sp.simplify(sp.integrate(sp.integrate(Y21c * Y21c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
N22c = sp.simplify(sp.integrate(sp.integrate(Y22c * Y22c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
```

Then at line 328 (within the `expect_zero` block), add immediately after the existing three cross-integral assertions and before the `Y20 norm - 1` assertion:

```python
expect_zero("Y00-Y21c cross integral", I00_21c)
expect_zero("Y00-Y22c cross integral", I00_22c)
expect_zero("Y21c-Y22c cross integral", I21c_22c)
expect_zero("Y00 norm - 1", N00 - 1)
expect_zero("Y21c norm - 1", N21c - 1)
expect_zero("Y22c norm - 1", N22c - 1)
```

Mirror in `.wl` between lines 185 and 187:

```
i0021c = FullSimplify[Integrate[Integrate[y00 y21c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i0022c = FullSimplify[Integrate[Integrate[y00 y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i21c22c = FullSimplify[Integrate[Integrate[y21c y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm00 = FullSimplify[Integrate[Integrate[y00 y00 dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm21c = FullSimplify[Integrate[Integrate[y21c y21c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm22c = FullSimplify[Integrate[Integrate[y22c y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
```

After line 189 (the existing Y20-Y22c assertion), insert (before the Y20-norm line 190):

```
expectZero["Y00-Y21c cross integral", i0021c];
expectZero["Y00-Y22c cross integral", i0022c];
expectZero["Y21c-Y22c cross integral", i21c22c];
expectZero["Y00 norm - 1", norm00 - 1];
expectZero["Y21c norm - 1", norm21c - 1];
expectZero["Y22c norm - 1", norm22c - 1];
```

(F4 will optionally replace this hand-picked list of six new checks with a single 4x4 overlap-matrix-equals-identity check in the .wl. Apply F3 first as listed here, then F4 will refactor it.)

**Verification command:**

After Codex applies, the SymPy output's Section IV should contain at least 10 lines of the form `... = 0` or `... = 0` (was 4); the Mathematica output likewise (was 4).

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- summary: Added the missing spherical harmonic cross-overlap and norm derivations.
- deviation: Mathematica scalar assertions were consolidated by F4 as directed.

## F4 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:136-169` (Section III)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:171-190` (Section IV)

**Issue:**

Sections III and IV of the Mathematica script are line-by-line transliterations of the SymPy script: same definitions of `D20, D21, D22`, same coefficient extraction at `omega^0, omega^2, omega^4`, same `d2bar, a2, b2` combinations, same isotropy substitution; same four spherical harmonics in the same order, same four integrals in the same order. The Mathematica engine is not providing an independent derivation.

**Required change:**

**For Section IV (apply after F3):** Replace the (now ten) per-pair integral checks in the .wl with a single overlap-matrix check.

After line 181 of the .wl (or wherever F3 inserted the new harmonic definitions), insert:

```
yList = {y00, y20, y21c, y22c};
overlap = FullSimplify[
  Table[Integrate[Integrate[yList[[i]] yList[[j]] dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}],
        {i, 1, 4}, {j, 1, 4}],
  Assumptions -> $Assumptions];
expectMatrixZero["spherical harmonic overlap matrix - identity",
                 overlap - IdentityMatrix[4]];
```

Then remove all individual `expectZero[..., i0020]`, `i2021c`, `i2022c`, `i0021c`, `i0022c`, `i21c22c`, `norm00`, `norm20 - 1`, `norm21c - 1`, `norm22c - 1` calls that F3 introduced — they're now subsumed by the matrix check. (Leave the integral definitions in place if you wish, but the assertions should be one matrix check, not ten scalar checks.)

This is structurally not a transliteration of the SymPy script (which checks individual pairs).

**For Section III:** Restructure so the second-engine path differs from SymPy. The minimal acceptable change:

Replace lines 143-156 (the `d20, d21, d22` definitions, the series extractions, and the `d220, d221, d222`, `d2Bar, a2, b2` definitions) with a derivation that first writes the *3-channel* self-energy as a single block-diagonal `3x3` matrix `dP2 = DiagonalMatrix[{d20, d21, d22}]`, then constructs the trace/anisotropy invariants via explicit projection onto representation-theoretic basis matrices `T0 = (1/Sqrt[5]) DiagonalMatrix[{1, Sqrt[2], Sqrt[2]}]` (trace direction, with l=2 channel weights), `Ta = (1/Sqrt[10]) DiagonalMatrix[{2, -1, -1}]` (uniaxial anisotropy), `Tb = (1/Sqrt[2]) DiagonalMatrix[{0, 1, -1}]` (biaxial anisotropy). Then:

```
dP2 = DiagonalMatrix[{d20, d21, d22}];
dP2s = Map[Normal[Series[#, {omega, 0, 4}]] &, dP2, {2}] // Expand;
d2coeffMat = Map[Coefficient[#, omega, 2] &, dP2s, {2}];

(* projections onto representation-theoretic basis *)
d2Bar = FullSimplify[Tr[d2coeffMat]/5, Assumptions -> $Assumptions];  (* trace weight 1+2+2 = 5 in this basis *)
a2 = FullSimplify[(2 d2coeffMat[[1,1]] - d2coeffMat[[2,2]] - d2coeffMat[[3,3]])/10,
                  Assumptions -> $Assumptions];
b2 = FullSimplify[(d2coeffMat[[2,2]] - d2coeffMat[[3,3]])/2,
                  Assumptions -> $Assumptions];
```

The defining formulas for `d2Bar, a2, b2` are still recognizable, but the derivation path goes through a diagonal-matrix construction rather than three independent scalar definitions, which makes the second engine's algebra structurally distinguishable from SymPy's.

If a more structurally independent rewrite cannot be made mechanically, append `## Blocked: F4` instead of applying — the orchestrator will route to a human reviewer.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 003` and confirm:
- Section IV output contains exactly one new line `spherical harmonic overlap matrix - identity = ...` (with the matrix simplifying to all zeros) and no per-pair integral PASS lines beyond what's needed (per F3, individual integrals defined are fine — the assertions consolidate).
- Section III defines `dP2` as a diagonal matrix and uses `Tr[d2coeffMat]/5` rather than `(d220 + 2 d221 + 2 d222)/5`.
- Script still exits 0; all final symbolic forms for `a2, b2, d2bar` still match SymPy's at PASS.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- summary: Reworked Mathematica Section III through a diagonal matrix trace path and consolidated Section IV into one overlap-matrix identity check.
- deviation: Weighted `d2coeffMat` before `Tr[d2coeffMat]/5` so the grouped-P2 trace invariant remains the same as the SymPy ledger.
