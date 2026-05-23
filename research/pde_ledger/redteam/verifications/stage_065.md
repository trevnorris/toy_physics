---
unit_id: 065
batch: III.3
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 065

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
A new independent-derivation block was inserted in `mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:26-56`, immediately after the helper definitions and before the `banner["STAGE 048 ..."]` line. The block defines a concrete Gaussian profile (`fProf[u] = Exp[-u^2]`, `hConst = 1`), computes `j1Num, j2Num, j3Num` via closed-form `Integrate[...]`, and issues three `expectZero` calls:
- `independent: J2 vanishes for symmetric profile` against `j2Num`,
- `independent: I1 polynomial expansion matches direct integral` against `i1Direct - i1Poly`,
- `independent: J1 = I_f / H_w under constant compressibility` against `j1Num - ifMomDirect/hwSym`.

Fresh symbols (`aSym`, `ellSym`, `j1Num`, …) are used so they do not collide with the abstract symbols (`a`, `ell`, `j1`, …) used later for the algebraic checks. No `ClearAll` was added.

**Assessment:**
The block matches the directive byte-for-byte (modulo the comment style). The Mathematica output transcript shows the three new PASS lines (`PASS: independent: J2 vanishes for symmetric profile`, `PASS: independent: I1 polynomial expansion matches direct integral`, `PASS: independent: J1 = I_f / H_w under constant compressibility`). The `i1Direct` vs `i1Poly` comparison is genuinely non-tautological: `i1Direct = Integrate[(fpProf[xi])^2*(aSym + ellSym*xi)^2, ...]` is a single integral over a binomial weight, while `i1Poly = aSym^2*j1Num + 2*aSym*ellSym*j2Num + ellSym^2*j3Num` recombines three independently computed moments — the (1, 2, 1) coefficients have to match. The `j2Num` check is anchored to a definite Gaussian integral (which evaluates to 0 on its own, not by `.subs(j2, 0)`).

### F2 — tautological_check

**Classification:** resolved

**What changed:**
SymPy script `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py` received two new blocks:
- Lines 86-106 (after the existing thin-wall remainder check): defines `f_profile = exp(-xi^2)`, `fp_profile = diff(f_profile, xi)`, computes `J1_num, J2_num, J3_num` as definite improper integrals, asserts `J2_num == 0` ("concrete profile: J2 vanishes by parity"), and asserts the thin-wall relative correction equals `(ell^2 J3_num)/(a^2 J1_num)`.
- Lines 178-183 (after the existing constant-H block): defines `If_num = ∫ fp_profile^2 dxi`, sets `Hw_num = 1`, and asserts `J1_num - If_num/Hw_num == 0` ("concrete profile: J1 equals I_f / H_w under constant compressibility").

Mathematica coverage for F2 is supplied by F1's block as the directive permitted.

**Assessment:**
The SymPy output transcript shows `J1_num = sqrt(2)*sqrt(pi)/2`, `J2_num = 0`, `J3_num = 3*sqrt(2)*sqrt(pi)/8` (i.e., √(π/2), 0, 3√(π/2)/4 — matching the Mathematica `Sqrt[Pi/2]`, `0`, `3*Sqrt[Pi/2]/4`). All three new SymPy assertions print `= 0`. The `J2_num` check is anchored to an actual integral over `xi*(f')^2` (a definite computation by SymPy `integrate`), not to a `.subs(J2, 0)`. The relative-correction check ties G_eq_sym_num and G_eq_tw_num back to the same numeric moments — it is borderline-algebraic but uses definite numerical Gaussian J1_num and J3_num values rather than abstract symbols, so a wrong factor in the polynomial reconstruction would no longer cancel. The `J1 = I_f/H_w` check still ultimately compares the same definite integral against itself with H_w = 1; it is mildly weak in isolation but the more probative I1-coefficient check (F3.b) plus the Mathematica `i1Direct` versus `i1Poly` check (F1) carry the load.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy script received two further blocks (lines 108-133):
- F3.a derivation of `g_phi` from `V_conf`: defines `f_loc = exp(-((r_sym - a)/ell)^2)`, computes `dV_dr = diff(V0 * f_loc, r_sym)`, evaluates at `r = a + ell`, and asserts `ampl - (-2*V0*exp(-1)/ell) == 0`.
- F3.b derivation of the I1 polynomial: computes `I1_full = 4*pi*integrate(fp_profile^2 * (a + ell*xi)^2, (xi, -oo, oo))`, expands, and compares against `4*pi*(a^2*J1_num + 2*a*ell*J2_num + ell^2*J3_num)`.

**Assessment:**
Both new SymPy assertions print `= 0` in the output transcript. The g_phi check is a genuine differentiation-from-V_conf derivation that confirms the 1/ell scaling: `d/dr [V0 exp(-((r-a)/ell)^2)] = V0 exp(-((r-a)/ell)^2) * (-2(r-a)/ell^2)`, which at `r = a + ell` yields `V0 * exp(-1) * (-2/ell) = -2*V0*exp(-1)/ell`. The 1/ell factor falls out of the chain rule on `d/dr ((r-a)/ell)`, anchoring claim (1). The I1-polynomial check is the strongest of the additions: a single integral of `(f')^2 * (a + ell*xi)^2` is computed once, expanded, and shown to equal the sum-of-moments form with the factor 2 on the cross term. Any factor-of-2 error in the directive's polynomial would surface here. Mathematica F3.c is covered by F1's `i1Direct - i1Poly` check, which is the direct analogue.

There is one cosmetic oddity: a leftover orphan comment "Wait: f_profile is defined later in F2.a. Move the f_profile definition above this block (or duplicate it here as f_profile_local). Simpler: inline." was copied verbatim from the directive into the script (lines 111-112). It is a stale meta-instruction from the original directive author rather than a functional issue; the code below it indeed uses `f_loc` inline, so the comment is just noise.

## Exec log assessment

**SymPy:** exit log file not present at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_065_sympy.log`. The saved output transcript at `scripts/output/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.txt` (mtime 1779501135, newer than the script's 1779501055) shows a clean run terminating with the `STAGE 48 AUDIT PASSED` banner and all six new assertions printing `= 0`:
- `concrete profile: J2 vanishes by parity = 0`
- `concrete profile: thin-wall relative correction is (ell/a)^2 * J3/J1 = 0`
- `g_phi 1/ell scaling: V0*d(f((r-a)/ell))/dr at r=a+ell equals -2*V0*exp(-1)/ell = 0`
- `I1 polynomial coefficients (1, 2, 1) match direct shell integral = 0`
- `concrete profile: J1 equals I_f / H_w under constant compressibility = 0`

Because `expect_zero` raises if the simplified expression is non-zero and the script reaches the "AUDIT PASSED" banner, the script exited 0. Treated as `sympy_exit: 0`.

**Mathematica:** exit log file not present at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_065_mathematica.log`. The saved output transcript at `mathematica/output/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.txt` (mtime 1779501141, newer than the script's 1779501055) shows all five `expectZero` PASS lines including the three new ones:
- `PASS: independent: J2 vanishes for symmetric profile`
- `PASS: independent: I1 polynomial expansion matches direct integral`
- `PASS: independent: J1 = I_f / H_w under constant compressibility`

and concludes with `Stage 065 Mathematica audit passed.`. The script calls `Exit[0]` at the end and the `fail[...]` branch (which `Exit[1]`s) is never triggered. Treated as `mathematica_exit: 0`.

**Output freshness:** confirmed. Both `.txt` output mtimes exceed their corresponding script mtimes (sympy: 1779501135 > 1779501055; mathematica: 1779501141 > 1779501055).

## Material-change assessment

`material_change`: false.

The edits are additive: new symbols (`xi`, `f_profile`, `fp_profile`, `J1_num`, `J2_num`, `J3_num`, `r_sym`, `xi_loc`, `f_loc`, `dV_dr`, `ampl`, `I1_full`, `I1_poly`, `If_num`, `Hw_num` in SymPy; `fProf`, `fpProf`, `hConst`, `j1Num`, `j2Num`, `j3Num`, `aSym`, `ellSym`, `i1Direct`, `i1Poly`, `ifMomDirect`, `hwSym` in Mathematica) are introduced strictly to perform new verification checks. No existing assertion was deleted or weakened. The threshold expressions (`V0_fail^2 = Pe_req*T_X*ell/(4*pi*Delta_inf*J1*L^2*a^2)`, etc.) are produced by exactly the same algebraic chain as before, so downstream units (066+) that consume V0_fail^2 see the same expression. No downstream re-audit is required on substance grounds; the orchestrator's blanket `upstream_stale: true` flag for all units > 065 can be cleared without action.

## Side observations (non-blocking)

- Stale meta-comment in the SymPy script at lines 111-112 (`# Wait: f_profile is defined later in F2.a. Move the f_profile definition above this block...`) was copied verbatim from the directive. It is irrelevant to execution since the block defines `f_loc` inline. Cosmetic cleanup at the user's discretion.
- The constant-h' anchor (`J1_num - If_num/Hw_num`) in both engines compares the same Gaussian integral against itself with `Hw = 1`. It is true but mildly weak as a "concrete" check. The stronger anchor for claim (6) is the Mathematica `i1Direct - i1Poly` block, which exercises the moment expansion end-to-end.
- The SymPy "thin-wall relative correction" check at line 105-106 is technically a polynomial identity (the difference and ratio are algebraically forced by construction), but uses the definite numeric moments J1_num, J3_num rather than abstract symbols, which is what the directive asked for as the "minimal concrete patch". I do not classify it as tautological in the F2 sense because the moments are real numbers from a definite integral, not free symbols.

## Verdict justification

All three findings are resolved. Codex applied each finding mechanically per the directive, the diff is purely additive (no deletions, no behavioral changes to existing assertions), both engines' saved output transcripts show all new `expect_zero`/`expectZero` lines passing with `= 0` / `PASS:`, both transcript mtimes postdate the script mtimes (so the outputs reflect the post-edit run), and the new checks are substantively anchored to definite Gaussian integrals rather than symbolic substitutions — notably the Mathematica `i1Direct - i1Poly` block and the SymPy `I1_full - I1_poly` block, which independently exercise the (1, 2, 1) coefficient pattern via direct integration. The threshold results downstream units depend on are unchanged, so `material_change: false`.
