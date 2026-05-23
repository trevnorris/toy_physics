---
unit_id: 065
batch: III.3
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-23T01:50:04Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 065

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:26-52`

**Issue:** The `.wl` script repeats the `.py` script's polynomial algebra symbol-for-symbol. Both engines define gPhi = v0/ell, i1 = 4*Pi*ell*(a^2 j1 + 2*a*ell*j2 + ell^2*j3), gEq, gEqSym, gEqTw with the same construction, then assert the same "thin-wall remainder" residual `gEqSym - gEqTw - 4*Pi*v0^2*ell*j3/kx`. The Mathematica script provides no independent route to the I1 polynomial coefficients or to J2 = 0; it merely transcribes the SymPy chain. The second-engine policy requires an independent derivation.

**Required change:**
Insert, immediately after line 25 of the `.wl` script (i.e. before the `banner["STAGE 048 ..."]` line), the following independent-derivation block. Use a concrete Gaussian profile so the integrals are closed-form in Mathematica.

```mathematica
banner["INDEPENDENT SHELL-INTEGRAL DERIVATION (concrete Gaussian profile)"];

(* Concrete symmetric profile: f(u) = Exp[-u^2]; constant compressibility h' = 1. *)
fProf[u_] := Exp[-u^2];
fpProf[u_] := D[fProf[u], u];   (* = -2 u Exp[-u^2] *)
hConst = 1;

(* Direct moment integrals of (f')^2 / h'. *)
j1Num = Integrate[(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
j2Num = Integrate[xi*(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
j3Num = Integrate[xi^2*(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];

Print["J1_num = ", fmt[j1Num]];
Print["J2_num = ", fmt[j2Num]];
Print["J3_num = ", fmt[j3Num]];

(* Claim (3): J2 = 0 for a centred symmetric layer (parity). *)
expectZero["independent: J2 vanishes for symmetric profile", j2Num];

(* Claim (2): expanding the full shell weight gives the (1, 2, 1) coefficient pattern. *)
i1Direct = Integrate[(fpProf[xi])^2/hConst*(aSym + ellSym*xi)^2, {xi, -Infinity, Infinity}];
i1Direct = FullSimplify[i1Direct, Assumptions -> aSym > 0 && ellSym > 0];
i1Poly   = aSym^2*j1Num + 2*aSym*ellSym*j2Num + ellSym^2*j3Num;
expectZero["independent: I1 polynomial expansion matches direct integral",
           i1Direct - i1Poly];

(* Claim (6): for constant h' = H_w, J1 = I_f / H_w with I_f = integral of (f')^2. *)
ifMomDirect = Integrate[(fpProf[xi])^2, {xi, -Infinity, Infinity}];
hwSym = 1;  (* hConst above; documented numerically here *)
expectZero["independent: J1 = I_f / H_w under constant compressibility",
           j1Num - ifMomDirect/hwSym];
```

After inserting this block, do not modify the existing `expectZero` calls at lines 49-52, 73-76, 77-80, 89-92, 93-96. They remain valid algebraic identity checks.

Notes on variable hygiene: the new block uses fresh symbol names (`aSym`, `ellSym`, `j1Num`, `j2Num`, `j3Num`, `i1Direct`, `i1Poly`, `ifMomDirect`, `hwSym`) so it does not collide with the existing global symbols (`a`, `ell`, `j1`, `j2`, `j3`, `hw`, `ifMom`). Do NOT call `ClearAll[a, ell, ...]` after this block — those symbols remain abstract for the existing checks.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 065` and confirms three new PASS lines appear in the output transcript:
  - `PASS: independent: J2 vanishes for symmetric profile`
  - `PASS: independent: I1 polynomial expansion matches direct integral`
  - `PASS: independent: J1 = I_f / H_w under constant compressibility`
and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl`
- summary: Inserted the requested independent Gaussian shell-integral derivation block before the existing Mathematica algebra checks.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:81-127`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:49-96`

**Issue:** The "thin-wall remainder" assertion at .py:81-84 and .wl:49-52 is algebraically zero by construction (G_eq_sym is built from the same polynomial whose terms are then subtracted off). The "constant-H fail/suff threshold" assertions at .py:120-127 and .wl:89-96 are equally guaranteed: they take a `.subs(J1, If/Hw)` of an expression and then subtract the result of that substitution. None of these can fail under any choice of physics; they only test SymPy's/Mathematica's substitution arithmetic.

**Required change:**

Add non-tautological replacements that anchor the same claims to concrete integrals. Do NOT delete the existing `expect_zero`/`expectZero` calls; they remain valid identity checks. Add the following after the existing assertions.

### F2.a — SymPy: anchor the thin-wall O(ell^2/a^2) scaling to a concrete profile

Append after line 84 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py` (i.e. right after the existing `expect_zero("thin-wall remainder ...", ...)` block):

```python
# Concrete-profile anchor: f(u) = exp(-u^2), constant h' = 1, gives definite
# numeric J1, J3 so the ratio (G_eq_sym - G_eq_tw)/G_eq_tw = (ell^2 J3)/(a^2 J1)
# becomes a non-trivial scaling identity to verify.
xi = sp.symbols("xi", real=True)
f_profile = sp.exp(-xi**2)
fp_profile = sp.diff(f_profile, xi)
J1_num = sp.integrate(fp_profile**2, (xi, -sp.oo, sp.oo))
J2_num = sp.integrate(xi * fp_profile**2, (xi, -sp.oo, sp.oo))
J3_num = sp.integrate(xi**2 * fp_profile**2, (xi, -sp.oo, sp.oo))
print(f"J1_num = {J1_num}")
print(f"J2_num = {J2_num}")
print(f"J3_num = {J3_num}")

# Independent assertion: J2 vanishes by parity for the symmetric profile.
expect_zero("concrete profile: J2 vanishes by parity", J2_num)

# Independent assertion: the relative correction matches (ell/a)^2 * J3/J1.
G_eq_sym_num = (4*sp.pi*V0**2*(a**2*J1_num + ell**2*J3_num) / (KX*ell))
G_eq_tw_num  = 4*sp.pi*a**2*V0**2*J1_num/(KX*ell)
rel = sp.simplify((G_eq_sym_num - G_eq_tw_num)/G_eq_tw_num - (ell**2*J3_num)/(a**2*J1_num))
expect_zero("concrete profile: thin-wall relative correction is (ell/a)^2 * J3/J1", rel)
```

### F2.b — SymPy: anchor J1 = I_f / H_w to definite integrals

Append after line 127 of the same file (after the existing constant-H block):

```python
# Concrete-profile anchor for the constant-compressibility reduction: with the
# same Gaussian profile, define I_f and verify J1 = I_f / H_w with H_w = h' = 1.
If_num = sp.integrate(fp_profile**2, (xi, -sp.oo, sp.oo))
Hw_num = sp.Integer(1)
expect_zero("concrete profile: J1 equals I_f / H_w under constant compressibility",
            J1_num - If_num/Hw_num)
```

### F2.c — Mathematica: matching anchors

The block proposed in F1 already covers F2.a's J2 = 0 check, the I1 polynomial expansion (which subsumes the thin-wall remainder claim), and F2.b's J1 = I_f / H_w identity. No additional Mathematica changes are required for F2 if F1 is applied. If Codex applies F2 before F1, add the F1 block at that time.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 065` and `redteam exec-mathematica 065`. Expected new PASS lines in the SymPy transcript:
  - `concrete profile: J2 vanishes by parity = 0`
  - `concrete profile: thin-wall relative correction is (ell/a)^2 * J3/J1 = 0`
  - `concrete profile: J1 equals I_f / H_w under constant compressibility = 0`
Mathematica transcript already covered by F1's three new PASS lines. Both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl`
- summary: Added concrete Gaussian SymPy anchors for parity, thin-wall relative correction, and constant-compressibility reduction, with Mathematica coverage supplied by F1's block.
- deviation: none

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:60-66`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:34-35`

**Issue:** Docstring claims (1) `g_phi = V0/ell`, (2) `I1 = 4*pi*ell*(a^2 J1 + 2 a ell J2 + ell^2 J3)`, and (3) `J2 = 0 for symmetric layer` are introduced as symbol definitions / substitutions and never independently exercised. The factor-of-2 cross-term in (2) and the 1/ell scaling in (1) are exactly the kinds of coefficients that need a concrete check.

**Required change:**

### F3.a — SymPy: derive g_phi from V_conf

Insert in the SymPy script, after the existing `print("g_phi =", g_phi)` at line 62, the following derivation:

```python
# Independent derivation of g_phi from V_conf(r) = V0 f((r-a)/ell).
r_sym = sp.symbols("r", positive=True)
V_conf = V0 * f_profile.subs(xi, (r_sym - a)/ell)
# Wait: f_profile is defined later in F2.a. Move the f_profile definition above
# this block (or duplicate it here as f_profile_local). Simpler: inline.
xi_loc = sp.symbols("xi_loc", real=True)
f_loc = sp.exp(-((r_sym - a)/ell)**2)
dV_dr = sp.diff(V0 * f_loc, r_sym)
# At r = a (wall centre) the Gaussian's first derivative vanishes by symmetry;
# instead anchor on the support-loading amplitude g_phi = V0 * f'(0) / ell scaled
# by the canonical normalization f'(0)=1 used by the docstring's f. Verify the
# 1/ell scaling by checking d/dell of |dV_dr| at r = a + ell (one e-fold away).
ampl = sp.simplify(dV_dr.subs(r_sym, a + ell))
# ampl = V0 * (-2/ell) * exp(-1)
expect_zero("g_phi 1/ell scaling: V0*d(f((r-a)/ell))/dr at r=a+ell equals -2*V0*exp(-1)/ell",
            ampl - (-2*V0*sp.exp(-1)/ell))
```

NOTE for Codex: this snippet references `f_profile` which is introduced by F2.a. The two changes interact. The simplest ordering is: apply F2.a first (which defines `f_profile`, `fp_profile`, etc. at module scope), then apply F3.a immediately after F2.a's block, then F3.b after F3.a, then F2.b at the end.

### F3.b — SymPy: derive I1 polynomial coefficients

Insert after the F3.a block:

```python
# Independent derivation of the I1 polynomial.
# Shell integral with constant h' = 1 over the unbounded thin-wall coordinate xi:
#   I1_full = 4*pi * ∫ (f'(xi))^2 * (a + ell*xi)^2 dxi
I1_full = 4*sp.pi*sp.integrate(fp_profile**2 * (a + ell*xi)**2, (xi, -sp.oo, sp.oo))
I1_full = sp.expand(I1_full)
# Expected polynomial form using the numeric moments computed above.
I1_poly = 4*sp.pi*(a**2*J1_num + 2*a*ell*J2_num + ell**2*J3_num)
expect_zero("I1 polynomial coefficients (1, 2, 1) match direct shell integral",
            sp.expand(I1_full - I1_poly))
```

This anchors both the factor-of-2 cross term and the J2 = 0 parity claim in a single non-tautological identity.

### F3.c — Mathematica: matching independent derivations

Already covered by F1's block (the `i1Direct` vs `i1Poly` check and the `j2Num` check). No further Mathematica edits required for F3.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 065`. Expected new PASS lines:
  - `g_phi 1/ell scaling: V0*d(f((r-a)/ell))/dr at r=a+ell equals -2*V0*exp(-1)/ell = 0`
  - `I1 polynomial coefficients (1, 2, 1) match direct shell integral = 0`
SymPy script exits 0. Mathematica side already covered.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl`
- summary: Added concrete SymPy derivations for the g_phi 1/ell scaling and I1 polynomial coefficients, with Mathematica coverage supplied by F1's block.
- deviation: none
