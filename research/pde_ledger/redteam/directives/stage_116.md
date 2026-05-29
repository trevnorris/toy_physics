---
unit_id: 116
batch: IV.3
created_at: 2026-05-28T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-28T20:15:14-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
f2_revised: 2026-05-28  # F2 revised by orchestrator+Codex consult (session 019e717e-ab7b-7db2-a2cc-0b891662bf3c): the audit's proposed "3*kappa0==9*gamma0" replacement was itself tautological (both operands hardcoded literals); agreed resolution is deletion + labeled prints.
---

# Codex directive — unit 116

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:38-46`

**Issue:** The `.wl` is a line-by-line port of the SymPy script: it posits the identical trial function `Sin[k x]` (line 38), runs the identical `First[Solve[...]]` (line 56) corresponding to SymPy's `sp.solve(...)[0]` (sympy line 52-55), and reads the identical polynomial coefficients (`Coefficient[dBare, z, 5]`, `Coefficient[dBare, z, 2]` at lines 73, 78) corresponding to `D_bare.coeff(z, 5)`, `D_bare.coeff(z, 2)`. To serve as an independent second engine, the Mathematica eigenvalue must be *derived*, not posited identically to the SymPy guess.

**Required change:**
Replace the posited-trial-function derivation of the D/N eigenvalue with an independent solve. Concretely, between the `$Assumptions` block (ends line 32) and the `OmegaW` definition (line 48), restructure lines 34-46 so that `kWValue` is obtained from solving the ODE+BCs rather than from a hand-written `Pi/(2 lW)`:

- Keep the comment block (lines 34-37) but update it to say the eigenvalue is now solved, not posited.
- Derive the general solution and characteristic equation independently. Suggested form (Codex may adapt syntax so the script runs and exits 0):
  ```
  (* General solution of q'' + k^2 q = 0 with q(0)=0 is q(x) = A Sin[k x].
     The Neumann BC q'(lW)=0 gives the characteristic equation Cos[k lW]==0.
     Solve for the smallest positive eigenvalue independently. *)
  gensol = DSolve[{q''[xv] + kSym^2 q[xv] == 0, q[0] == 0}, q, xv];
  qGen[xv_] = q[xv] /. First[gensol] /. {C[1] -> 0, C[2] -> 1};   (* pick the Sin branch *)
  charEq = FullSimplify[(D[qGen[xv], xv] /. xv -> lW), Assumptions -> $Assumptions];
  (* charEq is proportional to Cos[kSym lW]; smallest positive root: *)
  kWValue = First[Sort[Select[kSym /. Solve[Cos[kSym*lW] == 0 && 0 < kSym*lW < Pi, kSym], Positive]]];
  expectZero["D/N eigenvalue solves Cos[kW lW]==0", Cos[kWValue*lW]];
  expectZero["D/N eigenvalue kW = Pi/(2 lW)", kWValue - Pi/(2*lW)];
  ```
  If `DSolve` / `Solve` over the bounded range does not return cleanly in this environment, the acceptable fallback is to solve `Cos[kSym*lW] == 0` via `Reduce[Cos[kSym*lW] == 0 && 0 < kSym*lW < Pi, kSym]` and extract the unique root, OR to assert `Cos[kWValue*lW] == 0` for the closed form while explicitly solving the characteristic equation (the key requirement: `kWValue` must come from a solver acting on `Cos[k lW]==0`, NOT be typed as `Pi/(2 lW)` and then verified by substitution as the SymPy script does).
- Leave lines 48-87 unchanged; `kWValue` continues to feed `OmegaW` and downstream checks exactly as before.

Do not change the SymPy script for F1 (the SymPy derivation may remain the posit-and-verify form; only the second engine must be made independent).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 116` and confirm the new eigenvalue is obtained from a `Solve`/`Reduce`/`DSolve` on the characteristic equation (not a hand-typed `Pi/(2 lW)` later verified), that a `Cos[kW lW]==0` (or equivalent characteristic) check appears, and that the script prints all PASS lines and exits 0.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl`
- summary: Replaced the posited D/N trial eigenvalue with a DSolve-derived mode, characteristic Reduce solve, and Cos[kW lW] checks.
- deviation: none

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:64-94` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:62-81`

**Issue:** The renormalization/coefficient block is entirely unfalsifiable: (a) lines 64-69 round-trip the solved `L_W_required` back through `4 L_W^2/(pi^2 a^2)` to recover `(1+r_c)/3` (guaranteed by the solve); (b) lines 70-74 divide literals (`(1+r_c)/3`, `(1+r_c)/9`) by their own `(1+r_c)` factor; (c) lines 76-93 read the `z^2`/`z^5` coefficients back off the hand-typed ansatz `D_bare = (1+r_c)(1 - z^2/3 - i z^5/9)`. None can fail regardless of the physics.

The genuine, falsifiable physics — the D/N eigenvalue, the `kappa0` collapse, and the tube-length law — lives at lines 24-62 and must be kept exactly as-is. The renormalization to canonical coefficients (`kappa_c = 1/3`, `gamma_c = 1/9`) carries **no independent falsifiable content within this stage's scope**: `kappa0` is already established at 24-62, and `gamma0 = (1+r_c)/9` is an upstream-carried input (Stage 98 compensation requirement) with no in-stage derivation. Therefore any `expect_zero` on `kappa_c`/`gamma_c` (including a `3*kappa0 == 9*gamma0` relation between two literals) is tautological. The honest fix is to delete the falsifiable assertions at 64-93 and **report** the renormalized coefficients as labeled consequences (prints), not assert them. (This treatment was agreed between the orchestrator and Codex; do not re-introduce any `kappa_c`/`gamma_c`/`3*kappa0-9*gamma0`/`D_bare` `expect_zero` assertion.)

**Required change (SymPy):** Replace lines 64-94 (the block from the `# Round-trip…` comment through `kappa0_bare = sp.simplify((1 + r_c) / 3)`) with the following labeled-print block. Do NOT add any `expect_zero` in this block.

BEFORE (lines 64-94):
```python
# Round-trip the geometric kappa0 through L_W_required:
kappa0_bare_geom = sp.simplify(4 * L_W_required**2 / (sp.pi**2 * a**2))
expect_zero(
    "geometric kappa0 at L_W_required equals (1+r_c)/3",
    kappa0_bare_geom - (1 + r_c) / 3,
)
gamma0_bare = sp.simplify((1 + r_c) / 9)
kappa_c = sp.simplify(kappa0_bare_geom / (1 + r_c))
gamma_c = sp.simplify(gamma0_bare / (1 + r_c))
expect_zero("final kappa_c - 1/3", kappa_c - sp.Rational(1, 3))
expect_zero("final gamma_c - 1/9", gamma_c - sp.Rational(1, 9))

# Bare outgoing form as pure-scale deformation of canonical l=2 branch.
# Extract gamma0 from the z^5 coefficient and check it matches (1+r_c)/9.
D_bare = sp.expand((1 + r_c) * (1 - z**2 / 3 - sp.I * z**5 / 9))
gamma0_from_D = sp.I * D_bare.coeff(z, 5)
expect_zero(
    "gamma0 extracted from D_bare matches (1+r_c)/9",
    sp.simplify(gamma0_from_D - (1 + r_c) / 9),
)
D_final = sp.simplify(D_bare / (1 + r_c))
expect_zero(
    "bare scaled-canonical branch renormalizes to canonical",
    D_final - (1 - z**2 / 3 - sp.I * z**5 / 9),
)
kappa0_from_D = -D_bare.coeff(z, 2)
expect_zero(
    "kappa0_bare extracted from D_bare matches (1+r_c)/3",
    sp.simplify(kappa0_from_D - (1 + r_c) / 3),
)
kappa0_bare = sp.simplify((1 + r_c) / 3)
```

AFTER:
```python
# --- Renormalization to canonical coefficients (REPORTED, not asserted) ---
# The load-bearing, falsifiable physics is verified above (eigenvalue, kappa0
# collapse, tube-length law). The renormalization below carries no independent
# falsifiable content in this stage: kappa0 is already established above and
# gamma0 is an upstream-carried input (Stage 98 compensation requirement) with
# no in-stage derivation. Dividing out the common (1+r_c) hybridization factor
# is therefore a definitional consequence, so these values are PRINTED, not
# asserted (an expect_zero here would be tautological).
kappa0_bare = sp.simplify(kappa0_from_tube.subs(L_W, L_W_required))  # derived tube coeff at required length
gamma0_bare = sp.simplify((1 + r_c) / 9)                            # upstream-carried input (Stage 98), not derived in-stage
common_scale = 1 + r_c
kappa_c = sp.simplify(kappa0_bare / common_scale)
gamma_c = sp.simplify(gamma0_bare / common_scale)
print("Renormalization (definitional consequence, not an independent check):")
print("  kappa0_bare (derived tube coeff at L_W_required) =", kappa0_bare)
print("  gamma0_bare (upstream-carried input, Stage 98)   =", gamma0_bare)
print("  kappa_c = kappa0_bare/(1+r_c) =", kappa_c)
print("  gamma_c = gamma0_bare/(1+r_c) =", gamma_c)
```

Notes for Codex on the AFTER block:
- `kappa0_from_tube` and `L_W_required` are already defined above (lines 51-55); `.subs(L_W, L_W_required)` makes the provenance explicit. `kappa0_bare`, `gamma0_bare`, `kappa_c`, `gamma_c` remain defined so the existing summary block (current lines 96-99) still works — leave that summary block intact.
- Remove the `D_bare` ansatz and all its coefficient read-backs entirely. Do NOT re-add a `D_bare`/`coeff` read-back or any `expect_zero` for `kappa_c`/`gamma_c`.
- The `z` symbol declaration at the top is now unused by an assertion but may remain (do not refactor the symbol declarations).

**Required change (Mathematica):** mirror the SymPy change — delete the assertions at lines 62-81 and replace with labeled prints. Do NOT add any `expectZero` in this block.

BEFORE (lines 62-81):
```
kappa0BareGeom = FullSimplify[4*lWRequired^2/(Pi^2*a^2), Assumptions -> $Assumptions];
expectZero["geometric kappa0 at lWRequired equals (1+r_c)/3",
           kappa0BareGeom - (1 + rC)/3];
gamma0Bare = FullSimplify[(1 + rC)/9, Assumptions -> $Assumptions];
kappaC = FullSimplify[kappa0BareGeom/(1 + rC), Assumptions -> $Assumptions];
gammaC = FullSimplify[gamma0Bare/(1 + rC), Assumptions -> $Assumptions];

expectZero["final kappa_c - 1/3", kappaC - 1/3];
expectZero["final gamma_c - 1/9", gammaC - 1/9];

dBare = Expand[(1 + rC)*(1 - z^2/3 - I*z^5/9)];
gamma0FromD = FullSimplify[I * Coefficient[dBare, z, 5], Assumptions -> $Assumptions];
expectZero["gamma0 extracted from dBare matches (1+rC)/9",
           gamma0FromD - (1 + rC)/9];
dFinal = FullSimplify[dBare/(1 + rC), Assumptions -> $Assumptions];
expectZero["bare scaled-canonical branch renormalizes to canonical", dFinal - (1 - z^2/3 - I*z^5/9)];
kappa0FromD = FullSimplify[-Coefficient[dBare, z, 2], Assumptions -> $Assumptions];
expectZero["kappa0_bare extracted from dBare matches (1+rC)/3",
           kappa0FromD - (1 + rC)/3];
kappa0Bare = FullSimplify[(1 + rC)/3, Assumptions -> $Assumptions];
```

AFTER:
```
(* --- Renormalization to canonical coefficients (REPORTED, not asserted) --- *)
(* Load-bearing physics verified above; gamma0 is an upstream-carried input    *)
(* (Stage 98), so kappa_c/gamma_c here are definitional consequences, printed  *)
(* not asserted (an expectZero here would be tautological).                    *)
kappa0Bare = FullSimplify[4*lWRequired^2/(Pi^2*a^2), Assumptions -> $Assumptions];  (* derived tube coeff at required length *)
gamma0Bare = FullSimplify[(1 + rC)/9, Assumptions -> $Assumptions];                  (* upstream-carried input (Stage 98) *)
commonScale = 1 + rC;
kappaC = FullSimplify[kappa0Bare/commonScale, Assumptions -> $Assumptions];
gammaC = FullSimplify[gamma0Bare/commonScale, Assumptions -> $Assumptions];
Print["Renormalization (definitional consequence, not an independent check):"];
Print["  kappa0_bare (derived tube coeff at lWRequired) = ", kappa0Bare];
Print["  gamma0_bare (upstream-carried input, Stage 98) = ", gamma0Bare];
Print["  kappa_c = kappa0Bare/(1+rC) = ", kappaC];
Print["  gamma_c = gamma0Bare/(1+rC) = ", gammaC];
```

Notes for Codex: the summary block (current wl lines 83-90) references `lWRequired`, `kappa0Bare`, `gamma0Bare`, `kappaC`, `gammaC`, all still defined; leave it intact. Remove the `dBare` ansatz and its read-backs entirely. The `z` symbol declaration may remain.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 116` and `redteam exec-mathematica 116` and confirm: (1) no `kappa0_bare_geom`/`kappa0BareGeom` round-trip assertion from `L_W_required`/`lWRequired` remains; (2) no `D_bare`/`dBare` ansatz or coefficient read-back remains; (3) NO `expect_zero`/`expectZero` on `kappa_c`/`gamma_c` or a `3*kappa0-9*gamma0` relation remains (the renormalization is now labeled prints by design — this is the agreed non-tautological treatment, not an omission); (4) the genuine checks at lines 24-62 (eigenvalue, kappa0 collapse, tube-length law) are unchanged and still PASS; (5) both scripts print the renormalization values and exit 0.

## Applied: F2

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl`
- summary: Replaced tautological renormalization assertions and D_bare/dBare coefficient read-backs with labeled consequence prints in both scripts.
- deviation: none
