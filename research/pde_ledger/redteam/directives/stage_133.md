---
unit_id: 133
batch: IV.4
created_at: 2026-05-27T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 133

Apply the finding below. After applying, append an `## Applied: F1` block under
the finding with: `files_changed`, `summary` (one sentence), and `deviation`
(or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is
specified. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:33-52`

**Issue:**
The Mathematica audit duplicates the SymPy script's hand-ansatz algebra: it
introduces the same `cCoeff` and `aCoeff` derived in the notes and assembles
the same closed-form `u` as a literal expression. It then verifies the ODE and
BCs against that hand-supplied `u`. The second engine therefore offers no
independent algebraic verification of `u`; an error in the hand-derived `A` or
`C` would propagate identically to both engines and both would pass. Replace
the hand-ansatz construction with an independent derivation via `DSolveValue`
so Mathematica genuinely re-derives the closed-form solution from the PDE and
boundary conditions.

**Required change:**

Replace the block at lines 28-52 of
`mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl`
(everything from `Clear[x, piM, kappa, gSrc];` through the
`expectZero["mouth derivative kernel", sKernel - sTarget];` line) with the
following independent-derivation block.

Before:
```
Clear[x, piM, kappa, gSrc];
$Assumptions =
  Element[{x, piM, kappa, gSrc}, Reals] &&
  0 <= x <= 1 && piM > 0 && kappa > 0 && gSrc > 0 && kappa != piM;

sigma = piM*Exp[-piM*x]/(1 - Exp[-piM]);
cCoeff = FullSimplify[gSrc*piM/((1 - Exp[-piM])*(kappa^2 - piM^2)), Assumptions -> $Assumptions];
aCoeff = FullSimplify[cCoeff*(kappa*Sinh[kappa] + piM*Exp[-piM])/(kappa*Cosh[kappa]), Assumptions -> $Assumptions];
u = FullSimplify[aCoeff*Sinh[kappa*x] - cCoeff*Cosh[kappa*x] + cCoeff*Exp[-piM*x], Assumptions -> $Assumptions];

residual = FullSimplify[-D[u, {x, 2}] + kappa^2*u - gSrc*sigma, Assumptions -> $Assumptions];
bc0 = FullSimplify[u /. x -> 0, Assumptions -> $Assumptions];
bc1 = FullSimplify[D[u, x] /. x -> 1, Assumptions -> $Assumptions];

sKernel = FullSimplify[(D[u, x] /. x -> 0)/gSrc, Assumptions -> $Assumptions];
sTarget = FullSimplify[
  piM*(kappa*Tanh[kappa] + piM*(Exp[-piM]/Cosh[kappa] - 1))/((1 - Exp[-piM])*(kappa^2 - piM^2)),
  Assumptions -> $Assumptions
];

Print["u(x) = ", fmt[u]];
expectZero["ODE residual", residual];
expectZero["u(0)", bc0];
expectZero["u'(1)", bc1];
expectZero["mouth derivative kernel", sKernel - sTarget];
```

After:
```
Clear[x, piM, kappa, gSrc, uFun];
$Assumptions =
  Element[{x, piM, kappa, gSrc}, Reals] &&
  0 <= x <= 1 && piM > 0 && kappa > 0 && gSrc > 0 && kappa != piM;

sigma = piM*Exp[-piM*x]/(1 - Exp[-piM]);

(* Independent derivation: let DSolveValue solve the D/N problem from scratch. *)
uSol = DSolveValue[
  {-uFun''[x] + kappa^2*uFun[x] == gSrc*sigma, uFun[0] == 0, uFun'[1] == 0},
  uFun[x],
  x
];
u = FullSimplify[uSol, Assumptions -> $Assumptions];

residual = FullSimplify[-D[u, {x, 2}] + kappa^2*u - gSrc*sigma, Assumptions -> $Assumptions];
bc0 = FullSimplify[u /. x -> 0, Assumptions -> $Assumptions];
bc1 = FullSimplify[D[u, x] /. x -> 1, Assumptions -> $Assumptions];

sKernel = FullSimplify[(D[u, x] /. x -> 0)/gSrc, Assumptions -> $Assumptions];
sTarget = FullSimplify[
  piM*(kappa*Tanh[kappa] + piM*(Exp[-piM]/Cosh[kappa] - 1))/((1 - Exp[-piM])*(kappa^2 - piM^2)),
  Assumptions -> $Assumptions
];

Print["u(x) = ", fmt[u]];
expectZero["ODE residual", residual];
expectZero["u(0)", bc0];
expectZero["u'(1)", bc1];
expectZero["mouth derivative kernel", sKernel - sTarget];
```

Key differences:
- `cCoeff` and `aCoeff` are removed.
- `u` is now obtained by `DSolveValue[{ODE, BCs}, uFun[x], x]` — Mathematica
  derives the closed form from the PDE itself rather than copying the SymPy
  ansatz.
- All four downstream `expectZero` assertions are unchanged and now serve as
  independent verification that the `DSolve` answer satisfies the ODE, both
  BCs, and reproduces the paper's boxed `S(Pi, kappa)` kernel via `u'(0)/G`.

Do NOT change anything else in the file. Leave the static-shell limit block
(currently lines 54-58 in the original), the banner sections, the final two-
channel print block, and the `Exit[0]` line untouched.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 133` and
confirm:
1. The script exits 0.
2. The script no longer contains the identifiers `cCoeff` or `aCoeff`
   (grep both literals; both should be absent).
3. The script contains a `DSolveValue` (or `DSolve`) call with the inhomogeneous
   ODE and both BCs.
4. All four `expectZero` lines for `ODE residual`, `u(0)`, `u'(1)`, and
   `mouth derivative kernel` are present and pass.
5. The static-shell limit block still prints `S(Pi,0) = 1` and the
   `static-shell limit` `expectZero` still passes.
