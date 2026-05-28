---
unit_id: 133
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage133_coupled_mouth_fixedpoint.md
  paper_appendix: present
---

# Audit unit 133 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_133.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage133_coupled_mouth_fixedpoint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.txt`

## What the paper claims

Stage 133 derives the explicit coupled mouth-layer fixed-point law. The body
equation in the card states the two-channel boundary-layer problem reduces (via
diagonalisation of the stiffness `K`) to the fixed-point relation
"\(\Pi=\sum_\alpha M_\alpha\mathcal S(\Pi,\kappa_\alpha)\)" with the scalar D/N
mouth-response kernel `S(Pi, kappa)` given in closed form. The notes spell out
the derivation: (i) the dimensionless coupled mouth-layer system with D/N BCs
`U(0)=0`, `U'(1)=0`; (ii) exact diagonalisation into two scalar D/N problems;
(iii) the closed-form solution `u(x) = A sinh(kx) - C cosh(kx) + C e^{-Pi x}`
with explicit `A`, `C`; (iv) the boxed kernel
`S(Pi, kappa) = Pi[kappa tanh(kappa) + Pi(e^{-Pi} sech(kappa) - 1)] / ((1 - e^{-Pi})(kappa^2 - Pi^2))`;
and (v) the static-shell limit `S(Pi, 0) = 1`. The card's `Purpose` field
declares "its audit target is the verification output quoted below" — i.e., the
boxed `Pi = sum M_alpha S(Pi, kappa_alpha)` two-channel fixed-point law.

## What the script claims to verify

Both scripts solve the scalar D/N problem with exponential source `Sigma(x)` and
verify four assertions: (1) the ODE residual `(-d^2/dx^2 + kappa^2) u - G*Sigma`
vanishes for the closed-form `u`; (2) `u(0) = 0`; (3) `u'(1) = 0`; (4) the
computed mouth-derivative kernel `u'(0)/G` equals the boxed closed form
`S(Pi, kappa)`. They additionally check the static-shell limit `S(Pi, 0) = 1`.
The general two-channel fixed-point law is **printed** (banner block "GENERAL
TWO-CHANNEL FIXED-POINT LAW") but is **not** wrapped in an `assert` /
`expectZero` — both scripts treat it as the formal linear combination of the
already-verified per-channel kernel.

## Paper ↔ script cross-check

| Paper-side deliverable | Script coverage |
|---|---|
| Closed-form `u(x)` solves the scalar D/N problem (ODE + BCs) | match (SymPy lines 41-43, Mathematica lines 38-40) |
| Boxed `S(Pi, kappa)` kernel formula | match (SymPy lines 45-56, Mathematica lines 42-52) |
| Static-shell limit `S(Pi, 0) = 1` | match (SymPy lines 58-62, Mathematica lines 54-58) |
| Diagonalisation `K = R diag(kappa+^2, kappa-^2) R^T` | extra (treated as definition; not formally verified, but not required since it's by construction) |
| Boxed `Pi = sum_alpha M_alpha S(Pi, kappa_alpha)` (the stated audit target) | partial — printed only, not asserted (acceptable: linear combination of per-channel kernel after diagonalisation is mathematically trivial) |
| Card `Checks`: gain pair against outlet consistency, susceptibility closure, numerical fixed points | not covered here — but these are forward-looking checklist items, not in-stage verification targets per the card's `Purpose` field |

`paper_alignment: partial` — the per-channel kernel that the boxed two-channel
law is built from is fully verified; the two-channel assembly itself is not
asserted but is a trivial linear combination once `S` is verified. No identity
mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero("ODE residual", res)` where `res = -u'' + kappa^2 u - G Sigma` | scalar D/N ODE | yes |
| A2 | sympy | 54 | `expect_zero("u(0)", bc0)` | Dirichlet BC at x=0 | yes |
| A3 | sympy | 55 | `expect_zero("u'(1)", bc1)` | Neumann BC at x=1 | yes |
| A4 | sympy | 56 | `expect_zero("mouth derivative kernel", S - S_target)` | boxed `S(Pi, kappa)` formula | yes |
| A5 | sympy | 61-62 | `if sp.simplify(S0 - 1) != 0: raise` | static-shell limit `S(Pi,0)=1` | yes |
| A6 | mathematica | 49 | `expectZero["ODE residual", residual]` | scalar D/N ODE | yes |
| A7 | mathematica | 50 | `expectZero["u(0)", bc0]` | Dirichlet BC at x=0 | yes |
| A8 | mathematica | 51 | `expectZero["u'(1)", bc1]` | Neumann BC at x=1 | yes |
| A9 | mathematica | 52 | `expectZero["mouth derivative kernel", sKernel - sTarget]` | boxed `S(Pi, kappa)` formula | yes |
| A10 | mathematica | 58 | `expectZero["static-shell limit", s0 - 1]` | static-shell limit `S(Pi,0)=1` | yes |

All listed assertions are non-tautological: they verify the manually constructed
`u(x)` satisfies the ODE and BCs (a genuine algebraic identity in the
hyperbolic + exponential basis), and that the algebraic simplification of
`u'(0)/G` reproduces the boxed `S(Pi, kappa)` formula. Both could fail if `A`,
`C`, or `S_target` were misstated.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:33-52`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.py:34-56`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script's algebra, not an
independent re-derivation. Every load-bearing intermediate is the hand-derived
ansatz from the notes, copied identically across both engines:

SymPy (lines 34-39):
```
Sigma = Pi * sp.exp(-Pi*x) / (1 - sp.exp(-Pi))
C = G * Pi / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
A = sp.simplify(C * (kappa * sp.sinh(kappa) + Pi * sp.exp(-Pi)) / (kappa * sp.cosh(kappa)))
u = sp.simplify(A * sp.sinh(kappa*x) - C * sp.cosh(kappa*x) + C * sp.exp(-Pi*x))
```

Mathematica (lines 33-36):
```
sigma = piM*Exp[-piM*x]/(1 - Exp[-piM]);
cCoeff = FullSimplify[gSrc*piM/((1 - Exp[-piM])*(kappa^2 - piM^2)), ...];
aCoeff = FullSimplify[cCoeff*(kappa*Sinh[kappa] + piM*Exp[-piM])/(kappa*Cosh[kappa]), ...];
u = FullSimplify[aCoeff*Sinh[kappa*x] - cCoeff*Cosh[kappa*x] + cCoeff*Exp[-piM*x], ...];
```

Variable choreography, intermediate names, and the explicit closed form of `u`
are identical. The same `sTarget` literal is constructed by the same algebraic
expression. Mathematica never calls `DSolve` on
`-u''[x] + kappa^2 u[x] == gSrc*sigma[x]` with the D/N BCs to obtain `u`
independently, then test equality against the SymPy form. Both engines verify
that a hand-supplied `u` satisfies the ODE/BCs and that a hand-supplied
`sTarget` matches the manually computed `u'(0)/G` from that same `u`. The
second engine therefore contributes no independent algebraic check — it
reproduces the first engine's algebra in a different syntax.

**Why this matters:**
Per the second-engine policy, both engines must derive the result independently
from the physical premises so that a coding error or symbol-domain mistake in
one engine can be caught by the other. Here, an error in the hand-derived `A`
or `C` would propagate to both engines unchanged and both would "PASS." The
checkpoint trust audit treats a transliterated second engine as effectively
absent cross-checking.

**Required change:**
Refactor the Mathematica script so it derives the closed-form `u(x)` by an
independent path. The minimal acceptable approach is:

1. Set up the ODE symbolically:
   `ode = -u''[x] + kappa^2 u[x] == gSrc*piM*Exp[-piM*x]/(1 - Exp[-piM])`.
2. Use `DSolveValue[{ode, u[0] == 0, u'[1] == 0}, u[x], x]` (or `DSolve` with
   the same BCs) to obtain `u[x]` from first principles.
3. Compute `kernelDerived = FullSimplify[(D[uSol, x] /. x -> 0)/gSrc]`.
4. Assert `expectZero["kernel matches paper formula", kernelDerived - sTarget]`
   where `sTarget` is the boxed paper formula.
5. Independently verify ODE residual and BCs against `uSol` (these become
   sanity checks that `DSolve` returned the correct solution; not duplicates
   of the SymPy assertions because the source of `uSol` is now `DSolve` rather
   than the hand ansatz).

This removes all hand-imported `cCoeff`, `aCoeff`, and the manual `u`
construction from the Mathematica side; the SymPy script can retain its
hand-constructed `u` since it serves as the explicit-form witness.

**Verification:**
After Codex applies, the new Mathematica script must contain a `DSolve` (or
`DSolveValue`) call producing the closed-form `u(x)`, and the
"mouth derivative kernel" `expectZero` must still pass. The verifier runs
`redteam exec-mathematica 133` and confirms the script exits 0 and that the
`u(x)` printed line is structurally identical (up to `FullSimplify`) to the
current output. Grepping the script for `cCoeff` and `aCoeff` should return
no matches (they should be eliminated in favour of the `DSolve` solution).

## Independent-derivation check (Mathematica)

The `.wl` script is a transliteration of the `.py` script. Both build `u(x)`
from the same hand-supplied coefficients `A` and `C` from the notes; neither
engine derives `u(x)` from `DSolve` / `dsolve` / Green's function /
variation-of-parameters / Laplace transform. The Mathematica script's
`cCoeff` and `aCoeff` are direct ports of the SymPy `C` and `A`. See F1 for
the side-by-side comparison.

## Engine cross-check

Both engines report the same printed `u(x)` (modulo simplification form) and
the same printed `S(Pi, kappa)`. SymPy output prints `S(Pi,0) = 1`; Mathematica
prints `S(Pi,0) = 1` then `static-shell limit = 0`. All five corresponding
assertion pairs (`ODE residual`, `u(0)`, `u'(1)`, `mouth derivative kernel`,
`S(Pi, 0) = 1` / static-shell limit) report identical zeros. No
engine-level numerical or symbolic disagreement. (However, because of F1, this
agreement is not informative: both engines test the same hand-supplied
expression.)

## Verdict justification

The math holds up: I verified by hand that the proposed `u(x) = A sinh(kx) -
C cosh(kx) + C e^{-Pi x}` with `A = C(k sinh k + Pi e^{-Pi})/(k cosh k)` and
`C = G Pi / ((1-e^{-Pi})(k^2 - Pi^2))` satisfies the ODE
`(-d_x^2 + k^2) u = G Sigma`, satisfies `u(0) = 0` (since `-C + C = 0`),
satisfies `u'(1) = 0` (substitution gives `A k cosh(k) = C k sinh(k) +
C Pi e^{-Pi}`, exactly the definition of `A`), produces
`u'(0)/G = (Pi/((1-e^{-Pi})(k^2 - Pi^2))) [k tanh(k) + Pi(e^{-Pi} sech(k) - 1)]`
which matches the boxed `S`, and limits to `S(Pi, 0) = (-Pi^2(1-e^{-Pi}))/
(-Pi^2(1-e^{-Pi})) = 1`. The assertions are non-tautological and the assertion
inventory exercises every load-bearing piece of the paper's boxed kernel and
its static limit. The only finding is the second-engine transliteration: both
engines verify the same hand-supplied `u(x)` rather than independently
deriving it, so the Mathematica script does not provide the cross-engine
algebraic redundancy the policy requires. Verdict: `findings` (one
script-side finding), `stop_cold: null`.

## Self-test notes

I hand-traced (i) the ODE residual on each of the three building blocks
`A sinh(kx)`, `-C cosh(kx)`, `C e^{-Pi x}` and confirmed the first two annihilate
and the third produces `C(k^2 - Pi^2) e^{-Pi x} = G Pi e^{-Pi x}/(1-e^{-Pi})`,
matching `G Sigma`; (ii) the two boundary conditions; (iii) the kappa->0 limit
of the kernel (both numerator and denominator vanish like `Pi^2(1-e^{-Pi})`,
giving the unit limit); (iv) checked that `simplify` calls under
`positive=True` / `kappa > 0, piM > 0` assumptions cannot hide branch errors
because all hyperbolic functions are entire and the only pole is at
`kappa = piM` (excluded by `kappa != piM` in Mathematica's assumptions and
implicitly by the symbol declaration in SymPy). The directive's prescribed
`DSolve` approach was self-tested mentally: with the source `gSrc*sigma[x]`
and BCs `u[0] == 0, u'[1] == 0`, Mathematica's `DSolve` returns the same
analytic solution up to `FullSimplify`, so the equality check
`kernelDerived - sTarget == 0` will hold.
