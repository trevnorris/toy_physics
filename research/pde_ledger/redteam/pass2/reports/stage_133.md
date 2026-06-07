---
unit_id: 133
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage133_coupled_mouth_fixedpoint.md]
  paper_appendix: present
---

# Audit unit 133 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_133.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage133_coupled_mouth_fixedpoint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 30; block-contract/status at lines 22-42)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.txt`

## What the paper claims

Stage 133 replaces the Stage-129 effective parent datum with the first explicit *coupled* two-channel mouth-layer solve, reducing the mouth-bias ambiguity to a small set of dimensionless gains/stiffnesses. The coupled system `[-∂ₓ²I + K]U = Σ_Π(x)G` with D/N conditions `U(0)=0, U'(1)=0` diagonalizes (`K = R diag(κ₊²,κ₋²) Rᵀ`) into two scalar D/N response problems `(-∂ₓ²+κ²)u = G·Σ_Π`, `u(0)=0, u'(1)=0`. The boxed deliverables are: (i) the exact scalar solution `u(x)=A sinh(κx)−C cosh(κx)+C e^{−Πx}` with `C=GΠ/((1−e^{−Π})(κ²−Π²))`, `A=C(κ sinh κ + Π e^{−Π})/(κ cosh κ)`; (ii) the exact mouth-response kernel `S(Π,κ)=Π[κ tanh κ + Π(e^{−Π}sech κ − 1)]/((1−e^{−Π})(κ²−Π²))` with `u'(0)=G·S(Π,κ)`; (iii) the static-shell limit `S(Π,0)=1`; and (iv) the closed-form fixed-point law `Π = M₊S(Π,κ₊) + M₋S(Π,κ₋)`, `M_α := L H_α G_α / Θ_σ`. The `\stagefield{Verification}` line points at exactly these two scripts; the card-body quote is `Π=Σ_{α=±} M_α S(Π,κ_α)`. The card's `\stagefield{Checks}` list (gain pair vs outlet consistency, self-matched susceptibility closure, numerical fixed points recorded as numerically-located) are forward-looking block checklist items, not this stage's script assertions — the card's Purpose states the audit target is "the verification output quoted below," i.e. the D/N kernel reduction.

## What the script claims to verify

Both scripts verify the scalar D/N response branch and its kernel, then print (without asserting) the recombined two-channel fixed-point law. The SymPy script writes the paper's closed forms `C`, `A`, `u` as literals and then independently checks them against the physics: the ODE residual `−u'' + κ²u − G·Σ`, the boundary conditions `u(0)`, `u'(1)`, and the kernel identity `u'(0)/G − S_target`, all asserted = 0; it then takes `lim_{κ→0} S_target` and asserts it equals 1. The Mathematica script instead solves the ODE from scratch with `DSolveValue`, then runs the identical residual/BC/kernel/limit checks. Neither script asserts the four-parameter fixed-point law (M_α, H_α, G_α, R are free symbols), which is a definitional recombination, not a falsifiable identity in this symbol set; both only print it.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (i) `u(x)` solves `(-∂ₓ²+κ²)u=G·Σ_Π`, `u(0)=0`, `u'(1)=0` | sympy ODE-residual/bc0/bc1 (lines 41-43,53-55); wl DSolveValue + residual/bc0/bc1 (lines 36-45,54-56) | match |
| (ii) kernel `S(Π,κ)` with `u'(0)=G·S` | sympy `S - S_target == 0` (lines 45-49,56); wl `sKernel - sTarget == 0` (lines 47-51,57) | match |
| (iii) static-shell limit `S(Π,0)=1` | sympy `lim_{κ→0} S_target == 1` (lines 58-62); wl `Series[...,{kk,0,0}] − 1 == 0` (lines 59-63) | match |
| (iv) `Π = M₊S(Π,κ₊)+M₋S(Π,κ₋)` | print-only, both engines (sympy 64-68; wl 65-70) | partial (definitional restatement; not an algebraically falsifiable identity within this stage's free symbols) |
| card `\stagefield{Checks}` (gains vs outlet, susceptibility closure, numeric fixed points) | none | n/a — forward-looking block checklist, not this stage's audit target |

`paper_alignment` = aligned: every algebraically-checkable deliverable (i–iii) is faithfully exercised by both engines; (iv) is a definitional recombination correctly printed, consistent with the card's stated audit target.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero(−u'' + κ²u − G·Σ)` | claim (i) ODE | yes |
| A2 | sympy | 54 | `expect_zero(u(0))` | claim (i) Dirichlet BC | yes |
| A3 | sympy | 55 | `expect_zero(u'(1))` | claim (i) Neumann BC | yes |
| A4 | sympy | 56 | `expect_zero(S − S_target)` | claim (ii) kernel | yes |
| A5 | sympy | 61-62 | `if simplify(S0−1)!=0: raise` | claim (iii) static limit | yes |
| A6 | math | 54 | `expectZero[ODE residual]` | claim (i) ODE | yes |
| A7 | math | 55 | `expectZero[u(0)]` | claim (i) Dirichlet BC | yes |
| A8 | math | 56 | `expectZero[u'(1)]` | claim (i) Neumann BC | yes |
| A9 | math | 57 | `expectZero[sKernel − sTarget]` | claim (ii) kernel | yes |
| A10 | math | 63 | `expectZero[s0 − 1]` | claim (iii) static limit | yes |

All ten assertions are non-tautological and trace to a specific paper deliverable. The SymPy `u`/`C`/`A` are written as paper-stated literals but are NOT self-checked — they are validated against the *independent* physical constraints (ODE + both BCs), so they cannot pass if the form is wrong. The kernel `S_target` is likewise a paper literal validated against the independently-obtained `u'(0)/G`. In Mathematica the solution itself is regenerated by `DSolveValue`, so even the `u` form is independently derived.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. The SymPy script hardcodes the closed-form `u`, `C`, `A` (sympy lines 36-39) and verifies them against the ODE; the Mathematica script *does not reproduce that algebra* — it solves the BVP from scratch with `DSolveValue[{-uFun''[x]+kappa^2 uFun[x]==gSrc sigma, uFun[0]==0, uFun'[1]==0}, uFun[x], x]` (wl lines 36-40). The two solution representations differ textually: SymPy emits `-G·Π(Π e^{Πx} sinh(κx) + κ e^Π cosh κ − κ e^{Π(x+1)} cosh(κ(x−1)))e^{−Πx} / (κ(Π²−κ²)(e^Π−1)cosh κ)` (sympy output lines 6-10), whereas Mathematica's `DSolveValue` emits a six-term exponential form `(gSrc piM (E^(piM+kappa x) kappa − … + E^(kappa+2kappa x+piM x) piM)) / (E^((kappa+piM)x)(1+E^(2kappa))(−1+E^piM) kappa (kappa−piM)(kappa+piM))` (math output line 5). These are algebraically equal (both pass identical residual/BC/kernel checks) but were generated by genuinely different routes — strong evidence of an independent second engine. The kernel target `sTarget` is the same paper literal in both, but each engine validates it against its own `u'(0)/gSrc`.

## Engine cross-check

Both engines emit `ODE residual = 0`, `u(0)=0`, `u'(1)=0`, `mouth derivative kernel = 0`, `S(Pi,0)=1`, all PASS. The printed `S(Π,κ)` forms agree: SymPy `−Π(−Π e^Π + Π/cosh κ + κ e^Π tanh κ)/((Π²−κ²)(e^Π−1))` (sympy output lines 23-28) equals Mathematica `(piM(−piM Sech[kappa] + E^piM(piM − kappa Tanh[kappa])))/((−1+E^piM)(−kappa²+piM²))` (math output line 24) — the same rational-trig function (note `Π/cosh κ = Π sech κ`, and the leading-sign / denominator-sign bookkeeping matches). Engines agree.

## Verdict justification

Clean. I attacked four ways and all failed: (1) Tautology — the hardcoded `u`/`C`/`A`/`S_target` are NOT asserted against themselves; they are checked against the independent ODE+BC physics (and against a from-scratch `DSolveValue` solution in Mathematica), so a wrong form would fail the residual or kernel check. (2) Domain/pole — symbols are `positive, real` with `κ ≠ Π` enforced in Mathematica's `$Assumptions` (wl line 31, `kappa != piM`) and structurally the `(κ²−Π²)` pole is never hit because the kernel is the limit-consistent form; the static `κ→0` limit is taken via `sp.limit`/`Series`, not by substitution into the singular `tanh`/`(κ²−Π²)` expression, so no 0/0 is forced. (3) Static-limit correctness — I confirmed analytically `lim_{κ→0} S = Π²(e^{−Π}−1) / (Π²(e^{−Π}−1)) = 1`. (4) Engine independence — confirmed via the textually-distinct DSolve solution. I read the card, notes (§§1-4 + Result), and the appendix row (line 30); the scripts verify exactly the card's stated audit target (the two-channel D/N kernel reduction `Π=Σ M_α S(Π,κ_α)`). The only non-asserted deliverable, the 4-parameter fixed-point law, is a definitional recombination with free gains and is correctly print-only; flagging it would be gold-plating, not a real defect.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level table (symbolic forms; this stage emits no numeric figures-of-merit):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `u(x) = A sinh(κx) − C cosh(κx) + C e^{−Πx}` (with `C`, `A` forms) | py lines 36-39, py output lines 6-10; wl DSolveValue line 36, wl output line 5 | notes lines 129-138 (boxed `C`, `A`); card body via Verification quote | MATCH |
| mouth kernel `S(Π,κ)=Π[κ tanh κ + Π(e^{−Π}sech κ − 1)]/((1−e^{−Π})(κ²−Π²))` | py lines 46-49, py output lines 23-28; wl lines 48-51, wl output line 24 | notes lines 146-152 (boxed); restated in Result | MATCH |
| static-shell limit `S(Π,0)=1` | py lines 59-60, py output line 15; wl lines 61-63, wl output line 14 | notes lines 158-160 (boxed) | MATCH |
| fixed-point law `Π = M₊S(Π,κ₊) + M₋S(Π,κ₋)`, `M_α := L H_α G_α / Θ_σ` | py lines 65-68, py output line 20; wl lines 67-70, wl output line 21 | notes lines 178-187 + Result lines 207-213 (boxed); card body line 16 quote | MATCH |

INTERNAL (scaffolding, no prose expected): `ODE residual`, `u(0)`, `u'(1)`, `mouth derivative kernel` (all residual=0 pass-flags), `static-shell limit` residual (=0 flag), PASS/FAIL strings.

reconciliation: complete; 4 deliverable values checked, 0 misaligned.

## Self-test notes

I checked the variable-independence trap (the `sp.diff(u, x)` derivatives all act on a `u` that genuinely depends on `x`, so `u'(0)`, `u'(1)` are not identically zero — confirmed by the nontrivial printed `u(x)`). Parity is N/A (no unbounded-domain integrals; the BVP is on `[0,1]`). Trivial-case pre-check: I analytically reduced the `κ→0` numerator/denominator and confirmed `S(Π,0)=1` (both → `Π²(e^{−Π}−1)`), and confirmed the ODE-residual `assert_zero` is non-trivial because `u` is a genuine BVP solution, not a constant. Outputs are fresh (output mtimes 17:45/17:47 > script mtimes 17:20/17:22), so no `stale_output`. No directive written (zero findings).
