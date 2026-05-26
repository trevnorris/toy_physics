---
unit_id: 017
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 017 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_017.tex`
- notes: (none)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.txt`

## What the paper claims

Stage 017 transports the parent throat action into the weak-axisymmetric grouped P_2 lane structure. The paper card (`paper/stages/stage_017.tex:14-26`) states three load-bearing deliverables:

1. A Y_{20} wall perturbation yields the grouped-lane signature `(lambda_{20}, lambda_{21}, lambda_{22}) = (1, 1/2, -1)` (eq:stage017-wall-signature). Explicitly: "The audit derives \eqref{eq:stage017-wall-signature} from the actual triple-overlap coefficients."
2. The grouped trace/anomaly variables obey `b = 3a` (eq:stage017-b-equals-3a).
3. Wall-only obstruction: "Pure wall anisotropy closes the even gates only on the trivial branch delta K = delta M = 0. Thus parent promotion supplies the wall-side origin of the pattern, but does not by itself realize the full weak-axisymmetric branch."

The `\stagefield{Output}` packages these as: "Stage 017 exports the parent weak-axisymmetric lane law (eq:stage017-wall-signature)-(eq:stage017-b-equals-3a) and the wall-only gate obstruction." The Part I appendix row (`stage_appendix_part01.tex:56`) confirms the row text: "Weak-axisymmetric parent-action transport law." No notes file exists for stage 017.

## What the script claims to verify

The SymPy script (1) derives the lane factors `lambda_{20}, lambda_{21}, lambda_{22}` from `gaunt(2,2,2,0,m,-m)` ratios with a `(-1)^m` real-harmonic phase, including a same-sign cross-term vanishing pre-check for m=1,2; (2) constructs the grouped trace `Xbar` and anomaly variables `a, b` from `(X_{20}, X_{21}, X_{22}) = (eps lambda_A X_1)` for both wall inertia (M_1) and wall stiffness (K_{1w}); (3) assembles wall-only contributions to the live K_1 and H_even gates, cross-checks them lane-by-lane against generic formulas `K_{1,wall} = -delta M + delta K/9` and `H_{even,wall} = (2/3) delta M - delta K/27`, computes the 2x2 Jacobian determinant (= 1/27, with a mutated-determinant nonzero discriminator at line 110), and confirms the linear homogeneous system has only the trivial solution `{delta K = 0, delta M = 0}`; (4) extends the b=3a check to `Xi_load` and the prefactor shift `delta P_0`. The Mathematica script reproduces the same identities, but the lane-factor derivation goes through direct sphere integration of three `SphericalHarmonicY` factors (`tripleOnSphere`, lines 23-44) rather than through SymPy's Gaunt-symbol library — an independent path to the same overlap integrals.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `(lambda_{20}, lambda_{21}, lambda_{22}) = (1, 1/2, -1)` derived from triple overlap (eq:stage017-wall-signature) | SymPy lines 40-45 (`gaunt` ratios) + Mathematica lines 51-55 (sphere integration) | match |
| `b = 3a` for grouped trace/anomaly (eq:stage017-b-equals-3a) | SymPy lines 57-60, 125-128 + Mathematica lines 107-114 | match |
| Wall-only obstruction: only trivial branch `delta K = delta M = 0` | SymPy lines 99-110, 123-124 + Mathematica lines 61-69 | match |

`paper_alignment` = aligned. No paper deliverable is missing a script check; no script check is orphaned from the paper.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43 | `assert_zero('Y20 overlap lane 20', lam20 - 1)` | claim 1 (lambda_{20}=1) | yes |
| A2 | sympy | 44 | `assert_zero('Y20 overlap lane 21', lam21 - 1/2)` | claim 1 (lambda_{21}=1/2) | yes |
| A3 | sympy | 45 | `assert_zero('Y20 overlap lane 22', lam22 + 1)` | claim 1 (lambda_{22}=-1) | yes |
| A4 | sympy | 22-24 | `same_sign != 0` raises for m=1,2 (selection rule) | claim 1 selection-rule guard | yes |
| A5 | sympy | 57 | `assert_zero('wall inertia trace', Mbar)` | claim 2 (Mbar = 0) | yes |
| A6 | sympy | 58 | `assert_zero('wall inertia b=3a', bM - 3*aM)` | claim 2 (b=3a for M) | yes |
| A7 | sympy | 59 | `assert_zero('wall stiffness trace', Kbar)` | claim 2 (Kbar = 0) | yes |
| A8 | sympy | 60 | `assert_zero('wall stiffness b=3a', bK - 3*aK)` | claim 2 (b=3a for K) | yes |
| A9 | sympy | 92-94 | generic K1 vs lane gate, m=0,1,2 | claim 3 (lane structure of K1) | yes |
| A10 | sympy | 95-97 | generic Hev vs lane gate, m=0,1,2 | claim 3 (lane structure of Hev) | yes |
| A11 | sympy | 109 | `wall_matrix.det() - 1/27` | claim 3 (nondegenerate gate) | yes |
| A12 | sympy | 110 | `assert_nonzero` on mutated determinant | discrimination guard for A11 | yes |
| A13 | sympy | 123-124 | `sol_even == [{dKsym: 0, dMsym: 0}]` | claim 3 (trivial branch only) | yes |
| A14 | sympy | 125-126 | `Xibar`, `bXi - 3*aXi` | claim 2 extended to Xi | yes |
| A15 | sympy | 127-128 | `Pbar`, `bP - 3*aP` | claim 2 extended to prefactor | yes |
| M1 | math | 51 | `laneFactor[0] - 1` via direct sphere integration | claim 1 (lambda_{20}=1) | yes |
| M2 | math | 52 | `laneFactor[1] - 1/2` | claim 1 (lambda_{21}=1/2) | yes |
| M3 | math | 53 | `laneFactor[2] + 1` | claim 1 (lambda_{22}=-1) | yes |
| M4 | math | 54-55 | same-sign cross terms vanish | claim 1 selection rule | yes |
| M5 | math | 107-108 | wall inertia grouped trace=0 and b-3a=0 | claim 2 | yes |
| M6 | math | 109-110 | wall stiffness grouped trace=0 and b-3a=0 | claim 2 | yes |
| M7 | math | 79-85 | K_1 gate lane vs generic, m=0,1,2 | claim 3 (lane structure) | yes |
| M8 | math | 86-92 | H_even gate lane vs generic, m=0,1,2 | claim 3 (lane structure) | yes |
| M9 | math | 64 | `Det[wallJacobian] - 1/27` | claim 3 (nondegenerate gate) | yes |
| M10 | math | 67-69 | solution count and values trivial | claim 3 (trivial branch only) | yes |
| M11 | math | 111-112 | Xi load grouped trace=0 and b=3a | claim 2 extended | yes |
| M12 | math | 113-114 | prefactor grouped trace=0 and b=3a | claim 2 extended | yes |

All rows: Anchored = yes. No tautological assertions. The mutated-determinant nonzero check (A12, `wall_matrix.det() + 1/27`) supplies an explicit discriminator demonstrating the determinant check has real bite (mutated value evaluates to `2/27 != 0`).

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script derives the lane factors from first principles via

```
Integrate[Sin[theta] * SphericalHarmonicY[2,0,theta,phi] *
                       SphericalHarmonicY[2,m,theta,phi] *
                       SphericalHarmonicY[2,-m,theta,phi],
          {phi, 0, 2*Pi}, {theta, 0, Pi}]
```

(lines 23-44, called via `tripleOnSphere`). This is mechanically independent from SymPy's `gaunt(2,2,2,0,m,-m)` call, which uses tabulated Wigner 3-j coefficients internally. The two engines therefore arrive at `(1, 1/2, -1)` via genuinely distinct paths: closed-form 3-j symbol manipulation (SymPy) vs symbolic sphere integration (Mathematica).

The downstream algebra (grouped projector formulas at lines 94-96; wall Jacobian at lines 61-64; `Solve` on the linear system at lines 65-69) has the same structure as the SymPy side because that algebra IS the claim — there is no second derivation possible of "compute a 2x2 matrix determinant" or "weighted mean coefficients `(1, 2, 2)/5`". This is not `mathematica_transliteration`: the load-bearing physics input (the triple-overlap integral values) is independently derived, and the linear-algebra postprocessing has only one correct form.

## Engine cross-check

Both engines produce residual = 0 on every named check (sympy output `STATUS: PASS`; mathematica output 23 named PASS lines plus `STATUS: PASS`). The lane factors agree exactly: `(1, 1/2, -1)`. The wall-only Jacobian determinant evaluates to `1/27` in both. The linear system yields `{delta K = 0, delta M = 0}` as the sole solution in both. The b=3a tests evaluate cleanly in both. No engine disagreement.

Output freshness: sympy txt mtime is ~1.6 h after the .py mtime; mathematica txt mtime is ~1.6 h after the .wl mtime. `outputs_fresh: true`.

## Verdict justification

Verdict is `clean`. I attempted the following attacks and all failed: (a) verify lane signature `lambda_{22} = -1` from real-harmonic squared overlap — the `(-1)^m` phase convention is correct for m=2 (`(-1)^2 = +1`, and the bare `gaunt(2,2,2,0,2,-2)/gaunt(2,2,2,0,0,0)` integrates to `-1`, consistent with the independent sphere integration in Mathematica); (b) check the same-sign cross-term guard (sympy lines 22-24, math lines 54-55) is not tautological — `gaunt(2,2,2,0,m,m)` for m≠0 vanishes by the m1+m2+m3=0 selection rule, which is a real check on the engine's implementation; (c) compute the wall-only Jacobian by hand from `K1_wall = -dMsym + dKsym/9`, `Hev_wall = (2/3) dMsym - dKsym/27`: `det = (1/9)(2/3) - (-1)(-1/27) = 2/27 - 1/27 = 1/27` — matches; (d) check the mutated determinant assertion at sympy line 110 has discrimination power — `1/27 + 1/27 = 2/27 != 0`, so the assert_nonzero discriminator fires; (e) check that the sympy `Mbar = 0` is non-tautological — weighted mean is `(lambda_{20} + 2 lambda_{21} + 2 lambda_{22})/5 = (1 + 1 + (-2))/5 = 0`, a real cancellation of the `(1, 1/2, -1)` signature, not a definitional zero; (f) confirm both engines derive the lane signature via independent paths (Gaunt symbol vs sphere integration); (g) confirm the m=0 lane factor in the SymPy script genuinely passes through `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)` rather than being hard-coded (the previous v1 short-circuit `if m == 0: return sp.Integer(1)` has been removed — current line 22 only guards the `same_sign` check, not the return value); (h) confirm the previous v1 tautological "wall-only specialization" asserts have been replaced with the lane cross-checks at lines 92-97 — yes. The paper card's three load-bearing claims all have matching, non-tautological, well-anchored assertions in both engines, derived through independent paths where it matters (the Gaunt overlap values). No `paper_misalignment` items; no script-side findings; no user resolution needed.

## Self-test notes

I checked the variable-independence and parity traps. (1) `wall_matrix` takes derivatives of `K1_wall` and `Hev_wall` with respect to `dKsym, dMsym`, and both targets genuinely depend on both variables (Jacobian entries 1/9, -1, -1/27, 2/3 are all literal nonzero rationals), so no zero-derivative trap. (2) The lane factors come from sphere integrals over the bounded `[0, pi] x [0, 2 pi]` domain, so the symmetry/parity-on-unbounded-domain trap does not apply; selection rule m1+m2+m3=0 is the relevant constraint and is correctly handled. (3) Trivial-case substitution: with `dKsym=9, dMsym=1` the K1_wall row gives `-1 + 1 = 0` and the Hev_wall row gives `2/3 - 1/3 = 1/3 != 0`, confirming the two rows of the wall Jacobian are linearly independent and `det = 1/27` is non-vacuous. No additional traps triggered.
