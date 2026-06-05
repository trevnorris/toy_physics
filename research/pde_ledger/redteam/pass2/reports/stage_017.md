---
unit_id: 017
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T03:03:05Z
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
- notes: `(none)` (no `notes/stages/moving_throat_pde_stage017_*.md` exists)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 56 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.txt`

## What the paper claims

Stage 017 transports the parent throat action into the weak-axisymmetric grouped
\(P_2\) lane structure. The card states three deliverables. (1) A \(Y_{20}\) wall
perturbation gives the grouped-lane signature
\(\stagefield{(\lambda_{20},\lambda_{21},\lambda_{22})}=(1,\tfrac12,-1)\)
(eq:stage017-wall-signature, line 18), derived "from the actual triple-overlap
coefficients." (2) Consequently the grouped trace/anomaly variables obey
\(b=3a\) (eq:stage017-b-equals-3a, line 23). (3) A wall-only obstruction: "Pure
wall anisotropy closes the even gates only on the trivial branch
\(\delta K=\delta M=0\)" (line 30), so parent promotion supplies the wall-side
origin of the pattern but does not by itself realize the full weak-axisymmetric
branch. The `\stagefield{Output}` line (lines 34-37) explicitly bounds the exports
to eqs (1)-(2) plus the wall-only gate obstruction. The appendix row 56 summarizes
the stage as the "Weak-axisymmetric parent-action transport law," \StatusExactClosure{}.

## What the script claims to verify

The SymPy script derives the three lane factors from Gaunt (Wigner) overlap
coefficients (`real_y20_square_ratio`, lines 19-25) and asserts they equal
\(1, \tfrac12, -1\) (lines 43-45). It then forms wall inertia/stiffness shifts
\(X_A = \varepsilon\,\lambda_A\,X_1\), computes grouped trace/anomaly variables
(`grouped_trace_anomaly`, lines 28-32), and asserts the trace vanishes and
\(b=3a\) for both inertia and stiffness (lines 57-60). It builds the wall-only
even gates \(K_1\) and \(H_{\mathrm{ev}}\), asserts their Jacobian determinant is
\(1/27\) (line 109), confirms the determinant-mutation guard is nonzero (line 110),
and that the wall-only gate system solves only trivially \(\delta K=\delta M=0\)
(lines 123-124). It additionally checks \(b=3a\) for the \(\Xi_{\mathrm{load}}\)
and prefactor defects (lines 125-128) and cross-checks per-lane gate expressions
against the generic lane formulas (lines 92-97). The Mathematica script
independently re-derives the lane factors by integrating products of
`SphericalHarmonicY` over the sphere (lines 23-49) and mirrors the gate/trace/
determinant/trivial-branch checks (M1-M12).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) lane signature \((1,\tfrac12,-1)\), from triple overlap | py 43-45 (Gaunt); wl M1-M3 (sphere integral) | match |
| (2) grouped \(b=3a\) | py 57-60 (`bM-3aM`, `bK-3aK`); wl M5/M6 | match |
| (3) wall-only obstruction \(\delta K=\delta M=0\) | py 109 (det=1/27), 123-124 (`sol_even` trivial); wl M9/M10 | match |
| (— intermediate) \(\Xi_{\mathrm{load}}\), prefactor \(b=3a\) | py 125-128; wl M11/M12 | extra (supporting; intermediate, terse card legitimately omits) |
| (— intermediate) generic-formula factorization | py 92-97; wl M7/M8 | extra (redundant consistency, see note) |

`paper_alignment: aligned`. The three card-level deliverables each have a faithful,
non-tautological, independently-routed script-side check in both engines. The
"extra" rows are intermediate scaffolding that establish the obstruction and the
propagation of \(b=3a\); they are not separately-stated card deliverables, so their
absence from the terse card is legitimate per the reconciliation guards.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43-45 | `assert_zero(lam2m - target)` | claim 1 (lane signature) | yes |
| A2 | sympy | 24 | same-sign cross term vanishes (raises) | claim 1 (real-harmonic basis) | yes |
| A3 | sympy | 57,59 | `assert_zero(Mbar)`,`assert_zero(Kbar)` | claim 2 (trace) | yes |
| A4 | sympy | 58,60 | `assert_zero(bM-3aM)`,`(bK-3aK)` | claim 2 (b=3a) | yes |
| A5 | sympy | 92-97 | `assert_zero(gate - eps*lam*generic)` | claim 3 (support) | partial (see note) |
| A6 | sympy | 109 | `assert_zero(det - 1/27)` | claim 3 (obstruction) | yes |
| A7 | sympy | 110 | `assert_nonzero(det + 1/27)` | claim 3 (mutation guard) | yes |
| A8 | sympy | 123-124 | `sol_even == [{dKsym:0,dMsym:0}]` | claim 3 (trivial branch) | yes |
| A9 | sympy | 125-128 | `assert_zero(Xibar)`,`(bXi-3aXi)`,`(Pbar)`,`(bP-3aP)` | claim 2 (propagation) | yes |
| A10 | mathematica | 51-55 | `checkZero(laneFactor[m]-target)`, cross terms | claim 1 | yes |
| A11 | mathematica | 64 | `checkZero(Det-1/27)` | claim 3 | yes |
| A12 | mathematica | 67-69 | solution count + trivial \(dK,dM\) | claim 3 | yes |
| A13 | mathematica | 80-91 | M7/M8 gate-vs-generic | claim 3 (support) | partial (see note) |
| A14 | mathematica | 107-114 | M5/M6/M11/M12 trace + b=3a | claim 2 | yes |

Note on A5/A13: the per-lane gate is defined as `D21_2m + D01_2m/9` with
`D21_2m = -eps*lam2m*M1` and `D01_2m = eps*lam2m*K1w` (py 63-64, 47-53), so the
check `gate - eps*lam2m*(-M1+K1w/9) == 0` holds by substitution irrespective of the
value of `lam2m`. It is a redundant factorization/consistency check, not the
load-bearing physics assertion. It is harmless because the substantive content
(lane values, b=3a, determinant, trivial branch) is independently anchored by A1,
A3-A9, A11-A12. It does not change the verdict.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` derives the lane signature independently. The `.py` route is
algebraic Gaunt/Wigner-3j coefficients (`from sympy.physics.wigner import gaunt`,
`gaunt(2,2,2,0,m,-m)/gaunt(2,2,2,0,0,0)`, py 4, 19-25). The `.wl` route is a direct
spherical integral of the product of three `SphericalHarmonicY` functions over
\((\theta,\phi)\) (wl 23-49):
`Integrate[Sin[theta]*SphericalHarmonicY[2,q1,...]*...*SphericalHarmonicY[2,q3,...], {phi,0,2Pi},{theta,0,Pi}]`.
These are genuinely different computational paths to the same coefficients — not a
transliteration. The downstream gate/trace/determinant algebra (kWall, hWall,
Jacobian, Solve) operates on the same shared mathematical objects (the generic lane
formulas), which is expected agreement on the physics, not code-choreography
echoing. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit all-zero residuals and `STATUS: PASS`.
- SymPy output: lane signature `lambda_(20)=1, lambda_(21)=1/2, lambda_(22)=-1`;
  `Mbar=0, a_M=M1*eps/4, b_M=3*M1*eps/4`; `Kbar=0, a_K=K1w*eps/4, b_K=3*K1w*eps/4`;
  `sol_even = [{dKsym:0, dMsym:0}]`; `Xibar=0`, `Pbar=0` with `b=3a`.
- Mathematica output: M1-M3 lane ratios residual 0; M9 determinant residual 0; M10
  trivial dK/dM residual 0; M5/M6/M11/M12 trace + b=3a residual 0.
The determinant (1/27), trivial branch, lane values, and b=3a agree across engines.
No `engine_disagreement`.

## Verdict justification

`clean`. I read the paper card (and confirmed no `notes/stages/stage017_*` file
exists) and appendix row 56 first, built the model (three deliverables: lane
signature \((1,\tfrac12,-1)\); grouped \(b=3a\); wall-only obstruction
\(\delta K=\delta M=0\)), then attacked the scripts. Attacks tried and failed:
(a) tautology hunt — the lane checks are anchored to independent Gaunt/sphere-integral
derivations and the literal targets, so a wrong overlap coefficient would fail them;
the b=3a checks depend on the actual \((1,\tfrac12,-1)\) values (verified by hand:
\(a=\varepsilon X_1/4\), \(b=3\varepsilon X_1/4\), \(b-3a=0\)) and would fail for any
other signature; (b) determinant arithmetic verified by hand: Jacobian
\([[1/9,-1],[-1/27,2/3]]\), \(\det=2/27-1/27=1/27\) — matches assertion; the
trivial-branch claim follows correctly since \(\det\neq0\); (c) symbol domains: `D0,N0`
positive (justified, they are denominators/numerators kept nonzero), others real —
consistent with the physical setup; (d) the only borderline item (A5/A13 generic-formula
factorization) is redundant scaffolding, not a load-bearing assertion, so it does not
weaken any of the three card deliverables. Outputs are fresh (both `.txt` mtimes newer
than their script mtimes). All values reconcile (see below). Paper and script claims
match exactly.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\lambda_{20}=1\) | py 40/43, wl 51; sympy out L8 | stage_017.tex:18 (`\left(1,\frac12,-1\right)`) | MATCH |
| \(\lambda_{21}=1/2\) | py 41/44, wl 52; sympy out L8 | stage_017.tex:18 | MATCH |
| \(\lambda_{22}=-1\) | py 42/45, wl 53; sympy out L8 | stage_017.tex:18 | MATCH |
| grouped law \(b=3a\) (inertia) | py 58, wl 108; sympy out L11-13 (`a_M=M1*eps/4, b_M=3*M1*eps/4`) | stage_017.tex:23 (`b=3a`) | MATCH |
| grouped law \(b=3a\) (stiffness) | py 60, wl 110; sympy out L12-13 (`a_K=K1w*eps/4, b_K=3*K1w*eps/4`) | stage_017.tex:23 | MATCH |
| wall-only trivial branch \(\delta K=\delta M=0\) | py 123-124, wl 65-69; sympy out L21-22 (`[{dKsym:0, dMsym:0}]`) | stage_017.tex:30 (`\(\delta K=\delta M=0\)`) | MATCH |

INTERNAL items (genuine scaffolding, no card/notes counterpart expected, no finding):
- Jacobian determinant \(=1/27\) (py 109, wl 64) — establishes obstruction; intermediate.
- determinant-mutation guard \(+1/27\neq0\) (py 110) — verification scaffolding.
- \(\Xi_{\mathrm{load}}\) grouped defects `Xibar=0, a_Xi=-K1w*eps/(4*D0), b_Xi=-3*K1w*eps/(4*D0)` (py 116/125-126, wl 111-112) — supporting propagation of b=3a.
- prefactor grouped defects `Pbar=0, a_P=-K1w*N0*eps/(4*D0**2), b_P=-3*K1w*N0*eps/(4*D0**2)` (py 122/127-128, wl 113-114) — supporting propagation.
- per-lane explicit gate vectors \(K1_{(20,21,22)}\), \(H_{ev,(20,21,22)}\) (sympy out L16-17) — intermediate displayed forms.
- generic lane formulas `K1_wall=-dM+dK/9`, `H_even,wall=2dM/3-dK/27` (py 90-91, wl 61-62) — intermediate.

All six stage deliverable values reconcile against the card; the intermediate
quantities are legitimate scaffolding a terse `\StatusExactClosure{}` card need not
restate. No `value_mismatch` or `script_missing_paper_claim` raised.

## Self-test notes

Checked: (1) variable independence — the Jacobian `sp.diff(K1_wall, dKsym/dMsym)`
operands genuinely depend on both `dKsym` and `dMsym` (K1_full = dKsym - ... ,
D21_full = -(dMsym + ...)), so no identically-zero-derivative trap; det computed by
hand = 1/27 matches. (2) Trivial-case pre-check — substituting \((1,\tfrac12,-1)\)
into the grouped formulas by hand gives \(a=\varepsilon X_1/4\), \(b=3\varepsilon X_1/4\),
confirming `bM-3aM=0` and `Mbar=0` reduce to literal 0. (3) Tautology trap — flagged
A5/A13 as redundant-by-construction but verified the load-bearing assertions (A1, A3-A9)
are not, so no finding warranted. Conclusion: clean.
