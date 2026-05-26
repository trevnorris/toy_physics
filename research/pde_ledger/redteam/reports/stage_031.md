---
unit_id: 031
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
  paper_appendix: present
---

# Audit unit 031 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_031.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 52 references unit 031)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.txt`

## What the paper claims

Stage 031 proves that on the stable lower-mode branch of the rank-one loaded wall problem, the selected static prefactor `P_{0,-} = beta_0 s_- / lambda_-` is strictly monotone increasing in the loading parameter `alpha_0`, starts at a finite stable-side value, and diverges as `alpha_0 -> alpha_crit^-`. Verbatim, `\stagefield{Output}` reads: "Stage~031 outputs the derivative identities \eqref{eq:app-stage031-ds-dalpha}--\eqref{eq:app-stage031-dP-dalpha}, the refined threshold \eqref{eq:app-stage031-alpha-crit}, and the unique stable crossing theorem." The four distinct deliverables are: (i) exact derivative `ds_-/dalpha = 2 (Delta K_ax)^2 Pi_kappa / R^3 > 0`; (ii) exact derivative `dP_{0,-}/dalpha = beta_0 [(ds_-/dalpha) lambda_- + s_-^2] / lambda_-^2 > 0` whenever `lambda_- > 0` and `beta_0 > 0`; (iii) the refined softening threshold `alpha_crit = AB / (B kappa_0^2 + A kappa_1^2)` with `B = A + Delta K_ax`, together with the endpoint identities `P_{0,-}(0) = beta_0 kappa_0^2 / A` and `P_{0,-} -> +infinity` as `alpha_0 -> alpha_crit^-`; and (iv) the unique stable-side crossing theorem, i.e., every target `P_target > P_{0,-}(0)` is attained at a unique `alpha_* in (0, alpha_crit)`. The notes are slightly more explicit than the .tex on intermediate scaffolding: they spell out the Hellmann-Feynman identity `d lambda_- / d alpha_0 = -s_-` that connects (i) and (ii), and they spell out the determinant identity `lambda_- lambda_+ = AB - alpha_0 (B kappa_0^2 + A kappa_1^2) = T_0 (alpha_crit - alpha_0)` that the script uses in Part V.

## What the script claims to verify

The SymPy script (and its Mathematica mirror) treats the 2x2 loaded-wall eigenproblem analytically by defining `lam_minus`/`lam_plus` as the two algebraic roots of the secular polynomial (parametrized as `(2A + DK - alpha sigma ± R)/2` with `R^2 = (DK + alpha(x0-x1))^2 + 4 alpha^2 x0 x1`). It then uses the Hellmann-Feynman identity as a *definition*, `s_minus := -d lam_minus / d alpha`, and verifies five things: (I) the closed-form identity `ds_-/d alpha = 2 DK^2 x_0 x_1 / R^3` (with `x_0 = kappa_0^2`, `x_1 = kappa_1^2`, so `x_0 x_1 = Pi_kappa`); (II) the closed-form derivative `dP_{0,-}/d alpha = beta_0 (ds_-/d alpha * lam_- + s_-^2) / lam_-^2`, verified both as an abstract quotient/HF identity and as a direct identity on the physical expressions; (III) the initial values `lambda_-(0) = A`, `s_-(0) = x_0`, `P_{0,-}(0) = beta_0 x_0 / A`; (IV) the explicit form of `alpha_crit = A(A+DK)/T_0` with `T_0 = (A+DK) x_0 + A x_1`, together with the determinant factorization `lam_- lam_+ = A(A+DK) - alpha T_0` and the vanishing `lam_-(alpha_crit) = 0`; (V) the rewriting `P_{0,-} = beta_0 s_- lam_+ / (T_0 (alpha_crit - alpha))`, exposing the explicit `1/(alpha_crit - alpha)` pole.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (i) `ds_-/dalpha = 2 (DK_ax)^2 Pi_kappa / R^3` | Part I `expect_zero("ds_-/dalpha exact formula", ds_exact - ds_expected)` (sympy:50; wl:42) | match |
| (i) strict positivity | not asserted; manifest from closed form under declared positivity of `DK, x_0, x_1, R` | partial (formula matches; positivity is read off) |
| (ii) `dP_{0,-}/dalpha = beta_0 [(ds_-/dalpha) lam_- + s_-^2]/lam_-^2` | Part II `expect_zero("dP0_-/dalpha direct identity", dP_direct - dP_physical)` (sympy:73; wl:49) | match |
| (ii) strict positivity | not asserted; follows from positivity of (i) plus `beta_0 > 0`, `s_-^2 > 0`, `lam_-^2 > 0` | partial (formula matches; positivity is read off) |
| (iii) `P_{0,-}(0) = beta_0 kappa_0^2 / A` | Part III `expect_zero("P0_-(0)", ...)` (sympy:82; wl:55) | match |
| (iii) `alpha_crit = AB/(B kappa_0^2 + A kappa_1^2)` | Part IV `alpha_crit = A(A+DK)/T_0` with `T_0 = (A+DK)x_0 + A x_1` (sympy:87; wl:59); `expect_zero("lam_-(alpha_crit)", ...)` (sympy:94; wl:61) | match (algebraically identical: `B = A+DK`, `kappa_0^2 = x_0`, `kappa_1^2 = x_1`) |
| (iii) `P_{0,-} -> +infinity` as `alpha_0 -> alpha_crit^-` | Part V factorization exposes `1/(alpha_crit - alpha)` pole (sympy:117; wl:72); `lambda_+^(crit)` printed but its nonzero-ness not asserted | partial (pole shown; nonvanishing of numerator at `alpha_crit` only inspected, not asserted) |
| (iv) unique crossing theorem | not asserted — this is the continuity+IVT consequence of (i)/(ii)/(iii); ingredients are verified | partial (the theorem is a logical consequence; ingredients are verified) |

Dominant pattern is `match`. The `partial` rows reflect that strict positivity and the IVT-style crossing theorem are not directly assertable as algebraic identities — they are read off the verified closed forms. The script does not test any identity that the paper does not claim. Front-matter `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `expect_zero("ds_-/dalpha exact formula", ds_exact - 2 DK^2 x0 x1 / R^3)` | (i) closed form of `ds_-/dalpha` | yes |
| A2 | sympy | 69 | `expect_zero("generic quotient/HF identity", dP_generic - dP_expected)` | quotient-rule + HF substitution sanity | no — algebraically forced by chain rule with `L' = -S` substituted |
| A3 | sympy | 73 | `expect_zero("dP0_-/dalpha direct identity", dP_direct - dP_physical)` | (ii) closed form of `dP_{0,-}/dalpha` | yes |
| A4 | sympy | 80 | `expect_zero("lambda_-(0)", lam_-(0) - A)` | (iii) endpoint `lam_-(0) = A` | yes |
| A5 | sympy | 81 | `expect_zero("s_-(0)", s_-(0) - x_0)` | (iii) endpoint `s_-(0) = kappa_0^2` | yes |
| A6 | sympy | 82 | `expect_zero("P0_-(0)", P0_-(0) - beta_0 x_0/A)` | (iii) `P_{0,-}(0) = beta_0 kappa_0^2 / A` | yes |
| A7 | sympy | 88 | `expect_zero("det factorization", lam_- lam_+ - (A(A+DK) - alpha T_0))` | (iii) determinant identity → basis for `alpha_crit` | yes |
| A8 | sympy | 94 | `expect_zero("lam_-(alpha_crit)", lam_-(alpha_crit))` after Sqrt rewrite | (iii) `alpha_crit` is the actual softening point | yes |
| A9 | sympy | 98-101 | `expect_zero("threshold radical square identity", T_0^2 R^2(alpha_crit) - root_crit^2)` | internal lemma justifying A8's Sqrt rewrite | yes (it is a real polynomial identity, derived from `R^2 /. alpha->alpha_crit`) |
| A10 | sympy | 110-113 | `expect_zero("lambda_- * lambda_+ - T0(alpha_crit-alpha)", ...)` | identifies simple zero at `alpha_crit` | yes |
| A11 | sympy | 117 | `expect_zero("P0_- factorization", P_{0,-} - beta_0 s_- lam_+ / (T_0(alpha_crit-alpha)))` | (iii) divergence structure | yes |
| B1-B10 | mathematica | 42, 49, 53-55, 60, 61, 65, 70, 72 | mirror of A1, A3, A4-A6, A7, A8, A9, A10, A11 | mirror | mirror; A2 (generic quotient identity) is the one assertion the Mathematica file does *not* mirror |

A2 is the only assertion that is essentially forced by SymPy's chain-rule machinery rather than the physical setup. The companion A3 is the genuine load-bearing check (direct comparison on the physical `s_-` and `lam_-`); A2 is redundant. Notably the Mathematica file omits the A2 analog entirely, going straight to the direct check at line 49.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:26-73`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:32-119`

**What's wrong:**
The Mathematica script is a structural line-by-line port of the SymPy script. It uses the same five-part decomposition with identical banner strings (`PART I` – `PART V`, same English titles including the typo `PREFATOR`), the same intermediate variable choreography, the same closed forms typed in the same order, and the same assertion name strings. Three corresponding pairs:

SymPy Part I (lines 39-50):
```
sigma = x0 + x1
delta_kappa = x0 - x1
KappaProd = x0 * x1
R = sp.sqrt((DK + alpha * delta_kappa) ** 2 + 4 * alpha**2 * KappaProd)
lam_minus = sp.simplify((2 * A + DK - alpha * sigma - R) / 2)
lam_plus = sp.simplify((2 * A + DK - alpha * sigma + R) / 2)
s_minus = sp.simplify(-sp.diff(lam_minus, alpha))
ds_exact = sp.simplify(sp.diff(s_minus, alpha))
ds_expected = sp.simplify(2 * DK**2 * KappaProd / R**3)
expect_zero("ds_-/dalpha exact formula", ds_exact - ds_expected)
```

Mathematica Part I (lines 32-42):
```
sigma = x0 + x1;
deltaKappa = x0 - x1;
kappaProd = x0*x1;
r = Sqrt[(dK + alpha*deltaKappa)^2 + 4*alpha^2*kappaProd];
lamMinus = FullSimplify[(2*a + dK - alpha*sigma - r)/2, ...];
lamPlus  = FullSimplify[(2*a + dK - alpha*sigma + r)/2, ...];
sMinus   = FullSimplify[-D[lamMinus, alpha], ...];
dsExact  = FullSimplify[D[sMinus, alpha], ...];
dsExpected = FullSimplify[2*dK^2*kappaProd/r^3, ...];
expectZero["ds_-/dalpha exact formula", dsExact - dsExpected];
```

SymPy Part IV (lines 86-94) versus Mathematica Part IV (lines 58-65): both construct `T_0 = (A+DK) x_0 + A x_1` first, then back-derive `alpha_crit = A(A+DK)/T_0` from it, then verify the same `lam_- lam_+ - (A(A+DK) - alpha T_0)` factorization, then the same `T_0^2 R^2(alpha_crit) - root_crit^2` Sqrt-rewrite identity.

SymPy Part V (lines 110-117) versus Mathematica Part V (lines 70-72): both verify the same `lam_- lam_+ - T_0(alpha_crit - alpha)` identity and the same `P0_-` factorization.

The two engines are not independently re-deriving the result from physical premises (e.g., from `Eigenvalues[M]` of the loaded 2x2 wall matrix and `Eigenvectors[M]` for the overlap); the `.wl` is transliterating the `.py`'s algebraic choreography. The only substantive divergence is that the Mathematica script *omits* the redundant A2 ("generic quotient/HF identity") check and goes straight to the direct identity — but that is a single deletion, not an independent derivation.

**Why this matters:**
The second-engine policy exists to catch SymPy-specific algebra bugs and to provide independent confirmation. A transliteration cannot catch a shared modeling error: if the parametrization `(2A + DK - alpha sigma ± R)/2` were off by a sign convention, or if `s_- := -d lam_-/d alpha` were the wrong sign for the lower branch, both engines would silently echo the mistake. The practical risk for this unit is low because the algebra is standard quadratic-root manipulation, but the policy is uniformly violated.

**Required change:**
Restructure the Mathematica script so it derives the same five results via an *independent* route from the secular polynomial. Suggested independent route:
1. Build the loaded 2x2 wall matrix `M = {{a, 0}, {0, a + dK}} - alpha * Outer[Times, {Sqrt[x0], Sqrt[x1]}, {Sqrt[x0], Sqrt[x1]}]` (or its equivalent symmetric form). Use `Eigenvalues[M]` to obtain the two eigenvalues; assert that one equals SymPy's `lamMinus` parametrization, rather than constructing the parametrization directly.
2. Use `Eigenvectors[M]` to extract the lower-eigenvector `e_-`, compute `sMinusOverlap = (v . e_-)^2 / (e_- . e_-)` where `v = {Sqrt[x0], Sqrt[x1]}`, and assert the Hellmann-Feynman identity `D[lamMinus, alpha] + sMinusOverlap == 0` as a verified theorem. This converts the script's load-bearing *definition* `s_- := -d lam_-/d alpha` into something the Mathematica side checks independently.
3. Derive `alpha_crit` from `Det[M] == 0` rather than building `T_0` first and back-deriving.
4. Re-derive Parts III and V from the above without copying intermediate variable names from the SymPy script.

**Verification:**
The fixed `.wl` should still emit `PASS` for every existing check, but its variable names, top-level decomposition, and intermediate quantities should differ from the `.py`'s. In particular it should contain `Eigenvalues[...]` and `Eigenvectors[...]` calls that the SymPy script does not. A new HF identity check (item 2 above) should appear as an additional assertion in the Mathematica transcript.

---

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:60-69`

**What's wrong:**
The "generic quotient/HF identity" check at sympy:60-69 is algebraically forced by `sympy.diff` and the substitution map. It declares abstract symbols `Lsym, Ssym, DSsym`, computes `d/dalpha(beta_0 S(alpha)/L(alpha))` via SymPy's symbolic differentiation, substitutes `S -> Ssym`, `L -> Lsym`, `dS/dalpha -> DSsym`, `dL/dalpha -> -Ssym`, and asserts the result equals `beta_0 (DSsym Lsym + Ssym^2)/Lsym^2`. By the quotient rule, `d/dalpha(S/L) = (S' L - S L')/L^2 = (DS·L - S·(-S))/L^2 = (DS·L + S^2)/L^2`; the check is true for any `S` and `L` whatsoever, regardless of physics. It exercises SymPy's substitution machinery, not the moving-throat physics.

Quoted from the script (lines 60-69):
```
Lsym, Ssym, DSsym = sp.symbols("Lsym Ssym DSsym", nonzero=True, real=True)
dP_generic = sp.diff(beta0 * sp.Function("S")(alpha) / sp.Function("L")(alpha), alpha)
dP_generic = dP_generic.subs({
    sp.diff(sp.Function("S")(alpha), alpha): DSsym,
    sp.Function("S")(alpha): Ssym,
    sp.diff(sp.Function("L")(alpha), alpha): -Ssym,
    sp.Function("L")(alpha): Lsym,
})
dP_expected = sp.simplify(beta0 * (DSsym * Lsym + Ssym**2) / Lsym**2)
expect_zero("generic quotient/HF identity", sp.simplify(dP_generic - dP_expected))
```

The companion direct identity at line 73 (`dP_direct = sp.diff(P0_sel, alpha)` against `beta_0(ds_expected lam_- + s_-^2)/lam_-^2`) is the load-bearing physical check and is unaffected by this finding. The Mathematica script (line 47-49) skips the generic identity entirely and uses only the direct check — confirming that the generic check is not load-bearing.

**Why this matters:**
The check cannot fail except under a SymPy bug; it does not exercise anything about the moving-throat physics or the `lam_-`/`s_-` definitions. It is harmless (the next check is the genuine one), but its presence in the transcript creates the false impression that an additional physical identity was verified.

**Required change:**
Either (a) delete lines 60-69 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py` (the next assertion at lines 71-73 is the load-bearing check and is unaffected); or (b) re-purpose the block to verify the Hellmann-Feynman identity itself by constructing `s_-` independently as the eigenvector overlap `(v . e_-)^2 / (e_- . e_-)` from the 2x2 matrix and asserting `sp.simplify(sp.diff(lam_minus, alpha) + s_minus_via_overlap) == 0`. Option (b) would convert the script's load-bearing *definition* into a verified theorem; option (a) is purely cosmetic. Prefer (a) unless implementing (b) is straightforward.

**Verification:**
After option (a), the saved sympy output's "PART II" section should no longer contain the `generic quotient/HF identity = 0` line; the `dP0_-/dalpha direct identity = 0` line must still appear. After option (b), a new `HF identity` assertion appears and passes. Script must continue to exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is a structural transliteration of the SymPy script (see F1). It shares the same `(2A + DK - alpha sigma ± R)/2` parameterization of the eigenvalues, the same Hellmann-Feynman shortcut `s_- := -D[lam_-, alpha]`, the same banner labels (including the `PREFATOR` typo), the same five-part decomposition, and the same intermediate constructions (`T_0`, `radCritDerived`, `p0Factored`). The lone substantive divergence is that the Mathematica script omits the SymPy script's "generic quotient/HF identity" check (line 47-49 uses the direct check directly), which is an *improvement* but not an independent derivation. Three corresponding pairs were quoted under F1.

## Engine cross-check

Both engines emit `0` on every assertion and `PASS`/`EXIT 0` overall. Side-by-side of the printed final forms:

- `ds_-/dalpha`:
  - SymPy txt lines 5-13: `2 DK^2 x_0 x_1 / (4 alpha^2 x_0 x_1 + (DK + alpha(x_0 - x_1))^2)^(3/2)`
  - Mathematica txt line 7: `(2*dK^2*x0*x1)/((dK + alpha*(x0 - x1))^2 + 4*alpha^2*x0*x1)^(3/2)`
  - Algebraically identical.
- `alpha_crit`:
  - SymPy txt lines 73-76: `A(A + DK) / (A x_1 + x_0(A + DK))`
  - Mathematica txt line 35: `(a(a + dK)) / (dK x_0 + a(x_0 + x_1))`
  - Equal: `A x_1 + (A+DK) x_0 = A x_1 + A x_0 + DK x_0 = A(x_0 + x_1) + DK x_0`.
- `lambda_+^(crit)`:
  - Mathematica txt line 36 (cleaned): `((a + dK)^2 x_0 + a^2 x_1) / (dK x_0 + a(x_0 + x_1))`
  - SymPy txt lines 78-89 (not yet collapsed because `sp.simplify` does not denest `sqrt(perfect square)` here without a hint, but the symbolic form contains `sqrt((A(A+DK)(x_0-x_1) + DK T_0)^2) = (A+DK)^2 x_0 + A^2 x_1` which is the radicand identity that the script's A9 check verifies).
  - Same after simplification.
- `P_{0,-}` factored: both engines emit the same `1/(alpha_crit - alpha)` denominator (sympy txt lines 95-119; math txt line 45).

No `engine_disagreement` finding. Outputs are fresh (script mtimes 17:18 and 17:20; output mtimes 17:22 for both — outputs newer than scripts).

## Verdict justification

The script's bottom-line assertions match the paper card's four deliverables: closed form of `ds_-/dalpha`, closed form of `dP_{0,-}/dalpha`, endpoint values, `alpha_crit` (algebraically identical to `AB / (B kappa_0^2 + A kappa_1^2)`), and the pole-structure factorization that supports the divergence claim. Strict positivity of the two derivatives and the unique-crossing theorem are not symbolically asserted (positivity is read off the verified closed forms under declared positive symbols; the crossing theorem is a continuity+IVT consequence). Both engines pass all assertions and agree on all printed forms.

Attacks tried that did not break the audit: (a) verified that `s_-(0) = kappa_0^2` is not implicitly forced by the way `s_- := -d lam_-/d alpha` is defined — the algebra at `alpha=0` gives `-d lam_-/d alpha = x_0` only because `R(0) = DK` and the inner derivative of `R^2` evaluates to `2 DK (x_0 - x_1)`, so the check is a genuine algebraic test of the eigenvalue parametrization; (b) verified that the script's `alpha_crit` matches the paper's `AB/(B kappa_0^2 + A kappa_1^2)` under the substitution `B = A + DK`, `kappa_0^2 = x_0`, `kappa_1^2 = x_1`; (c) verified that `expect_zero("lam_-(alpha_crit)", ...)` is not faked via aggressive simplification — the script's manual Sqrt rewrite (lines 90-93) replaces `sqrt(radicand)` with `root_crit = A^2 x_1 + (A+DK)^2 x_0` only after asserting (line 98-101) that `T_0^2 R(alpha_crit)^2 = root_crit^2`, which is an exact polynomial identity verified by `sp.expand`; (d) verified that the `s_minus^2` term in the `dP_{0,-}/dalpha` formula (from the quotient rule with `dL/dalpha = -S` substituted) is correctly placed and not absorbed into an aggressive `simplify`.

The two findings are about *form*, not *substance*: F1 is the systemic transliteration issue, F2 is one redundant scaffold assertion that exists in the SymPy file only. Neither changes any quoted result; downstream units consuming `alpha_crit`, `ds_-/dalpha`, or `P_{0,-}` see identical values. No `UNFIXABLE`, no `CRITICAL_DOWNSTREAM`, no `paper_misalignment`. Verdict `findings`.

## Self-test notes

I checked: (1) `s_minus` is defined as `-d lam_-/d alpha` and not as the eigenvector overlap, so the HF identity is built-in rather than independently verified — flagged as the secondary option in F2's "option (b)"; (2) the `expect_zero("lam_-(alpha_crit)", ...)` Sqrt rewrite is itself anchored by the immediately following `threshold radical square identity` check (lines 97-101) which is a real polynomial identity, so A8 is not load-bearing on a hidden assumption; (3) no `assert_zero` here depends on a derivative with respect to a variable the expression doesn't contain — all `sp.diff(..., alpha)` calls act on expressions that genuinely depend on `alpha`; (4) `paper_misalignment` was carefully checked under `B = A + Delta K_ax`, `kappa_0^2 = x_0`, `kappa_1^2 = x_1`, `Pi_kappa = x_0 x_1` — every constant matches; (5) no integrals on unbounded domains are present, so parity traps do not apply.
