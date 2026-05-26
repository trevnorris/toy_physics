---
unit_id: 024
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md"]
  paper_appendix: present
---

# Audit unit 024 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_024.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 38: "Overlap integrals, O(3) isotropy, and first axisymmetric splitting" — \StatusExactClosure / \StatusReduced)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt`

## What the paper claims

The paper card's `\stagefield{Output}` is verbatim: "Stage 024 outputs the real-STF normalization (Y-orthonormal), the angular source-map identity (S_A^port = S_A), the isotropic overlap formulas (B_n, Z_n, N_n, D), the grouped isotropy collapse (D_{A,n} and N_{A,n} lane-independent), the weak-axisymmetric signature lambda = (1, 1/2, -1), the line b = 3a, and the normalization transport coefficient P_1 = (N_1 D_0 - N_0 D_1)/D_0^2." Concretely the card boxes seven results: (i) sphere orthonormality `int Y_A Y_B dOmega = delta_AB` from the exact fourth-moment identity; (ii) the angular projection `S_A^port = S_A`, hence `mhat_ang = 1`; (iii) the isotropic-branch closed forms for `B_n`, `Z_n`, `N_n` and the conservative operator `D(omega) = D_0 + D_2 omega^2 + D_4 omega^4`; (iv) the grouped isotropic collapse `D_{20,n} = D_{21,n} = D_{22,n}` and `N_{20,n} = N_{21,n} = N_{22,n}`; (v) the triple-overlap matrix `M^(20) = (sqrt(5)/(7 sqrt(pi))) diag(1, 1/2, 1/2, -1, -1)` from the sixth moment of the unit sphere; (vi) grouped signature `(lambda_20, lambda_21, lambda_22) = (1, 1/2, -1)` with the `b_x = 3 a_x` line via the `(1,2,2)` weighted trace/anomaly map; (vii) first-order transport `P_1 = (N_1 D_0 - N_0 D_1)/D_0^2` and the defect law `b_P = 3 a_P`. The notes confirm the same content and explicitly attribute the `(1,2,2)` grouped weighting to Stage 022. The card uses `\StatusExactClosure` for the angular identities and `\StatusReduced` for the separated radial/axial ansatz.

## What the script claims to verify

The SymPy script's docstring and section headers cover all seven of the paper-side deliverables. Concretely the assertions test: (I.1) the 5x5 Gram matrix of normalized real-STF harmonics is the identity (via the analytic delta-pairing form of `int n_i n_j n_k n_l dOmega`); (I.2) `gram * svec = svec`, identified with `mhat_ang = 1`; (II.1) substituting equal-lane data into the `(1,2,2)`-weighted trace/anomaly maps `xbar`, `a_x`, `b_x` gives `(x_0, 0, 0)`; (II.2) two scalar witnesses pinning the linear weights to `(2,-1,-1)/10` and `(0,1,-1)/2`; (III.1) a two-mode BdG witness whose Taylor coefficients of `sum_alpha C_alpha^2/(varpi_alpha^2 - omega^2)` reproduce `B_0`, `B_2`, `B_4`; (III.2) a one-pair Maxwell/mixed witness whose `Z(omega) = (Q - H omega^2)/(Delta - S omega^2 + omega^4)` and `N(omega) = (P - GW omega^2)^2/(...)^2` Taylor series reproduce the paper's `Z_0, Z_2, Z_4, N_0, N_2, N_4` and `D_0, D_2, D_4`; (IV.1) the 5x5 triple-overlap matrix `M_{AB}^{(20)} = int Y_A Y_{20} Y_B dOmega` equals `(sqrt(5)/(7 sqrt(pi))) diag(1, 1/2, 1/2, -1, -1)`; (IV.2) the `(1, 1/2, -1)` first-order ansatz yields `a_x = eps x_1/4`, `b_x = 3 eps x_1/4`, `b_x = 3 a_x`; (V) `lane_ratio(lambda) = (N_0 + eps lambda N_1)/(D_0 + eps lambda D_1)` expanded to first order in `eps` reproduces `P_0 + eps lambda P_1` with `P_1 = (N_1 D_0 - N_0 D_1)/D_0^2`, and the resulting `(a_P, b_P)` satisfy `b_P = 3 a_P`. The Mathematica `.wl` script mirrors all of these, with sphere-integral methods (`Integrate[n[i] n[j] ... Sin[theta], {theta, 0, Pi}, {phi, 0, 2 Pi}]`) for I.1 and IV.1 (genuinely independent re-derivation), and otherwise the same algebraic constructions and assertion forms as the SymPy script.

## Paper <-> script cross-check

| # | Paper-side deliverable | Script-side check | Status |
|---|---|---|---|
| P1 | `int Y_A Y_B dOmega = delta_AB` (boxed, Eq. \eqref{eq:app-stage024-Y-orthonormal}) | SymPy `Gram - I5` via `quad_overlap`; Mathematica `Gram - I5` via spherical `Integrate` | match |
| P2 | `S_A^port = S_A`, hence `mhat_ang = 1` (boxed, Eq. \eqref{eq:app-stage024-source-map}) | `gram * svec - svec` | match (trivial given P1) |
| P3 | `B_0, B_2, B_4` closed forms (boxed, Eq. \eqref{eq:app-stage024-B-moments}) | Two-mode Taylor of `sum C^2/(varpi^2-omega^2)`; checks each coefficient against `sum C^2/varpi^(2k+2)` | match (algebraic) |
| P4 | `Z_0, Z_2, Z_4` and `N_0, N_2, N_4` closed forms (Eq. \eqref{eq:app-stage024-ZN-moments}) | One-pair Taylor of `(Q - H omega^2)/(Delta - S omega^2 + omega^4)` and squared variant for N; coefficients matched against paper expressions | match (algebraic; see F2) |
| P5 | `D(omega) = D_0 + D_2 omega^2 + D_4 omega^4` with `D_0 = K - B_0 - Z_0`, `D_2 = -(M + B_2 + Z_2)`, `D_4 = -(B_4 + Z_4)` (boxed, Eq. \eqref{eq:app-stage024-D-isotropic}) | `D_0, D_2, D_4` Taylor coefficients of `K - M omega^2 - Bresp - Zresp` | match (algebraic) |
| P6 | Grouped isotropy collapse `D_{20,n} = D_{21,n} = D_{22,n}` (boxed, Eq. \eqref{eq:app-stage024-isotropic-collapse}) | Implicit from lane-independence of `C_alpha`, `G_U`, `G_W`, `R_r`; not separately asserted as such | partial (see F4) |
| P7 | `M^(20) = (sqrt(5)/(7 sqrt(pi))) diag(1, 1/2, 1/2, -1, -1)` (boxed, Eq. \eqref{eq:app-stage024-axisym-matrix}) | SymPy `M - M_target` via `triple_overlap`; Mathematica via spherical `Integrate` | match |
| P8 | `(lambda_20, lambda_21, lambda_22) = (1, 1/2, -1)` (boxed, Eq. \eqref{eq:app-stage024-lambda-signature}) | Lane weights pulled directly from M^(20) diagonal; signature used in IV.2 ansatz | match |
| P9 | `b_x = 3 a_x` (boxed, Eq. \eqref{eq:app-stage024-b-3a}) | `bx - 3*ax` zero check on first-order ansatz | match |
| P10 | `P_1 = (N_1 D_0 - N_0 D_1)/D_0^2` (boxed, Eq. \eqref{eq:app-stage024-PA-transport}) | Series expansion of `lane_ratio` to first order, plus three lane checks | match (algebraic) |
| P11 | `b_P = 3 a_P` (boxed, Eq. \eqref{eq:app-stage024-P-defects}) | `bP - 3*aP` zero check on first-order P_A | match |

Dominant pattern: aligned. All paper-side deliverables are tested; some are tested only by algebraic self-consistency rather than independent derivation (see F2). Paper alignment set to `aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 164 | `expect_zero("Gram - I5", gram - sp.eye(5))` | P1 | yes (sphere integral) |
| A2 | sympy | 171 | `expect_zero("projected - svec", projected - svec)` | P2 | partial (trivially follows from A1) |
| A3 | sympy | 193 | `expect_zero("xbar - x0", xbar.subs(equal_lane) - x0)` | P6 (carry to P9) | no (linear-combo definition substituted with equal values) |
| A4 | sympy | 194-195 | `expect_zero("a_x on equal lanes", ax.subs(equal_lane))`, same for bx | P6 | no |
| A5 | sympy | 200-203 | Four scalar witnesses on (a_x, b_x) maps | (weighting of P9 map) | no (arithmetic identity by construction) |
| A6 | sympy | 226-228 | `expect_zero("C_alpha - lambda_B,alpha I_etaalpha", ...)` | P3 helper | no (defines `C1 = lamB1*Ieta1`, then checks equal to itself — fully tautological) |
| A7 | sympy | 237-239 | B_0, B_2, B_4 Taylor coefficient checks (two-mode witness) | P3 | partial (algebraic Taylor of a chosen rational; cannot fail unless the paper's claimed simplification is wrong) |
| A8 | sympy | 275-283 | Z_0, Z_2, Z_4, N_0, N_2, N_4 Taylor coefficient checks (one-pair witness) | P4 | partial (algebraic Taylor of a chosen rational; see F2) |
| A9 | sympy | 284-286 | D_0, D_2, D_4 Taylor coefficient checks | P5 | partial (algebraic combination of A7/A8) |
| A10 | sympy | 305 | `expect_zero("M - M_target", M - M_target)` | P7 | yes (sphere integral; substantive) |
| A11 | sympy | 323-326 | xbar/a_x/b_x/b_x-3a_x checks on lambda=(1,1/2,-1) ansatz | P8, P9 | partial (algebraic; weights chosen from A10 result) |
| A12 | sympy | 359-361 | P_{20}/P_{21}/P_{22} lane expansion checks | P10 | partial (algebraic Taylor of constructed rational) |
| A13 | sympy | 368-371 | abar/a_P/b_P/b_P-3a_P checks on first-order P_A | P11 | partial (algebraic, follows from A12) |
| A14 | math | 73 | `expectZero["Gram - I5", gram - IdentityMatrix[5]]` (via spherical Integrate) | P1 | yes (independent sphere integral) |
| A15 | math | 78 | `expectZero["projected - svec", gram.svec - svec]` | P2 | partial (follows from A14) |
| A16 | math | 86-92 | xbar/ax/bx equal-lane and witness checks | P6, weights | no (same arithmetic identities as A3-A5) |
| A17 | math | 102 | `expectZero["C_alpha - lambda_B,alpha I_etaalpha", {c1, c2} - {lambdaB1*iEta1, lambdaB2*iEta2}]` | P3 helper | no (fully tautological) |
| A18 | math | 109-111 | B_0, B_2, B_4 Taylor coefficient checks | P3 | partial (same algebraic check as A7) |
| A19 | math | 142-150 | Z/N/D Taylor coefficient checks | P4, P5 | partial (same algebraic check as A8/A9; see F1) |
| A20 | math | 158 | `expectZero["M - M_target", m20 - mtarget]` (via spherical Integrate) | P7 | yes (independent sphere integral) |
| A21 | math | 169-172 | xbar/a_x/b_x/b_x-3a_x checks on lambda=(1,1/2,-1) ansatz | P8, P9 | partial (same as A11) |
| A22 | math | 183-192 | P_{20}/P_{21}/P_{22}/Pbar/a_P/b_P/b_P-3a_P checks | P10, P11 | partial (same as A12/A13) |

The "Exercises which paper claim?" column shows that every script-side check ties to a paper-side deliverable. No orphan assertions (no `paper_missing_script_claim`). However, only assertions A1, A10, A14, A20 are genuinely substantive (non-algebraic) checks of paper content — and A14/A20 are independent re-derivations of A1/A10 via a different integration method, so the cross-engine agreement carries real weight only for the angular sphere integrals (I.1 and IV.1).

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:80-92` (Section II)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:94-150` (Section III)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:160-192` (Sections IV.2 and V)

**What's wrong:**
Sections II, III, IV.2, and V of the Mathematica script are structurally line-by-line ports of the corresponding SymPy sections rather than independent re-derivations. Compare for example:

SymPy III.2 (lines 247-253):
```
GU = sp.simplify(lamU * IetaU)
GW = sp.simplify(lamW * IetaW)
Rr = sp.simplify(lamR * Iuw)
Delta = sp.simplify(OmegaU**2 * OmegaW**2 - Rr**2)
S = sp.simplify(OmegaU**2 + OmegaW**2)
Q = sp.simplify(GU**2 * OmegaW**2 + 2*GU*GW*Rr + GW**2*OmegaU**2)
H = sp.simplify(GU**2 + GW**2)
P = sp.simplify(OmegaU**2 * GW + Rr * GU)
```

Mathematica III (lines 117-123) does the same in the same order, with parallel variable names (`gU, gW, rPair`, etc.) and identical algebraic groupings. Both then build `Zresp = Series[(Q - H*omega^2)/(Delta - S*omega^2 + omega^4), omega, 0, 4/6]` and extract identical coefficients, comparing to identical hand-written closed forms (sympy lines 275-286 vs mathematica lines 142-150). Section V's `lane_ratio` function (sympy 347-349) is mirrored verbatim in Mathematica `laneRatio[lam_]` (line 177), with the same eps-expansion and the same three lane checks against `(P_0 + eps P_1)`, `(P_0 + eps P_1/2)`, `(P_0 - eps P_1)`.

Sections I.1 and IV.1 are NOT transliterated — Mathematica genuinely uses spherical `Integrate[..., Sin[theta], ...]` while SymPy uses the analytic delta-pairing form of the fourth/sixth moment. Those two sections constitute the only genuine independent cross-check. For the rest of the script, the two engines run the same algebra in two languages.

**Why this matters:**
The second-engine policy requires that the Mathematica side independently derive the result so that an algebra error in the SymPy side cannot silently propagate to a matching algebra error in Mathematica. For Sections II, III, IV.2, and V, both engines would arrive at the same answer for the same wrong reason. This is a checkpoint stage (is_checkpoint: true), so the bar for independent verification is higher.

**Required change:**
For Section III (Z, N, D coefficient closed forms), have the Mathematica engine derive `Z_0, Z_2, Z_4, N_0, N_2, N_4` via a different path. A natural independent path: construct the 2x2 coupled-pair system as a matrix
```
{{Omega_U^2 - omega^2, R}, {R, Omega_W^2 - omega^2}}
```
solve for `(U, W)` driven by `(G_U eta, G_W eta)` via `Inverse[]` or `LinearSolve[]`, contract back onto eta to obtain the per-pair Z response (and similarly the per-pair N), then Taylor expand to order omega^4 and assert the same coefficients. That gives `Z_n`/`N_n` from matrix inversion rather than from a hand-chosen rational form. SymPy continues to use the rational-form path.

For Section V (first-order P_A transport), the Mathematica side can take an alternative path by computing `P_A = N_A/D_A` directly with the closed-form `P_0 = N_0/D_0`, `P_1 = (N_1 D_0 - N_0 D_1)/D_0^2`, and then expanding `(N_0 + eps lambda N_1)/(D_0 + eps lambda D_1) - (P_0 + eps lambda P_1)` to first order via `Normal[Series[..., {eps, 0, 1}]]` — actually this is what it does, but it should at minimum NOT mirror SymPy's `lane_ratio` helper structure verbatim.

For Section II (grouped weighting), there is no obvious independent algebraic path because the `(1,2,2)` weighting is a definitional convention carried in from Stage 022. The minimum mitigation is to derive the weighting from the Stage 022 grouped metric inner product (which Stage 022's audit covers).

**Verification:**
After fixes, both engines should land on the same Z/N/D coefficient formulas but the Mathematica derivation chain should start from a 2x2 matrix inverse (or other physically anchored object) rather than from the same rational `(Q - H omega^2)/(Delta - S omega^2 + omega^4)` already used in SymPy. The new derivation should be visible as a `LinearSolve` or `Inverse` call in the Mathematica script, replacing the current direct `Series[(Q - H*omega^2)/...]` construction.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:230-239` (Section III.1, B_n)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:255-286` (Section III.2, Z_n, N_n, D_n)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:338-371` (Section V, P_A transport)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:104-150` (mirror of III)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:174-192` (mirror of V)

**What's wrong:**
The Section III and Section V assertions check that the paper's stated closed-form Taylor coefficients are algebraically consistent with chosen rational generating functions, but they do not check that the rational generating function arises from the right physics. Specifically:

For B_n: the script defines `Bresp = C1^2/(varpi1^2 - omega^2) + C2^2/(varpi2^2 - omega^2)` (line 230), Taylor-expands in `omega^2`, and verifies `B_0 = C1^2/varpi1^2 + C2^2/varpi2^2`. But this is `Bresp(0) = sum_alpha C_alpha^2/varpi_alpha^2` by direct substitution; the assertion cannot fail no matter what the underlying BdG physics demands.

For Z_n, N_n: the script defines `Zresp = (Q - H omega^2)/(Delta - S omega^2 + omega^4)` and `Nresp = (P - GW omega^2)^2/(Delta - S omega^2 + omega^4)^2`, Taylor-expands, and confirms coefficients against the paper's `Z_n = Q/Delta`, `(QS-HDelta)/Delta^2`, etc. The chosen rationals are precisely the closed forms whose series the paper's `Z_n, N_n` are stated to be — so the assertions are essentially "Taylor coefficient of f equals the paper's Taylor coefficient of f" which is true by symbolic algebra alone. They never test that the chosen rationals are the actual coupled Maxwell/mixed loop response.

For P_A: `lane_ratio(lambda) = (N_0 + eps lambda N_1)/(D_0 + eps lambda D_1)` Taylor-expanded to first order in `eps` and compared to `P_0 + eps lambda P_1` is a textbook quotient derivative identity that holds by construction.

**Why this matters:**
At checkpoint quality, these sections need to either (a) derive the rational generating functions from a physically anchored upstream construction (e.g., matrix inversion of the 2x2 conservative Maxwell/mixed sub-block), or (b) be explicitly labeled as "consistency of paper's algebraic simplifications" rather than "verification of physics." Currently the script transcripts read as full physics verification, but they are algebraic identities that cannot detect an error in the upstream construction. A typo of, say, `Q_r = G_U^2 Omega_W^2 + G_U G_W R_r + G_W^2 Omega_U^2` (missing the factor of 2) would NOT be caught because the same definition is used on both sides of the assertion.

This bites doubly because the rational forms ARE the load-bearing claims of Stage 024 (carried forward into stage 025+ as `B_n`, `Z_n`, `N_n`, `D_n`). If they are wrong, Stage 025's normalization is wrong. The script as written confirms only that the paper's stated closed forms can be put in a self-consistent rational presentation — it does not catch an upstream typo in the physical premise.

**Required change:**
Add a substantive anchor for the per-pair Maxwell/mixed response. Concretely, add a SymPy block (and have Mathematica derive independently per F1) that:

1. Constructs the per-pair 2x2 conservative matrix
   ```
   M_pair = sp.Matrix([[OmegaU**2 - omega**2, Rr], [Rr, OmegaW**2 - omega**2]])
   ```
   (these are the off-diagonal R-coupled `(U, W)` modes from the Stage-023 conservative sector — confirm sign of R_r against Stage 023 before committing).
2. Inverts and contracts against the eta-coupling vector `(GU, GW)`:
   ```
   coupling = sp.Matrix([GU, GW])
   Z_full = (coupling.T * M_pair.inv() * coupling)[0, 0]
   ```
3. Asserts that `Z_full == (Q - H omega^2)/(Delta - S omega^2 + omega^4)` after `sp.simplify`. This is the missing physics-to-rational anchor.
4. Then the existing Taylor-coefficient assertions (Z_0 = Q/Delta etc.) become real verifications of the paper's simplification of the matrix-inverse response.

For B_n, the existing two-mode witness is essentially `(1/(a-x)) Taylor` arithmetic and is not at risk; leave as is, but add a one-line comment that the check is the geometric-series identity, not an anchored physics check.

For Section V, no new physics is involved (it is pure first-order quotient calculus), but the existing checks should be tagged as algebraic. Alternatively, add a finite-difference numerical witness: pick concrete numeric values of N_0, N_1, D_0, D_1 and `eps = 1e-3`, evaluate `P_A` directly versus `P_0 + eps lambda P_1`, and assert the difference is `O(eps^2)`. That converts the symbolic identity into a falsifiable numerical claim.

**Verification:**
After fix, Section III gains an additional `expect_zero("Z_full matches rational form", Z_full - (Q - H*omega**2)/(Delta - S*omega**2 + omega**4))` (or equivalent) that exercises the 2x2-inverse derivation. If that assertion ever fails (e.g., due to a sign flip on R_r), it catches an upstream-physics error that the current script cannot. The Mathematica side gets the same matrix-inverse anchor, derived in Mathematica without copying SymPy's algebra (closes F1 for Section III simultaneously).

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:191-203` (Section II.1 and II.2)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:225-228` (Section III preamble C_alpha check)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:81-92` (mirror of II)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:100-102` (mirror of C_alpha check)

**What's wrong:**
The C_alpha assertion (sympy lines 223-228):
```
C1 = sp.simplify(lamB1 * Ieta1)
C2 = sp.simplify(lamB2 * Ieta2)
expect_zero(
    "C_alpha - lambda_B,alpha I_etaalpha",
    sp.Matrix([C1, C2]) - sp.Matrix([lamB1 * Ieta1, lamB2 * Ieta2]),
)
```
defines `C1` as exactly `lamB1 * Ieta1`, then asserts `C1 == lamB1 * Ieta1`. This is a definition-equals-itself check; it cannot fail. Same in the Mathematica version (line 102).

The Section II equal-lane checks (sympy lines 192-195):
```
xbar = sp.simplify((x20 + 2*x21 + 2*x22)/5)
...
equal_lane = {x20: x0, x21: x0, x22: x0}
expect_zero("xbar - x0", xbar.subs(equal_lane) - x0)
expect_zero("a_x on equal lanes", ax.subs(equal_lane))
expect_zero("b_x on equal lanes", bx.subs(equal_lane))
```
substitute `x_0` into linear maps that have row-sums `5/5=1`, `0/10=0`, `0/2=0` by construction. Cannot fail.

The Section II.2 witness checks (sympy lines 200-203) substitute single-coordinate offsets and check specific rationals (1/5, -1/10, 1/2). These are pure arithmetic identities of the chosen linear weights. They are slightly more useful in that they pin down the (1, 2, 2)/5 weighting, but cannot detect errors in the physics — only in the script-side definition of the weights, which then matches the assertion by construction.

**Why this matters:**
At checkpoint quality, definition-equals-itself assertions clutter the transcript with green PASSes that do not provide verification value. The transcript can mislead a reader into thinking that 7 things were checked when really 2 things were checked (the definitions and the linear-combination arithmetic). Lower severity than F1/F2 because the checks at least document the conventions and would catch a careless rename, but they should be either removed or relabeled as "definitional consistency."

**Required change:**
Either (a) remove the C_alpha self-equality assertion entirely, or (b) replace it with a check that adds verification value — e.g., assert that on the equal-coupling limit `lamB1 = lamB2`, `Ieta1 = Ieta2`, the result for `C_alpha` is the common scalar `C`. Similarly for the equal-lane II.1 assertions, prefer keeping only the II.2 witnesses (which at least pin specific arithmetic), and tag II.2 in the transcript banner as "definitional pinning of (1,2,2)/5 weighting carried from Stage 022."

Apply the same changes to the Mathematica script. Do NOT change the Section II semantics; just stop asserting tautologies and replace with witnesses that pin a specific weighting or convention.

**Verification:**
After fix, the C_alpha assertion is removed or rewritten. Section II.1 equal-lane assertions are removed or rewritten to test something falsifiable (e.g., a non-equal-lane configuration whose `xbar` value can be predicted only from the (1,2,2)/5 weighting). Section II.2 retains its arithmetic witnesses but is explicitly labeled as a definitional pinning, not a physics check.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py` (the script as a whole)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl` (the script as a whole)

**What's wrong:**
The paper's grouped isotropy collapse (Eq. \eqref{eq:app-stage024-isotropic-collapse}, boxed) reads
```
D_{20,n} = D_{21,n} = D_{22,n} = D_n,
N_{20,n} = N_{21,n} = N_{22,n} = N_n,
```
which the card describes as the consequence of an `O(3)`-invariant kernel forcing every grouped anisotropy defect to vanish. The script verifies the lane-independence of the constituent couplings `C_alpha`, `G_U`, `G_W`, `R_r` (by writing them with no lane subscript) but never explicitly forms `D_{20,n}, D_{21,n}, D_{22,n}` from per-lane data and asserts equality. The grouped collapse is implicit in the algebra — Stage 024's lane-free `D_n` is asserted to be the common value, but no test confronts a per-lane `D_{A,n}` against the common `D_n`.

The closest the script comes is Section II's equal-lane substitution into the `(1,2,2)`-weighted maps, which is not the same statement: equality of `xbar`/`a_x`/`b_x` on equal input is a property of the linear map, not a property of `D_n` arising from `O(3)`-invariant data.

**Why this matters:**
The grouped collapse is one of the paper's boxed outputs. Currently, an attacker who alters one lane's `G_U` (e.g., `G_{U,20} = G_U + delta`, `G_{U,21} = G_{U,22} = G_U`) and ran the per-lane D_n forward would observe `D_{20,n} != D_{21,n}` and the `a_D`/`b_D` defects becoming nonzero. The script does not exercise this scenario, so it does not directly verify the paper's claim that lane independence collapses to grouped isotropy.

This is lower severity than F1/F2 because the paper's mathematical content is sound and the lane-independence claim is structurally obvious (no lane index appears in `C_alpha`, `G_U`, etc.), but a checkpoint-quality script should make the collapse explicit.

**Required change:**
Add a short Section II.3 (or extend Section III) that:
1. Defines lane-decorated couplings `C_{alpha,A}`, `G_{U,r,A}`, `G_{W,r,A}`, `R_{r,A}` for `A in {20, 21, 22}`.
2. Forms the corresponding per-lane `B_{A,n}, Z_{A,n}, N_{A,n}, D_{A,n}` by substituting the lane-decorated couplings into the existing closed forms.
3. Applies the `O(3)`-invariance specialization `C_{alpha,20} = C_{alpha,21} = C_{alpha,22} = C_alpha` (and similarly for G, R).
4. Asserts that the per-lane `D_{A,n}` and `N_{A,n}` all collapse to the same scalar `D_n` and `N_n`, and that the `(a_D, b_D)` and `(a_N, b_N)` grouped defects vanish.
5. Additionally, perturb one lane with a symbolic `delta` and verify the defects are linear in `delta`, confirming the script can detect lane breaking when it occurs.

Mirror in Mathematica.

**Verification:**
After fix, the transcript contains a new `Section II.3 — Lane collapse under O(3) invariance` block with assertions `D_{20,n} - D_n = 0`, `D_{21,n} - D_n = 0`, `D_{22,n} - D_n = 0`, `a_{D,n} = 0`, `b_{D,n} = 0`, plus the perturbation witness `a_{D,n}|_{delta} = (specific linear function of delta)` to demonstrate non-tautology.

## Independent-derivation check (Mathematica)

The Mathematica script's Sections I.1 (Gram matrix orthonormality) and IV.1 (axisymmetric triple-overlap matrix) ARE genuine independent re-derivations. Compare SymPy:
```
def I4(i, j, k, l):
    s = 0
    for pr in PAIRINGS4:
        prod = 1
        for a, b in pr:
            prod *= delta[inds[a], inds[b]]
        s += prod
    return sp.simplify(4 * pi * s / 15)
```
vs Mathematica:
```
n[1] = Sin[theta] Cos[phi];
n[2] = Sin[theta] Sin[phi];
n[3] = Cos[theta];
i4[i_, j_, k_, l_] := Integrate[
  n[i] n[j] n[k] n[l] Sin[theta],
  {theta, 0, Pi}, {phi, 0, 2 Pi}
];
```
SymPy uses the analytic delta-pairing closed form of `int n_i n_j n_k n_l dOmega`. Mathematica computes the same integral by direct symbolic spherical integration. The two methods are mathematically equivalent but algorithmically independent; agreement (both produce the 5x5 identity for the Gram matrix, and both produce the diagonal `(sqrt(5)/(7*sqrt(pi))) * (1, 1/2, 1/2, -1, -1)` for M^(20)) is strong evidence the angular content is correct.

For Sections II, III, IV.2, V, the Mathematica script transliterates the SymPy algebra. See F1.

## Engine cross-check

Both engines pass every assertion. Specifically:

| Check | SymPy result | Mathematica result | Agree? |
|---|---|---|---|
| Gram - I5 | 5x5 zero matrix | 5x5 zero matrix | yes |
| projected - svec | 5-vec of zeros | 5-vec of zeros | yes |
| Section II equal-lane and witnesses | all zero | all zero | yes |
| C_alpha - lambda_B,alpha I_etaalpha | zero | zero | yes |
| B_0/B_2/B_4 sum formulas | zero | zero | yes |
| Z_0/Z_2/Z_4 formulas | zero | zero | yes |
| N_0/N_2/N_4 formulas | zero | zero | yes |
| D_0/D_2/D_4 formulas | zero | zero | yes |
| M^(20) - M_target | 5x5 zero matrix | 5x5 zero matrix | yes |
| Section IV.2 grouped weights | zero | zero | yes |
| Section V P-expansion and defects | zero | zero | yes |

No `engine_disagreement` finding. The agreement on Sections I.1 and IV.1 is informative because the two engines use genuinely different integration methods. The agreement on Sections II, III, IV.2, V is uninformative because both engines run the same algebra (see F1).

## Verdict justification

Verdict is `findings` (4 medium/low). Paper alignment is `aligned`: every paper-side boxed equation has a script-side check, no numeric or sign mismatches, no notes-vs-script convention disagreements. The angular sphere integrals (Sections I.1 and IV.1) are genuinely substantive and independently cross-checked between the two engines. However the closed-form coefficient assertions for `B_n, Z_n, N_n, D_n, P_n` (Sections III and V) are algebraic identities of self-defined rationals, and the Mathematica engine duplicates rather than independently re-derives those sections. For a checkpoint stage, this falls short of the required bar but is recoverable: F2 prescribes a matrix-inverse anchor that converts the algebraic identities into real physics checks, F1 routes that derivation through Mathematica independently, F3 trims definitional tautologies, and F4 adds the explicit lane-collapse witness. None of the findings are `CRITICAL_DOWNSTREAM` because the underlying mathematics is correct as far as the script reaches — the paper card's boxed results are reproduced by the (algebraic) checks; the issue is verification quality, not result correctness. No `UNFIXABLE` because the prescribed fixes are mechanical.

Attacks tried that failed: I tried to find a sign disagreement on R_r between paper and script (paper uses `R_r` symmetrically as `2 G_U G_W R_r` in Q_r and `R_r G_U` in P_r, script matches both); I tried to find a different `(1,2,2)` weighting in the notes vs script (notes do not state the weights explicitly but refer to Stage 022; the (1,2,2)/5 weighting in the script matches the standard grouped-metric convention and is consistent with `xbar`'s row-sum = 1 normalization); I tried to find an axisymmetric matrix off-diagonal that would break the diagonal claim (the triple integral for off-diagonal entries `<Y_{20} Y_{20} Y_{21c}>` etc. is zero by phi-integration parity, which the spherical integrate-by-numerical-methods Mathematica path independently confirms); I tried to find a missing factor on `kappa_* = sqrt(5)/(7 sqrt(pi))` by re-deriving from the standard 3j-symbol convention `<l=2,m;l=2,0|l=2,m'> = sqrt(5/(4 pi)) * <2,m;2,0|2,m'>` and got the same prefactor. The card's claims hold up against attack within the scope the script actually verifies; the verification depth is the issue, not the math.

## Self-test notes

I checked: (1) variable independence — the eps-Taylor series in Section V correctly tracks `lambda` (the lane weight enters `lane_ratio(lam)` as a real argument); no spurious zero-derivative traps. (2) Parity on the sphere integrals — the I4 fourth-moment integral and I6 sixth-moment integral are symmetric in their indices and produce the correct delta-pairing structure; the off-diagonal entries of M^(20) vanish by `cos/sin(m phi)` orthogonality, which I confirmed by inspection of the phi-integration. (3) Trivial-case pre-check — for F2's proposed `Z_full = coupling.T * M_pair.inv() * coupling`, setting `R_r = 0` reduces to `G_U^2/(Omega_U^2 - omega^2) + G_W^2/(Omega_W^2 - omega^2)`, whose constant term `G_U^2/Omega_U^2 + G_W^2/Omega_W^2` matches `Q/Delta` at `R_r=0` (where `Q = G_U^2 Omega_W^2 + G_W^2 Omega_U^2` and `Delta = Omega_U^2 Omega_W^2`, so `Q/Delta = G_U^2/Omega_U^2 + G_W^2/Omega_W^2`). The proposed anchor reduces correctly. (4) Path specifications in directive use `scripts/` for `.py` and `mathematica/` for `.wl`. (5) Paper round-trip — none of the proposed fixes introduce a constant or relation not already in the paper card or notes; the matrix-inverse anchor (F2) and the lane-collapse witness (F4) both derive their content from boxed paper equations.
