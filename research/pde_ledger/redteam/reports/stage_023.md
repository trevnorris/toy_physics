---
unit_id: 023
batch: I.2
auditor_model: claude-opus-4-7[1m]
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md"]
  paper_appendix: present
---

# Audit unit 023 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_023.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 68 + intro paragraph)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.txt`

## What the paper claims

Stage 023's `\stagefield{Output}` states: "Stage~023 outputs the projectors \eqref{eq:app-stage023-projectors}, the full coupled coefficient formulas \eqref{eq:app-stage023-b-moments}, \eqref{eq:app-stage023-z-moments-port}, \eqref{eq:app-stage023-n-moments-port}, and \eqref{eq:app-stage023-d-coeffs}, the isotropic normalization ratio \eqref{eq:app-stage023-normalization-test}, the constant-prefactor conditions \eqref{eq:app-stage023-constant-prefactor-conditions}, and the anisotropy transport laws \eqref{eq:app-stage023-u-anisotropy} and \eqref{eq:app-stage023-p-anisotropy}." Distinct deliverables: (1) the three exact `G_grp = diag(1,2,2)` projectors `P_bar, P_a, P_b` and their identities `P_iP_j = delta_ij P_i`, `P_bar + P_a + P_b = I`; (2) BdG moment definitions `B_{A0}, B_{A2}, B_{A4}`; (3) per-port conservative Maxwell/mixed moments `Z_{An}^{(r)}` for n=0,2,4 with `Z_{An} = sum_r Z_{An}^{(r)}`; (4) per-port outgoing-transfer moments `N_{An}^{(r)}` for n=0,2,4 with the per-port closed forms; (5) total operator coefficients `D_{A0}=K_A-B_{A0}-Z_{A0}`, `D_{A2}=-(M_A+B_{A2}+Z_{A2})`, `D_{A4}=-(B_{A4}+Z_{A4})`; (6) grouped-decomposition identities for `(Dbar, a_D, b_D)` and `(Nbar, a_N, b_N)`; (7) isotropic-branch closed forms for `u_2, u_4, P_0, P_2, P_4`; (8) the universal normalization product `mhat_0^2 N_0/(K-B_0-Z_0) = 54 G c_s^5/(5 a^5 c^5)`; (9) constant-prefactor conditions `N_2 = 2 D_2 N_0/D_0` and the `N_4` formula in eq:app-stage023-constant-prefactor-conditions; (10) first-order anisotropy transport laws for `u_2` and `P_0`. The notes additionally enumerate per-lane Lagrangian §2 (wall+BdG+U-W mixed sector) as the physical source from which Z,N moments arise.

## What the script claims to verify

The SymPy script verifies, in five sections: (I) `G_grp = diag(1,2,2)` projector orthogonality, norms, idempotency, completeness, and decomposition action `P_i x = x_i e_i` on an arbitrary grouped vector; (II) the per-port Z and N rational-function expansions match the closed forms by `series` then `coeff`, plus the grouped trace/anomaly decomposition `Dbar_0 = Kbar - Bbar_0 - Zbar_0` etc., plus a linearity additivity check on `groupedParts`; (III) the isotropic-branch series for `u_2, u_4, P_0, P_2, P_4` reduce to the closed forms; the constant-prefactor conditions are obtained by `solve(P2==0, N2)` then `solve(P4.subs(N2,N2_target)==0, N4)` and cross-checked against an independent closed form; the universal normalization is anchored by an independent Stage-5 Gamma5_port derivation from `j_2 + i y_2`; (IV) first-order anisotropy transport laws for `u_2` and `P_0` via `series` in `eps` to first order; (V) monotonicity derivatives of `P_0 = N_0/(K-B_0-Z_0)`. The Mathematica script mirrors structure with two additional independent paths: numerical substitution of `Z_n^{(r)}, N_n^{(r)}` at fixed `(Omega_U, Omega_W, R, g_U, g_W)` against `SeriesCoefficient`, and a direct Bessel small-z Taylor expansion of `h_2` to recover Gamma5_port = a^5/(27 c_s^5).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Projectors `P_bar, P_a, P_b` and identities | `grouped_projector_calculus()` / Section I | match |
| (2) BdG moments `B_{An}` | Symbolically declared; structure (`sum_alpha c^2/varpi^{2n+2}`) is not exercised — treated as Stage-3 carry-forward | partial |
| (3) Per-port `Z_{An}^{(r)}` closed forms | Section II.0 series expansion against closed forms | match |
| (4) Per-port `N_{An}^{(r)}` closed forms | Section II.0 series expansion against closed forms | match |
| (5) Total `D_{An} = ...` | Symbolically constructed in Section II and decomposed; the construction itself is the definition | match (by construction) |
| (6) Grouped trace/anomaly decomposition for `D_{An}` | Section II.1 explicit identity checks | match |
| (7) Isotropic-branch `u_2, u_4, P_0, P_2, P_4` | Section III.1 series→coeff→closed-form match | match |
| (8) Universal normalization `mhat^2 N_0/D_0 = 54Gc_s^5/(5a^5c^5)` | Section III.3 with Stage-5 Bessel `Gamma5_port` anchor | match |
| (9) Constant-prefactor `N_2, N_4` conditions | Section III.2 solver + closed-form cross-check | match |
| (10) Anisotropy transport laws | Section IV (generic + grouped-defect) | match |
| Lagrangian §2 → (Z_n, N_n) Schur derivation | Not exercised; rational form is taken as given input | partial (carry-forward; see F2) |

`paper_alignment: aligned` — every paper-side deliverable has a script-side counterpart that matches the stated identity. The two partial rows represent acceptable carry-forwards from upstream stages, not misalignments.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 75-80 | `expect_zero` on basis orthogonality & norms | (1) basis prep | yes |
| A2 | sympy | 94-100 | projector idempotency/orthogonality/completeness | (1) | yes |
| A3 | sympy | 109-112 | `P_i x - x_i e_i = 0`, decomposition completeness | (1) | yes |
| A4 | sympy | 140-149 | one-port Z_n series vs closed form | (3) | yes |
| A5 | sympy | 159-175 | one-port N_n series vs closed form | (4) | yes |
| A6 | sympy | 233-243 | grouped decomposition for D_{An} | (5)(6) | yes |
| A7 | sympy | 250-252 | groupedParts additivity (non-tautological linearity) | (6) | yes |
| A8 | sympy | 288-295 | isotropic u_2,u_4,P_0,P_2,P_4 vs closed forms | (7) | yes |
| A9 | sympy | 310-311 | N_2_target/N_4_target solver vs independent closed forms | (9) | yes |
| A10 | sympy | 317-320 | `mhat^2 N_0_target/D_0 = 54Gc_s^5/(5a^5c^5)` | (8) | partial — see F1 |
| A11 | sympy | 333 | Gamma5_port = a^5/(27 c_s^5) | (8) anchor | yes |
| A12 | sympy | 340-343 | ratio_target at mhat=1 reproduces universal target | (8) | yes |
| A13 | sympy | 376-377 | generic du_2, dP_0 first-order forms | (10) | yes |
| A14 | sympy | 392-395 | grouped-defect a_u2, b_u2, a_P0, b_P0 forms | (10) | yes |
| A15 | sympy | 424-427 | monotonicity derivatives of P_0 = N_0/(K-B_0-Z_0) | (5)(8) | yes |
| B1-B30 | mathematica | 48-248 | mirrors A1-A15 + numerical cross-check (97-117) + direct Bessel path (212-219) | all of (1)-(10) | yes |

A10 is the only borderline row: `N0_target` is defined as the solution of `mhat^2 N0/D0 = target` for `N0`, then the script asserts that substituting `N0_target` reproduces `target`. That is algebraically guaranteed.

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:314-320`

**What's wrong:**
The block

```python
N0_target = sp.simplify(sp.solve(sp.Eq(mhat**2 * (N0/D0), 54 * G * c_s**5 / (5 * a**5 * c**5)), N0)[0])
expect_zero(
    "N0_target reproduces universal normalization",
    (mhat**2 * (N0_target/D0)) - 54 * G * c_s**5 / (5 * a**5 * c**5),
)
```

is tautological by construction: `N0_target` is *defined* as the value of `N_0` that makes the LHS equal the RHS, so substituting it back is algebraically guaranteed to give zero. The check cannot fail regardless of whether the universal target value `54 G c_s^5/(5 a^5 c^5)` is correct.

**Why this matters:**
This block carries the appearance of independently verifying the universal normalization target, but it doesn't — the value `54 G c_s^5/(5 a^5 c^5)` is on both sides of the equation by construction. The real verification of the target is done a few lines below (lines 322-343) via the Stage-5 Gamma5_port = a^5/(27 c_s^5) derivation and `ratio_target = gamma_GR/(mhat^2 * Gamma5_port)` check, which IS substantive. The taut block is dead weight that risks giving a reader false confidence that the universal value was cross-checked twice.

**Required change:**
Replace the tautological block with a check that the equation `mhat^2 * P_0 = 54 G c_s^5/(5 a^5 c^5)` is equivalent to `mhat^2 * N_0/(K - B_0 - Z_0) = 54 G c_s^5/(5 a^5 c^5)` after substituting `P_0 = N_0/D_0` and `D_0 = K - B_0 - Z_0`. That exercises the paper's eq:app-stage023-normalization-test substitution (the equality between the abstract `P_0` form and the explicit `(K, B_0, Z_0)` form), which is a non-tautological identity the script does not currently test.

Concretely, after the existing `N0_target` line, insert symbolic substitution:

```python
K_sym, B0_sym, Z0_sym = sp.symbols("K_sym B0_sym Z0_sym", positive=True, real=True)
norm_abstract = mhat**2 * N0 / D0 - 54 * G * c_s**5 / (5 * a**5 * c**5)
norm_explicit = mhat**2 * N0 / (K_sym - B0_sym - Z0_sym) - 54 * G * c_s**5 / (5 * a**5 * c**5)
expect_zero(
    "normalization abstract == explicit under D0 = K - B0 - Z0",
    sp.simplify(norm_abstract.subs(D0, K_sym - B0_sym - Z0_sym) - norm_explicit),
)
```

Delete (or relabel as `# illustrative only`) the existing taut block on lines 317-320.

**Verification:**
The new `expect_zero` line should appear after line 314 and produce a fresh `... = 0` in the .txt output. Codex must also mirror the change in the Mathematica `.wl` script (add the analogous explicit/abstract equivalence assertion under Section III after line 191).

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:122-175`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:77-91`

**What's wrong:**
Paper §2 introduces the reduced per-lane Lagrangian with `(q_A, X_{A,alpha}, U_{A,r}, W_{A,r}, R_{A,r})` and mixing terms `R U W`, `g_U q U`, `g_W q W`. Paper §4 then *states* the per-port `(Delta, S, Q, G, P)` definitions and the resulting `Z_{An}^{(r)}` and `N_{An}^{(r)}` closed forms — but the derivation from the Lagrangian (Euler-Lagrange + Schur complement of the `(U,W)` block) is not shown in this stage's card; the notes treat them as carry-forwards from Stage 003 / Stage 021.

The script in Section II.0 verifies that the rational functions

- `(Q - H ω^2)/(Δ - S ω^2 + ω^4)` Taylor-expand to the `Z_n` formulas, and
- `(P - g_W ω^2)^2/(Δ - S ω^2 + ω^4)^2` Taylor-expand to the `N_n` formulas

but **never derives the rational functions themselves from the Lagrangian.** The denominator `Δ - S ω^2 + ω^4` is just typed in; nothing checks it is the Schur complement determinant of the `(U,W)` mass matrix at frequency `ω`. So the chain "Lagrangian → Schur complement → rational function → series → closed form" only has the last two links verified.

This is acceptable as a carry-forward (script comment at sympy line 181-182: `# Stage-003 carry-forward`), but the carry-forward is only loosely cited — there is no explicit reference to the upstream script's specific lines that derive the rational form from the Lagrangian, and the .wl makes no carry-forward comment at all in the analogous block.

**Why this matters:**
A reader checking the stage cannot trace where the rational-function ansatz comes from. If the upstream Schur complement was wrong (e.g., wrong sign on the `R` cross-term, or wrong S = Ω_U² + Ω_W² formula vs `S = Ω_U² + Ω_W² + 2R` or similar), this stage's audit would pass and the error would be invisible. The stage card's §2 Lagrangian → §4 `(Δ, S, Q, G, P)` step is the only un-anchored algebra in the whole stage.

**Required change:**
Add a brief in-script Schur-complement derivation under Section II.0 (sympy lines after 130; mathematica lines after 82). Specifically, construct the 2x2 frequency-dependent matrix

```
M(ω) = [[Ω_U² - ω², R], [R, Ω_W² - ω²]]
```

representing the `(U,W)` block, compute `det M(ω) = (Ω_U² - ω²)(Ω_W² - ω²) - R² = Δ - S ω² + ω^4` and assert this equals the script's `Delta_expr - S_expr*omega^2 + omega**4`. Similarly construct the `q`-to-`(U,W)` coupling vector `(g_U, g_W)` and verify the Schur reduction `(g_U, g_W) · M(ω)^{-1} · (g_U, g_W)^T` produces the `(Q - H ω²)/(Δ - S ω² + ω^4)` rational function used by the script.

Concretely insert before line 132 of the .py:

```python
Mblock = sp.Matrix([[OmegaU**2 - omega**2, Rmix], [Rmix, OmegaW**2 - omega**2]])
det_M = sp.expand(Mblock.det())
expect_zero("Schur denominator matches Delta - S omega^2 + omega^4",
            sp.expand(det_M - (Delta_expr - S_expr * omega**2 + omega**4)))
g_vec = sp.Matrix([gU, gW])
Z_schur = sp.simplify(((g_vec.T * Mblock.adjugate() * g_vec)[0]) / det_M)
expect_zero("Z rational form matches Schur (g_U,g_W) M^{-1} (g_U,g_W)^T",
            sp.simplify(Z_schur - (Q_expr - H_expr * omega**2) / det_M))
```

and mirror in the .wl with Mathematica equivalents (`Mblock = {{omegaU^2 - omega^2, rMix}, {rMix, omegaW^2 - omega^2}}; detM = Det[Mblock]; ...`).

**Why low and not medium:** the carry-forward citation, while loose, does point at the right upstream stages, and the rational-function form is a well-known coupled-oscillator Schur complement. The fix anchors the chain without breaking anything.

**Verification:**
Two new `expect_zero` lines in Section II.0 of each engine; new `= 0` lines in the .txt outputs. After the fix, the chain Lagrangian → Schur → rational → series → Z_n closed form is fully in-script.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a port of the SymPy script — same five sections in the same order, same `groupedParts` helper with same lane weights, same assertion sequence in Sections I, II.1, III.1, IV, V. **However**, it adds two genuinely independent verification paths:

1. Lines 93-117: numerical substitution of `Z_n, N_n` closed forms at concrete `(Omega_U=2, Omega_W=3, R=1, g_U=1, g_W=2)` cross-checked against `SeriesCoefficient` of the rational function at those values. This is a numerical verification that the closed forms agree with the rational expansion, independent of the symbolic `Series` path used by SymPy.

2. Lines 209-219: a direct small-z Taylor expansion of `j_2 + i y_2` (truncated at z^9 / z^8 for `h2` and `D[h2,z]`) feeding into `Λ_2`, then into `Y_2`, then extracting the `ω^5` coefficient. This bypasses the `omega*D[h2,z]/h2` symbolic ratio path and verifies `Gamma5_port = a^5/(27 c_s^5)` from a different algebraic route.

These two independent paths are sufficient to defeat a pure-transliteration finding. **No `mathematica_transliteration` finding is raised**, though the auditor notes the section/assertion skeleton is shared and a stricter reviewer might disagree.

## Engine cross-check

Both engines pass every assertion (`= 0` for sympy, `PASS:` for mathematica). The Stage-5 Gamma5_port anchor `a^5/(27 c_s^5)` appears in both. The universal normalization target `54 G c_s^5/(5 a^5 c^5)` appears in both. The constant-prefactor `N_2 = 2 D_2 N_0/D_0` and `N_4 = N_0(2 D_0 D_4 + D_2^2)/D_0^2` appear in both. The grouped projector formulas

- `Pbar = (1/5)·[[1,2,2],[1,2,2],[1,2,2]]`
- `Pa = (1/20)·[[16,-8,-8],[-4,2,2],[-4,2,2]]`
- `Pb = (1/4)·[[0,0,0],[0,2,-2],[0,-2,2]]`

are present and verified idempotent in both engines and match the paper's boxed eq:app-stage023-projectors. Engines agree.

## Verdict justification

`verdict: findings` with two low-severity issues. The stage's load-bearing math — projector calculus, grouped decomposition, isotropic-branch closed forms, Stage-5 Gamma5_port anchor, constant-prefactor conditions, anisotropy transport laws — is exercised by substantive, non-tautological assertions in both engines, and every paper-side deliverable maps to a matching script-side check. The two findings are: (F1) one block in §III.3 of the sympy script is a tautological "solve-then-substitute" against the universal target that adds no verification value (the real verification is the Stage-5 Bessel chain four lines below, which IS substantive); (F2) the Schur-complement step from the §2 Lagrangian to the rational function `(Q - H ω²)/(Δ - S ω² + ω⁴)` is taken as a carry-forward from Stages 003/021 rather than re-derived in-script. Neither finding undermines the stage's conclusion; both are upgrade-quality issues that strengthen the audit chain rather than fix wrong math. No `paper_misalignment`. No `engine_disagreement`. Outputs are fresh. Attacks tried: (a) verified the constant-prefactor `N_4` reduction `[2 D_0(D_2·(2 D_2 N_0/D_0) + D_4 N_0) - 3 D_2^2 N_0]/D_0^2 = N_0(D_2^2 + 2 D_0 D_4)/D_0^2` is algebraically correct — passes; (b) checked the projector matrix elements match the paper boxed form — match; (c) checked sign conventions on `du_2` and `dP_0` — match; (d) checked series truncation orders are sufficient for the extracted coefficients — sufficient (order 6 for ω^4 extraction in sympy, order 4 sufficient for ω^4 extraction in mathematica). The docstring still says "Stage 6" / "STAGE 006" rather than "Stage 023" (sympy line 3 docstring, mathematica line 37 banner) — cosmetic stale label, no math impact, not flagged.

## Self-test notes

- **Variable independence:** F2's proposed `Mblock = [[Ω_U²-ω², R], [R, Ω_W²-ω²]]` and `det_M` use symbols already declared as real in the sympy script (`omega, OmegaU, OmegaW, Rmix` at line 124); the `(g_U, g_W)` Schur quadratic form uses `gU, gW` real symbols at the same line. All symbols in the proposed expressions are independent and non-trivially present. No spurious zero-derivative trap.
- **Symmetry/parity:** N/A — no unbounded integrals proposed.
- **Trivial-case pre-check:** For F2's new `expect_zero("Schur denominator matches ...", det_M - (Delta_expr - S_expr*omega^2 + omega^4))`: substitute `R=0, OmegaU=1, OmegaW=1`: det_M = (1-ω²)² = 1 - 2ω² + ω⁴; RHS = (1·1 - 0) - (1+1)ω² + ω⁴ = 1 - 2ω² + ω⁴. Match. For the Schur quadratic-form check at the same substitution: M^{-1} = (1/(1-ω²)) · I; (g_U, g_W) M^{-1} (g_U, g_W)^T = (g_U² + g_W²)/(1-ω²). Paper's Q-form: Q = g_U²·1 + 2 g_U g_W·0 + g_W²·1 = g_U²+g_W²; H = g_U²+g_W². So (Q - H ω²)/det_M = (g_U²+g_W² - (g_U²+g_W²)ω²)/(1-ω²)² = (g_U²+g_W²)(1-ω²)/(1-ω²)² = (g_U²+g_W²)/(1-ω²). Matches. Trivial-case pre-check passes.
- **Path specifications:** F2 touches existing `.py` and `.wl` files at named line regions; no new files. F1 modifies existing `.py` script at named line region.
- **Paper round-trip:** F1's replacement check uses the same `54 G c_s^5/(5 a^5 c^5)` constant the paper states (eq:app-stage023-normalization-test), the same `D_0 = K - B_0 - Z_0` identity (eq:app-stage023-d-coeffs), and the same `P_0 = N_0/D_0` identity (eq:app-stage023-isotropic-pref). No new paper_misalignment introduced. F2's Schur derivation uses the Lagrangian §2 mixing terms and the §4 `(Δ, S, Q, G, P)` definitions exactly as the paper states them; no new constants introduced.
