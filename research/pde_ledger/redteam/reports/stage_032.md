---
unit_id: 032
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T22:39:50-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage032_source_map_from_mode_integrals.md"]
  paper_appendix: present
---

# Audit unit 032 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_032.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage032_source_map_from_mode_integrals.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 54; insertion at line 102)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.txt`

## What the paper claims

The stage card's `\stagefield{Output}` (line 65) is: "Stage~032 outputs the source-map factor [Eq. mhat-minus], the bound [Eq. mhat-bound], and the invariant selected product [Eq. Nminus]." The three boxed equations are:

1. (Eq. mhat-minus, lines 36-41) `mhat_- = (v . e_-)/kappa_0` and `mhat_-^2 = s_-/kappa_0^2`, derived by projecting the natural D/N source vector `J_src = g_Q Q_STF v` onto the selected eigenvector `e_-` and normalizing by `J_-(0) = g_Q Q_STF kappa_0`.
2. (Eq. mhat-bound, lines 44-47) `1 <= mhat_-^2 < sigma/kappa_0^2 = 11/9`, since `s_-` increases from `kappa_0^2` toward `sigma = kappa_0^2 + kappa_1^2`.
3. (Eq. Nminus, lines 53-57) `N_-(alpha_0) := mhat_-^2 P_{0,-} = beta_0 s_-^2 / (kappa_0^2 lambda_-)`, plus the target identification `N_-(alpha_0) = N_Q^target` (Eq. Nminus-target).

Inputs (line 9): selected eigenvector `e_-`, overlap vector `v`, and a source coupled through the D/N load profile. The notes file gives the same content with additional derivation context (Sections 1-7) and explicit `kappa_0 = 2 sqrt(2)/pi`, `kappa_1 = -4/(3 pi)`, `sigma = 88/(9 pi^2)`. The appendix row (line 54) states: "Natural D/N source map `mhat_-^2 = s_-/kappa_0^2`."

## What the script claims to verify

Per the SymPy docstring (lines 3-11), the script audits five things: (1) exact finite-throat axial integrals, (2) local-kernel reduction of the wall/internal/source couplings, (3) Schur-complement decomposition of the reduced wall operator, (4) the natural-D/N source map, and (5) elimination of the abstract source-map factor.

Concretely, the assertions verify: the values of `kappa_0`, `kappa_1`, `sigma`, `sigma/kappa_0^2 = 11/9` (15.1); that the local-kernel Lagrangian terms reduce to the `(v.q)` overlap structure including `L_src - g_Q Q (v.q) = 0` (15.2); that the Schur complement equals `Xi I + alpha vv^T` (15.3); that the formula `mhat_sq := s_minus_nat/kappa_0^2` equals 1 at alpha0=0 and has limit 11/9 as alpha0->infty (15.4); and that `Nprod_nat := (s_-/kappa_0^2)*(beta_0 s_-/lambda_-)` evaluates to `beta_0 kappa_0^2 / A` at alpha0=0 and to 0 as alpha0->infty (15.5).

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Eq. mhat-minus `mhat_-^2 = s_-/kappa_0^2` | `mhat_sq = s_minus_nat/kappa_0^2` (defined as a rename); endpoint values 1 at alpha=0 and 11/9 at alpha->infty | partial (the *projection step* `J_- = g_Q Q (v.e_-)`, hence `mhat_- = (v.e_-)/kappa_0`, is not exercised — `e_-` is never built and `(v.e_-)^2` is never compared to the imported `s_-` formula) |
| Eq. mhat-bound `1 <= mhat_-^2 < 11/9` | endpoint checks only (alpha0=0 and alpha0->infty) | partial (the strict inequality and monotonicity on the interior are not exercised; no check that mhat_-^2 stays in (1, 11/9) for finite alpha0 > 0) |
| Eq. Nminus `N_- = beta_0 s_-^2 / (kappa_0^2 lambda_-)` | `Nprod_nat = (s_-/kappa_0^2) * (beta_0 s_-/lambda_-)` — by construction; endpoint sanity checks | match (the identity is structural; endpoint checks are weak but consistent) |
| Eq. Nminus-target `N_-(alpha_0) = N_Q^target` | not present | n/a — paper card states this as a definitional target equation solved by downstream stages, not an identity Stage 032 is expected to prove. Not flagged. |
| Source vector form `J_src = g_Q Q_STF v` | `L_src - gQ*Qstf*(kappa_0 q_0 + kappa_1 q_1) = 0` | match |
| Local-kernel reductions (`L_etaU`, `L_etaphi`, `L_etaW`, `L_UW`) | five checks in 15.2 | match (script verifies more than the .tex spells out, but the notes file Section 2 explicitly enumerates these — extra coverage relative to .tex, supported by notes) |
| Schur decomposition `Sigma = Xi I + alpha vv^T` | `Sigma - [Xi I + alpha vv^T] = 0` in 15.3 | extra (not in stage_032.tex Output line; notes Section 3 covers this; was first proved at Stage 12/29; this stage 032 rederives in explicit kernel language) |

Dominant pattern: assertions exist for each paper claim, but two of them (Eq. mhat-minus and Eq. mhat-bound) are exercised only at endpoints or by definition. The Schur and kernel-reduction extras are sourced from the notes file's Sections 2-3, not the .tex Output line, so they are not orphan scaffolding.

`paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58-61 | basis orthonormality | input setup (kappa derivation) | yes |
| A2 | sympy | 62-65 | `kappa_0 - 2 sqrt(2)/pi = 0`, `kappa_1 + 4/(3 pi) = 0`, `sigma - 88/(9 pi^2) = 0`, `sigma/kappa_0^2 - 11/9 = 0` | inputs to all three deliverables | yes |
| A3 | sympy | 91-110 | local-kernel reductions including `L_src - gQ Q (v.q) = 0` | notes Section 2; supports source-vector derivation `J_src = g_Q Q_STF v` | yes |
| A4 | sympy | 143 | `Sigma - [Xi I + alpha vv^T] = 0` | notes Section 3 (extra, not in .tex Output) | yes (within notes scope) |
| A5 | sympy | 166 | `mhat_-^2(alpha=0) - 1 = 0` | Eq. mhat-minus endpoint + Eq. mhat-bound lower endpoint | yes (endpoint only) |
| A6 | sympy | 169-172 | `delta_kappa^2 + 4 Kprod - sigma^2 = 0` under natural subs | scaffolding for s_minus simplification (just (a-b)^2 + 4ab = (a+b)^2) | partial — algebraic identity, not a physics check |
| A7 | sympy | 180-190 | `s_minus_nat - s_minus_nat_simplified` and same at numeric point | redundant simplification check; same value computed two ways using the identity in A6 | partial — internal consistency, does not exercise the (v.e_-)^2 claim |
| A8 | sympy | 191 | `limit_{alpha->infty} mhat_-^2 = 11/9` | Eq. mhat-bound upper endpoint | yes (endpoint only) |
| A9 | sympy | 202-206 | `Nprod(alpha=0) - beta_0 kappa_0^2/A = 0` | Eq. Nminus at alpha=0 (consistency check) | partial (by construction, plus a substitution; consistent with paper) |
| A10 | sympy | 209-210 | `limit_{alpha->infty} Nprod_nat = 0` | Eq. Nminus asymptotic | partial (verifies a limit, not the closed-form identity, which is true by construction) |
| M1 | mathematica | 51-58 | parallel to A1, A2 | same | yes |
| M2 | mathematica | 81-85 | parallel to A3 | same | yes |
| M3 | mathematica | 121, 134, 135 | `Schur via Inverse vs LinearSolve = 0`, `Sigma - [Xi I + alpha vv^T] = 0`, `sigmaMatViaSolve - [...] = 0` | parallel to A4, with an *independent* second route via `LinearSolve` | yes (the LinearSolve route is a non-transliterated cross-check for 15.3 specifically) |
| M4 | mathematica | 151-176 | parallel to A5-A7 | same | yes (endpoint only) |
| M5 | mathematica | 177 | `limit_{alpha->infty} mhat_-^2 - 11/9 = 0` | parallel to A8 | yes |
| M6 | mathematica | 185, 189 | parallel to A9, A10 | same | partial |

Observation: A6 and A7 are essentially algebraic-identity scaffolding (showing two equivalent forms of `s_minus_nat` agree because `delta_kappa^2 + 4 Kprod = sigma^2`). They do not test that `s_- = (v.e_-)^2` — they test that two re-arrangements of the same closed-form expression for `s_-` agree, which they must by construction.

## Findings

### F1 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:145-191`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:137-177`

**What's wrong:**
The script's core claim (paper Eq. mhat-minus, paper card lines 36-41) is `mhat_- = (v . e_-)/kappa_0`, derived by projecting the source vector `J_src = g_Q Q_STF v` onto the selected eigenvector `e_-` of the loaded wall operator and dividing by `J_-(0) = g_Q Q_STF kappa_0`. The script never builds `e_-` and never computes `(v . e_-)^2`; it only imports a closed-form symbolic expression for `s_-` (sympy line 155) and then defines `mhat_sq := s_minus_nat / kappa_0^2` (sympy line 163), which makes the identity `mhat_-^2 = s_-/kappa_0^2` true by construction.

Quote (paper card, lines 36-41):
```
\boxed{
  \widehat m_-=\frac{v\cdot e_-}{\kappa_0},
  \qquad
  \widehat m_-^{\,2}=\frac{s_-}{\kappa_0^2}}.
```

Quote (sympy line 163):
```
mhat_sq = sp.simplify(s_minus_nat / kappa0**2)
```

The substantive verification of `s_- = (v . e_-)^2` would require constructing the loaded-wall 2x2 Sigma matrix at finite alpha0, diagonalizing it to obtain `e_-`, and confirming `(v . e_-)^2` matches the closed-form `s_-` expression. None of this is done. The script is consistent with the carry-forward formula coming from stages 028-031, but Stage 032's own headline claim ("mhat_-^2 = s_-/kappa_0^2") is verified only by the trivial rename.

Additionally, the bound `1 <= mhat_-^2 < 11/9` (paper Eq. mhat-bound, lines 44-47) is verified only at the endpoints alpha0=0 (gives 1) and alpha0->infty (gives 11/9). The script does not verify the *strict* inequality `mhat_-^2 < 11/9` for any finite alpha0 > 0, nor that `mhat_-^2 >= 1` strictly on the interior, nor monotonicity.

**Why this matters:**
The headline claim of the stage card is a derivation, not a definition. As written, the script accepts the derivation's output and confirms only its boundary values. A bug in the imported `s_-` formula (e.g., a sign or coefficient error in stage 028/029 carried forward) would not be caught here, even though Stage 032 is the unit the paper card and appendix row point to for verifying `mhat_-^2 = s_-/kappa_0^2`. The endpoint checks at alpha0=0 and alpha0->infty are *necessary but not sufficient* — there are many alternative `s_-` formulas with the same two endpoint values.

**Required change:**
Add an independent verification in Stage 15.4 (after the current `mhat_sq` definition) that constructs the loaded 2x2 wall operator
`M(alpha0) = diag(A, A+DK) - alpha0 * v v^T`,
takes its lower eigenvector `e_-`, computes `s_check := (v . e_-)^2`, and asserts `expect_zero("s_check - s_minus_nat", s_check - s_minus_nat)`. Mirror in the .wl using an independent route (e.g. `Eigensystem` or `Eigenvectors` with appropriate normalization). This is the missing step that ties the abstract `s_-` symbol back to the paper's derivation.

Optionally also add a sampling check that `1 <= mhat_sq <= 11/9` at a handful of finite alpha0 values to exercise the strict bound, though monotonicity is more naturally handled by stage 031.

**Verification:**
The new check should appear after sympy line 165 and after .wl line 150. The new assertion should evaluate to 0 symbolically (or at least at numeric alpha0 values). If the assertion passes, the carry-forward from upstream is confirmed and the headline claim is no longer accepted by construction.

### F2 — mathematica_transliteration

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl` (entire file)

**What's wrong:**
The Mathematica script is a near line-by-line port of the SymPy script. Compare:

SymPy (lines 121-131):
```
Kint = sp.Matrix([
    [AU, 0, 0, -gR * kappa0],
    [0, AU, 0, -gR * kappa1],
    [0, 0, Aphi, 0],
    [-gR * kappa0, -gR * kappa1, 0, AW],
])
Bmat = sp.Matrix([
    [gU, 0, gB * kappa0, gW * kappa0],
    [0, gU, gB * kappa1, gW * kappa1],
])
```

Mathematica (lines 94-102):
```
kInt = {
  {aU, 0, 0, -gR*kappa0},
  {0, aU, 0, -gR*kappa1},
  {0, 0, aPhi, 0},
  {-gR*kappa0, -gR*kappa1, 0, aW}
};
bMat = {
  {gU, 0, gB*kappa0, gW*kappa0},
  {0, gU, gB*kappa1, gW*kappa1}
};
```

SymPy (lines 153-155):
```
R = sp.sqrt((DK + alpha0 * delta_kappa)**2 + 4 * alpha0**2 * Kprod)
lambda_minus = sp.simplify((A + B - alpha0 * sigma_sym - R) / 2)
s_minus = sp.simplify(sp.Rational(1, 2) * (sigma_sym + ((DK + alpha0 * delta_kappa) * delta_kappa + 4 * alpha0 * Kprod) / R))
```

Mathematica (lines 143-145):
```
r = Sqrt[(dK + alpha0*deltaKappa)^2 + 4*alpha0^2*kappaProd];
lamMinus = FullSimplify[(2*a + dK - alpha0*sigmaSym - r)/2, Assumptions -> $Assumptions];
sMinus = FullSimplify[(sigmaSym + ((dK + alpha0*deltaKappa)*deltaKappa + 4*alpha0*kappaProd)/r)/2, Assumptions -> $Assumptions];
```

The same intermediate quantities, the same construction of `R`, the same `lambda_minus`/`s_minus` formulas in the same order, even the same interior-identity cross-check (sympy 176-190 vs. .wl 160-176) appear in lockstep. The variable names are renamed (camelCase) but the algebraic choreography is identical. The only genuinely independent computation in the .wl is the `LinearSolve` cross-check of the Schur complement in Stage 15.3 (.wl lines 104-121), which is a partial mitigation but does not cover Stages 15.4-15.5.

**Why this matters:**
Per the project's second-engine policy, both engines must derive results independently from the physical premises rather than echo each other's algebra. A transliterated .wl provides no protection against a logic error in the .py (and vice versa) — both would produce the same wrong answer in lockstep. For Stages 15.4-15.5, where the headline claim lives, the two engines are not independent.

**Required change:**
For Stage 15.4 specifically (the `mhat_-` derivation), have the Mathematica script take an independent route. One concrete option: compute `lamMinus` and `sMinus` directly from the loaded 2x2 wall matrix `M = {{a, 0}, {0, a+dK}} - alpha0 * v.Transpose[v]` via `Eigensystem[M]`, picking the lower eigenvalue/eigenvector and computing `(v.e_-)^2` directly. This both (a) gives an independent algebra path for the .wl, and (b) addresses F1 by constructing `e_-` explicitly. After computing `lamMinusIndep` and `sMinusIndep` this way, assert agreement with the closed-form `lamMinus`/`sMinus`. The existing checks may stay; this adds the independent path.

**Verification:**
The .wl should contain an `Eigensystem` (or `Eigenvectors`/`Eigenvalues`) call that does not appear in the .py, plus assertions that the independent eigenvalue/eigenvector results agree with the imported closed forms. After the change, the .wl's Stage 15.4 derives the same numbers from a different algebraic path.

### F3 — insufficient_verification

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:168-190`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:153-176`

**What's wrong:**
The "interior consistency" block (sympy 168-190, .wl 153-176) asserts:
1. `delta_kappa^2 + 4 Kprod - sigma^2 = 0` under the natural-D/N substitution.
2. `s_minus_nat - s_minus_nat_simplified = 0` (two re-arrangements of the same closed form agree after using identity (1)).
3. The same check at the numeric point `(alpha0=1, DK=1)`.

Identity (1) is the universal algebraic identity `(a-b)^2 + 4ab = (a+b)^2` evaluated at `a = kappa_0^2`, `b = kappa_1^2` — it holds for any reals; no D/N or wall physics is actually exercised. Identity (2) is then a consequence of (1) applied to two formally distinct simplifications of `R`. Identity (3) is the same check at a single numeric point, which is strictly weaker than (2).

Inline comments label this block as an "interior consistency" check, but it does not check anything about the interior of the alpha0 parameter range — it checks that two algebraic re-arrangements of the same expression agree.

**Why this matters:**
These three assertions consume audit attention but provide essentially no information about the physics: a bug in the s_minus formula would not show up here as long as the bug were carried consistently between the two algebraic re-arrangements. If F1 is addressed (independent construction of `s_-` from `e_-`), this entire block becomes redundant scaffolding.

**Required change:**
Either (a) remove the three assertions (sympy lines 169-190 and .wl lines 154-176) once F1's independent check is in place, or (b) leave them but rename their inline comments to reflect what they actually do ("algebraic re-arrangement consistency", not "interior consistency"). Recommendation: defer this until F1 is applied; if F1's independent check supersedes them, prefer removal.

**Verification:**
After F1 is applied, re-run sympy and mathematica and confirm the four removed assertions no longer appear in the output, or the comment text is updated. This is an optional cleanup; not a correctness gate.

## Independent-derivation check (Mathematica)

The .wl follows the .py step-for-step through Stages 15.1, 15.2, 15.4, and 15.5 (see F2 for line-paired excerpts). The single bona-fide independent move is in 15.3, where the .wl re-derives the Schur complement via `LinearSolve` and confirms agreement with the direct `Inverse` route (`expectZero["Schur via Inverse vs LinearSolve", sigmaMat - sigmaMatViaSolve]` at .wl line 121). That cross-check is genuine and useful, but it only covers Section 15.3. The headline derivation in 15.4 is parallel algebra, not parallel reasoning. F2 covers this.

## Engine cross-check

Both transcripts pass all assertions and produce equivalent algebraic forms:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `kappa_0` | `2*sqrt(2)/pi` | `(2*Sqrt[2])/Pi` |
| `kappa_1` | `-4/(3*pi)` | `-4/(3*Pi)` |
| `sigma` | `88/(9*pi**2)` | `88/(9*Pi^2)` |
| `sigma/kappa_0^2` | `11/9` | `11/9` |
| `Xi` | `g_U**2/A_U` | `gU^2/aU` |
| `alpha` | `g_B**2/A_phi + (A_U*g_W + g_R*g_U)**2/(A_U*(A_U*A_W - 88*g_R**2/(9*pi**2)))` | `gB^2/aPhi + (gR*gU + aU*gW)^2/(aU*(aU*aW - (88*gR^2)/(9*Pi^2)))` |
| `mhat_-^2` general | `(63*pi**2*DeltaK + 968*alpha0 + 11*sqrt(...))/(18*sqrt(...))` | `(11 + (968*alpha0 + 63*dK*Pi^2)/Sqrt[...])/18` (equivalent after factoring) |

All 16 SymPy assertions and all 18 Mathematica assertions print as `0`/`PASS`. No engine disagreement.

The `Limit::alimv` warnings in the .wl output (lines 75, 85) are benign: Mathematica's `Limit` does not honor positivity assumptions on the limit variable, but the limit values match SymPy.

## Stale output check

- SymPy script mtime: 2026-05-21 17:24:09
- SymPy output mtime: 2026-05-21 17:25:53 (newer — fresh)
- Mathematica script mtime: 2026-05-21 17:24:09
- Mathematica output mtime: 2026-05-21 17:26:05 (newer — fresh)

No `stale_output` finding.

## Verdict justification

Both scripts run end-to-end, produce equivalent results, and verify every input quantity (`kappa_0`, `kappa_1`, `sigma`, `sigma/kappa_0^2 = 11/9`), the local-kernel reductions, and the Schur structure of the wall self-energy. They confirm the endpoint values of the paper's mhat bound (1 at alpha0=0, 11/9 as alpha0->infty) and the asymptotic behavior of `Nprod_nat`. What they do *not* substantively verify is the paper's headline derivation step: that `s_- = (v . e_-)^2` for the eigenvector `e_-` of the loaded wall operator. The script accepts a closed-form `s_-` from upstream and defines `mhat_sq := s_-/kappa_0^2` as a rename, then checks endpoint consistency. The Mathematica script is structurally a port of the SymPy script for 15.4-15.5 (the headline-claim sections), with one genuine independent step in 15.3. Verdict is `findings` with `paper_alignment: partial`: nothing is wrong with the paper relative to the script, but the script does less work than the paper card promises for the headline derivation. Stop-cold is null — F1's fix is a constructive addition (build `e_-`, compute `(v.e_-)^2`, compare to `s_-`); it does not invalidate downstream stages.

## Self-test notes

Self-tests performed:
- **Variable independence**: F1's prescribed check builds `M(alpha0)` which genuinely depends on alpha0, A, DK, kappa_0, kappa_1. Computing `Eigensystem[M]` returns eigenvectors that depend on these symbols, so `(v.e_-)^2` is a non-trivial function of (alpha0, A, DK). Comparing to the closed-form `s_minus_nat` (same dependence) is a real check that can fail.
- **Symmetry/parity**: irrelevant — no integrals over symmetric domains in the new check.
- **Trivial-case pre-check**: at alpha0=0, the loaded matrix is `diag(A, A+DK)`, lower eigenvector is `(1, 0)` (assuming A < A+DK), so `(v.e_-)^2 = kappa_0^2`, matching `s_-(alpha0=0) = kappa_0^2`. At alpha0 large, lower eigenvector aligns with `v / sqrt(sigma)`, so `(v.e_-)^2 -> sigma`. Both endpoints match the known closed form. The check is non-trivial.
- **Paper round-trip**: F1's required new check uses the same kappa values and the same closed-form `s_-` already in the script; it does not introduce new constants the paper does not state. F2's required change (independent Eigensystem path) likewise stays inside the paper's symbol set. F3 is a comment cleanup with no paper impact.
