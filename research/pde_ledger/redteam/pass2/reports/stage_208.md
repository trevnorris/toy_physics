---
unit_id: 208
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md]
  paper_appendix: present
---

# Audit unit 208 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_208.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 47, 833-918)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.txt`

## What the paper claims

Stage 208 opens the first genuinely mixed (two-coordinate) search sector. `\stagefield{Output}`: "Pairwise mixed-ray cone, gradient and off-diagonal Hessian synergy laws, canonical two-ray screen, promotion/deferral rule, and minimal packet for the full ratio optimizer." The notes enumerate the distinct deliverables: (1) the pairwise mixed-ray cone `s_ij(r)=(e_i+r·e_j)/sqrt(1+r²)`, r≥0; (2) the gradient-synergy law `k_ij(r)=(k_i+r·k_j)/sqrt(1+r²)` with `dk_ij/dr=(k_j-k_i·r)/(1+r²)^{3/2}`, unique maximizer `r_grad=k_j/k_i`, max `sqrt(k_i²+k_j²)`; (3) the off-diagonal curvature-synergy law `H_{1,ij}=(h_ii+2r·h_ij+r²·h_jj)/(1+r²)` with cross-weight `w_x(r)=2r/(1+r²)` maximized at r=1, plus the diagonal-neutrality law when `h_ij=0`; (4) the mixed curvature envelopes `kappa_lo, kappa_hi` and the certified bracket `T(H_0,k;c)=2H_0/(k+sqrt(k²-2cH_0))`; (5) the two canonical screen rays (gradient-optimal `s_grad`, equal-mix `s_eq`) with their slopes and curvatures, coinciding iff `k_i=k_j`; (6) promotion/deferral rule and (7) the minimal packet `P_ij^mix`. The card states "Mathematica audit: none yet" — stale, since a `.wl` now exists (retrofit batch). This is a `\StatusExactClosure{}` mixed stage with a higher bar; no status-only carve-out.

## What the script claims to verify

The SymPy script (sections I-VI) verifies the symbolic identities for the slope law, its derivative, the gradient-optimal ratio/slope and the gradient gain, the mixed-curvature decomposition into `(w_i, w_x, w_j)` weights, diagonal neutrality, the cross-weight derivative/extremum at r=1, the two canonical screen rays' directions/slopes/curvatures, the coincidence condition, the entrywise curvature envelopes, and the fixed-ratio certified bracket via the quadratic root relation. The Mathematica script (M1-M9) verifies the same family of identities but additionally *solves* for the stationary ratio (`Solve`/`Reduce`) rather than positing it, and *extracts* the curvature weights via polynomial `Coefficient`/`Series` rather than positing them. Every `expect_zero`/`expectZero` compares a script-built object (matrix Rayleigh quotient, derivative, solved root) against the paper's closed form, so all are substantive (none tautological).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) pairwise cone `s_ij(r)` | sympy `s=[1,r]/sqrt(1+r²)`; wl `sRay[x]` | match |
| (2) slope law + derivative + r_grad + max | sympy A: "mixed slope law", "derivative law", "gradient-optimal slope/gain"; wl M1-M3 (Solve/Reduce) | match |
| (3) curvature decomp + w_x + r=1 max + neutrality | sympy "mixed curvature decomposition", "cross-weight...", "diagonal neutrality"; wl M4-M5 (Coefficient/Series) | match |
| (4) envelopes + certified bracket T | sympy IV-V; wl M8-M9 (rootMap) | match |
| (5) two canonical rays + coincidence iff k_i=k_j | sympy III, VI; wl M6-M7 | match |
| (6) promotion/deferral rule | — (interval-comparison logic, prose theorem) | not script-checked (acceptable: relational/prose, no closed-form identity to assert; the packet pieces it consumes are all verified) |
| (7) minimal packet `P_ij^mix` | — (definitional bundling of already-verified pieces) | not script-checked (acceptable: a tuple of verified quantities) |

Deliverables (6) and (7) are prose/relational consolidations whose constituent closed forms are all individually verified; they carry no new algebraic identity. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 55-58 | `expect_zero(K_ij + (k_i+r·k_j)/sqrt) ` | (2) slope law | yes |
| A2 | sympy | 63-66 | derivative law == 0 | (2) | yes |
| A3 | sympy | 74-76 | stationarity, slope=sqrt(k_i²+k_j²), gain | (2) | yes |
| A4 | sympy | 96-97 | curvature decomp; diagonal neutrality | (3) | yes |
| A5 | sympy | 102-104 | cross-weight deriv, stationarity r=1, value=1 | (3) | yes |
| A6 | sympy | 122-127 | grad ray dir/K/curvature | (5) | yes |
| A7 | sympy | 140-142 | equal-mix dir/slope/curvature | (5) | yes |
| A8 | sympy | 159-166 | envelope weighted forms | (4) | yes |
| A9 | sympy | 188-195 | quadratic root relation tau_lo/hi | (4) | yes |
| A10 | sympy | 202-217 | coincidence + grad cross-weight | (5) | yes |
| M1 | wl | 91-104 | oriented/positive slope + derivative law | (2) | yes |
| M2 | wl | 106-135 | Solve/Reduce stationary root, count==1, region==kj/ki, slope² | (2) | yes (DERIVES r_grad) |
| M3 | wl | 139-164 | curvature decomp via Coefficient/Series, neutrality, cross-weight | (3) | yes (DERIVES weights) |
| M4 | wl | 168-204 | canonical rays, coincidence | (5) | yes |
| M5 | wl | 211-225 | envelope Rayleigh + weighted forms | (4) | yes |
| M6 | wl | 230-245 | rootMap closure quadratic | (4) | yes |

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.txt:11` (mtime 2026-05-11) vs `scripts/...stage208..._sympy_audit.py` (mtime 2026-06-03)

**What's wrong:**
The committed SymPy `.txt` output banner reads `STAGE 191 — PAIRWISE MIXED RAYS...` (lines 11 and 206), but the current `.py` source banner reads `STAGE 208 — ...` (line 35 / 219). Git confirms the `.py` banner was renumbered 191→208 on 2026-06-03 (commit e2a4780), while the saved `.txt` predates that (last touched 2026-05-11). This is the expected P4-52 stale-banner artifact. The mathematical content of the captured run is otherwise faithful to the current source (all slope/weight/curvature/bracket forms in the `.txt` match the current `.py`); only the stage label is stale.

**Why this matters:**
Cosmetic only — the captured residuals/identities are unchanged; the divergence is purely the stage-number banner. The Mathematica `.txt` (mtime 2026-06-02) already carries the correct `STAGE 208` banner.

**Required change:**
None mandatory. A fresh `python3` re-run (which the verifier triggers) will regenerate the `.txt` with the correct `STAGE 208` banner. No script edit needed.

**Verification:**
After re-run, `scripts/output/...stage208..._sympy_audit.txt` line ~11 should read `STAGE 208 — PAIRWISE MIXED RAYS AND OFF-DIAGONAL HESSIAN SYNERGY`.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.** The `.wl` is NOT a transliteration of the `.py`. Two load-bearing objects are extracted by *genuinely different methods*:

1. **Gradient-optimal ratio `r_grad` (load-bearing for deliverable 2).**
   - SymPy *posits* the answer and checks stationarity at it:
     ```
     r_grad = sp.simplify(kj / ki)                          # py:68 — posited
     expect_zero("gradient-optimal ratio stationarity", kij_prime.subs(r, r_grad))  # py:74
     ```
   - Mathematica *solves* for it from the stationarity equation and proves uniqueness:
     ```
     stationarySolutions = Solve[{D[positiveSlope, r] == 0, r >= 0}, r, Reals]   # wl:106
     stationaryRegion = Reduce[D[positiveSlope, r] == 0 && r >= 0, r, Reals]      # wl:112
     expectZero["M3 unique stationary root count", Length[stationarySolutions] - 1]  # wl:126
     expectTrue["M3 stationarity region equals r == kj/ki", Equivalent[stationaryRegionOnDomain, r == kj/ki]]  # wl:127
     ```
   Derive (Solve/Reduce, with a uniqueness/count check absent on the SymPy side) vs posit-and-substitute. Different operation extracting the same object.

2. **Curvature weights `(w_i, w_x, w_j)` (load-bearing for deliverable 3).**
   - SymPy *posits* the weights and checks the decomposition:
     ```
     w_x = sp.simplify(2 * r / (1 + r**2))                  # py:89 — posited literal
     H1_expected = sp.simplify(w_i * hii + w_x * hij + w_j * hjj)   # py:90
     expect_zero("mixed curvature decomposition", H1ij - H1_expected)  # py:96
     ```
   - Mathematica *extracts* the weights from the polynomial numerator by coefficient projection (using `Series` to linearize the cross term in `hij`):
     ```
     curvatureNumerator = Expand[Together[mixedCurvature*(1 + r^2)]]            # wl:141
     linearHijNumerator = Normal[Series[curvatureNumerator, {hij, 0, 1}]]      # wl:142
     weightX = cleanScalar[Coefficient[linearHijNumerator, hij]/(1 + r^2)]      # wl:145
     expectZero["M4 coefficient-recovered cross weight", weightX - 2*r/(1 + r^2)]  # wl:156
     ```
   Derive-by-coefficient-extraction vs posit-and-check. Different operation.

Shared structure (the Rayleigh quotient `s.H.s`, the gradient `{-ki,-kj}`, the monomial definitions of the cone) is the legitimately-shared *physical premise*, which the policy permits. The METHOD that extracts the load-bearing objects (r_grad, the cross-weight) differs: SymPy posits, Mathematica solves/coefficient-extracts. The "each CAS runs its own simplifier" defense is not needed here — the discriminator is a real method difference, not a re-simplification of the same algebra. Some downstream checks (M5 cross-weight derivative, M6 rootMap closure) ARE the same operation as their SymPy counterparts, but the gate-defining objects of the two main theorems are independently derived, so the script as a whole clears the bar.

## Engine cross-check

Both engines agree. Side-by-side of the load-bearing finals (from the saved `.txt` outputs):
- slope `k_ij(r) = (ki + kj·r)/sqrt(1+r²)` — sympy txt:36, wl txt:11 — agree.
- derivative `(kj - ki·r)/(1+r²)^{3/2}` — sympy txt:43, wl M2 residual 0 — agree.
- `r_grad = kj/ki`, max `sqrt(ki²+kj²)` — sympy txt:50/54, wl txt:18/25-28 — agree.
- weights `{1/(1+r²), 2r/(1+r²), r²/(1+r²)}` — sympy txt:70-74, wl txt:34 — agree.
- cross-weight at r=1 equals 1 — both pass.
- equal-mix curvature `(hii+2hij+hjj)/2`, gradient curvature `(ki²hii+2kikj·hij+kj²hjj)/(ki²+kj²)` — sympy txt:107/130, wl M6 residuals 0 — agree.
- certified bracket `tau = 2H0·sqrt(1+r²)/(ki+kj·r+sqrt(...))` quadratic residual 0 — sympy txt:189-190, wl txt:87-89 — agree.
All residuals are exactly 0 in both engines; no engine_disagreement.

## Verdict justification

`findings` (one low-severity, informational `stale_output`). Attacks tried that failed: (a) hunted for tautology — every `expect_zero` compares an independently-built matrix/Rayleigh object against the paper's closed form, so none are construction-guaranteed; (b) probed the `.wl` for transliteration — it survives because r_grad is *solved* (Solve/Reduce + uniqueness) and the curvature weights are *coefficient-extracted* (Series/Coefficient), not posited as in the `.py`; (c) checked symbol domains — `ki,kj>0, r≥0, H0>0` match the paper's `k_i=|Γ_i|>0`, `r≥0`, `H_0>0`; (d) checked the `sp.diff`/`D` derivatives — `kij` and `positiveSlope` genuinely depend on `r`, so the derivatives are nonzero and stationarity is substantive (no zero-derivative trap); (e) no integrals → no parity trap. Paper↔script alignment is exact on all closed-form deliverables (1)-(5); deliverables (6)-(7) are prose/relational consolidations of already-verified pieces and legitimately carry no new identity. The only blemish is the stale SymPy `.txt` banner (191→208), which is the known P4-52 artifact and regenerates correctly on re-run.

## Value Reconciliation (pass-2 augmentation)

All Stage-208 deliverables are symbolic closed forms (no pinned numeric constants — the only literals are structural: the `1`, `2`, `1/2`, `sqrt(2)` inside the formulas, and r=1 as the cross-weight maximizer). Reconciling each emitted symbolic deliverable against the `.tex` card / appendix / `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| cone `s_ij(r)=(1,r)/sqrt(1+r²)` | py:44, wl:78; sympy txt:17, wl txt:9 | notes ln 84-89; appendix eq `pairwise-cone` (ln 839) | MATCH |
| slope `k_ij(r)=(k_i+r·k_j)/sqrt(1+r²)` | py:47, wl:97; sympy txt:36, wl txt:11 | notes ln 117-120; appendix eq `pairwise-slope` (ln 847) | MATCH |
| derivative `dk/dr=(k_j-k_i·r)/(1+r²)^{3/2}` | py:65, wl:103; sympy txt:43 | notes ln 124-128 | MATCH |
| `r_grad=k_j/k_i` | py:68, wl:131; sympy txt:50, wl txt:18 | notes ln 132-134, ln 328; appendix eq `pairwise-gradient-ratio` (ln 852) | MATCH |
| `k_grad=sqrt(k_i²+k_j²)` | py:75, wl:134; sympy txt:54 | notes ln 138-140, ln 341; appendix ln 854 | MATCH |
| curvature `H_1,ij=(h_ii+2r·h_ij+r²·h_jj)/(1+r²)` | py:85, wl:139; sympy txt:64, wl txt:33 | notes ln 178-182; appendix eq `pairwise-curvature` (ln 859) | MATCH |
| cross-weight `w_x(r)=2r/(1+r²)`, max=1 at r=1 | py:89/104, wl:156/164; sympy txt:73, wl txt:34 | notes ln 197, ln 209-211; appendix eq `pairwise-cross-weight` (ln 866) | MATCH |
| `dw_x/dr=2(1-r²)/(1+r²)²` | py:102, wl:161; sympy txt:78 | notes ln 204 | MATCH |
| `s_grad=(k_i·e_i+k_j·e_j)/sqrt(k_i²+k_j²)` | py:112, wl:169; sympy txt:90, wl txt:51 | notes ln 332-336; appendix eq `pairwise-canonical-screens` (ln 871-872) | MATCH |
| grad curvature `(k_i²h_ii+2k_ik_j·h_ij+k_j²h_jj)/(k_i²+k_j²)` | py:126, wl:191; sympy txt:107 | notes ln 346-350 | MATCH |
| `s_eq=(e_i+e_j)/sqrt(2)`, `k_eq=(k_i+k_j)/sqrt(2)` | py:130/141, wl:193/194; sympy txt:117/126, wl txt:52 | notes ln 359-367, ln 374; appendix eq `pairwise-canonical-screens` (ln 874-875) | MATCH |
| equal-mix curvature `(h_ii+2h_ij+h_jj)/2` | py:142, wl:198; sympy txt:130 | notes ln 372-375 | MATCH |
| coincidence iff `k_i=k_j` | py:203, wl:200; sympy txt:195 | notes ln 382-385 | MATCH |
| envelopes `kappa_lo/hi=(h*_lo/hi numerator)/(1+r²)` | py:152-153, wl:213-214; sympy txt:143-154 | notes ln 263-276 | MATCH |
| certified bracket `T(H_0,k;c)=2H_0/(k+sqrt(k²-2cH_0))` | py:177, wl:229; sympy txt:173, wl txt:85 | notes ln 281-285 | MATCH |
| grad-ray cross-weight `2k_ik_j/(k_i²+k_j²)` | py:216; sympy txt:197 | notes implicit (ln 229, ln 388-390 tradeoff); not explicitly boxed | MATCH (covered by tradeoff discussion / derivable) |

INTERNAL scaffolding (no prose expectation): `expect_zero`/`expectZero` residual flags, `pass`/`fail` markers, `cleanScalar`/`cleanTensor`/`stripCE`/`zeroTensorQ` helpers, `Solve`/`Reduce` intermediate region strings, `Series` linearization intermediate, banner/subbanner labels.

reconciliation: complete; 16 deliverable values checked, 0 misaligned.

## Self-test notes

Variable-independence: confirmed `kij` and `positiveSlope` depend on `r` (so `sp.diff`/`D` w.r.t. `r` are genuinely nonzero — no zero-derivative silent-pass); the `D[positiveSlope, r]==0` Solve/Reduce in the `.wl` therefore probes a real stationarity. Parity: no integrals in this unit, trap N/A. Trivial-case: each `expect_zero` is a closed-form residual identity that reduces to 0 by construction-vs-paper-form comparison, not by domain assumption collapse — e.g. substituting `hij=0` in the neutrality check genuinely drops the cross term. The single finding (stale_output) is informational and needs no directive — it regenerates on the verifier's re-run; no `## Resolve before fix_loop` and no Codex-applied script edit is warranted, so no directive file is written.
