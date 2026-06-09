---
unit_id: 216
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md]
  paper_appendix: present
---

# Audit unit 216 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_216.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (read sec. starting line 1093 plus the row at line 63)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.txt`

## What the paper claims

Stage 216 is a structural/algebraic gate, not a numeric one. The card's `\stagefield{Output}` reads verbatim: "Unique five-coordinate simplex, exact face reduction to support-\(\le4\), support-cardinality-five gate, and proof that no higher-support local search exists." The notes enumerate seven deliverables, of which the script-checkable algebraic ones are: (1) the unique positive spherical five-simplex `Delta_5^+` and its five codimension-one faces being the Stage-215 quadruple simplices; (2) the gradient-synergy theorem — gradient-optimal point `a_5^grad = (k_λ,k_c,k_γ,k_U,k_W)/‖k‖`, max slope `‖k‖`, and per-face first-order gap `(k_5^grad)^2 − (k_face^grad)^2 = k_axis^2 > 0`; (3) the ten-way cross-leverage law `w_Σ = 2Σ_{p<q}a_p a_q = (Σa_p)^2 − 1 ≤ 4` with equality at the equal-mix barycenter `a_5^eq=(1,1,1,1,1)/√5` (giving 4), and the Cauchy-slack identity `5Σa_p^2 − (Σa_p)^2 = Σ_{p<q}(a_p−a_q)^2`; (4) the fixed-simplex certified bracket `τ = 2H0/(k+√(k^2−2H0κ))`; (5) the seven-row canonical screen packet (5 imported faces + grad + eq); and (6) the support-cardinality-5 ceiling. The notes also report the leverage ladder 4/3/2/1 for five/quad/triple/pair barycenters. The numeric Bézout budget (162/324/750/1500) belongs to Stages 217/218, not 216.

## What the script claims to verify

The SymPy script verifies, mostly by posit-and-verify: M1/M2 that the posited `a_grad=k/‖k‖` is a unit vector and gives slope `‖k‖`; M3 that each face gap equals the missing `k_axis^2` and that `k_axis^2` is positive; M4 the cross-leverage identity and Cauchy-slack identity (both as expanded-polynomial equalities); M5 that `w_Σ(a_eq)=4` and prints the 3/2/1 ladder as imported constants; M6 that the posited bracket form satisfies the certified quadratic `½κτ^2 − kτ + H0 = 0`. The Mathematica script verifies the same six families but DERIVES the load-bearing objects independently: M1/M2 by solving the Lagrange stationarity system and selecting the positive-μ branch; M4/M5 by encoding leverage and slack as matrix quadratic forms `aᵀ(J−I)a` / `aᵀ(nI−J)a` and additionally proving the `≤4` bound via `Eigenvalues[J−I]={4,−1,−1,−1,−1}`; M6 by `Solve`-ing the quadratic and selecting the smaller positive root.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Unique five-simplex + 5 faces = Stage-215 quads | py/wl narrative only (faces named, count=5) | partial (combinatorial, not algebraically asserted — but it is a definitional statement, acceptable) |
| Gradient-optimal point + max slope ‖k‖ | py M1/M2 (posit-verify); wl M1/M2 (Lagrange-Solve) | match |
| Per-face first-order gap = k_axis^2 > 0 | py M3; wl M3 | match |
| Cross-leverage `w_Σ=(Σa)^2−1≤4`; Cauchy slack | py M4; wl M4 (+ eigenvalue bound in M5) | match |
| Barycenter leverage = 4 | py M5; wl M5 (+ spectral max) | match |
| Leverage ladder 4/3/2/1 | py prints 4/3/2/1 (3/2/1 are imported literals) | match (vs notes L271/276–278) |
| Fixed-simplex certified bracket τ form | py M6; wl M6 | match |
| Seven-row canonical screen packet | py narrative (count=7) | match (structural count) |
| Support-cardinality-5 ceiling | py narrative (dim=5) | match (structural) |
| Card states "Mathematica audit: none yet" | a passing `.wl` now exists | MISMATCH (card metadata stale) → F1 |

Set `paper_alignment: partial` — the mathematics is fully aligned; the only discrepancy is the stale `Verification` metadata line in the card.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero(grad_norm_sq − 1)` | gradient point unit-norm | yes |
| A2 | sympy | 40 | `expect_zero(k_grad − k_norm)` | max slope = ‖k‖ | yes |
| A3 | sympy | 56–57 | `expect_zero(diff − k_axis^2)`, `require k_axis^2 positive` | per-face gap | yes |
| A4 | sympy | 87–88 | `expect_zero` cross-leverage & Cauchy-slack identities | leverage law | yes |
| A5 | sympy | 96 | `expect_zero(w_eq − 4)` | barycenter leverage | yes |
| A6 | sympy | 113 | `expect_zero(½κτ^2 − kτ + H0)` | certified bracket | yes |
| M1 | math | 87 | `expectZero(lagrangePoint·lagrangePoint − 1)` | gradient point unit-norm | yes |
| M2 | math | 88,91 | `expectZero(lagrangeValue − √S)`, dominance | max slope = ‖k‖ | yes |
| M3 | math | 103–104 | `expectZero(S − faceSum − k^2)`, positivity | per-face gap | yes |
| M4 | math | 117–124 | `expectZero` matrix-form leverage & Laplacian slack | leverage law | yes |
| M5 | math | 131–133 | `expectZero` barycenter=4 AND `Max[eig]−4` | barycenter max leverage | yes |
| M6 | math | 153–156 | `expectZero` solved root = stated form, smaller positive | certified bracket | yes |

No tautological rows: in each posit-verify pair the assertion CAN fail (e.g. if the posited `a_grad` were wrong, `k_grad − k_norm` would be nonzero; the residual `½κτ^2 − kτ + H0` is not `0=0` by construction). The leverage identities expand to genuine polynomial cancellations.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (stale verification metadata)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_216.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl` (entire, present + passing)

**What's wrong:**
The card's `\stagefield{Verification}` line states: "SymPy audit: \StageFile{...stage216..._sympy_audit.py}.  Mathematica audit: none yet." But a Mathematica audit `.wl` now exists, is dated 2026-06-02, and passes all six check families (its saved output ends "All Stage 216 Mathematica audit checks passed."). The sibling stage 218 card already cites its `.wl` in this same field, so the field is meant to name the Mathematica audit when one exists. This is a documentation-metadata staleness from the pass-1 dual-engine retrofit, not a math disagreement.

**Why this matters:**
The card under-reports the verification coverage of the stage. Left alone, a reader (or a downstream coverage tracker) believes stage 216 is single-engine when it is in fact dual-engine. Direction (update the card vs. consider the `.wl` out of scope) is a user/paper-owner decision; Codex does not edit paper/.

**Required change:**
See `## Resolve before fix_loop` in the directive. No Codex-applied change pending user resolution.

**Verification:**
After user resolution, if direction (a): card line 11 names the `.wl` path exactly as stage 218 does.

## Independent-derivation check (Mathematica)

INDEPENDENT. This is the key question for a retrofit-heavy batch where sibling 211 was found to be a port. The `.wl` is NOT a transliteration — it uses a structurally different extraction method for every load-bearing object:

1. **Gradient optimum (M1/M2).** `.py` POSITS the closed form and verifies it: `a_grad = [sp.simplify(k / k_norm) for k in ks]` then `expect_zero("M2 ...", k_grad - k_norm)` (py:36,40). `.wl` DERIVES it by solving the constrained-optimization Lagrange system:
   ```
   lagrangian = kVec.aVec - mu (aVec.aVec - 1);
   stationarity = Thread[(D[lagrangian, #] & /@ aVec) == 0];
   lagrangeRules = Solve[Join[stationarity, {aVec.aVec == 1}], Append[aVec, mu], Reals];
   positiveRule = SelectFirst[lagrangeRules, TrueQ[FullSimplify[(mu /. #) > 0, ...]] &];
   ```
   (wl:72–82). Derive-via-Lagrange-Solve vs posit-the-answer = independent.

2. **Cross-leverage / Cauchy slack (M4).** `.py` writes the monomial sum explicitly: `w_sigma = 2*(a1*a2 + a1*a3 + ... + a4*a5)` and `pair_sum = sum((x-y)**2 ...)` then expands (py:80–86). `.wl` encodes them as matrix quadratic forms: `offDiagonalForm = ConstantArray[1,{n,n}] - IdentityMatrix[n]; wSigma = aVec.offDiagonalForm.aVec` and `laplacianForm = n IdentityMatrix[n] - Outer[Times,ones,ones]` (wl:111–124). Different construction (matrix vs explicit-monomial) of the same object.

3. **Bound `w_Σ ≤ 4` / bracket root (M5/M6).** `.py` substitutes one point (`w_sigma.subs(a_i → 1/√5)`, py:95) and posits the bracket then verifies the quadratic (py:112–113). `.wl` proves the bound spectrally via `Eigenvalues[offDiagonalForm]` and `Max[leverageEigenvalues] - 4` (wl:129,133) — a strictly stronger and different argument — and solves the quadratic with `Solve[quadraticResidual == 0, tau, Reals]` selecting the smaller positive root (wl:138–156).

In every case the METHOD that extracts the load-bearing object DIFFERS (derive-vs-posit, matrix-vs-monomial, spectral-vs-point). The "each CAS runs its own simplifier" defense is not needed here — the routes are genuinely different at the construction level. M3 (face gap = missing square) is the same trivial subtraction in both, but it is a corollary, not a load-bearing object, and sharing the monomial premise is permitted.

## Engine cross-check

Both engines agree. Common results: gradient point `k_i/√(Σk^2)`, max slope `√S`, face gaps `k_axis^2`, leverage identity `=0`, slack identity `=0`, barycenter leverage `4` (wl additionally: top eigenvalue `4`), bracket `2H0/(k+√(k^2−2H0κ))` with quadratic residual `0`. SymPy output (txt) and Mathematica output (txt) report the same forms; all checks PASS in both. No residual/sign/factor disagreement.

## Verdict justification

The mathematics holds up under attack. Attacks tried and failed: (i) tautology hunt — every posit-verify assertion can fail and the residual is not `0=0` by construction; (ii) transliteration attack (the batch's headline risk) — the `.wl` derives via Lagrange-Solve, matrix quadratic forms, eigenvalues, and quadratic Solve, none of which mirror the `.py`'s posit-and-verify choreography, so it is genuinely independent; (iii) domain/assumption attack — `k_i>0`, `H0,κ>0`, and `k^2>2H0κ` in the `.wl` `$Assumptions` exactly match the script's physical setup and the notes' positivity premises, so the Lagrange positive-branch selection and the real-distinct-positive-roots of M6 are justified; (iv) value reconciliation — all emitted deliverables (symbolic forms + the 4/3/2/1 ladder) reconcile with the notes and appendix. I read the paper card, the notes, and the appendix; the verified claims match. The single finding is a low-severity stale-metadata `paper_misalignment` (card says "Mathematica audit: none yet" while a passing `.wl` exists), which routes to the user.

## Value Reconciliation (pass-2 augmentation)

Stage 216 emits no free numeric constants of its own; its deliverables are symbolic. The only literals are structural integers (the leverage ladder and the support ceiling). The Bézout budget numbers (162/324/750/1500) belong to Stages 217/218 and are NOT emitted by this stage's scripts. Reconciliation against notes + appendix:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `a_5^grad = (k_λ,k_c,k_γ,k_U,k_W)/√(Σk^2)` | py M1 (out L9–10); wl M1 (out L9) | notes L161–166; appendix eq:app-part06-five-grad (L1138–1143) | MATCH |
| max slope `k_5 = √(Σk^2)` | py M2 (out L10); wl M2 (out L13) | notes L168–174 | MATCH |
| per-face gap `= k_axis^2` (>0) | py M3 (out L11–15); wl M3 (out L21–40) | notes L189–195 | MATCH |
| cross-leverage `w_Σ=(Σa)^2−1` | py M4 (out L24); wl M4 (out L45) | notes L240–244; appendix eq L1150–1152 | MATCH |
| Cauchy-slack identity `=Σ(a_p−a_q)^2` | py M4 (out L25); wl M4 (out L47) | notes L246–251 | MATCH |
| barycenter leverage `w_Σ(a_eq)=4` | py M5 (out L27); wl M5 (out L56) | notes L269–273; appendix L1152–1154 | MATCH |
| barycenter `a_5^eq=(1,1,1,1,1)/√5` | py M5 (out L26) | notes L262–266; appendix eq L1144–1148 | MATCH |
| leverage ladder quad=3 | py (out L28) imported literal | notes L276 | MATCH |
| leverage ladder triple=2 | py (out L29) imported literal | notes L277 | MATCH |
| leverage ladder pair=1 | py (out L30) imported literal | notes L278 | MATCH |
| certified bracket `τ=2H0/(k+√(k^2−2H0κ))` | py M6 (out L33); wl M6 (out L64) | notes L316–322; appendix eq:app-part06-five-root (L1166–1173) | MATCH |
| top leverage eigenvalue `=4` | wl M5 (out L53,58) | — (internal spectral proof of `w_Σ≤4`, supports notes L256) | MATCH (supports the `≤4` bound; not a standalone deliverable) |
| canonical screen rows `=7` | py (out L39) | notes L353–367 (5 faces + grad + eq); appendix L1129–1154 | MATCH |
| support ceiling `5` | py (out L43) | notes L471–481; appendix eq:app-part06-support-ceiling (L1222–1226) | MATCH |

INTERNAL (scaffolding, no finding): unit-norm residuals (`||a_grad||^2=1`), `Lagrange μ>0` branch selector, interior ratio prints `r/s/t/u = k_c/k_λ …` (these are the grad ratios, also stated notes L177–185 — MATCH if counted), Solve root list, eigenvalue list, pass/fail flags.

reconciliation: complete; 14 deliverable values checked, 0 misaligned (the lone finding is the stale `.tex` Verification-field metadata, F1, not a value mismatch).

## Self-test notes

Checked: (1) variable-independence — no zero-derivative trap; the only `D[...]` is the Lagrangian gradient `D[kVec.aVec − μ(a·a−1), a_i]`, and each `a_i` genuinely appears, so stationarity is non-degenerate. (2) Symmetry/parity — no unbounded integrals in this stage. (3) Trivial-case — substituting `a_i=1/√5` gives `w_Σ=2·10·(1/5)=4` ✓ and `Σa_i^2=1` ✓; the `J−I` top eigenvalue is `n−1=4` ✓; the bracket residual at `τ=2H0/(k+√(k^2−2H0κ))` cancels by rationalization ✓. No directive math is prescribed (the sole finding is a user-gated paper_misalignment), so the round-trip trap is N/A.
