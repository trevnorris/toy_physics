---
unit_id: 026
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md"]
  paper_appendix: present
---

# Audit unit 026 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_026.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 026 read; longtable line 42 and closure item 5)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.txt`

## What the paper claims

The stage's `\stagefield{Output}` states verbatim: "Stage~026 outputs the finite-throat D/N overlap \eqref{eq:app-stage026-kappa0}, the concrete couplings \eqref{eq:app-stage026-couplings}, the finite-throat normalization law \eqref{eq:app-stage026-normalization}, and the required stiffness formula \eqref{eq:app-stage026-Kreq}."

The four deliverables decompose as: (i) the overlap law `kappa_n = sqrt(2)/((n+1/2)*pi)` with the lowest-branch value `kappa = 2*sqrt(2)/pi`; (ii) the branch couplings `C = kappa*lambda_B, G_U = lambda_U, G_W = kappa*lambda_W, R = kappa*lambda_R`, plus the Stage-25 reductions `Delta = Omega_U^2*Omega_W^2 - kappa^2*lambda_R^2`, `Q = lambda_U^2*Omega_W^2 + 2*kappa^2*lambda_U*lambda_W*lambda_R + kappa^2*lambda_W^2*Omega_U^2`, `P = kappa*(Omega_U^2*lambda_W + lambda_R*lambda_U)`, and `X = kappa^2*lambda_B^2/varpi^2`; (iii) the normalization law `mhat_rad^2 * kappa^2*(Omega_U^2*lambda_W + lambda_R*lambda_U)^2 / [Delta*(K*Delta - Delta*kappa^2*lambda_B^2/varpi^2 - Q)] = N_Q^target`; (iv) the closed-form required stiffness `K_req = kappa^2*lambda_B^2/varpi^2 + Q/Delta + mhat_rad^2*kappa^2*(Omega_U^2*lambda_W + lambda_R*lambda_U)^2 / (N_Q^target * Delta^2)`. The `\stagefield{Checks}` block also itemizes "Solving (eq:app-stage026-normalization) for K gives (eq:app-stage026-Kreq)." The constant-wall geometric identity `K = K_eta + 6*T_Omega` (eq:app-stage026-Kgeo) is a carry-forward, not a new derivation. The notes file mirrors this content and identifies `N_Q^target = 54 G c_s^5 / (5 a^5 c^5)`.

## What the script claims to verify

The SymPy and Mathematica scripts each verify: (a) orthonormality of `u0(s) = 1/sqrt(L)` and `f0(s) = sqrt(2/L)*sin(pi*s/(2*L))`; (b) the general overlap law `kappa_n` (matched against the paper's closed form) and the lowest-branch value `kappa = 2*sqrt(2)/pi`; (c) the explicit branch overlap identifications `I_(eta,phi)=I_(eta,w)=I_(u,w)=kappa`, `I_(eta,u)=1`, with the duplicated identifications explained by a comment in the SymPy script noting that eta-with-u0 and (phi,w,u-leg)-with-f0 collapse three of the four overlaps to the same integral; (d) the printed symbolic forms of `C, G_U, G_W, R, Delta, Q, P, B0, Z0, N0, D0, P0` after substitution; (e) that solving `mhat^2 * N0 / D0 - target == 0` for K returns a `K_req` whose back-substitution into the residual yields zero. The Mathematica script also provides an independent path for `kappa_n` (indefinite integration + boundary substitution + analytic-short-form cross-check).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `kappa_n = sqrt(2)/((n+1/2)*pi)` (eq:app-stage026-kappa-n) | sympy A3: `kappa_n - kappa_n_expected == 0`; mathematica A3: indef-integral path vs short form `== 0` | match |
| `kappa = 2*sqrt(2)/pi` (eq:app-stage026-kappa0) | sympy A4 / mathematica A4: `kappa - 2*sqrt(2)/pi == 0` | match |
| `I_(eta,phi)=kappa, I_(eta,u)=1, I_(eta,w)=kappa, I_(u,w)=kappa` (eq:app-stage026-overlaps) | sympy A5: `overlap_u0_f0 - kappa == 0` and `overlap_u0_u0 - 1 == 0`; equivalent in mathematica | match (comment 124-128 documents the three-way identification, no false duplication of independent checks) |
| `C, G_U, G_W, R` couplings (eq:app-stage026-couplings) | constructed in `branch_substitution()`; output strings reduce to the paper forms after substituting `kappa = 2*sqrt(2)/pi` | match (constructive) |
| `Delta, Q, P` (eq:app-stage026-DeltaQP) | constructed; output strings match the paper form after substituting `kappa^2 = 8/pi^2` | match (constructive) |
| `X = kappa^2*lambda_B^2/varpi^2` (eq:app-stage026-X) | `B0 = C^2/varpi^2` constructed; output `8*lambda_B^2/(pi^2*varpi^2)` equals `kappa^2*lambda_B^2/varpi^2` | match (constructive) |
| Normalization law (eq:app-stage026-normalization) | residual `mhat^2*N0/D0 - target` displayed; algebraic identity `Delta*(K*Delta - Delta*B0 - Q) = Delta^2*(K - B0 - Q/Delta)` makes the script's `mhat^2 N0/D0` the same scalar as the paper's `mhat^2*kappa^2*(...)^2 / [Delta*(K*Delta - Delta*B0 - Q)]` | match |
| **`K_req` closed-form** (eq:app-stage026-Kreq) | sympy A6 / mathematica A6: `residual @ K_req == 0` only | partial — back-substitution is algebraically guaranteed by Solve and does not verify the paper's three-term decomposition `B0 + Q/Delta + mhat^2*kappa^2*(...)^2/(N_Q^target * Delta^2)` |
| `K = K_eta + 6*T_Omega` (eq:app-stage026-Kgeo) | printed only, no assertion | defensible — explicitly a carry-forward, not a new derivation |
| Branch condition `K_eta + 6*T_Omega = K_req` (eq:app-stage026-branch-condition) | not checked; physical condition, not algebraic | n/a |

Dominant pattern is `match`. The single `partial` is a verification-rigor concern, not a paper↔script disagreement. Setting `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 83 | `expect_zero("int u0^2 - 1", ...)` | normalization of nu_0 (eq:app-stage026-u0) | yes |
| A2 | sympy | 84 | `expect_zero("int f0^2 - 1", ...)` | normalization of f_0 (eq:app-stage026-DN-ladder) | yes |
| A3 | sympy | 105 | `expect_zero("kappa_n - expected", kappa_n - kappa_n_expected)` | overlap law (eq:app-stage026-kappa-n) | yes |
| A4 | sympy | 109 | `expect_zero("kappa - 2*sqrt(2)/pi", kappa - kappa_expected)` | minimal branch overlap (eq:app-stage026-kappa0) | yes |
| A5 | sympy | 145-146 | `expect_zero("overlap_u0_f0 - kappa", ...)`, `expect_zero("overlap_u0_u0 - 1", ...)` | overlap identifications (eq:app-stage026-overlaps) | yes |
| A6 | sympy | 202 | `expect_zero("residual @ K_req", residual.subs(K, K_req))` | normalization + K_req (eqs 026-normalization, 026-Kreq) | partial — verifies Solve gave a root, not that the root matches the paper's three-term form |
| M1 | mathematica | 48 | `expectZero["int u0^2 - 1", ...]` | normalization of nu_0 | yes |
| M2 | mathematica | 49 | `expectZero["int f0^2 - 1", ...]` | normalization of f_0 | yes |
| M3 | mathematica | 82 | `expectZero["kappa_n (analytic) - (fundamental thm)", kappaN - kappaNViaFundamentalThm]` | overlap law via independent path | yes |
| M4 | mathematica | 88 | `expectZero["kappa - 2*sqrt(2)/pi", kappa - kappaExpected]` | minimal branch overlap | yes |
| M5 | mathematica | 113-114 | `expectZero["overlap_u0_f0 - kappa", ...]`, `expectZero["overlap_u0_u0 - 1", ...]` | overlap identifications | yes |
| M6 | mathematica | 162 | `expectZero["residual @ K_req", residual /. k -> kReq]` | normalization + K_req | partial — same Solve-tautology issue as A6 |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:195-202`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:155-162`

**What's wrong:**
The paper card lists eq:app-stage026-Kreq as one of the four boxed `Output` deliverables and the `\stagefield{Checks}` block explicitly says "Solving (eq:app-stage026-normalization) for K gives (eq:app-stage026-Kreq)." The paper's closed form is

```
K_req = kappa^2*lambda_B^2/varpi^2
      + Q/Delta
      + mhat_rad^2 * kappa^2 * (Omega_U^2*lambda_W + lambda_R*lambda_U)^2
        / ( N_Q^target * Delta^2 ).
```

The script's verification of this deliverable, in full, is:

SymPy lines 197-202:
```python
K_req = sp.solve(sp.Eq(residual, 0), K)[0]
...
expect_zero("residual @ K_req", residual.subs(K, K_req))
```

Mathematica lines 157-162:
```mathematica
kReq = k /. First[Solve[residual == 0, k]];
kReq = FullSimplify[kReq, Assumptions -> $Assumptions];
...
expectZero["residual @ K_req", residual /. k -> kReq];
```

Substituting back into the equation you just solved is algebraically guaranteed by `sp.solve` (resp. `Solve`); the assertion cannot fail unless the CAS is broken. The check does not verify that `K_req` has the structure the paper claims (three terms `B0 + Q/Delta + mhat^2*kappa^2*(...)^2/(N_Q^target * Delta^2)`). A K_req formula that differed from the paper but still solved the equation would pass the current check silently.

**Why this matters:**
The K_req closed form is what downstream stages quote: stage 032 uses the natural D/N source map `mhat_-^2 = s_-/kappa_0^2` to convert the K_req-style target into the selected-branch target, and the Part-II appendix closure (line 124 of `stage_appendix_part02.tex`) lists this branch's `kappa` and its admissibility tests together. Without an independent structural check, the load-bearing form is currently dependent on `Solve`'s output formatting rather than a derived identity. The Mathematica script already FullSimplify's into the three-term form (output line 79), but no assertion enforces that, so a future refactor could regress silently.

**Required change:**
After computing K_req via Solve in both scripts, add an explicit structural check against the paper's closed form.

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`, immediately after line 197 (`K_req = sp.solve(...)`) and before line 202 (`expect_zero("residual @ K_req", ...)`), insert:

```python
    K_req_paper = (
        B0
        + Q / Delta
        + mhat**2 * kappa**2 * (Omega_U**2 * lambda_W + lambda_R * lambda_U)**2
          / (target * Delta**2)
    )
    expect_zero("K_req - K_req_paper", sp.simplify(K_req - K_req_paper))
```

(Use the existing in-scope symbols `B0, Q, Delta, mhat, kappa, Omega_U, lambda_W, lambda_R, lambda_U, target`.)

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`, immediately after line 158 (`kReq = FullSimplify[kReq, ...]`) and before line 162 (`expectZero["residual @ K_req", ...]`), insert:

```mathematica
  kReqPaper = FullSimplify[
    b0
    + q/delta
    + mhat^2 * kappa^2 * (omegaU^2*lambdaW + lambdaR*lambdaU)^2
      / (target * delta^2),
    Assumptions -> $Assumptions
  ];
  expectZero["K_req - K_req_paper", kReq - kReqPaper];
```

Both new checks are non-tautological: they compare Solve's symbolic output to the paper's structural decomposition. If Solve ever returned an algebraically-equivalent-but-structurally-different form that happened to equal `K_req_paper` after simplification, the check still passes; if it returned something that does not equal `K_req_paper`, the check fails — which is the desired behavior.

**Verification:**
After Codex applies, the SymPy output file should contain a new line `K_req - K_req_paper = 0` between sections IV.2's `K_req = ...` block and the closing `residual @ K_req = 0`. The Mathematica output should contain a `PASS: K_req - K_req_paper` line in the same neighborhood. Both scripts must continue to exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally close to a transliteration of the SymPy script (same module names — `concreteModes`, `overlapLaw`, `branchSubstitution`, `normalizationTest`; same variable choreography; same assertions in the same order; same I/II/III/IV banner pattern). The substantive divergence is in Section II.1 where Mathematica derives `kappa_n` via an indefinite-integration + boundary-substitution path (`indef = Integrate[u0*fN, s]; (indef /. s -> l) - (indef /. s -> 0)`) and cross-checks it against the analytic short form `Sqrt[2]/((n+1/2)*Pi)`, whereas the SymPy script uses a single `Integrate[u0*f_n, (s, 0, L)]` and checks against the same short form. That is a legitimately different code path for the load-bearing kappa identity. The remaining sections (Section III.1 `Integrate[u0*f0, {s,0,l}]`, Section III.2 constructive substitutions, Section IV residual + Solve) follow the same approach in both engines.

This is acceptable rather than flag-worthy: Section II.1 is the only nontrivial algebraic discovery in the stage (Sections III–IV are substitutions and rearrangements), and that section is independently rederived. The structural similarity in Sections III–IV reflects that there is essentially one canonical way to substitute a constant and rearrange a linear equation; no second-engine policy violation. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce the same final symbolic expressions, modulo CAS-formatting:

| Quantity | SymPy (output line) | Mathematica (output line) |
|---|---|---|
| `kappa_n` | `2*sqrt(2)/(pi*(2*n + 1))` (26) | `(2*Sqrt[2])/(Pi + 2*n*Pi)` (28) |
| `kappa` | `2*sqrt(2)/pi` (33) | `(2*Sqrt[2])/Pi` (36) |
| `kappa` numeric | `0.900316316157106` (35) | `0.90031631615710606955...` (39) |
| `Delta` | `Omega_U**2*Omega_W**2 - 8*lambda_R**2/pi**2` (56) | `omegaU^2*omegaW^2 - (8*lambdaR^2)/Pi^2` (62) |
| `Q` | `(8*Omega_U**2*lambda_W**2 + pi**2*Omega_W**2*lambda_U**2 + 16*lambda_R*lambda_U*lambda_W)/pi**2` (57) | `lambdaU^2*omegaW^2 + (8*lambdaW*(2*lambdaR*lambdaU + lambdaW*omegaU^2))/Pi^2` (63) |
| `P` | `2*sqrt(2)*(Omega_U**2*lambda_W + lambda_R*lambda_U)/pi` (58) | `(2*Sqrt[2]*(lambdaR*lambdaU + lambdaW*omegaU^2))/Pi` (64) |
| `B0` | `8*lambda_B**2/(pi**2*varpi**2)` (59) | `(8*lambdaB^2)/(Pi^2*varpi^2)` (65) |
| `K_req` form | single rational fraction (92-119) | three-term decomposition matching paper (79) |
| `residual @ K_req` | `0` (122) | `0` + PASS (80-81) |

The SymPy K_req is presented as one big rational; expanding `K_req_sympy` and dividing the numerator's three-term structure by the common denominator `pi^2*G*cs^5*varpi^2*(pi^2*Omega_U^2*Omega_W^2 - 8*lambda_R^2)^2` reproduces the Mathematica three-term decomposition. They are equivalent symbolic expressions. Both back-substitute to residual `= 0` as the respective scripts confirm. Engines agree. No `engine_disagreement`.

Stale-output check: SymPy script mtime is May 21 16:59; output mtime is May 21 17:01. Mathematica script mtime is May 21 17:00; output mtime is May 21 17:02. Both outputs are newer than their scripts. Outputs are fresh.

## Verdict justification

Verdict is `findings` with one medium-severity `insufficient_verification`. The paper's deliverables for stage 026 — the overlap law `kappa_n`, the branch overlaps, the explicit `C, G_U, G_W, R, Delta, Q, P, X` substitutions, the normalization law residual, and the recovered K_req — are all correctly reflected in the scripts; the SymPy and Mathematica engines agree on every symbolic output. The one defect is that the paper's K_req closed-form (eq:app-stage026-Kreq), explicitly named in `\stagefield{Output}` and `\stagefield{Checks}`, is currently anchored only by a back-substitution check `residual @ K_req == 0` that is algebraically guaranteed by Solve and cannot detect a divergence between the script's K_req and the paper's three-term structural form.

Attacks tried that failed: (a) probed whether `B0 = C^2/varpi^2` actually reduces to `kappa^2*lambda_B^2/varpi^2` after substitution — yes, numerically `8*lambda_B^2/(pi^2*varpi^2) = (2*sqrt(2)/pi)^2 * lambda_B^2 / varpi^2`; (b) probed whether the Mathematica `kappa_n` path is independent — Section II.1 uses an indefinite-integral + boundary-substitution path distinct from SymPy's definite integral; (c) probed whether the residual `mhat^2*N0/D0 - target` algebraically equals the paper's normalization expression — yes, `Delta*(K*Delta - Delta*B0 - Q) = Delta^2*(K - B0 - Z0)`, so script and paper formulations are identical scalars; (d) probed for hidden symbol-assumption errors — `K` is real (allowed to be sign-indefinite), `Delta` is allowed to be sign-indefinite, lambdas real, Omegas/varpi/G/c/cs/a/mhat positive, consistent with the physical setup; (e) probed for stale outputs — both outputs are newer than their scripts. No `paper_misalignment` (paper and script agree on all derivations), no `tautological_check` finding distinct from F1 (the Solve back-substitution is more accurately classified as `insufficient_verification` because it is the only check on a paper-named deliverable, not a redundant gloss on an otherwise-checked claim), no `engine_disagreement`, no `mathematica_transliteration`, no `stale_output`. The fix is mechanical and local to two files.

## Self-test notes

Walked through the proposed `K_req_paper` expression mentally. Substituting the script's symbolic values: `B0 = 8*lambda_B^2/(pi^2*varpi^2)`, `Q/Delta = (...)/(Omega_U^2*Omega_W^2 - 8*lambda_R^2/pi^2)`, and the third term becomes `mhat^2 * (8/pi^2) * (Omega_U^2*lambda_W + lambda_R*lambda_U)^2 / (target * Delta^2)`. The Mathematica output (line 79) already prints this exact three-term form, so the proposed assertion `kReq - kReqPaper == 0` will pass with the current Mathematica K_req; for SymPy, the single-fraction Solve output will simplify against the paper form to zero (sp.simplify is sufficient because the difference is identically zero — direct algebra). The variable-independence trap (sp.diff against a non-dependency) does not apply (no derivatives proposed); parity trap (odd integrand on symmetric domain) does not apply (no new integrals proposed); trivial-case substitution: at `mhat = 0` both `K_req` and `K_req_paper` collapse to `B0 + Q/Delta`, consistent. Paper round-trip (v2): the proposed assertion uses only symbols and constants the paper already states; no new numeric constants introduced; the eq:app-stage026-Kreq form is reproduced symbol-for-symbol. Path specifications are absolute; both files already exist and are the only edit targets.
