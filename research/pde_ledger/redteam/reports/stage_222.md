---
unit_id: 222
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 222 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_222.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows: line 56 table entry; lines 518-572 narrative for Stage 222; line 679 claim status)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The `\stagefield{Output}` (stage_222.tex:15) states: "Concrete branch pole census, residue/linewidth cancellation, and first branch-level dynamic survival test. The pole polynomial and survival gate are exact in the primitive closure; the displayed slice is a numerical placement." The notes (§Purpose, lines 21-36) and appendix (lines 518-572) enumerate the deliverables: (1) the exact primitive overlap constant `kappa = 2*sqrt(2)/pi`; (2) the conservative quartic pole polynomial `F(y)` with `D(omega) = F(omega^2)/((varpi^2-omega^2)*Delta(omega))`, degree 4; (3) the residue/linewidth cancellation `R_{Q,*} = |A_{Q,*}|/gamma_* = 27 c_s^5 / (a^5 omega_*^5 N(omega_*))` where the pole derivative `D_0'` cancels exactly, with `N(omega) = P(omega)^2/Delta(omega)^2`, `P(omega) = A(omega) G_W + R G_U`, `Gamma_5 = a^5/(27 c_s^5)`; (4) the exact low-loss survival threshold `R_{Q,*} >= 2 DeltaV_req (1+eta^2)/eta * x^6`; (5) an explicit numerical pole census on the admissible sample slice `(lam_B,lam_U,lam_W,lam_R,Omega_U,Omega_W,varpi,K,M) = (0.5,0.3,0.4,0.25,1.0,1.4,2.0,3.0,1.0)`, a=c_s=1, with stated static data and four-pole table; (6) the monotone static/dynamic tension under a `lambda_W` scan (P0 increases, upper-wall R_{Q,*} decreases). The card is checkpoint:False, status-only:False; both engines are nominally required.

## What the script claims to verify

The SymPy module's docstring/prints and assertions verify, in order: the overlap integral equals `2*sqrt(2)/pi` (line 56); the hand-built quartic `F_y` reconstructs the numerator of `D = K_B - Q/Delta` exactly, `quartic_relation == 0` and `Poly(F_y,y).degree() == 4` (lines 97-101); the residue/linewidth ratio collapses with `D0prime` cancelling to `27 c_s^5/(a^5 omega_*^5 N_star)` (line 117); the symbolic low-loss peak and survival threshold are formed (lines 119-124, printed but not asserted as identities); the static slice data `Delta0, D0, N0, P0` match the quoted numerics (lines 167-170); the uncoupled and coupled pole frequencies match the quoted census (lines 189-194); the four `R_{Q,*}` figures match (lines 204-208); the two illustrative low-loss thresholds match (lines 219-220); and the `lambda_W` scan rows match a hardcoded table and are monotone (lines 259-264). The script computes the poles and R_Q figures from the symbolic `F_y`/`N_omega` (not from raw literals), so the regression literals anchor an in-script derivation.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) kappa = 2*sqrt(2)/pi | line 56 `assert simplify(kappa - 2 sqrt2/pi)==0` | match |
| (2) quartic F(y), D=F/((varpi^2-w^2)Delta), deg 4 | lines 97-101 | match |
| (3) R_{Q,*}=27c_s^5/(a^5 w_*^5 N(w_*)), derivative cancels | line 117 (cancellation) + line 197 RQ_expr w/ N_omega line 153 | match |
| (4) low-loss survival threshold | lines 119-124 formed; lines 219-220 numeric thresholds asserted | partial (symbolic forms printed, not asserted as identities; numeric thresholds asserted) |
| (5) numerical pole census (static data + 4-pole table) | lines 167-208 | match |
| (6) lambda_W scan static/dynamic tension table + monotonicity | lines 259-264 | mismatch (lambda_W=0.2 R_{Q,*}: notes 213.483858657863 vs script 145.4838586578641) |

`paper_alignment: partial` — five of six deliverables are faithfully exercised; deliverable (6)'s table carries one row whose R_{Q,*} value disagrees between notes and script.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 56 | `simplify(kappa - 2 sqrt2/pi) == 0` | claim 1 | yes |
| A2 | sympy | 100 | `cancel(D - F(w^2)/((varpi^2-w^2)Delta)) == 0` | claim 2 | yes |
| A3 | sympy | 101 | `Poly(F_y,y).degree() == 4` | claim 2 | yes |
| A4 | sympy | 117 | `simplify(RQ_star - 27 c_s^5/(a^5 w_*^5 N_star)) == 0` | claim 3 (cancellation) | yes |
| A5 | sympy | 167-170 | `assert_close(Delta0/D0/N0/P0, lit)` | claim 5 (static data) | yes |
| A6 | sympy | 189-194 | `assert_close(roots, lit)` | claim 5 (pole census) | yes |
| A7 | sympy | 208 | `assert_close(RQ(pole), lit)` | claim 5 (R_Q figures) | yes |
| A8 | sympy | 219-220 | `assert_close(threshold_eta, lit)` | claim 4 (numeric thresholds) | yes (illustrative inputs) |
| A9 | sympy | 261 | `assert_close(scan_row, lit)` | claim 6 (scan table) | partial (row 0.2 disagrees w/ notes) |
| A10 | sympy | 263-264 | monotone increasing P0 / decreasing R_Q | claim 6 (tension) | yes (insensitive to the 213-vs-145 value; both decrease) |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md:395`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.py:253`

**What's wrong:**
The notes static/dynamic tension table (§8) gives, for `lambda_W = 0.2`:
> `| 0.2 | 0.00594740531769 | 2.82723442158450 | 2.04402272302752 | 213.483858657863 |`

The script's hardcoded expected scan row and the saved output give the same `P0` (0.005947405317693074), `D0` (2.8272344215844973), and upper-wall `omega_*` (2.04402272302752), but the upper-wall `R_{Q,*}` is **145.4838586578641**, not 213.483858657863:
> script line 253: `(0.2, 0.005947405317693074, 2.8272344215844973, 2.04402272302752, 145.4838586578641),`
> output line 56: `(0.2, 0.005947405317693074, 2.8272344215844973, 2.04402272302752, 145.4838586578641)`

All other scan rows (lambda_W = 0.4, 0.6, 0.8, 1.0) agree exactly between notes and script, including lambda_W=0.4 → 32.0025481088465. Only the lambda_W=0.2 row's R_{Q,*} disagrees. Because `R_{Q,*} = 27/(omega_*^5 N(omega_*))` depends only on the pole frequency and the transfer factor, and the script recomputes it from the symbolic `RQ_expr`/`N_omega` rather than from a literal, the script-side number is a derived value, whereas the notes table number appears to be a stale or independently-typed value. The script's `assert_monotone_decreasing` (line 264) cannot catch this: both 213→32 and 145→32 are strictly decreasing, so the assertion passes regardless.

**Why this matters:**
The notes table is the paper-facing artifact for deliverable (6) (the appendix narrative, lines 518-572, only describes the tension qualitatively and does not list per-row values, so the notes table is the only place this number is published). A reader citing the notes would carry 213.48; the verified figure is 145.48. This is a published-value discrepancy in a `\StatusNumerical{}` row.

**Required change:**
None by Codex — this is a paper/notes ↔ script disagreement and routes to the user. See `## Resolve before fix_loop` in the directive.

**Verification:**
After the user picks a direction: if the script value (145.4838586578641) is correct, the notes table line 395 R_{Q,*} cell is updated to 145.483858657864 (paper-side edit, Codex-applied only after explicit follow-up authorization); if the notes value is correct, the script's symbolic computation must be re-examined (but note all four poles at lambda_W=0.4 already match notes, so a script bug is unlikely).

### F2 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `moving_throat_pde_stage222_*` file present)

**What's wrong:**
No Mathematica audit script exists for this unit. The card (`stage_222.tex:11`) explicitly states "Mathematica audit: none yet." The dual-engine rule requires a `.wl` wherever Mathematica CAN independently verify the stage, and every deliverable here is squarely within Mathematica's native capability: the overlap is a closed-form trig integral (`Integrate`), the quartic identity and degree are polynomial-algebra facts (`Together`/`Numerator`/`Exponent`/`CoefficientList`), the residue/linewidth cancellation is a `FullSimplify` of a symbolic ratio, the pole census is `NSolve`/`Solve` on a quartic, and the scan is numeric substitution. There is no genuine impossibility, so this is a `missing_mathematica` finding, not a single-engine carve-out.

**Why this matters:**
The numerical pole census, R_{Q,*} figures, and the scan table (deliverable 5 and 6) are cross-engine-checkable; without a second engine the F1 discrepancy (213 vs 145) had no independent arbiter, and any single-engine root-finding/precision artifact goes unchallenged.

**Required change:**
Codex creates an independent Mathematica audit (see directive F2 with claim manifest M1-M6 and anti-transliteration guard). Target path is mandatory and exact.

**Verification:**
`redteam exec-mathematica 222` runs the new `.wl`, it exits 0, and its independently-derived numerics (kappa, quartic degree, R_{Q,*} cancellation, the four poles, the four R_{Q,*} figures, the two thresholds, the scan) agree with the SymPy figures (and surface whether 145 or 213 is the correct lambda_W=0.2 value).

## Independent-derivation check (Mathematica)

No `.wl` exists; transliteration cannot be assessed. The directive's claim manifest mandates native primitives and a different decomposition (Resultant/CoefficientList for the quartic, FullSimplify for the cancellation, NSolve over the quartic in y, Integrate for kappa) so the new engine does not echo the SymPy choreography.

## Engine cross-check

Only one engine present; cross-check deferred to the new `.wl`. The F1 discrepancy is precisely the kind of single-row error a second engine would catch.

## Verdict justification

`verdict: findings`, `stop_cold: null`. Five of six deliverables hold up under attack: the kappa integral, the quartic identity and degree-4 check, the residue/linewidth cancellation (genuinely non-tautological — `D0prime` symbolically cancels and the `Gamma_5` substitution produces the boxed `27 c_s^5/(a^5 w_*^5 N_star)` form), the static slice data, the pole census, the four R_{Q,*} figures, the two illustrative thresholds, and the monotonicity claims all verify; literals are regression anchors over in-script symbolic derivations, not standalone hardcoded answers; `assert_close` tol=1e-12 has adequate margin (observed residuals ~1e-13). The two open items are (F1) a published lambda_W=0.2 R_{Q,*} value in the notes (213.48) that the script does not reproduce (145.48), needing user resolution, and (F2) the absence of a Mathematica engine that the dual-engine rule requires here. I read the paper card, the full notes, and the relevant appendix rows before the script, and the script's verified claims match the paper for all deliverables except the one disputed scan-row value.

## Self-test notes

Checked: (1) Variable independence — `RQ_expr` (line 197) depends on `omega` and (through `N_omega`, line 153) on all coupling/frequency symbols, so the `.subs(omega, pole)` evaluations are non-trivial; no identically-zero derivative pattern is used. (2) Symmetry/parity — n/a, no unbounded-domain integrals; kappa is over [0,L] and evaluates to a nonzero constant 2√2/π. (3) Trivial-case — the residue/linewidth assertion (line 117) is the cancellation identity itself; substituting Gamma_5=a^5/(27c_s^5) gives exactly the RHS, confirmed by hand. (4) Paths — `.wl` target is in `mathematica/` with the exact mandated filename including `_mathematica_audit`. (5) Paper round-trip — F1 is left for user resolution (no Codex edit prescribed); the F2 manifest reuses only paper/notes-stated constants (Gamma_5, the sample slice, V_known) so it introduces no new constant.
