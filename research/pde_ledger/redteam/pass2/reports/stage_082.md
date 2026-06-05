---
unit_id: 082
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md"]
  paper_appendix: present
---

# Audit unit 082 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_082.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 142, definition preview lines 168–186)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.txt`

## What the paper claims

Stage 082 packages the entire reduced moving-throat PDE, on the surviving passive/outgoing quadrupole branch, into a single scalar diagnostic. The boxed `\stagefield{Output}` is "the master residual \eqref{eq:app-stage082-Rquad}", namely `R_quad = zeta_req(Pi_tr,C_mix,eps_blk) - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa)`, with the demand `zeta_req = (Pi_tr - C_mix)/(C_mix - eps_blk(2 C_mix - Pi_tr))` (eq. app-stage082-zeta-req) and supply `zeta_phys = Omega_Pe^2 (kappa+pi^2/4)/(kappa+y(eta)^2)` (eq. app-stage082-zeta-phys), plus the sign convention `R_quad<0` (support exceeds demand) / `=0` (saturation) / `>0` (support failure). The card's `\claimstatus` is explicit: `\StatusExactClosure{}` for the *residual definition*, `\StatusOpen{}` for the *actual PDE branch values* inserted into it — so the stage does NOT claim to solve the Stage-058 `Pe_*` fixed point; it claims the residual is correctly assembled and the surrounding exact algebra holds. The notes add four further exact deliverables that feed the residual: (1) the inverse map `Pi_tr = C_mix Q(zeta_req;eps_blk)` with `Q = [1+(1-2eps_blk)zeta]/[1-eps_blk zeta]` (§4), (2) the product thresholds `Pi_suff = C_mix Q(zeta_-)`, `Pi_fail = C_mix Q(zeta_+)` (§4), (3) the strict monotonicity of `zeta_req` in `Pi_tr` underlying that inversion (§4), and (4) the Family-1 specialization `(eta,kappa)=(37,12321/5)`, `Xi_F1 = W_wall = Upsilon_w Lambda_ell^2 = 1369 Upsilon_w = 136900 Theta_w` (§5, eq. app-stage082-XiF1). The appendix row 142 summarizes it as "Full reduced residual `R_quad = zeta_req - zeta_phys`".

## What the script claims to verify

The two engines verify the exact algebra around the residual (not the open PDE values, consistent with `claimstatus`). The SymPy docstring lists four checks: (1) the exact inverse map between `zeta_req` and `Pi_tr` via `Q`; (2) the exact product thresholds `Pi_suff`/`Pi_fail`; (3) the Family-1 strength identity `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w`; (4) the master residual definition. Concretely the scripts assert: the round-trip `zeta_req(C_mix*Q(zeta)) - zeta == 0`; `zeta_req(Pi_suff) - zeta_- == 0` and `zeta_req(Pi_fail) - zeta_+ == 0`; the residual collapses `R_quad(Pi_suff, zeta_phys=zeta_-) == 0` and `R_quad(Pi_fail, zeta_phys=zeta_+) == 0`; the monotonicity factorization (numerator `= C_mix(1-eps_blk)`, denominator `= (C_mix - eps_blk(2 C_mix - Pi_tr))^2`); the `zeta_phys` closed form's `Pe->oo` limit `Omega_Pe -> pi/2` and the resulting Family-1 ceiling reproducing the upstream `zeta_max^(F1) = 2.46752922945601` to `1e-10`; and the (explicitly "display only") arithmetic forms of `Xi_F1`.

## Paper ↔ script cross-check

| Paper / notes deliverable | Script-side check | Status |
|---|---|---|
| `R_quad = zeta_req - zeta_phys` definition (Output, eq. app-stage082-Rquad) | sympy L64 `R_quad = zeta_req - zeta_phys`; collapse checks L68–75; wl L59–69 | match |
| `zeta_req = (Pi_tr - C_mix)/(C_mix - eps_blk(2 C_mix - Pi_tr))` (eq. zeta-req) | sympy L38; wl derives it by `Solve` (L38–40) | match |
| Inverse map `Pi_tr = C_mix Q(zeta_req)` (notes §4) | sympy L46–49 round-trip; wl L45–48 round-trip | match |
| Thresholds `Pi_suff = C_mix Q(zeta_-)`, `Pi_fail = C_mix Q(zeta_+)` (notes §4) | sympy L51–60; wl L50–57 | match |
| Monotonicity of `zeta_req` in `Pi_tr` (notes §4, underlies inversion) | sympy L82–92 sign-controlled num/den factor; wl L76–86 | match |
| `zeta_phys = Omega_Pe^2(kappa+pi^2/4)/(kappa+y^2)` closed form (eq. zeta-phys, notes §1.4) | sympy L101–138 (Omega_Pe->pi/2 limit + F1 ceiling vs upstream `zeta_max^(F1)`); wl L95–113 | match |
| Family-1 `(eta,kappa)=(37,12321/5)` (eq. XiF1, notes §5) | sympy L123 `kappa_F1=12321/5`, root at eta=37 L118; wl L102,104 | match |
| `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w` (eq. XiF1) | sympy L147–161 ("display only"); wl L119–133 | match (honest display, anchored upstream) |
| Sign convention `R_quad<0/=0/>0` (eq. app-stage082-sign) | prose-only in script ledger print L164–166; not an executable claim | match (not a math identity; nothing to assert) |

`paper_alignment: aligned`. Every load-bearing identity in the card and notes has a faithful, non-tautological script check; the open PDE-value side is correctly left unverified, matching `\claimstatus`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46–49 | `simplify(zeta_req(C_mix*Q(zeta)) - zeta) == 0` | inverse map (notes §4) | yes |
| A2 | sympy | 59 | `zeta_req(Pi_suff) - zeta_- == 0` | threshold `Pi_suff` (§4) | yes |
| A3 | sympy | 60 | `zeta_req(Pi_fail) - zeta_+ == 0` | threshold `Pi_fail` (§4) | yes |
| A4 | sympy | 68–71 | `R_quad(Pi_suff, zeta_phys=zeta_-) == 0` | residual def (Output) | yes |
| A5 | sympy | 72–75 | `R_quad(Pi_fail, zeta_phys=zeta_+) == 0` | residual def (Output) | yes |
| A6 | sympy | 85–88 | `num(d zeta_req/dPi) - C_mix(1-eps_blk) == 0` | monotonicity (§4) | yes |
| A7 | sympy | 89–92 | `den(d zeta_req/dPi) - (C_mix-eps_blk(2C_mix-Pi))^2 == 0` | monotonicity (§4) | yes |
| A8 | sympy | 110 | `Omega_Pe(Pe->oo) - pi/2 == 0` | zeta_phys closed form (eq. zeta-phys) | yes |
| A9 | sympy | 135 | `|zeta_phys(F1,Pe->oo) - zeta_max^(F1)| < 1e-10` | F1 zeta_phys ceiling (§5, cross-stage 084) | yes |
| A10 | sympy | 160–161 | `Xi_F1` arithmetic forms (== 0, "display only") | F1 strength identity (eq. XiF1) | partial (display-only, honestly labeled) |
| B1 | wl | 45–48 | `FullSimplify[zeta_req(C_mix*Q(zeta)) - zeta] === 0`, zeta_req derived via `Solve` | inverse map | yes (independent inverse) |
| B2 | wl | 56–57 | `zeta_req(Pi_suff)-zeta_- === 0`, `(...)-zeta_+ === 0` | thresholds | yes |
| B3 | wl | 62–69 | `R_quad(Pi_suff,zeta_-) === 0`, `(Pi_fail,zeta_+) === 0` | residual def | yes |
| B4 | wl | 79–86 | monotonicity num/den factor | monotonicity | yes |
| B5 | wl | 101 | `Omega_Pe(Pe->Inf) - Pi/2 === 0` | zeta_phys closed form | yes |
| B6 | wl | 110–112 | `|zeta_phys(F1) - zeta_max^(F1)| < 1e-10` | F1 ceiling | yes |
| B7 | wl | 130–133 | `Xi_F1` arithmetic forms (display only) | F1 strength identity | partial (display-only, honestly labeled) |

A10/B7 are the only non-fully-anchored rows; both are explicitly self-labeled "(display only)" in the script and the in-script comments (sympy L155–159, wl L125–129) state plainly that they are arithmetic on hand-supplied integers and "not independent verifications" — the upstream stages 073/075 own those facts. This honest disclosure matches the card's own `claimstatus` (Open for inserted values) and the upstream anchoring (see Value Reconciliation), so it is not a finding.

## Findings

None.

I attacked the following and each held up:

- **Tautology on the inverse map (A1/B1).** `Q` and `zeta_req` are constructed as independent expressions; the round-trip `zeta_req(C_mix*Q(zeta)) - zeta` would be nonzero if `Q` were the wrong inverse. The Mathematica side strengthens this: it does NOT restate the SymPy closed form — it `Solve`s `PiTr == cMix*qMap(zetaSym)` for `zetaSym` (wl L38) and then checks the round-trip, an independent re-derivation. Not tautological.
- **Threshold redundancy (A2/A3).** These re-exercise the inverse at `zeta_±`; mildly redundant with A1 but not circular (they could still fail if `Pi_suff`/`Pi_fail` were mis-formed). They map to a distinct notes §4 deliverable. Acceptable.
- **Hardcoded F1 ceiling (A9/B6).** The `2.467529229455835` reference is NOT used as the answer — the script independently recomputes `zeta_phys(Pe->oo)` from the bisection root `y_F1` and the closed form, then cross-checks against the stage-084 reference to `1e-10`. The reported residual `1.77e-13` is genuine recomputation noise, well inside tolerance. Confirmed the reference equals stage 084 `.wl` L71 / output L26 / notes L144 (`2.46752922945601`). Not a hardcoded echo.
- **Delicate root stability.** The known tan-singularity root near `pi/2` is handled with `mpmath.findroot(..., solver='bisect')` on the bracket `(1.5, 1.55)` (sympy L118), exactly the robust route flagged in pass 1; Mathematica uses `FindRoot` seeded at `1.527` with `WorkingPrecision -> 30`. Both engines land on `1.5294824837146996499...` agreeing to 28+ digits. No stability defect.
- **Symbol assumptions.** `C_mix, Theta_w, Upsilon_w` positive; everything else real. The positivity is justified by the physical setup (mixing constant and wall-depth amplitudes are positive) and is needed only for the monotonicity sign argument; it does not over-constrain the inverse-map/threshold algebra. No `symbol_assumption_error`.
- **Constant alignment.** `eta=37`, `kappa=12321/5`, `Lambda_ell=37` (1369=37^2), `Upsilon_w=100 Theta_w` (136900=100*1369), `alpha_r=10` all reconcile with the card (L43/L46) and upstream stages 073 (L25/L74) and 075 (L7/L24). No `value_mismatch`.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. The decisive divergence is at the inverse map: SymPy hand-writes `Q = (1+(1-2eps_blk)zeta)/(1-eps_blk zeta)` and substitutes (sympy L39, L48), whereas Mathematica explicitly inverts — `zetaReqSolList = Solve[PiTr == cMix*((1+(1-2*epsBlk)*zetaSym)/(1-epsBlk*zetaSym)), zetaSym]` then strips the `ConditionalExpression` and renames (wl L37–40), with the comment "This forces Mathematica to find the inverse of qMap rather than restating the SymPy-side closed form." The F1 ceiling uses different primitives too: SymPy `mpmath.findroot(..., solver='bisect')` on `(1.5,1.55)` vs Mathematica `FindRoot[..., {ySym,1.527}, WorkingPrecision->30]`; both compute the same `Omega_Pe` limit via `sp.limit`/`Limit`. The shared identity-residual chassis (`expect_zero`/`expectZero` on the same comparisons) is acceptable: both arrive at the compared quantities by different routes. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at full precision:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| inverse map / thresholds / R_quad collapses | all `= 0` (lines 14,25,26,32,33) | all `=== 0`, PASS (lines 7,11,13,17,19) |
| monotonicity num/den | `= 0` (lines 36,37) | `=== 0`, PASS (lines 21,23) |
| `Omega_Pe -> pi/2` | `pi/2` (line 39) | `Pi/2` (line 26) |
| `y_F1` | `1.52948248371469964992710762240` (line 41) | `1.5294824837146996499271076224...` (line 28) |
| `zeta_phys(F1,Pe->oo)` | `2.4675292294560122333` (line 42) | `2.46752922945601223332958...` (line 29) |
| `|... - zeta_max^(F1)|` | `1.77e-13` (line 43) | `1.77e-13` (line 30) |
| `Xi_F1` forms | `1369 Upsilon_w`, `136900 Theta_w` (lines 47,49) | `1369*upsilonW`, `136900*thetaW` (lines 32,33) |

`engines_agree: true`.

## Verdict justification

`clean`. The paper card (with notes and appendix) claims the master residual `R_quad = zeta_req - zeta_phys` is exactly defined and that its surrounding algebra — the `Q`-inverse map, the `Pi_suff`/`Pi_fail` thresholds, the monotonicity of `zeta_req`, the `zeta_phys` closed-form ceiling, and the Family-1 strength identity — is exact, while explicitly leaving the inserted PDE branch values Open. Both scripts verify exactly that set of identities with non-tautological, well-anchored assertions; the Mathematica engine independently inverts via `Solve` rather than echoing SymPy; the delicate tan-singularity root is handled with the robust bisection/seeded solvers and the two engines agree to 28+ digits; every numeric constant reconciles with the card and the upstream 073/075 anchors; and the only display-only block (`Xi_F1` arithmetic) is honestly self-labeled and consistent with the card's `claimstatus`. I tried tautology, hardcoded-answer, symbol-domain, root-stability, transliteration, and constant-mismatch attacks; none produced a defect. Outputs are fresh (script mtimes `02:17:33`/`02:15:09` precede output mtimes `02:17:45`/`02:19:35`).

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_req = (Pi_tr - C_mix)/(C_mix - eps_blk(2 C_mix - Pi_tr))` | py L38 / out L5–8; wl out L5 | `stage_082.tex:16–18` (eq. zeta-req); notes §1.5 L112 | MATCH |
| `Q = [1+(1-2eps_blk)zeta]/[1-eps_blk zeta]` | py L39 / out L10–13; wl out L6 | notes §4 L185–186 | MATCH |
| `Pi_tr = C_mix Q(zeta_req)` (inverse map) | py L46–49 / out L14 | notes §4 L190 | MATCH |
| `Pi_suff = C_mix Q(zeta_-)` | py L51 / out L16–19 | notes §4 L200 | MATCH |
| `Pi_fail = C_mix Q(zeta_+)` | py L52 / out L21–24 | notes §4 L202 | MATCH |
| `R_quad = zeta_req - zeta_phys` | py L64 / out L28–31; wl out L15 | `stage_082.tex:29–32` (boxed Output); notes §2 L126–128; appendix L183 | MATCH |
| `dzeta_req/dPi_tr` num `= C_mix(1-eps_blk)`, den `= (C_mix-eps_blk(2C_mix-Pi))^2` | py L82–92 / out L35–37 | notes §4 (monotonicity, prose, L181) | MATCH (sign-content of §4) |
| `Omega_Pe -> pi/2` (Pe->oo) | py L108–110 / out L39–40 | notes §1.4 L101–102 (Omega_Pe form); limit implicit | MATCH (form), limit is INTERNAL derivation |
| `y_F1 = 1.529482483714699…` (root y tan y = 37) | py L118 / out L41; wl L102 / out L28 | notes §5 eta_F1=37 L214 (root not tabulated) | INTERNAL (intermediate; eta=37 anchored) |
| `zeta_phys(F1,Pe->oo) = 2.46752922945601…` | py L124–127 / out L42 | notes §5 (F1 residual, value not boxed here); = stage 084 notes L144 `zeta_max^(F1)≈2.46752922945601` | MATCH (cross-stage 084 reference) |
| `zeta_max^(F1)` reference `2.467529229455835` | py L132 / out L43; wl L107 | stage 084 `.wl:71`/out L26/notes L144 (`2.46752922945601`) | MATCH (agree to 1e-10; recomputed not pinned) |
| `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w` | py L148–166 / out L46–56; wl L119–139 | `stage_082.tex:46` (eq. XiF1); notes §5 L218 | MATCH |
| `(eta,kappa)_F1 = (37, 12321/5)` | py L123 / wl L104 | `stage_082.tex:43`; notes §5 L213–214 | MATCH |
| `Lambda_ell = 37` | py L147 / wl L119 | stage 073 `.tex:25`/notes L74 | MATCH (upstream anchor) |
| `Upsilon_w = 100 Theta_w` (alpha_r=10) | py L142–143,149 / wl L115–121 | stage 075 `.tex:7,24`/notes L108 | MATCH (upstream anchor) |

INTERNAL scaffolding (accounted for, no finding): `expect_zero`/`expectZero` residual values, all PASS flags, the `y_F1` root precision, `diff_to_reference`/`diffToReference` (1.77e-13), tolerance `1e-10`, `mpmath.mp.dps=30`/`WorkingPrecision->30`, the intermediate `zeta_phys_closed`/`OmegaPe` symbolic expressions, `Omega_Pe_limit` extraction.

reconciliation: complete; 15 deliverable values checked, 0 misaligned.

## Self-test notes

I checked variable-independence (the `sp.diff(zeta_req, Pi_tr)` / `D[zetaReq, PiTr]` genuinely depends on `Pi_tr` — both engines return a nonzero rational, output lines 35/20, so the monotonicity factorization is non-trivial), parity (no unbounded symmetric integrals in this stage), and the trivial-case substitution for each `expect_zero` (the inverse-map and threshold round-trips collapse to 0 only because `Q` is the true inverse, not by construction). Path/round-trip checks are moot — no directive is written since there are zero findings.
