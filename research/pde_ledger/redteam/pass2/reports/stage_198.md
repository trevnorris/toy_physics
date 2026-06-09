---
unit_id: 198
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage198_exact_finite_orbit_law.md]
  paper_appendix: present
---

# Audit unit 198 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_198.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage198_exact_finite_orbit_law.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 127, 265, 1329–1361, 1468)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage198_exact_finite_orbit_law_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage198_exact_finite_orbit_law_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` says the stage "Solves the dependent triple exactly, defines mismatch coordinates, and proves finite single-orbit lock." The notes enumerate five deliverables: (1) the exact finite single-orbit law for the dependent triple `(T_U, K_eta^(eff), mu_W)` from the three triangular coherent monomials `C_tr_*, C_nt_*, eps_eta_*`; (2) the exact dependent residual mismatch triple `(m_T, m_K, m_mu)` via the invariant ratios `R_tr = m_T^(1+chi0_*)`, `R_eta = 1/m_K`, `R_nt = m_mu/(m_K m_T^F_*)`; (3) the exact logarithmic chart `q_tr=(1+chi0_*)tau`, `q_eta=-kappa`, `q_nt=mu-kappa-F_* tau`, shown equal to the Stage-192 drift compiler `M_*` applied to a pure dependent mismatch vector; (4) the exact restoration map that returns a candidate branch to the orbit by changing only the dependent triple; (5) the sharp orbit-lock theorem `Delta_orbit=0 <=> m_T=m_K=m_mu=1 <=> q_tr=q_nt=q_eta=0`. The appendix (rows 1329–1361) carries the same free/dependent split, the same log chart, and the same orbit-lock equivalence verbatim.

## What the script claims to verify

The SymPy script defines the three exact monomials (`Ctr`, `Cnt`, `epsEta`), types the three closed-form orbit values (`Keta_orbit`, `TU_orbit`, `muW_orbit`) directly from the notes, and asserts (a) each orbit value reproduces its target invariant (`Ctr(TU_orbit)/Ctr_star-1=0`, etc.); (b) the three invariant ratios collapse exactly to `m_T^(1+chi0_*)`, `1/m_K`, `m_mu/(m_K m_T^F_*)`; (c) the explicit `M_*` matrix applied to the pure dependent mismatch vector `(0,0,0,0,kappa,0,mu,tau)` equals `(q_tr,q_nt,q_eta)`; (d) the restoration map returns each orbit value; (e) the inverse chart round-trips and `m=1` at `q=0`. The Mathematica script verifies the SAME five deliverables but builds them by independent routes: it SOLVES for the orbit by a log-linear `Coefficient`-matrix + `LinearSolve` (rather than substituting closed forms), reconstructs `M_*` as a finite-perturbation `D[...]` Jacobian of the monomial laws, performs restoration as a dependent-column `LinearSolve`, and inverts the chart with `Solve[...,Reals]`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) finite orbit law `(T_U^orbit, K_eta^orbit, mu_W^orbit)` | py L97–108 `Ctr/Cnt/epsEta(orbit)=star`; wl L155–162 log-linear solve + agreement | match |
| (2) mismatch triple / invariant ratios `R_tr,R_eta,R_nt` | py L134–136; wl L168–182 | match |
| (3) log chart `q_tr,q_eta,q_nt` = `M_*` on pure dep. mismatch | py L144–170; wl L191–229 (Jacobian compiler) | match |
| (4) restoration map (dependent triple only) | py L181–196; wl L233–265 (column LinearSolve) | match |
| (5) orbit-lock theorem `Delta_orbit=0 <=> m=1 <=> q=0` | py L227–238; wl L269–300 (`Equivalent` + `m=1@q=0`) | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological script-side counterpart on BOTH engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 97–100 | `simplify(epsEta(Keta_orbit)/epsEta_star -1)==0` | (1) | yes |
| A2 | sympy | 101–104 | `simplify(Ctr(TU_orbit)/Ctr_star -1)==0` | (1) | yes |
| A3 | sympy | 105–108 | `simplify(Cnt(orbit)/Cnt_star -1)==0` | (1) | yes |
| A4 | sympy | 134 | `Ctr_ratio - m_T^(1+chi0)==0` | (2) | yes |
| A5 | sympy | 135 | `epsEta_ratio - 1/m_K==0` | (2) | yes |
| A6 | sympy | 136 | `Cnt_ratio - m_mu/(m_K m_T^F)==0` | (2) | yes |
| A7 | sympy | 167–170 | `M_* Dx_mis - (q_tr,q_nt,q_eta)==0` | (3) | yes |
| A8 | sympy | 194–196 | `T/K/mu_restore - orbit ==0` | (4) | yes |
| A9 | sympy | 214–238 | inverse-chart round-trip + `m=1@q=0` | (5) | yes |
| B1 | mathematica | 142 | `Det(dep log matrix)+(1+chiS)==0` | (1) | yes |
| B2 | mathematica | 157–162 | `Log(orbitFromLogs/SymPy target)==0` ×3 + invariants | (1) | yes |
| B3 | mathematica | 168–182 | scaling ratios = `sT^(1+chiS)`, `1/sK`, `sM/(sK sT^fS)` | (2) | yes |
| B4 | mathematica | 196,225,228 | Jacobian compiler = `M_*` and `M_* . dep = q` | (3) | yes |
| B5 | mathematica | 242–265 | column `LinearSolve` restoration + returns orbit | (4) | yes |
| B6 | mathematica | 277–299 | inverse chart `Solve` + `m=1@q=0` + `Equivalent` lock | (5) | yes |

All rows anchored; no tautological or orphaned assertions. The SymPy `M_*` (L156–162) is a typed matrix compared against a typed mismatch vector, but its entries are non-trivially exercised by the matrix-vector product producing the independently-defined `qtr/qnt/qeta_expected` — A7 fails if any matrix entry is wrong, so it is not tautological.

## Findings

None.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT. The two engines solve the SAME deliverables by structurally different methods; the `.wl` is not a transliteration.

- **Orbit law (deliverable 1).** SymPy *posits* the closed forms and checks them: `Keta_orbit = sp.simplify(cEtaU**2/(KU*epsEta_star))`, `TU_orbit = ...` (py L79–88), then `expect_zero("Ctr(TU_orbit)/Ctr_star -1", ...)` (py L101–104). Mathematica instead *solves* for them from scratch: it linearizes the log-residuals, extracts the dependent coefficient matrix `depCoefficientMatrix = Table[Coefficient[invariantLogResiduals[[row]], depLogVars[[col]]], {row,3},{col,3}]` (wl L126–129) and `logOrbitSolution = LinearSolve[depCoefficientMatrix, -depOffset]` (wl L133–135), `orbitFromLogs = Exp/@ logOrbitSolution` (wl L136). It then *cross-checks* the solved values against the SymPy closed forms (`samePositive["K_eta orbit agrees with SymPy target", orbitFromLogs[[1]], kEtaSympy]`, wl L157). The SymPy formula is the comparison TARGET, not the method — the Mathematica orbit is produced by `Coefficient`+`LinearSolve`, an independent route. This is the textbook "solve vs. substitute-and-verify" split.

- **Drift compiler `M_*` (deliverable 3).** SymPy hard-types the matrix: `Mstar = sp.Matrix([[0,1+deltaUs,...],...])` (py L156–162). Mathematica *reconstructs* it as a Jacobian of the actual monomial laws under multiplicative perturbation: `driftRules` exponentiate each variable (wl L199–208), `finiteLogPacket = Log[(laws /. driftRules)/laws]` (wl L209–211), then `compilerByJacobian = Table[D[finiteLogPacket[[row]], driftVars[[col]]], {row,3},{col,8}]` (wl L212–215). The hard-typed matrix appears only as `compilerExpected` (wl L216–220), the comparison target. The Jacobian is genuinely computed from the laws, not copied.

- **Restoration (deliverable 4).** SymPy applies the explicit closed-form exponential restoration map and checks it returns the orbit (py L181–196). Mathematica solves restoration as a linear system on the dependent columns: `dependentColumns = compilerByJacobian[[All,{5,7,8}]]`, `restorationSolve = LinearSolve[dependentColumns, -packetExpected]` (wl L233–234). Different machinery, same result.

**Self-test trap (Jacobian zeros are real, not vacuous).** I confirmed each nonzero `D[...]` entry differentiates a variable that genuinely appears in `finiteLogPacket`, and each zero entry corresponds to a variable genuinely absent (output L42): row 1 `= dC+dC*deltaS+dG+deltaS*dG+dT+chiS*dT-(2+chiS+deltaS)*dU` depends only on {dC,dG,dT,dU}, so the zeros in cols {1,5,6,7} are real; row 2 depends on all of {dLam,dG,dU,dK,dW,dMu,dT} but NOT dC, so its col-2 zero is real; row 3 `= 2*dC-dK-dU` makes the other five columns genuinely zero. No entry passes vacuously through an absent variable.

## Engine cross-check

Both engines exit cleanly with all checks PASS. SymPy output shows every `expect_zero` residual as `0` (e.g. L72–74, L91–93, L110–115, L157–159, L190–192). Mathematica output shows every `PASS:` line (L13–79), including `det(dependent log matrix)+(1+chiS)=0` (L12), the three "orbit agrees with SymPy target = 0" lines (L14–19), the Jacobian agreement `{{0,...},{0,...},{0,...}}` (L44), the restoration `{0,0,0}` (L53,55), and the `Equivalent` lock `= True` (L78). The two engines agree on every shared deliverable.

## Verdict justification

CLEAN. I attacked the orbit solve (tried to find a hidden hardcode — the `.wl` genuinely solves via `LinearSolve`, the `.py` posits-then-verifies, and they cross-check), the `M_*` matrix (SymPy types it but A7's matrix-vector product non-trivially exercises every load-bearing entry; Mathematica reconstructs it as a real Jacobian), the Jacobian self-test (all zeros correspond to genuinely-absent variables), the restoration map (independent column-LinearSolve vs. closed-form), and the orbit-lock `Equivalent` claim (a real biconditional, not a tautology). I read the card, notes, and appendix first: the five paper deliverables map one-to-one to substantive, anchored checks on both engines, and the constants/exponents match the notes' boxed formulas exactly. The two engines are genuinely independent, not a transliteration. The only non-blocking issue is a stale SymPy output banner (`STAGE 181`, =198−17, see reconciliation) which the orchestrator's re-run refreshes, and a card-text lag (`Mathematica audit: none yet` despite a passing `.wl`) which is paper-side and out of red-team script scope.

## Self-test notes

Checked: (1) variable-independence of every `D[...]` Jacobian entry — all zeros trace to genuinely-absent drift variables (output L42), so the compiler reconstruction is not vacuous; (2) trivial-case for the orbit-lock biconditional — `Equivalent[packet==0, logs==0]` is a real two-way claim that fails if the chart degenerates; (3) no unbounded integrals in this stage, so parity trap N/A. No directive needed; no self-test failure.

## Value Reconciliation (pass-2 augmentation)

This stage is entirely symbolic — it emits no named numeric constants (no `Pi_star`, `gamma_0`, benchmark numbers). The deliverables are closed-form symbolic results, each of which is a boxed equation in the notes and (for the load-bearing ones) the appendix.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `K_eta^(orbit) = c_etaU^2/(K_U eps_eta_*)` | py L79, out L35–39; wl L144, out L14 | notes L139–143 (boxed) | MATCH |
| `T_U^(orbit) = (L^2 K_U/pi^2)[C_tr_*/(gamma c_etaU/K_U)^(1+deltaU_*)]^(1/(1+chi0_*))` | py L80–83, out L40–49; wl L145–148, out L16 | notes L150–158 (boxed) | MATCH |
| `mu_W^(orbit) = C_nt_* K_eta^orbit K_W^2/lambda_W^2 (gamma^2 lambda_W^2 sigma/(K_U K_W))^(-E_*)(pi^2 T_U^orbit/(L^2 K_U))^(F_*)` | py L84–88, out L50–71; wl L149–153, out L18 | notes L165–175 (boxed) | MATCH |
| `R_tr = m_T^(1+chi0_*)` | py L134, out L91; wl L168–172 | notes L220–223,246; appendix L1347 | MATCH |
| `R_eta = 1/m_K` | py L135, out L92; wl L173–177 | notes L226–229,250; appendix L1349 | MATCH |
| `R_nt = m_mu/(m_K m_T^F_*)` | py L136, out L93; wl L178–182 | notes L232–235,248 | MATCH |
| `q_tr = (1+chi0_*)tau`, `q_eta=-kappa`, `q_nt=mu-kappa-F_* tau` | py L144–146, out L98–103; wl L194 | notes L268–274; appendix L1347–1351 | MATCH |
| `M_*` drift matrix (3×8 rows) | py L156–162; wl L216–220, out L43 | notes L284–301 (matches the three encoded rows) | MATCH |
| restoration map (3 exp factors) | py L181–186, out L120–155; wl L235–239 | notes L314–333 (boxed) | MATCH |
| orbit-lock `Delta_orbit=0 <=> m=1 <=> q=0` | py L236–238; wl L295–300, out L78 | notes L357–369; appendix L1354–1361 | MATCH |

INTERNAL (scaffolding, no prose expected): `tau/kappa/mu` log symbols, `tau_inv/kappa_inv/mu_inv` inverse-chart expressions, `depCoefficientMatrix`, `depOffset`, `logOrbitSolution`, `finiteLogPacket`, `compilerByJacobian`, `dependentColumns`, `restorationSolve`, all `expect_zero`/`PASS` residuals.

reconciliation: complete; 10 deliverable values checked, 0 misaligned.

Note: the committed SymPy `.txt` carries a stale `STAGE 181` banner (out L3, L195) = 198−17 (known pre-renumber drift); the `.py` source itself prints `STAGE 198` (py L35). The `.txt` mtime (2026-06-01 12:18) predates the `.py` mtime (2026-06-03 15:59), so the output is stale — informational only; the orchestrator re-run refreshes it and the reconciliation above is based on the script source. (See `stale_output` note in the standard sections; not raised as a blocking finding because no value content disagrees.)
