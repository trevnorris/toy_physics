---
unit_id: 202
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
  notes_stage_files: [moving_throat_pde_stage202_free_quintuple_target_graph.md]
  paper_appendix: present
---

# Audit unit 202 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_202.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 35, 236, 387, 502; eqs. `app-part06-Ctr-direct` … `app-part06-graph-family-zero`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.txt`

## What the paper claims

Stage 202 removes the abstract target-orbit language of Stage 201 by solving the three carried monomial constraints (`\mathfrak C_{tr,*}`, `\mathfrak C_{nt,*}`, `\epsilon_\eta`) explicitly for the dependent triple `(T_U, K_\eta^{eff}, \mu_W)` as a *graph* over the five free microscopic coordinates `\mathbf y=(\lambda_W,c_{\eta U},\gamma,K_U,K_W^{eff})`. `\stagefield{Output}` reads: "Exact target graph `\mathbf x_*^{graph}(\mathbf y)`, graph errors `(E_T,E_K,E_\mu)`, and the first reduced-family test in graph-error coordinates." The notes enumerate six deliverables: (1) the graph map `\Phi_*` producing `T_U^{graph}, K_{\eta,*}^{graph}, \mu_W^{graph}`; (2) the theorem `\mathcal O_*=\{\mathbf x_*^{graph}(\mathbf y)\}` (substitution + uniqueness); (3) the graph-error packet `(E_T,E_K,E_\mu)`; (4) the exact identities `E_T=q_{tr}/(1+\chi_{0,*})`, `E_K=-q_\eta`, `E_\mu=q_{nt}-q_\eta+F_* q_{tr}/(1+\chi_{0,*})`; (5) the repair law `\Delta\mathbf x_{rep}=(0,0,0,0,-E_K,0,-E_\mu,-E_T)^T`; (6) the first reduced-family closure test `\Delta_{fam}^{graph}=0 \iff \mathbf x_{cand}\in\mathcal Z_*`. The notes also assert `\mu_W^{graph}` is independent of `L` and `\pi`, and that the graph projection equals the Stage 201 canonical projection. NOTE: the card's `\stagefield{Verification}` line says "Mathematica audit: none yet," which is now stale — a `.wl` exists (see F1).

## What the script claims to verify

The SymPy script (a) writes the three carried monomials, (b) *posits* the closed-form graph solutions `deltaU_graph, T_graph, Keta_graph, mu_graph` taken from the boxed paper formulas and verifies them by substituting back into the monomials and checking `log(reconstructed/target)==0`, (c) checks `mu_graph` is free of `L`, (d) checks the graph projection equals the canonical Stage 201 projection componentwise, (e) checks the three graph-error identities and the `m_T/m_K/m_\mu` rewrite, (f) checks the repair-vector rewrite, and (g) checks the reduced-family packet (and multiplicative chart) vanishes on the graph with `\chi_Q=1`. The Mathematica script *derives* the graph by `LinearSolve` of a 3x3 log-linear dependent-triple system, then runs the same family of substitution/identity residual checks (M1–M7). All assertions are residual-equals-zero or boolean-true checks.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Graph map `T_U^{graph},K_{\eta}^{graph},\mu_W^{graph}` | py II posits+substitutes; wl I LinearSolve + M1–M3 | match |
| `\mu_W^{graph}` free of `L,\pi` | py L-not-in-free_symbols; wl M4 `FreeQ`+`D[,L]` | match |
| (2) `\mathcal O_*=` graph (substitution leg) | py II `Ctr/Cnt/eps` reconstructions; wl M1–M3 | match (substitution leg only; uniqueness leg is by Stage-192 pivot, not re-proved — see note) |
| graph proj = canonical proj | py III 8-component residual; wl M6 | match |
| (3)(4) graph errors + identities | py IV; wl M5 | match |
| (5) repair-vector rewrite | py V; wl M7 | match |
| (6) reduced-family packet zero | py VI; wl M7 | match |
| `\Delta_{fam}=0 \iff \in\mathcal Z_*` (the iff) | py/wl check only the forward `\Rightarrow` (vanishing on graph) | partial (the converse rests on the substitution+projection theorem, which the scripts do exercise) |
| card "Mathematica audit: none yet" | `.wl` now present | mismatch (stale card line) → F1 |

Dominant pattern: aligned. The one mismatch is a stale `\stagefield{Verification}` line, routed to the user (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 115–118 | `log(Ctr[graph]/Ctr_tgt)==0` | (1)/(2) tracking | yes |
| A2 | sympy | 119–122 | `log(Cnt[graph]/Cnt_tgt)==0` | (1)/(2) nontracking | yes |
| A3 | sympy | 123–126 | `log(eps[graph]/eps_tgt)==0` | (1)/(2) dressing | yes |
| A4 | sympy | 129–131 | `L not in mu_graph.free_symbols` | `\mu_W^{graph}` L-free | yes |
| A5 | sympy | 162–165 | 8-comp `log(x_proj_can/x_graph)==0` | proj=canonical | yes |
| A6 | sympy | 187–192 | three graph-error identity residuals | (4) | yes |
| A7 | sympy | 198–200 | `E - log(m)` rewrites | (3)/(4) | yes |
| A8 | sympy | 227 | repair-vector residual | (5) | yes |
| A9 | sympy | 241–243 | family packet on graph ==0 | (6) | yes |
| A10 | sympy | 250–253 | multiplicative chart on graph ==1 | (6) | yes |
| M0 | mathematica | 142–145 | `depMatrix.sol - depRhs == 0` | LinearSolve consistency | yes |
| M1 | mathematica | 149–152 | `logCtr[logTGraph]-Log[Ctrtgt]==0` | (1)/(2) tracking | yes |
| M2 | mathematica | 154–157 | `logEta[logKetaGraph]-Log[epsEtatgt]==0` | (1)/(2) dressing | yes |
| M3 | mathematica | 159–162 | `logCnt[...]-Log[Cnttgt]==0` | (1)/(2) nontracking | yes |
| M4 | mathematica | 165–169 | `FreeQ[L,Pi]` + `D[logMu,L]==0` | `\mu_W^{graph}` L,π-free | yes |
| M5 | mathematica | 181–188 | three graph-error residuals | (4) | yes |
| M6 | mathematica | 217–218 | 8-comp proj residual + threaded eqns | proj=canonical | yes |
| M7 | mathematica | 236, 246–249 | repair residual + family packet ==0 | (5)/(6) | yes |

All rows non-tautological: each substitutes the *derived/posited* graph into an *independently-written* monomial/identity and checks a residual that would be nonzero if the graph formula were wrong.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (stale verification line)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_202.tex:11` quote: "SymPy audit: \StageFile{…stage202…sympy_audit.py}.  Mathematica audit: none yet."
- script side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl:1–252` (a full Mathematica audit exists and passes)
- notes side: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md:549–552` lists only the `.py` under "Supporting files" (no `.wl`).

**What's wrong:**
The card's `\stagefield{Verification}` says the Mathematica audit does not exist ("none yet"), and the notes' Supporting-files list omits the `.wl`. But a complete `.wl` is present and passes (output `…mathematica_audit.txt` shows all M1–M7 PASS, `Exit[0]`). This is a documentation/script reflection mismatch (the script side now does *more* than the prose records). Direction of resolution is the user's call — almost certainly "update the card/notes to record the Mathematica audit," but the auditor does not silently edit prose.

**Why this matters:**
A reader trusting the card would believe stage 202 is single-engine when it is now dual-engine. It is also the kind of stale-prose signal that should be reconciled rather than left to drift.

**Required change:**
`## Resolve before fix_loop` (see directive). No Codex script edit.

**Verification:**
After user resolution, card line 11 reads "Mathematica audit: \StageFile{mathematica/…stage202…_mathematica_audit.wl}." and notes Supporting-files lists the `.wl`.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.**

The single discriminating operation: **for the load-bearing graph solution, the `.py` POSITS the closed forms; the `.wl` DERIVES them by `LinearSolve` of a log-linear system.**

- SymPy (`…_sympy_audit.py:92–104`) hand-writes the answer:
  ```
  deltaU_graph = (Ctr_tgt / (gamma*cetaU/KU)**(1+deltaUs))**(1/(1+chi0s))
  T_graph    = L**2*KU*deltaU_graph/pi**2
  Keta_graph = cetaU**2/(KU*epsEta_tgt)
  mu_graph   = Cnt_tgt*Keta_graph*KW**2/lam**2 * (...)**(-Estar) * deltaU_graph**Fstar
  ```
  i.e. it transcribes the boxed paper formulas and then *verifies* them by substitution (lines 115–126). No solve is performed in the `.py`.

- Mathematica (`…_mathematica_audit.wl:112–135`) instead constructs the dependent-triple system in log variables and solves it:
  ```
  depMatrix = {{1+chi0s,0,0},{0,1,0},{-Fstar,-1,1}};
  depRhs    = { Log[Ctrtgt] - (1+deltaUs) logTrackingBase - (1+chi0s)(2 Log[Pi]-2 Log[L]-Log[KU]),
                2 Log[cetaU]-Log[KU]-Log[epsEtatgt],
                Log[Cnttgt]-2 Log[lamW]+2 Log[KW]-Estar logDriveBase+Fstar(2 Log[Pi]-2 Log[L]-Log[KU]) };
  {logTGraph, logKetaGraph, logMuGraph} = normalize[LinearSolve[depMatrix, depRhs]];
  ```
  The graph values are an *output* of `LinearSolve`, then exponentiated (`TGraph=Exp[logTGraph]`, etc.).

This is exactly derive-vs-posit. Note the prompt's specific worry: stage 202 had a first-pass iter-2 rework from a transcendental `Solve` to a log-linearized `LinearSolve`. That linearization lives ONLY in the `.wl`; the `.py` never solves at all (it posits), so the two engines do *not* share a common linearization+LinearSolve route — there is no shared solve to be a port of. The downstream identity sections (graph errors M5, projection M6, repair/family M7) do check the same target identities in both engines, but that is the *verification* layer; checking the same paper-stated identity in both CASs is expected and is not a port, because the load-bearing graph object reaching those sections was produced by two genuinely different operations. The `expectZero`/`Simplify` machinery in the `.wl` is the audit's own normalizer, not an echo of SymPy's `simplify_log` choreography (one works in raw expressions + `expand_log`; the other works in pre-formed log-linear combinations).

## Engine cross-check

Both engines reach all-zero residuals on every shared check. SymPy output: "log(target Ctr/Cnt/eps reconstructed by graph) = 0", 8-component projection vector all 0, all graph-error/repair/family residuals 0, `L in mu_graph.free_symbols? no`. Mathematica output: graph log-system residual `{0,0,0}`, M1–M3 `=0`, M4 `True` and `d log(mu_graph)/dL =0`, M5 `{0,0,0}`, M6 8-comp `{0,…,0}` + `True`, M7 repair `{0,…,0}` + family `{0,0,0,0}`. Agreement is exact. The `.wl` adds two genuinely-extra coverage items the `.py` lacks: an explicit `LinearSolve` consistency residual (M0, line 142) and an explicit `D[logMuGraph,L]==0` derivative check (M4, line 169) on top of the `FreeQ` membership test — strengthening, not echoing.

## Verdict justification

Verdict `findings` solely because of F1, a low-severity stale-prose `paper_misalignment` (card/notes say "Mathematica audit: none yet" but a passing `.wl` exists). The math is sound and dual-engine: I attacked the `δ_{U,*}` self-referential exponent (it is a CARRIED constant `deltaUs`/`deltaU_star`, NOT the solved `δ_U`; the paper, `.py`, and `.wl` all treat the `(1+δ_{U,*})` exponent on `(γc_ηU/K_U)` as a fixed symbol — internally consistent, no circularity); I attacked the substitution checks for tautology (they substitute a posited/derived graph into independently-written monomials, so they can fail); I attacked the projection/repair/family identities (anchored to paper eqs.); and I attacked the `.wl` for transliteration (it is not — derive-vs-posit on the load-bearing object). Both outputs PASS with exit 0. I confirmed I read the card, notes, and the part-06 appendix rows, and the script claims match the paper claims except the stale verification line.

## Self-test notes

(1) Variable independence: the only differentiation is `.wl:169` `D[logMuGraph, L]`; `logMuGraph` is the `LinearSolve` output, which (correctly, per the elimination) contains no `Log[L]` — so the derivative is identically 0 and the `expectZero` passes for the right reason, *corroborated* by the independent `FreeQ[muGraphFlat,L]&&FreeQ[muGraphFlat,Pi]` truth check, so it is not a vacuous-derivative trap. (2) No unbounded integrals — parity trap N/A. (3) Trivial-case: substituting the graph back collapses each monomial ratio to 1 (log 0), matching both outputs. F1 is a paper_misalignment, so no Codex script edit is prescribed; the directive carries only the `## Resolve before fix_loop` block.

## Value Reconciliation (pass-2 augmentation)

The scripts emit no benchmark numeric constants; every deliverable is a closed-form symbolic object. Note the SymPy `.txt` is STALE (mtime 2026-05-11, older than the `.py` mtime 2026-06-03) and carries a stale "STAGE 185" banner/renumber label (a known numbering-band deferral, not a content error); reconciliation below is grounded in the script *source* plus the fresh `.wl` output. All symbolic deliverables reconcile against the notes (`.md`) and/or appendix (`.tex`).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `\delta_{U,*}^{graph}(\mathbf y)` | py:92–94; wl:131–132 (out wl L9) | notes:177–185 / app `eq:app-part06-deltaU-graph` L408–416 | MATCH |
| `T_U^{graph}(\mathbf y)` | py:95; wl:129 (out wl L10) | notes:207–213 / app `eq:…TU-graph` L417–422 | MATCH |
| `K_{\eta,*}^{graph}(\mathbf y)=c_{\eta U}^2/(K_U\epsilon_{\eta,*}^{tgt})` | py:96; wl:129 (out wl L11) | notes:223–228 / app `eq:…Keta-graph` L423–428 | MATCH |
| `\mu_W^{graph}(\mathbf y)` | py:97–104; wl:129 (out wl L12) | notes:249–262 / app `eq:…muW-graph` L429–437 | MATCH |
| `\mu_W^{graph}` independent of `L` (and `\pi`) | py:129–131 (out py L90); wl:165–169 (out wl L26–29) | notes:264 ("independent of `L` and `\pi`") | MATCH |
| graph proj = canonical proj | py:162–165 (out py L253–268); wl:217 (out wl L42) | notes:347–351 / app `eq:…can-equals-graph` L463–468 | MATCH |
| `E_T=q_{tr}/(1+\chi_{0,*})` | py:187 (out py L320); wl:181–188 (out wl L35) | notes:392 / app row eq | MATCH |
| `E_K=-q_\eta` | py:188 (out py L321); wl:181–188 | notes:398 | MATCH |
| `E_\mu=q_{nt}-q_\eta+F_* q_{tr}/(1+\chi_{0,*})` | py:189–192 (out py L322); wl:181–188 | notes:402 | MATCH |
| repair vector `(0,0,0,0,-E_K,0,-E_\mu,-E_T)^T` | py:227 (out py L478–493); wl:236 (out wl L51) | notes:421–434 / app L420–434 | MATCH |
| reduced-family packet `(\chi_Q-1,E_T,E_K,E_\mu)=0` on graph | py:241–253 (out py L554–561,618–625); wl:246–249 (out wl L54) | notes:460–483 / app `eq:…graph-family-zero` L494–500 | MATCH |

INTERNAL scaffolding (no finding): `depMatrix`, `depRhs`, `logSplit/logTrackingBase/logDriveBase`, `Rtr_expr/Rnt_expr/Reta_expr`, `qtr/qnt/qeta`, `m_T/m_K/m_\mu`, `M_family`, all residual vectors, PASS/FAIL flags, the graph log-system residual (M0).

reconciliation: complete; 11 deliverable values checked, 0 misaligned. (The only documentation discrepancy is the stale "Mathematica audit: none yet" card line and the notes' Supporting-files omission of the `.wl`, captured as F1, not a value mismatch.)
