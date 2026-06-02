---
unit_id: 202
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage202_free_quintuple_target_graph.md]
  paper_appendix: present
---

# Audit unit 202 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_202.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 33–39, 91, 155–190, 236, 239, 384–502, 1304–1310 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 202 (anchor MTDC-T10.1, status \StatusExactClosure) removes the abstract target-orbit language by solving the three carried direct microscopic monomial constraints (`\mathfrak C_{{\rm tr},*}`, `\mathfrak C_{{\rm nt},*}`, `\epsilon_\eta`) explicitly for the dependent triple `(T_U, K_\eta^{(\rm eff)}, \mu_W)` as a graph over the five free coordinates `\mathbf y=(\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\rm eff)})`. The `\stagefield{Output}` is: "Exact target graph `\mathbf x_*^{\rm graph}(\mathbf y)`, graph errors `(E_T,E_K,E_\mu)`, and the first reduced-family test in graph-error coordinates." The distinct deliverables (notes §§1–7, appendix §§"Direct microscopic target graph" / "Graph-error realization coordinates"): (1) the exact graph solve `\delta_{U,*}^{\rm graph}, T_U^{\rm graph}, K_{\eta,*}^{\rm graph}, \mu_W^{\rm graph}` reconstructing the target monomials by direct substitution; (2) the graph theorem `\mathcal O_* = \{\mathbf x_*^{\rm graph}(\mathbf y)\}`; (2b) `\mu_W^{\rm graph}` independent of `L`; (3) graph projection equals the Stage-201 canonical orbit projection, `\Pi^{\rm graph}=\Pi^{\rm can}`; (4) the graph-error packet `(E_T,E_K,E_\mu)` and its exact identities `E_T=q_{\rm tr}/(1+\chi_{0,*})`, `E_K=-q_\eta`, `E_\mu=q_{\rm nt}-q_\eta+F_* q_{\rm tr}/(1+\chi_{0,*})`; (5) the repair vector `\Delta\mathbf x_{\rm rep}=(0,0,0,0,-E_K,0,-E_\mu,-E_T)^T`; (6) the first reduced-family closure packet `(\chi_Q-1,E_T,E_K,E_\mu)` vanishing on the target graph. The notes file is written under an older internal numbering ("Stage 253/252/243/237") but the physics maps one-to-one onto the Stage 202 card and the Part VI appendix.

## What the script claims to verify

The SymPy script defines the three direct monomials with the paper's exponents (`1+\delta_{U,*}`, `1+\chi_{0,*}`, `E_*`, `-F_*`), solves the graph triple (`deltaU_graph`, `T_graph`, `Keta_graph`, `mu_graph`) matching the appendix formulas verbatim, and then asserts: (II) substituting the graph solve back into each monomial reconstructs the target value (`log(...)=0`) and `mu_graph` is `L`-independent (`L not in free_symbols`); (III) an independently-built canonical projection `x_proj_can` (assembled from the monomial ratios Rtr/Rnt/Reta, *not* from the graph) equals `x_graph` componentwise; (IV) the three graph-error identities against the q-chart and against the multiplicative chart `m_T,m_K,m_\mu`; (V) the repair-vector rewrite in q-chart equals the graph-error chart; (VI) the reduced-family packet and multiplicative family chart vanish on the target graph at `\chi_Q=1`. Saved output exits 0 with all checks `= 0`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) graph solve reconstructs targets | Section II `expect_zero` ×3 (lines 115–126) | match |
| (1b) graph formulas match appendix eqs | lines 92–104 vs eqs deltaU/TU/Keta/muW-graph | match |
| (2b) `\mu_W^{\rm graph}` independent of `L` | lines 129–131 free_symbols check | match |
| (3) `\Pi^{\rm graph}=\Pi^{\rm can}` | Section III matrix `expect_zero` (lines 162–165) | match |
| (4) `E_T=q_{\rm tr}/(1+\chi_{0,*})`, `E_K=-q_\eta`, `E_\mu=...` | Section IV `expect_zero` ×3 (lines 187–192) + multiplicative ×3 (198–200) | match |
| (5) repair vector `(0,0,0,0,-E_K,0,-E_\mu,-E_T)` | Section V `expect_zero` (line 227) | match |
| (6) reduced-family packet vanishes on graph | Section VI `expect_zero` ×2 (lines 241–253) | match |
| (2) `\mathcal O_*={\mathbf x_*^{\rm graph}(\mathbf y)}` uniqueness | not separately asserted; the round-trip (II) + projection equality (III) establish containment + section consistency | partial (acceptable — uniqueness is the Stage-201/pivot result carried in, see notes §4.1) |
| — | no script check tests something the paper omits | no extra |

Dominant pattern is `match`; the single `partial` row is the uniqueness half of the graph theorem, which the notes (§4.1) explicitly route to the carried Stage-201 pivot block rather than re-proving in 202. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 115–118 | `log(Ctr[graph]/Ctr_tgt) == 0` | claim 1 (Ctr round-trip) | yes |
| A2 | sympy | 119–122 | `log(Cnt[graph]/Cnt_tgt) == 0` | claim 1 (Cnt round-trip) | yes |
| A3 | sympy | 123–126 | `log(eps_eta[graph]/epsEta_tgt) == 0` | claim 1 (eps round-trip) | yes |
| A4 | sympy | 130–131 | `L not in mu_graph.free_symbols` | claim 2b (L-independence) | yes |
| A5 | sympy | 162–165 | matrix `log(x_proj_can/x_graph) == 0` | claim 3 (graph=canonical proj) | yes |
| A6 | sympy | 187 | `E_T - q_tr/(1+chi0) == 0` | claim 4 | yes |
| A7 | sympy | 188 | `E_K + q_eta == 0` | claim 4 | yes |
| A8 | sympy | 189–192 | `E_mu - (q_nt - q_eta + F* q_tr/(1+chi0)) == 0` | claim 4 | yes |
| A9 | sympy | 198–200 | `E_T - log(m_T) == 0`, `E_K - log(m_K)`, `E_mu - log(m_mu)` | claim 4 (multiplicative chart) | yes |
| A10 | sympy | 227 | matrix `dx_rep - dx_rep_graph == 0` | claim 5 (repair vector) | yes |
| A11 | sympy | 241–244 | `Delta_family[graph,chiQ=1] == 0` | claim 6 | yes |
| A12 | sympy | 250–253 | `M_family[graph,chiQ=1] - (1,1,1,1) == 0` | claim 6 (multiplicative) | yes |

Every load-bearing assertion traces to a paper deliverable. No orphaned scaffolding. No assertion tests a claim the paper does not make.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 202 is not status-only (manifest `is_status_only_candidate: False`, checkpoint `False`) and performs substantive, fully symbolic algebra: solving three monomial constraints for a dependent triple, log-ratio identity verification, an 8-vector projection-equality, and a repair-vector rewrite. The stage card line 11 reads `Mathematica audit: none yet`, and no `.wl` exists anywhere (`find` confirms only the `.py`, output, and notes for stage 202; 197 sibling stages including the adjacent Stage 203 — which does the same class of symbolic algebra — do carry `.wl` scripts). Per the project dual-engine contract, a second-engine script is REQUIRED wherever Mathematica CAN independently verify the math. Every assertion in this stage (monomial round-trip, `E_T/E_K/E_\mu` identities, repair vector, family packet) is expressible in native Mathematica primitives (`Solve`/`Reduce`, `PowerExpand`, `Simplify`/`FullSimplify` under positive assumptions, `Log`, matrix ops). There is no impossibility carve-out, so the gap is a finding.

**Why this matters:**
The stage's exact-closure status rests on a single engine. The `\mu_W^{\rm graph}` `L`-independence claim, the `\Pi^{\rm graph}=\Pi^{\rm can}` equality, and the three quotient-equivalence identities are the load-bearing results that feed the Part VI realization compiler and finite mixed-ray search (downstream Stages 203–218). A SymPy-only verification leaves no cross-engine check that the log-ratio simplifications (which lean on `expand_power_base(..., force=True)` and `powsimp(..., force=True)`) did not silently mask a branch/sign error.

**Required change:**
Codex must author a NEW independent Mathematica script at the target path below. The directive states the requirement and claim manifest only; Codex designs and writes the implementation route. See directive F1.

**Verification:**
After Codex applies, `redteam exec-mathematica 202` runs the new `.wl`; it must exit 0 with each manifest check (M1–M7) printing a zero residual / `True`, and the new file must appear at `mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl`.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive's anti-transliteration guard (F1) requires the new script to derive the graph triple by an independent decomposition (e.g. `Solve`/`Reduce` of the monomial system, or `PowerExpand`+`Exponent` extraction) rather than echoing the SymPy hand-substitution `deltaU_graph = (Ctr_tgt/(...))**(1/(1+chi0))`.

## Engine cross-check

Only one engine present; not applicable. (`engines_agree: n/a`.)

## Verdict justification

The SymPy script holds up under attack against the paper. I checked: (a) constants — the exponents `1+\delta_{U,*}`, `1+\chi_{0,*}`, `E_*`, `-F_*` in lines 69–104 match the appendix monomials eq:app-part06-Ctr/Cnt-direct exactly, including the subtlety that `\delta_{U,*}` is carried as an independent symbol (`deltaU_star`), not as `T_U`-derived `delta_U`; (b) tautology — Section II solves the triple and then re-substitutes into the *original* monomials, a genuine round-trip, not `x==x`; Section III builds `x_proj_can` from monomial ratios independent of `x_graph`, so the equality is a real consistency test; (c) the graph-error identities (Section IV) genuinely exercise the paper's `E=q`-chart equivalence; (d) no `simplify` under physics-violating assumptions (all symbols `positive=True`, which is justified by the paper's "`\mathbf x>0`" target-orbit definition). The only attack that lands is the dual-engine gap: a substantive symbolic stage with no Mathematica second engine, which Mathematica can clearly verify. Verdict `findings` (one `missing_mathematica`); no `paper_misalignment`, no stop-cold.

## Self-test notes

Trap 1 (variable independence / `D[]`): the existing script has no derivatives, and the prescribed manifest has none — N/A. Trap 2 (integral parity): no integrals anywhere — N/A. Trap 3 (trivial-case): each manifest `assert_zero` is a `Log[ratio]` of a graph-reconstructed monomial over its target; substituting the graph solve collapses the ratio to 1 (residual 0), confirmed already by SymPy output lines 87–89 — the Mathematica re-derivation must reproduce these zeros. Trap 4 (paths): target `.wl` placed in `mathematica/` with the mirrored filename. Trap 5 (paper round-trip): the manifest uses only the paper's stated constants and introduces none, so no new `paper_misalignment` is created. Cosmetic note (not a finding): the script banners/comments say "STAGE 185 / Stage 201 / Stage 187" (internal-renumbering residue); no load-bearing identity differs, so this is not `paper_misalignment`.
