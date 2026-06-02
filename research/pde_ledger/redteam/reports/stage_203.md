---
unit_id: 203
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md]
  paper_appendix: present
---

# Audit unit 203 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_203.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (read intro + Stage 203 rows: §sec:app-part06-scalar-graph-slice-and-log-ray-compiler, eqs eq:app-part06-chi-hat-def … eq:app-part06-ray-closure; Theorem items 5–6)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Scalar graph slice \(\mathcal Z_*=\{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\chi_Q(\mathbf y)=1\}\), graph-family tangency, and one-parameter crossing theorem." The card (`\stagefield{Derivation ledger}`) and the Part VI appendix (§sec:app-part06-scalar-graph-slice-and-log-ray-compiler; Theorem item 5–6) name four distinct deliverables: (D1) the graph-closure scalar `\widehat\chi_Q(\mathbf y):=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y))` and the codimension-one slice `\mathcal Z_*`; (D2) the exact dependent-triple graph log-tangent formulas (notes §3); (D3) the kernel theorem `M_*\dot{\Delta\mathbf x}_{\rm graph}=0`, i.e. graph families lie in `\ker M_*` (notes §4); and (D4) the same-free-quintuple candidate decomposition and graph-error→quotient compiler `q_{\rm tr}=(1+\chi_{0,*})E_T`, `q_{\rm nt}=E_\mu-E_K-F_*E_T`, `q_\eta=-E_K`, with the inverse compiler (notes §5). The headline existence result is (D5) the one-parameter crossing theorem (notes §7, appendix eq:app-part06-IVT-crossing-condition/result): if `\widehat\Delta_Q(\tau_-)\widehat\Delta_Q(\tau_+)<0` then there is `\tau_*` with `\widehat\Delta_Q(\tau_*)=0`. The appendix additionally states (eq:app-part06-graph-lifted-ray) that the three target monomials are *invariant for all `\tau` on the graph lift*, which is the structural reason `\widehat\chi_Q` is a faithful scalar reduction. Status tag: `\StatusExactClosure{}`; primary anchor MTDC-T10.1.

Note (provenance, not a math defect): the notes body refers to this stage as "Stage 237" and to upstream sources as Stage 253/243/192/197, while the card/appendix consistently use 203/202/201 and the MTDC-T10.1 anchor. The scripts likewise banner "STAGE 186" and cite Stage 192/197/202/243. These are stale stage numbers from a renumbering; the algebra they carry is the same. Recorded as informational below; no number affects a verified identity.

## What the script claims to verify

Both scripts share six sections. §I builds the carried Stage 243 monomial-drift matrix `M_*` and verifies its canonical dependent-triple section `M_* S - I_3 = 0` (D3 scaffolding). §II derives the dependent graph log-tangents by differentiating the graph map and checks each against the carried closed forms (D2). §III assembles `dx_graph` from those tangents and verifies `M_* dx_graph = 0` (D3). §IV builds the three target monomials from the graph lift dressed by error exponents `E_T,E_K,E_mu` and verifies the graph-error→quotient relations (D4). §V solves the inverse compiler and verifies the repair vector `-S q = (0,0,0,0,-E_K,0,-E_mu,-E_T)` and `M_*\Delta x_{\rm rep}+q=0` (D4). §VI ("Stage 197 scalar closure composed with an explicit Stage 202 graph path") substitutes `\Sigma_0=0,\Sigma_5=0` into the carried Stage 197 closure `\chi_Q=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)`, reducing it to `\beta^5`, builds a path `\beta(\tau)=1+\rho(2\tau-1)/(1+\rho)`, and demonstrates the IVT crossing on `\widehat\Delta_Q=\beta^5-1` (D5).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 `\widehat\chi_Q`, slice `\mathcal Z_*` | §VI builds `hat_chi_graph`, `hat_delta_graph` (but as `\beta^5`, see F1) | partial |
| D2 dependent graph log-tangents | §II `dln delta`, `tau_1`, `kappa_eta`, `mu_1` vs carried forms | match |
| D3 kernel `M_*\dot{\Delta x}_{\rm graph}=0` | §III `M_* dot(Delta x)_graph` + §I `M_* S - I_3` | match |
| D4 graph-error→quotient + inverse compiler | §IV `q_tr/q_nt/q_eta`; §V inverse + repair vector | match |
| D5 one-parameter crossing theorem | §VI IVT on `\beta^5-1` over a γ-only path; four free coords inert, `\Sigma_0=\Sigma_5=0` hardcoded | partial |
| (appendix) three target monomials invariant on the graph lift | none — bypassed by setting `\Sigma_0=\Sigma_5=0` | missing |

Dominant pattern: D2–D4 match faithfully; D1 and D5 are only partially exercised, and the appendix's monomial-invariance claim has no check. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 106 | `expect_zero(M_* S - I_3)` | D3 scaffolding | yes |
| A2 | sympy | 158–178 | `expect_zero` ×4 graph log-tangents | D2 | yes |
| A3 | sympy | 185 | `expect_zero(M_* dx_graph)` | D3 | yes |
| A4 | sympy | 235–237 | `expect_zero` ×3 graph-error→quotient | D4 | yes |
| A5 | sympy | 257–271 | inverse compiler + repair vector + `M_* dx_rep + q` | D4 | yes |
| A6 | sympy | 307–311 | `beta_path` id + closure-numerator id | D5 (path setup) | partial |
| A7 | sympy | 312 | `expect_positive(denominator)` | D5 (continuity of `\widehat\Delta_Q`) | yes |
| A8 | sympy | 326–330 | `<0` at τ=0, `>0` at τ=1, `=0` at τ=1/2, unique real root | D5 (IVT crossing) | partial (object is `\beta^5`, not `\chi_Q\circ x^{graph}`) |
| A9 | math | 137 | `expectZero[M_* S - I_3]` | D3 scaffolding | yes |
| A10 | math | 169–185 | `expectZero` ×4 graph log-tangents (log-space) | D2 | yes |
| A11 | math | 192 | `expectZero[M_* dxGraph]` | D3 | yes |
| A12 | math | 228–230 | `expectZero` ×3 graph-error→quotient (log-space) | D4 | yes |
| A13 | math | 247–261 | inverse compiler + repair vector + kernel | D4 | yes |
| A14 | math | 283–306 | `beta_path` id, closure-numerator id, denom>0, signs, root, unique crossing | D5 | partial (object is `\beta^5`) |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py:273-330`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl:263-306`

**What's wrong:**
The headline deliverable D5 (one-parameter crossing theorem) and the supporting D1 object `\widehat\chi_Q(\mathbf y):=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y))` are exercised on a hand-reduced scalar, not on the paper's actual object. Section VI's own header reads "Stage 197 scalar closure composed with an explicit Stage 202 graph path," and the appendix states (eq:app-part06-graph-lifted-ray): "The three target monomials are invariant for all `\tau` on this lift. Hence the scalarized closure function `\Phi_{\mathbf s}(\tau):=\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau))` contains the entire graph-aligned realization question." But the script does this:

sympy lines 286, 289–291:
```
chi_from_stage197 = sp.simplify(3 * (Siso * beta**5 + 9 * Sigma5) / (3 * Siso - Sigma0))
hat_chi_graph = sp.simplify(chi_from_stage197.subs({beta: beta_path, Sigma0: 0, Sigma5: 0}))
```
mathematica lines 269, 272 do the same with `Sigma0 -> 0, Sigma5 -> 0`. So `\widehat\chi_Q` collapses to `\beta(\tau)^5`, and the free-quintuple path `y_graph_path = [lam_bar, cetaU_bar, gamma_path, KU_bar, KW_bar]` moves only `\gamma` (through `beta_path`); `lam_bar, cetaU_bar, KU_bar, KW_bar` are declared (sympy:278–280) but never enter `\widehat\chi_Q`. Two distinct shortcuts are taken without justification in-script:
1. `\Sigma_0=\Sigma_5=0` is asserted by hand as "on the graph," but it is not derived from `\mathbf x_*^{\rm graph}`; the appendix's stated reason the scalar reduces is *monomial invariance on the lift*, which is never checked.
2. `\chi_Q` is never composed with the Stage 202 graph map built in §IV (`Ctr, Cnt, epsEta` from `T_graph, Keta_graph, mu_graph`); the section reuses a bare upstream `\beta^5` instead, so D1's defining composition is not exercised at the slice level.

The IVT mechanics themselves (continuity via positive denominator, sign change `\Delta(0)<0`, `\Delta(1)>0`, root at τ=1/2, uniqueness) are internally correct and non-tautological — the defect is the *object* the theorem is run on, which is narrower than the paper's `\widehat\chi_Q` over the five-coordinate graph.

**Why this matters:**
On a checkpoint, the headline crossing theorem and the slice scalar must be exercised on the paper's actual object. As written, the section would still PASS if the Stage 202 graph map had a wrong dependent-coordinate exponent, because the graph map plays no role in §VI; the crossing is verified on a generic `\beta^5` toy. The appendix's load-bearing structural claim — that the three target monomials are invariant along the graph lift, which is what makes closure a *single* scalar — is not verified anywhere in the unit, even though it is directly checkable (it is the global/finite counterpart of §III's kernel identity).

**Required change:**
See directive F1. Requirement (claim-manifest, route left to Codex): within §VI, exercise the crossing/slice scalar on the genuine graph lift rather than a hand-zeroed `\beta^5`. Specifically: (i) add a check that the three target monomials evaluated on `\mathbf x_*^{\rm graph}(\mathbf y(\tau))` are identically equal to their target values (`q_tr=q_nt=q_eta=0`) along a non-trivial free-quintuple path that moves at least two of the five free coordinates — i.e. monomial invariance on the lift; and (ii) drive the crossing demonstration from `\chi_Q` *evaluated through the graph composition* (so the `\Sigma_0=\Sigma_5=0` reduction, if retained, is shown to follow from the lift rather than imposed by fiat). The IVT sign/root/uniqueness checks may remain. Anti-transliteration guard: the Mathematica check must form the monomial-invariance residual in native log-space (`Log`/`PowerExpand`/`Simplify`) independently of the SymPy power-space construction.

**Verification:**
A new `expect_zero`/`expectZero` block should appear in §VI of both scripts asserting the three target-monomial residuals vanish on the multi-coordinate graph path, and the crossing scalar should be sourced from the graph composition. Both scripts must exit 0 with the new checks passing; outputs refreshed.

## Independent-derivation check (Mathematica)

Genuinely independent for the load-bearing graph algebra. §II/§IV are derived in a *different representation* than SymPy: SymPy constructs `deltaU_graph_t` and the monomials as multiplicative power expressions and takes `sp.log(...)` (sympy:128–147, 213–226); Mathematica constructs the same quantities *directly in additive log-space* and differentiates/compares there:
- sympy:129–133 `deltaU_graph_t = (Ctr_tgt/(gamma_t*cetaU_t/KU_t)**(1+deltaUs))**(1/(1+chi0s))`; `T_graph_t = L**2*KU_t*deltaU_graph_t/pi**2` → then `sp.diff(sp.log(T_graph_t), t)`.
- wl:148–151 `logDeltaGraphT = (Log[CtrTgt] - (1+deltaUs)(Log[gammaT]+Log[cetaUT]-Log[KUT]))/(1+chi0s)`; `logTGraphT = 2 Log[L] + Log[KUT] + logDeltaGraphT - 2 Log[Pi]` → then `D[logTGraphT, t]`.
This is not a syntactic port; it sidesteps `PowerExpand` branch traps by working in logs from the outset, which is a legitimate second-engine check. Similarly §IV: sympy builds `Ctr` as a product of powers then `sp.log(Ctr/Ctr_tgt)` (sympy:213–224); wl builds `qtr` as a direct log-linear combination (wl:211–215).
Sections I, III, V, VI are structurally parallel (same `Mstar`, same pivot columns `{8,5,7}`, same `Solve` system, same `beta_path`, same `Sigma0->0,Sigma5->0`, same τ∈{0,1,1/2}). That parallelism is acceptable: §I/§III/§V are the same finite linear algebra any engine performs, and §VI shares the F1 defect in both engines rather than being an independent cross-check. On balance, **independent** — no `mathematica_transliteration` finding — but the F1 shortcut is duplicated in both engines, so the second engine does not currently catch it.

## Engine cross-check

Engines agree. Final forms match across both outputs:
- `\widehat\Delta_Q(0)` sympy `-rho*(rho^4+5rho^3+10rho^2+10rho+5)/(rho+1)^5` ↔ wl `-(rho*(5+10rho+10rho^2+5rho^3+rho^4))/(1+rho)^5` (identical).
- `\widehat\Delta_Q(1)` sympy `rho*(31rho^4+75rho^3+70rho^2+30rho+5)/(rho+1)^5` ↔ wl `(rho*(5+30rho+70rho^2+75rho^3+31rho^4))/(1+rho)^5` (identical).
- real crossing set: sympy `{1/2}` ↔ wl `2 tau == 1` (identical).
- `M_* S - I_3 = 0`, `M_* dx_graph = 0`, all §II/§IV/§V residuals = 0 in both. No `engine_disagreement`.

## Verdict justification

`findings`. Deliverables D2 (graph log-tangents), D3 (kernel identity + canonical section), and D4 (graph-error→quotient compiler + inverse + repair vector) hold up cleanly: every assertion is non-tautological (a wrong graph exponent or wrong matrix entry would surface a nonzero residual), well-anchored, and the two engines derive the graph algebra in genuinely different representations and agree. Attacks that failed: tried to find a tautology in §IV (the `q` residuals are built from independent monomial definitions, not the error-exponent expressions, so they can fail); tried to find a domain error in the `expect_positive`/`expect_negative` sign checks (they genuinely depend on `rho>0`, declared); tried to find an engine residual disagreement (none). The one substantive gap is the headline crossing theorem / slice scalar (D1, D5): §VI runs the IVT on a hand-reduced `\beta^5` with four free coordinates inert and `\Sigma_0=\Sigma_5=0` imposed by fiat, and never verifies the appendix's stated reason the scalar reduction is legitimate (target-monomial invariance on the graph lift). This clears the substantive bar for D2–D4 but not for D1/D5 at checkpoint strictness, hence `checkpoint_bar: concerns` and `paper_alignment: partial`. No `paper_misalignment` in the math sense — every formula the script *does* verify matches the card/appendix/notes exactly; the issue is coverage, not a wrong identity. No stop-cold: the fix is additive (a new check), does not change any derived constant, and does not propagate downstream.

## Self-test notes

Checked: (1) variable independence — proposed F1 monomial-invariance residual uses `Ctr,Cnt,epsEta` evaluated on `x_*^graph`, which by construction equal the target values; the check has no `diff`, so no identically-zero-derivative trap; the path must move ≥2 of the five free coordinates so the cancellation is non-trivial (a single-coordinate path could mask a wrong exponent in an untouched coordinate). (2) Parity — N/A, no unbounded integrals. (3) Trivial-case — substituting the existing γ-only `beta_path` already yields `q_tr=q_nt=q_eta=0` because §IV proves the graph monomials hit their targets; confirmed the residual reduces to 0. (5) Paper round-trip — the F1 fix reuses the same target constants `Ctr_target, Cnt_target, epsEta_target` the paper states (appendix eq:app-part06-target-values) and introduces no new constant, so it cannot create a new paper_misalignment.
