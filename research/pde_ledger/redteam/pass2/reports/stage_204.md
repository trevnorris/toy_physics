---
unit_id: 204
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
  notes_stage_files: [moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor.md]
  paper_appendix: present
---

# Audit unit 204 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_204.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 39, 236, 597–730, anchors 1304–1307)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Graph-invariant log-ray family, primitive free-direction table, scalar root equation \(\Phi_{\mathbf s}(\tau)=1\), monotone-ray uniqueness theorem, and first local root predictors." Stage 204 makes the abstract one-parameter family of Stage 203 explicit: it defines the free-quintuple log-ray \(\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s}\) with constant log-slopes \(s_i\); lifts it to the full graph point with four dependent exponents \(\sigma_\delta=-\mathfrak a_*(s_\gamma+s_c-s_U)\), \(\sigma_T=s_U+\sigma_\delta\), \(\sigma_{K_\eta}=2s_c-s_U\), \(\sigma_\mu=2s_c-s_U+2s_W-2s_\lambda-E_*(2s_\gamma+2s_\lambda-s_U-s_W)+F_*\sigma_\delta\), with \(\mathfrak a_*=(1+\delta_{U,*})/(1+\chi_{0,*})\); proves the three target monomials \(\mathfrak C_{\rm tr,*},\mathfrak C_{\rm nt,*},\epsilon_\eta\) are invariant along the lift for all \(\tau\) (equivalently \(M_*\dot{\Delta\mathbf x}=0\)); tabulates the five primitive-direction exponent rows; and gives the affine/log-linear root predictors \(\tau_{\rm aff}=(1-\Phi_0)/\Phi_1\), \(\tau_{\log}=-\ln\Phi_0/L_0\) with the first-order-agreement relation \(\tau_{\log}-\tau_{\rm aff}=-\varepsilon^2/(2L_0)+O(\varepsilon^3)\). The monotone-ray uniqueness theorem is a logical (sign-of-slope) statement, not a closed-form identity.

## What the script claims to verify

Both engines verify, in seven blocks: (I/M1) constant free log-slopes \(d\ln y_i/d\tau=s_i\); (II/M2,III/M3) the four dependent exponents — SymPy by equating a direct substitution against the posited exp-form \(\delta(0)e^{\sigma_\delta\tau}\), Mathematica by recovering the slope via \(\partial_\tau\ln(\cdot)\) and matching the posited formula plus confirming constancy; (III/M4) the three monomials reduce to their targets along the lift; (IV/M5) \(M_*\cdot\dot{\Delta\mathbf x}_{\rm ray}=0\); (V/M6) the five primitive-direction exponent rows; (VI/M7) first-order log-ray completeness via series; (VII/M7) the predictor difference \(\tau_{\log}-\tau_{\rm aff}\) starts at quadratic order with coefficient \(-1/(2L_0)\) (SymPy folds the eps^3 coefficient into one assertion; Mathematica checks coefficients 0..3 individually). The monotone-ray uniqueness theorem is not independently exercised by either script (it is a logical statement; the scripts verify only the differentiable scalar/predictor machinery it relies on).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Log-ray family / constant log-slopes | I (sympy 87–91), M1 (wl 86–90) | match |
| Four dependent exponents (\(\sigma_\delta,\sigma_T,\sigma_{K_\eta},\sigma_\mu\)) | II/III (sympy 109,116,123,138), M2/M3 (wl 100–128) | match |
| Finite monomial invariance | III (sympy 154–156), M4 (wl 156–158) | match |
| \(M_*\dot{\Delta\mathbf x}=0\) | IV (sympy 171), M5 (wl 167) | match |
| Primitive free-direction table | V (sympy 188–197), M6 (wl 178–181) | match |
| First-order log-ray completeness | VI (sympy 204), M7 (wl 195–196) | match |
| Root predictors \(\tau_{\rm aff},\tau_{\log}\) + quadratic-order agreement | VII (sympy 220–223), M7 (wl 190–193) | match |
| Scalar root equation \(\Phi_{\mathbf s}(\tau)=1\) | implicit (defines machinery; no symbolic \(\widehat\chi_Q\) available) | partial (legitimate — \(\widehat\chi_Q\) is PDE/candidate data, card flags this) |
| Monotone-ray uniqueness theorem | none (logical statement, not a closed-form identity) | partial (legitimate — not script-verifiable) |

`paper_alignment: aligned`. The two `partial` rows are expected: the card and appendix explicitly state \(\widehat\chi_Q\) values and derivatives are "PDE-returned or candidate-family data; the predictor formulas themselves are exact" (appendix line 729), and the uniqueness theorem is a sign-of-slope logical claim. The scripts verify exactly the exact-formula content; the non-symbolic / logical pieces are out of scope by the paper's own framing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 87–91 | `diff(log(y_i),tau) - s_i == 0` | constant log-slopes | yes |
| A2 | sympy | 109 | `delta_direct - delta_exp == 0` | sigma_delta | yes |
| A3 | sympy | 116/123/138 | `direct - exp == 0` (T,Keta,mu) | sigma_T/Keta/mu | yes |
| A4 | sympy | 154–156 | `monomial - target == 0` | invariance | yes |
| A5 | sympy | 171 | `M_* dx == 0` | quotient kernel | yes |
| A6 | sympy | 197 | `actual - expected == 0` (table) | primitive table | yes |
| A7 | sympy | 204 | `series - (Y0+Y1 eps) == 0` | first-order completeness | yes |
| A8 | sympy | 220–223 | `diff_series + eps^2/(2L0) - 2eps^3/(3L0) == 0` | predictor agreement | yes |
| M1 | wl | 86–90 | `logRate - s_i == 0` | constant log-slopes | yes |
| M2 | wl | 100–101 | `D[recovered,tau]==0` + `recovered - formula == 0` | sigma_delta | yes |
| M3 | wl | 123–128 | same for T,Keta,mu | sigma_T/Keta/mu | yes |
| M4 | wl | 156–158 | `monomial - target == 0` (via deltaFromTiming) | invariance | yes |
| M5 | wl | 167 | `mStar.dxRay == 0` | quotient kernel | yes |
| M6 | wl | 178–181 | `recovered/.subs - expected == 0` | primitive table | yes |
| M7 | wl | 190–196 | series coefficients 0..3 + completeness | predictor agreement / completeness | yes |

No tautological, hardcoded, or under-exercised assertions. Every exponent formula is checked against an independently-built quantity (direct substitution in SymPy; log-derivative recovery in Mathematica), so the equality can genuinely fail if a sign/factor were wrong.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.txt:3` (Date 2026-05-11T12:49:10) vs `.../scripts/...stage204_..._sympy_audit.py` (mtime 2026-06-03 15:59)

**What's wrong:**
The committed SymPy `.txt` output is older than the `.py` and was captured from a pre-renumber revision of the script. Two concrete symptoms: (a) the banner reads `STAGE 187 — …` (lines 11, 229) whereas the current `.py:35,225` prints `STAGE 204`; (b) the captured VII block shows `series[tau_log - tau_aff] = 2*eps^3/(3*L0) - eps^2/(2*L0)` (lines 221–225), which is consistent with the current assertion but is the only place the eps^3 term is visible — the current script's combined assertion on line 220–223 is not separately re-emitted. The Mathematica output (`mathematica/output/...204..._audit.txt`, mtime 2026-06-02) is fresh and correctly labelled STAGE 204.

**Why this matters:**
The captured SymPy transcript does not reflect the current script's identity (wrong stage number) and predates the Jun 3 edit, so it is not a faithful record. Content-wise the math shown still agrees with the current script (no numeric disagreement), so this is informational, not blocking.

**Required change:**
Re-run `python3 scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py` and overwrite the saved `.txt`. No source edit needed.

**Verification:**
New `.txt` header Date ≥ `.py` mtime; banner reads `STAGE 187` → `STAGE 204`; exit code 0; all checks still print `= 0`.

## Independent-derivation check (Mathematica)

Verdict: **INDEPENDENT.** The load-bearing objects are the four dependent exponents and the monomial invariance, and the two engines extract them by genuinely different methods.

1. **Exponents — derive-vs-posit/equate.** SymPy posits the closed form then proves it by equating two closed forms:
   ```
   # sympy 100, 109
   delta_graph_tau_exp = sp.simplify(delta_graph_0 * sp.exp(sigma_delta * tau))
   expect_zero("delta_U^graph direct - exp form", delta_graph_tau_direct - delta_graph_tau_exp)
   ```
   Mathematica instead RECOVERS the slope by differentiating the log of the graph quantity, and additionally proves constancy:
   ```
   # wl 99–101
   sigmaDeltaRecovered = logRate[deltaGraph[tau]];          (* = D[Log[deltaGraph[tau]], tau] *)
   expectZero["M2 d sigma_delta/dtau", D[sigmaDeltaRecovered, tau]];
   expectZero["M2 sigma_delta recovered - formula", sigmaDeltaRecovered - sigmaDelta];
   ```
   "build-direct-and-exp-form-then-subtract" vs "differentiate-the-log-to-recover-the-slope" are different extraction operations for the same exponent. Same for T/Keta/mu (sympy 113–138 vs wl 120–128).

2. **Monomial invariance — different reconstruction route.** SymPy rebuilds each monomial from the same `*_tau_direct` algebraic expressions (sympy 141–147). Mathematica instead defines abstract monomial functions and recovers \(\delta\) from the TIMING relation `deltaFromTiming[kU_,tU_] := Pi^2 tU/(L^2 kU)` (wl 131), applies them to the lifted ray vector `xGraphRay` whose dependent slots are `Exp[sigmaKeta tau]`, `Exp[sigmaMu tau]`, `Exp[sigmaT tau]` (wl 142–151) — i.e. it inverts the timing law rather than re-substituting the direct delta. Distinct choreography.

3. **Acknowledged shared/port-like sub-pieces (NOT load-bearing).** M5 reuses the same posited `M_*` matrix and `dx_ray` literal as the SymPy IV block (wl 161–166 ≈ sympy 159–166), and M7's predictor series both `Series`-expand the same `tau_log - tau_aff` closed form (wl 188 ≈ sympy 217). These are posit-vs-posit. But (a) the matrix is carried-forward Stage-192 data (sharing a definition is allowed), and its load-bearing inputs `sigmaKeta/sigmaMu/sigmaT` are independently recovered in M2/M3; (b) the predictor agreement is a secondary "first local predictor," not the stage's load-bearing exact closure object. Neither overturns the INDEPENDENT verdict on the exponents and invariance.

## Engine cross-check

Both engines pass all checks (sympy `.txt` and wl `.txt` both end "…AUDIT PASSED", exit 0). Concrete numeric/symbolic agreement on the one place a value is emitted: the predictor difference series. SymPy: `2*eps^3/(3*L0) - eps^2/(2*L0)` (sympy `.txt:221–225`). Mathematica: `-1/2*eps^2/L0 + (2*eps^3)/(3*L0) - (3*eps^4)/(4*L0)` (wl `.txt:83`). These agree through the orders both retain (eps^2 coeff `-1/(2L0)`, eps^3 coeff `+2/(3L0)`), and both match the paper's `-\varepsilon^2/(2L_0)+O(\varepsilon^3)` (notes §8.1; appendix line 670 family). `engines_agree: true`.

## Verdict justification

`findings` — one low-severity, informational `stale_output` (SymPy `.txt` predates the Jun 3 `.py` and carries the stale `STAGE 187` banner; content still agrees so it is non-blocking). I attacked: (i) the predictor self-test — hand-expanded `tau_log - tau_aff = (1/L0)[eps/(1+eps) - ln(1+eps)] = -eps^2/(2L0) + 2eps^3/(3L0) - 3eps^4/(4L0)`, which matches both the assertion and both engines, so it is non-tautological and correct; (ii) the exponent checks — each equates an independently-built quantity, so a wrong sign/factor would surface; (iii) transliteration — the exponents are recovered by log-derivative in Mathematica vs equate-closed-forms in SymPy, and the monomial route inverts the timing law rather than re-substituting, so the `.wl` is INDEPENDENT on the load-bearing objects despite sharing the `M_*` literal and the predictor closed form on two secondary checks; (iv) symbol domains — all base coordinates declared positive in both engines, matching the "positive base point \((\mathbb R_{>0})^5\)" of the notes; (v) paper alignment — every exact-formula deliverable maps to a check; the two unverified pieces (\(\widehat\chi_Q\) numeric values, monotone-ray uniqueness) are out of scope by the paper's own framing (appendix line 729). The stage holds up.

## Value Reconciliation (pass-2 augmentation)

The scripts are fully symbolic; their "results" are the derived exponent formulas, the invariance identities, and the predictor formulas — all of which are reported in the notes/appendix as boxed equations. No numeric figure-of-merit constants are emitted (only structural residuals = 0, which are scaffolding).

| value | source (py/wl + output line) | .tex / .md location | status |
|---|---|---|---|
| `aStar = (1+deltaU_star)/(1+chi0_star)` | py:57, wl:93 | notes line 167 (`\mathfrak a_*`); appendix eq line 607 | MATCH |
| `sigma_delta = -aStar*(s_gamma+s_c-s_U)` | py:58, wl:98; sympy.txt:48–51 | notes line 173; appendix line 611 | MATCH |
| `sigma_T = s_U + sigma_delta` | py:59, wl:112; sympy.txt:83–86 | notes line 181; appendix line 615 | MATCH |
| `sigma_Keta = 2 s_c - s_U` | py:60, wl:113; sympy.txt:88–89 | notes line 184; appendix line 617 | MATCH |
| `sigma_mu = 2s_c - s_U + 2s_W - 2s_lam - E_*(2s_gam+2s_lam-s_U-s_W) + F_* sigma_delta` | py:61–68, wl:114–118; sympy.txt:91–98 | notes lines 187–193; appendix lines 621–624 | MATCH |
| Monomial invariance `Ctr/Cnt/eps_eta = target` for all tau | py:154–156, wl:156–158 | notes lines 269–287; appendix line 641 | MATCH |
| `M_* dot(Delta x)_ray = 0` | py:171, wl:167 | notes line 302 | MATCH |
| Primitive table e_lambda `(0,0,0,-2-2E_*)` | py:182, wl:172; sympy.txt:179 | notes table line 334; appendix table (rows derive from σ formulas) | MATCH |
| Primitive table e_c `(-a,-a,2,2-F_* a)` | py:183, wl:173 | notes table line 335 | MATCH |
| Primitive table e_gamma `(-a,-a,0,-2E_*-F_* a)` | py:184, wl:174 | notes table line 336 | MATCH |
| Primitive table e_U `(a,1+a,-1,-1+E_*+F_* a)` | py:185, wl:175 | notes table line 337 | MATCH |
| Primitive table e_W `(0,0,0,2+E_*)` | py:186, wl:176 | notes table line 338 | MATCH |
| `tau_aff = (1-Phi0)/Phi1` | py:210, wl:186 | notes line 460; appendix line 670 | MATCH |
| `tau_log = -ln(Phi0)/L0` | py:211, wl:187 | notes line 466; appendix line 672 | MATCH |
| `tau_log - tau_aff = -eps^2/(2L0) + O(eps^3)` | py:222 / sympy.txt:221–225, wl.txt:83 | notes line 483; appendix line 670 family | MATCH |
| eps^3 coefficient `+2/(3L0)` (script refinement) | py:222, wl:193, wl.txt:83 | notes use `O(eps^3)` (line 484) → subsumes; not separately stated | MATCH (within `O(eps^3)`) |

Internal scaffolding (no prose expected): banner/subbanner strings, `expect_zero`/`expectZero`/`logRate`/`normalizeScalar`/`pretty`/`pass`/`fail` helpers, all `... = 0` residual prints, the intermediate `delta_graph_*`, `T_graph_*`, `Keta_graph_*`, `mu_graph_*`, `*_tau_direct`/`*_exp` working expressions, `deltaFromTiming`, `ray_component`/`localRayFirstOrder` series workspaces, `Phi0_eps`/`Phi1` placeholder symbols.

reconciliation: complete; 16 deliverable values checked, 0 misaligned.

## Self-test notes

Variable independence: every `sp.diff`/`D[...,tau]` differentiates an expression that genuinely depends on `tau` (each ray component carries `Exp[s_i tau]`), so no identically-zero-derivative trap; the `D[sigmaXRecovered, tau] == 0` constancy checks are nontrivial because `sigmaXRecovered` is the τ-derivative of a log that, pre-simplification, still contains τ. Trivial-case: substituting `s_gamma=s_c=s_U` gives `sigma_delta=0` consistently in both formula and recovery. Predictor parity/series: hand-expanded `tau_log-tau_aff` and confirmed the eps^2 coefficient `-1/(2L0)` and eps^3 `+2/(3L0)` match both engines and the paper. No `## Resolve before fix_loop` blocks needed (the lone finding is a script-side stale-output, not a paper_misalignment).
