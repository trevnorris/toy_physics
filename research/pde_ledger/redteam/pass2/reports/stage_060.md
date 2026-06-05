---
unit_id: 060
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage060_entropic_microclosure.md]
  paper_appendix: present
---

# Audit unit 060 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_060.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage060_entropic_microclosure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row line 98; `\input` line 238)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.txt`

## What the paper claims

Stage 060 (`\stagefield{Output}`: "A positive-density microclosure for the transport source family and the definition of the microscopic gain variable") shows the Stage-058 drift-diffusion source branch is *derivable* from a positive local source thermodynamics rather than postulated as a shape family. The notes enumerate the deliverables: (1) an explicit positive-density free energy `F[sigma,phi]`; (2) its Euler-Lagrange variations — chemical potential `mu_sigma = Theta_sigma log(sigma/sigma_*) - Lambda_phi phi` and support bulk equation `-T_X phi_ss + K_X phi = Lambda_phi sigma` with Robin/Neumann BCs; (3) the minimal Onsager current `J = -D_sigma sigma_s + M_sigma Lambda_phi sigma phi_s` with Einstein relation `D_sigma = M_sigma Theta_sigma`; (4) recovery of the Stage-056 exponential family `Sigma_Pe(x) = Pe e^{Pe x}/(e^{Pe}-1)` under the affine-drop reduction, with `Pe = (Lambda_phi/Theta_sigma) Delta phi`; (5) the exact microscopic coupling `Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X) = chi_sigma Lambda_phi^2 L^2/T_X = mu_sigma Lambda^2 L^2/(D_sigma T_X)`; (6) the passivity / Lyapunov identity `dF/dt = -int J^2/(M_sigma sigma) <= 0` under no-flux boundaries. The appendix row (line 98) summarizes it as a "Positive source/free-energy closure reproducing the transport source family."

## What the script claims to verify

Both scripts verify, by symbolic differentiation and simplification: (1) `mu = Theta log(sigma/sigma_*) - Lambda phi` as `delta F/delta sigma`; (2) the Onsager current decomposition into drift-diffusion form and the `D_sigma = M_sigma Theta` substitution; (3) that the exponential trial `C exp(a s)` with `a = Lambda Delta_phi/(Theta L)` zeroes the affine-drop current, that the normalization constant `C` integrates the family to 1 on `[0,L]`, that the rescaled `Sigma_Pe(x) = Pe e^{Pe x}/(e^{Pe}-1)` integrates to 1 on `[0,1]`, and that the rate `gamma` solved independently from `J=0` equals `Lambda Delta_phi/(Theta L)`; (3.5) a support BVP in the `K_X→0`, `K_m→∞` rigid-grounding limit yields end-to-end drop `Lambda L^2 sigma_0/(2 T_X)`, anchoring the `Lambda L^2/T_X` scale; (4) `Xi_micro = Lambda^2 L^2/(Theta T_X)` and its susceptibility and phenomenological (D/M) forms; (5) the dissipation density `M sigma mu_s^2 >= 0` is provably nonnegative (supporting `dF/dt <= 0`, the integral identity itself only printed as a calculus remark); (6) the full-free-energy support Euler-Lagrange bulk equation.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (2) `mu = Theta log(sigma/sigma_*) - Lambda phi` | py 42-43 / wl 48 | match |
| (2) support EL `-T_X phi_ss + K_X phi = Lambda sigma` | py 150-153 / wl 162-165 | match |
| (3) Onsager current + Einstein `D=M Theta` | py 46-52 / wl 52-62 | match |
| (4) exponential family `Sigma_Pe(x)` + normalization | py 56-78 / wl 64-96 | match |
| (4) `Pe = (Lambda/Theta) Delta phi` (rate identification) | py 80-89 / wl 98-107 | match |
| (5) `Xi_micro = Lambda^2 L^2/(Theta T_X)` (+ chi, D/M forms) | py 122-127 / wl 132-142 | match (wl headline posited, not derived — see F2) |
| (6) passivity `dF/dt <= 0` | py 135-145 / wl 145-156 (density positivity only; integral identity printed) | partial-by-design (acceptable) |

Paper alignment is `aligned`: every stated deliverable maps to a script-side check, and all numeric/symbolic forms agree with the notes. The card itself is terse (names the gain `G_micro`, defers the explicit form to the notes), which is legitimate. No `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero(mu - [...])` | (2) chem pot | yes (mu derived via `diff(f,sigma)`) |
| A2 | sympy | 48 | `expect_zero(J decomposition)` | (3) Onsager | yes (J from `diff(mu,s)`) |
| A3 | sympy | 51 | `expect_zero(D_sigma subst)` | (3) Einstein | yes |
| A4 | sympy | 62 | `expect_zero(J=0 by exp family)` | (4) exp branch | yes |
| A5 | sympy | 68 | `expect_zero(C normalizes)` | (4) normalization | yes |
| A6 | sympy | 76 | `expect_zero(Sigma_Pe rescale)` | (4) Sigma_Pe | yes |
| A7 | sympy | 77 | `expect_zero(Sigma_Pe integrates to 1)` | (4) normalization | yes |
| A8 | sympy | 88 | `expect_zero(gamma_derived - Lambda dphi/(Theta L))` | (4) Pe rate | yes (gamma from `solve(J=0)`) |
| A9 | sympy | 115 | `expect_zero(BVP rigid limit - L^2 Lam sig0/(2 T_X))` | scale anchor for (5) | yes |
| A10 | sympy | 125 | `expect_zero(Xi_micro - Lam^2 L^2/(Theta T_X))` | (5) gain | partial (Delta cancellation; weak but does exercise division) |
| A11 | sympy | 126-127 | `expect_zero(chi & D/M forms)` | (5) equivalences | yes |
| A12 | sympy | 141 | `assert ask(nonnegative(M sigma mu_s^2))` | (6) passivity | yes |
| A13 | sympy | 153 | `expect_zero(support bulk eq form)` | (2) EL | yes |
| B1 | math | 48 | `expectZero(mu - [...])` | (2) | yes |
| B2 | math | 54-62 | `expectZero(J decomp, D subst)` | (3) | yes |
| B3 | math | 75-96 | `expectZero(J=0, C norm, Sigma_Pe, integral)` | (4) | yes |
| B4 | math | 106 | `expectZero(gammaDerived - ...)` | (4) | yes (gamma from `Solve`) |
| B5 | math | 126 | `expectZero(BVP rigid limit)` | scale anchor for (5) | yes |
| B6 | math | 134, 137 | `expectZero(xiMicro via chi / via D/M)` | (5) equivalences | yes |
| B7 | math | 140 | `expectZero(xiMicro - lambdaPhi^2 ell^2/(theta tX))` | (5) gain headline | **no — tautological (F2)** |
| B8 | math | 141-142 | `expectZero(chi & D/M forms)` | (5) | yes |
| B9 | math | 151-156 | `Reduce[ForAll ... dissipationDensity>=0]` | (6) | yes |
| B10 | math | 165 | `expectZero(support bulk eq form)` | (2) | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.txt` (mtime 2026-05-22 17:59:13)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.txt` (mtime 2026-05-22 17:59:27)
- vs `scripts/...sympy_audit.py` (mtime 2026-06-03 15:59:11) and `mathematica/...mathematica_audit.wl` (mtime 2026-06-03 15:59:11)

**What's wrong:**
Both saved `.txt` outputs predate their scripts by ~12 days. The content is also visibly stale on numbering self-labels: the SymPy output banner reads `STAGE 43 — ENTROPIC SOURCE MICROCLOSURE` and `STAGE 43 THEOREM LEDGER` (txt lines 3, 37) and ledger item 3 reads "Stage-39 exponential family" (txt line 41), while the Mathematica output banner reads `STAGE 043 — ENTROPIC MICROCLOSURE` (txt line 4). The current SymPy source uses `banner("STAGE 60 ...")` (py:22) and `banner("STAGE 60 THEOREM LEDGER")` (py:156); the current Mathematica source uses `banner["STAGE 060 — ENTROPIC MICROCLOSURE"]` (wl:32). The numeric residuals/forms in the transcripts match what the current scripts would emit, so the staleness is confined to the banner/label text — but the transcript was nonetheless captured before the current script revision.

Related script-side stale numbering self-labels (low-severity `stale_output` class, flagged per the orchestrator's note; orchestrator resolves label scope):
- `scripts/...sympy_audit.py:3` docstring "Stage 43 SymPy audit" (should be Stage 60).
- `scripts/...sympy_audit.py:159` "exactly the Stage-39 exponential family" (canonical is Stage-056; cf. `.tex:21` "Stage~056" and notes:50 "Stage-056"; 39+17=56 offset).
- `scripts/...sympy_audit.py:22,156` banner text uses 2-digit "STAGE 60" (cosmetic vs 3-digit "STAGE 060" used by the `.wl`).

**Why this matters:**
The committed transcript is the authoritative record per the audit protocol; a stale transcript whose headers name a different stage misleads future readers and the reconciliation step. The residuals are unaffected, so this is informational, but a refresh is warranted and the SymPy docstring/ledger self-labels should be corrected.

**Required change:**
Correct the SymPy self-labels (py:3 "Stage 43"→"Stage 60"; py:159 "Stage-39"→"Stage-56"), then re-run both scripts to refresh the committed `.txt` outputs so the banners read STAGE 060 consistently. (Banner 2-digit→3-digit padding at py:22,156 is cosmetic — optional, orchestrator's call.)

**Verification:**
After refresh, both `.txt` mtimes postdate their scripts; the SymPy output banner reads "STAGE 60"/"STAGE 060" (no "43"); ledger line 3 reads "Stage-56"; all PASS lines remain.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:132,140`

**What's wrong:**
On the Mathematica side the headline gain is *posited*, not derived, and then "verified" against itself. Line 132 sets `xiMicro = lambdaPhi^2*ell^2/(theta*tX);`. Line 140 then asserts:
```
expectZero["Xi_micro - Lambda^2 L^2/(Theta T_X)", xiMicro - lambdaPhi^2*ell^2/(theta*tX)];
```
This is `lambdaPhi^2*ell^2/(theta*tX) - lambdaPhi^2*ell^2/(theta*tX)`, i.e. `expr - expr`, which is identically zero by construction and cannot fail no matter what the physics is. Contrast the SymPy side, where `Xi_micro` is *derived* by cancellation (`Xi_micro = simplify(Lam * phi_from_Phi / Theta / Delta)` with `phi_from_Phi = Lam*L**2*Delta/T_X`, py:122-123), so its matching check at py:125 at least exercises the `Delta` cancellation. The Mathematica script never forms `Xi_micro` from `phi_from_Phi/Delta`; the comment at wl:129-131 states the `deltaSupport cancellation is avoided`, so the gain is just declared.

The two genuinely-content-bearing Mathematica checks for this deliverable are the chi-substitution (wl:134-135) and D/M-substitution (wl:137-138) consistency checks, plus the susceptibility/phenomenological forms (wl:141-142). Line 140 adds nothing — it is dead scaffolding.

**Why this matters:**
On its own, line 140 contributes a false sense of independent confirmation of the headline value `Xi_micro = Lambda^2 L^2/(Theta T_X)` on the Mathematica engine, when in fact that value is only posited there. The deliverable IS still covered on the Mathematica side by B6/B8 (the chi and D/M consistency checks) and is independently derived on the SymPy side, so the overall claim is not unverified — but the specific assertion at line 140 is tautological and should either be removed or replaced with a check that ties `xiMicro` to the `Pe`/`Delta phi` route already present in the script.

**Required change:**
At `wl:140`, replace the self-referential assertion with one that derives the headline from quantities the script already computes rather than re-stating the literal. Concretely, build `Xi_micro` from the Peclet identification and the BVP-anchored support drop the script derives, e.g. set
`xiMicroDerived = FullSimplify[lambdaPhi/theta * (lambdaPhi*ell^2/tX), Assumptions -> $Assumptions];`
(the `Pe = (lambdaPhi/theta) Delta phi` rate times the `lambdaPhi ell^2/tX` per-unit-Delta support drop already used in the notes/SymPy route) and assert `expectZero["Xi_micro from Pe x drop", xiMicroDerived - xiMicro];`. This makes the headline value the *result* of a substitution chain instead of a restatement of itself. If a faithful in-script derivation is not practical, at minimum delete line 140 (the chi/D-M consistency checks at wl:134-142 already cover the deliverable non-tautologically).

**Verification:**
After the fix, `wl:140` no longer reads `xiMicro - lambdaPhi^2*ell^2/(theta*tX)`; the new/remaining `Xi_micro` checks each reference an expression formed from an independent route (Pe rate × support drop, or the chi/D-M substitutions), and the Mathematica script still exits 0 with all PASS lines. The refreshed Mathematica output reflects the new check label.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT a line-by-line transliteration. Three concrete divergences: (a) it forms `xiMicro` directly as a dimensional combination and explicitly avoids the SymPy `phi_from_Phi/Delta` cancellation (wl:129-132 vs py:122-123); (b) it proves dissipation positivity with `Reduce[ForAll[{muS,sigmaVal,mSigma}, Implies[...,>=0]]]` (wl:151-153), whereas SymPy uses `sp.ask(sp.Q.nonnegative(...))` (py:141) — genuinely different decision procedures; (c) the support BVP uses `DSolveValue` + `Solve` with `Cases[..., _C, Infinity]` constant extraction (wl:113-118) vs SymPy `dsolve` + free-symbol sorting (py:96-105). The choreography is parallel (same physics, same stage structure) but the routes to the load-bearing results differ. No `mathematica_transliteration` finding. The one place this independence backfires is exactly F2: by avoiding the cancellation route, the Mathematica `Xi_micro` ends up posited and its self-check tautological.

## Engine cross-check

The two engines agree on every load-bearing result. Side-by-side of the headline outputs:

- `mu_sigma^(chem)`: sympy `-Lambda_phi*phi(s) + Theta*log(sigma(s)/sigma_*)` (txt:6) ↔ math `-(lambdaPhi*phi) + theta*Log[sigma/sigmaStar]` (txt:6) — agree.
- `J expanded`: sympy `Lambda_phi*M_sigma*sigma(s)*phi'(s) - M_sigma*Theta*sigma'(s)` (txt:8) ↔ math `lambdaPhi*mSigma*sigma[s]*phi'[s] - mSigma*theta*sigma'[s]` (txt:9) — agree.
- `gamma_derived`: sympy `Delta_phi*Lambda_phi/(L*Theta)` (txt:19) ↔ math `(deltaDrop*lambdaPhi)/(ell*theta)` (txt:22) — agree.
- `Delta_derived` (BVP): sympy `L**2*Lambda_phi*sigma_0/(2*T_X)` (txt:23) ↔ math `(ell^2*lambdaPhi*sigma0)/(2*tX)` (txt:25) — agree.
- `Xi_micro`: sympy `L**2*Lambda_phi**2/(T_X*Theta)` (txt:25) ↔ math `(ell^2*lambdaPhi^2)/(theta*tX)` (txt:32) — agree.
- `dissipation density`: sympy `M_sigma*mu_s**2*sigma_val` (txt:29) ↔ math `mSigma*muS^2*sigmaVal` (txt:39) — agree.

All `expect_zero`/`expectZero` residuals are 0 in both transcripts. `engines_agree: true`.

## Verdict justification

Verdict is `findings` (2, both low-severity, both script-side, `material_change: false`). The physics is sound and fully aligned with the paper: every Stage-060 deliverable in the notes maps to a non-trivial check in at least one engine, both engines agree on all residuals, and no value disagrees with the notes/card/appendix. Attacks that failed: I tried to break the chemical-potential, Onsager, exponential-family, normalization, Pe-rate, BVP-scale, susceptibility/D-M, and EL checks for tautology or hidden assumptions — all are anchored to derived (differentiated/solved) quantities and could genuinely fail. The two findings are (F1) stale committed transcripts with leftover "Stage 43/043/39" self-labels and two stale SymPy docstring/ledger labels, and (F2) a single tautological assertion on the Mathematica side (`wl:140`, `expr - expr`) where the headline `Xi_micro` is posited rather than derived; the deliverable is still covered by the adjacent non-tautological chi/D-M consistency checks and by the SymPy derivation, so this is a cleanup, not a coverage gap. No `paper_misalignment`, no stop-cold.

## Self-test notes

- Variable-independence: F2's prescribed `xiMicroDerived = (lambdaPhi/theta)*(lambdaPhi*ell^2/tX)` depends on `lambdaPhi, theta, ell, tX` and reduces to `lambdaPhi^2*ell^2/(theta*tX)`, matching `xiMicro` — the `expectZero` residual is genuinely 0 only because the two routes coincide, not by construction (the LHS is built from the Pe×drop product, the RHS from line 132), so it is a real check. No spurious `diff(expr, var)` with var absent.
- Trivial-case: substituting any positive numeric profile (e.g. `lambdaPhi=theta=ell=tX=1`) gives `xiMicroDerived = 1 = xiMicro`, residual 0 — consistent.
- Parity/integral: the only unbounded-domain object is the BVP `dsolve`; integrals are over finite `[0,L]`/`[0,1]` and the normalization checks already return 0 in both transcripts — not affected by the fix.
- Paper round-trip: the F2 fix keeps the headline at `Lambda^2 L^2/(Theta T_X)`, identical to notes:203 — no new `paper_misalignment` introduced.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `mu_sigma^(chem) = Theta log(sigma/sigma_*) - Lambda_phi phi` | py:42 / wl:48 (sympy txt:6, math txt:6) | notes:31,114; eq form in .tex via EL closure | MATCH |
| support EL `-T_X phi_ss + K_X phi = Lambda_phi sigma` | py:153 / wl:165 (sympy txt:32, math txt:41) | notes:118 | MATCH |
| Onsager current `J = -D_sigma sigma_s + M_sigma Lambda sigma phi_s` | py:48 / wl:56 (sympy txt:8, math txt:9) | notes:145 | MATCH |
| Einstein relation `D_sigma = M_sigma Theta` | py:51 / wl:60 (sympy txt:11, math txt:13) | notes:149 | MATCH |
| `Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1)` | py:74 / wl:93 (sympy txt:16) | notes:56,176; eq:app-stage060 family | MATCH |
| `Pe = (Lambda_phi/Theta) Delta phi` / `gamma = Lambda Delta_phi/(Theta L)` | py:87 / wl:104 (sympy txt:19, math txt:22) | notes:59,180 | MATCH |
| `Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X)` | py:125 / wl:132,140 (sympy txt:25, math txt:32) | notes:203 | MATCH |
| `Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X` (chi form) | py:126 / wl:141 (sympy txt:27) | notes:211,266,89 | MATCH |
| `Xi_micro = M_sigma Lambda^2 L^2/(D_sigma T_X)` (D/M form) | py:127 / wl:142 (sympy txt:28) | notes:215,77 | MATCH |
| passivity `dF/dt = -int J^2/(M sigma) <= 0` | py:144-145 (density>=0) / wl:151-156 | notes:248 | MATCH (qualitative; density positivity asserted) |

INTERNAL (scaffolding, no prose expected): `f_sigma` free-energy density (py:38); `J on affine-drop branch` intermediate (txt:12); normalization constant `C = Delta_phi Lambda_phi/(L Theta (exp(...)-1))` (txt:14); BVP `phi_BVP(0)`, `phi_BVP(L)`, and `Delta_derived = L^2 Lambda sigma_0/(2 T_X)` rigid-grounding scale anchor (txt:21-23 — sanity check for the `Lambda L^2/T_X` scale, not a stated deliverable; the notes' Stage-058 normalization route gives the headline drop with coefficient 1, the BVP route gives 1/2 under different BCs, and the script comment py:120-122 flags the distinction); `dissipation density = M sigma mu_s^2` (txt:29); pass/fail flags and residual-zero check values.

reconciliation: complete; 10 values checked, 0 misaligned
