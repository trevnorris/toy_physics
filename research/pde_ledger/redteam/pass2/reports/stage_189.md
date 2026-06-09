---
unit_id: 189
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
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md]
  paper_appendix: present
---

# Audit unit 189 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_189.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 109, 1466 for MTDC-T9.4)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.txt`

## What the paper claims

Stage 189 is the "transfer-shape / outgoing-prefactor compiler" (anchor MTDC-T9.4). The card `\stagefield{Output}` says it "Equates \(\Xi_1\) with the logarithmic transfer-shape and outgoing-prefactor slope, including corrected nontracking form." The notes enumerate five deliverables: (1) the exact first-order transfer-shape identities `δln T² = Ξ₁`, `δln(1−ε_η) = R₁+Ξ₁`, `δln R_target = R₁`, with the selected-branch compatibility `δln R_target + δln T² − δln(1−ε_η) = 0` and `R_target·T² = Λ₀(1−ε_η)`, `Λ₀ = 27π²G c_s⁵/(20 a⁵ c⁵)`; (2) the isotropic grouped compiler `(D₀,D₂,D₄,N₀,N₂,N₄) ↦ (u₂,u₄,P₀,P₂,P₄)` with the closed forms `u₂=−D₂/D₀`, `u₄=(D₂²−D₀D₄)/D₀²`, `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³`; (3) the outgoing-branch coefficients `K₀=P₀`, `K₂=P₂+A P₀`, `K₄=P₄+A P₂+B P₀`, `Γ₅=G₅ P₀=(a⁵/27c_s⁵)P₀` with `A=a²/9c_s²`, `B=4a⁴/81c_s⁴`, `G₅=a⁵/27c_s⁵`; (4) the normalization equivalence `m₀²Γ₅ = 2G/5c⁵ ⟺ m₀²P₀ = 54G c_s⁵/(5 a⁵ c⁵)`; (5) the carry-forward bridge `Ξ₁ = δln T_A²/(ελ_A) = P₁/P₀` with the direct-slope law `δln T_A² = ε λ_A Ξ₁` and `P₁=(N₁D₀−N₀D₁)/D₀²`. The card's Checks also demand: keep `Ξ₁` separate from selected-branch dressing residuals, do not collapse Packet A and Packet B prematurely, and use the corrected nontracking form. The card's `\stagefield{Verification}` line still says "Mathematica audit: none yet" — but that is a paper-prose lag, not a script defect (a `.wl` does exist; flagged below as a non-blocking note, paper-side, deferred).

## What the script claims to verify

The SymPy script verifies all five notes-deliverables as `expect_zero` residual checks: the observable→transfer-packet compiler matrix and the three defect identities (Sec I), the one-port continuum transfer shape + coherent specialization + the direct-slope bridge `δln T_A² − ελ_A Ξ₁ = 0` (Sec II), the `u`/`P` series compiler and the `P₀ = (K_bl/D₀)T_eff²` bridge (Sec III), the weak-axisymmetric prefactor slope `P₁` and `P₁/P₀ = N₁/N₀ − D₁/D₀` (Sec IV), the `K`/`Γ₅` outgoing fingerprint compiler (Sec V), and the normalization product `m₀²Γ₅ = 2G/5c⁵` at `P₀ = 54G c_s⁵/(5a⁵c⁵ m₀²)` plus the constant-prefactor branch (Sec VI). The Mathematica script verifies the same five deliverables via independent routes (finite-log Jacobian, exponential-path slopes, `D[]` autodiff of the rational series, `Series`, `Coefficient`, `Solve`).

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `δln T² = Ξ₁`, `δln R_target = R₁`, `δln(1−ε_η) = R₁+Ξ₁` | py 76,78,80 / wl 94,95 | match |
| compatibility `δln R_target + δln T² − δln(1−ε_η) = 0` | py 84-87 / wl 96-99 | match |
| `C_obs→trf` compiler matrix, rank 2 | py 62-63 / wl 82-83 | match |
| `R_target·T² = Λ₀(1−ε_η)`, `Λ₀ = 27π²G c_s⁵/(20a⁵c⁵)` | py 98,101 (printed) / wl 103,105 (printed) | match (definition, not assertion — see note) |
| `δln T_A² = ε λ_A Ξ₁` | py 136 / wl 127-130 | match |
| `u₂,u₄` compiler | py 149,152 / wl 141-147 | match |
| `P₀,P₂,P₄` compiler | py 155-160 / wl 152-160 | match |
| `P₀ = (K_bl/D₀)T_eff²` | py 164 / wl 163 | match |
| `P₁`, `P₁/P₀ = N₁/N₀ − D₁/D₀` | py 172-176 / wl 168-174 | match |
| `K₀,K₂,K₄,Γ₅` + `A,B,G₅` | py 192-198 / wl 184-197 | match |
| `m₀²Γ₅ = 2G/5c⁵`, `Γ₅=a⁵P₀/27c_s⁵` | py 207-214 / wl 202-206 | match |
| `m₀²P₀ = 54G c_s⁵/(5a⁵c⁵)` (as `P₀_target`) | py 206,209 / wl 201,204 | match |
| constant-prefactor branch `P₂=0,P₄=0` ⟹ `K₂=A P₀`, `K₄=B P₀` | py 216-221 / wl 208-216 | match |

All notes deliverables have a faithful, non-tautological script-side check in both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `expect_zero(C·obs − trf)` | deliverable 1 (compiler matrix) | yes |
| A2 | sympy | 76 | `dlnT2.subs − Xi1 == 0` | `δln T² = Ξ₁` | yes |
| A3 | sympy | 80 | `dlnOneMinus.subs − (R1+Xi1) == 0` | `δln(1−ε_η) = R₁+Ξ₁` | yes |
| A4 | sympy | 84-87 | compatibility residual | rank-2 compatibility | yes |
| A5 | sympy | 136 | `direct_slope − ελ_A Xi1_closed == 0` | `δln T_A² = ελ_A Ξ₁` | yes |
| A6 | sympy | 152 | `Y_series − (1+u₂w²+u₄w⁴) == 0` | `u₂,u₄` | yes |
| A7 | sympy | 160 | `Pref_series − (P₀+P₂w²+P₄w⁴) == 0` | `P₀,P₂,P₄` | yes |
| A8 | sympy | 164 | `P0 − (K_bl/D0)T_eff² == 0` | static prefactor bridge | yes |
| A9 | sympy | 173-176 | `P₁`, `P₁/P₀` checks | weak-axisym slope | yes |
| A10 | sympy | 196-198 | outgoing-branch expansion | `K₀,K₂,K₄,Γ₅` | yes |
| A11 | sympy | 207-214 | normalization product | `m₀²Γ₅`, `Γ₅` | yes |
| A12 | sympy | 218-221 | constant-prefactor branch | `K₂=A P₀`, `K₄=B P₀` | yes |
| B1 | mathematica | 82 | `expectZero(cObsToTrfNative − expected)` | compiler matrix (Jacobian route) | yes |
| B2 | mathematica | 83 | `MatrixRank − 2 == 0` | rank-2 compatibility | yes |
| B3 | mathematica | 94 | `xiFromObservable − xi1 == 0` | `δln T² = Ξ₁` (path-slope route) | yes |
| B4 | mathematica | 95 | `oneMinusSlope − (rOneSlope+xi) == 0` | `δln(1−ε_η) = R₁+Ξ₁` | yes |
| B5 | mathematica | 96-99 | compatibility | rank-2 compatibility | yes |
| B6 | mathematica | 127-130 | `directSlope − ε... xi1Closed == 0` | `δln T_A² = ελ_A Ξ₁` | yes |
| B7 | mathematica | 144-147 | `u₂,u₄` from `D[]` autodiff | `u₂,u₄` | yes |
| B8 | mathematica | 156-160 | `P₀,P₂,P₄` from `D[]` autodiff | `P₀,P₂,P₄` | yes |
| B9 | mathematica | 163 | `P0 − (K_bl/D0)T_eff² == 0` | static bridge | yes |
| B10 | mathematica | 172-174 | `P₁` / `P₁/P₀` from `D[]` | weak-axisym slope | yes |
| B11 | mathematica | 194-197 | outgoing expansion via `Series`/`Coefficient` | `K₀,K₂,K₄,Γ₅` | yes |
| B12 | mathematica | 202-206 | normalization product | `m₀²Γ₅`, `Γ₅` | yes |
| B13 | mathematica | 208-216 | constant-prefactor branch via `Solve` | `K₂=A P₀`, `K₄=B P₀` | yes |

No tautological rows: every check substitutes an independently-supplied form into a residual that genuinely could fail (a wrong closed-form would leave a nonzero residual).

## Findings

None. No script-side or paper-side (value) findings.

## Independent-derivation check (Mathematica)

This is the heart of the V.3 retrofit audit. For each deliverable I state the method on both sides:

**Deliverable 1 — observable→transfer compiler matrix `C_obs→trf`.**
- `.py:45-62`: posits the linear symbolic expressions `dlnT2 = dn − Bstar*dr`, `dlnOneMinus = −epsetas/(1−epsetas)*de`, `dlnRtarget = dlnOneMinus − dlnT2`, hand-writes the 3×3 matrix `C_obs_to_trf`, and checks `C·obs − trf == 0`.
- `.wl:64-83`: builds *nonlinear finite-log functions* `transferLogFunctions = {xN − bStar·xR, Log[(1 − etaStar·Exp[xEta])/(1 − etaStar)], Log[...] − (xN − bStar·xR)}`, then recovers the matrix as the **Jacobian** `D[transferLogFunctions[[i]], obsVars[[j]]] /. zeroObs` (autodiff at the origin) and compares to the expected matrix.
- Verdict: **GENUINELY INDEPENDENT.** The `.py` posits the linearized matrix directly; the `.wl` derives it by linearizing nonlinear finite-log functions around `obs = 0`. Different derivation route (autodiff of a nonlinear map vs. hand-posited linear matrix).

**Deliverable 1b — the three defect identities `δln T² = Ξ₁` etc.**
- `.py:73-87`: substitutes `obs_sub = {dr:Theta1, dn:Xi1+Bstar*Theta1, de:Sigma_eta}` into the *already-linearized* symbolic `dlnT2`/`dlnRtarget`/`dlnOneMinus`.
- `.wl:85-99`: builds *finite exponential paths* `rTrPath = r0Obs·Exp[s·theta1]`, `nStarPath = n0Obs·Exp[s·(xi1+bStar·theta1)]`, `etaPath = etaStar·Exp[s·sigmaEta]`, `tSqPath = nStarPath/rTrPath^bStar`, then takes `D[Log[tSqPath], s] /. s→0` — a finite-perturbation log-slope. `oneMinusSlope` is taken from the *nonlinear* `Log[1 − etaPath]`, not from the pre-linearized `−ε/(1−ε)·δ` form the `.py` uses.
- Verdict: **GENUINELY INDEPENDENT.** Linearized-symbolic substitution (py) vs. differentiation of finite exponential paths (wl) is a different route to the same slopes.

**Deliverable 5 / the prompt's iteration-2 concern — direct-slope bridge `δln T_A² = ελ_A Ξ₁`.**
- `Ξ₁` is supplied INDEPENDENTLY on both sides from the raw input perturbation amplitudes: `.py:119-125` `eps1 = diff(epssplit,epsW)*epsW1 + diff(epssplit,deltaU)*deltaU1; Xi1_closed = zetaZ − omegaW + 2*chi1/(1+chi0) + 2*eps1/(1−epssplit)`; `.wl:116-120` identical. It is NOT back-derived from `direct_slope`.
- The slope itself: `.py:130-135` perturbs each variable of `T2_coh` by `s·epsA·λ_A·<amplitude>` and takes `diff(log(T2_coh_pert/T2_coh), s).subs(s,0)`; `.wl:121-126` does the identical perturbation construction and `D[Log[...], s]/.s→0`. The check `direct_slope − ε λ_A Xi1_closed == 0` then genuinely tests that the perturbed-transfer-shape log-slope reproduces the hand-built `Xi1_closed` (a wrong amplitude combination would leave a nonzero residual).
- The residual back-definition tautology flagged in iteration-2 (`R_target·T² − Λ₀(1−ε_η) ≡ 0`) is confirmed DEMOTED in both engines: `.py:101,105` only `sp.pprint`s `Rtarget_definition`; `.wl:105,107` only `Print`s `rTargetDefinition`. Neither is fed to an `expect_zero`/`expectZero`. The demotion held.
- Verdict on this section: **GENUINELY INDEPENDENT for the Jacobian/defect/slope derivation routes** (Sec I), and the demoted definition is correctly non-asserted in both engines. NOTE: the `Xi1_closed` perturbation construction (Sec II, py 119-136 vs wl 116-130) is the SAME variable choreography on both sides (same amplitude list, same `T2_coh_pert` build, same `d/ds log` at 0) — i.e. Sec II is a port of the perturbation algebra. But Sec I's recovery of the *same* `δln T² = Ξ₁` slope is fully independent (finite-log Jacobian + exponential paths vs. linearized substitution), so the deliverable `δln T_A² = ελ_A Ξ₁` is independently cross-validated by a genuinely different route in Sec I. Net: BORDERLINE-on-Sec-II-alone but INDEPENDENT at the deliverable level — not a `mathematica_transliteration` finding.

**Deliverables 2,3,4 — `u`/`P` compiler, `K`/`Γ₅`, normalization, constant-prefactor branch.**
- `.py` Sec III/V: uses `sp.series(D0/Dcons, ω, 0, 6).removeO()` and `sp.series(D0*Nseries/Dcons**2, ...)` then compares to hand-written closed forms `u2,u4,P0,P2,P4`; outgoing via `sp.Poly(Pref_trunc*Yhat, ω).as_dict()`.
- `.wl` Sec III/V: extracts the same coefficients via `D[response[w], {w,2}]/2 /. w→0` (Taylor coefficient by autodiff) and `D[prefactor[w], {w,4}]/24`, and the outgoing via `Series[...,{w,0,5}]` + `Coefficient[outTrunc, w, k]`.
- Constant-prefactor branch: `.py:216-217` hand-writes `N2_const = 2 D2 N0/D0`, `N4_const = (D2²+2D0D4)N0/D0²`; `.wl:208-209` instead **`Solve`s** `p2Expected == 0` for `n2` and `(p4Expected/.n2→n2Const)==0` for `n4`. Genuinely different (algebraic Solve vs. hand-substituted form).
- Verdict: **GENUINELY INDEPENDENT.** `sp.series`/`Poly.as_dict` (py) vs. `D[]`-autodiff Taylor coefficients + `Series`/`Coefficient` + `Solve` (wl) are different mechanizations of the series-coefficient extraction. The closed forms are the targets, not the routes; the residual checks would fail if a closed form were wrong.

Overall: the `.wl` is a genuine second engine. The only section with shared variable choreography is the Sec II `Xi1_closed` perturbation build, but that deliverable is independently confirmed by Sec I's distinct route, so no transliteration finding is warranted.

## Engine cross-check

Both engines run the same 12-13 residual checks and all report exactly `= 0` / `PASS` (sympy output lines 42-44, 73, 84, 94-95, 100-102, 107, 112-117; mathematica output lines 11-19, 27-30, 35-39, 45-49, 55, 61-72). The printed intermediate forms agree: `Λ₀ = 27π²G c_s⁵/(20a⁵c⁵)` (sympy 61-64, math 26), `T_A² = Z_W(1+ρ)²/(Ω_W²(1−ε_W)²)` (sympy 49-54, math 24). No residual/sign/factor disagreement.

## Verdict justification

`clean`. Attacks tried and failed: (1) I checked each `expect_zero`/`expectZero` for tautology — none defines `x=expr` then asserts `x==expr`; each substitutes an independently-built closed form into a residual that could be nonzero. (2) I checked the iteration-2 concern: the back-definition `R_target·T² − Λ₀(1−ε_η)` is confirmed demoted to a printed-only definition in both engines, and `Ξ₁`/`Xi1_closed` is built from raw input amplitudes, not back-derived from the slope. (3) I checked the V.3 retrofit independence deliverable-by-deliverable: the `.wl` uses finite-log Jacobian, exponential-path slopes, `D[]`-autodiff Taylor coefficients, `Series`/`Coefficient`, and `Solve` — genuinely different mechanizations from the `.py`'s `sp.diff`/`sp.series`/`Poly` choreography. (4) Numbering: banner reads "STAGE 189" in both scripts and both outputs; no stale 172/206/240/239 labels in the card or scripts (the card carries no `+17=206` self-label and no legacy "Stage 240"/239 prose — that legacy is in the notes prose only, which is out of red-team scope). (5) All notes deliverables map to a script check on both sides; all emitted constants reconcile with card+notes (see reconciliation section). One non-blocking observation: the card's `\stagefield{Verification}` line says "Mathematica audit: none yet", but a `.wl` does exist and passes — this is paper-prose lag, paper-side, deferred (not a script defect). Outputs are fresh (both .txt mtimes 11:53 > script mtimes 11:48). Verdict stands as clean.

## Self-test notes

I checked variable-independence on the `D[]`/`sp.diff` derivatives: every differentiated expression genuinely depends on the differentiation variable (`tSqPath` depends on `s` via the path exponentials; `response[w]`/`prefactor[w]` depend on `w`; `pLane[eps]` depends on `eps`; `epssplit_pert`/`t2CoherentPert` depend on `s`), so no derivative is identically zero — the `expect_zero` checks are non-trivially satisfied. No unbounded-domain integrals appear (no parity trap). Trivial-case spot check: setting `P₂=P₄=0` (constant-prefactor branch) reduces `K₂ → A P₀`, `K₄ → B P₀` exactly as both engines assert. No directive is written (zero findings).

## Value Reconciliation (pass-2 augmentation)

Reconciliation of every deliverable value the scripts emit against the `.tex` card and the `.md` notes. Outputs are fresh; reconciliation is based on script source + saved outputs.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δln T² = Ξ₁` | py:76 / wl:94 (sympy out 42, math out 14) | notes:128 (boxed), card Output line 15 | MATCH |
| `δln R_target = R₁`, `δln(1−ε_η) = R₁+Ξ₁` | py:78,80 / wl:95 (sympy out 43) | notes:166-173 (boxed) | MATCH |
| compatibility `δln R_target + δln T² − δln(1−ε_η) = 0` | py:84-87 / wl:96-99 (sympy out 44) | notes:216-224 (boxed) | MATCH |
| `Λ₀ = 27π²G c_s⁵/(20a⁵c⁵)` | py:98 / wl:103 (sympy out 61-64, math out 26) | notes:144 (boxed) | MATCH |
| `T_A² = Z_W(1+ρ)²/(Ω_W²(1−ε_W)²)` | py:100 / wl:104 (sympy out 49-54, math out 24) | notes:238-243 (boxed) | MATCH |
| `δln T_A² = ε λ_A Ξ₁` | py:136 / wl:127 (sympy out 73, math out 29) | notes:291-295 (boxed) | MATCH |
| `u₂=−D₂/D₀`, `u₄=(D₂²−D₀D₄)/D₀²` | py:149,152 / wl:141-142 (sympy out 84) | notes:330-339 (boxed) | MATCH |
| `P₀=N₀/D₀` | py:155 / wl:152 (sympy out 94) | notes:355-357 (boxed) | MATCH |
| `P₂=(D₀N₂−2D₂N₀)/D₀²` | py:156 / wl:153 | notes:390-392 (boxed) | MATCH |
| `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` | py:157 / wl:154 | notes:395-398 (boxed) | MATCH |
| `P₀=(K_bl/D₀)T_eff²` | py:164 / wl:163 (sympy out 95) | notes:361-364 (boxed) | MATCH |
| `P₁=(N₁D₀−N₀D₁)/D₀²`, `P₁/P₀=N₁/N₀−D₁/D₀` | py:172,176 / wl:169,174 (sympy out 102) | notes:418-422 (boxed) | MATCH |
| `Ξ₁ = δln T_A²/(ελ_A) = P₁/P₀` | (cross-checked via A5+A9/B6+B10) | notes:432-439 (boxed), appendix:47,1466 | MATCH |
| `A=a²/9c_s²`, `B=4a⁴/81c_s⁴`, `G₅=a⁵/27c_s⁵` | py:183-185 / wl:178-180 | notes:466-472 (boxed) | MATCH |
| `K₀=P₀`, `K₂=P₂+A P₀`, `K₄=P₄+A P₂+B P₀`, `Γ₅=G₅P₀` | py:192-198 / wl:189-197 (sympy out 107) | notes:491-501 (boxed) | MATCH |
| `Γ₅ = a⁵P₀/(27c_s⁵)` | py:213 / wl:206 (sympy out 113) | notes:500,520 (boxed) | MATCH |
| `m₀²Γ₅ = 2G/(5c⁵)` | py:207-210 / wl:202-205 (sympy out 112) | notes:514 (boxed) | MATCH |
| `m₀²P₀ = 54G c_s⁵/(5a⁵c⁵)` (`P₀_target`) | py:206 / wl:201 (ledger sympy out 133) | notes:524-526 (boxed) | MATCH |
| constant-branch `K₂=A P₀`, `K₄=B P₀` | py:220-221 / wl:213-216 (sympy out 116-117) | notes:562-564 (boxed) | MATCH |
| `N₂=2D₂N₀/D₀` (P₂=0 cond), `N₄=(D₂²+2D₀D₄)N₀/D₀²` (P₄=0 cond) | py:216-217 / wl:208-209 | notes:547,555 (boxed) | MATCH |

INTERNAL scaffolding (accounted for, no finding): `C_obs_to_trf` matrix entries, `epssplit`/`eps1` perturbation amplitudes (`zetaZ, omegaW, chi1, epsW1, deltaU1`), `Xi1_closed` intermediate, `T2_coh`/`T2_coh_pert`, `Y_series`/`Pref_series` intermediate series, `rank=2`, `T_eff²=N₀/K_bl`, all pass/fail flags and residual `= 0` print lines.

reconciliation: complete; 21 deliverable values checked, 0 misaligned.
