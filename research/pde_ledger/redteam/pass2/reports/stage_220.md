---
unit_id: 220
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 220 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_220.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows: line 52 table entry; lines 366-513 the "Linear dynamic mixed-bundle kernel" narrative; theorem `thm:app-part07-dynamic-no-free-lunch` line 508)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.txt`

## What the paper claims

The card (`\stagefield{Output}`, line 15) states: *"Dynamic product-family theorem and phase-lag no-go: linear dynamic mixing preserves the static spatial families, and the first-order passive/outgoing insertion is a phase-lag term rather than a new real barrier-lowering kernel off resonance."* The derivation ledger (line 13) enumerates: build the dynamic 3×3 susceptibility; prove `det K_dyn = Delta_Pi D_Pi`; repeat the product-family factorization at finite frequency; differentiate the outgoing-port insertion; separate in-phase conservative reshaping from quadrature pumping. The notes (Section 10) and the Part-VII appendix (eqs `app-part07-dyn-coefficients` … `app-part07-phase-lag-no-real`) give the explicit deliverables: (1) the dynamic coefficients `K_B,A,W` and invariants `Delta_Pi,Q_Pi,D_Pi`; (2) `det K_dyn = Delta_Pi D_Pi`; (3) static reduction to Stage 219; (4) the six exact inverse entries `chi_*`; (5) `V_mix = -1/2 J^T K_dyn^{-1} J`; (6) the collinear factorization `chi_s = N_s/(Delta_Pi D_Pi)`; (7) the product-family theorem preserving exactly `x^-6, e^{-2κx}/x^4, e^{-4κx}/x^2`; (8) the outgoing-port derivative identity `∂_Π V_mix = -1/2 T_J^2`; (9) the linear correction `δV^{(1)} = -i/2 Γ T_J^2`; (10) the phase-lag no-go `Re δV^{(1)}=0`, `P_abs^{(1)} = (ωΓ/2) T_J^2 ≥ 0`.

## What the script claims to verify

Both engines verify, symbolically and over fully free symbols, the ten deliverables above: the determinant identity (M1/det), the six-residual static reduction to Stage 219 (M2), the six native-inverse-vs-closed-form entries (M3/`inverse_checks`), the quadratic susceptibility law (M4), the collinear factorization (M5), the primitive product-family form (M6/`deltaV_primitive`), the outgoing-port derivative identity (M7/`dVdPi`), the linear outgoing correction (M8), and the phase-lag no-go with the perfect-square absorbed-power form (M9). The `.wl` additionally proves the *closure* of the spatial-family support set via a `CoefficientRules`-based Laurent-support extraction (`laurentSupportSet === {{-6,0},{-4,1},{-2,2}}`), and carries an M7 self-test that `V_mix` genuinely depends on `Pi` (guarding the identically-zero-derivative trap). A small numeric off-pole slice confirms `Re=0`, `Im≠0`, `P_abs>0` are simultaneously realizable.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `K_B,A,W,Delta_Pi,Q_Pi,D_Pi` definitions (app eqs 372-388) | py 42-48 / wl 97-109; printed | match |
| `det K_dyn = Delta_Pi D_Pi` (app eq 403) | py 64-65 / wl 121 | match |
| static reduction to Stage 219 (notes §2.3) | py 69-74 / wl 130-138 | match |
| six exact inverse entries (notes §3.1) | py 85-100 / wl 144-161 | match |
| `V_mix = -1/2 J^T K^{-1} J` (app eq 409) | py 114-123 / wl 166-175 | match |
| collinear `chi_s = N_s/(Delta_Pi D_Pi)` (notes §3.2) | py 128-139 / wl 179-192 | match |
| product family `x^-6,e^{-2κx}/x^4,e^{-4κx}/x^2` only (app eq 416) | py 146-167 / wl 196-229 (+Laurent closure) | match |
| `∂_Π V_mix = -1/2 T_J^2` (app eq 436) | py 177-179 / wl 233-239 | match |
| `δV^{(1)} = -i/2 Γ T_J^2` (app eq 449) | py 188-191 / wl 243-246 | match |
| `Re δV^{(1)}=0`, `P_abs=(ωΓ/2)T_J^2≥0` (app eq 457) | py 193-198 / wl 250-264 | match |
| card `\stagefield{Verification}`: "Mathematica audit: none yet" | a passing `.wl` exists and verifies | **mismatch (paper stale)** |

`paper_alignment: partial` — every *symbolic deliverable* matches both engines; the single divergence is doc-side: the card and notes still assert SymPy-only coverage while a verified Mathematica engine is present.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65 | `simplify(det - Delta_Pi*D_Pi)==0` | det identity | yes |
| A2 | sympy | 69-74 | 6× `simplify(.subs(ω=0,Π=0) - static)==0` | static reduction | yes |
| A3 | sympy | 100 | `all(Kinv[i,j]-chi==0)` | inverse entries | yes |
| A4 | sympy | 123 | `simplify(deltaV_matrix - deltaV_formula)==0` | susceptibility law | yes |
| A5 | sympy | 139 | `simplify(deltaV_col + 1/2 chi_s S^2)==0` | collinear thm | yes |
| A6 | sympy | 167 | `simplify(deltaV_primitive - expected)==0` | product family | yes (form match only) |
| A7 | sympy | 179 | `simplify(dVdPi + 1/2 T_J^2)==0` | outgoing-port deriv | yes |
| A8 | sympy | 191 | `simplify(deltaV_linear_out - expected)==0` | linear correction | yes |
| A9 | sympy | 194 | `real_part == 0` | Re δV^{(1)}=0 | yes |
| A10 | sympy | 198 | `simplify(P_abs - ωΓ/2 T_J0^2)==0` | perfect-square P_abs | yes (load-bearing; F2 fix) |
| A11 | sympy | 238-242 | sample `Δ≠0,D≠0,|Re|<1e-12,|Im|>0,P_abs>0` | off-pole realizability | yes (sample) |
| M1 | math | 121 | `expectZero[detNative - DeltaPi DPi]` | det identity | yes |
| M2 | math | 138 | `expectZero[6 static residuals]` | static reduction | yes |
| M3 | math | 151-161 | 6× `expectZero[chi - Kinv[[i,j]]]` | inverse entries | yes |
| M4 | math | 175 | `expectZero[Vmix - VmixExpected]` | susceptibility law | yes |
| M5 | math | 189-192 | `expectZero[-1/2 Jcol.Kinv.Jcol + 1/2 chiS S^2]` | collinear thm | yes |
| M6a | math | 212 | `expectZero[Vprim - VprimExpected]` | product family (form) | yes |
| M6b | math | 216-221 | `expectZero[3 extracted family coeffs]` | product family (coeff extraction) | yes |
| M6c | math | 229 | `expectTrue[laurentSupportSet === {{-6,0},{-4,1},{-2,2}}]` | **family closure (no extra families)** | yes (independent route) |
| M7a | math | 235-238 | `expectTrue[(dVdPi/.cons)!=0]` | self-test: V_mix depends on Pi | yes (trap guard) |
| M7b | math | 239 | `expectZero[dVdPi + 1/2 TJ^2]` | outgoing-port deriv | yes |
| M8 | math | 246 | `expectZero[linearOutgoing - linearExpected]` | linear correction | yes |
| M9a | math | 255 | `expectTrue[phaseReal === 0]` | Re δV^{(1)}=0 | yes |
| M9b | math | 256 | `expectZero[Pabs - ω γ/2 TJ0^2]` | perfect-square P_abs | yes (load-bearing; F2 fix) |
| M9c | math | 257-264 | `expectTrue[Δ≠0 && D≠0 && Re=0 && Im≠0 && Pabs>0]` | off-pole realizability | yes (sample) |

No tautological rows. A6/M6a are form-matches but are backed by the independent extraction (M6b) and closure (M6c) checks, so the product-family deliverable is exercised beyond a mere construct-then-assert.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (doc stale relative to script: a verified second engine exists but the card/notes deny it)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_220.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md:560-592` (Section 10 + "Supporting file" list)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl` (the existing, passing engine)

**What's wrong:**
The card's `\stagefield{Verification}` line reads (stage_220.tex:11):
> `SymPy audit: \StageFile{...stage220...sympy_audit.py}.  Mathematica audit: none yet.`
But a Mathematica audit `.wl` exists, is committed, and passes all nine M-blocks (committed output, all `PASS`). The neighbouring stage 221 card already cites its `.wl` correctly:
> `Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage221_..._mathematica_audit.wl}.` (stage_221.tex:11)
The notes file is also stale: Section 10 opens "The accompanying **SymPy** audit verifies:" and the "Supporting file" list (line 592) names only the `.py`. Neither doc mentions the Mathematica engine.

This is a documentation-vs-script discrepancy. Per the v2 prompt and pass-2 rule, Codex must not silently edit paper/ or notes/; the direction is the user's call.

**Why this matters:**
The card understates the verification coverage of a checkpoint-bar stage. A reader citing this stage would believe it is single-engine when it is in fact dual-engine verified. It is the inverse of the usual `value_mismatch` (the docs claim *less* than the scripts deliver), but it is still a paper↔script misalignment that should be reconciled before the card is published.

**Required change:**
See `## Resolve before fix_loop` in the directive. No script change. The likely resolution is to update the card line 11 to cite the `.wl` via `\StageFile{...}` (mirroring stage 221) and refresh the notes Section 10 wording to "SymPy and Mathematica audits verify" with the `.wl` added to the supporting-file list — but that is a paper/notes edit the user must authorize.

**Verification:**
After user-authorized doc edit, stage_220.tex:11 should name the `.wl` path with `\StageFile{}` and drop "none yet"; notes Section 10 should list the `.wl`. No script re-run required (scripts already pass).

## Independent-derivation check (Mathematica)

**INDEPENDENCE CALL: independent (not a transliteration), with one acknowledged shared-operation signature.**

The two engines share the same high-level route on the *core* identities — both build the same 3×3 `K_dyn`, take the *native* `Det`/`Inverse`, and compare against the *same* hand-supplied closed forms (`chi_qq=1/D_Pi`, `chi_qW=(A G_W+R G_U)/(Delta_Pi D_Pi)`, etc.). This shared signature is unavoidable and *correct* second-engine behaviour: neither script "derives" the chi forms procedurally; each independently confirms `native-inverse == paper-published-form` (notes §3.1 / app inverse entries). Confirming the same paper claim with each engine's own algebra is exactly what the dual-engine policy wants — it is not a line-by-line port of derivation steps.

Crucially, the `.wl` carries genuinely *independent* verification machinery the `.py` lacks, which is what distinguishes it from a transliteration:
1. **M6 Laurent-support closure** (wl 223-229): `supportPolynomial = -2 VprimY DeltaPi DPi x^6`, then `CoefficientRules[..., {x,y}]` → `laurentSupportSet === {{-6,0},{-4,1},{-2,2}}`. This *independently proves the spatial-family set is closed* (no fourth family appears), a strictly stronger statement than the `.py`'s form-match `deltaV_primitive - expected == 0` (py 167). The `.py` never establishes closure; the `.wl` does, via a different extraction (`CoefficientRules` rather than subtract-and-simplify).
2. **M6b coefficient extraction** (wl 216-221): independently re-extracts `C6,C4,C2` via `Coefficient[VprimYCollected, x^-6]` etc. — a second, route-distinct confirmation.
3. **M7a self-test** (wl 235-238): `expectTrue[(dVdPi/.cons) != 0]` guards the identically-zero-derivative trap; the `.py` has no such guard.

Quoted correspondence to justify "independent, not port":
- Core identity (shared route, expected): py 64-65 `det_identity = K_dyn.det(); assert simplify(det_identity - Delta_Pi*D_Pi)==0` ↔ wl 120-121 `detNative = Det[Kdyn]; expectZero["M1...", detNative - DeltaPi DPi]`.
- Divergent route (only in `.wl`): wl 223-229 `supportPolynomial = ...; CoefficientRules[...]; laurentSupportSet === expectedSupportSet` — has **no `.py` counterpart**, so this block cannot be a transliteration of anything in the `.py`.

Conclusion: the `.wl` is an independent re-verification, not a `mathematica_transliteration` finding. The shared-operation signature is the native-inverse-vs-published-closed-form comparison, which is legitimate.

## Engine cross-check

Both engines agree at the level they claim. Every shared identity (det, static reduction, six inverse entries, susceptibility law, collinear theorem, product family, outgoing-port derivative, linear correction, phase-lag) reduces to `0` / `True` in both committed outputs:
- SymPy output: "All Stage 220 symbolic checks passed." with the off-pole slice `Delta_Pi=133.8125`, `D_Pi=10.0205…`, `V_mix=-0.40159…`, `deltaV_linear=-0.000675…*I` (pure imaginary, Re=0), `P_abs=0.000337…>0`.
- Mathematica output: all M1-M9 `PASS`, `M6 Laurent support = {{-6,0},{-4,1},{-2,2}}`, "All Stage 220 Mathematica checks passed."
The numeric off-pole slice uses *identical* sample rules in both engines (K=11,M=2,C=1,ϖ=5,Ω_U=3,Ω_W=4,R=2,G_U=1,G_W=2,ω=1/2,J=(1,2,1),Γ=1/10), and both report Re=0, Im≠0, P_abs>0 consistently. No `engine_disagreement`.

## Verdict justification

`verdict: findings`, one `paper_misalignment` (low). The math is sound and dual-engine verified: I tried to break the determinant identity, the inverse closed forms, the product-family closure, and the phase-lag no-go, and all held up against both the paper appendix equations and the notes. The first-pass F2 fix (the P_abs perfect-square assert) is genuinely load-bearing and non-tautological in BOTH engines (py 198 / wl 256): each independently confirms `ComplexExpand/expand_complex` of `Im[-i/2 Γ T_J0^2]` collapses to `-1/2 Γ T_J0^2` — which only holds because `T_J0` is real on the conservative branch — and therefore `P_abs = (ωΓ/2) T_J0^2 ≥ 0` is a real square, not a sample coincidence. The prior fix still holds. The `.wl` is an independent route, not a transliteration (it adds the Laurent-support closure proof and the M7 self-test that the `.py` lacks). The only defect is doc-side: the card (line 11) and notes (Section 10) still say "Mathematica audit: none yet" / "SymPy audit verifies" while a passing `.wl` is committed — a `paper_misalignment` routed to the user, not auto-fixed by Codex.

## Self-test notes

Checked: (1) variable-independence trap — the `.wl` explicitly self-tests `dVdPi /. cons != 0` (wl 237), confirming `V_mix` truly depends on `Pi` before the `∂_Π V_mix = -1/2 T_J^2` assert; not identically zero. (2) Parity/perfect-square — `P_abs` is `ωΓ` times `T_J0^2`, a genuine square; with `ω` free (real) the sample `ω=1/2>0` gives `P_abs>0`, and the symbolic assert proves the perfect-square form for all real symbols (≥0 contingent on ω,Γ≥0, which `Gamma` carries via `nonnegative=True`/`gammaOut>=0`; `omega` is only declared real, so the *symbolic* P_abs sign is square-times-ω, correctly ≥0 only for ω≥0 — the no-go statement is about magnitude, consistent with the paper's `ω Γ/2 T_J^2`). (3) Trivial-case — substituting the sample profile gives Re=0, Im=-0.000675, P_abs=+0.000337, matching both outputs. No new finding from the self-test; the only issue is the doc staleness in F1.

## Value Reconciliation (pass-2 augmentation)

All Stage 220 deliverables are **symbolic closed forms**; the only numeric values emitted are the off-pole sample-slice numbers, which are pure verification scaffolding (INTERNAL). Reconciliation is therefore symbolic-form matching against the appendix equations and notes.

| value (deliverable) | source (py / wl + output line) | .tex / .md location | status |
|---|---|---|---|
| `K_B = K - Mω² - C²/(ϖ²-ω²)` | py 42 / wl 97; sympy.txt:11 | app eq 372 `app-part07-dyn-coefficients` | MATCH |
| `A = Ω_U²-ω²`, `W = Ω_W²-ω²-Π` | py 43-44 / wl 98-99; sympy.txt:12-13 | app eq 372 | MATCH |
| `Delta_Pi = A·W - R²` | py 46 / wl 107; sympy.txt:14 | app eq 381 | MATCH |
| `Q_Pi = G_U²W + 2G_UG_WR + G_W²A` | py 47 / wl 108; sympy.txt:15 | app eq 383 | MATCH |
| `D_Pi = K_B - Q_Pi/Delta_Pi` | py 48 / wl 109; sympy.txt:16 | app eq 388 / notes eq 154 (boxed) | MATCH |
| `det K_dyn = Delta_Pi D_Pi` | py 65 / wl 121 (=0) | app eq 403 / notes eq 174 (boxed) | MATCH |
| `chi_qq = 1/D_Pi` | py 85 / wl 144; sympy.txt:22 | notes §3.1 eq 226 | MATCH |
| `chi_qU = (G_U W + R G_W)/(Δ_Π D_Π)` | py 82,86 / wl 145; sympy.txt:23 | notes §3.1 eq 228-231 | MATCH |
| `chi_qW = (A G_W + R G_U)/(Δ_Π D_Π)` | py 83,87 / wl 146; sympy.txt:24 | notes §3.1 eq 233-236 | MATCH |
| `chi_UU = (K_B W - G_W²)/(Δ_Π D_Π)` | py 88 / wl 147; sympy.txt:25 | notes §3.1 eq 238 | MATCH |
| `chi_UW = (K_B R + G_U G_W)/(Δ_Π D_Π)` | py 89 / wl 148; sympy.txt:26 | notes §3.1 eq 241 | MATCH |
| `chi_WW = (K_B A - G_U²)/(Δ_Π D_Π)` | py 90 / wl 149; sympy.txt:27 | notes §3.1 eq 244 | MATCH |
| `V_mix = -1/2 J^T K_dyn^{-1} J` | py 114-123 / wl 166-175 | app eq 409 / notes eq 213 | MATCH |
| `chi_s = N_s/(Δ_Π D_Π)` (collinear) | py 137 / wl 188; sympy.txt:33 | notes §3.2 eq 265 (boxed) + N_s eq 270 | MATCH |
| product family `{x^-6, e^{-2κx}/x^4, e^{-4κx}/x^2}` | py 160-167 / wl 206-229; math.txt:63 | app eq 416 / notes eq 304 (boxed) | MATCH |
| `C6,C4,C2` susceptibility coeffs | py 156-158 / wl 203-205; sympy.txt:36-38 | notes eqs 318-328 | MATCH |
| `∂_Π V_mix = -1/2 T_J²` | py 177-179 / wl 233-239 | app eq 436 / notes eq 371 (boxed) | MATCH |
| `T_J = e_W^T K_cons^{-1} J` | py 177 / wl 233; sympy.txt:42 | notes eq 391 | MATCH |
| `δV^{(1)} = -i/2 Γ T_J²` | py 190 / wl 245; sympy.txt:45 | app eq 449 / notes eq 414 (boxed) | MATCH |
| `Re δV^{(1)} = 0` | py 194 / wl 255; sympy.txt:46 | app eq 454 / notes eq 420 (boxed) | MATCH |
| `P_abs^{(1)} = (ωΓ/2) T_J² ≥ 0` | py 198 / wl 256; sympy.txt:48 | app eq 457 / notes eq 441 (boxed) | MATCH |
| static reduction → Stage 219 bundle | py 69-74 / wl 138; math.txt:25 | notes §2.3 (eqs 194-201) | MATCH |
| Mathematica engine present + passing | wl all M1-M9; math.txt:91 | card line 11 says "none yet"; notes §10 says "SymPy audit" only | **MISMATCH (doc stale)** → F1 |

INTERNAL items (verification scaffolding, no doc reflection expected): off-pole sample slice `Delta_Pi=133.8125`, `D_Pi=10.0205…`, `V_mix=-0.40159…`, `deltaV_linear=-0.000675…I`, `P_abs=0.000337…`; the `expectZero` residuals (all 0); `M6 Laurent support` set (internal closure check); the M7 self-test boolean; the M9 dissipative-sample boolean; tolerances `1e-12`.

**reconciliation: 22 symbolic deliverables checked, 21 MATCH, 1 MISMATCH** (the Mathematica-engine-present-vs-card-says-"none yet" doc staleness, folded into F1 as the sole `paper_misalignment`). Every *numeric* emission is INTERNAL scaffolding. No symbolic-value disagreement exists between scripts and docs; the lone misalignment is a verification-coverage label that understates the scripts.
