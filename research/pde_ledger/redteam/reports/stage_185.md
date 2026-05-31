---
unit_id: 185
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage185_microscopic_monomials.md]
  paper_appendix: present
---

# Audit unit 185 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_185.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage185_microscopic_monomials.md` (internal heading "Stage 253"; this is the original-stage-order source for renumbered Stage 185)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (row 185 at line 101; full §"Microscopic monomial invariants" lines 918-1047)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.txt`

## What the paper claims

Stage 185 (`\stagefield{Output}`: "Pushes the branch composites to direct monomials \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\).") sharpens the Stage 184 branch composites into *direct microscopic monomials* of the coherent-kernel couplings. The load-bearing deliverables, drawn from the notes (authoritative derivation source) and the appendix §"Microscopic monomial invariants":
1. The three normal-form coordinates are exactly the first logarithmic drifts of three direct monomials: `δln C_tr,* = Σ_tr`, `δln C_nt,* = Σ_nt`, `δln ε_η = Σ_η`, with `C_tr,* = χ0^(1+δU*) δU^(1+χ0*)`, `C_nt,* = (Z_W/Ω_W²) ε_W^{E*} δU^{-F*}`, `ε_η = c_ηU²/(K_U K_η^eff)`, and `E_*`, `F_*` as in eq. (app-part05-EF-defs).
2. The explicit microscopic forms of `C_tr,*` and `C_nt,*` in the eight positive primitives `(λ_W, c_ηU, γ, K_U, K_η^eff, K_W^eff, μ_W, T_U)` (notes §2,§3; appendix eqs. Ctr-micro, Cnt-micro).
3. The triangular observable law restated in monomials: `Θ1 = -C_tr,* δlnC_tr,*`, `Ξ1 = A_tr,* δlnC_tr,* + δlnC_nt,*`, `R1+Ξ1 = -ε_η*/(1-ε_η*) δlnε_η`.
4. The first-order zero-defect compatibility ledger solved for the dependent triple `(τ1, κ_η, μ1)` (notes §6.1-6.3).
5. (Appendix only, lines 1026-1037) the 3×8 exponent matrix `M_*` whose dependent-(τ,κη,μ) minor has `det = 1+χ0* > 0`, hence rank 3 / 5-dim kernel. The notes do NOT carry this matrix; it is the lowest-authority source and is the setup consumed by Stage 186.

## What the script claims to verify

Both engines (a) reconstruct the five microscopic-ratio drifts `Σ_χ,Σ_δ,Σ_ε,Σ_Z,Σ_η` by building independent multiplicative perturbations `x → x·exp(εΛ·d_x)` of each kernel primitive, forming the kernel-definition ratios, differentiating at ε=0, and matching against the Stage-182 symbolic Σ definitions; (b) build the tracking/nontracking/dressing monomial ratios two ways (from the composite ratios and from a primitive-exponent compiler), confirm the two forms are equal, and confirm each monomial's first log-drift equals Σ_tr/Σ_nt/Σ_η; (c) restate the triangular observable law; (d) solve the zero-defect system for τ1, κ_η, μ1 and back-substitute to confirm the three Σ vanish. The script docstring matches deliverables 1-4 of the paper. It does not assemble `M_*` or its determinant minor.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Verdict |
|---|---|---|
| 1. `δlnC_tr,*=Σ_tr`, `δlnC_nt,*=Σ_nt`, `δlnε_η=Σ_η` | `Sigma_tr_direct - Sigma_tr` (sympy:167), `Sigma_nt_direct - Sigma_nt` (184), `Sigma_eta_direct - Sigma_eta` (194) | match |
| 1b. monomial defs (exponents) | `Ctr_ratio` (156), `Cnt_ratio` (171), `epseta0` (105) built from paper exponents | match |
| 2. explicit microscopic (primitive) forms | `Ctr_ratio_primitive` (157), `Cnt_ratio_primitive` (172), `epseta_ratio_primitive` (187) cross-checked to composite forms (165,182,191) | match |
| 3. triangular observable law (coefficients C_tr,*, A_tr,*) | `Theta_1/Xi_1 monomial law` (202,203), `R1+Xi1 complement` (206) | **partial** — coefficient factors not independently exercised (F1) |
| 4. zero-defect solve for (τ1,κη,μ1) | `sp.solve` (212,213) + hand form (214) + back-subst (222-233) | match |
| 5. (appendix) `M_*`, `det M_*^(τ,κη,μ) = 1+χ0*` | none explicit; rows verified via Σ defs, nonsingularity implied by unique solve | **partial/missing** — explicit minor determinant unverified (F2) |

Dominant pattern: the math the scripts verify matches the paper's stated claims (deliverables 1-4 fully; 5 substantively but not explicitly). No identity verified differs from a paper identity → `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 124-131 | `d ln x - d_x == 0` (8 primitives) | building block of 1/2 | yes |
| A2 | sympy | 132-145 | composite ratio = primitive product | kernel-def correctness for 1/2 | yes |
| A3 | sympy | 148-152 | `Σ_*_direct - Σ_* == 0` (5 ratios) | claim 1 (Σ defs) | yes |
| A4 | sympy | 165-167 | `Ctr` composite=primitive; `δlnCtr - Σ_tr` | claim 1, 1b, 2 | yes |
| A5 | sympy | 182-184 | `Cnt` composite=primitive; `δlnCnt - Σ_nt` | claim 1, 1b, 2 | yes |
| A6 | sympy | 189-194 | `ε_η` composite=primitive; `δlnε_η - Σ_η` | claim 1, 2 | yes |
| A7 | sympy | 202-203 | `Theta1 - (-C_tr* Σ_tr_direct)`, `Xi1 - (...)` | claim 3 | **partial** (coeff on zero) |
| A8 | sympy | 206 | `Rcombo_direct - Rcombo` | claim 3 | yes |
| A9 | sympy | 222-233 | back-subst of solve → Σ_tr/η/nt = 0 | claim 4 | yes |
| A10 | mathematica | 116-223 | identical check set; μ1 via `Solve[sigmaNt==0]` (207) | claims 1-4 | yes (A7 mirror partial) |

## Findings

### F1 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py:199-203`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl:191-198`

**What's wrong:**
The checks named "Theta_1 monomial law" and "Xi_1 monomial law" do not exercise the observable coefficients `C_tr,*` and `A_tr,*` that they appear to validate. In SymPy:
```
Theta1 = sp.simplify(-C_tr_star * Sigma_tr)            # uses symbolic Sigma_tr
expect_zero("Theta_1 monomial law", Theta1 - (-C_tr_star * Sigma_tr_direct))
```
This residual is `-C_tr_star * (Sigma_tr - Sigma_tr_direct)`. Because `Sigma_tr == Sigma_tr_direct` was already established at line 167, the residual is `C_tr_star · 0 ≡ 0` for *any* value of `C_tr_star`. The same holds for `Xi_1 monomial law` (`A_tr_star` multiplies `Sigma_tr - Sigma_tr_direct` plus `Sigma_nt - Sigma_nt_direct`, both already zero). The Mathematica mirror (lines 197-198) has the identical structure. So the paper-deliverable-3 coefficients `C_tr,* = χ0δU/[(1+χ0)(1+δU)(1+χ0+δU)]` and `A_tr,* = 2χ0/[(1+χ0)(1+δU)]` are correct only insofar as they were transcribed correctly from the paper; no assertion would fail if they were wrong.

**Why this matters:**
On a checkpoint stage the triangular observable law is one of the four headline deliverables. The current checks give false confidence that the coefficient algebra is verified when it is not exercised at all. A transcription slip in `C_tr_star`/`A_tr_star` would pass silently.

**Required change:**
Anchor the coefficients to an independent in-script source. The Stage 182/183 source for `Θ1` is `Θ1 = -[χ0δU/((1+χ0)(1+δU)(1+χ0+δU))]Σ_tr` and the Stage-251/183 source for `Ξ1` is `Ξ1 = A_tr,*Σ_tr + Σ_nt`. Replace the multiply-onto-zero checks with checks that the coefficients themselves equal their defining ratios, independently of the already-proven `Σ_*_direct == Σ_*` identity. Concretely (SymPy), add explicit coefficient identities that cannot collapse to `0·anything`:
```
expect_zero("C_tr,* coefficient", C_tr_star - chi0s*deltaUs/((1+chi0s)*(1+deltaUs)*(1+chi0s+deltaUs)))
expect_zero("A_tr,* coefficient", A_tr_star - 2*chi0s/((1+chi0s)*(1+deltaUs)))
```
That alone is still self-referential (definition vs definition). The substantive fix is to verify the *normalization* that fixes these coefficients: confirm `C_tr,*` is the reciprocal of the Stage-184 `C_*` exponent used in `𝔗_*`, i.e. that `C_tr,* · [(1+χ0)(1+δU)(1+χ0+δU)/(χ0δU)] - 1 == 0`, and that `A_tr,*` matches the appendix `A_tr` definition. If no independent in-script source for these factors exists, at minimum restructure the existing checks so the coefficient is NOT multiplied onto a quantity already known to be zero — e.g. compare `Theta1` against `-C_tr_star*Sigma_tr` recomputed through the *primitive* drift `Sigma_tr_compiled` (line 164) rather than `Sigma_tr_direct`, so the check fails on a wrong primitive-exponent compiler. See directive for the precise minimal edit.

**Verification:**
New `expect_zero`/`expectZero` line(s) for the coefficient anchoring appear in both scripts and exit 0; the existing "Theta_1 monomial law"/"Xi_1 monomial law" residuals are rewritten so they are not identically `coeff·0`.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py:211-233`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl:204-223`

**What's wrong:**
The appendix (the source the stage card points to) attributes to Stage 185 the exponent matrix `M_*` and the explicit claim (lines 1034-1038)
```
\det M_*^{(\tau,\kappa_\eta,\mu_1)} = 1+\chi_{0,*} > 0,
```
which is what gives `M_*` rank three and a five-dimensional kernel (the input Stage 186 consumes). The scripts verify the *rows* of `M_*` (each Σ drift law) and demonstrate nonsingularity *operationally* by solving the 3-equation system uniquely for `(τ1,κη,μ1)` and back-substituting, but they never assemble the dependent 3×3 minor nor assert its determinant equals `1+χ0*`. On the checkpoint "exact paper alignment" bar this is a gap: the specific determinant value (a load-bearing positivity statement) is unverified.

**Why this matters:**
The rank-3 / `det = 1+χ0* > 0` fact is the bridge from this stage to Stage 186's similarity-orbit tangent theorem. Verifying it directly in this stage closes the carry-forward rather than relying on Stage 186 to re-derive it. The determinant is cheap to check from the already-defined drift expressions.

**Required change:**
Add a determinant check that extracts the dependent-block coefficients from the already-defined `Sigma_tr`, `Sigma_nt`, `Sigma_eta` (do NOT hardcode `M_*` entries) and verifies the minor determinant equals `1+χ0*`. See directive for the precise edit and self-tested coefficient values.

**Verification:**
A new check `det M_*^(tau,keta,mu) - (1+chi0s) == 0` appears in both scripts and exits 0.

### F3 — paper_misalignment (NOT — cosmetic only; see note)

This is intentionally NOT filed as a finding requiring resolution; recorded here as a low-severity cosmetic observation. Both scripts open with `banner("STAGE 168 — DIRECT MICROSCOPIC MONOMIALS")` (sympy:41, wl:26) and the saved outputs echo "STAGE 168", while this is Stage 185 (the notes carry the original-stage label "Stage 253"). "168" is a stale copy-paste (Stage 168 is "Off-bundle slippage decomposition" per appendix line 67). It is a print label only — no assertion depends on it — so it is not a math or paper-alignment defect. Fixing the label is a trivial hygiene edit and is folded into the directive as a non-blocking item.

## Independent-derivation check (Mathematica)

The `.wl` is structurally close to the `.py`: identical Σ definitions, identical variable choreography, the same `firstRatioDrift` (D[ratio,ε]/Λ at ε=0) method, and the same ordered check list. This is a parallel implementation of the *same valid method* applied to the *paper's* premises (the Σ defs and monomial exponents are the paper's, not SymPy's invention), executed on a genuinely different CAS kernel (`FullSimplify`/`Solve` vs `simplify`/`cancel`/`solve`). There is one real divergence in derivation: SymPy hand-constructs `mu_sol` from the notes §6.3 closed form (sympy:214) then validates it against `Sigma_nt=0`, whereas Mathematica solves `Solve[sigmaNt==0, mu1]` (wl:207) directly — these arrive at the same μ1 by different routes. I judge this borderline-acceptable rather than `mathematica_transliteration`: the shared structure reflects a canonical drift-by-differentiation method and shared paper premises, not an echo of SymPy-specific intermediate algebra. Recorded as a note, not a finding.

## Engine cross-check

Both saved transcripts (exit 0) report every `expect_zero`/`expectZero` residual `= 0` and every `PASS`. The final printed forms agree algebraically: e.g. SymPy `tau1 = (-c1*deltaUs - c1 + chi0s*kU - deltaUs*gam1 + deltaUs*kU - gam1 + 2*kU)/(chi0s + 1)` equals Mathematica `tau1 = (-((1 + deltaUs)*(c1 + gam1)) + (2 + chi0s + deltaUs)*kU)/(1 + chi0s)` (verified by hand-expansion); `kappa_eta = 2*c1 - kU` in both; `R1+Xi1` matches (`c_etaU²(-2c1+kU+keta)/(K_U K_eta - c_etaU²)` vs `cetaURef²(2c1-keta-kU)/(cetaURef²-ketaRef kuRef)`, same up to sign cancellation in numerator/denominator). Engines agree.

## Verdict justification

`findings` (not clean, not stop-cold). I read the stage card, the notes (the authoritative derivation, internally "Stage 253"), and the appendix §"Microscopic monomial invariants". The math the scripts verify matches the paper's deliverables 1-4 exactly: I checked the kernel definitions (χ0, δU, ε_W, Z_W/Ω_W², ε_η — the dropped constant prefactors π²/L² and σ cancel in the log-drift ratios, so they are harmless), the Σ definitions, the monomial exponents and their primitive-compiler forms, the E_*/F_* coefficients, and the τ1/κη/μ1 compatibility solve, all against the notes — all aligned. Attacks that failed: tautology hunt on the ratio-drift checks (genuine — they differentiate kernel-built exponentials and compare to symbolic Σ; a wrong kernel definition would fail), the primitive-vs-composite cross-checks (substantive), the back-substitution solve (genuine), and engine agreement. Attacks that landed: the two "monomial law" checks multiply the observable coefficients onto an already-zero quantity (F1), and the appendix's `det M_*^(τ,κη,μ)=1+χ0*` claim is exercised only implicitly via the solve, not asserted (F2). Both are low-severity, script-side, mechanically fixable, and introduce no paper disagreement. `paper_alignment: aligned`; `engines_agree: true`; `outputs_fresh: true` (sympy script 2026-05-11 11:58 < output 12:48; wl script 2026-04-21 < output 2026-05-11 13:24).

## Self-test notes

(1) Variable-independence: confirmed `diff(Sigma_tr, tau1)=1+chi0s`, `diff(Sigma_nt, {tau1,keta,mu1})=(-F*,-1,1)`, `diff(Sigma_eta,keta)=-1`, all other dependent-block derivatives are 0 — so the F2 determinant minor `[[1+chi0s,0,0],[-F*,-1,1],[0,-1,0]]` has det `=1+chi0s` (cofactor along row 1), matching the appendix `P_(T,Kη,μ)` det (line 1165). No spurious-zero derivatives. (2) No unbounded integrals in this stage — parity trap N/A. (3) Trivial-case: with all drifts zero the new det check gives literal `1+chi0s ≠ 0` (nonzero as required), and the F1 coefficient-anchoring residuals reduce to literal 0. (4) Paths: edits target existing `.py` in `scripts/` and `.wl` in `mathematica/` — no new files. (5) Paper round-trip: F1 anchors to the paper's own `C_tr,*`/`A_tr,*` defs (appendix Ctr-Atr-defs) and F2 to `det = 1+χ0*` (appendix Mstar-minor); neither introduces a constant absent from the paper.
