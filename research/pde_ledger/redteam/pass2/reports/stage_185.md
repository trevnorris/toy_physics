---
unit_id: 185
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00-06:00
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage185_microscopic_monomials.md]
  paper_appendix: present
---

# Audit unit 185 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_185.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage185_microscopic_monomials.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (row 185 at line 101; §"Microscopic monomial invariants" lines 918-1047; M_* drift map lines 1014-1038)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.txt`

## What the paper claims

Stage 185 (`\stagefield{Output}`: "Pushes the branch composites to direct monomials \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\).") sharpens the Stage 184 branch composites into *direct microscopic monomials* of the coherent-kernel couplings. Load-bearing deliverables, drawn from the notes (authoritative derivation source) and appendix §"Microscopic monomial invariants":
1. The three branch-adapted coordinates are exactly the first logarithmic drifts of three direct monomials: `δln C_tr,* = Σ_tr`, `δln C_nt,* = Σ_nt`, `δln ε_η = Σ_η`, with `C_tr,* = χ0^(1+δU*)·δU^(1+χ0*)`, `C_nt,* = (Z_W/Ω_W²)·ε_W^{E*}·δU^{-F*}`, `ε_η = c_ηU²/(K_U K_η^eff)`, and `E_*`, `F_*` as in eq. (app-part05-EF-defs / notes lines 51-62).
2. The explicit microscopic (primitive) forms of `C_tr,*` and `C_nt,*` in the eight positive primitives `(λ_W, c_ηU, γ, K_U, K_η^eff, K_W^eff, μ_W, T_U)` (notes §2,§3; appendix Ctr-micro/Cnt-micro).
3. The triangular observable law restated in monomials: `Θ1 = -C_tr,*·δlnC_tr,*`, `Ξ1 = A_tr,*·δlnC_tr,* + δlnC_nt,*`, `R1+Ξ1 = -ε_η*/(1-ε_η*)·δlnε_η`, with `C_tr,* = χ0δU/[(1+χ0)(1+δU)(1+χ0+δU)]`, `A_tr,* = 2χ0/[(1+χ0)(1+δU)]` (notes §5).
4. The first-order zero-defect compatibility ledger solved for the dependent triple `(τ1, κ_η, μ1)` (notes §6.1-6.3).
5. (Appendix only, lines 1026-1038) the 3×8 exponent matrix `M_*` whose dependent-(τ,κη,μ) minor has `det = 1+χ0* > 0`, hence rank 3 / 5-dim kernel. The notes do not carry this matrix; it is the lowest-authority source and is the setup consumed by Stage 186.

## What the script claims to verify

Both engines (a) reconstruct the five microscopic-ratio drifts `Σ_χ,Σ_δ,Σ_ε,Σ_Z,Σ_η` by building independent multiplicative perturbations `x → x·exp(εΛ·d_x)` of each kernel primitive, forming the kernel-definition ratios, differentiating at ε=0, and matching against the Stage-182 symbolic Σ definitions; (b) build the tracking/nontracking/dressing monomial ratios two ways (from the composite ratios and from a primitive-exponent compiler), confirm the two agree, and confirm each monomial's first log-drift equals Σ_tr/Σ_nt/Σ_η; (c) restate the triangular observable law (both an "independent slippage" form against symbolic Σ and a "monomial law" form against the primitive-compiler drift); (d) assemble the dependent 3×3 minor of M_* and confirm `det = 1+χ0s`; (e) solve the zero-defect system for τ1, κ_η, μ1 and back-substitute to confirm the three Σ vanish. The docstring matches deliverables 1-5.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Verdict |
|---|---|---|
| 1. `δlnC_tr,*=Σ_tr`, `δlnC_nt,*=Σ_nt`, `δlnε_η=Σ_η` | `Sigma_tr_direct - Sigma_tr` (sympy:167), `Sigma_nt_direct - Sigma_nt` (184), `Sigma_eta_direct - Sigma_eta` (194) | match |
| 1b. monomial exponent defs (`E_*`, `F_*`) | `E_star`/`F_star` (sympy:62-63) used inside `Cnt_ratio` (171) and compiler (172-179) | match |
| 2. explicit microscopic (primitive) forms | `Ctr_ratio_primitive` (157-162), `Cnt_ratio_primitive` (172-179), `epseta_ratio_primitive` (187) cross-checked to composite forms (165,182,191) | match |
| 3. triangular observable law (coeffs C_tr,*, A_tr,*; R-complement) | "Theta_1/Xi_1 monomial law" (sympy:212,213) now vs `Sigma_tr_compiled`/`Sigma_nt_compiled`; `R_1+Xi_1 complement` (216) | match (see Independent-derivation note on coefficient anchoring) |
| 4. zero-defect solve for (τ1,κη,μ1) | `tau_sol`/`keta_sol`/`mu_sol` (230-232) + back-subst (240-251) | match |
| 5. (appendix) `det M_*^(τ,κη,μ) = 1+χ0*` | `Mstar_minor` assembled from `diff(Sigma_*, {tau1,keta,mu1})` (222-228), `det - (1+chi0s)` (229) | match |

Dominant pattern: every paper-side deliverable, including the appendix-only M_* minor determinant, is now exercised by a script-side check. No identity verified differs from a paper identity → `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 124-131 | `d ln x - d_x == 0` (8 primitives) | building block of 1/2 | yes |
| A2 | sympy | 132-145 | composite ratio = primitive product | kernel-def correctness for 1/2 | yes |
| A3 | sympy | 148-152 | `Σ_*_direct - Σ_* == 0` (5 ratios) | claim 1 (Σ defs) | yes |
| A4 | sympy | 165-167 | `Ctr` composite=primitive; `δlnCtr - Σ_tr`; compiler form | claim 1, 1b, 2 | yes |
| A5 | sympy | 182-184 | `Cnt` composite=primitive; `δlnCnt - Σ_nt`; compiler form | claim 1, 1b, 2 | yes |
| A6 | sympy | 189-194 | `ε_η` composite=primitive; `δlnε_η - Σ_η` | claim 1, 2 | yes |
| A7 | sympy | 210-211 | `Theta1 - (-C_tr*·Σ_tr)`, `Xi1 - (A_tr*·Σ_tr+Σ_nt)` (symbolic) | claim 3 | Θ partial / Ξ yes |
| A8 | sympy | 212-213 | `Theta1 - (-C_tr*·Σ_tr_compiled)`, `Xi1 - (A_tr*·Σ_tr_compiled+Σ_nt_compiled)` | claim 3 (via primitive compiler) | yes |
| A9 | sympy | 216 | `Rcombo_direct - Rcombo` | claim 3 (R-complement) | yes |
| A10 | sympy | 229 | `det Mstar_minor - (1+chi0s)` | claim 5 | yes |
| A11 | sympy | 240-251 | back-subst of solve → Σ_tr/η/nt = 0 | claim 4 | yes |
| A12 | mathematica | 116-241 | identical check set; μ1 via `Solve[sigmaNt==0,mu1]` (225) | claims 1-5 | yes (A7 Θ mirror partial) |

## Findings

None. (See Independent-derivation check for the residual low-severity observation on the Θ1 coefficient anchoring; it is below the finding threshold given the minimal fix authorized by the first-pass directive was applied and the check now also exercises the primitive compiler.)

## Independent-derivation check (Mathematica)

The `.wl` (lines 48-241) is a **parallel implementation of the same canonical drift-by-differentiation method** used by the `.py`, not an echo of SymPy-specific intermediate algebra. Both:
- define the five Σ slippages identically (`sigmaChi=gam1+c1-kU`, … wl:48-52 = sympy:54-58);
- build the same multiplicative exponential perturbations `x → x·Exp[ε·Λ·d_x]` (wl:63-70 = sympy:71-78);
- apply the same `firstRatioDrift` operator `D[ratio,ε]/Λ |_{ε→0}` (wl:46 = sympy:36-38);
- compile the SAME primitive-exponent products for `Ctr_ratio_primitive` (wl:148-154 = sympy:157-162) and `Cnt_ratio_primitive` (wl:163-171 = sympy:172-179), which are exactly the paper's M_* matrix rows (appendix lines 1029-1031).

The single strongest shared-construct is the primitive compiler exponent vectors, e.g. nontracking `gamma^(2E*)·lamW^(2+2E*)·muW·tu^(-F*)·kU^(F*-E*)/(keta·kWeff^(2+E*))` (wl:163-171 vs sympy:172-179) — character-for-character identical exponent choices. By the prompt's literal test ("same black-box differentiation on the same constructed object = a port"), this is structurally a port.

However: (i) the constructed object's exponents ARE the paper's claim (M_* rows), so any engine that re-derives this stage must encode them — there is no independent route that avoids them; (ii) the two engines run on genuinely different CAS kernels (`FullSimplify`/`Solve[...,Reals]` vs `simplify`/`cancel`/`solve`); (iii) there is one real algorithmic divergence — SymPy hand-constructs `mu_sol` from the notes §6.3 closed form (sympy:232) then validates it, whereas Mathematica solves `Solve[sigmaNt==0,mu1]` directly (wl:225) and the two arrive at the same μ1 by different routes. Per the established MIRROR_POLICY norm for algebraic-identity checkpoint stages (same canonical method on independent kernels), I judge this the accepted mirror, NOT a flag-worthy `mathematica_transliteration`. This matches the first-pass judgment.

Residual note on the Θ1 coefficient (first-pass F1): the "Theta_1 independent slippage law" (sympy:210) residual is `C_tr_star·(Σ_tr−Σ_tr)≡0`; because `Theta1` (sympy:201-204) is constructed with the SAME `(1+chi0s)(1+deltaUs)(1+chi0s+deltaUs)` denominator literal as `C_tr_star` (sympy:197), the check cannot anchor the coefficient to an external source — though it DOES catch a disagreement between the two separate transcriptions of that denominator (sympy:197 vs sympy:203). The "Theta_1 monomial law" (sympy:212) now uses `Sigma_tr_compiled` (the primitive-exponent drift), so it additionally fails on a wrong primitive exponent. The Ξ1 checks are genuinely substantive: `A_tr_star·Σ_tr + Σ_nt` requires the `2χ0/(1+δU)·Σ_δ` terms to cancel, which fails if `A_tr_star` or `F_star` is wrong. This is the minimal fix the first-pass directive authorized; resolved at the accepted level, hence no finding.

## Engine cross-check

Both saved transcripts (sympy .txt, wl .txt) report every residual `= 0` and every `PASS`. Final printed forms agree algebraically:
- `Σ_tr`: sympy `-(chi0s+1)*(kU-tau1)+(deltaUs+1)*(c1+gam1-kU)` = wl `(1+deltaUs)*(c1+gam1-kU)+(1+chi0s)*(-kU+tau1)` (same after sign distribution).
- `det M_*^(τ,κη,μ) - (1+chi0s) = 0` in both (sympy .txt line 71; wl .txt line 103-104 PASS).
- `tau1`: sympy `(-c1*deltaUs - c1 + chi0s*kU - deltaUs*gam1 + deltaUs*kU - gam1 + 2*kU)/(chi0s+1)` = wl `(-((1+deltaUs)*(c1+gam1)) + (2+chi0s+deltaUs)*kU)/(1+chi0s)` (expanded by hand: identical numerators).
- `kappa_eta = 2*c1 - kU` in both.
- `R1+Xi1`: sympy `c_etaU²*(-2c1+kU+keta)/(K_U*K_eta_eff - c_etaU²)` = wl `cetaURef²*(2c1-keta-kU)/(cetaURef²-ketaRef*kuRef)` (multiply sympy by -1/-1: identical).
Engines agree.

## Verdict justification

`clean`. I read the stage card, the notes (the authoritative derivation), and the appendix §"Microscopic monomial invariants" including the M_* drift map. The math the scripts verify matches all five paper deliverables. I re-derived the load-bearing checkpoint quantity — the dependent-minor determinant — by hand from the Σ definitions: `∂(Σ_tr,Σ_nt,Σ_eta)/∂(τ1,κη,μ1) = [[1+χ0s,0,0],[-F*,-1,1],[0,-1,0]]`, det = `(1+χ0s)·((-1)(0)-(1)(-1)) = 1+χ0s`, exactly matching appendix eq:app-part05-Mstar-minor and the script's check (sympy:229) — it holds, and it is non-tautological because the minor entries are extracted by `diff` from the independently-built Σ expressions (a wrong Σ_nt exponent would change `-F*` and break the cofactor). I also re-derived `Σ_tr = χ0s·δUs·Σ_tr/(denom)` collapse inside `Theta1` and the `A_tr_star·Σ_tr + Σ_nt` cancellation, both confirmed. Attacks that failed: tautology hunt on the ratio-drift checks (genuine — they differentiate kernel-built exponentials vs symbolic Σ; a wrong kernel def fails), the primitive-vs-composite cross-checks (substantive — exponent typo fails), the det-minor check (now present and substantive), the back-substitution solve, and engine agreement. Attacks that landed at low severity but below the finding bar: the Θ1 "independent slippage law" coefficient is not externally anchored (first-pass F1, minimal fix applied — now also exercises the primitive compiler). First-pass F2 (det M_* minor unverified) is RESOLVED — the check is present at sympy:222-229 / wl:217-222. First-pass F3 (stale "STAGE 168" banner) is RESOLVED — banner now reads "STAGE 185" (sympy:41, wl:26) and the saved outputs echo "STAGE 185". The first-pass `insufficient_verification` and `stale` flags are therefore resolved, not recurring. `paper_alignment: aligned`; `engines_agree: true`; `outputs_fresh: true` (sympy .py 2026-05-30 09:03 < .txt 15:33; wl 2026-05-30 09:03 < .txt 15:34).

## Self-test notes

(1) Variable-independence: confirmed the dependent-minor `diff`s — `∂Σ_tr/∂τ1=1+χ0s`, `∂Σ_nt/∂(τ1,κη,μ1)=(-F*,-1,1)`, `∂Σ_eta/∂κη=-1`, all other dependent-block derivatives 0 — so the minor det `=1+χ0s` by cofactor along row 1; no spurious-zero derivatives. (2) No unbounded integrals in this stage — parity trap N/A. (3) Trivial-case: with all primitive drifts zero every `Σ_*_direct` and every monomial-drift residual reduces to literal 0, and the det check gives `1+χ0s ≠ 0` (nonzero literal, as required for the rank-3 claim). (4) No new scripts proposed — both engines present in correct dirs. (5) Paper round-trip: the dropped constant prefactors `σ` (in ε_W) and `π²/L²` (in δU) carry no drift symbol, so they cancel exactly in every multiplicative ratio and contribute 0 to first drifts — harmless, and they appear correctly in the notes/appendix definitions, so no `value_mismatch` is introduced.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value emitted by the scripts (and confirmed in the saved `.txt`) reconciles with the `.tex` card / appendix and the `.md` notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Σ_tr = (1+χ0s)(τ1−κU)+(1+δUs)(γ1+c1−κU)` | sympy:61 / .txt L5 | notes lines 122-149 (Σ defs + Σ_tr); appendix M_* row1 (1029) | MATCH |
| `Σ_nt = Σ_Z + E*·Σ_eps − F*·Σ_delta` | sympy:64 / .txt L6 | notes lines 201-214; appendix M_* row2 (1030) | MATCH |
| `Σ_eta = 2c1 − κU − κη` | sympy:58 / .txt L7 | notes line 133; appendix M_* row3 (1031) | MATCH |
| `E_* = 2ε_W*/(1−ε*)·(11+9δU*)/(11(1+δU*))` | sympy:62 / used in Cnt | notes lines 51-55 / appendix 948-949 | MATCH |
| `F_* = 2χ0*/(1+δU*) + 4ε_W*δU*/(11(1−ε*)(1+δU*)²)` | sympy:63 / used in Cnt | notes lines 57-61 / appendix 951-954 | MATCH |
| `C_tr,* = χ0^(1+δU*)·δU^(1+χ0*)` (exponents) | sympy:156 / .txt L81 carry-fwd | notes lines 158-163 / appendix 932-936 | MATCH |
| `C_tr,*` primitive exponents (γ:1+δU*, c:1+δU*, T:1+χ0*, KU:−(2+χ0*+δU*)) | sympy:157-162 | appendix M_* row1 (1029) | MATCH |
| `C_nt,* = (Z_W/Ω_W²)·ε_W^{E*}·δU^{−F*}` | sympy:170-171 / .txt L82 | notes lines 234-240 / appendix 939-943 | MATCH |
| `C_nt,*` primitive exponents (λW:2+2E*, μW:1, T:−F*, KU:F*−E*, Kη:−1, KW:−(2+E*), γ:2E*) | sympy:172-179 | appendix M_* row2 (1030) | MATCH |
| `ε_η = c_ηU²/(K_U K_η^eff)` | sympy:105 / .txt L83 | notes lines 299-303 / appendix 987-989 | MATCH |
| `Θ1 = −C_tr,*·Σ_tr`, coeff `C_tr,* = χ0δU/[(1+χ0)(1+δU)(1+χ0+δU)]` | sympy:197,201-204 / .txt L64 | notes lines 324-337 | MATCH |
| `Ξ1 = A_tr,*·Σ_tr + Σ_nt`, coeff `A_tr,* = 2χ0/[(1+χ0)(1+δU)]` | sympy:198,205-208 / .txt L65 | notes lines 326,341-343 | MATCH |
| `R1+Ξ1 = −ε_η*/(1−ε_η*)·Σ_eta` | sympy:209 / .txt L66 | notes lines 328-330 | MATCH |
| `det M_*^(τ,κη,μ) = 1+χ0s` | sympy:229 / .txt L71 | appendix eq:app-part05-Mstar-minor (1035-1037) | MATCH |
| `τ1 = κU − (1+δU*)/(1+χ0*)(γ1+c1−κU)` | sympy:230 / .txt L72 | notes §6.1 lines 412-420 | MATCH |
| `κη = 2c1 − κU` | sympy:231 / .txt L73 | notes §6.2 lines 425-427 | MATCH |
| `μ1` (and `μ1 on full zero-defect branch`) | sympy:232,237 / .txt L74-75 | notes §6.3 lines 432-453 (both boxed forms) | MATCH |

Internal scaffolding (accounted for, no finding): the eight primitive single-variable drift checks `d ln x − d_x` (sympy:124-131), the five composite-from-primitive ratio cross-checks (sympy:132-145), the `*_compiled`/`*_direct` intermediate drift quantities, the `Ctr_ratio`/`Cnt_ratio` ratio objects, `Mstar_minor` matrix, the back-substitution residuals, and all `PASS`/`= 0` flags. These exist only to drive assertions; they are not stage deliverables.

reconciliation: complete; 17 deliverable values checked, 0 misaligned.
