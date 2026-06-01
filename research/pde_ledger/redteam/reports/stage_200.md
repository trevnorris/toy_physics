---
unit_id: 200
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage200_reference_free_home_stretch_theorem.md]
  paper_appendix: present
---

# Audit unit 200 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_200.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage200_reference_free_home_stretch_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 131, 1379–1456, 1468)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.txt`

## What the paper claims

Stage 200 (anchor MTDC-T9.6) is the reference-free home-stretch theorem: it combines the Packet-A branch scalar with the Packet-B orbit-quotient triple into the single four-scalar final verdict packet. The card's `\stagefield{Output}` states verbatim: "Combines Packet A and Packet B into the four-scalar final verdict packet \((\Delta_Q,q_{\rm tr},q_{\rm nt},q_\eta)\)." The appendix (eq:app-part05-full-packet-additive/-multiplicative/-mismatch) defines the additive packet `\Delta_{\rm full}=(\Delta_Q,q_{tr},q_{nt},q_\eta)` with `\Delta_Q:=\chi_Q-1`, the multiplicative chart `(\chi_Q, R_{tr},R_{nt},R_\eta)`, and the mismatch chart `(\chi_Q,m_T,m_K,m_\mu)`, with the exact chart conversions `R_{tr}=e^{q_{tr}}=(m_T)^{1+\chi_{0,*}}`, `R_\eta=e^{q_\eta}=1/m_K`, `R_{nt}=e^{q_{nt}}=m_\mu/(m_K m_T^{F_*})`. The notes add the load-bearing deliverables: (1) the Packet-B compiler matrix derived from primitive monomial log-ratios reproduces the carried Stage 192 matrix `M_*`; (2) orbit-representative (witness) independence; (3) the mismatch chart re-derived from the dependent-triple orbit solve "rather than posited"; (4) the additive cocycle law; (5) the Packet-A linearization `\Delta_Q^{lin}`. Note: the notes body and both scripts' banners carry stale stage labels ("Stage 251" in notes, "STAGE 183" in script banners) inherited from a pre-renumbering revision; the mathematical content is the Stage 200 / MTDC-T9.6 packet.

## What the script claims to verify

The SymPy script verifies, in five sections: (I) the four-scalar Packet-B compiler matrix `M_*` obtained as the Jacobian of the log of the primitive monomial ratios (Ctr, Cnt, epsEta) equals the hardcoded carried Stage 192 matrix, and `q=M_*·Δx`; (II) witness invariance — applying the orbit scaling with the derived dependent-triple multipliers `Φ_T, Φ_K, Φ_μ` leaves all three monomials invariant, hence the pairwise-witness packet equals the intrinsic orbit packet; (III) the mismatch chart `q=((1+χ₀)ln m_T, ln m_μ−ln m_K−F_*ln m_T, −ln m_K)`; (IV) the additive cocycle `Δx^31=Δx^32+Δx^21` and `q^31=q^32+q^21`; (V) `Δ_Q^lin=ε(5ε_β+δΣ₀/(3S)+9δΣ₅/S)` from a series expansion of `χ_Q=3(Sβ⁵+9Σ₅)/(3S−Σ₀)`. The Mathematica script asserts the same five sections with the same checks. Both scripts exit 0.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Four-scalar packet `(\Delta_Q,q_{tr},q_{nt},q_\eta)`, `\Delta_Q=\chi_Q-1` | sympy V `Delta_H_linear`, II `Delta_H_intrinsic`; math equivalents | match |
| Packet-B compiler `M_*` from monomials = carried Stage 192 matrix | sympy I `Mderived-Mexpected==0`; math I | match |
| Orbit-representative independence (witness) | sympy II witness invariance + pairwise=intrinsic; math II | match |
| Mismatch chart conversions `R_{tr}=m_T^{1+\chi_0}`, `R_\eta=1/m_K`, `R_{nt}=m_\mu/(m_K m_T^{F_*})` | sympy III (substantive); math III (tautological — see F2) | match (sympy) / partial (math) |
| Additive cocycle law | sympy IV; math IV | match |
| Packet-A linearization `\Delta_Q^{lin}` | sympy V series; math V series | match |

`paper_alignment: aligned` — every paper deliverable maps to a script-side check, and the symbolic forms agree exactly. The stale internal stage labels (notes "251", banners "183") are cosmetic prose-side artifacts; they do not change any verified identity and are not a math `paper_misalignment`. Codex must not touch them (notes/banners on the prose side are out of scope; the banner strings are in-script but purely decorative, so I do not direct an edit).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 167 | `expect_zero(Mderived - Mexpected)` | M_* compiler (claim 1) | yes |
| A2 | sympy | 168 | `expect_zero(q_pair - Mderived*Dvec)` | M_* linearity | yes |
| A3 | sympy | 227–238 | `ln[Ctr/Cnt/epsEta(witness)/(*)]==0` | witness invariance (claim 2) | yes |
| A4 | sympy | 269 | `expect_zero(pair_witness - intrinsic)` | orbit-rep independence | yes |
| A5 | sympy | 313–315 | `Ctr/Cnt/epsEta(actual) ratio - m-form` | mismatch chart (claim 3) | yes (rebuilds monomials) |
| A6 | sympy | 339–340 | `q_mismatch - carried; M_*Δx_mis - q` | mismatch chart | yes |
| A7 | sympy | 396–397 | cocycle `D31-D32-D21`, `q31-q32-q21` | cocycle (claim 4) | partial (near-definitional) |
| A8 | sympy | 423 | `DeltaQ_linear - expected` (series) | Packet-A lin (claim 5) | yes |
| A9 | math | 167–168 | `Mderived - Mexpected`, `qPair - Mderived.Dvec` | M_* compiler | yes |
| A10 | math | 199–201 | witness invariance | witness invariance | yes |
| A11 | math | 218 | pairwise - intrinsic | orbit-rep independence | yes |
| A12 | math | 240–245 | mismatch `Log[Ctr/Cnt/epsEtaActualRatio]` | mismatch chart | NO — tautological (F2) |
| A13 | math | 253–254 | `qMismatch - expected`, `M_*.Dmis - qMismatch` | mismatch chart | partial |
| A14 | math | 274–275 | cocycle | cocycle | partial |
| A15 | math | 294 | `DeltaQLinear - expected` (series) | Packet-A lin | yes |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl:116-319`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py` rather than an independent native re-derivation. Concretely, the two scripts share the same Section I–V choreography, the same intermediate variable names, the same hardcoded target matrix, the same substitution dictionaries, and the same ledger prose. Three corresponding excerpts:

1. Hardcoded target matrix is identical, not independently reconstructed:
   - sympy:155–161 `Mexpected = sp.Matrix([[0, 1 + deltaUs, ...], ...])`
   - math:157–161 `Mexpected = {{0, 1 + deltaUs, ...}, ...}`
2. Same `ratio→exp(D)` substitution dictionary and same `Dvec` ordering:
   - sympy:133–142 `ratio_to_logs = {rla: sp.exp(Dl), ...}`; sympy:153 `Dvec = sp.Matrix([Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT])`
   - math:139–148 `ratioToLogs = {rla -> Exp[Dl], ...}`; math:155 `Dvec = {Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT}`
3. Section III mismatch uses the same `Dmis` literal vector and same expected chart:
   - sympy:326 `Dmis = sp.Matrix([0, 0, 0, 0, sp.log(mK), 0, sp.log(mMu), sp.log(mT)])`
   - math:248 `Dmis = {0, 0, 0, 0, Log[mK], 0, Log[mMu], Log[mT]}`

The Mathematica script does use `D[]` and `Series[]` (native primitives) in Sections I and V, so it is not a pure token-rewrite there; but the *derivation route* (which ratios, in which collapsed form, fed to which Jacobian, compared to which hardcoded literal) mirrors SymPy step-for-step. As a checkpoint (higher bar), the second engine must derive the result by an independent route to confer independent confidence; this one echoes SymPy's algebra. This finding overlaps with F2 (the most damaging consequence is the Section III tautology).

**Why this matters:**
A transliterated second engine provides no independent cross-check. If the SymPy route contains a hidden error in how a monomial ratio is set up, the mirrored Mathematica route reproduces the same error and still "agrees," giving false dual-engine confidence at a trust-anchor checkpoint.

**Required change:**
Rewrite the Mathematica Sections so the load-bearing results are obtained by Mathematica-native, independent routes rather than mirroring the SymPy pre-collapsed forms:
- Section I: build `CtrRatio`, `CntRatio`, `epsRatio` by dividing two `ctrMonomial[...]` / `cntMonomial[...]` / `epsEtaMonomial[...]` calls evaluated at scaled vs. unscaled arguments (the same way SymPy builds them from the helper functions and `ratio_subs`), then take the Jacobian with `D[]`. Do NOT hardcode the pre-reduced `(rg rc/rU)^(1+deltaUs)(rT/rU)^(1+chi0s)` form at math:118-130.
- Section III: see F2 — reconstruct the actual-vs-target ratios from the monomial helper functions with `TActual/KetaActual/muActual` substituted, dividing by an independent `CtrTarget/CntTarget/epsEtaTarget` symbol, exactly as SymPy does at sympy:304-311. This both removes the tautology and supplies the independent route.

**Verification:**
After the rewrite, math Section I should show `CtrRatio` built from `cntMonomial`/`ctrMonomial` quotients (not the hand-collapsed literal), and math Section III should show `CtrActualRatio = normalizeExpr[ctrMonomial[gf, cf, KUf, TActual, chi0s, deltaUs, L]/CtrTarget]` (not `(TActual/TOrbit)^(1+chi0s)`). All five sections still exit 0 with identical printed residuals (all-zero), and the engine cross-check still passes.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl:236-245`

**What's wrong:**
In Mathematica Section III the three actual/target ratios are written in *already-collapsed* form built directly from the `m·orbit` definitions, so the three `expectZero` checks are algebraic identities about `Log` of a power, not tests of the monomials' physics:

```
math:232  TActual    = normalizeExpr[mT TOrbit];
math:233  KetaActual = normalizeExpr[mK KetaOrbit];
math:234  muActual   = normalizeExpr[mMu muOrbit];
math:236  CtrActualRatio = normalizeExpr[(TActual/TOrbit)^(1 + chi0s)];
math:238  epsEtaActualRatio = normalizeExpr[KetaOrbit/KetaActual];
math:240  expectZero["...", Log[CtrActualRatio] - (1 + chi0s) Log[mT]];
math:245  expectZero["...", Log[epsEtaActualRatio] + Log[mK]];
```

Since `TActual = mT·TOrbit`, `CtrActualRatio` reduces to `mT^(1+chi0s)` by construction, so `Log[CtrActualRatio] - (1+chi0s)Log[mT] = (1+chi0s)Log[mT] - (1+chi0s)Log[mT] = 0` independently of any monomial structure. Likewise `epsEtaActualRatio = KetaOrbit/KetaActual = 1/mK` by the `KetaActual=mK·KetaOrbit` definition, so the check is `Log[1/mK]+Log[mK]=0` — trivially true. `CntActualRatio` (math:237) is assembled the same way from the collapsed pieces. None of these checks exercises the chart-conversion claim (paper eq:app-part05-chart-conversions-1/2) that the *Ctr/Cnt/epsEta monomials' exponent structure* yields `R_{tr}=m_T^{1+\chi_0}`, `R_\eta=1/m_K`, `R_{nt}=m_\mu/(m_K m_T^{F_*})`.

Contrast SymPy Section III (sympy:304-315), which rebuilds the full monomial with the actual values substituted and divides by an *independent* target symbol:
```
sympy:304  Ctr_actual_ratio = simplify_expr(ctr_monomial(gf, cf, KUf, T_actual, chi0s, deltaUs, L) / Ctr_target)
sympy:313  expect_zero("...", Ctr_actual_ratio - mT ** (1 + chi0s))
```
Here `T_orbit` is solved (sympy:287-289) so that `ctr_monomial(...,T_orbit,...)=Ctr_target`, so the ratio genuinely tests that the Ctr monomial's T-exponent is `(1+chi0s)`. The SymPy check is substantive; the Mathematica check is tautological. Thus the paper's claim that the mismatch chart is "re-derived from the dependent-triple orbit solve rather than posited" (Ledger point 3) is verified by only one of the two engines.

**Why this matters:**
At a trust-anchor checkpoint the dual-engine requirement exists precisely so that the mismatch chart conversion laws (which downstream Part V citations rely on, MTDC-T9.6) are confirmed independently. As written, the Mathematica engine confirms a `Log[a^b]=b·Log[a]` identity, not the chart law, so a sign or exponent error in the underlying Cnt/Ctr monomial conversion would not be caught by the second engine.

**Required change:**
Rewrite math Section III to reconstruct the ratios from the monomial helper functions, exactly mirroring SymPy's substantive route (this simultaneously discharges F1 for Section III):
- Replace math:236-238 with:
  - `CtrActualRatio = normalizeExpr[ctrMonomial[gf, cf, KUf, TActual, chi0s, deltaUs, L] / CtrTarget];`
  - `CntActualRatio = normalizeExpr[cntMonomial[lamf, gf, KUf, KetaActual, KWf, muActual, TActual, eStar, fStar, L, sigma] / CntTarget];`
  - `epsEtaActualRatio = normalizeExpr[epsEtaMonomial[cf, KUf, KetaActual] / epsEtaTarget];`
- Keep `TOrbit`, `KetaOrbit`, `muOrbit` (math:222-230) as the orbit solve, and keep `TActual/KetaActual/muActual` (math:232-234). The `CtrTarget/CntTarget/epsEtaTarget` symbols already exist in `$Assumptions` (math:110) so no new symbol declaration is needed.
- Leave the three `expectZero` targets (math:240-245) unchanged in form; they will now test the genuine monomial conversion.

**Verification:**
After the change, math:236 reads `CtrActualRatio = normalizeExpr[ctrMonomial[gf, cf, KUf, TActual, chi0s, deltaUs, L]/CtrTarget]` (monomial-derived, not pre-collapsed), and the three Section III `expectZero` checks still print `= 0` / `PASS`. The script exits 0 and the printed `q(mismatch)` matrix is unchanged from the committed output.

## Independent-derivation check (Mathematica)

The `.wl` is NOT an independent native re-derivation; it is a transliteration of the `.py` (F1). It does use Mathematica-native `D[]` (math:156) and `Series[]` (math:280-291) in Sections I and V, but every section follows SymPy's exact derivation route, reuses SymPy's intermediate variable names, and compares against the same hardcoded `Mexpected`/`Dmis` literals. Section III goes further and degrades into tautology (F2) by writing the actual/target ratios in pre-collapsed form. Corresponding-section evidence:
- math:118-130 hardcodes `CtrRatio = (rg rc/rU)^(1+deltaUs)(rT/rU)^(1+chi0s)` vs sympy:118-119 which divides two `ctr_monomial` calls.
- math:236-238 pre-collapses the Section III ratios vs sympy:304-311 which rebuilds monomials and divides by independent target symbols.
- math:157-161 `Mexpected` is byte-for-byte the same literal as sympy:155-161.

## Engine cross-check

Both engines exit 0 and produce matching final residuals where the checks are genuine:
- Section I: both print `derived M_* - carried Stage 192 matrix = 0` (sympy output L63-68; math output L27-28) and `q^(2<-1) - M_* Δx = 0` (sympy L69-74; math L29-30). Agree.
- Section II: both print three witness `... = 0` and `pairwise - intrinsic = 0` (sympy L79-185; math L35-46). Agree.
- Section III: both report `= 0` (sympy L190-213; math L51-62), but the Mathematica zeros are tautological (F2), so the *agreement* in Section III is not independent confirmation.
- Section IV: both print all-zero cocycle residuals (sympy L218-239; math L67-70). Agree.
- Section V: both print `Delta_Q^lin - ... = 0` (sympy L244; math L75-76) and matching `Delta_H^lin`. Agree.
No `engine_disagreement`: where both engines compute the same quantity, the residuals match. The defect is not disagreement but that the Mathematica Section III agreement is vacuous.

## Verdict justification

`findings`. The SymPy script is substantive throughout: I hand-verified the Section I Jacobian rows (Ctr row gives exactly `[0,1+δ,1+δ,−(2+χ₀+δ),0,0,0,1+χ₀]`), the Section II derived dependent-triple multipliers `Φ_K=rc²/rU`, `Φ_T=rU(rU/(rg rc))^α` with `α=(1+δ)/(1+χ₀)` (these genuinely enforce monomial invariance — not tuned tautologies), the Section III monomial re-derivation through an independent target symbol, and the Section V series expansion of `χ_Q=3(Sβ⁵+9Σ₅)/(3S−Σ₀)` reproducing `ε(5ε_β+δΣ₀/(3S)+9δΣ₅/S)`. All of SymPy holds up and aligns exactly with the paper card and appendix. The two findings are Mathematica-side: (F1) the `.wl` is a transliteration, not an independent route — disqualifying at a checkpoint's higher bar; and (F2) its Section III checks are tautological (`Log[mT^{1+χ₀}]−(1+χ₀)Log[mT]≡0`, `Log[1/mK]+Log[mK]≡0`), so the second engine does not independently confirm the chart-conversion law that MTDC-T9.6 hands downstream. Paper alignment is exact; the stale "251"/"183" labels are cosmetic and out of scope for Codex. No stop-cold: the SymPy proof stands, the fix is a localized Mathematica rewrite that does not change any derived constant, so nothing downstream is invalidated.

## Self-test notes

- Variable independence (trap #1): the only `D[]`/`diff` ops are the Section I/V Jacobian over `Dvec` and the Section V `Series` in `eps`; in both engines the differentiated expression genuinely depends on those variables (the log-ratios are linear in the `D`'s; `chi_from_def` depends on `eps` after `.subs`). No vacuous-derivative trap. My proposed F2 fix introduces no new `D[]`.
- Trivial-case / hardcoded check (trap #3): I confirmed F2 is tautological by substituting the definitions `TActual=mT·TOrbit` etc. and seeing the residual collapse to `0` for any monomial, and confirmed the SymPy counterpart is non-trivial because `Ctr_target` is an independent symbol satisfying the orbit solve. The proposed math Section III fix uses only pre-existing symbols (`CtrTarget` etc. already in `$Assumptions` at math:110), so it will run without new declarations.
- Paper round-trip (trap #5): the F2 fix reproduces the same `R_{tr}=m_T^{1+\chi_0}`, `R_\eta=1/m_K`, `R_{nt}=m_\mu/(m_K m_T^{F_*})` conversions the paper states (eq:app-part05-chart-conversions-1/2), introducing no new constant and no new `paper_misalignment`.
