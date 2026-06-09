---
unit_id: 199
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T18:51:04Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage199_pairwise_orbit_transport_law.md]
  paper_appendix: present
---

# Audit unit 199 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_199.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage199_pairwise_orbit_transport_law.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 129, 1329–1377, 1468)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage199_pairwise_orbit_transport_law_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage199_pairwise_orbit_transport_law_mathematica_audit.txt`

## What the paper claims

Stage 199 (`MTDC-T9.6`) gives the **exact pairwise / reference-independent orbit-transport law**, removing the last orbit-side base-point privilege from Part V. The card `\stagefield{Output}` states it "Removes representative privilege by giving exact two-point transport, mismatch, cocycle, and orbit-lock laws." The notes enumerate the deliverables precisely: (1) the two-point quotient packet `q^(2<-1) = M_* Δx^(2<-1)` with `q_tr=ln(C_tr,2/C_tr,1)`, `q_nt=ln(C_nt,2/C_nt,1)`, `q_η=ln(ε_η,2/ε_η,1)`, plus closed forms for the three coherent-monomial ratios; (2) the exact transport factors `Φ_T = r_U (r_U/(r_γ r_c))^α_*`, `Φ_K = r_c²/r_U`, `Φ_μ` (both factored and expanded-monomial forms), with `α_* = (1+δ_U,*)/(1+χ_0,*)`; (3) the reference-independent mismatch triple `m_T=r_T/Φ_T, m_K=r_K/Φ_K, m_μ=r_μ/Φ_μ` and the invariant-ratio collapse `C_tr ratio = m_T^(1+χ0), ε ratio = 1/m_K, C_nt ratio = m_μ/(m_K m_T^F)`; (4) the projector split `Δx = O_orb Δx + Q_quot Δx` with the explicit 8-vectors; (5) the cocycle/composition laws `Φ^31=Φ^32Φ^21`, `m^31=m^32m^21`, `q^31=q^32+q^21`; (6) the two-point orbit-lock theorem (5-way equivalence). The notes §4.2 add an exact restoration map, and §7 the reduction-to-Stage-198 specialization (free ratios = 1 ⇒ Φ=1 ⇒ m=raw ratio). The card's `\stagefield{Verification}` line says **"Mathematica audit: none yet"** even though a passing `.wl` exists (paper-side card-text lag — see F2).

## What the script claims to verify

The SymPy script (banner line 40) verifies the full deliverable set: §I checks the three coherent-monomial pairwise ratios match the closed forms (`.py:123-125`); §II builds `Φ_T,Φ_K,Φ_μ`, checks the factored vs expanded-monomial `Φ_μ` agree and that all three monomial ratios are 1 on the same orbit (`.py:154-174`); §III checks the mismatch-triple invariant-ratio collapse (`.py:192-203`); §IV builds `M_*` and the Stage-192 projectors (`Pdep.inv()` route) and checks the q-chart and the `O_orb/Q_quot` split against expected 8-vectors plus `M_* O_orb = 0` (`.py:227-287`); §V the restoration map (`.py:298-308`); §VI the cocycle laws for Φ, m, Δx, q (`.py:350-388`); §VII the reduction to Stage 198 (`.py:396-401`). The Mathematica script asserts the same deliverables but via independent routes: `compilerRows` built from primitive weight vectors (`.wl:98-120`), `Solve` of the same-orbit constraint for the transport factors (`.wl:144-176`), `LinearSolve`-of-augmented-system projectors with idempotency checks (`.wl:213-245`), plus an extra mismatch-to-q determinant/orbit-lock solve block (`.wl:303-317`).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| q-packet `q=M_*Δx`, three monomial-ratio closed forms (notes §1) | py §I `.py:123-125`; wl §I `.wl:120,133-135` | match |
| Transport factors Φ_T,Φ_K,Φ_μ + α_* (notes §2) | py §II `.py:154-174`; wl §II Solve `.wl:173-176,183` | match |
| Mismatch triple + invariant-ratio collapse (notes §3) | py §III `.py:192-203`; wl §III `.wl:195-200` | match |
| q-chart `q_tr=(1+χ0)ln m_T`, `q_η=-ln m_K`, `q_nt=...` (notes §4) | py §IV `.py:227-237`; wl §IV `.wl:209-211` | match |
| Projector split O_orb/Q_quot 8-vectors (notes §4.1) | py §IV `.py:283-287`; wl §IV `.wl:237-245` | match |
| Restoration map (notes §4.2) | py §V `.py:298-308`; wl §V `.wl:253-255` | match |
| Cocycle/composition Φ,m,q (notes §5) | py §VI `.py:350-388`; wl §VI `.wl:282-285` | match |
| Two-point orbit-lock theorem (notes §6) | py (via q-chart + reduction); wl §VII `.wl:296-317` | match |
| Reduction to Stage 198 (notes §7) | py §VII `.py:396-401`; wl §VIII `.wl:323-328` | match |

All paper-side deliverables have a faithful, non-tautological script-side check on both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 123-125 | `expect_zero(ratio_pairs - expected)` | monomial-ratio closed forms (§1) | yes |
| A2 | sympy | 154-157 | `expect_zero(ln Φ_μ - ln expanded)` | Φ_μ two forms (§2) | yes |
| A3 | sympy | 163-174 | `expect_zero(ln[ratio on same orbit])` | transport solves same-orbit (§2) | yes |
| A4 | sympy | 192-203 | `expect_zero(ratio - m-form)` | mismatch collapse (§3) | yes |
| A5 | sympy | 227-237 | `expect_zero(q_* - chart form)` | q-chart (§4) | yes |
| A6 | sympy | 283-287 | `expect_zero(O/Q vs expected, M_*O=0)` | projector split (§4.1) | yes |
| A7 | sympy | 298-308 | `expect_zero(ln restore - ln ΦT₁)` | restoration map (§4.2) | yes |
| A8 | sympy | 350-388 | `expect_zero(Φ/m/q cocycle)` | composition laws (§5) | yes |
| A9 | sympy | 396-401 | `expect_zero(Φ@free=1, m@free=raw)` | reduction to 198 (§7) | yes |
| B1 | mathematica | 120 | `expectZero[compilerRows - compilerTarget]` | M_* built from primitive weights (§1) | yes |
| B2 | mathematica | 133-135 | `expectZero[rowLog - ln closed]` | monomial ratios (§1) | yes |
| B3 | mathematica | 173-176,183 | `expectZero[Solve Φ - ln closed]` | transport from native Solve (§2) | yes |
| B4 | mathematica | 195-200 | `expectZero[ratio - m-form]` | mismatch collapse (§3) | yes |
| B5 | mathematica | 237-245 | `expectZero[O/Q vs expected, Q²-Q, QO]` | projector split + idempotency (§4.1) | yes |
| B6 | mathematica | 253-255 | `expectZero[restore - ΦT₁]` | restoration map (§4.2) | yes |
| B7 | mathematica | 282-285 | `expectZero[Φ/m/q cocycle]` | composition (§5) | yes |
| B8 | mathematica | 296-317 | `expectZero[lock solves agree, det, q=0⇒m=0]` | orbit-lock theorem (§6) | yes |
| B9 | mathematica | 323-328 | `expectZero[Φ@free=1, m@free=raw]` | reduction to 198 (§7) | yes |

No tautological rows. Every assertion compares an independently-derived object against the paper-stated closed form, or checks a structural identity (idempotency, kernel agreement) that can genuinely fail.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.txt:3,420` (mtime 2026-06-01 12:17:54)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:**
The committed SymPy `.txt` is older than the `.py` (script mtime 2026-06-03 > output mtime 2026-06-01), and its captured banner carries the stale pre-renumber stage label: `STAGE 182 — EXACT PAIRWISE ORBIT-TRANSPORT LAW…` (line 3) and `STAGE 182 LEDGER` (line 420), where 182 = 199 − 17. The **current** `.py:40` banner correctly reads `STAGE 199 …`, so the saved transcript predates the banner fix. This is the known V.3 batch pattern. The body math (all `= 0` residuals) is otherwise consistent with the current script — every check still passes and the printed forms match the current code. The Mathematica `.txt` is fresh and correctly labeled `STAGE 199` throughout.

**Why this matters:**
Cosmetic/provenance only; the orchestrator's independent re-run refreshes the transcript and the banner will then read 199. No math is affected.

**Required change:**
No code edit. Orchestrator re-runs `python3 scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py` to refresh `scripts/output/...sympy_audit.txt`; confirm the new transcript banner reads `STAGE 199`.

**Verification:**
After re-run, the output's line 3 reads `STAGE 199 — EXACT PAIRWISE ORBIT-TRANSPORT LAW…` and the ledger banner reads `STAGE 199 LEDGER`; all residuals remain `0`.

### F2 — paper-side card-text lag (informational; NOT a script finding)

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_199.tex:11`

**What's wrong:**
The card's `\stagefield{Verification}` reads "SymPy audit: `…stage199…sympy_audit.py`. **Mathematica audit: none yet.**" A passing, independent `.wl` (`mathematica/...stage199...mathematica_audit.wl`) and its fresh transcript now exist, so the "none yet" text is out of date.

**Why this matters:**
Prose-only currency lag; no script defect. Noted per the V.3 heads-up. Codex does not edit paper/; this is recorded for the user/orchestrator to update the card text (not a `paper_misalignment` requiring user resolution — the script and paper math agree; only the bookkeeping sentence lags).

**Required change:**
None for Codex. Orchestrator/user may update `paper/stages/stage_199.tex:11` to cite the Mathematica audit file. No directive entry is written for this item (paper-side, no math discrepancy).

**Verification:**
Card line 11 eventually names the `.wl` path; not gated by this audit.

## Independent-derivation check (Mathematica)

The `.wl` is **GENUINELY INDEPENDENT** of the `.py`, not a transliteration. Three load-bearing deliverables, each derived by a structurally different route:

1. **Compiler matrix M_*.** SymPy hardcodes the 3×8 exponent matrix as a literal (`.py:211-217`: `Mstar = sp.Matrix([[0, 1+deltaUs, 1+deltaUs, -(2+chi0s+deltaUs), 0,0,0, 1+chi0s], …])`). Mathematica instead **constructs** the rows by additively combining named primitive monomial weight vectors — `compilerRows = {(1+delU)*chiCoreWeights + (1+chi)*thermalWeights, nontrackPrefactorWeights + eStar*epsWWeights - fStar*thermalWeights, epsEtaWeights}` (`.wl:104-108`) with `chiCoreWeights={0,1,1,-1,0,0,0,0}` etc. (`.wl:98-102`) — then asserts `compilerRows - compilerTarget == 0` (`.wl:120`). The `.py` does not have this primitive-weight decomposition at all; the `.wl` builds M_* from physics-meaning monomial exponents.

2. **Transport factors Φ_T,Φ_K,Φ_μ.** SymPy writes the closed forms directly and verifies them: `PhiT = sp.simplify(rU*(rU/(rg*rc))**alpha_star)` (`.py:133`), `PhiK = rc**2/rU` (`.py:134`). Mathematica derives them from a **native same-orbit linear Solve**: `sameOrbitSolve = First[Solve[Thread[compilerRows.logVars == ConstantArray[0,3]], depLogVars, Reals]]` (`.wl:144-146`), `phiLogTemplate = depLogVars /. sameOrbitSolve` (`.wl:148`), then compares the *solved* `phiTLog/phiKLog/phiMuLog` against the SymPy closed forms (`.wl:173-176`). Different route: solve the constraint system for the dependent logs vs. assert the answer.

3. **Stage-192 projectors O_orb/Q_quot.** SymPy builds the dependent submatrix and inverts it: `Pdep = sp.Matrix.hstack(Mstar[:,T_idx], Mstar[:,Keta_idx], Mstar[:,mu_idx]); Pdep_inv = Pdep.inv(); Sdep = Edep*Pdep_inv; Qquot = Sdep*Mstar` (`.py:244-251`). Mathematica builds an **augmented constrained system and uses `LinearSolve` with unit-vector RHS**: `constrainedSystem = Join[compilerRows, freeSelector]; sectionFromSolve = Transpose[Table[LinearSolve[constrainedSystem, Join[UnitVector[3,col], ConstantArray[0,Length[freeSlots]]]], {col,3}]]; quotientProjector = sectionFromSolve.compilerRows` (`.wl:213-226`). Explicit submatrix inverse vs. augmented-system LinearSolve — different linear-algebra routes — and the `.wl` adds idempotency/orthogonality checks (`Q²-Q`, `O²-O`, `QO`, `OQ`, `.wl:242-245`) absent in the `.py`. The `.wl` §VII also adds a determinant of the mismatch-to-q chart (`Det[mismatchToQ]==1+chi`, `.wl:308-309`) and two cross-solve agreement checks (`.wl:296-301`) that have no SymPy counterpart.

The two engines share the same *naming and ordered basis* and ultimately compare against the same closed forms (which is correct — both must land on the paper's result), but the **derivation machinery differs in every load-bearing block**. This is the intended second-engine design, not a port. No `mathematica_transliteration` finding.

## Engine cross-check

Both transcripts report every check as zero / PASS. The SymPy `.txt` prints the symbolic forms (e.g. `Phi_K = r_c²/r_U`, `m_K = r_K r_U / r_c²`) and all `= 0` residuals; the Mathematica `.txt` prints the InputForm of the solved logs (`ln Phi_K from Solve = Log[rC^2/rU]`, `ln m_K = Log[(rK*rU)/rC^2]`) and `PASS` on every assertion including the extra idempotency block. The independently-derived transport factors, mismatch triple, projector 8-vectors, cocycle laws, and reduction-to-198 results agree between engines. `engines_agree: true`.

## Verdict justification

`findings` (both informational): one `stale_output` (the SymPy `.txt` predates the banner fix and carries the pre-renumber `STAGE 182` label; orchestrator re-run refreshes) and one paper-side card-text lag (`Verification` field still says "Mathematica audit: none yet" — not a math discrepancy, no Codex action). No script-side math defect was found. Attacks that failed: (i) tried to find a tautology — every assertion compares an independently-built object to the paper closed form or checks a structurally-failable identity (idempotency, kernel agreement, determinant), none are construction-guaranteed; (ii) tried to break independence — the `.wl` builds M_* from primitive weights, derives Φ from a native Solve, and builds projectors via augmented-system LinearSolve, so it is not a transliteration; (iii) checked symbol domains — all ratios and exponents declared `positive, real`, justified by the "positive microscopic states" premise (notes §1), and the `force=True`/`PowerExpand` log-expansions are valid under positivity; (iv) checked the q-chart sign/coefficient conventions (`q_tr=(1+χ0)ln m_T`, `q_η=-ln m_K`, `q_nt=ln m_μ-ln m_K-F_* ln m_T`) against notes §4 and appendix eq. `app-part05-log-chart-mismatches` — all match including the negative `q_η`; (v) checked α_* and the Φ_μ expanded-monomial exponents against notes §2.3 — exact match. I read the card, notes, and appendix; the script's verified claim matches the paper's claim exactly.

## Value Reconciliation (pass-2 augmentation)

Stage 199 emits **no numeric figures-of-merit**; every deliverable is a closed-form symbolic law in the exponent symbols `(χ_0,*, δ_U,*, E_*, F_*)` and the ratio symbols. The reconciliation is therefore symbolic-form matching against the notes (the natural carrier; the `.tex` card is deliberately terse and the appendix carries the chart/lock forms).

| value (symbolic deliverable) | source (py / wl) | .tex / .md location | status |
|---|---|---|---|
| `q = M_* Δx`, q_tr/q_nt/q_η = ln of monomial ratios | py:222-223, ledger 408-410; wl:204-207 | notes §1 (`q_tr=ln(C_tr,2/C_tr,1)…`), appendix:1366-1368 | MATCH |
| `C_tr ratio = (r_γ r_c/r_U)^(1+δ_U)(r_T/r_U)^(1+χ0)` | py:108, out:9-11; wl:122 | notes §1 boxed (lines 154-160) | MATCH |
| `C_nt ratio = (r_λ²r_μ/(r_K r_W²))(r_γ²r_λ²/(r_U r_W))^E (r_T/r_U)^-F` | py:109-113; wl:123-126 | notes §1 boxed (lines 169-176) | MATCH |
| `ε_η ratio = r_c²/(r_U r_K)` | py:114, out:21-25; wl:127 | notes §1 boxed (lines 162-167) | MATCH |
| `α_* = (1+δ_U,*)/(1+χ_0,*)` | py:132; wl:154 | notes §2 boxed (lines 195-198) | MATCH |
| `Φ_T = r_U(r_U/(r_γ r_c))^α` | py:133, out:33-39; wl:156 | notes §2.1 boxed (lines 204-210) | MATCH |
| `Φ_K = r_c²/r_U` | py:134, out:40-44; wl:157 | notes §2.2 boxed (lines 219-225) | MATCH |
| `Φ_μ` (factored + expanded monomial form) | py:135-146, out:45-66; wl:158-167 | notes §2.3 boxed (lines 234-253) | MATCH |
| `m_T=r_T/Φ_T, m_K=r_K/Φ_K, m_μ=r_μ/Φ_μ` | py:181-183, out:75-104; wl:187-189 | notes §3 boxed (lines 284-291), appendix:1342 | MATCH |
| invariant-ratio collapse `C_tr=m_T^(1+χ0), ε=1/m_K, C_nt=m_μ/(m_K m_T^F)` | py:192-203; wl:195-200 | notes §3 boxed (lines 296-314), appendix:1426-1432 | MATCH |
| q-chart `q_tr=(1+χ0)ln m_T, q_η=-ln m_K, q_nt=ln m_μ-ln m_K-F ln m_T` | py:227-237; wl:209-211 | notes §4 boxed (lines 346-352), appendix:1347-1351 | MATCH |
| O_orb / Q_quot 8-vectors | py:262-281, out:309-340; wl:231-232 | notes §4.1 boxed (lines 371-399) | MATCH |
| restoration map (T,K,μ restore = Φ·base) | py:294-308; wl:249-255 | notes §4.2 boxed (lines 407-432) | MATCH |
| cocycle `Φ^31=Φ^32Φ^21, m^31=m^32m^21, q^31=q^32+q^21` | py:350-388; wl:282-285 | notes §5 boxed (lines 454-481), appendix:1377 | MATCH |
| reduction-to-198 `Φ=1 ⇒ m=raw ratio` at free=1 | py:396-401; wl:323-328 | notes §7 boxed (lines 558-570) | MATCH |
| `det(mismatch-to-q chart) = 1+χ0` | wl:308-309, out:112 | implicit in notes §4 chart (Jacobian of the 3×3 chart map) | MATCH (notes carry the chart; determinant is its consistency invariant) |

INTERNAL scaffolding (no prose expected): pass/fail flags, `expect_zero/expectZero` residuals (all 0), the `L, sigma, pi` placeholder constants in the `.py` monomial definitions (cancel in every ratio), the `baseT/baseK/baseMu` placeholders in the `.wl` restoration block, intermediate `Pdep/Sdep/Edep` matrices, `constrainedSystem/sectionFromSolve`, idempotency residuals (`Q²-Q` etc.).

`reconciliation: complete; 16 deliverables checked, 0 misaligned.` Every emitted symbolic deliverable reconciles against the notes (and, where applicable, the appendix). The card `.tex` is terse by design and omits the closed forms (legitimately — they live in the notes/appendix), so no MISSING-DELIVERABLE finding. The only freshness caveat is the stale SymPy `.txt` (F1), which still prints forms consistent with the current `.py`.

## Self-test notes

Checked: (1) variable-independence — there are no `sp.diff`/`D[]` derivatives in either script; all checks are algebraic log-identity residuals, so the "zero-derivative trivial-pass" trap does not apply. (2) No unbounded integrals → no parity trap. (3) Trivial-case pre-check — substituted free ratios = 1 (the script's own §VII/§VIII): `Φ_T,Φ_K,Φ_μ → 1` and `m_T,m_K,m_μ → r_T,r_K,r_μ`, matching notes §7; and the same-orbit substitution drives all three monomial ratios to 1 as asserted. (4) No new `missing_verification_script` directive (both engines present), so no path-spec concern. (5) Paper round-trip — the only directive item (F1) is a no-op re-run, introducing no new constant; the paper-side card-text lag (F2) is informational with no Codex edit. Conclusion: independence holds, alignment is exact, the two findings are both informational; a directive is written carrying only the F1 stale-output re-run note.
