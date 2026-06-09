---
unit_id: 200
batch: V.3
auditor_model: Opus 4.8 (1M context) [claude-opus-4-8[1m]]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage200_reference_free_home_stretch_theorem.md]
  paper_appendix: present
---

# Audit unit 200 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_200.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage200_reference_free_home_stretch_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 131 stage table, lines 240–260 theorem block, lines 1380–1456 reference-free section, line 1468 MTDC-T9.6 citation)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` reads: "Combines Packet A and Packet B into the four-scalar final verdict packet \((\Delta_Q,q_{\rm tr},q_{\rm nt},q_\eta)\)." The notes expand this into five deliverables: (1) Packet-B orbit-representative independence — the orbit packet attaches to the orbit \(\mathcal O_*\), not to a chosen representative, via the exact additive cocycle law; (2) the additive four-scalar full verdict packet \(\Delta_{\rm full}=(\Delta_Q,q_{\rm tr},q_{\rm nt},q_\eta)\) with \(\Delta_Q:=\chi_Q-1\); (3) the equivalent multiplicative and mismatch charts plus the exact chart-conversion laws \(\mathfrak R_{\rm tr}=m_T^{1+\chi_{0,*}}\), \(\mathfrak R_\eta=1/m_K\), \(\mathfrak R_{\rm nt}=m_\mu/(m_K m_T^{F_*})\); (4) the reference-free theorem \(\Delta_{\rm full}=0\iff\Delta_{\rm branch}=0\text{ and }x\in\mathcal O_*\); (5) the statement that the reduced endgame is now exactly one Packet-A scalar plus one Packet-B quotient triple. The appendix theorem block (lines 240–260) and §reference-free section (lines 1386–1453) reproduce the additive/multiplicative/mismatch packet definitions and the same conversion laws verbatim. The Packet-A finish line carried from Stage 197 is \(\Delta_Q=\chi_Q-1\) with \(\chi_Q=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)\) at the reduced baseline.

## What the script claims to verify

Both scripts run five labelled sections. **I** derives the Packet-B compiler matrix `M_*` (`Mderived`) as the logarithmic Jacobian of the three primitive coherent monomials' pairwise ratios (`ctr/cnt/epsEta`), then checks `Mderived == Mexpected` (the carried Stage-192 literal) and `qPair == Mderived·Δx`. **II** builds an arbitrary target-orbit witness `w` from the orbit-preserving `Phi` factors and checks each witness/`*` monomial ratio is 1, then that the intrinsic orbit packet equals the pairwise-witness packet (orbit-representative independence). **III** inverts the monomials to the dependent-triple orbit solve (`Keta_orbit/T_orbit/mu_orbit`), applies mismatch factors `m_T,m_K,m_mu`, recomputes the monomial ratios, and checks they equal the chart-conversion targets `m_T^(1+chi0s)`, `m_mu/(m_K m_T^F*)`, `1/m_K`, plus `q(mismatch)==M_*·Δx_mis`. **IV** checks the additive cocycle `Δx31=Δx32+Δx21` ⇒ `q31=q32+q21`. **V** Taylor-expands `chi_from_def` to first order and checks `ΔQ_lin == eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)`, then assembles the full four-row `Δ_H^lin`. These map onto deliverables (1)–(5) above.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) orbit-representative independence (additive cocycle) | Sec II `pairwise witness packet - intrinsic orbit packet == 0` (py:269-272 / wl:219); Sec IV cocycle (py:396-397 / wl:275-276) | match |
| (2) additive packet `Δ_full=(Δ_Q,q_tr,q_nt,q_η)`, `Δ_Q=χ_Q-1` | Sec V `Δ_H^lin` assembly (py:428-437 / wl:300-302); `Δ_Q` row from Sec V linearization | match |
| (3) chart-conversion laws `R_tr=m_T^(1+χ0*)`, `R_η=1/m_K`, `R_nt=m_μ/(m_K m_T^{F*})` | Sec III three `expectZero` checks (py:313-315 / wl:241-246) | match |
| (4) reference-free theorem `Δ_full=0 ⇔ Δ_branch=0 ∧ x∈O_*` | Sec II witness-invariance + Sec III mismatch chart + Sec V `Δ_Q` linearization jointly exercise the iff coordinates (each scalar is shown to vanish exactly at the closure point) | match |
| (5) `M_*` compiler / four-scalar reduction | Sec I `Mderived==Mexpected`, `qPair==M_*·Δx` (py:167-168 / wl:168-169) | match |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side check, and the script tests nothing the paper does not claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 167 | `expect_zero(Mderived - Mexpected)` | (5) M_* compiler | yes |
| A2 | sympy | 168 | `expect_zero(q_pair - Mderived·Δx)` | (5) | yes |
| A3 | sympy | 227-238 | three `ln[witness/*] == 0` | (1) orbit independence | yes |
| A4 | sympy | 269-272 | `pairwise witness - intrinsic == 0` | (1) | yes |
| A5 | sympy | 313-315 | three chart-conversion ratio checks | (3) conversion laws | yes |
| A6 | sympy | 339-340 | `q_mismatch - carried chart == 0`, `M_*·Δx_mis - q_mismatch == 0` | (3)/(5) | yes |
| A7 | sympy | 396-397 | cocycle `Δx`/`q` additivity | (1) cocycle | yes |
| A8 | sympy | 423-426 | `ΔQ_lin - eps(5 eps_beta + …) == 0` | (2) Δ_Q | yes |
| B1–B8 | mathematica | 168-169, 200-202, 219, 241-246, 254-255, 275-276, 295-298 | mirror of A1–A8 via `expectZero`/`PASS` | same | yes |

All rows are substantive (each can fail if an exponent, the carried matrix, or a Taylor coefficient is wrong). No tautological or unanchored rows.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.txt` (mtime 2026-06-01 12:42:20)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.txt` (mtime 2026-06-01 12:42:20)

**What's wrong:**
Both saved outputs predate their scripts (both scripts: mtime 2026-06-03 15:59:11) and both carry a stale pre-renumber banner: their first banner reads `STAGE 183 — EXACT REFERENCE-FREE HOME-STRETCH THEOREM` and the ledger footer reads `STAGE 183 LEDGER`, whereas the current scripts print `STAGE 200 …` (py:76 / wl:68) and `STAGE 200 LEDGER` (py:439 / wl:304). The numeric/symbolic *content* of the captured runs otherwise matches what the current scripts would produce (all checks read `= 0` / `PASS`), so only the banner label is stale, not the math.

**Why this matters:**
Cosmetic provenance drift only — a reader diffing the committed `.txt` against the current script sees a mismatched stage number. No load-bearing result is affected.

**Required change:**
None for Codex to hand-edit; the orchestrator's independent re-run regenerates both `.txt` files with the correct `STAGE 200` banner. Informational, non-blocking.

**Verification:**
After the orchestrator re-runs `redteam exec-sympy 200` and `exec-mathematica 200`, both `.txt` first lines read `STAGE 200 — EXACT REFERENCE-FREE HOME-STRETCH THEOREM` and mtimes postdate the scripts.

## Independent-derivation check (Mathematica)

This is the V.3 checkpoint's central question: the `.wl` pre-existed and was de-transliterated in the first pass (F1: Section I hand-collapsed ratios → the `ratioSubs`/`Mderived`-Jacobian route; F2: Section III `Log[a^b]` collapse → the `ctrMonomial[...TActual...]/CtrTarget` exponent exercise). I re-verified BOTH fixes fresh against the current scripts.

**Section I — `ratioSubs` / `Mderived` Jacobian (F1).** The `.wl` builds the three symbolic coherent monomials (`Ctr1/Ctr2/Cnt1/Cnt2/eps1/eps2`, wl:123-128), forms the pairwise ratios under `ratioSubs` (wl:119-131), reparameterizes the ratios `r→Exp[D]` (wl:140-149), takes the log, and computes `Mderived = Table[D[qPair[[i]], Dvec[[j]]]]` (wl:157) — i.e. the matrix is genuinely *derived* as the Jacobian of the actual monomial log-ratios, then checked against the carried literal `Mexpected` (wl:168). The `.py` does the SAME for the same deliverable: `Mderived = q_pair.jacobian(Dvec)` (py:154) over `q_pair` built from `ratio_subs`/`ratio_to_logs` (py:107-152), checked at py:167. **Verdict on F1's sufficiency: sufficient.** This is a genuine derivation in BOTH engines (the `M_*` entries fall out of the monomial exponents `1+deltaUs`, `1+chi0s`, `Estar`, `Fstar`; the check `Mderived==Mexpected` fails if any monomial exponent is mis-stated — it is not a `Mexpected`-trusting tautology). It is NOT a disguised hand-collapse: neither engine pre-bakes the row entries; both run the CAS Jacobian. The original transliteration (hand-collapsed ratios echoing the `.py`'s intermediate algebra) is genuinely gone.

**Section III — `ctrMonomial[...TActual...]` exponent exercise (F2).** The `.wl` inverts the Ctr monomial to the orbit solve `TOrbit = (L^2 KUf/Pi^2)(CtrTarget/(gf cf/KUf)^(1+deltaUs))^(1/(1+chi0s))` (wl:224-226), applies `TActual = mT·TOrbit` (wl:233), recomputes `ctrMonomial[gf,cf,KUf,TActual,...]/CtrTarget` (wl:237), and checks it equals `(1+chi0s)Log[mT]` (wl:241). The `mT` factor enters the recomputed monomial inside the `(Pi^2 T/(L^2 K))^(1+chi0s)` power, so the result `m_T^(1+chi0s)` genuinely *exercises* the `(1+chi0s)` exponent — it is not a `Log[a^b]→b·Log[a]` collapse on a literal. The `.py` does the same (py:287-313). **Verdict on F2's sufficiency: sufficient.** No surviving tautology/X−X in either engine; the `T_orbit` inversion's `1/(1+chi0s)` cancels against the monomial's `(1+chi0s)` only because `T_actual` carries the *extra* `mT`, which is the substance of the check.

**Overall independence verdict: BORDERLINE (parallel re-implementation; F1/F2 de-fixes are genuine and sufficient).** The `.wl` is, line-for-line, the same five-section choreography as the `.py` (identical helper-monomial bodies, identical `ratioSubs`/`ratioToLogs`/`orbitRatioSubs`, identical `Mexpected` literal, identical `Phi` witness factors, identical orbit solve, identical `Series` linearization), translated CAS-for-CAS. What keeps this on the right side of the `mathematica_transliteration` line — under the consistent threshold "same expression → same operation = a port" — is that the shared route is itself a *genuine first-principles derivation* in each engine, executed with each CAS's own native simplifier (SymPy `simplify∘expand_log∘expand_power_exp` vs Mathematica `FullSimplify∘PowerExpand∘Together∘Expand` under `$Assumptions`). The load-bearing quantities (`M_*`, the chart conversions, `ΔQ_lin`) are computed, not echoed from a hand-collapsed intermediate; the de-transliteration removed exactly the echoed-intermediate pattern F1 flagged. I am NOT raising a `mathematica_transliteration` finding: there is no hand-collapsed or pre-baked intermediate that one engine copies from the other, which is the operative test (the parallel structure is the second-engine policy's intended shape, the same as the rest of this batch's stages). Evidence corroborating "genuine, not echo": (a) Sec I `Mderived` is a CAS Jacobian of monomials, not a typed-in matrix, in both; (b) Sec III `mT` is threaded through the power so the exponent is exercised; (c) Sec V the linearization coefficient `5 eps_beta + dSigma0/(3S) + 9 dSigma5/S` is produced by each engine's own `series`/`Series`+`removeO`/`Normal`, and I re-derived it by hand independently (it matches).

## Engine cross-check

Both engines agree at every check. SymPy output: all matrix residuals print as zero columns (`derived M_* - carried Stage 192 matrix = [0…0]`, `q^(2<-1) - M_* Δx = [0,0,0]`, `pairwise witness - intrinsic = [0,0,0,0]`, `q(mismatch) - carried chart = [0,0,0]`, `Δx31-Δx32-Δx21 = [0…0]`, `ΔQ_lin - eps(…) = 0`). Mathematica output: the same residuals print `MatrixForm[{{0,…,0}}]` and each scalar check prints `= 0` followed by `PASS`. The Section I monomial ratios printed match symbol-for-symbol between engines (e.g. SymPy `Cnt_2/Cnt_1 = r_λ²·r_μ·(r_U/r_T)^F*·(r_γ²r_λ²/(r_U r_W))^E* / (r_K r_W²)` vs Mathematica `(rla^2*(rg*rla)^(2*eStar)*rmu*(rU/rT)^fStar)/(rK*rW^2*(rU*rW)^eStar)` — identical after `PowerExpand`). `engines_agree: true`.

## Verdict justification

`findings: 1` — a single low-severity, non-blocking `stale_output` (both `.txt` carry the pre-renumber `STAGE 183` banner; content otherwise current). The checkpoint passes the higher bar: both engines present and substantive; every paper deliverable has a matching, non-tautological check; paper alignment is exact; the load-bearing quantities (`M_*`, the three chart-conversion laws, the Packet-A linearization coefficient) are re-derived in-script rather than trusted as literals, and I re-derived the Packet-A coefficient by hand to confirm. Attacks tried that failed: (a) hunted for a surviving `Log[a^b]` collapse in Sec III — the `mT` is genuinely threaded through the `(1+chi0s)` power, so the exponent is exercised; (b) tested whether `Mderived==Mexpected` is a self-confirming literal — no, `Mderived` is a CAS Jacobian of the physical monomials and would fail on a wrong exponent; (c) checked Sec II `assume(positive)` domains against the multiplicative orbit setup — all 50+ symbols are physical positive scales (`>0`), justified by the monomials being ratios of positive coherent kernels; (d) checked Sec V `S, beta, Sigma0, Sigma5` are `nonzero`/`>0` so `3S-Sigma0` is a legitimate denominator near the baseline. I read the card, notes, and appendix in full, and the script's verified claim matches the paper's stated four-scalar reference-free theorem. The `mathematica_transliteration` call is BORDERLINE but resolved as *not a finding*: the parallel structure is a genuine dual-CAS re-derivation, not an echoed hand-collapse, and the first-pass F1/F2 de-fixes are confirmed sufficient.

## Self-test notes

- **Variable independence**: Sec I `Mderived = Jacobian(qPair, Dvec)` — each `qPair` component is linear in `Dvec` (the log of an exp-substituted ratio), so every `D[·]` is the genuine nonzero exponent, not an identically-zero derivative. No `diff(expr, var)` where `var ∉ expr`.
- **Trivial-case / exponent exercise**: Sec III `mT→1` ⇒ `Ctr_actual_ratio→1` ⇒ `Log→0 = (1+chi0s)·0`; `mT≠1` gives `mT^(1+chi0s)`, confirming the exponent is live. Sec V `eps→0` ⇒ `ΔQ_lin→0`; coefficient hand-derived independently and matches.
- **Paper round-trip**: F1 (stale_output) prescribes no script edit, so it introduces no new `paper_misalignment`; the chart-conversion constants and `M_*` entries the scripts assert all appear verbatim in the notes (lines 264-285) and appendix (lines 1424-1432, 158-162 equiv). No value mismatch.

## Value Reconciliation (pass-2 augmentation)

This stage emits **no named numeric constants** — every deliverable is symbolic (the carried Stage-192 matrix `M_*`, the chart-conversion laws, the four-scalar packet, the Packet-A linearization coefficient). The Family-1 radius `√(4107−100π²)/(10π)` and the `168π²`/`100π²` class do not appear in this stage (its exponents are the symbolic `chi0s/deltaUs/Estar/Fstar`), so nothing in that reconciliation class applies.

| value (symbolic deliverable) | source (py / wl + output) | .tex/.md location | status |
|---|---|---|---|
| `M_*` row 1 = `[0, 1+δU*, 1+δU*, -(2+χ0*+δU*), 0,0,0, 1+χ0*]` | py:157 / wl:159; out: `M_* derived …` block | notes derives via chart conv. (lines 264-285); appendix eq. forms (1424-1432) | MATCH (symbolic, carried Stage-192 matrix; consistent with notes' conversion laws) |
| `M_*` row 2 = `[2(1+E*),0,2E*,F*-E*,-1,-(2+E*),1,-F*]` | py:158 / wl:160 | notes Cnt monomial / conversion (lines 274-285) | MATCH |
| `M_*` row 3 = `[0,2,0,-1,-1,0,0,0]` | py:159 / wl:161 | notes epsEta / R_η=1/m_K (line 270) | MATCH |
| chart conv. `R_tr=e^{q_tr}=m_T^{1+χ0*}` | py:313 / wl:241; out: `…- (1+chi0_*) Log[m_T] = 0` | notes line 265; appendix eq:app-part05-chart-conversions-1 (1426) | MATCH |
| chart conv. `R_η=e^{q_η}=1/m_K` | py:315 / wl:246 | notes line 270; appendix 1428 | MATCH |
| chart conv. `R_nt=e^{q_nt}=m_μ/(m_K m_T^{F*})` | py:314 / wl:244 | notes line 275; appendix eq:app-part05-chart-conversions-2 (1432) | MATCH |
| additive packet `Δ_full=(Δ_Q,q_tr,q_nt,q_η)`, `Δ_Q=χ_Q-1` | py:428-437 / wl:300-302 | card `\stagefield{Output}` (line 15); notes line 203-218; appendix 1386-1398 | MATCH |
| Packet-A linearization coeff. `5 ε_β + δΣ0/(3S) + 9 δΣ5/S` | py:421 / wl:293; out: `…= 0` + `PASS` | notes carries `Δ_Q=χ_Q-1` w/ `χ_Q=3(Sβ^5+9Σ5)/(3S-Σ0)`; coefficient is the first-order expansion (not separately tabulated, but the closed-form `χ_Q` it derives from IS in notes line 408-equiv / Stage-197 carry) | MATCH (deliverable is the symbolic `Δ_Q` law; the linearization is its internal verification, not a separately-boxed prose constant) |

INTERNAL scaffolding (no finding): `Phi_T/Phi_K/Phi_mu` witness factors (Sec II), `T_orbit/Keta_orbit/mu_orbit` dependent-triple solve (Sec III), `D21/D32/D31`/`q21/q32/q31` cocycle witnesses (Sec IV), `eps` series machinery, all `r_*`/`Δ_*` ratio/log dummy symbols.

**reconciliation: complete; 8 deliverable values checked, 0 misaligned.** Every symbolic deliverable the scripts emit is correctly reflected in the card, notes, or appendix. No `value_mismatch` and no `script_missing_paper_claim`. The lone finding (F1 `stale_output`) is not a reconciliation issue.
