---
unit_id: 199
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage199_pairwise_orbit_transport_law.md]
  paper_appendix: present
---

# Audit unit 199 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_199.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage199_pairwise_orbit_transport_law.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 129, 1329-1377, 1468 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The stage card's `\stagefield{Output}` states: "Removes representative privilege by giving exact two-point transport, mismatch, cocycle, and orbit-lock laws." The notes enumerate six deliverables for the pairwise (state-1 → state-2) setting: (1) the exact two-point quotient packet `q^(2<-1) = M_* Delta x^(2<-1)` with `q_tr=ln(Ctr2/Ctr1)`, `q_nt=ln(Cnt2/Cnt1)`, `q_eta=ln(eps2/eps1)`; (2) the exact pairwise transport factors `Phi_T = r_U (r_U/(r_g r_c))^alpha`, `Phi_K = r_c^2/r_U`, `Phi_mu = Phi_K r_W^2/r_la^2 (r_g^2 r_la^2/(r_U r_W))^{-E_*}(Phi_T/r_U)^{F_*}` with `alpha_* = (1+deltaU_*)/(1+chi0_*)`; (3) the reference-independent mismatch triple `m_T=r_T/Phi_T`, `m_K=r_K/Phi_K`, `m_mu=r_mu/Phi_mu` and the invariant-ratio collapse `Ctr2/Ctr1=m_T^{1+chi0_*}`, `eps2/eps1=1/m_K`, `Cnt2/Cnt1=m_mu/(m_K m_T^{F_*})`; (4) the projector split `Delta x = O_orb Delta x + Q_quot Delta x` (with `Q_quot Delta x` supported only on the dependent triple); (5) the multiplicative/additive cocycle laws `Phi^31=Phi^32 Phi^21`, `m^31=m^32 m^21`, `q^31=q^32+q^21`; (6) the two-point orbit-lock theorem `x^(2) in G_*.x^(1) <=> m_T=m_K=m_mu=1 <=> q_tr=q_nt=q_eta=0` (appendix eq. \ref{eq:app-part05-two-point-orbit-lock}, lines 1372-1376), plus the Stage-198 reduction (free ratios all 1 => `Phi=1`, `m=r`). Status is `\StatusExactClosure{}`; the stage introduces no new constitutive law.

## What the script claims to verify

The SymPy script verifies, over positive real symbols, the full pairwise-transport algebra in seven sections matching the notes one-to-one. Section I builds `Ctr/Cnt/eps` from the primitive coordinate monomials, substitutes `state2 = ratio * state1`, and checks the resulting pairwise ratios equal the hand-written "expected" pairwise monomial forms. Section II defines the three transport factors and checks (a) `Phi_mu`'s factored and expanded monomial forms agree, and (b) substituting each transport factor back into the corresponding ratio yields 1 (the same-orbit condition). Section III defines the mismatch triple and verifies the three invariant-ratio collapse identities. Section IV applies the linear map `M_*` to the log-ratio vector, cross-checks the q-chart relations against the transport-derived mismatches, builds the orbit/quotient projectors from `M_*`'s dependent-column inverse, and checks the projector split, the support pattern of `Q_quot Delta x`, and `M_* O_orb Delta x = 0`. Section V checks the restoration map; Section VI checks the additive/multiplicative cocycle laws with independent 2<-1 and 3<-2 ratio symbols; Section VII checks the Stage-198 reduction.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Quotient packet q = M_* Delta x; q_tr/q_nt/q_eta as log invariant-ratios | Section I (lines 116-125), Section IV (lines 219-238) | match |
| (2) Transport factors Phi_T, Phi_K, Phi_mu (incl. expanded monomial form) | Section II (lines 132-157) | match |
| (3) Mismatch triple + invariant-ratio collapse | Section III (lines 181-203) | match |
| (4) Projector split O_orb + Q_quot, support on dependent triple | Section IV (lines 244-287) | match |
| (5) Cocycle laws Phi^31, m^31, q^31 | Section VI (lines 350-388) | match |
| (6) Two-point orbit-lock theorem (q=0 <=> m=1 <=> invariant ratios=1) | Section IV q-chart + Q_quot support (lines 227-285), Section II same-orbit (lines 163-174) | match (algebraic core fully covered; "in G_*" is the definitional invariant-ratio equality, exercised in Section II) |
| (Stage-198 reduction, in notes §7 and card "Inputs") | Section VII (lines 395-401) | match |
| (restoration map, notes §4.2) | Section V (lines 294-309) | match (extra detail beyond card, fully aligned with notes) |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side check; no script check tests anything the paper/notes do not assert.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 123 | `Ctr_ratio_pairs - Ctr_ratio_expected == 0` | claim 1 (Ctr pairwise form) | yes |
| A2 | sympy | 124 | `Cnt_ratio_pairs - Cnt_ratio_expected == 0` | claim 1 (Cnt pairwise form) | yes |
| A3 | sympy | 125 | `eps_ratio_pairs - eps_ratio_expected == 0` | claim 1 (eps pairwise form) | yes |
| A4 | sympy | 154 | `ln Phi_mu - ln expanded == 0` | claim 2 (Phi_mu factored=expanded) | yes |
| A5 | sympy | 163 | `ln[Ctr ratio on same orbit] == 0` | claim 6 fwd / claim 2 (Phi_T solves same-orbit) | yes |
| A6 | sympy | 167 | `ln[eps ratio on same orbit] == 0` | claim 2 (Phi_K solves same-orbit) | yes |
| A7 | sympy | 171 | `ln[Cnt ratio on same orbit] == 0` | claim 2 (Phi_mu solves same-orbit) | yes |
| A8 | sympy | 192 | `ln Ctr_ratio - (1+chi0)ln m_T == 0` | claim 3 (collapse) | yes |
| A9 | sympy | 197 | `ln eps_ratio - ln(1/m_K) == 0` | claim 3 (collapse) | yes |
| A10 | sympy | 200 | `ln Cnt_ratio - (ln m_mu - ln m_K - F ln m_T) == 0` | claim 3 (collapse) | yes |
| A11 | sympy | 227 | `q_tr - (1+chi0)ln m_T == 0` | claim 6 (chart) | yes |
| A12 | sympy | 231 | `q_eta + ln m_K == 0` | claim 6 (chart) | yes |
| A13 | sympy | 235 | `q_nt - (ln m_mu - ln m_K - F ln m_T) == 0` | claim 6 (chart) | yes |
| A14 | sympy | 283 | `Q_quot Delta x - expected == 0` | claim 4 (quotient support = m-triple) | yes |
| A15 | sympy | 284 | `O_orb Delta x - expected == 0` | claim 4 (orbit part = free + Phi) | yes |
| A16 | sympy | 285 | `Q + O - Delta x == 0` | claim 4 (split is a partition) | yes |
| A17 | sympy | 286 | `M_* O_orb Delta x == 0` | claim 4/6 (O_orb in kernel of M_*) | yes |
| A18 | sympy | 287 | `M_* Q_quot Delta x - qpair == 0` | claim 4 (quotient carries full q) | yes |
| A19 | sympy | 298-308 | `ln T2/K2/mu2_restore - ln(Phi.x1) == 0` | restoration map (notes §4.2) | yes |
| A20 | sympy | 350-352 | `Phi^31 - Phi^32 Phi^21 == 0` (T,K,mu) | claim 5 (transport cocycle) | yes |
| A21 | sympy | 364-374 | `ln m^31 - ln m^32 - ln m^21 == 0` (T,K,mu) | claim 5 (mismatch cocycle) | yes |
| A22 | sympy | 387-388 | `Delta x^31 - ... == 0`, `q^31 - q^32 - q^21 == 0` | claim 5 (additive cocycle) | yes |
| A23 | sympy | 396-401 | `Phi=1`, `m=r` at equal free ratios | Stage-198 reduction | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

No `.wl` script exists for this unit, so the transliteration check does not apply. See the engine cross-check and the single-engine justification below.

## Engine cross-check

Only the SymPy engine is present. No cross-engine comparison is possible; `engines_agree: n/a`.

### Single-engine sufficiency (line-114 judgment, not a reflexive line-118 demand)

This unit is `is_status_only_candidate: False` and checkpoint `False`, but it is SymPy-only by design (card `\stagefield{Verification}`: "Mathematica audit: none yet"). Per the prompt's own line-114 guidance and established pipeline precedent (stages 121/122/123 verified SymPy-only; ~56 SymPy-only stages pipeline-wide), single-engine is acceptable when SymPy genuinely settles the claimed result. Here every assertion is an exact symbolic identity in a commutative monomial/log-linear algebra over strictly positive real symbols:
- The invariant ratios are products of integer/symbolic powers of positive ratios; equality reduces to comparing exponent vectors. SymPy's `expand_log(..., force=True)` on positive bases is exact (no branch ambiguity), and `simplify` on monomials is decisive.
- The transport-factor checks substitute an independently-defined `Phi` back into a ratio built from the primitive monomials and confirm the result is unity — a genuine "is `Phi` the same-orbit solution?" test, not a re-statement.
- The projector algebra is exact linear algebra: `M_*`'s dependent 3x3 block has determinant `1+chi0_*` (nonzero since chi0_* > 0), so `Pdep_inv` exists and the kernel/range split is well-posed; `M_* O_orb Delta x = 0` is a real test of that split.
There is no transcendental closed form, no numeric tolerance, no integral, and no branch-sensitive simplification that a second CAS would adjudicate differently. A Mathematica port would re-run the identical commutative-algebra reduction and add no independent confidence. **`missing_mathematica` is NOT flagged**: I cannot point to any claimed result that SymPy fails to genuinely verify, which is the prompt's stated bar for the finding.

## Verdict justification

Clean. Attacks tried that failed: (1) Tautology hunt — A5/A6/A7 are not "value substituted into its own defining equation": `Phi_T` (line 133) is defined from `r_U,(r_U/(r_g r_c))^alpha`, then substituted into `Ctr_ratio_expected` (independently built in Section I from the primitive monomials) to confirm it gives 1; A11/A12/A13 cross-link the linear-map output `qpair` (line 222) against the transport-factor-derived mismatches `m_T,m_K,m_mu` (lines 181-183), two independently constructed objects. (2) Hardcoded-result hunt — the "expected" pairwise forms in Section I are themselves cross-checked against the substitution of the primitive monomials (A1-A3), so no literal answer is asserted against itself; the `M_*` entries match the monomial exponents I re-derived by hand (row 1 r_U coeff `-(2+chi0+deltaU)` = `-(1+deltaU)-(1+chi0)`; row 2 r_U coeff `Fstar-Estar`; row 3 `[0,2,0,-1,-1,0,0,0]` = `ln(r_c^2/(r_U r_K))`). (3) Symbol-domain hunt — all symbols `positive=True, real=True`, justified by the notes' "two positive microscopic states"; `expand_log(force=True)` is exact on positive bases. (4) Missing-branch hunt — the cocycle section (VI) uses fully independent 21- and 32-ratio symbols and Phi_mu_of genuinely carries the r_lambda dependence, so the homomorphism property is exercised over the real parameter space, not a special case; the equal-free-ratio reduction (VII) is a legitimate limiting check that does not stand in for the general claim (the general claim is carried by Sections I-VI). I confirm I read the card, the notes, and the part-05 appendix rows: the script's verified claims match the paper's `Output` and all six notes deliverables.

Non-blocking observation (not a finding in any of the ten categories; recorded for transparency): the script's print strings carry stale stage labels — banner "STAGE 182" (lines 40, 403), "Stage 192 projectors" (line 208), "Stage 198"/"Stage 182 LEDGER" — and the notes file internally calls this "Stage 250" (with upstream "Stage 249/248/243"), whereas the card/appendix use 199/198/192-equivalent numbering. These are cosmetic renumbering artifacts in comments/print output only; they touch no assertion and do not alter any verified identity, so they neither block nor warrant a directive (and the script-label text is not a paper-side claim, so this is not `paper_misalignment` of a mathematical result).

## Self-test notes

Variable-independence trap: N/A — the script contains no `sp.diff`; every check is an algebraic/log-linear identity, so no derivative can be vacuously zero. Symmetry/parity trap: N/A — no integrals over any domain. Trivial-case pre-check: I hand-verified the load-bearing identities reduce as claimed (Phi_T solves the same-orbit Ctr=1 condition; Ctr_ratio = m_T^{1+chi0}; M_* dependent-block determinant = 1+chi0 != 0 so the projector split is well-posed; equal-free-ratio limit gives Phi=1). All passed, consistent with the committed output's `EXIT_CODE: 0`, which is fresh (output mtime 2026-05-11 12:48 > script mtime 2026-05-11 11:58).
