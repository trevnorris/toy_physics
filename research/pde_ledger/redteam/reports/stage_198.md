---
unit_id: 198
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
  notes_stage_files: [moving_throat_pde_stage198_exact_finite_orbit_law.md]
  paper_appendix: present
---

# Audit unit 198 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_198.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage198_exact_finite_orbit_law.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 127, 1329–1361, 1468)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` reads verbatim: "Solves the dependent triple exactly, defines mismatch coordinates, and proves finite single-orbit lock." The notes file (titled "Stage 249" internally, the original stage number; the audit-card number is 198) enumerates five deliverables: (1) the exact finite single-orbit law giving the dependent triple `(T_U, K_eta^(eff), mu_W)` as closed forms of the five free coordinates and the invariant triple `(C_tr,*, C_nt,*, eps_eta,*)` — formulas in notes §2.1–2.3; (2) the exact dependent residual mismatch triple `(m_T, m_K, m_mu)` and the invariant ratios `R_tr = m_T^(1+chi0*)`, `R_eta = 1/m_K`, `R_nt = m_mu/(m_K m_T^(F*))` — §3; (3) the exact logarithmic chart `q_tr=(1+chi0*)ln m_T`, `q_eta=-ln m_K`, `q_nt=ln m_mu - ln m_K - F* ln m_T`, together with the statement that applying the upstream drift compiler `M_*` to a pure dependent-mismatch vector reproduces this chart exactly — §4 and §4.1; (4) the exact algebraic restoration map that returns a candidate branch to the orbit by changing only the dependent triple — §5; (5) the sharp orbit-lock theorem `Delta_orbit=0 <=> m_T=m_K=m_mu=1 <=> q_tr=q_nt=q_eta=0` — §6. The appendix (lines 1339–1361) restates the same chart and lock condition and assigns the result to anchor `MTDC-T9.6`.

## What the script claims to verify

The SymPy script defines the three exact coherent monomials `Ctr, Cnt, epsEta` symbolically, writes down the proposed closed-form orbit values `Keta_orbit, TU_orbit, muW_orbit`, and then verifies by substitution that each monomial evaluated at the orbit point reproduces the corresponding invariant constant (`epsEta(Keta_orbit)/epsEta_* = 1`, `Ctr(TU_orbit)/Ctr_* = 1`, `Cnt(...)/Cnt_* = 1`). It then introduces actual values `m_T·orbit`, `m_K·orbit`, `m_mu·orbit` and verifies the invariant ratios collapse exactly to `m_T^(1+chi0*)`, `1/m_K`, `m_mu/(m_K m_T^(F*))`. It defines the log chart, verifies the hardcoded drift matrix `M_*` times the pure dependent-mismatch vector `(0,0,0,0,kappa,0,mu,tau)` equals `(q_tr,q_nt,q_eta)`, verifies the restoration map returns the orbit values, and verifies the inverse chart is consistent with the forward chart and that `q=0 => m=1`. All checks are `expect_zero` (raise on nonzero); the saved output shows every residual equal to 0 and exit code 0.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Exact orbit triple formulas §2.1–2.3 | lines 79–88 define the same forms; lines 97–108 verify each monomial at the orbit point returns its invariant | match |
| (2) Mismatch ratios / invariant ratios §3 | lines 121–136 verify `Ctr_ratio=m_T^(1+chi0*)`, `epsEta_ratio=1/m_K`, `Cnt_ratio=m_mu/(m_K m_T^F*)` | match |
| (3a) Log chart §4 | lines 144–146 define exactly the paper's chart | match |
| (3b) `M_* Δx_mis = (q_tr,q_nt,q_eta)` §4.1 | lines 156–170 verify the matrix-vector product for the pure dependent-mismatch vector | match |
| (4) Restoration map §5 | lines 177–196 verify `T/K/mu_restore = orbit` using the paper's exponents | match |
| (5) Orbit-lock theorem §6 | lines 204–238 verify forward/inverse chart consistency and `q=0 => m=1` | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 97–100 | `epsEta(Keta_orbit)/epsEta_* - 1 == 0` | claim 1 (K_eta orbit) | yes |
| A2 | sympy | 101–104 | `Ctr(TU_orbit)/Ctr_* - 1 == 0` | claim 1 (T_U orbit) | yes |
| A3 | sympy | 105–108 | `Cnt(muW_orbit,Keta_orbit,TU_orbit)/Cnt_* - 1 == 0` | claim 1 (mu_W orbit) | yes |
| A4 | sympy | 134 | `Ctr_ratio - m_T^(1+chi0*) == 0` | claim 2 (R_tr) | yes |
| A5 | sympy | 135 | `epsEta_ratio - 1/m_K == 0` | claim 2 (R_eta) | yes |
| A6 | sympy | 136 | `Cnt_ratio - m_mu/(m_K m_T^F*) == 0` | claim 2 (R_nt) | yes |
| A7 | sympy | 167–170 | `M_* Δx_mis - (q_tr,q_nt,q_eta) == 0` | claim 3b (drift map) | yes (cols 4,6,7) |
| A8 | sympy | 194 | `T_restore - T_orbit == 0` | claim 4 | yes |
| A9 | sympy | 195 | `K_restore - K_orbit == 0` | claim 4 | yes |
| A10 | sympy | 196 | `mu_restore - mu_orbit == 0` | claim 4 | yes |
| A11 | sympy | 214–225 | three forward/inverse chart consistency checks | claim 3a / 5 | yes |
| A12 | sympy | 236–238 | `m_T(q),m_K(q),m_mu(q)` at `q=0` minus 1 == 0 | claim 5 (lock) | yes |

All rows trace to a paper deliverable; no orphaned assertions. No row is tautological in the bad sense (each substitutes a proposed closed-form solution into an independently-defined expression and confirms a nontrivial collapse, or multiplies an independently-built matrix against a vector).

## Findings

None. (See "Observation" below for a non-blocking cosmetic note.)

### Observation (non-finding) — stage-number labels in print strings

The script's banner (line 35), ledger header (line 240), and ledger prose (lines 257, 264) print "STAGE 181", "Stage 192 drift compiler", and "Stage 197 Packet-A finish-line theorem". The audit card is numbered 198 and the notes use the original numbers 249 / 243 / 248. These are decorative `print()` labels only; they touch no `expect_zero` assertion and do not affect any verified identity. This is a labeling inconsistency, not a math defect and not a `paper_misalignment` of a verified claim. I deliberately do NOT raise it as a formal finding or write a directive: the only change would be editing print-string literals, which the directive template explicitly forbids as a stylistic change, and the pipeline is known to carry multiple parallel renumbering conventions. Recorded here for traceability.

## Independent-derivation check (Mathematica)

No `.wl` exists, so `mathematica_transliteration` does not apply.

## Engine cross-check

Only one engine present; `engines_agree: n/a`.

## Verdict justification

`clean`. I read the card, the notes, and the appendix rows first and built the model of all five deliverables, then attacked the script. Attacks tried and failed: (1) Tautology — the orbit-solve checks (A1–A3) are NOT `X==X`; each substitutes the author's proposed closed form into an independently-defined monomial and confirms it collapses to the invariant constant, which would fail under any algebra slip (e.g., a missing exponent or factor). The critical simplification `((X)^(1/(1+chi0*)))^(1+chi0*) = X` is valid under the declared positivity and is confirmed collapsing to 0 in the saved output. (2) Homogeneity (A4–A6) genuinely extracts the monomial exponents `(1+chi0*)`, `-1`, `(1,-1,-F*)`; a wrong exponent in the chart would surface here. (3) The drift-map check (A7) exercises the load-bearing columns 4/6/7 of `M_*` against the pure dependent-mismatch vector exactly as the paper's §4.1 claim requires; the free-coordinate columns are Stage 243's responsibility, not this stage's, so leaving them unexercised is correct scoping, not `insufficient_verification`. (4) The restoration (A8–A10) and inverse-chart (A11–A12) checks confirm the exp-exponents cancel and the linear chart inverts consistently, supporting the lock theorem in both directions via the verified bijection. Symbol assumptions (`positive=True` on all coordinates and on `chi0*, deltaU*, E*, F*`) match the notes' "positive microscopic state" and are, if anything, conservative; they do not make any identity pass for the wrong reason since the ratio identities hold for all real exponents. On the single-engine question: every claim is exact commuting-real symbolic algebra with no numerical approximation, no branch ambiguity, and only one elementary power-law simplification that SymPy handles deterministically — so SymPy fully and genuinely settles the stage and a second engine would merely re-derive the same trivial algebra. I therefore do not flag `missing_mathematica`; I cannot name any claimed result SymPy fails to verify. The output is fresh (output mtime 12:48:54 > script mtime 11:58:53, same day), so no `stale_output`.

## Self-test notes

Variable-independence trap: there is no `sp.diff` in this script, so no zero-derivative-against-unwired-symbol risk; the one matrix product (A7) was hand-checked to pick out columns 4/6/7, which carry the correct coefficients `(0,0,1+chi0*)`, `(-1,1,-F*)`, `(-1,0,0)`. Parity/integral trap: not applicable (no integrals). Trivial-case pre-check: I hand-substituted each orbit form back into its monomial and confirmed the residual reduces to 0 for non-degenerate symbols, and confirmed the restoration exp-exponents algebraically cancel to `T_orbit/K_orbit/mu_orbit`; the `q=0 => m=1` checks reduce to `exp(0)-1=0`, consistent with the verified bijective chart. No directive is written (zero script-side findings).
