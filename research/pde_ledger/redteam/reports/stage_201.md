---
unit_id: 201
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection.md]
  paper_appendix: present
---

# Audit unit 201 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_201.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows: line 33 status row; lines 129–237 theorem chain items 1–2 and eq:app-part06-main-graph-quotient; lines 239–313 intrinsic realization packet narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` reads verbatim: "Intrinsic realization packet, canonical orbit projection \(\Pi_{\mathcal O_*}^{\rm can}\), fixed-point criterion \(\mathbf x\in\mathcal Z_*\iff \chi_Q=1\) and \(\Pi_{\mathcal O_*}^{\rm can}(\mathbf x)=\mathbf x\)." Stage 201 turns the reference-free four-scalar home-stretch theorem into an executable realization audit: it defines intrinsic target monomial ratios \((\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta)\) and their logs \((q_{\rm tr},q_{\rm nt},q_\eta)\); forms the four-scalar packet \(\Delta_{\rm real}^{\rm int}=(\chi_Q-1,q_{\rm tr},q_{\rm nt},q_\eta)\); solves the dependent-triple repair equation \(M_*\,\Delta\mathbf x_{\rm rep}=-\mathbf q\) using the canonical section \(S_{(T,K_\eta,\mu)}\); exponentiates the repair to the canonical orbit projection \(\Pi^{\rm can}_{\mathcal O_*}\); and proves the projected state is the unique target-orbit point with the same free quintuple \((\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\rm eff)})\). Distinct deliverables (from card + notes §2–§7 + appendix theorem items 1–2 and eq:app-part06-main-graph-quotient): (D1) the intrinsic packet form; (D2) the mismatch chart \(m_T=\mathfrak R_{\rm tr}^{1/(1+\chi_{0,*})}\), \(m_K=\mathfrak R_\eta^{-1}\), \(m_\mu=\mathfrak R_{\rm nt}\mathfrak R_\eta^{-1}\mathfrak R_{\rm tr}^{F_*/(1+\chi_{0,*})}\); (D3) the exact repair vector supported only on \((T,K_\eta^{(\rm eff)},\mu_W)\); (D4) the cancellation \(M_*\,\Delta\mathbf x_{\rm rep}=-\mathbf q\) via \(M_*S=I_3\); (D5) the canonical orbit projection \(\Pi^{\rm can}_{\mathcal O_*}\); (D6) same-free-quintuple uniqueness; (D7) the fixed-point criterion and equivalence to \(\Delta_{\rm real}^{\rm int}=0\); (D8) the first-order linearized compiler \(M_*(\Delta\mathbf x+\Delta\mathbf x_{\rm rep}^{\rm lin})=0\).

## What the script claims to verify

The SymPy script builds the explicit \(3\times8\) Stage-192 quotient map \(M_*\) (parameterized symbolically by \(\chi_{0,*},\delta U_*,E_*,F_*\)), extracts the dependent pivot block \(P_{(T,K_\eta,\mu)}\) from columns \(\{T,K_\eta,\mu\}\), inverts it, builds the canonical section \(S_{(T,K_\eta,\mu)}=E_{\rm dep}P^{-1}\), and asserts \(M_*S-I_3=0\) (D4). It then asserts the mismatch-chart identities \(m_T,m_K,m_\mu\) (D2), the repair vector \(\Delta\mathbf x_{\rm rep}=-S\mathbf q\) equals the closed form on the dependent triple (D3), the cancellation \(M_*\Delta\mathbf x_{\rm rep}+\mathbf q=0\), an independent `sp.solve` of \(M_*\,\Delta x_{\rm dep}=-\mathbf q\) over \((dT,dK_\eta,d\mu)\) reproducing the same triple (D6 uniqueness), the orbit projection's log-ratio equals \(\Delta\mathbf x_{\rm rep}\) (D5), the pairwise-witness packet equals the intrinsic packet under \(\mathbf x_1\to\mathbf x_*\), the fixed-point reduction \(\Pi^{\rm can}(x)=x\) at \(\mathfrak R=1\) (D7), and the linearized \(M_*(\Delta x+\Delta x_{\rm rep}^{\rm lin})=0\) (D8). All checks pass (exit 0).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 intrinsic packet \((\chi_Q-1,q_{\rm tr},q_{\rm nt},q_\eta)\) | §II `Delta_real_int`; §VI pairwise=intrinsic | match |
| D2 mismatch chart \(m_T,m_K,m_\mu\) | §II `expect_zero` on `mT`,`mK`,`mMu` (lines 116–121) | match |
| D3 repair vector on \((T,K_\eta,\mu)\) | §III `repair vector - expected vector` (line 147) | match |
| D4 cancellation via \(M_*S=I_3\) | §I `M_* S - I_3` (line 85); §III `M_* Δx_rep + q` (line 148) | match |
| D5 orbit projection \(\Pi^{\rm can}_{\mathcal O_*}\) | §V `log(x_proj/x) - Δx_rep` (line 219) | match |
| D6 same-free-quintuple uniqueness | §IV independent `sp.solve` + three `expect_zero` (lines 162–186) | match |
| D7 fixed-point criterion \(\Pi^{\rm can}(x)=x\iff\mathfrak R=1\) | §VII line 246 (forward direction at \(\mathfrak R=1\)) | match (see note) |
| D8 linearized compiler | §VII `M_*(Δx+Δx_rep^lin)` (line 261) | match |

`paper_alignment: aligned`. Every paper-side deliverable maps to a non-tautological script check that exercises it; no `mismatch`, no `extra`. (Note on D7: §VII verifies the forward direction — \(\mathfrak R_{\rm tr}=\mathfrak R_{\rm nt}=\mathfrak R_\eta=1\Rightarrow\Pi^{\rm can}(x)=x\). The converse \(\Pi^{\rm can}(x)=x\Rightarrow\mathfrak R=1\) follows trivially from §V because \(\log(x_{\rm proj}/x)=\Delta\mathbf x_{\rm rep}\) and the repair vector's nonzero entries vanish only when all three \(q\)'s are zero; this is structurally established, so I do not raise it as an `insufficient_verification` finding. The Mathematica claim manifest below makes the iff explicit so the second engine closes the gap.)

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 85 | `expect_zero(M_* S - I_3)` | D4 | yes |
| A2 | sympy | 113–115 | `exp(q_*) - R_* == 0` | D1 bookkeeping | partial (definitional; exp∘log) |
| A3 | sympy | 116 | `m_T - R_tr^(1/(1+chi0))` | D2 | yes |
| A4 | sympy | 117 | `m_K - R_eta^(-1)` | D2 | yes |
| A5 | sympy | 118–121 | `m_mu - R_nt R_eta^-1 R_tr^(F/(1+chi0))` | D2 | yes |
| A6 | sympy | 147 | `repair vector - expected vector` | D3 | yes |
| A7 | sympy | 148 | `M_* Δx_rep + q == 0` | D4 | yes |
| A8 | sympy | 175–186 | `unique dT,dKeta,dmu` from `sp.solve` | D6 | yes |
| A9 | sympy | 219 | `log(x_proj/x) - Δx_rep == 0` | D5 | yes |
| A10 | sympy | 236–239 | `pairwise witness packet - intrinsic packet` | D1 | yes |
| A11 | sympy | 246 | `Pi^can(x)=x at R=1` | D7 (forward) | yes |
| A12 | sympy | 261 | `M_*(Δx + Δx_rep^lin) == 0` | D8 | yes |

A2 (lines 113–115) is definitional: with `qtr=sp.log(Rtr)`, `exp(qtr)-Rtr` is algebraically guaranteed zero. It is harmless bookkeeping that sits beside the substantive A3–A12 checks, so I do not raise a standalone `tautological_check` finding — the load-bearing verification (D2–D8) is carried by genuine matrix algebra (\(M_*S=I_3\), the inverted pivot block, the independent `sp.solve`).

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `(missing)` — target `mathematica/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_audit.wl`

**What's wrong:**
Stage 201 is non-status-only (`is_status_only_candidate: False`) and not a checkpoint, yet only the SymPy engine exists. The paper card line 11 itself records: "Mathematica audit: none yet." The stage's content — building the explicit \(3\times8\) integer/symbol matrix \(M_*\), inverting the \(3\times3\) dependent pivot block, verifying \(M_*S=I_3\), the closed-form repair vector, an independent linear solve, and the log/exp orbit-projection identities — is exactly the kind of finite symbolic linear algebra Mathematica can verify independently (matrix `Inverse`/`LinearSolve`, `Solve`, `FullSimplify`, `Series`/`PowerExpand` for the log-ratio identities). There is no impossibility carve-out: the second engine is genuinely possible here, so its absence is a finding under the dual-engine contract (the only legitimate omission is a pure status/label carry-forward, which this stage is not).

**Why this matters:**
The repair vector, mismatch chart, and orbit projection are the load-bearing forward results cited downstream (the card's "Downstream use" routes \resultanchor{MTDC-T10.1} to Part VI). A single-engine verification cannot catch a SymPy-specific simplification artifact in the \(P^{-1}\) inversion or the `expand_log(..., force=True)` steps used throughout `expect_zero`. A second engine deriving the same results by a different decomposition is the only guard against an engine-local algebra error.

**Required change:**
Add an independent Mathematica `.wl` audit per the directive's claim manifest (M1–M8). It must derive the results natively, not transliterate the SymPy choreography.

**Verification:**
`redteam exec-mathematica 201` produces a saved output, the script exits 0, and each manifest claim M1–M8 appears as an explicit `expectZero`/`If[...Exit[1]]` check whose residual is reported zero.

## Independent-derivation check (Mathematica)

No `.wl` exists — not applicable. The directive requires the new `.wl` to use a different decomposition than the SymPy script (e.g., `LinearSolve` for the dependent triple rather than mirroring `Edep . Inverse[Pdep]`, and `PowerExpand`/`Simplify` of explicit exponent arithmetic rather than the `expand_log(force=True)` path), with an anti-transliteration guard.

## Engine cross-check

Only one engine present; `engines_agree: n/a`. The directive's manifest is the cross-check target.

## Verdict justification

Verdict `findings`, one finding (`missing_mathematica`). The SymPy script is substantive and paper-aligned: I attacked the pivot-block inversion (det \(P_{(T,K_\eta,\mu)}=1+\chi_{0,*}\), inverse matches output lines 35–44, genuine computation not a literal), checked the matrix `M_*` column conventions against the appendix repair vector eq:app-part06-main-repair-vector (lines 146–153) and the mismatch chart eq:app-part06-main-graph-quotient (lines 174–181) — all signs and the \(F_*/(1+\chi_{0,*})\) factor match; the §IV uniqueness solve is genuinely independent of the §III closed form (it re-derives the triple via `sp.solve` rather than reusing \(S\)); the §VII linearized check is non-trivial because `q_from_Dx = M_* Dx` actually depends on all eight drift components. I confirmed I read the paper card, the notes, and the appendix theorem chain, and the script's verified claim matches the paper's stated claim. The only defect is the absent second engine, which the dual-engine contract requires.

## Self-test notes

Variable-independence trap: the §VII linearized residual `M_*(Dx + dx_rep_lin)` is built from `q_from_Dx = Mstar*Dx` where `Dx` carries all eight independent drift symbols, so the result is not identically zero by construction — it vanishes only because `S` is a genuine right-section (`M_*S=I_3`); the prescribed Mathematica M8 will reproduce this non-trivially. Symmetry/parity: no unbounded-domain integrals in this stage, so the parity trap does not apply. Trivial-case pre-check: substituting \(\mathfrak R_{\rm tr}=\mathfrak R_{\rm nt}=\mathfrak R_\eta=1\) (all \(q=0\)) into the repair vector and into \(\Pi^{\rm can}\) gives the zero vector / identity map respectively, matching §VII A11, and substituting a single \(\mathfrak R_\eta\neq1\) gives a nonzero \(K_\eta\)-entry only — the support claim D3 holds. Path spec: target `.wl` placed in `mathematica/`. Paper round-trip: the manifest constants \(\chi_{0,*},F_*\) and the exponent \(1/(1+\chi_{0,*})\), \(F_*/(1+\chi_{0,*})\) match eq:app-part06-main-graph-quotient and notes §4 exactly; no new constant introduced.
