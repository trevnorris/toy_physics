---
unit_id: 178
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage178_outgoing_port_coloading_theorem.md]
  paper_appendix: present
---

# Audit unit 178 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_178.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage178_outgoing_port_coloading_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (read intro, source-status row for 178 at line 87, and the dedicated stage-178 block at lines 654–669)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.txt`

## What the paper claims

The stage card's `\stagefield{Output}` is: "Rewrites \(\Xi_1\) as the weighted mismatch between actual static outgoing-port loading and wall-baseline loading." The appendix block (lines 654–669, eqs `app-part05-port-N-static`, `app-part05-nu-r-def`, `app-part05-Xi1-port-coloading`) makes this concrete: for each outgoing port the static transfer coefficient is \(N_{A,0}^{(r)}=P_{A,r}^2/\Delta_{A,r}^2\), its weak-axisymmetric log-slope is \(\delta\ln N_{A,0}^{(r)}=\epsilon\lambda_A\nu_r\), and the remaining grouped defect is \(\Xi_1=\bar\nu_N-\kappa_1\) (outgoing-weighted port slope minus wall-baseline slope). The notes add the full derivation: (i) \(\nu_r:=2(\mathfrak p_r-\mathfrak d_r)\); (ii) the weighted collapse \(\Xi_1=\sum_r\rho_r^{(N)}(\nu_r-\kappa_1)=\bar\nu_N-\kappa_1\) using \(\sum_r\rho_r^{(N)}=1\); (iii) closed forms \(\mathfrak p_r=\alpha_r(\mathfrak o_U+\mathfrak g_W)+\beta_r(\mathfrak r_r+\mathfrak g_U)\) and \(\mathfrak d_r=\chi_r(\mathfrak o_U+\mathfrak o_W)-2\zeta_r\mathfrak r_r\) with \(\alpha_r+\beta_r=1,\ \chi_r-\zeta_r=1\); and (iv) the exact equivalence with the Stage-177 slippage language, \(\nu_r=\kappa_1+2\mathfrak m_r+\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\), giving \(\sigma_r=\nu_r-\kappa_1\). The corollaries in notes §5 (zero-defect \(\bar\nu_N=\kappa_1\), per-port sufficient condition, dominant-port limit, naive-rigidity \(\Xi_1=-\kappa_1\)) are immediate substitutions into those identities, not independent deliverables.

## What the script claims to verify

The SymPy docstring lists four checks that map one-to-one onto the four notes deliverables: (1) the static outgoing-transfer slope \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) from series-expanding \(P_{A,r}^2/\Delta_{A,r}^2\); (2) the weighted collapse \(\Xi_1=\bar\nu_N-\kappa_1\) under \(\sum\rho=1\); (3) the closed port-variable forms for \(\mathfrak p_r,\mathfrak d_r\) plus the weight identities \(\alpha+\beta=1\) and \(\chi-\zeta=1\); (4) the equivalence \(2\mathfrak p_r-2\mathfrak d_r=\kappa_1+2\mathfrak m_r+\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\), i.e. \(\sigma_r=\nu_r-\kappa_1\). The Mathematica `.wl` performs the identical four checks. All checks are residual-equals-zero (`expect_zero` / `expectZero`), so the bottom line is genuinely "these residuals vanish symbolically over free slope symbols."

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\), slope of \(N_{A,0}^{(r)}=P^2/\Delta^2\) | sympy §1 / math §1 | match |
| \(\Xi_1=\bar\nu_N-\kappa_1\) via \(\sum\rho_r^{(N)}=1\) | sympy §2 / math §2 | match |
| \(\mathfrak p_r,\mathfrak d_r\) port-variable forms; \(\alpha+\beta=1\), \(\chi-\zeta=1\) | sympy §3 / math §3 | match |
| equivalence with Stage-177 slippage \(\sigma_r=\nu_r-\kappa_1\) | sympy §4 / math §4 | match |
| notes §5 corollaries (zero-defect, per-port, dominant-port, naive-rigidity) | — | not separately checked, but each is an algebraic substitution into §1+§2 identities (not an independent deliverable) |

Dominant pattern: every load-bearing paper deliverable is faithfully exercised. `paper_alignment: aligned`. No `extra` checks (nothing the paper does not mention). No numeric constant mismatch — the script carries the generic lane factor `eps*lam` (= \(\epsilon\lambda_A\)) exactly as the paper defines \(\delta\ln N_{A,0}^{(r)}=\epsilon\lambda_A\nu_r\); the signature values \(\lambda_{20}=1,\lambda_{21}=1/2,\lambda_{22}=-1\) are not needed because every identity is lane-independent (the \(\epsilon\lambda_A\) factor cancels).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_zero(nu_r - 2(p_r-d_r))` | \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) | yes |
| A2 | sympy | 69–72 | `expect_zero(Xi_1 - (nubar_N - kappa_1))` (rho3→1−rho1−rho2) | \(\Xi_1=\bar\nu_N-\kappa_1\) | yes |
| A3 | sympy | 103 | `expect_zero(p_r formula)` | \(\mathfrak p_r\) form | yes |
| A4 | sympy | 104 | `expect_zero(d_r formula)` | \(\mathfrak d_r\) form | yes |
| A5 | sympy | 105 | `expect_zero(alpha+beta-1)` | weight identity | yes |
| A6 | sympy | 106 | `expect_zero(chi-zeta-1)` | weight identity | yes |
| A7 | sympy | 135 | `expect_zero(nu_r - [kappa1 + sigma_r])` | Stage-177 equivalence | yes |
| A8 | math | 37 | `expectZero[nu_r - 2(p_r-d_r)]` | \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) | yes |
| A9 | math | 47–50 | `expectZero[Xi_1 - (nubar_N - kappa_1)]` | \(\Xi_1=\bar\nu_N-\kappa_1\) | yes |
| A10 | math | 76–79 | `expectZero[p_r/d_r/alpha+beta-1/chi-zeta-1]` | \(\mathfrak p_r,\mathfrak d_r\), weights | yes |
| A11 | math | 101 | `expectZero[nu_r - [kappa1 + sigma_r]]` | Stage-177 equivalence | yes |

All assertions are non-tautological: each compares a quantity built one way (series expansion of the actual squared port ratio, or the slippage combination) against an independently-constructed target. I worked the §4 algebra by hand and confirmed \(2\mathfrak p_r-2\mathfrak d_r\) collapses to \(2\mathfrak g_W-2\mathfrak o_W+2\beta\mathfrak i_r+2\zeta\mathfrak h_r\) using \(\alpha+\beta=1,\ \chi=1+\zeta,\ 2\mathcal I_r/(1+\mathcal I_r)=2\beta,\ 2\mathcal H_r/(1-\mathcal H_r)=2\zeta\), matching `nu_expected` term-by-term. The identity is correct and genuinely exercised.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl:29–101`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py:42–135`

**What's wrong:**
The `.wl` is a line-by-line transliteration of the `.py`, not an independent re-derivation. Each section mirrors the SymPy choreography with the same variable roles, the same intermediate quantities, and the same target expressions rewritten in Mathematica syntax:

- §1 SymPy `NAr = sp.expand(PAr**2 / DAr**2)` … `nu_from_series = sp.simplify((NAr_lin / N0r - 1) / (eps * lam))` ↔ `.wl` `nAr = Expand[pAr^2/dAr^2]` … `nuFromSeries = FullSimplify[(Normal[Series[nAr,{eps,0,1}]]/n0r - 1)/(eps*lam), ...]` — identical construction and the identical comparison `… - 2*(pR - dR)`.
- §3 SymPy `alpha = OU2*GW/P; beta = R*GU/P; chi = OU2*OW2/Delta; zeta = R**2/Delta` with `p_expected = alpha*(oU+gW)+beta*(rr+gU)`, `d_expected = chi*(oU+oW)-2*zeta*rr` ↔ `.wl` `alpha = ou2*gw/p; beta = r*gu/p; chi = ou2*ow2/delta; zeta = r^2/delta` with `pExpected = alpha*(oU+gW)+beta*(rr+gU)`, `dExpected = chi*(oU+oW)-2*zeta*rr` — same weights, same combination, same order.
- §4 SymPy `nu_expected = kappa1 + 2*m_r + (2*Ir_expr/(1+Ir_expr))*i_r + (2*Hr_expr/(1-Hr_expr))*h_r` ↔ `.wl` `nuExpected = kappa1 + 2*mR + 2*iRExpr*iR/(1+iRExpr) + 2*hRExpr*hR/(1-hRExpr)` — the same slippage formula assembled in the same way.

Per category 6 this violates the second-engine policy: both engines must derive the result independently from the physical premises rather than echo each other's algebra. The algebra here is sound (verified by hand), but the Mathematica run provides essentially no independent confirmation — a construction error in the shared choreography would be reproduced identically, not caught.

**Why this matters:**
The second engine exists to catch a derivation mistake the first engine might make. A verbatim port reduces that to a CAS-kernel comparison only; the structural-independence guarantee is lost.

**Required change:**
Re-derive the central §4 equivalence in the `.wl` through a route that does not reuse the SymPy variable choreography, then cross-check it against the §1–§3 outputs. Concretely, build \(\nu_r\) from the port data the *other* direction: form the static transfer coefficient \(N_{A,0}^{(r)}=P^2/\Delta^2\) symbolically with the component log-drifts substituted directly (`Ω_U^2→ou2*(1+eps*lam*oU)`, `G_W→gw*(1+eps*lam*gW)`, `R→r*(1+eps*lam*rr)`, `G_U→gu*(1+eps*lam*gU)`, `Ω_W^2→ow2*(1+eps*lam*oW)`), extract `nuFromData = Coefficient[Normal[Series[Log[nAr], {eps,0,1}]], eps*lam]` (or equivalent), and assert `expectZero["nu via log-data vs slippage", nuFromData - nuExpected]`. That gives the `.wl` an independent path from premises to \(\nu_r\) that does not retrace the `.py`'s `p_expected`/`d_expected` factoring, while still confirming the same physical identity. Keep the existing §1–§3 checks. (See directive for the bounded edit. This is a policy-driven change; if an established MATHEMATICA_MIRROR_POLICY explicitly sanctions mirrored second engines for this paper, the orchestrator should down-weight or waive this finding — the auditor cannot see that tracker.)

**Verification:**
After the edit the `.wl` contains a new `expectZero` whose left expression is built from the log of \(N_{A,0}^{(r)}\) (not from `pExpected`/`dExpected`), residual prints `0`, the script exits 0, and the new check line appears in the refreshed transcript.

### F2 — stale_output (label defect, informational)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py:34`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl:26`

**What's wrong:**
Both scripts print `banner("STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM")`, and this wrong stage number propagates into both saved transcripts (`.../scripts/output/...sympy_audit.txt:11` and `.../mathematica/output/...mathematica_audit.txt:11`: "STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM"). This unit is stage 178. The SymPy module docstring header and the Mathematica final line ("Stage 178 Mathematica audit passed.") correctly say 178, so the banner is an internal-inconsistency mislabel. No assertion depends on the banner string, so the math is unaffected; this is a provenance/transcript-attribution defect only. (Filed under `stale_output` as the closest transcript-content category; mtimes themselves are fine — outputs are newer than scripts.)

**Why this matters:**
A transcript whose header misattributes the stage number can mislead provenance audits (CHECKPOINT_TRUST_AUDIT-style) and makes the saved output look like it belongs to a different unit.

**Required change:**
In both scripts change the banner title literal from `STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM` to `STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM` (sympy line 34, math line 26). No other change.

**Verification:**
Refreshed transcripts' first banner reads `STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM`; scripts still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py`, not an independent derivation — see F1 for three corresponding-section quotes. The same variable choreography (`pAr/dAr/n0r/nAr`, `alpha/beta/chi/zeta`, the §4 `kappa1 + 2*mR + …` assembly) appears in both engines with only syntax differences. The underlying algebra is correct (hand-verified), but structural independence is absent.

## Engine cross-check

Both engines report all checks passing with matching residuals (all `0` / `PASS`). Final symbolic forms agree:

- SymPy §1 `nu_r = -2*d_r + 2*p_r`; Mathematica §1 `nu_r = 2*(-dR + pR)` — identical.
- SymPy §3 weights `alpha_r = GW*OU2/(GU*R + GW*OU2)`, `chi_r = OU2*OW2/(OU2*OW2 - R**2)`; Mathematica `alpha_r = (gw*ou2)/(gw*ou2 + gu*r)`, `chi_r = (ou2*ow2)/(ou2*ow2 - r^2)` — identical up to variable-name casing.
- SymPy §4 `nu_r` and `sigma_r` rational forms equal the Mathematica forms (different surface arrangement, same expression; both reduce to the verified `2*p - 2*d` identity since both engines assert the residual is 0).

No `engine_disagreement`.

## Verdict justification

Every load-bearing identity the paper card and notes require for stage 178 is present and non-tautologically exercised in both engines, the constants match (lane-independent `eps*lam` factor, no numeric mismatch), the symbol domains are benign, and I confirmed the most intricate check (§4 slippage-equivalence) by hand. So the math holds up and `paper_alignment` is `aligned`. Two script-side defects remain: the Mathematica engine is a verbatim transliteration of the SymPy engine (F1, medium, policy-driven), and both scripts/transcripts carry a wrong stage-number banner "STAGE 161" (F2, low, provenance label). Neither is a paper_misalignment, neither is stop-cold, and neither propagates to downstream units (the proven identities are unchanged). Verdict: `findings`.

## Self-test notes

(1) Variable independence: the only `diff`-like operations are `series`/`Series` in `eps`; every series target genuinely depends on `eps` (via `eps*lam` in PAr/DAr/P_A/D_A), so no identically-zero-derivative trap. (2) Symmetry/parity: no integrals — N/A. (3) Trivial-case pre-check: I hand-substituted the §4 identity using `alpha+beta=1`, `chi=1+zeta`, `2I/(1+I)=2*beta`, `2H/(1-H)=2*zeta` and confirmed `2p-2d` matches `nu_expected` term-by-term (residual 0), and confirmed §2's `Xi_1-(nubar-kappa1)` vanishes only because `rho3→1-rho1-rho2` is substituted (genuine use of Σρ=1). (4) Path specs: F1 names the `.wl` under `mathematica/`; no missing-script targets. (5) Paper round-trip: the F1 fix adds an independent route to the *same* `nu_r` identity and the F2 fix is a label only — neither introduces a new constant or alters any paper-stated value, so no new paper_misalignment is created.
