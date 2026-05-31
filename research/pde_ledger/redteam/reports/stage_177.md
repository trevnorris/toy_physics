---
unit_id: 177
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T06:42:53Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md]
  paper_appendix: present
---

# Audit unit 177 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_177.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/section read: line 85 status row; §"Outgoing-load scalar" lines 556-715; weak-axisymmetric signature lines 386-393)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.txt`

## What the paper claims

The stage card `\stagefield{Output}` states verbatim: "Shows all weak-axisymmetric outgoing slippages inherit the same grouped signature and collapse to \(\Xi_1=P_1/P_0\)." The notes (the authoritative intent here) enumerate the deliverables explicitly: (D1) the three microscopic outgoing slippages transport with weak-axisymmetric log-slopes \(\mathfrak m_r=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\tfrac12\kappa_1\), \(\mathfrak i_r=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r}\), \(\mathfrak h_r=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r}\); (D2) each port collapses to one amplitude \(\sigma_r=2\mathfrak m_r+\tfrac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\tfrac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\), so \(\Sigma_{A,r}^{(N)}=\epsilon\lambda_A\sigma_r\); (D3) the full grouped defect collapses to one scalar \(\Xi_{\rm load}^{(A)}=\epsilon\lambda_A\Xi_1\) with \(\Xi_1=\sum_r\rho_r^{(N)}\sigma_r\) (the boxed "main structural theorem"); (D4) the fixed grouped signature \(\lambda_{20}=1,\lambda_{21}=\tfrac12,\lambda_{22}=-1\) and the lane pattern \(\Xi^{(20)}=\epsilon\Xi_1,\Xi^{(21)}=\tfrac{\epsilon}{2}\Xi_1,\Xi^{(22)}=-\epsilon\Xi_1\); (D5) anisotropy variables \(\bar\Xi=0,a_\Xi=\tfrac{\epsilon}{4}\Xi_1,b_\Xi=\tfrac{3\epsilon}{4}\Xi_1\), hence \(b_\Xi=3a_\Xi\); (D6) the physical identification \(\Xi_1=P_1/P_0\) via the Stage-241 relation \(P_1/P_0=N_{01}/N_0-D_{01}/D_0=\Xi_{\rm load}\); (D7) the quadrupole-normalization lane pattern \(\Delta_Q^{(2m)}\); (D8) rigidity (\(\mathfrak i_r=\mathfrak h_r=0\Rightarrow\Xi_1=2\sum_r\rho_r^{(N)}\mathfrak m_r\)) and dominant-port corollaries. The appendix §"Outgoing-load scalar" anchors the load-factor factorization \(\Lambda_r^2/K=\mathcal M_r^2(1+\mathcal I_r)^2/(1-\mathcal H_r)^2\) (eq:app-part05-load-factor-factorization) as the bridge between raw load factor and the \((\mathcal M,\mathcal I,\mathcal H)\) decomposition.

## What the script claims to verify

The SymPy docstring lists four checks: (1) the weak-axisymmetric grouped slopes of \(\mathcal M_r,\mathcal I_r,\mathcal H_r\); (2) the per-port collapse \(\Sigma_{A,r}=\epsilon\lambda_A\sigma_r\); (3) the grouped trace/anomaly law \(\Xi_A=\epsilon\lambda_A\Xi_1\); (4) the outgoing-prefactor slope identification. The assertions match: block 1 builds a multiplicative log-linear perturbation `exp(eps*lam*slope)` of each microscopic variable and confirms `dln M/I/H` equals `lam*m_r/i_r/h_r`; block 2 confirms the order-eps slope of `log(Lambda_p^2/Kp / (Lambda^2/K))` equals `lam*sigma_r` with `Lambda=(OU2*GW+R*GU)/(OU2*OW2-R^2)`; block 3 hardcodes the lane pattern `(1,1/2,-1)` and confirms the standard grouped trace/anomaly projections give `(0, eps*Xi1/4, 3*eps*Xi1/4)` and `b=3a`; block 4 confirms `P1/P0 = n1-d1 = Xi_load` for `P=N/D` with log-linear `N,D`. The Mathematica `.wl` mirrors all four blocks identically (differing only by using `D[...,eps]/.eps->0` where SymPy uses `series/eps`).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 slopes m_r,i_r,h_r | `dln M/I/H == lam*m_r/i_r/h_r` (sympy 67-69 / wl 57-59) | match |
| D2 per-port amplitude sigma_r, Sigma=eps lam sigma_r | `Sigma_exact - lam*sigma_r == 0` (sympy 89 / wl 64) | match |
| D3 grouped collapse Xi_1 = sum_r rho_r^(N) sigma_r | none — `Xi1` is an abstract free symbol; the weighted port sum is never formed | missing (trivial-by-linearity) |
| D4 lane pattern (1,1/2,-1) | hardcoded `lam20/21/22` (sympy 94-96 / wl 69-71) | match (value-consistent with appendix signature) |
| D5 anisotropy bar=0,a=eps Xi1/4,b=3 eps Xi1/4, b=3a | sympy 106-109 / wl 80-83 | match |
| D6 Xi_1 = P1/P0 (microscopic tie) | `P1/P0 - Xi_load == 0` with Xi_load:=n1-d1, abstract (sympy 121-122 / wl 94-95) | partial (verifies only the definitional log-quotient identity; never ties to sigma_r/Xi_1) |
| D7 quadrupole-normalization Delta_Q lane pattern | none (structurally identical to D4) | missing (redundant with D4) |
| D8 rigidity / dominant-port corollaries | none | missing (trivial specializations of sigma_r) |
| Appendix load-factor factorization Lambda^2/K = M^2(1+I)^2/(1-H)^2 | none verified standalone; assumed in both scripts | missing |

Dominant pattern: the nontrivial core identities (D1, D2, D4, D5) are faithfully and non-tautologically exercised; D3/D6/D7/D8 are either trivial linear consequences or unverified specializations. No identity is verified with a *wrong* constant or a *different* target — so this is `partial` alignment, not `misaligned`. No `paper_misalignment` finding.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 67-69 | `series(log(M_p/M))/eps - lam*m_r == 0` (×3) | D1 | yes |
| A2 | sympy | 89 | `series(log(Lambda_p^2/Kp / (Lambda^2/K)))/eps - lam*sigma_r == 0` | D2 | yes |
| A3 | sympy | 106 | `(Xi20+2Xi21+2Xi22)/5 == 0` | D5 (trace) | yes (given hardcoded D4) |
| A4 | sympy | 107-109 | `a_Xi==eps Xi1/4`, `b_Xi==3eps Xi1/4`, `b-3a==0` | D5 | yes (given hardcoded D4) |
| A5 | sympy | 121-122 | `P_slope - P0*(n1-d1)==0`, `P_slope/P0 - (n1-d1)==0` | D6 | partial (definitional log-quotient only) |
| A6 | mathematica | 57-59 | `D[log(M_p/M),eps]/.eps->0 - lam*mR == 0` (×3) | D1 | yes |
| A7 | mathematica | 64 | `D[log(lambdaP^2/kP / (lambda0^2/k)),eps]/.eps->0 - lam*sigmaR == 0` | D2 | yes |
| A8 | mathematica | 80-83 | trace/anomaly projections | D5 | yes (given hardcoded D4) |
| A9 | mathematica | 94-95 | prefactor-slope identity | D6 | partial |

Every "yes" assertion would change value if the corresponding transport law / amplitude / signature were perturbed away from the paper's stated form, so A1-A4 and A6-A8 are non-tautological. A5/A9 cannot fail (log of a quotient = difference of logs by construction) — see F3.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl:26-104`
- (corresponding) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py:32-130`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:
- Background invariants — wl 32-35 `mCal=gw/(Sqrt[k]*ow2); iCal=r*gu/(ou2*gw); hCal=r^2/(ou2*ow2); lambda0=(ou2*gw+r*gu)/(ou2*ow2-r^2)` vs py 39-42 `Mcal=GW/(sp.sqrt(K)*OW2); Ical=R*GU/(OU2*GW); Hcal=R**2/(OU2*OW2); Lambda=(OU2*GW+R*GU)/(OU2*OW2-R**2)` — identical variable choreography.
- Perturbation — wl 37-42 `kP=k*Exp[eps*lam*kappa1]; ...` vs py 45-50 `Kp=K*sp.exp(eps*lam*kappa1); ...` — same multiplicative `exp` construction in the same order.
- Slopes & amplitude — wl 49-51/63 `mR=gW1-oW1-kappa1/2; ...; sigmaR=2*mR+2*iCal*iR/(1+iCal)+2*hCal*hR/(1-hCal)` vs py 58-60/84-88 — identical.
- Lane pattern / projections / prefactor block / carry-forward prints — wl 69-104 vs py 94-130 — identical, down to byte-for-byte print strings.
The only difference is the slope-extraction primitive (`D[...,eps]/.eps->0` vs `series(...,eps,0,2)/eps`), a trivial syntactic swap. Per the second-engine policy both engines must derive the result independently from the physical premises; here they echo the same algebra. (Bundled cleanup: both scripts mislabel the banner as "STAGE 160" — sympy:32, wl:26 — and both transcripts echo it; this is the wrong stage number.)

**Why this matters:**
The second engine is supposed to be a genuine cross-check. A transliteration cannot catch a setup error in the shared algebra (e.g., a wrong `Lambda` definition or a wrong `sigma_r` weight) because both scripts encode the identical expression; an error would pass in both. The independence guarantee is the whole point of running two engines.

**Required change:**
Differentiate the Mathematica derivation by adding an independent, non-shared verification and using the factored route the appendix actually states. Concretely (see directive F1 for exact edits): in the `.wl`, add a standalone check of the load-factor factorization `lambda0^2/k == mCal^2 (1+iCal)^2/(1-hCal)^2` (appendix eq:app-part05-load-factor-factorization), which neither script currently verifies, and re-derive the per-port amplitude slope from that *factored* form rather than from the raw `lambda0`. Also correct the banner string in both scripts to "STAGE 177".

**Verification:**
A new `expectZero["load-factor factorization", lambda0^2/k - mCal^2*(1+iCal)^2/(1-hCal)^2]` line should appear in the `.wl` output reading `= 0` / `PASS`, the `.wl` should still exit 0, and both transcripts should print the "STAGE 177" banner.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py:91-100`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl:66-74`

**What's wrong:**
The stage's headline (card `Output`, stage title, boxed notes line 186-192) is the *grouped collapse to one scalar* \(\Xi_{\rm load}^{(A)}=\epsilon\lambda_A\Xi_1\) with \(\Xi_1=\sum_r\rho_r^{(N)}\sigma_r\). The script never forms this weighted port sum: in block 3 `Xi1 = sp.symbols("Xi1")` (py:93 / wl:67) is an abstract free symbol, disconnected from the per-port `sigma_r` verified in block 2. The corollaries D7 (\(\Delta_Q\) lane pattern) and D8 (rigidity / dominant-port limits) are also untested. The per-port amplitude `Sigma_{A,r}=eps*lam*sigma_r` (block 2) IS verified, and `Xi_load=sum_r rho_r Sigma_r` is a definitional Stage-244 carry-forward, so the collapse is a trivial linearity step — which is why no check was written and why no non-tautological check *can* be written for it (see Self-test notes).

**Why this matters:**
A reader trusting the transcript would believe the literal "collapse to one scalar" is machine-verified, when only the per-port amplitude is. The gap is low-risk because the missing step is provably trivial (linear in `sigma_r`), but it leaves the card's bottom-line phrasing one inferential step ahead of the script.

**Required change:**
None applied by Codex (any added "sum_r rho_r sigma_r == Xi_1" check would be tautological — it reduces to linearity of a sum, see Self-test trap 3). This finding is informational: it documents that the headline collapse rests on a trivial linear step downstream of the verified per-port amplitude. Do NOT add a tautological collapse assertion. No directive edit for F2.

**Verification:**
No code change; the report records that D3/D7/D8 are trivial consequences of the verified D2/D4.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py:111-122`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl:85-95`

**What's wrong:**
The "Outgoing-prefactor slope identification" block defines `Xi_load = n1 - d1` (py:120 / wl:93) and `PA = NA/DA` with `NA=N0*exp(eps*lam*n1)`, `DA=D0*exp(eps*lam*d1)`, then asserts `P_slope/P0 - Xi_load == 0`. Since the order-eps log-slope of `N0 e^{ελn1}/(D0 e^{ελd1})` is `n1-d1` by construction, the assertion cannot fail — it confirms only that "the log of a quotient is the difference of the logs." The paper's deliverable D6 (\(\Xi_1=P_1/P_0\)) requires tying the *physical* prefactor slope to the *microscopic* `Xi_1 = sum_r rho_r sigma_r`; with `n1,d1` left as free abstract symbols this connection is never exercised.

**Why this matters:**
The check passes regardless of the stage's physics, so it provides no adversarial signal on D6. It encodes the Stage-241 carry-forward `P1/P0 = N01/N0 - D01/D0` in log-slope form (which is fine as documentation), but does not verify the nontrivial identification with the collapsed scalar.

**Required change:**
None applied by Codex. The nontrivial content of D6 lives in the (trivial-by-linearity) collapse of F2 and in Stage 241's carry-forward, not in an identity available within this stage's symbol scope; strengthening it would require either re-deriving Stage 241 (out of scope) or a tautological sum (F2). Documented as a known weak check; no directive edit for F3.

**Verification:**
No code change; report notes that A5/A9 are non-adversarial definitional identities.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py` (see F1): identical background invariants (wl 32-35 ↔ py 39-42), identical multiplicative perturbation (wl 37-42 ↔ py 45-50), identical `sigma_r` (wl 63 ↔ py 84-88), identical hardcoded lane pattern and projection weights (wl 69-78 ↔ py 94-104), and byte-identical carry-forward print blocks (wl 99-104 ↔ py 125-130). Only the slope-extraction primitive differs (`D[...,eps]/.eps->0` vs `series/eps`). This is `mathematica_transliteration`, F1.

## Engine cross-check

Both transcripts report every check `= 0` and (Mathematica) `PASS`, with exit code 0:
- `weak-axisymmetric d ln M/I/H = 0` (sympy out 13-15 / wl out 13-18)
- `Sigma_{A,r} = lambda_A sigma_r = 0` (sympy out 20 / wl out 23)
- `grouped trace vanishes = 0`, `a_Xi - eps Xi1/4 = 0`, `b_Xi - 3 eps Xi1/4 = 0`, `b_Xi - 3 a_Xi = 0` (sympy out 25-28 / wl out 29-36)
- `P_A slope = P0 * Xi_load = 0`, `(P1/P0) - Xi_load = 0` (sympy out 33-34 / wl out 41-44)
The engines agree everywhere. (This agreement is weakened as evidence by F1: they agree because they encode the same algebra.)

## Verdict justification

Within the scope the scripts actually exercise, the core nontrivial identities hold up under attack: the multiplicative-perturbation slope transport (D1) genuinely tests `m_r/i_r/h_r` against an independently parameterized background; the per-port amplitude check (D2) is non-tautological because `Lambda_p` is built from independent perturbed microscopics while `sigma_r` is built from `(M,I,H,m,i,h)`, and I confirmed by hand that `Lambda^2/K = M^2(1+I)^2/(1-H)^2` so the slope match is a real cross-check; the anisotropy projections (D5) reproduce the boxed notes values `(0, eps Xi1/4, 3 eps Xi1/4)` and `b=3a` given the paper's hardcoded signature (D4), which itself matches the appendix exactly. I tried to break the constants (lane pattern, slope coefficients, amplitude weights) against the notes/appendix and found no mismatch — so the verdict is `findings`, not `paper_misalignment` or `stop_cold`. The findings are: F1 (the Mathematica engine is a transliteration, plus a wrong "STAGE 160" banner) — the one fixable, non-tautological script-side item; and F2/F3 (the headline grouped collapse and the prefactor-slope identification rest on trivial-by-linearity / definitional steps that are documented but not — and cannot non-tautologically be — machine-checked within this stage's scope). Alignment is `partial`: every verified identity is correct and faithful, but several paper deliverables (D3/D6/D7/D8) are trivial consequences left unverified. Outputs are fresh (both `.txt` mtimes newer than their scripts).

## Self-test notes

Trap 1 (variable independence): the new factorization check `lambda0^2/k - mCal^2(1+iCal)^2/(1-hCal)^2` contains no derivatives; I verified the algebraic identity by hand (GW^2 cancels, OU2^2 cancels, OW2^2 cancels → `(ou2 gw + r gu)^2/(k(ou2 ow2 - r^2)^2) = lambda0^2/k`), so `expectZero` returns 0 non-trivially and would fail if `lambda0` were defined inconsistently with `(M,I,H)`. Trap 3 (trivial-case / anti-tautology): I confirmed that the F2 "collapse" check `sum_r rho_r (eps lam sigma_r) - eps lam sum_r rho_r sigma_r` reduces to 0 by linearity for any profile, i.e. it is tautological — so I deliberately did NOT prescribe it, and F2/F3 carry no Codex edits. Trap 5 (paper round-trip): the only prescribed edit (F1 factorization check + banner) introduces no new constant — `M,I,H,Lambda` are exactly the script's existing definitions and match the appendix factorization equation, so no new `paper_misalignment` is created.
