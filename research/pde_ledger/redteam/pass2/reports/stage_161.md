---
unit_id: 161
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage161_dn_similarity_slippage.md]
  paper_appendix: present
---

# Audit unit 161 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_161.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage161_dn_similarity_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read subsec `app-part04-bare-slippage`, lines ~1064–1138, and the MTDC-T8 anchor/firewall rows)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.txt`

## What the paper claims

Stage 161 decomposes the tangential mixed-port slippage susceptibility `\Upsilon_\Pi` into its
microscopic D/N-similarity content. The card's audit target (the boxed body equation) is
"`\Delta_Q=-\sigma_*\Xi_{\rm slip}\delta\Pi_{\rm tan}/(1-\sigma_*)`, with
`\Xi_{\rm slip}=\Xi_\gamma-2\Xi_L`" (stage_161.tex:16). The notes enumerate the supporting
deliverables: (1) the exact bare-slippage decomposition `\mathfrak B_W=((1+r_c)/9)(\varepsilon_\gamma-\varepsilon_\kappa)`
(notes §1, boxed); (2) its linearization `\delta\mathfrak B_W=((1+r_{c,*})/9)(\delta\varepsilon_\gamma-\delta\varepsilon_\kappa)`;
(3) the exact D/N-tube even defect `\varepsilon_\kappa=12L_W^2/(\pi^2a^2(1+r_c))-1` and its first variation
`\delta\varepsilon_\kappa=2\,\delta\ln(L_W/a)-\delta\ln(1+r_c)` (§2); (4) the odd defect
`\varepsilon_\gamma=9\gamma_0/(1+r_c)-1` with `\delta\varepsilon_\gamma=\delta\ln\gamma_0-\delta\ln(1+r_c)` (§3);
(5) the `\delta r_c`-cancellation difference identity
`\delta\varepsilon_\gamma-\delta\varepsilon_\kappa=\delta\ln\gamma_0-2\,\delta\ln(L_W/a)` (§3);
(6) `\Upsilon_\Pi=((1+r_{c,*})/9)(\Xi_\gamma-2\Xi_L)` with the numeric prefactor
`(1+r_{c,*})/9\approx0.462362334687869` (§4); (7) the collapse of the Stage-160 defect law to
`\Delta_Q=-(\sigma_*/(1-\sigma_*))\Xi_{\rm slip}\delta\Pi_{\rm tan}` after the `9` and `(1+r_c)`
prefactors cancel (§5); and (8) the D/N similarity-preservation theorem
`\Xi_\gamma=2\Xi_L\Rightarrow\Xi_{\rm slip}=0\Rightarrow\Delta_Q=0` (§6). The card's `\stagefield{Checks}`
list maps to (5)/(6)/(8). The appendix (eq. `app-part04-DeltaQ-Xislip`, `app-part04-Xislip-def`,
`app-part04-deltaPi-tan`) restates the same closed forms and the `\delta\Pi_{\rm tan}` transport
coefficients `0.832409471081635`, `1.16275838754222`.

## What the script claims to verify

The SymPy docstring lists six checks that mirror notes deliverables (1)–(6)/(8): the exact
similarity-defect decomposition, the linearized law, the D/N-tube `\varepsilon_\kappa` and its
variation, the `\delta r_c`-cancellation difference identity, the collapse of the final defect law,
and the preservation theorem `\Xi_\gamma=2\Xi_L\Rightarrow\Delta_Q=0`. Each is enforced by an
`expect_zero` (SymPy) / `expectZero` (Mathematica) that simplifies a residual and aborts on nonzero.
The scripts also PRINT the numeric prefactor `(1+r_{F1}^2)/9` and the `\delta\Pi_{\rm tan}`→`\delta\widehat T_m`
expanded `\Delta_Q` coefficients as labeled (non-asserted) carry-forward results.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `\mathfrak B_W=((1+r_c)/9)(\varepsilon_\gamma-\varepsilon_\kappa)` | py:44 / wl:35 `expect_zero("exact similarity-defect decomposition", …)` | match |
| (2) linearized `\delta\mathfrak B_W` | py:52 / wl:43 | match |
| (3) `\varepsilon_\kappa` + `\delta\varepsilon_\kappa=2\delta\ln(L_W/a)-\delta\ln(1+r_c)` | py:58–68 / wl:50–63 `"d eps_kappa identity"` | match |
| (4) `\varepsilon_\gamma` + `\delta\varepsilon_\gamma=\delta\ln\gamma_0-\delta\ln(1+r_c)` | py:73–85 / wl:67–76 | match |
| (5) difference identity (`\delta r_c` cancels) | py:86–93 / wl:77–89 `"difference identity"` | match |
| (6) `\Upsilon_\Pi=((1+r_{c,*})/9)\Xi_{\rm slip}`, prefactor `0.462362…` | py:100/135–137 / wl:96/122–124 (printed) | match (printed, not asserted) |
| (7) collapse `\Delta_Q=-(\sigma_*/(1-\sigma_*))\Xi_{\rm slip}\delta\Pi_{\rm tan}` | py:107–114 / wl:103–110 `"collapsed Delta_Q law"`/`"collapsed N_Q-1 law"` | match (but construction-guaranteed — see Verdict) |
| (8) `\Xi_\gamma=2\Xi_L\Rightarrow\Delta_Q=0` | py:126–133 / wl:118–120 | match |
| `\delta\Pi_{\rm tan}` transport coeffs `0.83240…`, `1.16275…`, `6.42981…`→`5.35223…` | py:117–123 / wl:112–116 (printed) | match (printed) |

`paper_alignment: aligned`. Every paper-side deliverable has a corresponding script-side check;
the values reconcile (see Value Reconciliation). The only defect is engine-independence, not
content.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero(B_W - (1+rc)(eps_g-eps_k)/9)` | deliverable 1 | yes |
| A2 | sympy | 52 | `expect_zero(dB_W - (1+rc*)(deps_g-deps_k)/9)` | deliverable 2 | yes |
| A3 | sympy | 68 | `expect_zero(d eps_kappa identity)` | deliverable 3 | yes |
| A4 | sympy | 82 | `expect_zero(d eps_gamma = dln g0 - dln(1+rc))` | deliverable 4 | yes |
| A5 | sympy | 93 | `expect_zero(difference identity)` | deliverable 5 | yes |
| A6 | sympy | 107 | `expect_zero(collapsed Delta_Q law)` | deliverable 7 | partial (construction-guaranteed) |
| A7 | sympy | 111 | `expect_zero(collapsed N_Q-1 law)` | deliverable 7 | partial (construction-guaranteed) |
| A8 | sympy | 126 | `expect_zero(Xi_g=2Xi_L => Delta_Q=0)` | deliverable 8 | yes |
| A9 | sympy | 130 | `expect_zero(Xi_g=2Xi_L => N_Q-1=0)` | deliverable 8 | yes |
| A10 | mathematica | 35 | `expectZero(bW - (1+rc)(epsG-epsK)/9)` | deliverable 1 | yes |
| A11 | mathematica | 43 | `expectZero(dBW - (1+rcStar)(depsG-depsK)/9)` | deliverable 2 | yes |
| A12 | mathematica | 63 | `expectZero(d eps_kappa identity)` | deliverable 3 | yes |
| A13 | mathematica | 73 | `expectZero(d eps_gamma = dln g0 - dln(1+rc))` | deliverable 4 | yes |
| A14 | mathematica | 89 | `expectZero(difference identity)` | deliverable 5 | yes |
| A15 | mathematica | 103 | `expectZero(collapsed Delta_Q law)` | deliverable 7 | partial (construction-guaranteed) |
| A16 | mathematica | 107 | `expectZero(collapsed N_Q-1 law)` | deliverable 7 | partial (construction-guaranteed) |
| A17 | mathematica | 119 | `expectZero(Xi_g=2Xi_L => Delta_Q=0)` | deliverable 8 | yes |
| A18 | mathematica | 120 | `expectZero(Xi_g=2Xi_L => N_Q-1=0)` | deliverable 8 | yes |

A6/A7 (and the Mathematica mirror A15/A16) are flagged "partial": `DeltaQ` is *defined*
(py:103) as `-9*sigma*UpsilonPi*dPi/((1-sigma)*(1+rc_star))` with `UpsilonPi=(1+rc_star)(Xi_g-2Xi_L)/9`,
so the assertion `DeltaQ + sigma*(Xi_g-2Xi_L)*dPi/(1-sigma) == 0` is the algebraic cancellation of the
`9` and `(1+rc_star)` factors the script itself introduced one line earlier. It cannot fail and cannot
detect an error in the carried Stage-160 prefactor — but the `9·(1+r_c)/9/(1+r_c)=1` cancellation it
confirms *is* the stage's stated collapse claim, so it is weak rather than wrong (no separate finding).

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:26-135`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Stage 161 sits in
the 105–175 orchestrator-direct band where this defect recurs. The choreography, intermediate
quantities, print labels, and even the inline comments are identical between the two engines:

- Same symbol setup and definitions, transliterated 1:1:
  py:40-41 `kappa0 = (1 + rc)*(1 + eps_kappa)/3; gamma0 = (1 + rc)*(1 + eps_gamma)/9`
  wl:31-32 `kappa0 = (1 + rc)*(1 + epsKappa)/3; gamma0 = (1 + rc)*(1 + epsGamma)/9`
- Same perturbation trick, same variable choreography:
  py:49-50 `BW_pert = BW.subs({eps_kappa: eps*deps_k, eps_gamma: eps*deps_g, rc: rc_star + eps*drc}); dBW = sp.diff(BW_pert, eps).subs(eps,0)`
  wl:40-41 `bWPert = bW /. {epsKappa -> eps*depsK, epsGamma -> eps*depsG, rc -> rcStar + eps*drc}; dBW = D[bWPert, eps] /. eps -> 0`
- Identical inline comment ported verbatim:
  py:81 `# On the branch, d gamma_0 = gamma_{0,*} d ln gamma_0 = (1+r_c) d ln gamma_0 / 9.`
  wl:72 `(* On the branch d gamma_0 = (1+r_c)/9 * d ln gamma_0. *)`
- Identical literal re-typing of the Stage-159 transport coefficients rather than re-deriving them:
  py:117 `0.832409471081635*dSigma0 - 1.16275838754222*dS`; py:122 `6.42981496203006*dThat`
  wl:112 `0.832409471081635*dSigma0 - 1.16275838754222*dS`; wl:115 `6.42981496203006*dThat`
- Same `UpsilonPi`/`deltaQ` construction so the "collapsed" assertion cancels by identical construction
  in both engines (py:100-110 ↔ wl:96-110).

The only genuine algebraic divergence is the branch reduction in the two D/N-tube identities — SymPy
uses `expr.subs(12*LW**2 - …, 0)` (py:67,92) while Mathematica uses `PolynomialRemainder[…, -12*lW^2 + a^2*Pi^2*(1+rc), lW]`
(wl:57-61,77-87). That is a different reduction mechanism for ~2 of the 9 checks, but the surrounding
derivation, the symbolic targets, and all printed/carried results are echoed line-for-line. Per the
second-engine policy both engines must derive the result independently from the physical premises; here
the second engine echoes the first engine's algebra.

**Why this matters:**
A transliterated second engine cannot catch a transcription error in the first — both would carry the
same mistake. The independence claim of the dual-engine pipeline is the whole point of the second
engine; a port provides false redundancy. (It does not change the math result; the verdict is `findings`,
not stop-cold.)

**Required change:**
Re-author the `.wl` so the load-bearing identities are reached by a route that does not mirror the
SymPy variable choreography step-for-step. Concretely: (a) for the bare-slippage decomposition and its
linearization, drive the linearization from a Taylor/Series expansion of `bW` rather than the
`eps -> 0` derivative substitution; (b) for the `\delta\varepsilon_\kappa`/`\delta\varepsilon_\gamma`
identities, reach the branch-reduced differentials via an explicit logarithmic-derivative substitution
(`D[Log[…], …]`) instead of re-typing the same `epsKExact`/`epsGExact` closed forms and the same
`12 lW^2 -> Pi^2 a^2 (1+rc)` replacement; (c) for the collapse, build `\Upsilon_\Pi` and the Stage-160
prefactor independently and confirm `\Delta_Q` reduces — not re-type the SymPy `deltaQ` definition.
The set of asserted identities (and their expected zeros) must stay identical; only the derivation
route changes so the engine is genuinely independent. The `.wl` must still `Exit[0]` with all
`expectZero` passing and the printed numeric `(1+r_{F1}^2)/9 = 0.46236233468786880…` and the
`Delta_Q in (dThat,dS)` coefficients unchanged.

**Verification:**
The re-authored `.wl` runs to `Stage 161 Mathematica audit passed.` with all `PASS:` lines, the
printed prefactor and `Delta_Q` coefficients match the prior transcript (engine agreement preserved),
and the script no longer mirrors the SymPy line structure (no `eps -> 0` derivative trick, no verbatim
`# On the branch …` comment port).

## Independent-derivation check (Mathematica)

Transliteration, as detailed in F1. The `.wl` reproduces the `.py`'s exact step sequence, variable
naming map (snake_case → camelCase), intermediate prints, and even comments, with the sole exception
of the `PolynomialRemainder` branch reduction (vs SymPy's `.subs(poly,0)`) used in two of the nine
checks. That single divergence is not enough to call the second engine independent across the stage;
the asserted targets and all carried literals are echoed. → `mathematica_transliteration` finding.

## Engine cross-check

Both engines are present and their final outputs agree:

- `exact similarity-defect decomposition = 0`, `linearized slippage law = 0`, `d eps_kappa identity = 0`,
  `d eps_gamma = … = 0`, `difference identity = 0`, `collapsed Delta_Q law = 0`, `collapsed N_Q-1 law = 0`,
  `Xi_gamma = 2 Xi_L => Delta_Q = 0`/`N_Q-1 = 0` — all zero in both transcripts.
- Prefactor: SymPy `(1+r_F1^2)/9 = 0.462362334687868748`; Mathematica `0.46236233468786880105…` — agree
  to 16 digits.
- `Delta_Q in (dThat,dS)` coefficients: SymPy `5.35223887169623`, `10.7044777433925`; Mathematica
  `5.352238871696225`, `10.70447774339245` — agree.

`engines_agree: true`.

## Verdict justification

The scripts correctly and non-tautologically verify deliverables (1)–(5) and (8) of the stage, and they
match the paper card, notes, and appendix exactly (no `168π²`/`100π²` or Family-1-radius content is in
play in this sub-topic, so that recurring stale-constant watch item is N/A here). I attacked: (a) the
`d eps_kappa`/`d eps_gamma`/`difference` identities for a zero-derivative trap — all derivatives are wrt
symbols the expressions genuinely depend on, so the `expect_zero`s are substantive, not trivially-true;
(b) the positivity assumptions (`LW,a` positive) — they do not invalidate the exact polynomial identities;
(c) the collapse assertion A6/A7 — it is construction-guaranteed (DeltaQ is defined from the same
UpsilonPi and prefactor it then cancels against), so it confirms the stated `9·(1+r_c)/9/(1+r_c)=1`
collapse but cannot independently validate the carried Stage-160 prefactor — weak but not a standalone
finding. The one real defect is engine independence: the `.wl` is a line-by-line transliteration of the
`.py` (F1, medium). Hence `verdict: findings`, not clean and not stop-cold.

## Self-test notes

Checked variable-independence: every `sp.diff`/`D` in the D/N-tube identities is taken wrt a symbol the
expression actually depends on (LW, a, rc, gamma0Sym), so no identically-zero-derivative trap (the
known prior failure mode) — the `expect_zero`s are real. No unbounded integrals, so no parity trap.
Trivial-case: the difference identity reduces to 0 only after the branch substitution `12 LW^2 → π²a²(1+rc)`,
which both engines apply, and the F1 re-authoring does not change any asserted target, so it introduces
no new `paper_misalignment` (paper round-trip clean).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `B_W = ((1+r_c)/9)(eps_gamma-eps_kappa)` | py:44 / wl:35; sympy out L5 | notes:81-88 (boxed §1) | MATCH |
| `dB_W = ((1+r_c*)/9)(deps_g-deps_k)` | py:52 / wl:43; out L7 | notes:93-100 (boxed §1) | MATCH |
| `eps_kappa = 12 LW^2/(π²a²(1+r_c)) - 1` | py:58 / wl:50; out L13 | notes:124-130 (boxed §2) | MATCH |
| `d eps_kappa = 2 dln(L_W/a) - dln(1+r_c)` | py:65 target / wl:56; out L14 | notes:143-151 (boxed §2) | MATCH |
| `eps_gamma = 9γ0/(1+r_c) - 1` | py:73 / wl:67; out L16 | notes:160-166 (boxed §3) | MATCH |
| `d eps_gamma = dln γ0 - dln(1+r_c)` | py:82 / wl:73; out L18 | notes:181-189 (boxed §3) | MATCH |
| `difference id: deps_g - deps_k = dln γ0 - 2 dln(L_W/a)` | py:88/93 / wl:81/89; out L19 | notes:193-200 (boxed §3) | MATCH |
| `Upsilon_Pi = ((1+r_c*)/9)(Xi_g - 2 Xi_L)` | py:100 / wl:96; out L24 | notes:236-242 (boxed §4); tex:16 | MATCH |
| `(1+r_F1^2)/9 = 0.462362334687869` | py:136 / wl:123; out L37 | notes:260 (`0.462362334687869`) | MATCH |
| `Delta_Q = -(σ*/(1-σ*)) Ξ_slip dΠ_tan` (collapse) | py:103/107 / wl:99/103; out L25/27 | tex:16 (boxed); notes:301-308; appendix eq DeltaQ-Xislip | MATCH |
| `Delta_Q in (dThat,dS)`: coeffs `5.35223887169622`, `−1.16275838754222` (×σ*/(1−σ*)) | py:122-123 / wl:115-116; sympy out L30 | notes:377/379 (`5.35223887169622`, `1.16275838754222`) | MATCH |

Sub-constants `0.832409471081635`, `1.16275838754222`, `6.42981496203006` (the Stage-159/160 transport
coefficients re-typed into `dPi_tan_expr`) also reconcile to notes:332/334/367. The `\delta\Pi_{\rm tan}`
coefficients additionally match appendix eq. `app-part04-deltaPi-tan` (lines 1116-1117).

INTERNAL (scaffolding, no finding expected in prose): `eps`/`drc`/`deps_k`/`deps_g` perturbation
helpers, `BW_pert`/`bWPert`, `depsk_direct`/`depsKDirect`, `depsk_target`/`depsKTarget`, the polynomial
reduction residuals, the `PASS:`/`FAIL:` flags, `$Assumptions`, and the printed carry-forward formula
list (1)–(6).

All emitted deliverable values reflect correctly in the notes (the natural carrier; the `.tex` card is
terse by design and carries the load-bearing `\Delta_Q`/`\Xi_{\rm slip}` forms). No MISMATCH or
MISSING-DELIVERABLE. The single audit finding (F1) is engine-independence, not a value discrepancy, so
this augmentation adds no `paper_misalignment` and does not change the verdict.
</content>
</invoke>
