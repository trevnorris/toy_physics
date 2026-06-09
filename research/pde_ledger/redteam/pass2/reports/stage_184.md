---
unit_id: 184
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage184_branch_invariant_coordinates.md]
  paper_appendix: present
---

# Audit unit 184 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_184.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage184_branch_invariant_coordinates.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 99 ledger row; subsection "Branch-observable coordinates" lines 881-916)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.txt`

## What the paper claims

Stage 184 (`\stagefield{Output}`: "Identifies exact branch composites \(\mathfrak T_*\), \(\mathfrak N_*\), and \(\epsilon_\eta\) whose logarithmic drifts are the three normal-form coordinates.") repackages the three Stage-183 branch-adapted defect coordinates \((\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)\) as first-order logarithmic drifts of three exact coherent-branch composites. The deliverables (notes §§2-6, appendix eqs. `app-part05-Rtr-def`, `app-part05-Nstar-def`, `app-part05-Tstar-def`, `app-part05-branch-coordinate-laws`) are: (i) the positive constants \(B_*=2(1+\chi_{0,*}+\delta_{U,*})/\delta_{U,*}\) and \(C_*=(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})/(\chi_{0,*}\delta_{U,*})\); (ii) the three branch-coordinate laws \(\delta\ln\mathfrak T_*=\Sigma_{\rm tr}\), \(\delta\ln\mathfrak N_*=\Sigma_{\rm nt}\), \(\delta\ln\epsilon_\eta=\Sigma_\eta\) with \(\mathfrak T_*:=R_{\rm tr}^{-C_*}\), \(\mathfrak N_*:=\mathcal T^2R_{\rm tr}^{B_*}\); (iii) the exact selected-branch demand identity \(R_{\rm target}\mathcal T^2=\Lambda_0(1-\epsilon_\eta)\), \(\Lambda_0=27\pi^2Gc_s^5/(20a^5c^5)\), with its complementary drift \(\delta\ln[(R_{\rm target}\mathcal T^2)/\Lambda_0]=-\tfrac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta\); and (iv) the rewritten triple-rigidity equivalence. The notes carry the full derivation; the .tex card is terse and routes detail to the appendix and source.

## What the script claims to verify

The SymPy docstring (lines 8-14) enumerates five checks: the product identity \(R_{\rm target}\mathcal T^2=\Lambda_0(1-\epsilon_\eta)\); \(\delta\ln\mathfrak T_*=\Sigma_{\rm tr}\); \(\delta\ln\mathfrak N_*=\Sigma_{\rm nt}\); \(\delta\ln\epsilon_\eta=\Sigma_\eta\); and the triple-invariance theorem. The script pins `Bstar`/`Cstar` symbolically (lines 39-40), takes `SigmaTr = -Cstar*Theta1` and `SigmaNT = Xi1 + Bstar*Theta1` as Stage-183 carry-forwards (lines 52-53), builds first-order-perturbed `Rtr`, `T2`, `eps_eta_var` (lines 56-58), and verifies that the order-`small` coefficient of `log(composite/composite0)` matches each carry-forward coordinate (lines 67-81). It also checks the selected-branch complement (lines 83-89), three redundant restatements (lines 94-96), and three zero-map substitutions (lines 100-105). The Mathematica `.wl` mirrors this exactly.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `B_* = 2(1+χ0+δU)/δU` | py:39 `Bstar`, output shows `2*(chi0+deltaU+1)/deltaU` | match |
| `C_* = (1+χ0)(1+δU)(1+χ0+δU)/(χ0 δU)` | py:40 `Cstar`, output matches | match |
| `δln T_* = Σ_tr`, `T_* = R_tr^{-C_*}` | py:65-69 | match |
| `δln N_* = Σ_nt`, `N_* = T² R_tr^{B_*}` | py:72-76 | match |
| `δln ε_η = Σ_η` | py:79-81 | match |
| `R_target T² = Λ0(1-ε_η)` (exact branch identity) | py:62 | match-but-tautological (see F2) |
| `δln[(R_target T²)/Λ0] = -(ε_η/(1-ε_η))Σ_η` | py:83-89 | match |
| triple-rigidity equivalence | py:94-105 (forward zero-map only) | partial (forward direction only; iff not shown, but acceptable — see verdict) |
| `Λ0 = 27π²Gc_s⁵/(20a⁵c⁵)` explicit form | not used; `Lam0` is opaque positive symbol | n/a (inert at linear order, correctly elided) |

`paper_alignment: aligned` — every load-bearing deliverable has a faithful, value-correct script counterpart. No `value_mismatch`, no `target_mismatch`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `Rtarget*T2 - Lam0*(1-eps_eta_var) == 0` | product identity (iii) | no (tautological — `Rtarget` defined as `Lam0*(1-eps_eta_var)/T2` at line 59) |
| A2 | sympy | 69 | `dln_Ttr - SigmaTr == 0` | law `δln T_*=Σ_tr` | yes |
| A3 | sympy | 76 | `dln_Ntr - SigmaNT == 0` | law `δln N_*=Σ_nt` | yes |
| A4 | sympy | 81 | `dln_eps_eta - SigmaEta == 0` | law `δln ε_η=Σ_η` | yes |
| A5 | sympy | 87-90 | `dln_Ecomp + eps_eta/(1-eps_eta)*SigmaEta == 0` | complement drift (iii) | yes |
| A6-A8 | sympy | 94-96 | `SigmaTr - dln_Ttr`, etc. == 0 | same as A2-A4 | redundant (sign-flipped duplicates of A2-A4) |
| A9-A11 | sympy | 100-105 | zero-map substitutions == 0 | forward direction of triple-rigidity | weak (follows trivially from linearity already established by A2-A4) |
| B1-B11 | mathematica | 48,55,62,67,73,76-78,83-85 | same checks, `expectZero`/`SeriesCoefficient` | same as A1-A11 | same anchoring as A-row counterparts |

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.txt` (mtime May 30 01:41)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.txt` (mtime May 30 01:41)
- both scripts mtime Jun 3 15:59

**What's wrong:**
Both saved `.txt` outputs predate their scripts by ~4 days. The staleness is not cosmetic: both transcripts open with `STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES` (sympy .txt:3, math .txt:3), whereas the *current* script banner is `STAGE 184 — EXACT BRANCH-INVARIANT COORDINATES` (py:35, wl:26). The committed outputs were therefore produced by an earlier revision whose banner still read "167" (the pre-renumber +17 label, 184−17=167). The Mathematica footer (`Stage 184 Mathematica audit passed.`, .txt:65) already says 184, confirming the banner line is the only stale fragment of the header and the body math is otherwise identical to what the current script would emit.

**Why this matters:**
The committed transcript does not reflect the current script's header. The body residuals are unaffected (all checks still resolve to 0 / PASS), so the math conclusion stands, but a fresh re-run is needed so the committed output banner matches the source. This is the verifier's standard fresh-run trigger.

**Required change:**
None for Codex to author. The orchestrator's independent re-run will regenerate both `.txt` with the corrected `STAGE 184` banner. (NUMBERING NOTE for orchestrator only — do not fix in this audit: the stale "STAGE 167" banner in both committed outputs is a +17 pre-renumber artifact; it is already corrected in the live `.py`/`.wl` source.)

**Verification:**
After fresh run, both `.txt` line 3 should read `STAGE 184 — EXACT BRANCH-INVARIANT COORDINATES`; all PASS lines unchanged.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py:59-62`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl:47-48`

**What's wrong:**
The product identity `R_target * T^2 = Lambda0 * (1 - eps_eta)` (paper deliverable iii, notes §1.3, appendix eq. line ~728) is checked at py:62 / wl:48, but `Rtarget` is *defined* one line earlier as exactly that quotient:
```
59  Rtarget = Lam0 * (1 - eps_eta_var) / T2
62  expect_zero("R_target * T^2 - Lambda0 * (1 - eps_eta)", sp.simplify(Rtarget * T2 - Lam0 * (1 - eps_eta_var)))
```
`Rtarget * T2` is therefore `Lam0*(1-eps_eta_var)` by construction; the subtraction is `Lam0*(1-eps_eta_var) - Lam0*(1-eps_eta_var) ≡ 0`. The assertion cannot fail regardless of physics. The paper presents this identity as a *premise* of the selected-branch description, not as something derived here — so a script check that re-derives it from its own definition verifies nothing. (Note: the *downstream* complement check A5 at py:87-90 is genuinely substantive — it linearizes `log((1-eps_eta_var)/(1-eps_eta))` — so the tautology is confined to the standalone product-identity line.)

**Why this matters:**
This is one of the five docstring-enumerated checks (line 9). As written it gives false assurance that the script independently confirms the branch identity. It does not. The other four checks are sound, so this is a quality/credibility defect rather than a math error.

**Required change:**
Replace the round-trip with a check anchored to an *independent* construction of `R_target`. Concretely, do not derive `Rtarget` from `Lam0*(1-eps_eta_var)/T2`; instead introduce `Rtarget` as its own perturbed symbol `Rtarget = Rtarget0 * exp(small * R1)` with `Rtarget0` a free positive symbol, and verify the *drift* relation the paper actually asserts — `δln(R_target T²) = δln(1-eps_eta)` i.e. the order-`small` coefficient of `log((Rtarget*T2)/(Rtarget0*T20))` equals that of `log((1-eps_eta_var)/(1-eps_eta))` — which is the genuine content of notes eq. lines 120-127. If instead the intent is only the algebraic identity at the *background* level, state explicitly in a comment that `Rtarget` is *defined* by the identity and remove the self-referential `expect_zero`, since a definition needs no verification.

**Verification:**
The new check at py:62 (and wl:48) must reference a `Rtarget` symbol that is NOT algebraically `Lam0*(1-eps_eta_var)/T2`; the residual it drives to zero must be a non-trivial linearization (the verifier confirms the residual expression printed is a function of `R1`/`Xi1`/`SigmaEta`, not a literal `0` produced by cancellation of identical terms).

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl` (whole file)
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:

1. Same variable choreography and the *identical* construction of the perturbed objects:
   - py:56-59 `Rtr = Rtr0*exp(small*Theta1)` / `T2 = T20*exp(small*Xi1)` / `eps_eta_var = eps_eta*(1+small*SigmaEta)` / `Rtarget = Lam0*(1-eps_eta_var)/T2`
   - wl:44-47 `rTr = rTr0*Exp[small*theta1]` / `t2 = t20*Exp[small*xi1]` / `epsEtaVar = epsEta*(1+small*sigmaEta)` / `rTarget = lam0*(1-epsEtaVar)/t2`

2. Same drift extraction by the *same* black-box on the *same* object:
   - py:67 `sp.series(sp.log(Ttr/Ttr0), small, 0, 2).removeO().coeff(small,1)`
   - wl:53 `SeriesCoefficient[Log[tTr/tTr0], {small,0,1}]`
   (identical: log of the same ratio, extract the order-1 coefficient in `small`)

3. Identical assertion list in identical order with identical names: py:62,69,76,81,87-90,94-96,100-105 map one-to-one to wl:48,55,62,67,73,76-78,83-85, including the same redundant restatements (A6-A8) and the same forward-only zero-maps (A9-A11). `Cstar`/`Bstar` are even built by the same factored expression (py:39-40 ↔ wl:31-32).

There is no second, structurally-distinct route (e.g. Mathematica could close the same identities via `D[..., small]` / `Limit`, or via `PowerExpand`/`Logarithmic` manipulation of the exponents, or by directly differentiating the closed-form composites symbolically). Instead it echoes SymPy's exponential-ansatz-plus-`Series`-coefficient algebra verbatim. Per the second-engine policy this is `mathematica_transliteration`.

**Why this matters:**
A transliterated second engine cannot independently catch an algebra error in the first — both will make the same `Series`/`coeff` choice on the same constructed object. The stage's value laws are believable, but the "two engines" do not constitute two independent confirmations.

**Required change:**
Re-author the `.wl` to derive at least the three branch-coordinate laws by an independent Mathematica route. Suggested independent route: define the composites in closed form (`tTr = rTr^(-cStar)`, `nTr = t2*rTr^bStar`) and compute `δln = D[Log[composite], small] /. small -> 0` (a derivative-based linearization) rather than `SeriesCoefficient[..., {small,0,1}]` of the ratio; or use `Normal@Series` with `PowerExpand`/`Logarithmic`-`Collect` so the exponent `bStar`/`cStar` is brought down by log-rules rather than by the same series machinery SymPy used. The final residuals must still be `0` but reached by a path that does not mirror py:67/74/79 token-for-token.

**Verification:**
The verifier confirms the `.wl` no longer uses `SeriesCoefficient[Log[ratio], {small,0,1}]` as the sole drift extractor for all three laws, and that at least one law is closed via a distinct construct (`D[...]`/`Limit`/`PowerExpand`), while all `expectZero` checks still PASS.

## Independent-derivation check (Mathematica)

Not independent — see F3. The `.wl` reproduces the SymPy script's exponential-ansatz construction (wl:44-47 ↔ py:56-59) and extracts every drift with `SeriesCoefficient[Log[ratio], {small,0,1}]` (wl:53,60,65,71), the direct Mathematica spelling of SymPy's `series(log(ratio), small, 0, 2).coeff(small,1)` (py:67,74,79,85). Assertion order, names, the redundant restatements, and the forward-only zero-maps are all 1:1. Strongest single piece of evidence: the constructed object `Rtarget`/`rTarget` is built by the identical quotient `lam0*(1-epsEtaVar)/t2` (wl:47) ≡ `Lam0*(1-eps_eta_var)/T2` (py:59), so both engines inherit the *same* tautology in F2 simultaneously — a hallmark of a port, not an independent derivation.

## Engine cross-check

Both engines AGREE. SymPy output and Mathematica output report identical residuals for every check:
- `B_*`, `C_*` identical (modulo print formatting): `2*(chi0+deltaU+1)/deltaU` vs `(2*(1+chi0+deltaU))/deltaU`; same `C_*`.
- `δln T_*` = `-((1+chi0)(1+deltaU)(1+chi0+deltaU)theta1)/(chi0 deltaU)` (math .txt:17) = SymPy's expanded `Theta_1*(-chi0**2*deltaU-...)/(chi0*deltaU)` (sympy .txt:16) — same rational function, one expanded, one factored.
- `δln N_*` matches (sympy .txt:22 vs math .txt:24).
- All `... = 0` / `PASS` lines coincide.
No `engine_disagreement`. (Agreement here is expected and low-value precisely because the engines are not independent — see F3.)

## Verdict justification

The stage's load-bearing physics — the three branch-coordinate laws `δln T_*=Σ_tr`, `δln N_*=Σ_nt`, `δln ε_η=Σ_η`, the constants `B_*`/`C_*`, and the selected-branch complement drift — is correctly and non-tautologically verified by checks A2-A5 (and their Mathematica twins), and every value reconciles exactly with the notes and appendix. Paper alignment is exact (no value or target mismatch). Three defects keep this from `clean`: (F1) both committed transcripts are stale and still carry the pre-renumber "STAGE 167" banner; (F2) the standalone product-identity check is tautological because `R_target` is defined by the very identity it then "verifies"; and (F3) the Mathematica script is a line-by-line port of the SymPy script rather than an independent second engine, so the two-engine guarantee is not actually met. None of these affect the truth of the stage's results, but F2 and F3 are genuine verification-quality defects and F1 needs a fresh run. Verdict: `findings`.

## Self-test notes

Variable-independence trap: checked every `series(log(ratio), small, ...)` — each composite genuinely depends on `small` via the exponential ansatz, so the order-1 coefficients are non-trivial (not identically zero); A2-A5 are real. Trivial-case trap: confirmed F2's residual cancels to literal `0` by identical-term subtraction (the tautology) while A5's residual is a genuine function of `eps_eta`/`SigmaEta`. Round-trip trap: confirmed `Rtarget*T2` reduces to its own definition at py:62 (the F2 finding). No symmetry/parity integrals in this stage. No symbol-assumption error: all positivity assumptions (`chi0,deltaU,eps_eta,Rtr0,T20,Lam0 > 0`) match the paper's positive-coupling state-space premise (card `\stagefield{Inputs}`).

## Value Reconciliation (pass-2 augmentation)

**Deliverable-level table:**

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `B_* = 2(1+χ0+δU)/δU` | py:39, wl:31; sympy .txt:5, math .txt:5 | notes:218-220 `B_*:=2(1+χ_{0,*}+δ_{U,*})/δ_{U,*}`; appendix:896 eq.`app-part05-Nstar-def` | MATCH |
| `C_* = (1+χ0)(1+δU)(1+χ0+δU)/(χ0 δU)` | py:40, wl:32; sympy .txt:6, math .txt:6 | notes:159-161 `C_*:=...`; appendix:903 eq.`app-part05-Tstar-def` | MATCH |
| `T_* = R_tr^{-C_*}` | py:65,108; wl:51,89 | notes:166 `\mathfrak T_*:=R_{rm tr}^{-C_*}`; appendix:901 | MATCH |
| `N_* = T² R_tr^{B_*}` | py:72,109; wl:58,90 | notes:223-224 `\mathfrak N_*:=\mathcal T² R_{rm tr}^{B_*}`; appendix:894 | MATCH |
| `δln T_* = Σ_tr` | py:69; sympy .txt:17, math .txt:18 | notes:181-183, appendix:908 eq.`app-part05-branch-coordinate-laws` | MATCH |
| `δln N_* = Σ_nt` | py:76; sympy .txt:23, math .txt:25 | notes:239-241, appendix:910 | MATCH |
| `δln ε_η = Σ_η` | py:81; sympy .txt:29, math .txt:32 | notes:290-292 (`δln\mathfrak D=Σ_η`), appendix:912 | MATCH |
| `R_target T² = Λ0(1-ε_η)` | py:62, wl:48; sympy .txt:11 | notes:108-114 (§1.3 boxed); appendix:~728 | MATCH (value), but check is tautological → F2 (verification quality, not a value misalignment) |
| `δln[(R_target T²)/Λ0] = -(ε_η/(1-ε_η))Σ_η` | py:86-89; sympy .txt:30-31, math .txt:34 | notes:303-309 (`δln\mathfrak E`); appendix mentions complement | MATCH |
| `Σ_tr = -C_* Θ1` (Stage-183 input) | py:52; sympy .txt:16 form | notes:153-155 | MATCH |
| `Σ_nt = Ξ1 + B_* Θ1` (Stage-183 input) | py:53; sympy .txt:22 form | notes:210-213 | MATCH |
| `Λ0 = 27π²Gc_s⁵/(20a⁵c⁵)` | NOT emitted as a value (kept opaque as `Lam0`) | notes:114; appendix:710 | MATCH (correctly elided; inert at linear order — script uses opaque positive symbol, no value to reconcile) |

**INTERNAL scaffolding (no finding):** perturbation bookkeeping symbol `small`; background base symbols `Rtr0`/`rTr0`, `T20`/`t20`, `Lam0`/`lam0`; perturbed objects `Rtr`, `T2`, `eps_eta_var`, `Rtarget`, `Ttr`/`Ttr0`, `Ntr`/`Ntr0`, `Ecomp`/`Ecomp0`; the redundant restatement residuals (A6-A8) and the forward zero-map residuals (A9-A11), all of which are pass/fail scaffolding equal to 0.

**reconciliation: complete; 12 deliverable values checked, 0 misaligned.** Every emitted deliverable value reflects correctly in both the notes `.md` and the part-05 appendix `.tex`. The F2 (tautological) and F3 (transliteration) findings are *verification-quality* defects, not value mismatches — no `value_mismatch` or `script_missing_paper_claim` arises from this reconciliation.
