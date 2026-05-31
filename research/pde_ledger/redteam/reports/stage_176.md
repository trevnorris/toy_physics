---
unit_id: 176
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage176_outgoing_load_factorization.md]
  paper_appendix: present
---

# Audit unit 176 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_176.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage176_outgoing_load_factorization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 83, 611-641)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.txt`

## What the paper claims

Stage 176 (anchor MTDC-T9.2, status ExactClosure) factors the outgoing load factor `\(\Lambda_r = P_r/\Delta_r\)` (with `\(P_r = \Omega_{U,r}^2 G_{W,r} + R_r G_{U,r}\)`, `\(\Delta_r = \Omega_{U,r}^2\Omega_{W,r}^2 - R_r^2\)`) into three microscopic slippages. The `\stagefield{Output}` states verbatim: "Factors the outgoing load into mixed-leg, interference, and hybridization slippages; under rigidity, zero defect requires \(G_W/\Omega_W^2\propto\sqrt K\)." The notes enumerate four distinct deliverables: (D1) the exact factorization `\(\Lambda_r^2/K = \mathcal M_r^2 (1+\mathcal I_r)^2/(1-\mathcal H_r)^2\)` with `\(\mathcal M_r = G_{W,r}/(\Omega_{W,r}^2\sqrt K)\)`, `\(\mathcal I_r = R_r G_{U,r}/(\Omega_{U,r}^2 G_{W,r})\)`, `\(\mathcal H_r = R_r^2/(\Omega_{U,r}^2\Omega_{W,r}^2)\)` (§2, appendix eq:app-part05-load-factor-factorization); (D2) the factored first-order defect `\(\Sigma_r^{(N)} = 2\delta\ln\mathcal M_r + \frac{2\mathcal I_r}{1+\mathcal I_r}\delta\ln\mathcal I_r + \frac{2\mathcal H_r}{1-\mathcal H_r}\delta\ln\mathcal H_r\)` with the three primitive `\(\delta\ln\)` formulas of §4 (appendix eq:app-part05-SigmaN-factorized); (D3) the fully expanded primitive-variable transport law of §5; and (D4) the rigidity corollary §6/§7 that setting `\(\delta\ln\mathcal I_r=\delta\ln\mathcal H_r=0\)` collapses `\(\Sigma_r^{(N)}\)` to `\(2\delta\ln\mathcal M_r\)`, hence zero defect requires the square-root mixed-leg law `\(G_W/\Omega_W^2\propto\sqrt K\)`.

## What the script claims to verify

The SymPy docstring lists exactly the four deliverables (factorization; first-order defect decomposition; expanded primitive-variable transport; rigidity corollary → square-root law). The four `expect_zero` checks test: (A1) the algebraic factorization identity `\(\Lambda^2/K - \mathcal M^2(1+\mathcal I)^2/(1-\mathcal H)^2 = 0\)`; (A2) that an independently-computed first-order drift `Sigma_exact` (from an exponential perturbation `X \to X e^{\epsilon\,\delta\ln X}` of the primitive variables, log-ratio, series in `\(\epsilon\)`, coefficient of `\(\epsilon^1\)`) equals the hand-built `Sigma_factored`; (A3) the same `Sigma_exact` equals the expanded primitive-variable form `Sigma_expanded`; and (A4) that `Sigma_factored` with `dlnI, dlnH` set to 0 equals `2*dlnM`. The Mathematica script mirrors all four checks with `expectZero`. Both engines exercise a single one-port case (subscript `r` dropped).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Verdict |
|---|---|---|
| D1 factorization `\(\Lambda^2/K = M^2(1+I)^2/(1-H)^2\)` | A1 (sympy:41-44, math:37-40) | match — non-tautological, independently anchored |
| D2 factored defect `\(\Sigma^{(N)}\)` + §4 `\(\delta\ln\)` formulas | A2 (sympy:62-82, math:60-72) | match — checked vs independent `Sigma_exact` |
| D3 expanded primitive transport (§5) | A3 (sympy:84-92, math:74-83) | match — checked vs independent `Sigma_exact` |
| D4 rigidity corollary `\(\Sigma^{(N)}\to 2\delta\ln M\)` | A4 (sympy:102-103, math:93-94) | partial — assertion is true by construction (tautological); the *substantive* corollary content is only implicitly covered via A2 |

The identity content (D1–D3) is faithfully and non-tautologically verified by both engines; the constants, definitions, and `\(\delta\ln\)` formulas all match the notes and appendix exactly. `paper_alignment` is `aligned` — the only defect (D4) is a verification-construction weakness, not a claim mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-44 | `expect_zero(Lambda**2/K - Mcal**2*(1+Ical)**2/(1-Hcal)**2)` | D1 | yes |
| A2 | sympy | 82 | `expect_zero(Sigma_exact - Sigma_factored)` | D2 | yes |
| A3 | sympy | 92 | `expect_zero(Sigma_exact - Sigma_expanded)` | D3 | yes |
| A4 | sympy | 103 | `expect_zero(Sigma_rigid - 2*dlnM)` where `Sigma_rigid = Sigma_factored.subs({dlnI:0,dlnH:0})` | D4 | no (tautological) |
| A1' | mathematica | 37-40 | `expectZero[lambda^2/k - mCal^2*(1+iCal)^2/(1-hCal)^2]` | D1 | yes |
| A2' | mathematica | 72 | `expectZero[sigmaExact - sigmaFactored]` | D2 | yes |
| A3' | mathematica | 83 | `expectZero[sigmaExact - sigmaExpanded]` | D3 | yes |
| A4' | mathematica | 94 | `expectZero[sigmaRigid - 2*dlnM]` where `sigmaRigid = sigmaFactoredForm /. {dlnI->0,dlnH->0}` | D4 | no (tautological) |

A2/A3 are the load-bearing checks: `Sigma_exact` is derived by differentiating the *actual* `\(\Lambda\)` definition through an independent exponential-perturbation route, while `Sigma_factored`/`Sigma_expanded` are the paper's hand-written decompositions, so equality genuinely validates the §4/§5 weights and `\(\delta\ln\)` formulas at once. A1 is a real algebraic identity (verified by hand below). A4 is the only "no" row.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py:102-103`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl:93-94`

**What's wrong:**
The rigidity-corollary check (paper deliverable D4, notes §6/§7, appendix eq:app-part05-square-root-mixed-leg) is verified against the *constructed* `Sigma_factored`, not against the independent `Sigma_exact`. SymPy:
```python
# line 79-81
Sigma_factored = sp.simplify(2 * dlnM + 2 * Ical / (1 + Ical) * dlnI + 2 * Hcal / (1 - Hcal) * dlnH)
# line 102-103
Sigma_rigid = sp.simplify(Sigma_factored.subs({dlnI: 0, dlnH: 0}))
expect_zero("rigidity reduction to 2 d ln M", Sigma_rigid - 2 * dlnM)
```
Because `Sigma_factored` is *built* as `2*dlnM + (I-weight)*dlnI + (H-weight)*dlnH`, substituting the expressions `dlnI -> 0` and `dlnH -> 0` and asserting the result equals `2*dlnM` is algebraically guaranteed by construction — it is the statement `a + b*0 + c*0 == a`. It cannot fail for any choice of the weights, the `\(\delta\ln\)` formulas, or the underlying `\(\Lambda\)` definition; it contains no physics. The Mathematica A4' (`sigmaRigid = sigmaFactoredForm /. {dlnI->0, dlnH->0}` then `expectZero[sigmaRigid - 2*dlnM]`) has the identical defect. The *substantive* corollary — that the actual first-order drift of `\(\Lambda^2/K\)` reduces to `\(2\delta\ln\mathcal M\)` once `\(\delta\ln\mathcal I=\delta\ln\mathcal H=0\)` — is only implicitly covered because A2 already proved `Sigma_exact == Sigma_factored`; A4 itself adds nothing.

**Why this matters:**
D4 is the stage's headline positive theorem (the square-root mixed-leg law and the `\stagefield{Output}` "under rigidity, zero defect requires \(G_W/\Omega_W^2\propto\sqrt K\)"). A tautological check leaves that deliverable without an independent verification: if a future edit broke the relationship between `Sigma_exact` and the rigidity reduction, A4 would still pass silently.

**Required change:**
Re-derive the rigidity reduction from the *independent* `Sigma_exact` rather than from the constructed `Sigma_factored`. In SymPy, after `Sigma_exact` is computed (line 72), build the rigidity-constrained drift by substituting the rigidity conditions `dlnIExpr = 0`, `dlnHExpr = 0` into `Sigma_exact` and confirm it equals `2*dlnM`. Concretely, solve the two linear rigidity constraints for two of the primitive log-drifts and substitute into `Sigma_exact`, then assert `(that) - 2*dlnM == 0`. See directive for the exact substitution. Mirror the same independent reduction in the `.wl`.

**Verification:**
After the fix, the rigidity check's residual must derive from `Sigma_exact` (not `Sigma_factored`). The verifier confirms: (a) the new check references `Sigma_exact`/`sigmaExact`; (b) both scripts still print the rigidity line `= 0` and exit 0; (c) deliberately perturbing one weight in `Sigma_factored` would no longer leave A4 passing (the check now depends on the true drift).

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl:32-94`
- (compare) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py:36-103`

**What's wrong:**
The `.wl` is a near line-by-line port of the `.py` with identical variable choreography. Corresponding sections:
- factorization: SymPy `Lambda = (OU2*GW + R*GU)/(OU2*OW2 - R**2)` / `Mcal = GW/(sp.sqrt(K)*OW2)` / `Ical = R*GU/(OU2*GW)` / `Hcal = R**2/(OU2*OW2)` (lines 36-39) vs Mathematica `lambda = (ou2*gw + r*gu)/(ou2*ow2 - r^2)` / `mCal = gw/(Sqrt[k]*ow2)` / `iCal = r*gu/(ou2*gw)` / `hCal = r^2/(ou2*ow2)` (lines 32-35) — same names, same order, same forms.
- perturbation: SymPy `Kp = K*sp.exp(eps*dK)` … `Lambdap = (OU2p*GWp + Rp*GUp)/(OU2p*OW2p - Rp**2)` (54-61) vs Mathematica `kP = k*Exp[eps*dK]` … `lambdaP = (ou2P*gwP + rP*guP)/(ou2P*ow2P - rP^2)` (52-59) — identical.
- decomposition coefficients: SymPy `Sigma_factored = ... 2*Ical/(1+Ical)*dlnI + 2*Hcal/(1-Hcal)*dlnH` (79-81) vs Mathematica `sigmaFactoredForm = 2*dlnM + 2*iCal*dlnI/(1+iCal) + 2*hCal*dlnH/(1-hCal)` (70) — identical.

The one genuine difference is the first-order extraction: SymPy uses `sp.series(...).removeO()/eps` (Taylor truncation) while Mathematica uses `D[Log[...], eps] /. eps -> 0` (symbolic differentiation). That single mechanistic difference is the only thing keeping this from being a verbatim transliteration. Per the second-engine policy the structural port is still flagged; severity is low because (i) the deliverables are pure algebraic identities for which there is essentially one correct symbolic form, and (ii) the `\(\epsilon^1\)`-coefficient extraction does differ between engines.

**Why this matters:**
A line-by-line port provides weaker cross-engine assurance than an independent re-derivation: a transcription error in the shared algebra would be reproduced in both engines. For algebraic identities the risk is modest, hence low severity, but the policy requires the flag.

**Required change:**
No structural rewrite is mandated for a low-severity transliteration of pure-algebra identities; the minimum acceptable remediation is to make the Mathematica first-order section derive `sigmaExact` by a route that does not mirror the SymPy variable-by-variable construction — it already uses `D[Log[...], eps]`, which differs from SymPy's `series`. Document this divergence (a one-line comment) so the independence is explicit, and keep the distinct extraction. If the orchestrator's policy treats any structural port as blocking, escalate; otherwise this is informational. Do not introduce new physics.

**Verification:**
Confirm the `.wl` retains a first-order extraction (`D[Log[...], eps]`) distinct from the SymPy `series` route and that a clarifying comment notes the engines use different mechanisms. Both scripts still exit 0.

### F3 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py:32` (banner literal "STAGE 159")
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl:26` (banner literal "STAGE 159")
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.txt:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.txt:11`

**What's wrong:**
Both scripts' top banner reads `banner("STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION")` (sympy:32, math:26) and both saved transcripts echo "STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION", but this is Stage 176. The docstring (sympy:5) and the Mathematica tail (`Print["Stage 176 Mathematica audit passed."]`, math:104) correctly say 176; only the banner string is wrong. This is a stale/copy-pasted label, not a math error.

The saved `.txt` mtimes (sympy 12:47, math 13:22) are newer than the script mtimes (11:56), so the outputs are otherwise fresh and content-faithful; the only mismatch they carry is the same mislabeled banner the script emits. I file this under `stale_output` (the closest of the ten categories) as informational; it is not blocking and does not change any verified result.

**Why this matters:**
A wrong stage number in the banner/transcript hurts traceability when the transcript is read in isolation and could mask a genuinely misfiled script in a later audit. Pure cosmetics — no numeric result is affected.

**Required change:**
Change the banner literal in both scripts from `STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION` to `STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION` (sympy:32, math:26). Re-run both scripts so the transcripts pick up the corrected banner.

**Verification:**
Both scripts print `STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION`; both transcripts' line 11 reads `STAGE 176 …`; both exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is a structural port of the `.py` (see F2 for the corresponding-section quotes): identical variable names (`lambda`/`mCal`/`iCal`/`hCal`), identical exponential-perturbation choreography, identical decomposition coefficients. The lone genuine divergence is the first-order extraction — SymPy `sp.series(log_ratio, eps, 0, 2).removeO()/eps` vs Mathematica `D[Log[log_ratio], eps] /. eps -> 0`. This is enough to keep the engines from echoing each other's `\(\epsilon\)`-truncation algebra but not enough to make the rest an independent re-derivation. Flagged as low-severity `mathematica_transliteration` (F2).

## Engine cross-check

Both engines produce identical residuals and pass states. SymPy transcript (lines 13, 21-22, 33): `exact factorization of Lambda^2/K = 0`; `factored first-order defect formula = 0`; `expanded primitive-variable transport = 0`; `rigidity reduction to 2 d ln M = 0`; `EXIT_CODE: 0`. Mathematica transcript (lines 14, 23, 25, 37): same four residuals `= 0` each followed by `PASS:`; `EXIT_CODE: 0`. The two engines agree exactly on all four checks. (The agreement on A4/A4' is unsurprising since both are tautological — see F1.)

## Verdict justification

Verdict is `findings`. The core algebra is sound and paper-aligned: I verified the factorization identity D1 by hand (`1+I = P/(Ω_U² G_W)`, `1-H = Δ/(Ω_U² Ω_W²)`, `M² = G_W²/(K Ω_W⁴)`, giving `RHS = P²/(K Δ²) = Λ²/K`), and the load-bearing first-order checks A2/A3 genuinely test the paper's §4/§5 decomposition against an *independently* differentiated `Sigma_exact`, so D2 and D3 hold up under attack. I tried to break A1 (non-tautological), the symbol assumptions (`positive=True` does not mask any branch in these pure identities or in the `\(\epsilon^1\)`-coefficient extraction, since the log branch offset drops under differentiation), and the `½δ_K` factor (matches `\(\ln M = \ln G_W - \tfrac12\ln K - \ln\Omega_W^2\)`) — all held. The genuine defects: (F1, medium) the rigidity-corollary check A4/A4' is tautological — it substitutes into the *constructed* `Sigma_factored` rather than the independent `Sigma_exact`, leaving the headline square-root-law deliverable D4 without a falsifiable check; (F2, low) the `.wl` is a structural transliteration of the `.py`; (F3, low) both banners mislabel the stage as "159". No stop-cold condition: the math is fully reconcilable, nothing downstream is invalidated by the fixes (D4's true content is already established by A2, so the fix strengthens verification without changing any result). No `paper_misalignment`: every paper deliverable maps to a script check with matching constants and forms.

## Self-test notes

Trap 1 (variable independence): no `sp.diff`/`D` against a variable the expression lacks — the only differentiation is `D[Log[...], eps]` / `series(..., eps)`, and the perturbed quantities all depend on `eps` by construction, so the `\(\epsilon^1\)` coefficient is genuinely nonzero. For the F1 fix I confirmed `Sigma_exact` depends on all six primitive log-drifts, so substituting the two rigidity constraints leaves a nonzero `2*dlnM` residual to test against. Trap 3 (trivial-case): the F1-prescribed substitution of `dlnIExpr=0, dlnHExpr=0` into `Sigma_exact` reduces to `2*(dGW - dOW - dK/2)` = `2*dlnM`, a nonzero literal, so the corrected check is non-vacuous; I confirmed the rigidity conditions zero exactly the I- and H-channels and not the M-channel. Trap 5 (paper round-trip): the F1 fix introduces no new constant — it re-uses the existing `dlnMExpr`/`dlnM` and the rigidity conditions exactly as the notes §6 state them, so no new `paper_misalignment` is created. Traps 2 (parity) and 4 (path) are not applicable (no unbounded integrals; no missing-script findings).
