---
unit_id: 115
batch: IV.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage115_core_balance_compensation.md]
  paper_appendix: present
---

# Audit unit 115 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_115.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage115_core_balance_compensation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (relevant subsection lines 511-573)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.txt`

## What the paper claims

The stage proves that the concrete two-channel core model lands exactly on the nontrivial compensated canonical branch. The paper card body equation states the coupling balance "`g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2` realizes the compensated branch". The notes expand this into multiple deliverables: (i) the boxed canonical conditions `rho_c=4 sigma_c, kappa_c=1/3, gamma_c=1/9`; (ii) the boxed coupling-balance law identical to the card; (iii) the boxed bare-channel normalization `kappa_0=(1+r_c)/3, gamma_0=(1+r_c)/9`; (iv) the boxed exact collapse `delta Lambda_core(z) = 4 sigma_* - sigma_*/(1 - z^2/3 - i z^5/9)` with `sigma_* = g_s^2/(4 K_s)`; (v) the normalized outgoing fingerprint `hat Y_2(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)`. The Part-IV appendix (lines 511-573) restates (i)-(ii) and adds a parent-action overlap track via `mathfrak r = lambda/sqrt(K_sK_q)`, `mathfrak g = g_q sqrt(K_s)/(g_s sqrt(K_q))`, with the equivalent parent compensation family `1 + mathfrak r^2 = 4(mathfrak g - mathfrak r)^2`, `mathfrak g_+- = mathfrak r +- (1/2) sqrt(1 + mathfrak r^2)`.

## What the script claims to verify

Both scripts (SymPy and Mathematica) assert the same three identities: (A1) on the balance surface defined by solving `rho_c - 4 sigma_c = 0` for `g_q`, the substituted `sigma_c` equals `sigma_* = g_s^2/(4 K_s)`; (A2) the constructed `delta_core = rho_c - sigma_c / (1 - (kappa_0/(1+r_c)) z^2 - i (gamma_0/(1+r_c)) z^5)` with `g_q -> gq_branch`, `kappa_0 -> (1+r_c)/3`, `gamma_0 -> (1+r_c)/9` is identically equal to the target `4 sigma_* - sigma_*/(1 - z^2/3 - i z^5/9)`; (A3) the normalized response `Y_eff = Lambda_eff(0)/Lambda_eff(z)` (with `Lambda_eff = Lambda_out + target_delta`, `Lambda_out = -3 + z^2/3 + z^4/9 + i z^5/9`) agrees with the target fingerprint `1 + z^2/9 + 4 z^4/81 + i z^5/27` through order z^5. Both scripts also print `sigma_*`, `kappa_0`, `gamma_0` for visual inspection.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (i) `rho_c = 4 sigma_c`, `kappa_c=1/3`, `gamma_c=1/9` canonical conditions | Implicit: A2 substitutes `kappa_0/(1+r_c) = 1/3` and `gamma_0/(1+r_c) = 1/9`; A1 enforces `rho_c = 4 sigma_c` via `solve(balance_eq=0, g_q)` | match |
| (ii) coupling balance `g_s^2(K_sK_q+lambda^2) = 4(K_sg_q - lambda g_s)^2` | Encoded in `balance_eq = rho_c - 4 sigma_c` (sympy line 28; math line 37); printed as Factor[]/expanded form | match (printed, not asserted, but A1 exercises it via the substitution) |
| (iii) `kappa_0 = (1+r_c)/3`, `gamma_0 = (1+r_c)/9` | A2 substitutes these explicit values (`kappa0_can`, `gamma0_can`) before the simplify-to-zero check | match |
| (iv) exact collapse `delta Lambda_core = 4 sigma_* - sigma_*/(1 - z^2/3 - i z^5/9)`, `sigma_* = g_s^2/(4K_s)` | A1 + A2 together | match |
| (v) `hat Y_2(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)` | A3 series check | match |
| Appendix `mathfrak r, mathfrak g` parent compensation family | not exercised — only the unscaled form of the balance equation is solved | extra-paper (parent reparametrization could provide an independent algebraic route for Mathematica) |

Set `paper_alignment: aligned` — every paper/notes deliverable has a corresponding script check, and the scripts do not test extra content the paper doesn't claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-44 | `expect_zero("sigma_c on balance surface", sigma_c.subs(g_q, gq_branch) - sigma_star)` | (i) `rho_c=4 sigma_c` ⇒ `sigma_c = sigma_*`; (iv) sigma_* value | yes |
| A2 | sympy | 49-61 | `expect_zero("exact collapse to compensated branch", delta_core - target_delta)` | (i) kappa_c=1/3, gamma_c=1/9; (iii) bare kappa_0, gamma_0; (iv) exact collapse | yes |
| A3 | sympy | 66-70 | `expect_zero("normalized outgoing fingerprint preserved", series(Y_eff,z,0,6) - Y_target)` | (v) Y_2 fingerprint | yes (but Lambda_eff factors as (1-sigma_*) * canonical-form, so the (1-sigma_*) factor cancels in Y_eff = Lambda_eff(0)/Lambda_eff — the check verifies a relation independent of sigma_*; this is still a non-trivial check of the canonical-branch shape, just narrower than its label suggests) |
| B1 | mathematica | 46 | `expectZero["sigma_c on balance surface", (sigmaC /. gQ -> gQBranch) - sigmaStar]` | (i)+(iv) | yes |
| B2 | mathematica | 51-57 | `expectZero["exact collapse to compensated branch", deltaCore - targetDelta]` | (i)+(iii)+(iv) | yes |
| B3 | mathematica | 60-63 | `expectZero["normalized outgoing fingerprint preserved", series(yEff,z,0,5) - yTarget]` | (v) | yes (same caveat as A3) |

No tautological assertions; all three engine pairs exit 0 with residual 0 in the saved outputs. The scripts also implicitly require `Length[gQSolutions] === 2` in Mathematica (line 41), which is a valid non-trivial structural check on the balance-equation branches.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl:33-63`
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.py:24-70`

**What's wrong:**
The `.wl` script is a near-mechanical transliteration of the `.py` script. Corresponding sections, with names rewritten in CamelCase but no structural divergence:

- sympy 24-26 (`r_c = lam**2/(K_s*K_q)`, `rho_c = g_s**2/K_s`, `sigma_c = (K_s*g_q - lam*g_s)**2/(K_s**2*K_q*(1+r_c))`)
  ↔ math 33-35 (`rC = FullSimplify[lam^2/(kS*kQ), …]`, `rhoC = FullSimplify[gS^2/kS, …]`, `sigmaC = FullSimplify[(kS*gQ - lam*gS)^2/(kS^2*kQ*(1 + rC)), …]`)
- sympy 28-32 (`balance_eq = sp.expand(rho_c - 4*sigma_c)`; print; `gq_solutions = sp.solve(sp.Eq(balance_eq,0), g_q)`)
  ↔ math 37-42 (`balanceEq = Expand[rhoC - 4*sigmaC]`; Print; `gQSolutions = Solve[balanceEq == 0, gQ, Reals]`)
- sympy 38-44 (`gq_branch = sp.simplify(gq_solutions[0])`; `sigma_star = (g_s**2)/(4*K_s)`; expect_zero on sigma_c)
  ↔ math 44-46 (`gQBranch = FullSimplify[gQ /. First[gQSolutions], …]`; `sigmaStar = FullSimplify[gS^2/(4*kS), …]`; expectZero on sigma_c)
- sympy 46-61 (kappa0_can, gamma0_can, delta_core, target_delta, expect_zero on collapse)
  ↔ math 48-57 (kappa0Can, gamma0Can, deltaCore, targetDelta, expectZero on collapse)
- sympy 63-70 (Lambda_out, Lambda_eff via series, Y_eff = Lambda_eff(0)/Lambda_eff, expect_zero on fingerprint)
  ↔ math 59-63 (lambdaOut, lambdaEff via Series, yEff = lambdaEff(z=0)/lambdaEff, expectZero on fingerprint)

Every algebraic step, the choice of branch (`gq_solutions[0]` ↔ `First[gQSolutions]`), the series truncation order (5/6), the explicit substitution sequence, and the final identity targets are line-for-line identical. The `.wl` script therefore does not provide an independent re-derivation of the result from the physical premises — it simply rewrites the SymPy algebra in Mathematica notation. The notes file (and the Part-IV appendix lines 530-543) provide a genuine independent algebraic route via the parent reparametrization `mathfrak r = lambda/sqrt(K_s K_q)`, `mathfrak g = g_q sqrt(K_s)/(g_s sqrt(K_q))`, yielding the compensation family `1 + mathfrak r^2 = 4(mathfrak g - mathfrak r)^2` and `mathfrak g_+- = mathfrak r +- (1/2) sqrt(1 + mathfrak r^2)`. Reformulating the balance check through `mathfrak r`/`mathfrak g` (and showing equivalence to the `g_q`-solve route) would supply genuine engine independence.

**Why this matters:**
The second-engine policy exists to catch SymPy-specific bugs (e.g., a bad `simplify` branch, a wrong canonical form, an unnoticed branch choice). A transliterated `.wl` script will reproduce any such SymPy bug. With the present `.wl` script, a SymPy `simplify` that happened to collapse a non-vanishing residual to zero under unjustified positivity assumptions on `K_s, K_q, lambda` would also be returned as zero by the Mathematica script, because the Mathematica script uses the same assumption set (`kS > 0 && kQ > 0`) over the same intermediate forms. The two outputs agreeing is therefore much weaker evidence than independent agreement would be.

**Required change:**
Add (alongside the existing `g_q`-solve derivation, not replacing it) an independent block in the `.wl` script that derives the balance surface and the σ_* identity through the parent-overlap reparametrization `mathfrak r, mathfrak g` of equations `eq:app-part04-r-g-parent-ratios` and `eq:app-part04-parent-compensation-family` (Part-IV appendix lines 530-543). Specifically, after line 46 (the existing sigma_c-on-balance check) and before line 48 (`kappa0Can = …`), insert a block that:

1. defines `frakR = lam/Sqrt[kS*kQ]` and `frakG = gQ*Sqrt[kS]/(gS*Sqrt[kQ])`;
2. forms `parentFamilyResidual = 1 + frakR^2 - 4*(frakG - frakR)^2` (the appendix family condition);
3. exhibits `Solve[parentFamilyResidual == 0, frakG]` returning the two roots `frakR +- (1/2)*Sqrt[1 + frakR^2]` and asserts each root, after substituting back to `g_q = frakG * gS * Sqrt[kQ]/Sqrt[kS]`, satisfies `balanceEq == 0` (i.e. the parent reparametrization is *equivalent* to the original solve route — `expectZero` on the substitution residual, NOT on the literal `Solve` output);
4. independently shows that under `parentFamilyResidual == 0` the ratio `sigma_c / sigma_*` equals 1, by substituting `g_q = frakG * gS * Sqrt[kQ]/Sqrt[kS]` with `frakG -> frakR - (1/2) Sqrt[1 + frakR^2]` into `sigma_c` and `FullSimplify`ing the difference from `sigma_*`.

The load-bearing new check is the final `expectZero["independent: sigma_c = sigma_* via parent reparametrization", ...]`. Leave the existing `g_q`-solve assertions in place — the two derivations are now complementary engine paths, not echoes.

**Verification:**
After Codex applies, the verifier will run `redteam exec-mathematica 115` and confirm:
- the new `independent: sigma_c = sigma_* via parent reparametrization` line appears in the transcript with residual `0`;
- the script exits 0;
- the existing three `expectZero` lines still pass with residual `0`.

## Independent-derivation check (Mathematica)

The `.wl` script is a transliteration of the `.py` script — see F1 for the line-by-line correspondence. The Mathematica `Solve[..., gQ, Reals]` matches SymPy `sp.solve(..., g_q)`, the `First[gQSolutions]` branch pick matches `gq_solutions[0]`, and the subsequent `FullSimplify` calls take the same intermediate forms. No independent reformulation (e.g. via `mathfrak r, mathfrak g`) is present. The appendix supplies a ready-made independent route that the script does not use.

## Engine cross-check

Both engines agree:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| sigma_c on balance surface − sigma_* | `0` (line 31 of `.txt`) | `0` PASS (lines 15-16) |
| delta_core − target_delta | `0` (line 32) | `0` PASS (lines 17-18) |
| series(Y_eff) − Y_target | `0` (line 33) | `0` PASS (lines 19-20) |
| sigma_* | `g_s**2/(4*K_s)` | `gS^2/(4*kS)` |
| kappa_0 | `(K_q*K_s + lam**2)/(3*K_q*K_s)` | `(1 + lam^2/(kQ*kS))/3` (algebraically equivalent) |
| gamma_0 | `(K_q*K_s + lam**2)/(9*K_q*K_s)` | `(1 + lam^2/(kQ*kS))/9` (algebraically equivalent) |

Outputs are not stale: sympy script mtime 2026-04-01, sympy `.txt` mtime 2026-05-11; mathematica script mtime 2026-05-11 11:56, mathematica `.txt` mtime 2026-05-11 13:08.

## Verdict justification

The audit verdict is `findings` with a single `medium`-severity `mathematica_transliteration` issue. Within the scope of what each script verifies, the assertions are non-tautological, exercise the load-bearing identities of the paper card and notes (coupling balance, exact collapse to the canonical δΛ form, normalized Y_2 fingerprint through z^5), and produce zero residuals in both engines. The scripts faithfully cover every paper-side deliverable in the notes and appendix subsection. The single defect is that the Mathematica script does not derive the result independently from the physical premises — it is a notation rewrite of the SymPy algebra and reproduces SymPy's branch choices and intermediate forms. The appendix's `mathfrak r, mathfrak g` parent reparametrization supplies a genuine independent algebraic route; once added, the second-engine policy is honored. No `paper_misalignment` was found; the paper card's body equation and the notes-level deliverables all map to script-side checks. No `CRITICAL_DOWNSTREAM` flag is warranted — the math is correct, only the second-engine independence is at issue.

Attempted attacks that failed:
- substituted `gq_branch = g_s(2 lam - sqrt(K_qK_s + lam^2))/(2 K_s)` back into `(K_s g_q - lam g_s)^2 / (K_s^2 K_q (1 + r_c))` by hand and recovered `g_s^2 / (4 K_s) = sigma_*` exactly — A1 is non-tautological and correct;
- expanded `delta_core - target_delta` symbolically against the balance solution and found exact cancellation (the 1/(1 - z^2/3 - i z^5/9) denominators match by construction once kappa_0/(1+r_c) = 1/3 and gamma_0/(1+r_c) = 1/9 substitute through) — A2 is correct;
- factored `Lambda_eff = (1 - sigma_*) * (-3 + z^2/3 + z^4/9 + i z^5/9)` so `Y_eff = -3 / (-3 + z^2/3 + z^4/9 + i z^5/9) = 1/(1 - z^2/9 - z^4/27 - i z^5/27)`, expanded to z^5 and recovered `1 + z^2/9 + 4 z^4/81 + i z^5/27` — A3 is correct (but its `(1 - sigma_*)` cancellation makes it weaker than its label suggests; this is noted in the assertion inventory but is not itself a finding because the check still pins down the canonical-branch shape);
- checked z's `positive=True` declaration in sympy for risk on series/simplify steps; the zero checks are on polynomial/rational expressions whose vanishing is independent of z's domain (no branch cuts, no log/sqrt(z) involved), so this is benign — not flagged.

## Self-test notes

- Variable independence: not applicable (no `diff` calls in either script).
- Symmetry/parity: not applicable (no integrals).
- Trivial-case pre-check: substituted gq_branch by hand for A1 and recovered sigma_* analytically; factored (1 - sigma_*) out of Lambda_eff for A3 and verified the residual series matches Y_target. Both pass on paper.
- Path specifications: F1 targets the existing `.wl` (insertion site between lines 46 and 48); no missing-script issue.
- Paper round-trip: F1's prescribed insertion uses the parent-overlap reparametrization that is already stated in appendix lines 530-543, so it does not introduce a new paper-side claim. The mfrak r and mfrak g definitions and the family equation `1 + mfrak r^2 = 4(mfrak g - mfrak r)^2` are paper-anchored.
