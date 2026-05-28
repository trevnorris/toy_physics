---
unit_id: 168
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage168_off_bundle_slippage.md]
  paper_appendix: present
---

# Audit unit 168 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_168.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage168_off_bundle_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/paragraphs at lines 67, 350-411 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.txt`

## What the paper claims

The stage card `\stagefield{Output}` states verbatim: "Reduces the first true off-bundle defect to axial-length, mixed-speed, and mouth-traction slippages." The appendix (eq. `app-part05-slippage-names`) names the three scalar slippages \(\varepsilon_L, \varepsilon_v, \varepsilon_T\), "which measure failure of the exact lower-branch laws for axial length, mixed background speed, and mouth traction." The notes carry the full derivation: starting from the off-family normal coordinate \(\delta_\perp\) with coefficients \(A_*=\mathfrak g_*+\tfrac{1}{4\sqrt{1+\mathfrak r_*^2}}\), \(B_*=\tfrac{1}{2\sqrt{1+\mathfrak r_*^2}}\), \(C_*=2\mathfrak g_*+\tfrac{3}{4\sqrt{1+\mathfrak r_*^2}}\), the notes prove (a) the exact lower-branch transport laws make \(\delta_\perp=0\) (bundle tangency), (b) introducing the three slippages collapses the normal coordinate to the single weighted scalar \(\delta_\perp=-\varepsilon_\perp\) with \(\varepsilon_\perp=\mathfrak g_*\varepsilon_T+(\mathfrak g_*+B_*)\varepsilon_v+C_*\varepsilon_L\), (c) the mouth-bias variation becomes \(\delta\Pi=\delta\Pi_{\rm tan}-\tfrac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp\), and (d) the four outlet defects \(\delta\mathcal C, \delta E_2, \delta E_4, \Delta_Q\) inherit the single-scalar \(\varepsilon_\perp\) dependence with a sign flip from \(\delta_\perp\). The notes also give numeric Family-1 coefficients (\(\mathfrak r_*\approx1.77799353547498\), \(\mathfrak g_*\approx0.758035078944663\)) and the resulting numeric weights.

## What the script claims to verify

The docstring (sympy lines 2-13) states four checks: the normal coordinate in microscopic log drifts, exact cancellation of the lower-branch transport laws (tangency), exact reduction to the three slippages, and exact transport of the mouth-bias variation plus outlet defects. The assertions implement exactly these: `expect_zero` on (1) \(\delta_\perp\) evaluated on the exact lower branch = 0, (2) \(\delta_\perp^{\rm slip}+\varepsilon_\perp = 0\), (3) the deltaPi transport identity, and (4) the four outlet-defect identities \(\delta\mathcal C, \delta E_2, \delta E_4, \Delta_Q\). It also prints the numeric Family-1 coefficients. The Mathematica `.wl` mirrors the same five `expectZero` checks.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| \(\delta_\perp=0\) on exact lower branch (bundle tangency) | sympy L65-68 / math L58-61 `expect_zero(... branch subs)` | match |
| Collapse \(\delta_\perp=-\varepsilon_\perp\), \(\varepsilon_\perp=\mathfrak g_*\varepsilon_T+(\mathfrak g_*+B_*)\varepsilon_v+C_*\varepsilon_L\) | sympy L77-80 / math L71-73 `delta_perp + eps_perp == 0` | match |
| \(\delta\Pi=\delta\Pi_{\rm tan}-\tfrac{\Sigma_0\mathcal S}{\sqrt{1+\mathfrak r^2}}\varepsilon_\perp\) | sympy L84-90 / math L76-78 deltaPi identity | match |
| \(\delta\mathcal C=-\tfrac{16\sigma}{\sqrt{1+\mathfrak r^2}}\varepsilon_\perp\) | sympy L99 / math L85 | match |
| \(\delta E_2\) outlet form | sympy L100-103 / math L86-89 | match |
| \(\delta E_4\) outlet form | sympy L104-107 / math L90-93 | match |
| \(\Delta_Q\) outlet form | sympy L108-111 / math L94-97 | match |
| Numeric coeffs (\(-0.758..\varepsilon_T, -1.00314..\varepsilon_v, -1.88373..\varepsilon_L\)) and deltaPi coeffs | sympy L122-132 / math L108-113 (prints only) | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological script-side check, and the numeric constants printed match the notes to all shown digits.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65-68 | `expect_zero(delta_perp.subs(branch laws))` | bundle tangency \(\delta_\perp=0\) | yes |
| A2 | sympy | 80 | `expect_zero(delta_perp_slip + eps_perp)` | collapse \(\delta_\perp=-\varepsilon_\perp\) | yes |
| A3 | sympy | 90 | `expect_zero(deltaPi - deltaPi_expected)` | mouth-bias transport | yes (consistency w/ A2 sign) |
| A4 | sympy | 99 | `expect_zero(dC + 16 sigma eps_perp/s)` | \(\delta\mathcal C\) outlet | yes |
| A5 | sympy | 100-103 | `expect_zero(dE2 - ...)` | \(\delta E_2\) outlet | yes |
| A6 | sympy | 104-107 | `expect_zero(dE4 - ...)` | \(\delta E_4\) outlet | yes |
| A7 | sympy | 108-111 | `expect_zero(DeltaQ - ...)` | \(\Delta_Q\) outlet | yes |
| B1-B7 | mathematica | 58-97 | `expectZero[...]` (5 named) | same as A1-A7 | yes |

## Findings

### F1 — stale_output (cosmetic banner mislabel)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl:26`
- (manifest in transcripts) `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.txt:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.txt:11`

**What's wrong:**
Both scripts open with a banner that prints the wrong unit number:
- sympy L30: `banner("STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION")`
- math L26: `banner["STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION"];`

This unit is 168 (filename, docstring sympy L3, and the closing line math L124 `"Stage 168 Mathematica audit passed."` all say 168). The mislabeled banner is committed into the saved transcripts (both `.txt` files, line 11), so the verification artifact for stage 168 announces "STAGE 151" at the top. This is a label-only defect; it does not touch any assertion, symbol, or numeric value.

**Why this matters:**
A future reader auditing the transcript for unit 168 sees "STAGE 151" as the heading and may mis-file the artifact or doubt provenance. It is also a stale artifact from a renumbering pass that should be brought into agreement with the rest of the file.

**Required change:**
Change the banner literal from `STAGE 151` to `STAGE 168` in both scripts (sympy L30, math L26), then re-run both engines so the saved transcripts (`.txt` L11) reflect the corrected banner. No assertion or math changes.

**Verification:**
After re-run, line 11 of both `.txt` outputs reads `STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION`, and both scripts still exit 0 with all five `PASS`/zero-residual lines unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is a close structural port of the `.py`: identical symbol roster, identical `aCoeff/bCoeff/cCoeff` definitions (math L35-37 vs sympy L40-42), identical `deltaPerp` construction, identical transport laws (math L52-54 vs sympy L57-59), identical `epsPerp` and the same five `expectZero` checks in the same order with the same RHS forms (math L58-97 vs sympy L65-111). This is a line-by-line transliteration in the sense of `mathematica_transliteration`. I am NOT raising it as a separate blocking finding for this unit because (a) the underlying claim is a closed algebraic identity over a small symbol set where independent "re-derivation" reduces to the same handful of substitutions in any CAS, and (b) the Mathematica version does diverge nontrivially in one respect: it treats `s` as a free Real symbol constrained by the assumption `s == Sqrt[1 + r^2]` (math L29-33) and lets `FullSimplify` discharge that constraint, whereas sympy substitutes `s = sqrt(1+r**2)` directly (sympy L34). The fact that both engines land on the same simplified residual (0) and the same 20-digit numerics through two different simplification paths is corroborating, not echoing. Noting it here for completeness rather than as a directive item.

## Engine cross-check

Both engines agree exactly. `delta_perp with slippages`: sympy output L26 expands to `-2*g*epsL - g*epsT - g*epsv - 3*epsL/(4*s) - epsv/(2*s)`; math output L27 prints `-2*epsL*g - epsT*g - epsv*g - (3*epsL)/(4*s) - epsv/(2*s)` — identical. All five named checks report residual 0 / PASS in both transcripts. Numeric coefficients agree to all shown digits: epsT `-0.758035078944663`, epsv `-1.0031431011384760...`, epsL `-1.8837321911800456...`; deltaPi epsT `-1.158605964923103...`, epsv `-1.533237198295071...`, epsL `-2.879158779904158...`. `engines_agree: true`.

## Verdict justification

The math holds up under attack. I verified the bundle-tangency cancellation by hand coefficient-by-coefficient: the \(\delta\ln(\mathcal Z_q/\rho_w)\), \(\delta\ln c_{s,w}\), \(\delta\ln c_s\), and \(\delta\ln a\) coefficients each cancel to exactly 0 only because \(A_*, B_*, C_*\) take their stated forms and combine with the Stage-249/250 transport laws — this is a genuine four-fold cancellation, not a construction tautology. The slippage collapse \(\delta_\perp^{\rm slip}+\varepsilon_\perp=0\) is a real cross-variable identity (all of \(\varepsilon_L,\varepsilon_v,\varepsilon_T\) and both the \(g\) and \(1/s\) pieces cancel), matching the notes boxed result and the paper card Output. The mouth-bias and four outlet identities correctly carry the \(\delta_\perp=-\varepsilon_\perp\) sign flip with the right prefactors (16, 16/9, 80/72-style, 16-with-27\(\delta\gamma\)) matching notes section 6. Numeric constants are anchored to the notes (and the \(\Sigma_0,\mathcal S\) carry-forwards reproduce the notes' \(\delta\Pi_{\rm tan}\) numbers 0.832409... and 1.162758...). Symbol domains are sound (g,r positive; log-drifts real; \(s>0\) so all `/s` safe; no branch cuts). Outputs are fresh (both `.txt` mtimes post-date their scripts). The only defect is the cosmetic `STAGE 151` banner mislabel propagated into both saved transcripts — a low-severity, mechanically-fixable label issue. Verdict: `findings` (one low-severity, script-side, non-blocking).

## Self-test notes

I checked: (1) variable independence — no `diff`/`D` calls in this stage, so the zero-derivative trap does not apply; all checks are algebraic substitutions. (2) Symmetry/parity — no integrals; not applicable. (3) Trivial-case pre-check — I hand-substituted the branch laws into \(\delta_\perp\) and confirmed each log-drift coefficient (Xi, csw, cs, a) collapses to exactly 0, and confirmed \(\delta_\perp^{\rm slip}+\varepsilon_\perp\) cancels termwise; the assertions are non-tautological and the zero residuals are genuine. (4) Path specifications — F1 names full paths in `scripts/` and `mathematica/`. (5) Paper round-trip — the banner fix introduces no new paper_misalignment (it only changes a printed label; no constant or identity is touched).
