---
unit_id: 182
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage182_microscopic_coherent_slippage.md]
  paper_appendix: present
---

# Audit unit 182 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_182.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage182_microscopic_coherent_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 95, 768-817 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` states: "Reduces the defect to \((\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta)\) plus dressing slippage \(\Sigma_\eta\), and isolates \(\Sigma_{\rm tr}\)." The notes and appendix expand this into distinct deliverables: (1) the exact four-slippage grouped-defect law for `\Xi_1` (appendix eq:app-part05-Xi1-slippage-law / notes §4 boxed eq); (2) the selected-branch demand law `\mathcal R_1` carrying one extra dressing slippage `\Sigma_\eta` (notes §5); (3) the tracking-factor drift `\Theta_1` factorizing through the single combination `\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi` (appendix eq:app-part05-Sigma-tr-def and eq:app-part05-Theta-Sigma-tr / notes §6); (4) the exact tracking/nontracking split of `\Xi_1` (notes §6, lines 400-412); (5) the Stage-30 physical-branch drift definitions (`\zeta_Z,\omega_W,\chi_1,\eta_1,\varepsilon_W,\delta_{U,1}`, notes §2.1); and (6) microscopic support-blindness — the support drifts `\delta\ln\lambda_\phi,\delta\ln K_\phi` do not enter the defect (notes §2.2, §8 item 1; appendix line 764 states `\partial_\zeta\Xi_1=0`). Crucially, notes §2.2 delegates the support-blindness mechanism to Stage 249 ("Stage 249 already proved that `\zeta` drops out … so the underlying support drifts do not enter the defect at all"), so at the Stage-182 level support-blindness is a carry-forward corollary, not an independent re-derivation.

## What the script claims to verify

The docstring (sympy lines 6-13) enumerates six checks mirroring the six paper deliverables: the Stage-181/30 physical variables, the four-slippage law, the `\mathcal R_1` dressing slippage, the `\Theta_1` factorization through `\Sigma_{\rm tr}`, the tracking/nontracking split, and "Support microscopic drifts do not appear." The load-bearing assertions are `expect_zero`/`expectZero` calls comparing a microscopic-log "direct" form (built from `lam1,c1,gam1,kU,keta,kW,mu1,tau1`) against the abstract-Sigma "slippage" form with the Sigma symbols substituted back to their microscopic definitions. Four of these checks (lines 93, 100, 114, 130) are substantive and faithfully match the paper's identities. The six "physical branch drifts" checks (lines 74-79) are restatements of the Stage-30 definitions, and the six "support-blindness" checks (lines 139-144) differentiate the abstract-Sigma forms with respect to two declared-but-never-wired support symbols.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) four-slippage `\Xi_1` law | sympy L93 / math L69 `Xi_1 direct - slippage form` | match |
| (2) `\mathcal R_1` + extra `\Sigma_\eta` | sympy L100 / math L75 `R_1 direct - slippage form` | match |
| (3) `\Theta_1` factorizes via `\Sigma_{\rm tr}` | sympy L114 / math L90 `Theta_1 factorization` | match |
| (4) tracking/nontracking split | sympy L130 / math L102 `Xi_1 split - slippage form` | match |
| (5) Stage-30 drift definitions | sympy L74-79 / math L54-59 | match but tautological (see F3) |
| (6) support-blindness | sympy L139-144 / math L107-112 `d.../dlamphi1`, `d.../dkphi` | partial — vacuous (see F1) |

Deliverables (1)-(4) are exercised faithfully and non-tautologically; the identities, constants (the `11`, `9`, factor-2 coefficients), and parameter forms all match the appendix/notes verbatim. Deliverable (6) is only nominally covered: the checks pass trivially because the support symbols are never present in the differentiated expressions. `paper_alignment: partial` reflects (6).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74-79 | `simplify(def - def) == 0` (×6) | claim 5 | no (tautological) |
| A2 | sympy | 93 | `simplify(Xi1_direct - Xi1_slip.subs) == 0` | claim 1 | yes |
| A3 | sympy | 100 | `simplify(R1_direct - R1_slip.subs) == 0` | claim 2 | yes |
| A4 | sympy | 114 | `simplify(Theta1_direct - Theta1_fact.subs) == 0` | claim 3 | yes |
| A5 | sympy | 130-133 | `simplify(Xi1_split.subs - Xi1_slip) == 0` | claim 4 | yes |
| A6 | sympy | 139-144 | `diff(Xi1_slip/R1_slip/Theta1_fact, lamphi1/kphi) == 0` (×6) | claim 6 | no (vacuous) |
| B1 | math | 54-59 | `FullSimplify[def - def] === 0` (×6) | claim 5 | no (tautological) |
| B2 | math | 69 | `FullSimplify[xi1Direct - (xi1Slip /. slipSubs)] === 0` | claim 1 | yes |
| B3 | math | 75 | `FullSimplify[r1Direct - (r1Slip /. slipSubs)] === 0` | claim 2 | yes |
| B4 | math | 90 | `FullSimplify[theta1Direct - (theta1Fact /. ... /. slipSubs)] === 0` | claim 3 | yes |
| B5 | math | 102 | `FullSimplify[(xi1Split /. ...) - xi1Slip] === 0` | claim 4 | yes |
| B6 | math | 107-112 | `D[xi1Slip/r1Slip/theta1Fact, lamphi1/kphi] === 0` (×6) | claim 6 | no (vacuous) |

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py:40,139-144`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl:28-29,107-112`

**What's wrong:**
The support-blindness deliverable (paper notes §2.2, §8 item 1; appendix line 764 `\partial_\zeta\Xi_1=0`) is "verified" by six checks per engine that are vacuous by construction. The support symbols `lamphi1, kphi` are declared (sympy line 40; math lines 28-29) but **never wired into any expression**. The checks then differentiate the *abstract-Sigma* forms `Xi1_slip`, `R1_slip`, `Theta1_fact` — which by definition contain only `{epsW, chi0, Sigma_eps, Sigma_chi, deltaU, Sigma_Z, Sigma_del}` (confirmed by the script's own output line 97: `free symbols of Xi_1: {epsW, chi0, Sigma_eps, Sigma_chi, deltaU, Sigma_Z, Sigma_del}`) — with respect to symbols those expressions do not contain. `sp.diff(Xi1_slip, lamphi1)` is identically `0` because `lamphi1` is not a free symbol of `Xi1_slip`; the `expect_zero` passes no matter what the physics is. The same is true for all six sympy checks (lines 139-144) and all six Mathematica checks (lines 107-112). This is exactly the "variable independence" failure mode the audit self-test warns about (derivative w.r.t. an absent symbol is trivially zero).

Note also that even the microscopic-log "direct" forms (`Xi1_direct` etc.) do not contain `lamphi1`/`kphi`, because the script's construction starts from the post-`\zeta`-cancellation Stage-249 form. So differentiation w.r.t. these symbols cannot be made non-vacuous within this stage's scope; the honest in-scope check is a *structural free-symbol-absence* assertion on the microscopic-log forms.

**Why this matters:**
The transcript and paper card report support-blindness as a verified deliverable, but nothing is actually verified — if a future edit accidentally introduced a support-lane drift into the defect construction, these checks would still pass. The check provides false assurance for one of the six enumerated deliverables.

**Required change:**
Replace the six vacuous derivative checks in each engine with non-vacuous structural assertions that the support symbols are absent from the **microscopic-log "direct" forms** (`Xi1_direct`, `R1_direct`, `Theta1_direct`), which are the forms actually built from microscopic logs. See directive F1 for the exact replacement. (Do not attempt to verify the full `\zeta`-cancellation here — notes §2.2 delegates that mechanism to Stage 249; reconstructing it would be a scope extension.)

**Verification:**
After the fix, sympy lines (was 139-144) should contain assertions of the form `assert lamphi1 not in Xi1_direct.free_symbols` (and `kphi`, and the same for `R1_direct`, `Theta1_direct`); the Mathematica counterpart should use `FreeQ[xi1Direct, lamphi1]` etc. Both scripts must still exit 0. The transcript's support-blindness section should report the structural absence rather than six `= 0` derivative lines.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl` (whole file)
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:

1. Definitions — sympy lines 46-56 (`zetaZ = 2*lam1 - keta - kW; omegaW = kW - mu1; chi1 = chi0*(gam1+c1-kU); eta1 = eps_eta*(2*c1-kU-keta); varepsW = epsW*(2*gam1+2*lam1-kU-kW); deltaU1 = deltaU*(tau1-kU); eps = ...; eps1 = ...`) are reproduced verbatim in math lines 32-42 with the same right-hand sides and the same `eps`/`eps1` intermediate definitions.
2. The slippage substitution map — sympy `slip_subs` (lines 65-71) is identical to math `slipSubs` (lines 45-51), same five entries in the same order.
3. The defect construction — sympy `Xi1_direct = zetaZ - omegaW + 2*chi1/(1+chi0) + 2*eps1/(1-eps)` (lines 82-84) is the same expression as math `xi1Direct` (line 62); likewise `R1_direct`, `Theta1_direct`, `Xi1_split` and the assertion order are identical.

Both engines start from the identical algebraic ansatz and the identical intermediate forms, so they cannot independently corroborate the result — they echo each other's algebra. The second engine adds no derivation independence.

**Why this matters:**
The second-engine policy requires independent re-derivation from the physical premises so a wrong intermediate cannot pass both engines. A transliteration would reproduce any algebra error in the `.py`. The agreement in the transcripts is therefore weak evidence.

**Required change:**
See directive F2. Have the Mathematica script build at least one of the load-bearing forms by an independent route — e.g. construct `xi1Direct` from the Stage-249 `eps1`/`varepsW` expressions symbolically and then verify the four-slippage law by an alternative reduction (collect on the Sigma symbols and compare coefficients) rather than re-typing the same `Xi1_split`/`Xi1_slip` literals. At minimum, the Mathematica `xi1Slip`, `r1Slip`, `theta1Fact`, `xi1Split` target forms should be re-derived by Mathematica algebra (e.g. `Collect`/`Coefficient` on the direct form) instead of copied from the SymPy script's hand-typed coefficients.

**Verification:**
The verifier confirms the `.wl` no longer mirrors the `.py` line-for-line in the defect construction, and that the engine still exits 0 with all checks passing.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py:74-79`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl:54-59`

**What's wrong:**
The six "Physical branch drifts" checks restate the Stage-30 definitions and then subtract them from themselves. E.g. sympy line 46 defines `zetaZ = 2*lam1 - keta - kW`, then line 74 asserts `zetaZ - (2*lam1 - keta - kW) == 0` — `X - X == 0` by construction. The same pattern holds for `omega_W`, `chi_1`, `eta_1`, `varepsilon_W`, `delta_U,1` (lines 75-79) and the Mathematica mirror (lines 54-59). These can never fail regardless of the physics; they exercise no claim beyond the definition that was just typed.

**Why this matters:**
Low: they are mislabeled as verification of deliverable (5) but verify nothing. They inflate the apparent check count without adding assurance. They are harmless but should not be counted as substantive checks.

**Required change:**
See directive F3. Either remove the six tautological checks, or — if they are meant to document the Stage-30 inputs — convert at least one into a substantive cross-check (e.g. assert the slip definitions reproduce the Stage-30 ratios via `\Sigma_\chi == chi1/chi0`, `\Sigma_\eta == eta1/eps_eta`, `\Sigma_\delta == deltaU1/deltaU`, which the notes §3 state as equalities and which are NOT tautological since `chi1`, `eta1`, `deltaU1` are independently defined). The latter is preferred because it anchors deliverable (5) to the notes §3 boxed equalities.

**Verification:**
Either the six tautological lines are gone, or they are replaced by `\Sigma_\chi = chi1/chi0`-style checks that reference independently-built quantities. Script still exits 0.

## Independent-derivation check (Mathematica)

See F2. The `.wl` is a structural transliteration of the `.py`: identical variable definitions (math 32-42 vs sympy 46-56), identical `slipSubs` map (math 45-51 vs sympy 65-71), identical defect/`R_1`/`Theta_1`/split constructions and assertion ordering. There is no independent derivation path; the second engine re-types the same algebra.

## Engine cross-check

Both engines pass all checks and agree on the final simplified forms. SymPy output (lines 27, 44, 65, 76, 98-103) and Mathematica output (lines 33-34, 40-41, 47-48, 55-56, 63-74) report `0`/`PASS` for every check. The printed closed forms agree up to normalization: SymPy's `Theta_1 = -Σ_tr·χ₀·deltaU / ((χ₀+1)(deltaU+1)(χ₀+deltaU+1))` (output line 69-71) equals Mathematica's `-((chi0*deltaU*sigmaTr)/((1+chi0)*(1+deltaU)*(1+chi0+deltaU)))` (output line 50). The agreement is real but weak evidence given F2 (transliteration). The Mathematica "free symbols of Xi_1" print (output line 62) is garbled formatting of `List @@ xi1Slip`, but it is a `Print`, not an assertion, so it is not load-bearing.

## Verdict justification

The four substantive identities (four-slippage law, `\mathcal R_1`, `\Theta_1` factorization, tracking/nontracking split) hold up under attack: I checked the coefficients (`11`, `9`, factors of 2, signs of `omegaW` and the `\Sigma_\delta` term), confirmed `Xi1_direct` is genuinely built from microscopic logs (not pre-baked), confirmed `slip_subs`/`slipSubs` map the Sigma symbols to the notes §3 definitions verbatim, and confirmed `Xi1_split` vs `Xi1_slip` are independently typed forms so the equality is non-tautological. Those four match the appendix and notes exactly. The verdict is `findings` because (F1) the support-blindness deliverable is verified only vacuously — the support symbols are declared but never enter the differentiated expressions, so all twelve support-blindness checks pass trivially; (F2) the Mathematica script is a line-by-line transliteration rather than an independent re-derivation; and (F3) the six "physical branch drifts" checks are tautological restatements of the Stage-30 definitions. No `paper_misalignment`: the paper's identities are correct and the script's substantive checks match them; the failures are script insufficiencies, not paper-vs-script disagreements. No stop-cold: the fixes are local and do not propagate downstream (the carried-forward results `\Sigma_{\rm tr}`, the four-slippage law are unchanged; only the verification rigor improves).

Minor (not a numbered finding): both scripts' opening `banner` prints "STAGE 165 …" (sympy line 34, math line 26) although this is Stage 182. This is a cosmetic label error in a transcript cited by the paper card; the directive folds the one-line fix into F1's edits.

## Self-test notes

Checked traps: (1) Variable independence — confirmed `lamphi1`/`kphi` are not free symbols of `Xi1_slip`/`R1_slip`/`Theta1_fact` NOR of `Xi1_direct`/`R1_direct`/`Theta1_direct` (the script never wires them in), so the existing derivative checks are trivially zero; my prescribed replacement uses free-symbol-membership assertions, which are non-vacuous (they would fail if a support log were introduced) and do not introduce a new `sp.diff` over an absent symbol. (2) No unbounded integrals appear in this unit, so the parity trap is n/a. (3) Trivial-case: `lamphi1 not in Xi1_direct.free_symbols` evaluates to a definite boolean (True given current construction), so the assertion is well-formed. (5) Paper round-trip: the F1/F3 replacements introduce no new constants — they reuse `chi1/chi0`, `eta1/eps_eta`, `deltaU1/deltaU` which are the notes §3 boxed equalities — so no new `paper_misalignment` is created.
