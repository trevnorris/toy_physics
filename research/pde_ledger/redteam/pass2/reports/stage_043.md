---
unit_id: 043
batch: III.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage043_support_direction_extraction.md]
  paper_appendix: present
---

# Audit unit 043 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_043.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage043_support_direction_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 64; `\input` row 204)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage043_support_direction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage043_support_direction_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` is verbatim: "A physical support-direction insertion rule \eqref{eq:app-stage043-replacements}", i.e. the replacement rule `q=tR_U, r=tR_phi, m=M_mix, n=M_supp` with `t=kappa_1/kappa_0`, to be inserted into the Stage-041 rank-2 determinant. The card body frames the stage as resolving whether `R_phi = R_U` (tracking / source-tied) or `R_phi != R_U` (a new rank-2 branch). The `.tex` card is terse; the notes file is the authoritative derivation and enumerates seven exact deliverables: (1) the rank-1 effective support vector `y = gB v + gU gS D_U v`; (2) the support-direction factor `R_phi = [1 + sigma0/(1+delta_U)]/(1+sigma0)`; (3) the splitting invariant `D_phi = -kappa0 kappa1 gB sigma0 delta_U/(1+delta_U)` vanishing iff `sigma0=0` or `delta_U=0`; (4) the split support-blocking ratio `eps_phi^split = eps_phi[1 - (2/11) delta_U/(1+delta_U)]`; (5) the physical support baseline `M_supp = 8 Z_phi (1+sigma0)^2 / [pi^2 (1-eps_eta)(1-eps_phi^split)]`; (6) the tracking condition `g_B g_R = g_W g_S` (equivalently `sigma0=rho0`) via `D_(phi z)=0`; (7) the mismatch `R_phi - R_U = delta_U(rho0-sigma0)/[(1+delta_U)(1+rho0)(1+sigma0)]`. The notes also pin the overlap constants `sigma = kappa0^2+kappa1^2 = 88/(9 pi^2)` and `kappa1^2/sigma = 2/11` (so `kappa0^2/sigma = 9/11`, `kappa0^2 = 8/pi^2`). The appendix row (64) is a one-line status: "Physical support direction extraction ... Continuum support direction data needed for the tracking/source-tied decision."

## What the script claims to verify

Both engines verify the seven notes deliverables as exact symbolic identities. They construct `y = gB v + gU gS D_U v` from the diagonal eliminator `D_U = diag(1/K_U, 1/(K_U(1+delta_U)))`, then assert (i) `R_phi` (computed from `y`) equals the closed form in `sigma0`; (ii) `D_phi` equals the expected splitting invariant; (iii) the overlap contraction `v.D_U.v` equals `(sigma/K_U)[1-(2/11)delta_U/(1+delta_U)]` after substituting `kappa1^2=(2/11)sigma`; (iv) the renormalized pole `A_phi^eff` and the split/minimal overlap ratio; (v) `M_supp` is independent of the bare masses `mu_eta, mu_phi`, matches a hand-written structural form with a free baseline `B`, and evaluates at `B=8/pi^2`; (vi) the tracking determinant `D_(phi z)` and `R_phi=R_U` under `g_B g_R=g_W g_S`; (vii) the mismatch closed form. The Mathematica script adds independent cross-checks the SymPy lacks: `delta_U->0` and `delta_U->infinity` endpoint limits of the overlap, a `Limit`-based minimal pole, a leading-order series check of the mismatch, and three numeric test-point sign checks of `R_phi-R_U`.

## Paper ↔ script cross-check

| paper-side deliverable (notes §) | script-side check | status |
|---|---|---|
| (1) effective support vector `y` | `y = gB*v + gU*gS*DU*v` (sympy 64; wl 45) | match |
| (2) `R_phi` factor | `R_phi - Rphi_expected == 0` (sympy 74; wl 58) | match |
| (3) splitting invariant `D_phi` | `D_phi - expected == 0` (sympy 81; wl 59) | match |
| (4) `eps_phi^split` / overlap / pole | overlap contraction, `A_phi^eff`, split-vs-min ratio (sympy 93/101/108/115; wl 87/88/93/98) | match |
| (5) baseline `M_supp = 8 Z_phi(1+sigma0)^2/...` | mu-independence + structural form + `B=8/pi^2` eval (sympy 127-148; wl 110-129) | partial (structural form exact; the literal `8` is substituted via `B=8/pi^2`, never derived from `kappa0^2`) |
| (6) tracking condition `g_B g_R=g_W g_S` | `D_(phi z)` + tracking via `gS=gB gR/gW` (sympy 158/162; wl 144/145) | match |
| (7) mismatch `R_phi-R_U` | `mismatch - expected == 0` (sympy 170; wl 160) + wl sign/series checks | match |
| Output insertion rule `(R_U,R_phi,M_mix,M_supp)` | covered piecewise: `R_phi` (2), `R_U`/`rho0` (6), `M_supp` (5); `M_mix`/`q=tR_U`,`r=tR_phi` are carry-forwards stated only in the card | match (the insertion rule is assembled from the verified pieces; `t`, `q`, `r`, `M_mix` are framing/carry-forward, not new claims) |

`paper_alignment: aligned` — every notes deliverable has a faithful, non-tautological script check; the only soft spot is the redundant baseline-value evaluation (F2). The terse `.tex` card omits the numeric constants, but the notes carry them correctly, so per the augmentation guards those are MATCH.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74 | `simplify(R_phi - Rphi_expected) == 0` | deliverable 2 | yes |
| A2 | sympy | 81 | `D_phi - Dphi_expected == 0` | deliverable 3 | yes |
| A3 | sympy | 93 | overlap contraction `== 0` | deliverable 4 | yes |
| A4 | sympy | 101 | `A_phi^eff - expected == 0` | deliverable 4 | yes |
| A5 | sympy | 108 | `A_phi^eff(delta_U=0) - expected == 0` | deliverable 4 (limit) | yes |
| A6 | sympy | 115 | split-vs-minimal overlap ratio `== 0` | deliverable 4 | yes |
| A7 | sympy | 127-128 | `diff(M_supp, mu_eta/mu_phi) == 0` | deliverable 5 (mass-cancellation) | yes |
| A8 | sympy | 141 | `M_supp_in_B - struct_expected == 0` | deliverable 5 (structural form) | yes |
| A9 | sympy | 148 | `M_supp(B=8/pi^2) - expected == 0` | deliverable 5 (numeric baseline) | no (redundant after A8; see F2) |
| A10 | sympy | 158 | `D_(phi z) - expected == 0` | deliverable 6 | yes |
| A11 | sympy | 162 | `(Rphi_exp - RU).subs(gS=gB gR/gW) == 0` | deliverable 6 (tracking) | yes |
| A12 | sympy | 170 | `mismatch - expected == 0` | deliverable 7 | yes |
| A13 | wl | 58-59 | `R_phi`, `D_phi` `=== 0` | deliverables 2,3 | yes |
| A14 | wl | 85-86 | overlap endpoint limits `=== 0` | deliverable 4 (independent endpoints) | yes |
| A15 | wl | 87-88,93,98 | overlap/pole/ratio `=== 0` | deliverable 4 | yes |
| A16 | wl | 110-129 | M_supp mass-cancel + structural + `B=8/pi^2` | deliverable 5 | yes (123); wl 129 redundant like A9 |
| A17 | wl | 144-145 | `D_(phi z)`, tracking `=== 0` | deliverable 6 | yes |
| A18 | wl | 157,160 | mismatch leading-order + closed form `=== 0` | deliverable 7 | yes |
| A19 | wl | 171,178,181 | numeric sign test-points `=== 0` | deliverable 7 (independent sign) | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage043_support_direction_sympy_audit.txt` (mtime May 26 02:03) vs script (mtime Jun 3 15:59)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage043_support_direction_mathematica_audit.txt` (mtime May 26 02:04) vs script (mtime Jun 3 15:59)

**What's wrong:**
Both saved transcripts predate their scripts. The staleness is visible in content, not just mtime: the SymPy transcript banner reads `STAGE 26 — CONTINUUM EXTRACTION...` (line 3) and the ledger header reads `STAGE 26 THEOREM LEDGER` (line 106), whereas the current `.py` prints `STAGE 43` (script line 55) and `STAGE 43 THEOREM LEDGER` (line 172). Likewise the Mathematica transcript banner reads `STAGE 026 ...` (line 3) while the current `.wl` prints `STAGE 043` (script line 33). The numerical/symbolic results in the transcripts are otherwise consistent with the current scripts (every `R_phi`, `D_phi`, `A_phi^eff`, `M_supp`, `D_(phi z)`, mismatch form, and PASS line matches what the current source would emit), so this is a label/banner-only drift, not a result disagreement.

**Why this matters:**
The committed transcript misrepresents which stage produced it; a reader cross-referencing the transcript banner against the card would see "Stage 26" for a Stage-043 audit. The underlying math is unaffected.

**Required change:**
Re-run both scripts and recommit the refreshed `.txt` outputs so the banners read STAGE 43 / STAGE 043. (Note: the in-source subbanner labels `26.1`–`26.5` and the docstring "Stage 26 SymPy audit" on sympy lines 4/57/83/117/150/164 are the same SCRIPT/OUTPUT-band stale-label class that the project has explicitly deferred to a dedicated content-keyed numbering pass — they are NOT to be edited by this red-team directive.)

**Verification:**
After a fresh exec, the regenerated sympy `.txt` line 3 reads `STAGE 43 — ...` and the mathematica `.txt` line 3 reads `STAGE 043 — ...`; all PASS lines remain.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:143-148`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:125-129`

**What's wrong:**
The "baseline value identification" block (sympy 144-148, wl 126-129) computes `Msupp_cont_eval = Msupp_cont_in_B.subs(B, 8/pi^2)` and `Msupp_expected = Msupp_struct_expected.subs(B, 8/pi^2)`, then asserts their difference is zero. But line 141 (`expect_zero("M_supp structural form (free baseline)", Msupp_cont_in_B - Msupp_struct_expected)`) has ALREADY proven `Msupp_cont_in_B == Msupp_struct_expected` as an identity in the free symbol `B`. Substituting the same value `B=8/pi^2` into both sides of a proven identity cannot make them differ, so the line-148 assertion is guaranteed-true by construction and exercises nothing new. Critically, it does NOT verify the stage's actual numeric claim: that the support baseline `B = kappa0^2` really equals `8/pi^2`. The notes derive that (`kappa0^2/sigma = 9/11`, `sigma = 88/(9 pi^2)` => `kappa0^2 = (9/11)(88/(9 pi^2)) = 8/pi^2`), but the script never connects `8/pi^2` to `sigma` or `kappa0^2` — it just substitutes the literal. The `8` in deliverable 5's headline `M_supp = 8 Z_phi (1+sigma0)^2/...` is therefore asserted by labeling, not exercised.

**Why this matters:**
The redundant assertion gives false confidence that the numeric baseline `8/pi^2` has been checked, when in fact only the structural form (which is independent of that number) is verified. If `kappa0^2` were mis-stated, this stage would not catch it.

**Required change:**
Replace the redundant line-148 (sympy) / line-129 (wl) assertion with a check that actually exercises the baseline value: assert `kappa0^2 == 8/pi^2` derived from the stage's own frozen overlap constants, e.g. assert `(sp.Rational(9,11)*sp.Rational(88,1)/(sp.Integer(9)*sp.pi**2)) - sp.Rational(8,1)/sp.pi**2 == 0`, or assert `B_value == (9/11)*sigma` with `sigma -> 88/(9*pi^2)` so that the literal `8/pi^2` is produced, not assumed. Keep the structural-form check (line 141) unchanged. (Self-test caveat: this fix only rearranges which identity is asserted; it introduces no new paper_misalignment because `8/pi^2 = (9/11)*88/(9 pi^2)` matches the notes exactly.)

**Verification:**
The regenerated transcript shows a new check whose label references `kappa0^2` or `8/pi^2 = (9/11) sigma` and prints `= 0`; the structural-form PASS remains.

## Independent-derivation check (Mathematica)

The `.wl` is structurally parallel to the `.py` (same five subbanner sections, same deliverable order), which is expected when both engines derive the same physical result. It is NOT a transliteration: (a) it constructs `y = gB v + gU gS (dU.v)` directly from the physical eliminator rather than echoing the SymPy intermediate-variable choreography; (b) it uses Mathematica-native constructs the SymPy never uses — `Det[{{kappa0,kappa1},{y0,y1}}]` for `D_phi` (wl 52), `Limit[..., deltaU->Infinity]` endpoint checks (wl 81-86), `Series[mismatch,{deltaU,0,1}]` leading-order verification (wl 154-157); and (c) it adds three numeric test-point sign checks of `R_phi-R_U` (wl 166-181) that have no SymPy counterpart and independently pin the sign of the mismatch. These are genuine independent cross-checks of the same claims, so no `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit identical final forms:
- `y0 = kappa0(gB + gS gU/kU)`, `y1 = kappa1(gB + gS gU/(kU(1+deltaU)))` — sympy txt 10-16, wl txt 9.
- `R_phi - expected = 0`, `D_phi - expected = 0` — sympy txt 21/26, wl txt 12-15.
- overlap `sigma(9 deltaU+11)/(11 kU(1+deltaU))` — sympy txt 32, wl txt 20 (same after normalization).
- `A_phi^eff`, split ratio, `M_supp` structural and `B=8/pi^2` forms — sympy txt 36-84, wl txt 21-55: identical.
- `D_(phi z) = -deltaU gU kappa0 kappa1(gB gR - gS gW)/(kU(1+deltaU))` — sympy txt 90, wl txt 60.
- mismatch `deltaU gU(gB gR - gS gW) kU/[(1+deltaU)(gS gU+gB kU)(gR gU+gW kU)]` — sympy txt 100-102, wl txt 72.
All PASS. The Mathematica's extra endpoint/series/sign checks all pass (wl txt 27-29, 71, 75-80). `engines_agree: true`.

## Verdict justification

`findings`: two low-severity findings, both informational/non-blocking. I attacked every assertion: the `R_phi`, `D_phi`, overlap, `A_phi^eff`, `D_(phi z)`, and mismatch checks all compare a construction-derived quantity against an independently-written closed form, so none is tautological — I hand-verified each algebraically (e.g. `D_phi = -kappa0 kappa1 gU gS deltaU/(kU(1+deltaU))` matches the `gB sigma0` form, and the overlap `(9 deltaU+11)/(11(1+deltaU))` matches `1-(2/11)deltaU/(1+deltaU)`). Symbol domains are sane (`kU,Kphi_eff,deltaU,mu>0`; `gB,gW,kappa!=0`), and the `simplify`/`FullSimplify` assumptions don't mask a branch. The mu-independence check is genuinely structural (the masses cancel between numerator and the two pole denominators). The only weak check is the redundant `B=8/pi^2` baseline evaluation (F2), which is guaranteed-true after the structural-form check and does not exercise the numeric baseline; it is low severity because the structural claim is fully verified and the `8/pi^2` value reconciles with the notes' frozen `sigma`. The transcripts are stale at the banner level (F1, STAGE 26/026 vs current STAGE 43/043) but otherwise content-consistent. Paper alignment is exact against the authoritative notes; the `.tex` card is terse but consistent. No paper_misalignment, no stop_cold.

## Self-test notes

Checked: (1) Variable independence — the only derivatives are `diff(M_supp, mu_eta/mu_phi)` (sympy 127-128, wl 110-111); `M_supp` genuinely depends on `mu_eta, mu_phi` (they appear in numerator `1/(mu_eta mu_phi)` and in both pole factors `1/mu`), and they cancel, so the `==0` result is a real cancellation, not a vacuous zero-derivative. (2) Symmetry/parity — no unbounded integrals in this stage; the only "limit" operations are `delta_U->0/infinity` endpoints, which I verified give `sigma/kU` and `(9/11)sigma/kU`. (3) Trivial-case — substituting `sigma0=0` (gS=0) into `R_phi` gives 1 (source-tied) and into `D_phi` gives 0, matching deliverable 3; `sigma0=rho0` (gS=gB gR/gW) makes the mismatch vanish, matching deliverable 6. My F2 fix only re-points an existing assertion and re-derives `8/pi^2 = (9/11)*88/(9 pi^2)`, introducing no new misalignment.

## Value Reconciliation (pass-2 augmentation)

Authoritative carrier is the notes file `moving_throat_pde_stage043_support_direction_extraction.md` (the `.tex` card is intentionally terse and states only the symbolic insertion rule, not the constants — per the augmentation guards, notes-coverage counts as MATCH).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `y0 = kappa0 gB(1+sigma0)`, `y1 = kappa1 gB[1+sigma0/(1+delta_U)]` | sympy txt 10-16; wl txt 9 | notes:93-95 | MATCH |
| `R_phi = [1+sigma0/(1+delta_U)]/(1+sigma0)` | sympy txt 17-20; wl txt 10 | notes:103 (`.tex` eq. uses `R_phi` symbolically, 13) | MATCH |
| `sigma0 := gU gS/(KU gB)` | sympy 69; wl 49 | notes:85 | MATCH |
| `D_phi = -kappa0 kappa1 gB sigma0 delta_U/(1+delta_U)` | sympy txt 22-25; wl txt 11 | notes:121 | MATCH |
| overlap `v.D_U.v = (sigma/KU)[1-(2/11)delta_U/(1+delta_U)]` | sympy txt 31-34; wl txt 20 | notes:153 | MATCH |
| `kappa1^2/sigma = 2/11` (so `kappa0^2/sigma = 9/11`) | sympy 89; wl 65 | notes:151 (`2/11`); `9/11` implied | MATCH |
| `sigma = kappa0^2+kappa1^2 = 88/(9 pi^2)` | used implicitly (script keeps `sigma` symbolic; not emitted as a number) | notes:149 | MATCH (notes-only; script keeps symbolic — see F2) |
| `eps_phi^split = eps_phi[1-(2/11)delta_U/(1+delta_U)]` | sympy 96; wl 70; A_phi^eff txt | notes:161 | MATCH |
| `A_phi^eff = K_phi^eff - g_S^2 v.D_U.v` | sympy txt 36-40; wl txt 21 | notes:141-143 | MATCH |
| `M_supp = 8 Z_phi(1+sigma0)^2/[pi^2(1-eps_eta)(1-eps_phi^split)]` | sympy txt 72-84; wl txt 53 | notes:165-171 | MATCH (structural form exact; numeric `8` via `B=8/pi^2` substitution, see F2) |
| `Z_phi := c_etaphi^2/(K_eta^eff K_phi^eff)` | embedded in `Msupp_struct` (`cB^2/(Keta_eff Kphi_eff)`), sympy 135; wl 118 | notes:171 | MATCH |
| `D_(phi z) = -delta_U kappa0 kappa1 gU(gB gR - gW gS)/(KU(1+delta_U))` | sympy txt 90; wl txt 60 | notes:191 | MATCH |
| tracking condition `g_B g_R = g_W g_S` (`sigma0=rho0`) | sympy 162; wl 145 | notes:195-199 | MATCH |
| `rho0 := gU gR/(KU gW)`, `R_U = [1+rho0/(1+delta_U)]/(1+rho0)` | sympy 160-161; wl 140-141 | notes:107 (`R_U`), 199 (`rho0`) | MATCH |
| mismatch `R_phi-R_U = delta_U(rho0-sigma0)/[(1+delta_U)(1+rho0)(1+sigma0)]` | sympy txt 99-102; wl txt 72 | notes:211-214 | MATCH |
| baseline literal `8/pi^2` (= `kappa0^2`) | sympy 144; wl 126 | notes:149+151 (derivable: `(9/11)*88/(9 pi^2)=8/pi^2`) | MATCH (consistent; not independently re-derived in-script — F2) |

INTERNAL (scaffolding, no prose expected): `expect_zero`/`expectZero` PASS flags, all `... - expected = 0` residuals, the free baseline symbol `B`/`bBaseline`, the `mu_eta/mu_phi` mass-independence residuals, the `delta_U->0`/`delta_U->infinity` endpoint residuals, the `Series` leading-coefficient residual, and the three numeric sign test-point residuals (`1/4`, `-1/4`, `0`) — these are sign/consistency cross-checks, not deliverable values.

`reconciliation: complete; 16 deliverable values checked, 0 misaligned`. Every emitted deliverable value reconciles against the notes (the natural carrier). The only soft point is the in-script numeric baseline `8/pi^2`, which is consistent with the notes but asserted by substitution rather than derived — captured as F2 (insufficient_verification), NOT a value_mismatch, since both sides agree.
