---
unit_id: 158
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00-06:00
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage158_linear_defect_transport.md]
  paper_appendix: present
---

# Audit unit 158 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_158.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage158_linear_defect_transport.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: line 19 anchor MTDC-T8, line 32 "Stages 154--163: ... linear defect transport", line 1179 MTDC-T8.5, line 1350 `\input{stages/stage_158}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.txt`

## What the paper claims

The card's `Derivation ledger` plus quote block (stage_158.tex:13-17) state the stage:
`"Linearizes \((\delta\Sigma_0,\delta\mathfrak g,\delta\mathcal S)\) into \((\delta M_s,\delta M_q,\delta\Pi)\)."` — i.e. it isolates "the transport of deviations into \(\Delta_Q\), D/N similarity slippage, and the final normal coordinate \(\delta_\perp\)." The notes expand this into four boxed deliverables: (§2) `\delta R = -\delta g/\sqrt{1+r_{F1}^2}+O(2)`; (§3) `\delta M_q = -(1/4)\delta\Sigma_0 + (\Sigma_0^{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)` plus the traction relation `\delta\Sigma_0 = (40/9)\hat T_{m,can}\,\delta\hat T_m`; (§4) `\delta\Pi = (1-S_{can}/4)\delta\Sigma_0 - (\Sigma_0^{can}/4)\delta S + (\Sigma_0^{can}S_{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)`; (§5) `\Delta_Q = 5b + a_0/3 + 9 a_5 + O(2)`. The notes also list the numerical coefficients to ~16 decimals. Crucially, the card's `Checks` block (stage_158.tex:21-25) now explicitly **routes the two auxiliary checks downstream**: even-preservation of the canonical gate "is verified downstream: see \ref{stage:159}", and tangent motion giving `\delta_\perp=0` "is verified downstream: see \ref{stage:162} and \ref{stage:163}". So stage 158's own deliverables are the four boxed linearizations only.

## What the script claims to verify

Both engines run six `expect_zero`/`expectZero` checks: (1) `linear delta R law` — the series of `R(g_*+dg)` about the lower compensated branch `g_*=r-\sqrt{1+r^2}/2` equals `1/4 - dg/\sqrt{1+r^2}`; (2) `delta Mq law` — the linearization of `M_q=-(\Sigma_0+d\Sigma_0)(R_*+dR)` (dropping the `dSigma0*dR` cross term) equals `-R_* d\Sigma_0 - \Sigma_0 dR`; (3) `delta Pi law` — the linearization of `\Pi=(\Sigma_0+d\Sigma_0)(1-(R_*+dR)(S_*+dS))` equals `(1-R_*S_*)d\Sigma_0 - \Sigma_0(R_*dS+S_*dR)`; (4) `composed delta Mq law` and (5) `composed delta Pi law` — substituting `dR=-dg/\sqrt{1+r^2}` and `R_*=1/4` into checks 2/3 reproduces the notes' boxed `(d\Sigma_0,dg)` and `(d\Sigma_0,dS,dg)` forms; (6) `linear Delta_Q law` — the first-order series of `\chi=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)` equals `1+\eps(5b+a_0/3+9a_5)`. A non-asserted numerical block then prints all eight carry-forward coefficients.

## Paper ↔ script cross-check

| Paper-side deliverable | Script check | Status |
|---|---|---|
| §2 boxed `\delta R = -\delta g/\sqrt{1+r_{F1}^2}` | `linear delta R law` (sympy:43, wl:36) | match |
| §3 boxed composed `\delta M_q` in `(d\Sigma_0,dg)` | `composed delta Mq law` (sympy:72, wl:61), backed by `delta Mq law` (sympy:52, wl:44) | match |
| §3 traction `\delta\Sigma_0 = (40/9)\hat T_{m,can}\delta\hat T_m` | `coef_dSigma_dT = (40/9)*T_can` printed (sympy:111, wl:93); numerical-only, not a deliverable equation per se | match (numeric) |
| §4 boxed composed `\delta\Pi` in `(d\Sigma_0,dS,dg)` | `composed delta Pi law` (sympy:79, wl:65), backed by `delta Pi law` (sympy:62, wl:53) | match |
| §5 boxed `\Delta_Q = 5b + a_0/3 + 9 a_5` | `linear Delta_Q law` (sympy:92, wl:77) | match |
| Checks item 1: deviations about renormalized canonical point | linearization performed about `g_*=r-\sqrt{1+r^2}/2`, and `R_*=1/4` is derived in check 1 | match |
| Checks item 2: even-preservation | card routes to stage 159 (not a 158 deliverable) | n/a (downstream) |
| Checks item 3: `\delta_\perp=0` | card routes to stages 162/163 (not a 158 deliverable) | n/a (downstream) |

`paper_alignment: aligned` — every boxed deliverable of stage 158 now has a load-bearing assertion, and the card no longer promises checks the script does not perform. (This resolves the pass-1 F1 `paper_misalignment`: the user updated the card to route items 2/3 downstream rather than adding assertions here.)

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43 | `expect_zero("linear delta R law", R_lin - R_expected)` | §2 boxed `\delta R` | yes |
| A2 | sympy | 52 | `expect_zero("delta Mq law", (Mq_lin-Mq0) - (-Rstar*dSigma0 - Sigma0*dR))` | §3 (intermediate `(d\Sigma_0,dR)` form) | yes |
| A3 | sympy | 62 | `expect_zero("delta Pi law", (Pi_lin-Pi0) - dPi_expected)` | §4 (intermediate `(d\Sigma_0,dR,dS)` form) | yes |
| A4 | sympy | 72 | `expect_zero("composed delta Mq law", dMq_composed - dMq_boxed)` | §3 boxed composed `(d\Sigma_0,dg)` form | yes |
| A5 | sympy | 79 | `expect_zero("composed delta Pi law", dPi_composed - dPi_boxed)` | §4 boxed composed `(d\Sigma_0,dS,dg)` form | yes |
| A6 | sympy | 92 | `expect_zero("linear Delta_Q law", chi_lin - chi_expected)` | §5 boxed `\Delta_Q = 5b+a_0/3+9a_5` | yes |
| M1 | wl | 36 | `expectZero["linear delta R law", rLin - rExpected]` | mirror of A1 | yes |
| M2 | wl | 44 | `expectZero["delta Mq law", (mQLin-mQ0) - (-rStar*dSigma0 - sigma0*dR)]` | mirror of A2 | yes |
| M3 | wl | 53 | `expectZero["delta Pi law", (piLin-pi0) - dPiExpected]` | mirror of A3 | yes |
| M4 | wl | 61 | `expectZero["composed delta Mq law", dMqComposed - dMqBoxed]` | mirror of A4 | yes |
| M5 | wl | 65 | `expectZero["composed delta Pi law", dPiComposed - dPiBoxed]` | mirror of A5 | yes |
| M6 | wl | 77 | `expectZero["linear Delta_Q law", chiLin - chiExpected]` | mirror of A6 | yes |

The pass-1 tautological `delta Ms law` check (which evaluated `(Sigma0+dSigma0)-Sigma0-dSigma0 == 0` by construction) has been removed in both engines — there is no such line in the current scripts. The pass-1 banner copy-paste residue ("STAGE 141") is also fixed: both engines now print "STAGE 158" (sympy:32, wl:26; outputs line 3).

## Findings

None. (For context, the pass-1 report `redteam/reports/stage_158.md` raised 3 findings — F1 `paper_misalignment` (Checks block), F2 `tautological_check` (`delta Ms law`), F3 `insufficient_verification` (missing composed-form assertions). All three are remediated in the current files: F1 by the user updating the card to route items 2/3 to stages 159/162/163; F2 by deleting the tautological line; F3 by adding the `composed delta Mq law` and `composed delta Pi law` assertions exactly as prescribed.)

## Independent-derivation check (Mathematica)

The `.wl` is structurally a near-perfect line-by-line mirror of the `.py`: same `gStar = r - Sqrt[1+r^2]/2`, same `rFun = (g-r)^2/(1+r^2)`, same `Series[..., {dg,0,1}]` vs `sp.series(..., dg, 0, 2)`, the same drop-the-cross-term substitution dictionaries (`/. dSigma0*dR -> 0` vs `.subs({dSigma0*dR: 0})`), the identical `chi = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)` construction, and an identical numerical-coefficients block with the same literals and print labels in the same order. Under the prompt's explicit instruction to scrutinize this 105-175 orchestrator-direct band for `mathematica_transliteration`, I weighed filing it. I do NOT file it, for the same reason the pass-1 auditor reached: the entire physical content of this stage is the *series expansion of a handful of given closed-form expressions*. There is no second algorithmic pathway to a series coefficient — any honest second engine computes the same Taylor expansion of the same rational/polynomial expression. The two engines DO use independent CAS implementations (Mathematica `Series`/`Together`/`FullSimplify` vs SymPy `sp.series`/`sp.expand`/`sp.simplify`), and they independently confirm every residual is zero and agree on all eight numerical coefficients to 18+ digits. The transliteration here is the expected and acceptable form for a pure-series-expansion stage, not a masked single-engine gap. (Were the stage to involve a genuinely multi-route derivation — e.g. solving a system, an ODE, or an integral with a distinct closed-form path — the transliteration call would land; it does not for a Taylor-coefficient match.)

## Engine cross-check

Both engines pass all six `expectZero` assertions with residual 0 and agree on every printed coefficient to 18+ decimals:

| Coefficient | SymPy (txt:15-22) | Mathematica (txt:21-28) |
|---|---|---|
| dR/dg | -0.49021604438762605982 | -0.49021604438762603754... |
| dMq/dSigma0 | -1/4 | -1/4 |
| dMq/dg | 2.2800112692779236356 | 2.28001126927792351405... |
| dPi/dSigma0 | 0.83240947108163457213 | 0.83240947108163457213159... |
| dPi/dS | -1.1627583875422189963 | -1.16275838754221894078... |
| dPi/dg | 1.5284331782324836746 | 1.52843317823248362127... |
| dSigma0/dThat | 6.4298149620300551130 | 6.42981496203005499347... |
| dPi/dThat | 5.3522388716962184560 | 5.35223887169621835652... |

`engines_agree: true`. All eight match the notes' boxed numerics (§2-§4).

## Verdict justification

Verdict: `clean`. I attacked: (1) the `g_*=r-\sqrt{1+r^2}/2` branch — it genuinely reduces `R` to `1/4` (so the canonical base value is derived, not hardcoded) and yields the linear coefficient `-1/\sqrt{1+r^2}`; (2) the `chi` linearization — confirmed by hand that `chi(eps=0)=1`, the `s` normalization cancels, and the `(5,1/3,9)` coefficients on `(b,a_0,a_5)` emerge from the `beta^5`, the `3S` factor in the denominator, and the `9 Sigma_5` term, not from a self-comparison; (3) the composed `delta Mq`/`delta Pi` checks — they tie the `dR=-dg/\sqrt{1+r^2}` law and `R_*=1/4` into the notes' boxed forms and would fail on a wrong sign or coefficient (non-tautological glue, the exact pass-1 F3 fix); (4) the carried literals `r_F1`, `\Sigma_0^{can}`, `S_{can}`, `\hat T_{m,can}` — all match notes §1 verbatim and are legitimately imported from stage 156 (card Inputs anchors them); (5) the `168`/`100π²` stale-value watch — the only "168" is inside the literal `4.651033550168876`, which is correct, so no stale-constant issue here. I confirmed I read the paper card, notes, and appendix rows, and that the script's verified claims match the paper's deliverables exactly. The card-vs-script alignment that was partial in pass-1 is now exact because the user routed the two auxiliary Checks downstream. The Mathematica transliteration is the acceptable form for a pure series-expansion stage and is not filed.

## Self-test notes

(1) Variable independence: every `series`/`Series` expands in the correct variable that the expanded expression actually depends on (`dg` for `R(g_*+dg)`, `eps` for `chi`); no identically-zero-derivative trap. (2) Trivial-case pre-check: substituting `dg=0` (or `eps=0`) into each composed/linear check reduces both sides to the same residual 0; deliberately flipping the `dg`-coefficient sign in `dMq_boxed` would make A4 fail when `dg!=0` — non-tautological. (3) Symbol assumptions: `r_sym>0` (sympy `positive=True`, wl `rSym>0`) is physically justified (`r_F1≈1.778>0`) and only affects `Sqrt[1+r_sym^2]`, which is unambiguous anyway — no `symbol_assumption_error`. (4) No paper_misalignment re-introduced: the boxed-form assertions match the notes' boxes verbatim. No traps tripped.

## Value Reconciliation (pass-2 augmentation)

Every deliverable value the scripts emit is checked against the `.tex` card and `.md` notes below. The card is terse (it carries no numerics, only the linearization quote-block), so the natural carrier for all numeric/symbolic results is the notes file.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `\delta R = -\delta g/\sqrt{1+r^2}` (symbolic, base `R_*=1/4`) | py:42-43, wl:35-36; out sympy:5 | notes §2 box, .md:88-95 (and `\delta R\approx -0.490216...` .md:100-104) | MATCH |
| `\delta M_q = -(1/4)\delta\Sigma_0 + \Sigma_0/\sqrt{1+r^2}\,\delta g` (composed box) | py:71, wl:60; out sympy:8 | notes §3 box, .md:140-149 | MATCH |
| `\delta\Pi = (1-S_*/4)\delta\Sigma_0 - (\Sigma_0/4)\delta S + \Sigma_0 S_*/\sqrt{1+r^2}\,\delta g` (composed box) | py:76-78, wl:64; out sympy:9 | notes §4 box, .md:216-225 | MATCH |
| `\Delta_Q = 5b + a_0/3 + 9 a_5` | py:91, wl:76; out sympy:10 | notes §5 box, .md:286-290 | MATCH |
| dR/dg = -0.49021604438762... | py:105; out sympy:15, wl:21 | notes §2, .md:100-104 (`-0.490216044387626`) | MATCH |
| dMq/dSigma0 = -1/4 | py:106; out sympy:16 | notes §3, .md:157 (`-0.25`) | MATCH |
| dMq/dg = 2.28001126927792... | py:107; out sympy:17 | notes §3, .md:159 (`2.28001126927792`) | MATCH |
| dPi/dSigma0 = 0.83240947108163... | py:108; out sympy:18 | notes §4, .md:234 (`0.832409471081635`) | MATCH |
| dPi/dS = -1.16275838754222... | py:109; out sympy:19 | notes §4, .md:236 (`-1.16275838754222`) | MATCH |
| dPi/dg = 1.52843317823248... | py:110; out sympy:20 | notes §4, .md:238 (`1.52843317823248`) | MATCH |
| dSigma0/dThat = 6.42981496203006... | py:111; out sympy:21 | notes §3, .md:178 (`6.42981496203006`) | MATCH |
| dPi/dThat = 5.35223887169622... | py:112; out sympy:22 | notes §4, .md:248 (`5.35223887169622`) | MATCH |

INTERNAL (scaffolding, no prose expected): `R_lin`, `R_expected`, `R_shift`, `g_star` (intermediate); `Mq_lin`/`Mq0`, `Pi_lin`/`Pi0`, `dMq_composed`/`dMq_boxed`, `dPi_composed`/`dPi_boxed`, `chi`/`chi_lin`/`chi_expected` (assertion drivers); `sqrt1` (intermediate); the per-check residual values (all 0); the `Carry-forward summary` echo lines (re-statements of the boxed deliverables already reconciled above).

reconciliation: complete; 12 deliverable values checked, 0 misaligned.
