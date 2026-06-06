---
unit_id: 085
batch: III.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage085_quadrupole_demand_cancellation.md]
  paper_appendix: present
---

# Audit unit 085 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_085.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage085_quadrupole_demand_cancellation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 148; also `\input` at line 288)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage085_quadrupole_demand_cancellation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage085_quadrupole_demand_cancellation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage085_quadrupole_demand_cancellation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage085_quadrupole_demand_cancellation_mathematica_audit.txt`

## What the paper claims

Stage 085 proves that once the selected outgoing quadrupole branch is normalized, the absolute outgoing-normalization factors cancel out of the support theorem, leaving a dependence only on the loading ratio `rho_alpha = alpha_req/alpha_mix`. The card's boxed deliverables are: (i) the product identities `\Pi_{tr} = (N_Q^{target}/\beta_0)\alpha_{req}` and `C_{mix} = (N_Q^{target}/\beta_0)\alpha_{mix}` (eq:app-stage085-products); (ii) the cancellation theorem `\Pi_{tr}/C_{mix} = \alpha_{req}/\alpha_{mix} =: \rho_\alpha` (eq:app-stage085-ratio); (iii) the reduced demand law `\zeta_{req}(\rho_\alpha,\epsilon_{blk}) = (\rho_\alpha-1)/[1-\epsilon_{blk}(2-\rho_\alpha)]` with unblocked limit `\zeta_{req}=\rho_\alpha-1` (eq:app-stage085-zeta-rho). `\stagefield{Output}` reads verbatim: "The cancellation theorem \eqref{eq:app-stage085-ratio} and reduced demand formula \eqref{eq:app-stage085-zeta-rho}." The notes add detail: the products come from `R_target = N_Q A/(beta_0 kappa_0^2)` with `kappa_0^2 = 8/pi^2`, `G_tr = 8 alpha_req/(pi^2 A)`, `M_mix = 8 alpha_mix/(pi^2 A)`, `Pi_tr=R_target G_tr`, `C_mix=R_target M_mix` (§1); the spectral forms via the Stage-030 relation `N_Q^{target} = mhat_-^2 beta_0 s_-/lam_-` ⇒ `Pi_tr = mhat_-^2 (s_-/lam_-) alpha_req`, likewise for `C_mix` (§2); and the demand law derived from the Stage-052 form `zeta_req = (Pi_tr-C_mix)/[C_mix - eps_blk(2 C_mix - Pi_tr)]` (§3). The appendix row (line 148) summarizes: "Selected support theorem depends only on `rho_alpha = alpha_req/alpha_mix`."

## What the script claims to verify

Both scripts construct `R_target = NQ*A/(beta0*kappa0^2)` with `kappa0^2 = 8/pi^2`, and `G_tr = 8*alpha_req/(pi^2*A)`, `M_mix = 8*alpha_mix/(pi^2*A)`; form `Pi_tr = R_target*G_tr` and `C_mix = R_target*M_mix`; then assert that the `pi^2`/`A`/`8` factors cancel so `Pi_tr - NQ*alpha_req/beta0 == 0`, `C_mix - NQ*alpha_mix/beta0 == 0`, and `Pi_tr/C_mix - alpha_req/alpha_mix == 0`. They then substitute the Stage-030 spectral relation `NQ -> mhat^2*beta0*s_-/lam_-` and assert the spectral forms `Pi_sel - mhat^2 s_- alpha_req/lam_- == 0` and `C_sel - mhat^2 s_- alpha_mix/lam_- == 0`. Finally they build the demand law `zeta_req = (Pi_tr-C_mix)/[C_mix - eps_blk(2 C_mix - Pi_tr)]` and assert: it equals the pure-loading form `(alpha_req-alpha_mix)/[alpha_mix - eps_blk(2 alpha_mix - alpha_req)]`; that loading form under `alpha_req -> rho_alpha*alpha_mix` equals `(rho_alpha-1)/[1-eps_blk(2-rho_alpha)]`; and that the unblocked limit (`eps_blk=0`) equals `alpha_req/alpha_mix - 1`. This matches all three paper deliverables plus the spectral substitution.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (i) Product identities `Pi_tr=(N_Q/beta0)alpha_req`, `C_mix=(N_Q/beta0)alpha_mix` (eq:products) | sympy L39–40 / wl L48–49 `expect_zero(Pi_tr - NQ*alpha_req/beta0)`, `(C_mix - NQ*alpha_mix/beta0)` | match |
| (ii) Cancellation theorem `Pi_tr/C_mix = alpha_req/alpha_mix` (eq:ratio) | sympy L41 / wl L50 `expect_zero(Pi_tr/C_mix - alpha_req/alpha_mix)` | match |
| (iii) Reduced demand `zeta_req=(rho_alpha-1)/[1-eps_blk(2-rho_alpha)]` (eq:zeta-rho) | sympy L60,L61–64 / wl L66–67 product↔loading + loading↔rho_alpha checks | match |
| unblocked limit `zeta_req=rho_alpha-1` | sympy L65–68 / wl L68 `expect_zero(zeta_expected.subs(eps_blk,0) - (alpha_req/alpha_mix-1))` | match |
| spectral forms `Pi_tr=mhat^2(s_-/lam_-)alpha_req` (notes §2) | sympy L52–53 / wl L59–60 `Pi_sel - mhat^2 s_- alpha_req/lam_-` etc. | match |

`paper_alignment: aligned` — every paper-side deliverable and notes-level intermediate has a faithful, non-tautological script-side check; no extras and no mismatches.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero(Pi_tr - NQ*alpha_req/beta0)` | (i) product identity | yes |
| A2 | sympy | 40 | `expect_zero(C_mix - NQ*alpha_mix/beta0)` | (i) product identity | yes |
| A3 | sympy | 41 | `expect_zero(Pi_tr/C_mix - alpha_req/alpha_mix)` | (ii) cancellation theorem | yes |
| A4 | sympy | 52 | `expect_zero(Pi_sel - mhat^2 s_- alpha_req/lam_-)` | spectral form (notes §2) | yes |
| A5 | sympy | 53 | `expect_zero(C_sel - mhat^2 s_- alpha_mix/lam_-)` | spectral form (notes §2) | yes |
| A6 | sympy | 60 | `expect_zero(zeta_req - zeta_expected)` | (iii) demand law (product↔loading) | yes |
| A7 | sympy | 61–64 | `expect_zero(zeta_expected[alpha_req→rho*alpha_mix] - rho-form)` | (iii) rho_alpha form | yes |
| A8 | sympy | 65–68 | `expect_zero(zeta_expected[eps_blk→0] - (alpha_req/alpha_mix-1))` | (iii) unblocked limit | yes |
| B1 | mathematica | 48 | `expectZero[piTr - nQ*alphaReq/beta0]` | (i) | yes |
| B2 | mathematica | 49 | `expectZero[cMix - nQ*alphaMix/beta0]` | (i) | yes |
| B3 | mathematica | 50 | `expectZero[piTr/cMix - alphaReq/alphaMix]` | (ii) | yes |
| B4 | mathematica | 59 | `expectZero[piSel - mhat^2*sMinus*alphaReq/lamMinus]` | spectral | yes |
| B5 | mathematica | 60 | `expectZero[cSel - mhat^2*sMinus*alphaMix/lamMinus]` | spectral | yes |
| B6 | mathematica | 66 | `expectZero[zetaReq - zetaExpected]` | (iii) product↔loading | yes |
| B7 | mathematica | 67 | `expectZero[(zetaExpected/.alphaReq→rhoAlpha*alphaMix) - zetaRho]` | (iii) rho form | yes |
| B8 | mathematica | 68 | `expectZero[(zetaExpected/.epsBlk→0) - (alphaReq/alphaMix-1)]` | (iii) unblocked limit | yes |

All rows "yes": each assertion non-tautologically exercises a paper deliverable. See Verdict justification for the why-not-tautological argument on the load-bearing A1/A3 checks.

## Findings

None. (zero findings)

The adversarial probes below all failed to break the scripts:

**Tautology probe (A1/A2 — the load-bearing cancellation).** `Pi_tr` is *not* defined as `NQ*alpha_req/beta0`; it is defined as `R_target*G_tr` where `R_target = NQ*A/(beta0*(8/pi^2))` and `G_tr = 8*alpha_req/(pi^2*A)`. The assertion `Pi_tr - NQ*alpha_req/beta0 == 0` therefore tests a *genuine* cancellation of `pi^2`, `A`, and the factor `8`: the residual is zero only because `kappa0^2 = 8/pi^2` exactly matches the `8/(pi^2 A)` numeric factor in `G_tr`. If `kappa0^2` were, say, `4/pi^2`, the residual would be `NQ*alpha_req/(2*beta0)` ≠ 0 and the check would fail. This is the physics content of the stage and the check is well-anchored, not tautological.

**Tautology probe (A6 — zeta product↔loading).** `zeta_req` is built from the *derived* `Pi_tr`/`C_mix` (in `NQ, beta0, alpha_*` variables), while `zeta_expected` is built directly from `alpha_*`. They agree only because the common factor `NQ/beta0` cancels between numerator and denominator (each homogeneous degree-1 in the product `NQ/beta0`). Not a defined-then-asserted round trip — it confirms the cancellation propagates through the full rational demand expression.

**Tautology probe (A7/A8).** A7 requires `alpha_mix` to cancel after the `alpha_req → rho_alpha*alpha_mix` substitution (non-trivial). A8's `eps_blk=0` reduction `(alpha_req-alpha_mix)/alpha_mix = alpha_req/alpha_mix - 1` is a real simplification, not an identity-by-construction.

**Symbol-assumption probe.** `eps_blk` is declared `real` (not positive) — correct, since the blocking fraction needs no sign restriction for these rational identities. All other symbols (`A, beta0, NQ, alpha_req, alpha_mix, mhat, s_minus, lam_minus, rho_alpha`) are `positive, real`, consistent with their physical roles (normalizations, overlaps, eigenvalues, loadings) and with the paper's setup. No positivity is doing illegitimate work in any `simplify`. The Mathematica `$Assumptions` (L29–32) mirror these exactly. No `symbol_assumption_error`.

**Hardcoded-result probe.** The only literal is `kappa0_sq = 8/pi^2` (sympy L24, wl L34), which the notes §1 explicitly state as "the exact D/N constant `kappa_0^2 = 8/pi^2`" carried from Stages 035–036. It is anchored in the notes, not an unexplained magic number. The `8/(pi^2 A)` factors in `G_tr`/`M_mix` likewise match notes §1. No `hardcoded_result`.

**Missing-branch probe.** The card's only branch distinction is the unblocked limit `eps_blk=0`, which is checked explicitly (A8). The blocked general form is checked symbolically (A6/A7), covering all `eps_blk`. No `missing_branch`.

**Stale-output probe.** Both `.txt` outputs are newer than their scripts (sympy: script 10:15:47, output 10:24:50; mathematica: script 10:15:49, output 10:26:08, all 2026-05-27) and their contents match what the current scripts would emit. Outputs fresh; no `stale_output`.

**Self-label/numbering probe.** Both scripts and both outputs label themselves "STAGE 085" correctly (sympy L19, wl L26, outputs L3). No stale "Stage NN" self-labels. The notes reference upstream Stages 030, 035–036, 052, 084 — these are cross-references (deferred, out of scope) and they are correct anyway. No numbering finding.

## Independent-derivation check (Mathematica)

The `.wl` is structurally parallel to the `.py`: same three building blocks (`rTarget`, `gTr`, `mMix`), same products, same eight assertions in the same order. For this stage that is acceptable, not a `mathematica_transliteration` finding. The stage's entire content is a short algebraic cancellation identity (`R_target*G_tr` collapsing to `NQ*alpha_req/beta0` once `kappa0^2=8/pi^2`); there is no meaningfully distinct "second derivation route" — both engines must build the same expression tree from the same physical premises (notes §1: the `R_target`/`G_tr`/`M_mix` definitions). The two engines run independent simplifiers (`sympy.simplify∘expand` vs `FullSimplify∘Together∘Expand`) on that tree and both reach 0; that is genuine cross-engine corroboration of the same identity rather than one echoing the other's intermediate algebra. Corresponding excerpts: `R_target = sp.simplify(NQ*A/(beta0*kappa0_sq))` (py L26) vs `rTarget = FullSimplify[nQ*aNorm/(beta0*kappa0Sq), ...]` (wl L35); `expect_zero("Pi_tr/C_mix - alpha_req/alpha_mix", sp.simplify(Pi_tr/C_mix - alpha_req/alpha_mix))` (py L41) vs `expectZero["Pi_tr/C_mix - alpha_req/alpha_mix", piTr/cMix - alphaReq/alphaMix]` (wl L50).

## Engine cross-check

Both engines emit identical final forms (modulo symbol spelling): `R_target = pi^2*A*NQ/(8*beta0)` (py output L5) = `(aNorm*nQ*Pi^2)/(8*beta0)` (wl output L5); `Pi_tr = NQ*alpha_req/beta0` (py L8) = `(alphaReq*nQ)/beta0` (wl L8); `NQ_selected = beta0*mhat^2*s_minus/lam_minus` (py L13) = `(beta0*mhat^2*sMinus)/lamMinus` (wl L16); `Pi_sel = alpha_req*mhat^2*s_minus/lam_minus` (py L14) = `(alphaReq*mhat^2*sMinus)/lamMinus` (wl L17). All eight residual checks report `= 0` in both transcripts and the Mathematica run prints `PASS:` for each plus "Stage 085 Mathematica audit passed." `engines_agree: true`.

## Verdict justification

Verdict is **clean**. I read the paper card, the notes, and the appendix row first, built the model (three boxed deliverables: product identities, ratio cancellation, reduced demand law, plus the spectral-form intermediate), then attacked the scripts. Every paper deliverable maps to a non-tautological, well-anchored assertion in both engines. The load-bearing cancellation (A1/A2/A3) genuinely tests the `8/pi^2` constant against the `8/(pi^2 A)` loading factors — it fails if `kappa0^2` is wrong — so it is not an identity-by-construction. The `zeta_req` checks confirm the `NQ/beta0` factor cancels through the full rational demand expression, the `rho_alpha` reduction, and the unblocked limit. Symbol assumptions are physically justified; the sole literal `kappa0^2=8/pi^2` is anchored in the notes; outputs are fresh; self-labels are correct. The Mathematica script is a parallel encoding of an irreducibly short algebraic identity (no distinct second route exists) and corroborates via an independent simplifier, so it is not a transliteration finding. Nothing broke.

## Self-test notes

Checked: (1) no `sp.diff`/`D[]` in this stage (pure algebraic identities), so the zero-derivative trap is N/A. (2) No integrals/parity considerations — N/A. (3) Trivial-case substitution: the `8/pi^2` cancellation reduces to exact 0 only with the correct constant (verified the residual would be nonzero for a wrong `kappa0^2`), and the unblocked-limit reduction `(alpha_req-alpha_mix)/alpha_mix = alpha_req/alpha_mix-1` is a genuine simplification. (4) Path specs N/A (no missing-script findings). (5) Paper round-trip: no fix prescribed, so no new misalignment introduced. Conclusion: all probes failed to break the scripts; clean.

## Value Reconciliation (pass-2 augmentation)

Enumerating every RESULT/deliverable value the scripts emit (from the `.py`/`.wl` source and the committed `.txt` outputs), excluding pass/fail flags and residual-zero scaffolding. All emitted deliverables are symbolic closed forms (this stage produces no numeric figures of merit beyond the literal constant `kappa0^2`).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `R_target = pi^2*A*NQ/(8*beta0)` | py L33 / py-out L5; wl L42 / wl-out L5 | notes §1 L52 `R_target = N_Q A/(beta_0 kappa_0^2)` (with kappa0^2=8/pi^2) | MATCH |
| `G_tr = 8*alpha_req/(pi^2*A)` | py L34 / py-out L6; wl L43 / wl-out L6 | notes §1 L60 `G_tr = 8 alpha_req/(pi^2 A)` | MATCH |
| `M_mix = 8*alpha_mix/(pi^2*A)` | py L35 / py-out L7; wl L44 / wl-out L7 | notes §1 L62 `M_mix = 8 alpha_mix/(pi^2 A)` | MATCH |
| `Pi_tr = NQ*alpha_req/beta0` | py L36 / py-out L8; wl L45 / wl-out L8 | tex eq:products L16–17; notes §1 L72 | MATCH |
| `C_mix = NQ*alpha_mix/beta0` | py L37 / py-out L9; wl L46 / wl-out L9 | tex eq:products L18–19; notes §1 L74 | MATCH |
| `Pi_tr/C_mix = alpha_req/alpha_mix` | py-out L12; wl-out L14 (asserted A3/B3) | tex eq:ratio L24–27; notes §1 L78 | MATCH |
| `NQ_selected = beta0*mhat^2*s_minus/lam_minus` | py L48 / py-out L13; wl L56 / wl-out L16 | notes §2 L88–92 `mhat_-^2 P_{0,-}=N_Q^{target}`, `P_{0,-}=beta_0 s_-/lam_-` | MATCH |
| `Pi_sel = alpha_req*mhat^2*s_minus/lam_minus` | py L49 / py-out L14; wl L57 / wl-out L17 | notes §2 L100 `Pi_tr = mhat_-^2 (s_-/lam_-) alpha_req` | MATCH |
| `C_sel = alpha_mix*mhat^2*s_minus/lam_minus` | py L50 / py-out L15; wl L58 / wl-out L18 | notes §2 L102 `C_mix = mhat_-^2 (s_-/lam_-) alpha_mix` | MATCH |
| `zeta_req = (alpha_req-alpha_mix)/[alpha_mix-eps_blk(2 alpha_mix-alpha_req)]` | py L74 (final ledger) / py-out L26 | notes §3 L119–120; tex eq:zeta-rho L32–34 (rho form) | MATCH |
| `zeta_req(rho,eps) = (rho_alpha-1)/[1-eps_blk(2-rho_alpha)]` | py L63 / wl L64 (asserted A7/B7) | tex eq:zeta-rho L32–34; notes §5 L181–182 | MATCH |
| unblocked: `zeta_req = rho_alpha-1` | py L67 / wl L68 (asserted A8/B8) | tex L36 "In the unblocked limit, zeta_req=rho_alpha-1"; notes §3 L134 | MATCH |
| `kappa0^2 = 8/pi^2` (literal constant) | py L24; wl L34 | notes §1 L57 "kappa_0^2 = 8/pi^2" | MATCH |

INTERNAL scaffolding (accounted for, no finding expected in prose): `expect_zero`/`expectZero` residual values (all `0`), the `banner`/`pass`/`fail`/`fmt` helpers, and the final-ledger print lines (py L70–74) which are restatements of the boxed results already in the table.

reconciliation: complete; 13 values checked, 0 misaligned
