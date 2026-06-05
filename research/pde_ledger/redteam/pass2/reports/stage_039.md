---
unit_id: 039
batch: III.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage039_split_u_sector.md]
  paper_appendix: present
---

# Audit unit 039 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_039.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage039_split_u_sector.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 039 at line 56)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage039_split_u_sector_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage039_split_u_sector_mathematica_audit.txt`

## What the paper claims

Stage 039 turns on the first axial stiffness of the internal `U` sector (`K_U1 = K_U(1+delta_U)`, `delta_U = pi^2 T_U/(L^2 K_U)`) and asks what survives. `\stagefield{Output}` states the deliverables verbatim: "The split placement law \eqref{eq:app-stage039-split-placement}, the direction factor \eqref{eq:app-stage039-RU}, and the collinearity theorem \eqref{eq:app-stage039-collinearity}." The boxed equations give four closed forms: (1) the split placement pair `delta_split = [delta_0 + eps_eta delta_U/(1+delta_U)]/(1-eps_eta)` and `eps_W,split = eps_W[1 - (2/11) delta_U/(1+delta_U)]`; (2) the direction factor `R_U = [1+rho_0/(1+delta_U)]/(1+rho_0)`; (3) the direction-splitting invariant `D_dir = -kappa_0 kappa_1 g_W rho_0 delta_U/(1+delta_U)`; (4) the collinearity theorem `D_dir = 0 ⟺ delta_U = 0 or rho_0 = 0`. The notes add two further exact deliverables: the split placement map (`M_mix^(splitU)`, `R_target^(splitU)`) and the **surviving exact product law** `R_target^(splitU) M_mix^(splitU) = 8 Lambda (1 - eps_W,split)/pi^2` (notes §5, lines 207–218), plus the overlap constants `sigma = 88/(9 pi^2)`, `lambda_0 = 2/9` (notes lines 67, 69). The appendix row (line 56) summarizes: "Split placement, direction-splitting invariant, and collinearity iff condition."

## What the script claims to verify

Both scripts (a) pin the overlap data `kappa_0, kappa_1, sigma, lambda_0`; (b) check the direct softenings `A0, A1` reduce to `A0 = K_eta_eff(1-eps_eta)/mu_eta`, `A1 = A0(1+delta_split)` with the postulated `delta_split` (sympy postulates and checks `A1`; `.wl` derives `delta_split` from `A1/A0-1` and checks it equals the postulated form); (c) check `eps_W,split` derived from `S_U = kappa_0^2/K_U + kappa_1^2/K_U1` equals the `2/11` split formula; (d) compute the loading vector `z0, z1, R_U`, check the `z1/z0` ratio identity, derive `D_dir = kappa_0 z1 - kappa_1 z0` and assert it equals the postulated closed form, then prove the collinearity theorem (both if-legs by substitution + the only-if leg by factoring `rho_0 delta_U` out of the numerator and asserting the residual is nonzero and free of `rho_0, delta_U`); (e) **print** (no assertion) the placement map `M_mix, R_target` and the product; (f) print the small-`delta_U` series.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `delta_split` (split-placement box) | sympy A1 check (l.87), `.wl` derived==postulated (l.71) | match |
| `eps_W,split` (split-placement box) | both engines assert derived==split formula (sympy l.98, `.wl` l.83) | match |
| `R_U` direction factor (eq-RU, boxed Output) | only printed; nearest assertion is z1/z0 ratio (sympy l.114, `.wl` l.95) which is **tautological** (see F1) | partial→missing |
| `D_dir` invariant (eq-Ddir, boxed Output) | both assert `D_dir - postulated == 0` (sympy l.122, `.wl` l.105) | match |
| collinearity theorem (eq-collinearity, boxed Output) | if-legs + only-if factoring, both engines (sympy l.125–136, `.wl` l.108–119) | match |
| split placement map `M_mix, R_target` (notes §5) | only printed, no assertion (sympy l.149–151, `.wl` l.129–131) | partial |
| **exact product law** `R_t M_mix = 8 Lambda(1-eps_W,split)/pi^2` (notes §5, l.217) | only `print("product = …")`; **never asserted** (sympy l.151, `.wl` l.131) | missing |
| `sigma = 88/(9pi^2)`, `lambda_0 = 2/9` (notes) | computed and printed, both engines | match |

Dominant pattern: the core split/direction/collinearity claims are faithfully and non-tautologically checked, but two stated deliverables (`R_U` as a falsifiable check, and the surviving product law) are only printed → `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 86 | `A0.subs - A0_expected == 0` | A0 reduction (delta_split scaffold) | yes |
| A2 | sympy | 87 | `A1.subs - A1_expected == 0` | delta_split | yes |
| A3 | sympy | 98 | `eps_W_direct.subs - eps_W_split == 0` | eps_W,split (2/11) | yes |
| A4 | sympy | 114 | `z1(1+rho0) - (k1/k0)z0(1+rho0/(1+dU)) == 0` | R_U (intended) | **no — tautological** |
| A5 | sympy | 122 | `D_dir - D_dir_expected == 0` | D_dir invariant | yes |
| A6 | sympy | 125 | `D_dir(deltaU=0) == 0` | collinearity if-leg | yes |
| A7 | sympy | 126 | `D_dir(rho0=0) == 0` | collinearity if-leg | yes |
| A8 | sympy | 131–136 | numerator/`(rho0 deltaU)` free of rho0,deltaU and nonzero | collinearity only-if | yes |
| B1 | math | 69 | `(a0/.) - a0Expected == 0` | A0 reduction | yes |
| B2 | math | 70 | `(a1/.) - a1Expected == 0` | delta_split | yes |
| B3 | math | 71 | `deltaSplitDerived - deltaSplitPostulated == 0` | delta_split | yes |
| B4 | math | 83 | `epsWSplitDerived - postulated == 0` | eps_W,split (2/11) | yes |
| B5 | math | 95 | z1/z0 ratio identity | R_U (intended) | **no — tautological** |
| B6 | math | 105 | `dDirDerived - dDirPostulated == 0` | D_dir invariant | yes |
| B7 | math | 108 | `dDir(deltaU=0) == 0` | collinearity if-leg | yes |
| B8 | math | 109 | `dDir(rho0=0) == 0` | collinearity if-leg | yes |
| B9 | math | 113–119 | numerator/`(rho0 deltaU)` free + nonzero | collinearity only-if | yes |
| — | both | M_mix/R_target/product | Print only, no assertion | placement map + product law | n/a |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:107-117`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:88-98`

**What's wrong:**
The `z1/z0` ratio assertion is the only check intended to exercise the boxed direction factor `R_U`. But `z0` and `z1` are *defined* with exactly `kappa0` and `kappa1` as prefactors:
- sympy l.107–108: `z0 = kappa0*g_W*(1+rho0)`, `z1 = kappa1*g_W*(1+rho0/(1+deltaU))`
- the assertion l.114–117: `z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU)) == 0`.

Substituting the definitions: `z1*(1+rho0) = kappa1 g_W (1+rho0/(1+deltaU))(1+rho0)` and `(kappa1/kappa0) z0 (1+rho0/(1+deltaU)) = (kappa1/kappa0) kappa0 g_W (1+rho0)(1+rho0/(1+deltaU)) = kappa1 g_W (1+rho0)(1+rho0/(1+deltaU))`. The two sides are **identically equal by construction** — the `kappa0` in `z0` cancels the `kappa0` in the denominator, and the `(1+rho0)`/`(1+rho0/(1+deltaU))` factors match term-for-term. The assertion cannot fail for any value of `rho0, deltaU, g_W`, so it does not verify the paper's claim `z1/z0 = (kappa1/kappa0) R_U` — it merely re-states the definitions of `z0, z1`. The `.wl` (l.95–98) is the identical tautology. `R_U` is printed (sympy l.113, `.wl` l.90) but never falsifiably checked.

**Why this matters:**
`R_U` is a boxed deliverable named explicitly in `\stagefield{Output}` ("the direction factor \eqref{eq:app-stage039-RU}"). The intended check is that the loading vector's component ratio equals `(kappa1/kappa0) R_U` with `R_U = [1+rho0/(1+deltaU)]/(1+rho0)` — i.e., that `R_U` correctly captures the direction rotation. The current assertion proves nothing about that; a wrong `R_U` formula (e.g., a sign error or wrong denominator) would not be caught because `R_U` is not referenced in the assertion at all.

**Required change:**
Replace the tautological ratio check with a falsifiable one that ties the component ratio to `R_U`. Add (both engines): assert `simplify(z1/z0 - (kappa1/kappa0)*R_U) == 0` where `R_U` is the independently-defined `(1+rho0/(1+deltaU))/(1+rho0)` (sympy l.109 / `.wl` l.90). This is non-tautological because `R_U` is defined from `rho0, deltaU` directly (not via `z0, z1`), so the check confirms the *paper's* closed form for `R_U` reproduces the loading-vector ratio. (Keep the existing print of `R_U`.)

**Verification:**
After fix, sympy l.~114 and `.wl` l.~95 should show a check named e.g. `"z1/z0 - (kappa1/kappa0)*R_U"` evaluating to `0`; the residual must reduce to 0 only because `R_U` equals the boxed form, and substituting a deliberately-wrong `R_U` must make it nonzero.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:147,151`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:127,131`

**What's wrong:**
The notes (§5, lines 215–218) and the appendix framing make the *surviving exact product law* a stated deliverable: "the exact product law survives: `R_target^(splitU) M_mix^(splitU) = 8 Lambda (1 - eps_W,split)/pi^2`." The scripts **compute** the product (sympy l.147 `product = sp.simplify(M_mix_split * R_target_split)`; `.wl` l.127) and **print** it (sympy l.151, `.wl` l.131), but there is **no assertion** that it equals `8 Lambda (1 - eps_W,split)/pi^2`. The printed value (sympy output l.45) is `8*Lambda*(11*deltaU - eps_W*(9*deltaU+11) + 11)/(11*pi**2*(deltaU+1))`, which does equal `8 Lambda(1-eps_W,split)/pi^2` — but nothing checks that, so a regression would pass silently. Likewise `M_mix_split` and `R_target_split` themselves are only printed (no check against the notes' closed forms).

**Why this matters:**
"The factorization survives" is the headline of notes §5 and §6 (the whole point that scalar placement is unbroken while direction splits). It is currently unverified — a falsifiable claim reduced to a `print`. The audit's verdict cannot certify this deliverable.

**Required change:**
Add to both engines a falsifiable check: `expect_zero("product law survives", product - 8*Lambda*(1 - eps_W_split)/pi**2)` (sympy after l.151; `.wl` `expectZero["product law survives", product - 8 lambda (1 - epsWSplit)/Pi^2]` after l.131), using the already-defined `eps_W_split`/`epsWSplit`. Optionally also assert `M_mix_split` and `R_target_split` equal the notes' boxed forms, but the product law is the load-bearing one.

**Verification:**
After fix, both outputs show `product law survives = 0` / `PASS`. The check is non-tautological: `product` is built from independently-`simplify`d `M_mix_split` and `R_target_split`, while the RHS is the notes' closed form expressed through `eps_W_split`.

### F3 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage039_split_u_sector_sympy_audit.txt` (mtime 2026-05-26 02:03:02)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage039_split_u_sector_mathematica_audit.txt` (mtime 2026-05-26 02:03:07)
- scripts mtime: both 2026-06-03 15:59:11 (newer than outputs)

**What's wrong:**
Both saved transcripts predate their scripts by ~8 days, so they are stale. The staleness is visible in label drift: the sympy transcript banner reads "STAGE 22" (l.3) / "STAGE 22 THEOREM LEDGER" (l.57) and sub-banners "22.1…22.5", while the current sympy *script* banner emits "STAGE 39" (l.62) — i.e., the script's top banner was already updated but the committed output was not regenerated. The Mathematica transcript banner reads "STAGE 022" (l.4) while the current `.wl` emits "STAGE 039" (l.33). The numeric *results* in the stale transcripts still match what the current scripts produce (verified by re-deriving `S_U`, `eps_W_split`, `D_dir`, product by hand — all agree), so the staleness is label/banner-level only, not a result discrepancy. This is the known deferred SCRIPT/OUTPUT-band numbering drift (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING) — the *source* `.py` also still carries stale self-labels in its docstring (l.3 "Stage 22"), sub-banners (l.68/89/103/140/153 "22.x"), and "Stage-21" cross-refs (l.15, 142, 178); the `.wl` source labels are already clean (sections "1."–"5.", banner "STAGE 039").

**Why this matters:**
The orchestrator's independent re-run refreshes these transcripts; once F1/F2 land, both outputs regenerate. The label drift itself is cosmetic and routed to the dedicated numbering pass, not fixed ad hoc here (per project policy: content-keyed, never offset-sweep). Flagging for completeness and to confirm result content is unaffected.

**Required change:**
No content edit required for the staleness per se — the orchestrator re-run regenerates the `.txt`. The stale self-labels in the `.py` source ("Stage 22"/"22.x"/"Stage-21") are left to the deferred SCRIPT/OUTPUT-band numbering pass, not corrected piecemeal in this directive.

**Verification:**
After the orchestrator re-runs sympy+mathematica for unit 039, both transcript mtimes exceed the script mtimes and the banners read "STAGE 39"/"STAGE 039".

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. The two engines take genuinely different routes on the load-bearing `delta_split`: SymPy *postulates* `delta_split` (l.78), builds `A1_expected = A0_expected(1+delta_split)`, and checks the physics-defined `A1` matches (l.87) — i.e., it validates the postulated form by forward substitution. Mathematica instead *derives* `delta_split` from the physics: `deltaSplitDerived = a1Direct/a0Expected - 1` (l.60) and then asserts the derived value equals the postulated closed form (l.71). These are opposite directions (postulate-and-confirm vs. derive-and-match), not a line-by-line port. The `.wl` also adds explicit `deltaSplit derived matches postulated` and `eps_W_split derived matches postulated` checks (l.71, 83) that the SymPy script folds into its A1/eps_W substitutions. The shared scaffolding (`kappa0=2√2/π`, `kappa1=-4/(3π)`, banner helpers) is unavoidable common premise, not echoed algebra. No `mathematica_transliteration` finding.

Note (non-finding): the `.wl` has a benign forward-reference ordering quirk — `deltaSplitDerived` (l.60) and `a1Direct` (l.59) reference `a0Expected` before it is assigned (l.63). Because WL `Set` stores the still-symbolic `a0Expected` and re-evaluates it once defined, the printed/asserted results resolve correctly (output l.15 shows `a0Expected` correctly substituted). It produces the right answer here and the assertions pass, so it is not a current defect — but it is fragile style. Not flagged as a finding since no assertion is weakened or made to pass falsely by it.

## Engine cross-check

The engines agree on every shared quantity (comparing sympy output vs. mathematica output):
- `sigma = 88/(9π²)`, `lambda0 = 2/9` — agree (sympy l.7–8, math l.7–8).
- `S_U = 8(9deltaU+11)/(9π²K_U(1+deltaU))` — agree (sympy l.22, math l.26).
- `eps_W_split = eps_W(9deltaU+11)/(11(1+deltaU))` — agree (sympy l.23, math l.27).
- `D_dir = 8√2 c_etaW deltaU rho0/(3π² √(muW muEta)(1+deltaU))` — agree (sympy l.33, math l.39).
- `Numerator(D_dir)/(rho0 deltaU) = 8√2 c_etaW` — agree (sympy l.37, math l.46).
- product, M_mix, R_target, all series — agree term-for-term (sympy l.43–54, math l.52–63).
All `expect_zero`/`expectZero` checks report `0`/`PASS` in both transcripts. No `engine_disagreement`.

## Verdict justification

The core physics of Stage 039 — the split anisotropy `delta_split`, the `2/11` split blocking ratio `eps_W,split`, the direction-splitting invariant `D_dir`, and the collinearity iff theorem (both if-legs plus a genuinely non-trivial only-if leg) — is faithfully and non-tautologically verified by both engines, which agree throughout and derive `delta_split` by independent routes. I tried to break the `2/11` factor (re-derived `S_U` and the `eps_W,split` substitution by hand: exact), the `D_dir` cancellation (exact), and the only-if factoring (genuinely confirms `D_dir ∝ rho0·deltaU`): all hold. The verdict is `findings`, not `clean`, for two real gaps: (F1) the only check intended to exercise the boxed direction factor `R_U` is tautological — `z0, z1` are defined with `kappa0, kappa1` prefactors so the ratio identity holds by construction and never references `R_U`; and (F2) the surviving exact product law (a stated notes deliverable, the headline "factorization survives" claim) is computed and printed but never asserted. Plus (F3) both transcripts are stale (label-only; results still match). Both F1 and F2 are script-side fixable (add falsifiable assertions tying the checks to the paper's `R_U` and product-law forms); no `paper_misalignment` (every value the script emits matches the paper/notes), no `UNFIXABLE`, no `CRITICAL_DOWNSTREAM`. I confirm I read the stage card, the notes, and the appendix row before auditing the scripts.

## Self-test notes

I checked: (1) variable independence — F1's proposed `z1/z0 - (kappa1/kappa0)*R_U` references `R_U` defined from `rho0, deltaU` (not from `z0, z1`), so it is non-tautological and reduces to 0 only when `R_U` equals the boxed form (verified by substituting the definitions: `z1/z0 = (kappa1/kappa0)(1+rho0/(1+deltaU))/(1+rho0) = (kappa1/kappa0)R_U`, exact). (2) Trivial-case: F2's `product - 8Lambda(1-eps_W_split)/pi^2` — I substituted `deltaU→0` mentally: `product→8Lambda(1-eps_W)/pi^2` and RHS→same, residual 0; for `deltaU≠0` the printed `product` (output l.45) equals the RHS expanded, so the assert_zero is satisfiable and non-trivial. (3) No new derivatives or symmetric-domain integrals introduced. (4) Path specs: F1/F2 edit existing files in `scripts/` and `mathematica/` — no new-script paths. (5) Paper round-trip: both fixes assert the paper/notes' own closed forms (`R_U`, product law), introducing no new constant, so no new `paper_misalignment`.

## Value Reconciliation (pass-2 augmentation)

Every labeled RESULT the scripts emit, located in `.tex`/`.md`:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `kappa0 = 2√2/π` | py l.57, out l.5; wl l.43, out l.5 | notes l.57–60 (basis), implied | MATCH |
| `kappa1 = -4/(3π)` | py l.58, out l.5; wl l.44, out l.5 | notes (basis §1) | MATCH |
| `sigma = 88/(9π²)` | py l.59, out l.7; wl l.45, out l.6 | notes l.67 `sigma = … = 88/(9 pi^2)` | MATCH |
| `lambda_0 = 2/9` | py l.60, out l.8; wl l.46, out l.8 | notes l.69 `lambda_0 := … = 2/9` | MATCH |
| `delta_split = [delta0+eps_eta deltaU/(1+deltaU)]/(1-eps_eta)` | py l.78/84, out l.15; wl out l.15 | tex l.26–27 (boxed); notes l.97 | MATCH |
| `eps_W_split = eps_W[1-(2/11)deltaU/(1+deltaU)]` | py l.94, out l.23; wl out l.27 | tex l.29–30 (boxed); notes l.123 | MATCH |
| `R_U = [1+rho0/(1+deltaU)]/(1+rho0)` | py l.109, out l.31; wl out l.36 | tex l.35–36 (boxed); notes l.157 | MATCH |
| `D_dir = -kappa0 kappa1 g_W rho0 deltaU/(1+deltaU)` | py l.120, out l.33; wl out l.39 | tex l.41–42 (boxed); notes l.172 | MATCH |
| `M_mix^(splitU) = 8 Z_W(1+rho0)²/[π²(1-eps_eta)(1-eps_W_split)]` | py l.145, out l.43; wl out l.52 | notes l.207–209 | MATCH |
| `R_target^(splitU) = Lambda(1-eps_eta)(1-eps_W_split)²/[Z_W(1+rho0)²]` | py l.146, out l.44; wl out l.53 | notes l.211–213 | MATCH |
| product law `R_t M_mix = 8 Lambda(1-eps_W_split)/π²` | py l.147, out l.45; wl out l.54 | notes l.217–218 | MATCH (value); **see F2 — printed not asserted** |
| series `delta_split`, `eps_W_split`, `R_U`, M/R ratios | py l.163–167, out l.50–54; wl out l.59–63 | notes l.103, 129, 187, 225–227 | MATCH |

INTERNAL (scaffolding, no prose expected, no finding): `A0`, `A1`, `A0_expected`, `A1_expected`, `S_U` (intermediate), `eps_W_direct`, `eps_eta_def`, `g_W`, `z0`, `z1`, `K_U1`, `a0`/`a1`/`a1Direct`/`deltaSplitDerived`/`deltaSplitPostulated` (`.wl` intermediates), `Numerator(D_dir)/(rho0 deltaU) = 8√2 c_etaW` (only-if residual), all `expect_zero`/`expectZero` residual flags.

Note: `S_U`, `M_mix`, `R_target` are emitted as labeled results AND carried in the notes, so they are listed in the main table as MATCH; their *use* inside the unasserted product is the substance of F2.

reconciliation: complete; 12 deliverable values checked, 0 misaligned (every emitted value matches the .tex card and/or notes). The two findings F1/F2 are NOT value mismatches — they are missing/tautological *assertions* over values that are themselves correct.
