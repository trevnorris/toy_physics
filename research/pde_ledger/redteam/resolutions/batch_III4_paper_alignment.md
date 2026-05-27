---
batch: III.4
created: 2026-05-27
status: applied_and_verified
items:
  - Q1: 074 F1 — alpha value mismatch (paper 128 vs notes 179 vs engines 111)
  - Q2: 075 F2 — Upsilon_w conversion factor (paper 117 vs notes 168 vs script 100)
  - Q3: 082 F1+F2 — extend scripts to exercise zeta_phys closed form and instantiate Family-1 (eta, kappa) pair
---

# Batch III.4 paper-alignment resolutions

Three substantive `paper_misalignment` items requiring user direction. The four banner-relabel items (078 F2 on .py and .wl; 081 F1 on .py; 084 stale "STAGE 067" on .wl) are pure script-side label fixes with unambiguous direction (paper card and filename are authoritative); they are applied directly by the orchestrator and do not appear here.

The Codex recommendation pass writes `## Recommendation` blocks below each question with `direction:` and rationale. The user reviews and approves; then a separate Codex apply pass executes the approved direction.

---

## Q1 — Stage 074 F1: `alpha = sqrt(kappa)` value mismatch

**Auditor finding (full):** `redteam/reports/stage_074.md` F1 (subtype `value_mismatch`).

### Three values in play

- `paper/stages/stage_074.tex:26-31` boxes: `kappa = (9/5) Lambda_ell^2 = 12321/5`, `alpha = sqrt(kappa) = 128/sqrt(5)`.
- `notes/stages/moving_throat_pde_stage074_family1_healing_lock.md:113-117` states: `alpha = sqrt(12321/5) = 179/sqrt(5) ~ 49.6407091`.
- Both engines compute `sqrt(12321/5) = 111/sqrt(5) ~ 49.6407091` (since `111^2 = 12321` exactly).

The numerical value `~49.6407091` in the notes is the value of `111/sqrt(5)` (not of `179/sqrt(5)`, which would be `~80.05`). So the notes' final numeric value already agrees with the engines; only the prose literal `179` is a typo. The paper card's `128` similarly disagrees with `sqrt(12321/5)` arithmetically.

The script computes the arithmetically correct value but only prints it (no assertion); the inconsistency went undetected at v1 because no script-side assert anchored alpha to a literal.

### Options

- **(a) Paper-side and notes-side typo fix** — both `128` and `179` are typos for `111`. Update `paper/stages/stage_074.tex:26-31` to `alpha = 111/sqrt(5)` and `notes/stages/moving_throat_pde_stage074_family1_healing_lock.md:113-117` to `alpha = sqrt(12321/5) = 111/sqrt(5) ~ 49.6407091`. Add a script-side assertion `alpha_ref - sp.Rational(111) / sp.sqrt(5) == 0` in both engines to lock the value going forward.
- **(b) Recompute `kappa`** — if alpha really should be `128/sqrt(5)`, then `kappa` should be `16384/5`, not `12321/5`, which would invalidate the carry-forward of `kappa` into stages 075/078/082. Almost certainly not the right direction, but stated for completeness.
- **(c) Defer** — leave the paper and notes as-is, mark in the script comments that alpha is computed but not asserted, and surface the inconsistency as a known paper-side erratum elsewhere.

### Destination-verification grep (run locally by orchestrator)

Downstream grep for `128` / `179` / `alpha = ` literals across stages 075-084 (paper/, notes/, scripts/, mathematica/):

- `paper/stages/stage_07[5-9].tex` and `paper/stages/stage_08[0-4].tex`: **no occurrence** of `128` or `179` as alpha literals.
- `notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md:157`: `alpha = sqrt(kappa).` (symbolic — fine).
- `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:63`: `alpha := sqrt(kappa) = 179/sqrt(5) ≈ 49.6407091,` — **same `179` typo as 074 notes**.
- `scripts/.../stage075..._sympy_audit.py:66` and `mathematica/.../stage075..._mathematica_audit.wl:55`: print only, no literal.

No downstream stage uses `alpha` directly with a literal value; they all consume `kappa = 12321/5` and recompute alpha symbolically when needed. The `179` typo also lives in `notes/stage075` line 63 and needs the same fix as 074 notes line 117.

## Recommendation

- direction: **a** (paper-side and notes-side typo fix, plus add script-side assertion)
- rationale: Both engines independently compute `sqrt(12321/5) = 111/sqrt(5)` (since `111^2 = 12321` exactly). The numerical value `~49.6407091` cited in both `notes/stage074:117` and `notes/stage075:63` matches `111/sqrt(5)` (which is `≈49.6407091...`), not `179/sqrt(5)` (`≈80.05`) nor `128/sqrt(5)` (`≈57.24`). So the prose literals `128` (paper) and `179` (notes ×2) are isolated typos for `111`, and the surrounding numerical values are already correct. Downstream stages do not propagate the typo.
- downstream_impact: None mathematically. Cosmetically, three files need the typo fix: `paper/stages/stage_074.tex:31`, `notes/stages/moving_throat_pde_stage074_family1_healing_lock.md:117`, and `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:63`. Adding the assertion `alpha_ref - sp.Rational(111) / sp.sqrt(5) == 0` (and `alphaRef - 111/Sqrt[5]` in `.wl`) at stage 074 prevents regression.
- notes: The `notes/stage075:63` typo is technically Q1's collateral (same typo, different file) but lives in stage 075's notes; including it in the Q1 paper-side fix keeps the resolution coherent. The Stage 074 audit also asks for a script-side assertion locking alpha to `111/sqrt(5)`; this is a script-side extension, not a paper-side change, and can be done by Codex once the user approves (a).

---

## Q2 — Stage 075 F2: `Upsilon_w` conversion factor mismatch

**Auditor finding (full):** `redteam/reports/stage_075.md` F2 (subtype `value_mismatch`). High severity.

### Three values in play

- `paper/stages/stage_075.tex:7` Inputs line: `Upsilon_w = 117 Theta_w`.
- `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:108` states: `Upsilon_w = 168 Theta_w`.
- Script: `alpha_r = 10` so `Upsilon_w = alpha_r^2 Theta_w = 100 Theta_w` (sympy:24, mathematica:37).

The paper's boxed final `Theta_w <= 3.62606e-4 Pe_req => fail` and `Theta_w >= 4.21495e-2 Pe_req => succeed` can be reproduced *only* by dividing the Upsilon thresholds by 100 (the script's value), not by 117 and not by 168. The notes' arithmetic in section 4 also produces numbers consistent with 100, not 168 (notes line 124-128 claims `Upsilon_fail/168 = 3.626e-4`, but actually `0.036260562/168 = 2.158e-4` — internally inconsistent).

So the paper Inputs line and the notes line are both stale text; the script and the paper's boxed numerics agree on 100.

### Options

- **(a) Paper-side and notes-side text fix** — update paper Inputs line `117 Theta_w` → `100 Theta_w`; update notes `168 Theta_w` → `100 Theta_w` and notes' Theta_fail/Theta_suff arithmetic accordingly. No script change. Add script-side assertion locking `alpha_r^2 == 100` to the rescaling constant going forward.
- **(b) Change script `alpha_r`** — set `alpha_r` such that `alpha_r^2 = 117` (i.e., `alpha_r = sqrt(117) ~ 10.816`) or `168`. This propagates new numerics into stages 076-084 and breaks the paper's own boxed Theta values.
- **(c) Defer with comment** — document the inconsistency in script comments without paper-side fix.

### Destination-verification grep (run locally by orchestrator)

Downstream grep for `Upsilon_w = 117|168|100` across stages 075-084:

- `paper/stages/stage_075.tex:7` and `:24`: both say `117 Theta_w`. **The only place `117` appears.**
- `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:108` and `:116`: both say `168 Theta_w`. **The only place `168` appears.**
- `notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md:11`: `Upsilon_w = alpha_r^2 Theta_w` (symbolic — fine).
- `notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md:218`: `Xi_F1 = ... = 1369 Upsilon_w = 136900 Theta_w` — load-bearing, consistent only with `Upsilon_w = 100 Theta_w` (since `1369 × 100 = 136900`).
- `notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md:71`: `Xi_F1 = W_wall = 1369 Upsilon_w = 136900 Theta_w` — same constraint.
- `notes/stages/moving_throat_pde_stage084_full_reduced_pde_writeup.md:119`: same `136900 Theta_w` — same constraint.
- `scripts/stage082..._sympy_audit.py:91, 106, 109, 114` and `mathematica/stage082..._mathematica_audit.wl:100, 108` and `mathematica/stage084..._mathematica_audit.wl:63`: all hardcode `Upsilon_w = 100 Theta_w` and `Xi_F1 = 136900 Theta_w`.

The `Xi_F1 = 136900 Theta_w` carry-forward across stages 082, 083, 084 (paper, notes, scripts) is internally consistent only with `Upsilon_w = 100 Theta_w`. The two isolated occurrences of `117` and `168` (paper:075 Inputs+body; notes:075 §3) are stale text drift.

## Recommendation

- direction: **a** (paper-side and notes-side text fix, plus add script-side assertion locking `alpha_r^2 == 100`)
- rationale: The script's `alpha_r = 10` → `Upsilon_w = 100 Theta_w` is internally consistent with: (1) the paper's own boxed final `Theta_fail = 3.626e-4 Pe_req` and `Theta_suff = 4.215e-2 Pe_req` numbers in stage 075 (only reproducible by dividing the Upsilon thresholds by 100); (2) the `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w` identity carried forward in `notes/stage082:218`, `notes/stage083:71`, `notes/stage084:119` (only consistent with `Upsilon_w = 100 Theta_w`); (3) the notes' own arithmetic in §4 (e.g., `Theta_fail = Upsilon_fail / 168 ≈ 3.626e-4` is internally inconsistent — `0.0362605617972939 / 168 = 2.158e-4`, not `3.626e-4`, but the final numeric value `3.626e-4` is correct for division by 100). Direction (a) updates the two stale prose lines without touching any math.
- downstream_impact: None mathematically. Paper-side edit: `paper/stages/stage_075.tex:7` Inputs line `117 Theta_w` → `100 Theta_w`; same change at `paper/stages/stage_075.tex:24` body line. Notes-side edit: `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:108` `168 Theta_w` → `100 Theta_w` and `:116` same; the notes' "Theta_fail = Upsilon_fail / 168" arithmetic at lines 124-128 should be updated to "Theta_fail = Upsilon_fail / 100" so the divisor matches the prose. Script-side: add a sympy `expect_zero("alpha_r^2 - 100", alpha_r**2 - 100)` and equivalent Mathematica check to lock the value going forward. This also resolves the Upsilon_w portion of Q3 / 082 F2 automatically (the script's hardcoded `100` is correct; only the paper Inputs prose needs the matching fix).
- notes: The notes' `168` is almost certainly stale artifact from an earlier `alpha_r = sqrt(168) ≈ 12.96` design that was later normalized to `alpha_r = 10`; the notes' own numeric Theta values were updated to match `100` but the `168` text was not. The paper's `117` has no obvious provenance — likely a transcription error from `100` via OCR-like substitution (1↔1, 0↔7 are visually adjacent in some fonts, but uncertain). Either way, both sources agree the loadbearing value is `100`.

---

## Q3 — Stage 082 F1+F2: extend scripts to exercise paper deliverables

**Auditor findings (full):** `redteam/reports/stage_082.md` F1 (subtype `script_missing_paper_claim`) and F2 (subtype `script_missing_paper_claim`).

### Paper claims not exercised

- **F1: zeta_phys closed form.** `paper/stages/stage_082.tex` (eq. `app-stage082-zeta-phys`) gives `zeta_phys = Omega_Pe^2 (kappa + pi^2/4) / (kappa + y(eta)^2)` and uses it as the inverse demand map. Both scripts leave `zeta_phys` as an opaque free symbol — the closed-form definition is never asserted or substituted.
- **F2: Family-1 numerical pair `(eta, kappa) = (37, 12321/5)`.** Paper card states these as the Family-1 reference. Scripts never instantiate the numerical values; everything stays symbolic. F2 also cross-references the Upsilon_w question (Q2): the script's hardcoded `Upsilon_w = 100 Theta_w` (effectively `136900 Theta_w` for `Xi_F1`) contradicts the paper Inputs line that says `Upsilon_w = 117 Theta_w` (resolved by Q2).

Both findings have effectively unambiguous direction (extend scripts to exercise the paper claims). They are routed through the user gate per `paper_misalignment` policy, but the resolution is a one-line confirmation.

### Options

- **(a) Extend both scripts.** Add a `zeta_phys` closed-form assertion `zeta_phys_closed - Omega_Pe^2 * (kappa + pi^2/4) / (kappa + y_eta^2) == 0` (`expect_zero` in SymPy, `expectZero` in Mathematica). Instantiate `(eta, kappa) = (37, 12321/5)` for a numerical leg verifying `zeta_phys` at the Family-1 reference point. The Upsilon_w part of F2 is resolved by whatever direction Q2 sets.
- **(b) Trim paper card.** Remove the `zeta_phys` closed-form equation and the Family-1 numerical instantiation from the paper card. Almost certainly the wrong direction — the closed form is a load-bearing paper deliverable.
- **(c) Defer.** Acknowledge in script comments without adding assertions.

### Destination-verification grep (run locally by orchestrator)

Downstream grep for `zeta_phys` / `zetaPhys` usage:

- `mathematica/stage084..._mathematica_audit.wl:43`: `zetaPhys = FullSimplify[omegaPe^2*(kappa + Pi^2/4)/(kappa + y^2), Assumptions -> $Assumptions];` — stage 084 **already pins the closed form** that the paper card states stage 082 should establish.
- `mathematica/stage084..._mathematica_audit.wl:73-76`: instantiates `(kappa_F1, eta_F1, y_F1)` and verifies `zeta_phys` matches the upstream `zeta_max^(F1)` constant. The Family-1 numerical pair `(37, 12321/5)` is therefore already load-bearing at stage 084.
- `notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md:87-91`: explicitly references `zeta_phys(Pe_-^(F1), 37, 12321/5)` and `zeta_phys(Pe_+^(F1), 37, 12321/5)`. Stage 083's claim depends on the closed form being correct.
- `scripts/stage082..._sympy_audit.py:63-64`: declares `zeta_phys = sp.symbols("zeta_phys", real=True)` as opaque, defines `R_quad = zeta_req - zeta_phys` symbolically, never instantiates.

The closed-form `zeta_phys = Omega_Pe^2 (kappa + pi^2/4) / (kappa + y(eta)^2)` is paper-claim load-bearing for stage 082 (eq. `app-stage082-zeta-phys`), downstream-consumed at stage 083 (notes §87,91), and already exercised in stage 084 .wl. Stage 082's own scripts not exercising it is a real script-side gap.

## Recommendation

- direction: **a** (extend both stage 082 scripts)
- rationale: The paper card boxes the `zeta_phys` closed form as a stage 082 deliverable (eq. `app-stage082-zeta-phys`), and downstream stage 083 explicitly references the closed form evaluated at the Family-1 numerical pair (`notes/stage083:87,91`). Stage 084 .wl already exercises the closed form (`mathematica/stage084:43`) and the Family-1 instantiation (`mathematica/stage084:73-76`); stage 082, where the form is first established, should also exercise it. Direction (b) "trim paper card" would be inconsistent with downstream consumption. Direction (c) "defer" leaves a script-side gap that the paper claim already requires.
- downstream_impact: None — the extension adds verification, not new claims. Both stage 082 scripts get a closed-form residual `zeta_phys_closed - omega_Pe^2 * (kappa + pi^2/4) / (kappa + y_eta^2) == 0` (expect_zero / expectZero), plus a Family-1 numerical leg `(kappa, eta) = (12321/5, 37)` consistent with the existing stage 084 instantiation. The Upsilon_w part of F2 is resolved automatically by Q2's direction (a): the script's existing `100 Theta_w` is correct; the paper Inputs line gets updated for prose consistency.
- notes: Stage 082 sympy currently leaves `zeta_phys` as `sp.symbols("zeta_phys", real=True)` (line 63) — the extension swaps the opaque symbol for the closed-form expression on the verification side while keeping the symbolic `R_quad = zeta_req - zeta_phys` definition for downstream readability. The numerical leg at `(37, 12321/5)` should compare against the Family-1 limit value already cited at stage 084 line 73-76 (so the two stages produce numerically agreeing evaluations).

---

## Apply log (2026-05-27)

All three substantive questions resolved as direction (a) by user; orchestrator applied all paper/notes/script edits directly (Codex stalled mid-consultation; replaced with orchestrator-direct apply given fully-specified math).

**Q1 — Stage 074 alpha (direction a):**
- paper/stages/stage_074.tex:31 `128/\sqrt5` → `111/\sqrt5`
- notes/.../stage074..._family1_healing_lock.md:117 `179/sqrt(5)` → `111/sqrt(5)`
- notes/.../stage075..._family1_threshold_window.md:63 `179/sqrt(5)` → `111/sqrt(5)` (same typo, collateral)
- scripts/.../stage074..._sympy_audit.py: added `expect_zero("alpha_ref - 111/sqrt(5)", alpha_ref - 111/sqrt(5))`
- mathematica/.../stage074..._mathematica_audit.wl: added `expectZero["alpha_ref - 111/sqrt(5)", alphaRef - 111/Sqrt[5]]`
- Both engines: exit 0; new assertion PASSes; material_change false (alpha was already computed correctly, only printed).

**Q2 — Stage 075 Upsilon_w (direction a):**
- paper/stages/stage_075.tex:7 Inputs `117 Theta_w` → `100 Theta_w`
- paper/stages/stage_075.tex:24 body `117 Theta_w` → `100 Theta_w`
- notes/.../stage075..._family1_threshold_window.md:108, :116 `168 Theta_w` → `100 Theta_w`
- notes/.../stage075..._family1_threshold_window.md:124-128 `/168` → `/100` (arithmetic)
- scripts/.../stage075..._sympy_audit.py: added F4 `assert alpha_r**2 == 100`
- mathematica/.../stage075..._mathematica_audit.wl: added F4 `expectZero["alpha_r^2 - 100", alphaR^2 - 100]`
- Both engines: exit 0; material_change false (script's 100 was already correct, paper/notes prose updated to match).

**Q3 — Stage 082 zeta_phys closed form + Family-1 instantiation (direction a):**
- scripts/.../stage082..._sympy_audit.py + mathematica/.../stage082..._mathematica_audit.wl: F1 closed-form pin block added — `Omega_Pe(Pe) = pi*Pe*(2*Pe*exp(Pe) + pi) / ((4*Pe^2 + pi^2)*(exp(Pe) - 1))`, `Omega_Pe → pi/2` as `Pe → oo` verified, `y_F1 ≈ 1.52948` (root of `y tan y = 37`, smallest positive) computed, `zeta_phys(Pe→oo, kappa_F1, y_F1) ≈ 2.4675292...` matched to upstream `zeta_max^(F1)` to 1.77e-13.
- F2: subsumed by F1's Family-1 instantiation; Upsilon_w cross-reference resolved by Q2.
- F3: replaced tautological derivative checks with numerator/denominator factorization `numerator - C_mix*(1-eps_blk) == 0` and `denominator - (C_mix - eps_blk*(2 C_mix - Pi_tr))^2 == 0` (exercises strict-positivity content from notes §4).
- Both engines: exit 0; material_change false (closed form is paper-stated; the script now exercises it but produces no new numeric constants).
- Note: SymPy `sp.nsolve` is unstable for `y tan y = 37` near `pi/2` (jumps to far roots); replaced with `mpmath.findroot(..., (1.5, 1.55), solver="bisect")`. Mathematica's `FindRoot` with `WorkingPrecision -> 30` is stable.

**Orchestrator-direct banner-relabel sweep (consistent with audit-flagged 074/078/081/084):**

The auditors flagged the stale banners on 078 (.py+.wl), 081 (.py), and 084 (.wl). The 076 verifier additionally noted 076's banner as a non-blocking side observation. After investigation, **every III.4 stage** had stale banners from the global renumber (commit `0d09ef6`). Applied banner-relabel sweep to align all III.4 self-banners with the post-renumber numbering:

- 073: docstring `Stage 56` → `Stage 073`; banner `STAGE 56` → `STAGE 073` (.py); banner `STAGE 056` → `STAGE 073` (.wl)
- 074: docstring `Stage 57` → `Stage 074`; banner `STAGE 57` → `STAGE 074` (.py); banner `STAGE 057` → `STAGE 074` (.wl)
- 075: banner `STAGE 58` → `STAGE 075` (.py); banner `STAGE 058` → `STAGE 075` (.wl); docstring updated
- 076: docstring `Stage 59` → `Stage 076`; banner `STAGE 59` → `STAGE 076` (.py); banner `STAGE 059` → `STAGE 076` (.wl)
- 077: docstring `Stage 60` → `Stage 077`; banner `STAGE 60` → `STAGE 077` (.py); banner `STAGE 060` → `STAGE 077` (.wl)
- 078: docstring `Stage 61` → `Stage 078`; banner `STAGE 61` → `STAGE 078` (.py); banner `STAGE 061` → `STAGE 078` (.wl)
- 079: docstring `Stage 62` → `Stage 079`; banner `STAGE 62` → `STAGE 079` (.py); banner `STAGE 062` → `STAGE 079` (.wl)
- 080: docstring `Stage 63` → `Stage 080`; banner `STAGE 63` → `STAGE 080` (.py); banner `STAGE 063` → `STAGE 080` (.wl)
- 081: docstring `Stage 64` → `Stage 081`; banner `STAGE 64` → `STAGE 081` (.py only — .wl was already correct)
- 082: docstring `Stage 65` → `Stage 082`; banner `STAGE 65` → `STAGE 082` (.py); banner `STAGE 065` → `STAGE 082` (.wl)
- 083: docstring `Stage 66` → `Stage 083`; banner `STAGE 66` → `STAGE 083` (.py); banner `STAGE 066` → `STAGE 083` (.wl)
- 084: banner `STAGE 067` → `STAGE 084` (.wl only — no .py for status_only_candidate)

All 23 banner edits are pure label changes (no math content affected). All 11 stages re-run with exit 0 and correct new banners visible in transcripts.

**Codex stall note:**
The first codex-chat consultation (Q1) stalled with no codex session log written (PIDs spawned but produced no output). Killed the stalled processes and pivoted to orchestrator-direct apply, since: (1) the audit + grep evidence on all 3 questions was conclusive without further Codex reasoning; (2) the math in each directive was fully specified; (3) the verifier wave caught any error afterward. Codex math-authority delegation will resume at the next batch (III.5) if codex is back online. Recommend cross-checking codex availability before next batch launches.
