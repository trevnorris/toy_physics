---
unit_id: 037
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md"]
  paper_appendix: present
---

# Audit unit 037 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_037.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 52 + summary rows 192/307/345/359)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage037_continuum_kernel_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt`

## What the paper claims

Stage 037 replaces the abstract Part-II selected-mode branch data by closed formulas extracted from an explicit finite-throat continuum kernel. The `\stagefield{Output}` reads: *"Closed continuum formulas \eqref{eq:app-stage037-A-DK-continuum}--\eqref{eq:app-stage037-M-delta-continuum} and the stability gate \eqref{eq:app-stage037-stability}."* The distinct deliverables are the boxed closed forms for the baseline scalar stiffness `A = (K_U K_eta^eff - c_etaU^2)/(mu_eta K_U)`, the axial gap `DeltaK_ax = pi^2 T_w/(mu_eta L^2)`, the mixed loading `alpha_mix = (K_U c_etaW + c_UW c_etaU)^2 / [mu_eta K_U (K_U K_W^eff - c_UW^2 sigma)]`, the outgoing normalization `beta_0 = (mu_W/mu_eta)(K_U c_etaW + c_UW c_etaU)^2/(K_U K_W^eff - c_UW^2 sigma)^2`, the dimensionless mixed baseline `M_mix = 8(...)^2/[pi^2 (K_U K_eta^eff - c_etaU^2)(K_U K_W^eff - c_UW^2 sigma)]`, and the anisotropy ratio `delta = pi^2 T_w K_U/[L^2(K_U K_eta^eff - c_etaU^2)]`, plus the two-part stability gate `K_U K_eta^eff > c_etaU^2` and `K_U K_W^eff > c_UW^2 sigma`. The card's `\stagefield{Inputs}` pins the overlap data `kappa_0 = 2 sqrt2/pi`, `kappa_1 = -4/(3 pi)`, `sigma = 88/(9 pi^2)`. The notes file additionally documents the intermediate normalized-coordinate route (`K_0`, `Omega_U^2`, `Omega_W^2 = K_W^eff/mu_W`, `Delta_0`, `Chi`) and the exact rank-1 Schur factorization `Sigma_wall = Xi I_2 + alpha v v^T`, which the card states only structurally in \eqref{eq:app-stage037-schur-form}.

## What the script claims to verify

The SymPy script (a) constructs the exact N/N modes `u_0,u_1` and the D/N half-wave `f_0`, verifies their normalization, orthogonality, and boundary conditions, and computes the overlaps to confirm `kappa_0 = 2 sqrt2/pi`, `kappa_1 = -4/(3 pi)`, `sigma = 88/(9 pi^2)`; (b) builds the 4x4 internal block `B` and the 2x4 coupling `C`, forms the exact Schur complement `Sigma = C B^{-1} C^T`, and asserts it equals `Xi I + alpha v v^T` with the closed-form `Xi`, `alpha`; (c) substitutes the continuum couplings `g_* = c_*/sqrt(...)` into `A`, `Delta0`, `Chi`, `beta0`, `alpha_mix`, `M_mix`, `delta`, and asserts each equals the boxed closed form; (d) re-asserts the `A` and `Delta0` numerators for the stability gates. The Mathematica script is an independent re-derivation that adds two stronger numerator-cross-multiplication checks (`A numerator matches Schur form`, `delta numerator matches closed form`) and an extra rank-1 self-consistency check (it solves `alpha`,`Xi` from the (1,1)/(1,2) entries and verifies the (2,2) entry). Every assertion is an `expect_zero`/`expectZero` against a derived (not hardcoded) target.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `A = (K_U K_eta^eff - c_etaU^2)/(mu_eta K_U)` | py L186/208–211; wl L176/208–211 | match |
| `DeltaK_ax = pi^2 T_w/(mu_eta L^2)` | py L102 (derived from premise, printed); wl L84; used in `delta` check | match (derived, asserted via `delta`) |
| `alpha_mix = (...)^2/[mu_eta K_U(K_U K_W^eff - c_UW^2 sigma)]` | py L190; wl L187 | match |
| `beta_0 = (mu_W/mu_eta)(...)^2/(K_U K_W^eff - c_UW^2 sigma)^2` | py L189; wl L186 | match |
| `M_mix = 8(...)^2/[pi^2(K_U K_eta^eff - c_etaU^2)(K_U K_W^eff - c_UW^2 sigma)]` | py L191; wl L188 | match |
| `delta = pi^2 T_w K_U/[L^2(K_U K_eta^eff - c_etaU^2)]` | py L192; wl L189–196 | match |
| Stability `K_U K_eta^eff > c_etaU^2`, `K_U K_W^eff > c_UW^2 sigma` | py L208–215 (re-asserts the two numerators that gate positivity); wl L208–215 | match (gate identities verified; inequality itself is a printed corollary of the identity sign) |
| Schur form `Sigma_wall = Xi I_2 + alpha v v^T` (card eq, notes §4) | py L144; wl L123,134 | match |
| Overlap inputs `kappa_0,kappa_1,sigma` | py L83–85; wl L68–70 | match |

All paper deliverables map to a non-tautological script-side check that derives the LHS from the continuum premises and subtracts the boxed closed form. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65–73 | `expect_zero(int/boundary)` | mode-basis premise for overlaps | yes |
| A2 | sympy | 83–85 | `kappa0/kappa1/sigma - target == 0` | Inputs `kappa,sigma` | yes |
| A3 | sympy | 144 | `Sigma - (Xi I+alpha vv^T) == 0` | Schur form eq. (app-stage037-schur-form) | yes |
| A4 | sympy | 186 | `A - A_expected == 0` | `A` closed form | yes |
| A5 | sympy | 187–192 | `Delta0/Chi/beta0/alpha_mix/M_mix/delta - expected == 0` | the five+1 closed forms | yes |
| A6 | sympy | 208–215 | `A`/`Delta0` numerator identities | stability gate | yes |
| B1 | math | 51–70 | `expectZero(modes,overlaps)` | mode basis + Inputs | yes |
| B2 | math | 123 | `expectMatrixZero(Sigma-...)` | Schur form | yes |
| B3 | math | 134 | `(2,2) entry consistency` | Schur rank-1 form (independent fit) | yes |
| B4 | math | 176 | `a - aDerived == 0` | `A` reduces to closed form | yes |
| B5 | math | 177–183 | `A numerator * den - ... == 0` (cross-mult) | `A` closed form (strong) | yes |
| B6 | math | 184–188 | `Delta0/Chi/beta0/alphaMix/mMix - expected == 0` | the closed forms | yes |
| B7 | math | 189–196 | `delta numerator cross-mult == 0` | `delta` closed form (strong) | yes |
| B8 | math | 208–215 | `A`/`Delta0` gate numerators | stability gate | yes |

No assertion is tautological-by-construction: in every case the LHS is assembled from the mode integrals + Schur inverse + continuum coupling substitutions, and the RHS is the independently typed boxed formula, so a transcription error in either route would make the residual nonzero. The Schur-complement check (A3/B2) is the genuine load-bearing one — it confirms the abstract rank-1 loading law of Part II is the exact first Schur complement, which is the stage's actual physics content.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage037_continuum_kernel_sympy_audit.txt` (mtime 2026-05-26 02:02) vs script mtime 2026-06-03 15:59
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt` (mtime 2026-05-26 02:03) vs script mtime 2026-06-03 15:59

**What's wrong:**
Both committed `.txt` transcripts predate their scripts by ~8 days, and the content confirms the staleness rather than being a mere timestamp artifact. The SymPy script's current banner prints `STAGE 37 — CONTINUUM-KERNEL EXTRACTION AUDIT` (py L51) and footer `STAGE 37 AUDIT COMPLETE` (L222), but the committed SymPy output still reads `STAGE 20 — ...` (txt L3) and `STAGE 20 AUDIT COMPLETE` (txt L71). The Mathematica script's current banner prints `STAGE 037 — CONTINUUM-KERNEL EXTRACTION` (wl L40), but the committed Mathematica output header reads `STAGE 020 — ...` (txt L3) while its own footer line `Stage 037 Mathematica audit passed.` (txt L95) is current — i.e. the transcript was produced from a mixed/older state. The numeric/symbolic result lines themselves agree with the current scripts (all residuals `= 0`, all closed forms identical), so the staleness is in the banner labels only and does not change any verified value.

**Why this matters:**
The committed transcript is the audit's evidentiary record. A reader cross-checking the saved output against the current script sees a banner-number mismatch (20 vs 37/037) that looks like the wrong stage was run. Refreshing the outputs removes the ambiguity.

**Required change:**
Re-run both scripts and overwrite the two `.txt` files with the current output. No script-source change is required for this finding (the result lines already match); this is the orchestrator's standard fresh-run step. Informational, non-blocking.

**Verification:**
After a fresh run, `scripts/output/..._sympy_audit.txt` line 3 reads `STAGE 37 — ...` (matching py L51) and line 71 reads `STAGE 37 AUDIT COMPLETE`; `mathematica/output/..._mathematica_audit.txt` line 3 reads `STAGE 037 — ...` (matching wl L40).

### F2 — paper_misalignment (stale self-label, deferred class)

**Subtype:** notes_contradicts_script (label-level only; no math impact)
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py:3` — `moving_throat_pde_stage20_continuum_kernel_sympy_audit.py`
- `.../sympy_audit.py:5` — `Stage 20 SymPy audit:`
- `.../sympy_audit.py:6` — `derive the Stage-17/19 reduced branch data ...`
- `.../sympy_audit.py:224` — `print("Stage-17/19 branch data exactly.")`
- (banner) `.../sympy_audit.py:51,222` — `STAGE 37` (uses `37`, not the canonical zero-padded `037`)

**What's wrong:**
The SymPy docstring still self-identifies as `stage20` / `Stage 20` and refers to its forward-data source as `Stage-17/19`, the pre-renumber labels from before the EM-extension realignment (canonical is stage 037; per [[numbering-drift-root-cause]] the +17 chain maps 20→037, 17→034, 19→036). The Mathematica `.wl` is already clean (banner `STAGE 037`, footer `Stage 037`), so this is an internal inconsistency between the two engines' self-labels as well. None of these labels are load-bearing for any assertion — every computed quantity and every boxed target is correct, the engines agree, and the cross-check against the paper is exact.

**Why this matters:**
Per the user-approved plan (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, commit 1fa4f7a) the SCRIPT/OUTPUT-band stale self-labels (bare `Print`/docstring/banner stage numbers, including hyphenated cross-refs like `Stage-17/19`) are a known, pervasive class deferred to a DEDICATED content-keyed pass, NOT to be offset-swept. This finding is recorded so the dedicated pass picks up stage 037's SymPy docstring/footer; it is explicitly **not** a blocking script-math defect and **not** for ad-hoc Codex repair in this fix_loop.

**Required change:**
None in this fix_loop. Routed to the deferred numbering pass: in the dedicated content-keyed cleanup, retag py L3 `stage20`→`stage037`, L5 `Stage 20`→`Stage 037`, L6/L224 `Stage-17/19`→`Stage 034/036` (the canonical sources of the reduced branch data), and normalize the banners L51/L222 `STAGE 37`→`STAGE 037` to match the `.wl`. Content-keyed, never offset-swept.

**Verification:**
Handled by the dedicated numbering pass; no assertion changes, scripts still exit 0 unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent** re-derivation, not a transliteration. Evidence: (1) it uses `LinearSolve[bMat, Transpose[cMat]]` for the Schur complement (wl L114) where the `.py` uses an explicit `B.inv()` matrix product (py L134) — a different linear-algebra route; (2) it adds checks the SymPy script does not have — an *independent rank-1 recovery* that `Solve`s `alpha` and `Xi` from the (1,1)/(1,2) entries and then tests the (2,2) entry (wl L124–134), and two numerator cross-multiplication identities (`A numerator matches Schur form` wl L177–183, `delta numerator matches closed form` wl L189–196) that confirm the closed form without dividing, guarding against a `Together`/`simplify` masking a denominator error; (3) symbol-domain handling differs — the `.wl` carries explicit `aU != 0 && ... && (aU*aW - gR^2 sigma) != 0` nonzero assumptions for the inverse (wl L99–100) whereas SymPy relies on the generic symbolic inverse. These are genuine second-engine checks, not echoed algebra.

## Engine cross-check

Final closed forms agree across engines. `A`: py `(K_U*(K_eta + 6*T_Omega) - c_etaU**2)/(K_U*mu_eta)` (sympy txt L52) = wl `(-cEtaU^2 + kU*(kEta + 6*tOmega))/(kU*muEta)` (math txt L75) — identical. `delta`: py `pi**2*K_U*T_w/(L**2*(K_U*(K_eta + 6*T_Omega) - c_etaU**2))` = wl `(kU*Pi^2*tw)/(ell^2*(-cEtaU^2 + kU*(kEta + 6*tOmega)))` — identical. `M_mix`: py `288*L**2*(...)^2/((K_U*(K_eta+6*T_Omega)-c_etaU**2)*(9*pi**2*K_U*(4*K_W*L**2+pi**2*T_W)-352*L**2*c_UW**2))` (sympy txt L57) vs wl `(-288*ell^2*(...)^2)/((cEtaU^2 - kU*(kEta+6*tOmega))*(4*ell^2*(-88*cUW^2+9*kU*kW*Pi^2)+9*kU*Pi^4*tW))` (math txt L80); the leading `-` on the wl form cancels against `(cEtaU^2 - kU*kEtaEff) = -(kU*kEtaEff - cEtaU^2)`, and the denominator factor `4*ell^2*(-88*cUW^2+9*kU*kW*Pi^2)+9*kU*Pi^4*tW = 9*pi^2*K_U*(4*K_W*L^2+pi^2*T_W) - 352*L^2*c_UW^2` (since `352 = 4*88` and `9*Pi^2*K_U*4*K_W*L^2 + 9*Pi^2*K_U*Pi^2*T_W = 36*Pi^2 K_U K_W L^2 + 9 Pi^4 K_U T_W`, and `9 Pi^2 K_U Pi^2 T_W = 9 Pi^4 K_U T_W` matches the `+9*kU*Pi^4*tW` term) — same rational function. All 7 deliverables, both stability numerators, the Schur matrix, and the overlaps match. No `engine_disagreement`.

## Verdict justification

The script math is sound and aligned with the paper. Attacks tried and failed: (1) tautology hunt — every assertion subtracts an independently-typed boxed target from a quantity assembled through the mode-integral/Schur/coupling route, so none is guaranteed-zero by construction; the Schur (1,1)/(1,2)→(2,2) recovery in the `.wl` is a real falsifiable consistency check. (2) Symbol-domain attack — the only inverse/`B.inv()` is over symbols where the `.wl` explicitly assumes the determinant nonzero and the `.py` uses the generic symbolic inverse; positivity assumptions (`mu_*>0`, `K_*>0`, `L>0`) match both the physical setup and the stated stability gates, and do not pre-bake the inequalities (the gates are verified as numerator *identities*, with the inequality printed as a sign corollary). (3) Engine-disagreement attack — all 7 closed forms and both gate numerators match between SymPy and Mathematica after accounting for cosmetic sign/grouping. (4) Paper-misalignment attack on values — every emitted deliverable reconciles with the boxed `.tex` formulas and the notes (see Value Reconciliation). The only findings are F1 (`stale_output`, informational — committed transcripts predate the banner update and show `STAGE 20/020` while result lines are current) and F2 (deferred SCRIPT/OUTPUT-band stale self-labels in the SymPy docstring/banner, explicitly routed to the dedicated numbering pass, not this fix_loop). Neither touches a verified value, so the substantive verdict is effectively clean; `verdict: findings` only because two label/freshness items are recorded. I confirm I read the paper card, the notes, and the appendix row, and the script's verified claim matches the paper's stated claim exactly.

## Value Reconciliation (pass-2 augmentation)

Every deliverable value the scripts emit, located in the `.tex` card and/or the `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `kappa_0 = 2 sqrt2/pi` | py L83 / sympy txt L18; wl L68 / math txt L27 | tex L7 (`\kappa_0=2\sqrt2/\pi`); md L101 | MATCH |
| `kappa_1 = -4/(3 pi)` | py L84 / sympy txt L19; wl L69 / math txt L27 | tex L7 (`\kappa_1=-4/(3\pi)`); md L103 | MATCH |
| `sigma = 88/(9 pi^2)` | py L85 / sympy txt L20; wl L70 / math txt L29 | tex L7 (`\sigma=\dots=88/(9\pi^2)`); md L107 | MATCH |
| `A = (K_U K_eta^eff - c_etaU^2)/(mu_eta K_U)` | py L194 / sympy txt L52; wl L198 / math txt L75 | tex L65 (eq app-stage037-A-DK-continuum); md L233 | MATCH |
| `DeltaK_ax = pi^2 T_w/(mu_eta L^2)` | py L102 / sympy txt L29; wl L84 / math txt L41 | tex L67; md L139/L235 | MATCH |
| `Delta0 = (K_U K_W^eff - c_UW^2 sigma)/(mu_U mu_W)` | py L195 / sympy txt L53; wl L199 / math txt L76 | md L237 (Delta0 is the internal Schur quantity; card states it structurally) | MATCH (notes) |
| `Chi = (K_U c_etaW + c_UW c_etaU)/(mu_U sqrt(mu_eta mu_W))` | py L196 / sympy txt L54; wl L200 / math txt L77 | md L239–240 | MATCH (notes) |
| `beta_0 = (mu_W/mu_eta)(K_U c_etaW + c_UW c_etaU)^2/(K_U K_W^eff - c_UW^2 sigma)^2` | py L197 / sympy txt L55; wl L201 / math txt L78 | tex L79–82 (eq app-stage037-beta-continuum); md L251–253 | MATCH |
| `alpha_mix = (K_U c_etaW + c_UW c_etaU)^2/[mu_eta K_U(K_U K_W^eff - c_UW^2 sigma)]` | py L198 / sympy txt L56; wl L202 / math txt L79 | tex L72–74 (eq app-stage037-alpha-beta-continuum); md L244–246 | MATCH |
| `M_mix = 8(K_U c_etaW + c_UW c_etaU)^2/[pi^2(K_U K_eta^eff - c_etaU^2)(K_U K_W^eff - c_UW^2 sigma)]` | py L199 / sympy txt L57; wl L203 / math txt L80 | tex L86–89 (eq app-stage037-M-delta-continuum); md L257–260 | MATCH |
| `delta = pi^2 T_w K_U/[L^2(K_U K_eta^eff - c_etaU^2)]` | py L200 / sympy txt L58; wl L204 / math txt L81 | tex L90–91; md L264–266 | MATCH |
| Stability gate `K_U K_eta^eff > c_etaU^2`; `K_U K_W^eff > c_UW^2 sigma` | py L218–219 / sympy txt L66–67; wl L218–219 / math txt L66–67 | tex L96–99 (eq app-stage037-stability); md L281/L290 | MATCH |

Internal scaffolding (accounted for, no finding expected): `K_0 = K_eta^eff/mu_eta` (intermediate, in notes L137); `varpi^2` (the support/BdG field frequency, printed but not a Stage-037 deliverable — it carries the `phi`-sector kernel, in notes L141); `Omega_U^2 = K_U/mu_U`, `Omega_W^2 = K_W^eff/mu_W` (intermediate, notes L143–146); `Xi (recovered) = gU^2/aU`, `alpha (recovered)` (frequency-space Schur entries printed by the `.wl` recovery step, intermediate); `K_eta^eff = K_eta + 6 T_Omega`, `K_W^eff = K_W + pi^2 T_W/(4 L^2)` (named combinations, both in notes L227–229 and tex L20/L61). All are present in the notes and are intermediate, so even the printed ones are MATCH-or-INTERNAL, not MISSING.

`reconciliation: complete; 12 deliverable values checked, 0 misaligned`

Note: the `varpi^2` printed by both scripts is the only emitted symbolic quantity that is neither a Stage-037 boxed deliverable nor a stability gate; it is the `phi`-sector frequency carried for completeness of the continuum kernel and is documented in the notes (L141), so it is INTERNAL, not a MISSING-DELIVERABLE.

## Self-test notes

Variable-independence trap: no `sp.diff`/`D[...]` derivative-of-an-independent-symbol checks are present that could vanish identically — the only derivatives (`D[u0,s]`, `D[u1,s]`, `D[f0,s]`, py L69–73 / wl L55–59) are mode-slope boundary conditions where the function genuinely depends on `s`, so the BC checks are non-trivial. Symmetry/parity: the overlap integrals `int u_i f_0` are over the bounded `[0,L]`, not an unbounded symmetric domain, so no odd/even-vanishing trap applies; the nonzero `kappa_0,kappa_1` results are consistent with the products `u_i f_0` not integrating to zero. Trivial-case: the load-bearing `expect_zero` residuals are differences of two independently-built rational functions; a sign or factor error in either the Schur route or the boxed target would leave a nonzero residual (confirmed by the `M_mix`/`delta` numerator cross-multiplication checks in the `.wl`, which avoid denominator-masking). No directive with Codex-applied edits is written: F1 is an orchestrator fresh-run, F2 is a deferred-pass numbering label, and neither is a fix_loop script-math change.
