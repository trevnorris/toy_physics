---
batch_id: III.2-v2
created: 2026-05-26
auditor_version: v2 (paper-grounded)
total_user_gate_items: 2
codex_role: math-authority + cross-stage research; provide direction (a/b/c/skip) + rationale per Q
user_role: review Codex recommendations, approve or redirect before any apply session
---

# Batch III.2 v2 — user-gate item consolidation

The v2 audit of stages 049–060 (Tracking, zeta thresholds, asymmetry, boost) produced **2 paper_misalignment items** routed to user resolution. Script-side findings (covered separately via fix_loop) are NOT in this file — they go through Codex directly.

The two items are opposite directions of the same alignment axis:
- **Q1 (050 F2)** — scripts contain MORE than the paper card advertises (notes Section 5 has it, paper body does not).
- **Q2 (057 F1)** — paper card claims MORE than the scripts verify (paper says "monotone in Pe," scripts only check `partial_kappa` and `partial_y`).

For each item Codex must fill in the `## Recommendation` block with `direction:` (a/b/c/skip) plus rationale citing file:line evidence. Codex MUST NOT apply any change in this session — recommendations only.

---

## Q1 — Stage 050 F2 (paper_misalignment / paper_missing_script_claim, medium)

**Item:** Both scripts verify an enhancement-ceiling theorem `S_n^(twin)(x) < S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` (SymPy `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:90,94,104-107`; Mathematica `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:83`). The notes describe `S_n^(max)` in Section 5 (`notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md:172-174`). However, the paper card `paper/stages/stage_050.tex` body (lines 11-42) contains no boxed equation for `S_n^(max)`, and the `\stagefield{Output}` (paper:44) reads only `"Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success} and higher-harmonic exclusion/softness thresholds."` — no mention of an enhancement ceiling.

**Audit finding location:** `redteam/reports/stage_050.md` § F2; `redteam/directives/stage_050.md` § F2 (`## Resolve before fix_loop` block)

**Context:** Stage 050's paper card boxed equations are the lowest-twin sufficiency criterion `zeta_req <= 1` and the higher-harmonic exclusion thresholds `x <= x_max(n; zeta_req)`. The ceiling `S_n^(max)` is a derived quantity — it bounds how large the higher-harmonic enhancement can grow, given the bracket structure of `zeta_n^(twin)`. Whether it belongs in stage 050 at all (vs. a downstream stage that uses it) is the question.

**Options:**
- **(a) Paper card should advertise the ceiling** — extend `paper/stages/stage_050.tex` with a fifth boxed equation `S_n^(twin) < S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` and update `\stagefield{Output}` to mention "...and higher-harmonic enhancement ceiling." No script change.
- **(b) Ceiling belongs to a different stage** — strip lines `sympy:88-112` and `wl:82-95` from stage 050; relocate the ceiling check (and its notes paragraph) to whichever stage the paper card identifies as its owner. Inspect notes Section 5 + downstream stages to find the rightful home.
- **(c) Both stage 050's paper card and notes Section 5 are derived from a third source** that contradicts both — flag for deeper review (e.g., the ceiling may belong to a tower-construction stage rather than stage 050).
- **skip**

## Recommendation

direction: a
rationale: Both audits make the ceiling a Stage-050 checked claim: SymPy defines `S_n_max` and verifies saturation/factored positivity at `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:88-107`, and Mathematica does the same at `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:82-95`. The Stage-050 notes already own the result as "Exact enhancement bounds for the higher harmonics" and state `S_n^(twin) < S_n^(max)` at `notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md:163-173`. The paper card currently boxes only `S_0`, lowest-twin success, higher-harmonic exclusion, and `x_max`, then outputs only exclusion/softness thresholds at `paper/stages/stage_050.tex:13-44`; its downstream-use field names Stage 051's branch-product rewrite and Stage 052's non-twin regime, not a separate ceiling owner, at `paper/stages/stage_050.tex:45`. The requested grep found no downstream Stage 051-059 paper/notes/SymPy consumer of `S_n^(max)` or the named ceiling variants, so no downstream stage appears to be the natural relocation target.

---

## Q2 — Stage 057 F1 (paper_misalignment / script_missing_paper_claim, medium)

**Item:** Paper card `paper/stages/stage_057.tex:23` states *"It is monotone in the constructive Peclet direction, with zero-bias value"*, and notes `notes/stages/moving_throat_pde_stage057_physical_parameter_map.md:140-148,303` claim `partial_Pe zeta_0^(Pe+R) > 0` on the constructive branch `Pe >= 0`. The notes mark this as carry-forward from Stage 056 ("From Stage 39, `dOmega_Pe/dPe > 0` on the constructive branch"). However:
- SymPy `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:62-73` only computes `partial_kappa` and `partial_y`; no `sp.diff(zeta_phys, Pe)`.
- Mathematica `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:56-67` only computes `D[zetaPhys, kappa]` and `D[zetaPhys, y]`; no `D[zetaPhys, pe]`.

So the paper claims monotonicity in Pe; neither script confirms it nor explicitly references the upstream Stage-056 source.

**Audit finding location:** `redteam/reports/stage_057.md` § F1; `redteam/directives/stage_057.md` § F1 (`## Resolve before fix_loop` block)

**Cross-stage research needed:** Stage 056 ("Construction-side throat support") is the natural carry-forward source — verify that `dOmega_Pe/dPe > 0` is actually proved there. Look at:
- `paper/stages/stage_056.tex`, `notes/stages/moving_throat_pde_stage056_*.md`
- `scripts/moving_throat_pde_stage056_*_sympy_audit.py` — does it have an explicit `D[Omega_Pe, Pe]` sign check?
- Stage 057's scripts (already audited): is `Omega_Pe` imported from Stage 056 or recomputed locally?

The user-resolution direction depends on (i) whether Stage 056 actually carries the monotonicity, and (ii) whether the carry-forward is sound enough to trust without re-checking at Stage 057.

**Context:** Same character as II.1 v2's stage 035 (paper polynomial coefficient fix) and III.1 v2's Q2 (Stage-044 residual import into 045) — a downstream stage relies on an upstream property that may or may not be properly anchored in the upstream stage's scripts. The destination-verification guardrail applies.

**Options:**
- **(a) Stage-057 scripts must check Pe-monotonicity locally** — in both scripts, after the existing `partial_kappa`/`partial_y` blocks, add a numerical sweep (e.g. `Pe in {0.1, 0.5, 1, 2, 5, 10}`, `kappa = 1`, `y = pi/4`) confirming `D[zeta_phys, Pe] > 0`, plus a comment explaining the carry-forward from Stage 056. Re-run sympy+mathematica.
- **(b) Pe-monotonicity stays as Stage 056 carry-forward** — follow-up directive authorizes a paper/notes edit adding a single line ("Pe-monotonicity established at Stage 056; not re-verified here"). No stage-057 script change. The Codex apply pass does NOT include this edit; the user must issue a separate directive after choosing direction (b). Codex should first verify Stage 056 actually proves the monotonicity before recommending (b).
- **(c) Both** — user adds the carry-forward note in paper/notes AND a confirming numerical spot-check in scripts.
- **skip**

## Recommendation

direction: a
rationale: Stage 057's paper claims monotonicity in the constructive Peclet direction at `paper/stages/stage_057.tex:23`, and the notes convert the upstream claim `dOmega_Pe/dPe > 0` into `partial_Pe zeta_0^(Pe+R) > 0` at `notes/stages/moving_throat_pde_stage057_physical_parameter_map.md:140-148`. Stage 057 recomputes `Omega_Pe` locally at `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:37-40`, forms `zeta_phys` at `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:54-59`, and then audits only the `kappa` and `y` derivatives at `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:61-73`; the Mathematica audit likewise only differentiates in `kappa` and `y` at `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:56-67`. Stage 056 has a prose monotonicity argument in the notes, deriving `dOmega_Pe/dPe = Cov_Pe(chi_0,s)/I_W` and asserting positivity at `notes/stages/moving_throat_pde_stage056_transport_source_asymmetry.md:137-151`, but the Stage-056 paper output advertises only the transport interpretation and physical overlap factor at `paper/stages/stage_056.tex:39-40`. Its audits verify the derivative-covariance identity at `scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py:73-79` and `mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl:60-65`, then move on to asymptotics without an audited sign check at `scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py:81-84` and `mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl:67-70`. Because the upstream executable checks stop at an identity rather than `dOmega_Pe/dPe > 0`, Stage 057 should verify its own Pe-monotonicity claim locally.

---

## End of questions

After Codex fills in `## Recommendation` blocks above, halt and present to user for review. Apply session is a separate Codex invocation with explicit per-question scope.

## Apply log

### Q1 applied
- direction: a
- files modified: `paper/stages/stage_050.tex:43-51`
- destination_verified: yes — no duplicate found before edit; post-edit label at `paper/stages/stage_050.tex:45`
- post-edit checks: n/a — paper-only
- notes: post-edit grep found only the new label definition and the Output reference.

### Q2 applied
- direction: a
- files modified: `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:75-83`; `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:69-80`
- destination_verified: yes — no pre-existing `dPe`/`D[zetaPhys, pe]` derivative or Pe sign-check block found before edit
- post-edit checks: `python3 scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py` exit 0 with `partial_Pe zeta > 0 on constructive branch (numerical sweep): PASS`; `math -script mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl` exit 0 with `PASS: partial_Pe zeta > 0 on constructive branch (numerical sweep)`
- notes: Mathematica still emits the existing `Limit::alimv` warnings during the pre-existing closure-ceiling limit check, but exits 0.
