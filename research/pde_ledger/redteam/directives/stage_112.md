---
unit_id: 112
batch: IV.2
created_at: 2026-05-27T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-29T11:31:49-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 112

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_112.tex:1` quote: `\section[Stage~129]{Stage~129: Exact Robin–Mixed Compensation Law}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_112.tex:2` quote: `\label{stage:112}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:3` quote: `Stage 95 SymPy audit.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:54` quote: `print('stage95: PASS')`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:26` quote: `banner["STAGE 095 — EXACT ROBIN-MIXED COMPENSATION LAW"];`

## Resolve before fix_loop

The SymPy docstring/final print and the Mathematica banner refer to this unit as "Stage 95 / STAGE 095", while the paper card heading prints "Stage 129" and the LaTeX label is `stage:112`. The actual mathematics in both scripts matches the paper card; only the identifying strings disagree. Which label is canonical for script transcripts in this project?

Possible directions (the user picks one):
- (a) Internal unit number ("Stage 112") is canonical for scripts → update sympy:3 docstring to "Stage 112 SymPy audit.", sympy:54 to `print('stage112: PASS')`, and math:26 banner to `banner["STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW"];`. (math:70 already says "Stage 112 Mathematica audit passed." — leave it.)
- (b) Display section number ("Stage 129", per `\section[Stage~129]{...}`) is canonical for scripts → update the three strings above to "Stage 129" / "STAGE 129", and also update math:70 to match.
- (c) The legacy "Stage 95" identifier has provenance value that should be retained → leave script labels alone, but add a comment in each script noting the renumber chain ("formerly Stage 95; re-anchored to unit 112 / display Stage 129 in batch IV.2").

## RESOLVED — direction (a) (Claude+Codex consult 019e748e, 2026-05-29)

Per the canonical-internal-stage-number convention (IV.4/IV.5) and Codex CONCUR: **direction (a)** — the math already matches the `stage:112` card; only the ID strings were stale. This paper_misalignment now has an APPROVED direction; the "do nothing for paper_misalignment" default is OVERRIDDEN for this finding.

**Codex: apply these exact edits as part of this fix loop, then append `## Applied: F1`:**
- `.py:3` → `Stage 112 SymPy audit.`
- `.py:54` → `print('stage112: PASS')`
- `.wl:26` → `banner["STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW"];`
- `.wl:70` already says "Stage 112 Mathematica audit passed." — leave it.

No paper change.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl`
- summary: Confirmed the audit transcript labels use the canonical internal Stage 112 identifier.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:33-67`

**Issue:**
The Mathematica script follows the SymPy script's algebraic choreography line-for-line: identical series construction, identical coefficient extraction with `/I` for the L5 normalization, identical Solve target, identical chi formula `(-L5/L0)/(1/27)` with identical branch substitution, identical "scaled identity" construction at branch B with `gamma -> 1/9`. The notes file (`notes/stages/moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md`, "Branch-selection data" section) provides an explicit *independent* derivation route — the Stage-92 linearized cross-check via `(b, a_0, a_5) = (0, 3 sigma_W, -sigma_W gamma_W)` and the condition `a_0/3 + 9 a_5 = sigma_W (1 - 9 gamma_W) = 0` — which the Mathematica script does not exercise. Using that route would give genuine engine independence.

**Required change:**
Insert (alongside the existing chi_Q-based derivation, not replacing it) an independent block in the `.wl` script that derives `gamma_W = 1/9` by the Stage-92 linearized cross-check. Specifically, after line 53 (after the four branch-identity `expectZero` calls) and before the `chiA` definition at line 55, add a block such as:

```mathematica
(* Independent derivation: Stage-92 linearized (b, a0, a5) cross-check. *)
(* On the nontrivial compensated branch (solB), the net deformation data are *)
(*   b   = 0,                                                                  *)
(*   a0  = 3 sigma_W,                                                          *)
(*   a5  = -sigma_W gamma_W,                                                   *)
(* and the linearized preservation condition is a0/3 + 9 a5 == 0.              *)
bDef = (rho - sigma) /. solB;                  (* expected 0 from rho->4 sigma minus a sigma term in the trivial-cancellation seam; *)
                                                (* the constant piece of (Lambda_hyb - Lambda_out) on solB is exactly the linearized b *)
a0Def = FullSimplify[Coefficient[(lambdaHyb /. solB) - lambdaOut, z, 2] * 3, Assumptions -> $Assumptions];
                                                (* a0 = 3 * (z^2 coefficient of the deformation), matching the notes' normalization *)
a5Def = FullSimplify[Coefficient[(lambdaHyb /. solB) - lambdaOut, z, 5] / I, Assumptions -> $Assumptions];
                                                (* a5 = (z^5 imaginary coefficient of the deformation) *)
expectZero["independent: b = 0 on solB", bDef];
expectZero["independent: a0 - 3 sigma on solB", a0Def - 3*sigma];
expectZero["independent: a5 + sigma*gamma on solB", a5Def + sigma*gamma];
(* Preservation condition FACTORIZES: a0/3 + 9 a5 = sigma (1 - 9 gamma).        *)
presCond = FullSimplify[a0Def/3 + 9*a5Def, Assumptions -> $Assumptions];
expectZero["independent: preservation cond = sigma (1 - 9 gamma)", presCond - sigma*(1 - 9*gamma)];
(* NONTRIVIAL branch (sigma_W != 0) forces gamma_W = 1/9. Use Reduce WITH the    *)
(* sigma != 0 assumption so the result is unconditional (no ConditionalExpr).    *)
gammaReduce = Reduce[presCond == 0 && sigma != 0, gamma, Reals];
expectZero["independent: gamma_W = 1/9 on nontrivial branch (sigma != 0)", (gamma /. ToRules[gammaReduce]) - 1/9];
(* DEGENERATE branch (sigma_W = 0): chi_B = 1 for ANY gamma — preservation is    *)
(* trivial there, so the "iff gamma_W=1/9" claim IS the sigma != 0 statement.    *)
chiBgen = (1 - 9*sigma*gamma)/(1 - sigma);
expectZero["degenerate sigma=0: chi_B = 1 for any gamma", (chiBgen /. sigma -> 0) - 1];
```

Then leave the existing `chiA`/`chiB` block in place — the two derivations are now complementary engine paths to `gamma_W = 1/9` rather than echoes of the SymPy procedure.

Note: the auditor sketched the `(b, a0, a5)` extraction above by reading the notes' "Branch-selection data" block; Codex should verify the normalization (factor of 3 on a0, sign on a5) against the notes file directly before committing the patch. If the normalization differs, re-derive from the notes' exact wording and adjust the `expectZero` targets accordingly. The load-bearing check is the final `gammaFromLinear - 1/9 == 0`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 112` and confirm:
- the new `independent: gamma_W = 1/9 on nontrivial branch (sigma != 0)` line appears and reports `0` (PASS),
- the `preservation cond = sigma (1 - 9 gamma)` and `degenerate sigma=0: chi_B = 1 for any gamma` lines appear and PASS,
- the script exits 0,
- and the existing chi_Q-based checks still pass unchanged.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl`
- summary: Added the independent Stage-92 linearized branch-data check, including nontrivial sigma != 0 reduction and the sigma = 0 degenerate case.
- deviation: Used coefficient-shift extraction for `a0`/`a5` and exact even-ratio extraction for `b` to match the notes' `(b,a0,a5)` normalization.

## F3 — symbol_assumption_error  (sigma_W != 0 qualifier; folded in 2026-05-29)

**Source:** codex_review R1 + the 108-F1 consult side-finding — the "preservation iff gamma_W = 1/9" claim holds only on the NONTRIVIAL branch (sigma_W != 0); at sigma_W = 0, chi_B = 1 for any gamma_W.

**.wl side:** already handled by the amended F2 linearized block above (factorization `a0/3 + 9 a5 = sigma (1 - 9 gamma)` + `Reduce[... && sigma != 0]` for gamma_W = 1/9 + the degenerate-case assertion).

**.py side — required change:** the SymPy script asserts `chi_B(gamma=1/9) = 1` (the gamma=1/9 ⟹ chi_B=1 direction, true unconditionally). Add the converse qualifier so the "iff" is honest. Using the branch-B normalization `chi_B = (1 - 9*sigma_W*gamma_W)/(1 - sigma_W)` symbolic in (sigma_W, gamma_W), add immediately after the existing gamma=1/9 assertion:
- a comment: preservation (chi_B = 1) holds iff `sigma_W*(1 - 9*gamma_W) = 0`, i.e. on the nontrivial branch (sigma_W != 0) iff gamma_W = 1/9; at sigma_W = 0, chi_B = 1 for any gamma_W;
- `expect_zero('chi_B - 1 factorizes as sigma(1-9 gamma)', sp.together(chi_B_general - 1).as_numer_denom()[0] - (sign-matched) sigma_W*(1 - 9*gamma_W))` — match SymPy's numerator sign normalization (sympy may return `-sigma_W*(9*gamma_W - 1)`; adjust to land at 0);
- `expect_zero('degenerate sigma=0: chi_B = 1 for any gamma', chi_B_general.subs(sigma_W, 0) - 1)`.

Use the script's actual symbol names; define a symbolic `chi_B_general = (1 - 9*sigma_W*gamma_W)/(1 - sigma_W)` if the script's `chi_B` already has gamma=1/9 substituted. Load-bearing: the factorization (chi_B = 1 ⟺ sigma_W(1 - 9 gamma_W) = 0) and the degenerate-case assertion.

**Verification command:** `redteam exec-sympy 112` shows the new factorization + degenerate assertions PASS, exit 0, existing checks unchanged.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py`
- summary: Added the branch-B factorization and sigma = 0 degenerate-case checks after the gamma = 1/9 assertion.
- deviation: none
