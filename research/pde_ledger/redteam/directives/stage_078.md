---
unit_id: 078
batch: III.4
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 078 (v2)

Apply each non-`paper_misalignment` finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-`paper_misalignment` finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:45-47` (comment block)

**Issue:**
The v1 directive (F3) required the Mathematica script to derive its coefficients from symbolic closed forms. The applied patch correctly derives `thetaFailSym` from the `Sinh/Cosh` closed form, but bootstraps `thetaSuffSym` from `thetaFailSym` times a literal decimal ratio rather than its own Stage-75 closed form. This leaves part of the engine independence unfulfilled. Additionally, the inline comment at lines 45-47 claims the script "verify[s] their ratio chi:J matches the Stage-77 ratio", but no such ratio check exists.

**Required change:**

(a) In `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`, replace line 44

```mathematica
thetaSuffSym = thetaFailSym * (4.21495341569977*^-2 / 3.62605617972939*^-4);
```

with the explicit symbolic closed form from Stage-75 (sympy output line 21):

```mathematica
(* Independent closed form for Theta_suff from Stage-75 sympy output line 21:
   Theta_suff/Pe_req = -(45 cosh(alpha) + 27 sqrt(5) sinh(alpha))
                        / (2500 - 2500 cosh(alpha)),  alpha = 111 sqrt(5)/5. *)
thetaSuffSym = -(45 Cosh[111 Sqrt[5]/5] + 27 Sqrt[5] Sinh[111 Sqrt[5]/5])
               / (2500 - 2500 Cosh[111 Sqrt[5]/5]);
```

(b) In the same file, at lines 45-47, the existing comment reads:

```mathematica
(* The chi^2 and Jensen-floor Theta values are recorded numerically in
   stage077 output; we adopt them at high precision but verify their
   ratio chi:J matches the Stage-77 ratio.                              *)
```

Remove the misleading "verify their ratio" claim by replacing those three lines with:

```mathematica
(* The chi^2 and Jensen-floor Theta values are adopted from Stage-77 at
   extended precision; their independent re-derivation is the subject of
   the Stage-077 audit, not this one.                                    *)
```

(No ratio check is being added; the goal is to remove a comment that overpromises what the script does.)

**Verification:**
After Codex applies:
- Line 44 contains `Cosh[111 Sqrt[5]/5]` and `Sinh[111 Sqrt[5]/5]`, no literal decimals on its RHS.
- Lines 45-47 no longer claim a ratio check.
- `redteam exec-mathematica 078` exits 0, all `expectApprox` and `expectTrue` PASSes remain.
- The numeric value of `thetaSuffCoeff = N[thetaSuffSym, 30]` is unchanged at the working precision (still `0.042149534156997728721`).

**Self-test mental walk-through:**
The Stage-75 closed form for `Theta_suff/Pe_req` is `-(45 cosh(α) + 27√5 sinh(α)) / (2500 - 2500 cosh(α))` with α = 111√5/5 ≈ 49.62. At α ≈ 49.62, cosh ≈ sinh, both very large positive. Numerator → -(45 + 27√5) × cosh(α) ≈ -105.37 × cosh(α); denominator → -2500 × cosh(α). Ratio → 105.37/2500 ≈ 0.04215, matching the Stage-75 numeric 0.042149534156997728721 — confirmed.

The `thetaSuff = thetaSuffCoeff*peReq` line (currently line 58) still consumes `thetaSuffCoeff = N[thetaSuffSym, 30]` at line 53, which is unaffected by this change in form. The downstream `peSuffChi = thetaChiCoeff/thetaSuffCoeff` etc. arithmetic is unchanged in value. All four `expectApprox` PASSes hold.

## F2 — paper_misalignment

**Subtype:** notes_contradicts_script (and partially target_mismatch on the script banners)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_078.tex:1` quote: `\section[Stage 078]{Stage 078: Explicit Family--1 Branch Comparison and Closing Verdict for This Subprogram}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md:1` quote: `# Moving-Throat PDE — Stage 078: ...` (filename and H1 use "078"; body uses legacy "Stage 58 / 60 / 61")

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:3` quote: `"""Stage 61 SymPy audit."""`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:23` quote: `banner("STAGE 61 — FAMILY-1 BRANCH VERDICT")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32` quote: `banner["STAGE 061 — FAMILY-1 BRANCH VERDICT"];`

## Approved by user (2026-05-27)

- direction: (b) — script-side banner relabel only; notes/ legacy numbering preserved as historical record.
- applied_by: orchestrator on 2026-05-27 (script-side):
  - `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:3` docstring `Stage 61` -> `Stage 078`
  - `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:23` banner `STAGE 61` -> `STAGE 078`
  - `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32` banner `STAGE 061` -> `STAGE 078`
- Codex applies F1 below (mathematica_transliteration); F2 is closed.
