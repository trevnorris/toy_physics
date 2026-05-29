---
unit_id: 045
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-26T07:44:00Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 045

Apply each non-paper_misalignment finding below in order. After applying each, append an `## Applied: F<n>` block under that finding with `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For F3, the orchestrator is holding for user resolution; do NOT apply F3 until the user has identified the Stage-27 normalization residual expression to import. The `## Resolve before fix_loop` block below states the question.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:34-69` and `:78-91` and `:125-136`

**Issue:** The Mathematica script is a line-by-line port of the SymPy script. Variable names, channel-list construction, branch-equation form, and `fTrack` hand-written expression all mirror the SymPy script step for step. This violates the second-engine independence policy: the wl script must derive at least one of the load-bearing identities by a route the py script does not use, so that a hidden error in the py algebra cannot survive into the wl path unchanged.

**Required change:**

Replace the wl algebra in three places with independent routes while preserving the same set of `expectZero` assertions and the same final result. F1 depends on F3's resolution because the primary independence fix is at the F_tr check; treat F1 as deferred until F3 is resolved. If F3 is applied as described below, F1 is satisfied at the F_tr check; the remaining wl independence work is to vary the routes for the coherent identity and the branch equation:

- wl lines 36–47 (coherent identity): replace the manual coefficient extraction with a direct check using the reference forms. After F2's cleanup, the substantive check becomes:

  Before:
  ```
  cWeta = Coefficient[Coefficient[couplingDensity, Wsym], etasym];
  cWU = Coefficient[Coefficient[couplingDensity, Wsym], Usym];
  cPhiEta = Coefficient[Coefficient[couplingDensity, phisym], etasym];
  cPhiU = Coefficient[Coefficient[couplingDensity, phisym], Usym];
  gWext = cWeta/Sqrt[muEta*muW];
  gRext = -cWU/Sqrt[muU*muW];
  gBext = cPhiEta/Sqrt[muEta*muPhi];
  gSext = -cPhiU/Sqrt[muU*muPhi];
  gW = lamW/Sqrt[muEta*muW];
  gR = gamma*lamW/Sqrt[muU*muW];
  gB = lamPhi/Sqrt[muEta*muPhi];
  gS = gamma*lamPhi/Sqrt[muU*muPhi];
  expectZero["g_W extracted - reference", gWext - gW];
  expectZero["g_R extracted - reference", gRext - gR];
  expectZero["g_B extracted - reference", gBext - gB];
  expectZero["g_S extracted - reference", gSext - gS];
  ```

  After:
  ```
  (* Independent route: derive g_* from the coupling density via partial derivatives,
     not via the SymPy-style `.coeff(...).coeff(...)` chain. The bilinear coefficients
     of W*eta, W*U, phi*eta, phi*U appear as second cross-derivatives. *)
  gWext = D[D[couplingDensity, Wsym], etasym]/Sqrt[muEta*muW];
  gRext = -D[D[couplingDensity, Wsym], Usym]/Sqrt[muU*muW];
  gBext = D[D[couplingDensity, phisym], etasym]/Sqrt[muEta*muPhi];
  gSext = -D[D[couplingDensity, phisym], Usym]/Sqrt[muU*muPhi];
  gW = lamW/Sqrt[muEta*muW];
  gR = gamma*lamW/Sqrt[muU*muW];
  gB = lamPhi/Sqrt[muEta*muPhi];
  gS = gamma*lamPhi/Sqrt[muU*muPhi];
  expectZero["g_W extracted - reference", gWext - gW];
  expectZero["g_R extracted - reference", gRext - gR];
  expectZero["g_B extracted - reference", gBext - gB];
  expectZero["g_S extracted - reference", gSext - gS];
  ```

  The `D[D[expr, var1], var2]` route is genuinely different from the SymPy `coeff(...).coeff(...)` chain; the two would diverge if the coupling density were nonlinear in W, phi, eta, or U.

- wl lines 98–112 (branch equation collapse): the SymPy script does `branch_eq.subs(R_phi, R_U)`, while the current wl uses `Series[branchNumRaw, {rPhi, rU, 0}]` (already somewhat divergent — leave this part alone, it is the existing independence).

**Verification:** After Codex applies, the verifier will run `redteam exec-mathematica 045` (after F3 is also applied) and confirm:
- the new `D[D[...]]` route still produces the same `g_*` values (the four `expectZero` assertions still pass),
- the wl transcript no longer textually mirrors the py extraction step,
- the wl script still exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- summary: Replaced the WL coefficient-chain extraction with second partial derivatives of the coupling density.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:114-129`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:78-91`

**Issue:** `M_tr - channel_sum` is built from the same template on both sides; the assertion cannot fail and provides no verification of the total-baseline formula. The substantive total-baseline claim is downstream of Stages 022 and 026, not Stage 045.

**Required change:**

In the SymPy script, edit lines 114–129 from:
```python
channels = [
    ("W", Z_W, eps_W_split),
    ("phi", Z_phi, eps_phi_split),
]
prefactor = 8 * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta))
M_tr_channel_sum = sp.simplify(
    sum(prefactor * Z_i / (1 - eps_i) for (_, Z_i, eps_i) in channels)
)
M_mix = sp.simplify(prefactor * Z_W / (1 - eps_W_split))
M_supp = sp.simplify(prefactor * Z_phi / (1 - eps_phi_split))
M_tr = sp.simplify(M_mix + M_supp)
print("M_mix  =", M_mix)
print("M_supp =", M_supp)
print("M_tr   =", M_tr)
print("M_tr_channel_sum =", M_tr_channel_sum)
expect_zero("M_tr - channel_sum", M_tr - M_tr_channel_sum)
```
to:
```python
prefactor = 8 * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta))
M_mix = sp.simplify(prefactor * Z_W / (1 - eps_W_split))
M_supp = sp.simplify(prefactor * Z_phi / (1 - eps_phi_split))
M_tr = sp.simplify(M_mix + M_supp)
print("M_mix  =", M_mix)
print("M_supp =", M_supp)
print("M_tr   =", M_tr)
# M_mix and M_supp are carried forward from Stages 022 and 026 in symbolic form;
# the substantive verification of the prefactor structure lives in those stages.
```

In the Mathematica script, edit lines 78–91 from:
```
channels = {{ZW, epsWSplit}, {ZPhi, epsPhiSplit}};
prefactor = 8*(1 + chi0)^2/(Pi^2*(1 - epsEta));
mTrChannelSum = FullSimplify[
  Total[prefactor*#[[1]]/(1 - #[[2]]) & /@ channels],
  Assumptions -> $Assumptions
];
mMix = FullSimplify[prefactor*ZW/(1 - epsWSplit), Assumptions -> $Assumptions];
mSupp = FullSimplify[prefactor*ZPhi/(1 - epsPhiSplit), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["M_tr = ", fmt[mTr]];
Print["M_tr_channel_sum = ", fmt[mTrChannelSum]];
expectZero["M_tr - channel_sum", mTr - mTrChannelSum];
```
to:
```
prefactor = 8*(1 + chi0)^2/(Pi^2*(1 - epsEta));
mMix = FullSimplify[prefactor*ZW/(1 - epsWSplit), Assumptions -> $Assumptions];
mSupp = FullSimplify[prefactor*ZPhi/(1 - epsPhiSplit), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["M_tr = ", fmt[mTr]];
(* M_mix and M_supp are carried forward from Stages 022 and 026 in symbolic form;
   the substantive verification of the prefactor structure lives in those stages. *)
```

**Verification:** After Codex applies, the verifier will run `redteam exec-sympy 045` and `redteam exec-mathematica 045` and confirm:
- the transcripts no longer contain the tautological `M_tr - channel_sum = 0` line,
- the `M_mix`, `M_supp`, `M_tr` print lines remain,
- both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- summary: Removed the tautological channel-sum total-baseline assertions while preserving the M_mix, M_supp, and M_tr printouts.
- deviation: none

## F3 — paper_misalignment

**Subtype:** notes_contradicts_script (the script's F_tr check fails to exercise the notes-stated deliverable)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_045.tex:38-43` quote (eq:app-stage045-tracking-laws): "M_tr = G_tr(xi, delta; R_tr), R_target = F_tr(xi, delta; R_tr)"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:210-226` quote: "Because the coherent local kernel lands on the tracking surface, the Stage-27 normalization residual also collapses to the Stage-23 tracking form: `R_target = F_tr(xi, delta; R_tr)`, with F_tr(xi, delta; R) = [9 delta + (9+2 R^2) xi]^2 [9 delta + (9+2 R) xi]^2 / [81 (1-xi) (9 delta^2 + 18 delta xi + (9+2 R^2) xi^2)^2]."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:178-189` quote: "F_track = sp.simplify((delta + (1 + lam0*R_U**2)*xi)**2 * (delta + (1 + lam0*R_U)*xi)**2 / ((1 - xi) * ((delta + xi)**2 + lam0*R_U**2 * xi**2)**2)) ... expect_zero('F_tr normalization law', sp.simplify(F_track.subs(lam0, lam0_dn)) - F_tr_expected)"
- The script never imports or restates the Stage-27 normalization residual; the check only verifies that an algebraically equivalent generic-lam0 form specializes to the lam0=2/9 form at lam0=2/9.

## Resolve before fix_loop

The Stage-27 normalization residual (the upstream expression that Stage 045 claims collapses to F_tr on the tracking branch) is not currently restated or imported by this script. To turn the F_tr check into a substantive verification, the new assertion must compute `<Stage-27 normalization residual>.subs(R_phi, R_U).subs(lam0, 2/9)` and compare against the notes' closed-form F_tr.

The question for the user:

> Which Stage-27 normalization-residual expression should the Stage-045 audit import as the upstream object that collapses to F_tr? Two possibilities:
>
> (a) The expression labeled "R_target = ..." in the Stage-27 or Stage-044 audit script (please specify which file and which assertion); the Stage-045 audit imports it as a symbolic expression (e.g., via shared module or by restating the expression once with a comment citing the upstream source) and then evaluates the collapse at `R_phi -> R_U`, `lam0 -> 2/9`.
>
> (b) The expression stated in notes/stages/moving_throat_pde_stage027_*.md (please specify file and section); the Stage-045 audit restates it inline once with a `# from notes/stages/moving_throat_pde_stage027_*.md sec N` comment.
>
> (c) Both Stage-27 sources differ on the residual form — flag for a deeper review of Stage 027.

The orchestrator will not invoke Codex on F3 until the user has chosen a direction and pointed to the canonical source.

Once resolved, the F3 patch instructions are:

In the SymPy script, replace lines 178–189 with:
```python
# Stage-27 normalization residual on the rank-2 closure, carried forward from
# <USER-SPECIFIED SOURCE>.  On the coherent tracking branch (R_phi -> R_U) at
# the D/N value lam0 = 2/9, the residual collapses to R_target - F_tr.
R_target_residual_stage27 = <expression from user-specified source>
collapsed_residual = sp.simplify(
    R_target_residual_stage27.subs(R_phi, R_U).subs(lam0, sp.Rational(2, 9))
)
F_tr_expected = sp.simplify(
    (9 * delta + (9 + 2 * R_U ** 2) * xi) ** 2
    * (9 * delta + (9 + 2 * R_U) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta ** 2 + 18 * delta * xi + (9 + 2 * R_U ** 2) * xi ** 2) ** 2)
)
expect_zero(
    "F_tr collapse from Stage-27 residual",
    collapsed_residual - (R_target - F_tr_expected),
)
print("coherent normalization residual =", sp.simplify(R_target - F_tr_expected))
```

In the Mathematica script, replace lines 125–136 with the corresponding `Together`/`FullSimplify`-based derivation route — concretely, the wl version should derive the collapsed residual via a different algebraic manipulation (e.g., factoring `Numerator[Together[stage27Residual /. rPhi -> rU /. lambda0 -> 2/9]]` and comparing) so that F1's independence requirement is also satisfied at this check.

**Verification command:** After Codex applies (post-resolution), the verifier will run `redteam exec-sympy 045` and `redteam exec-mathematica 045` and confirm:
- a new assertion line `F_tr collapse from Stage-27 residual = 0` appears in both transcripts,
- the upstream Stage-27 residual is referenced (not hand-written from scratch),
- both scripts exit 0.

## F4 — paper_misalignment

**Subtype:** notes_contradicts_script (labeling drift between paper card and script docstring/banners; script side is fixable without user gate)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_045.tex:1` quote: "Stage 045: Coherent Local D/N Support Kernel and Exact Tracking Reduction"
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex:68` quote: "045 & Coherent local D/N support kernel"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3` quote: `Stage 28 SymPy audit.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:31` quote: `banner("STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT")`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:191` quote: `print("\nAll Stage-28 symbolic checks passed.")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26` quote: `banner["STAGE 028 — COHERENT LOCAL TRACKING"]`

Direction of resolution: paper is authoritative; the script side has stale "Stage 28" / "Stage 028" labels from before the renumbering. Apply mechanically without user gate.

**Required change:**

Apply the following four single-line edits (no math content changes):

1. `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3`
   - Before: `Stage 28 SymPy audit.`
   - After: `Stage 045 SymPy audit.`

2. `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:31`
   - Before: `banner("STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT")`
   - After: `banner("STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT")`

3. `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:191`
   - Before: `print("\nAll Stage-28 symbolic checks passed.")`
   - After: `print("\nAll Stage-045 symbolic checks passed.")`

4. `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26`
   - Before: `banner["STAGE 028 — COHERENT LOCAL TRACKING"];`
   - After: `banner["STAGE 045 — COHERENT LOCAL TRACKING"];`

The wl footer at line 139 (`Print["Stage 045 Mathematica audit passed."]`) is already correct; do not change it.

**Verification command:** After Codex applies, the verifier will run `redteam exec-sympy 045` and `redteam exec-mathematica 045` and confirm:
- both `.txt` outputs banner with "STAGE 045",
- the SymPy footer reads "All Stage-045 symbolic checks passed.",
- the wl footer continues to read "Stage 045 Mathematica audit passed.",
- both scripts exit 0.

## Applied: F3

files_changed: scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:178-207; mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:125-156
summary: Imported Stage-044 F_cont/R_target residual into stage 045. Substituting tracking (R_phi → R_U) + D/N (lambda_0 → 2/9) yields notes' F_tr form. Generic tracking-collapse check kept as subsidiary anchor. Per user-approved Q2 (a) in batch III.1 v2. Stage 044 source destination_verified: scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:83,87,140,146.
deviation: none

## Applied: F4

files_changed: scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3,31,209; mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26; notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:232
summary: Relabeled docstring, banner, final-print, and notes heading from "Stage 28/028" to "Stage 045" per user-approved Q3 (b) in batch III.1 v2.
deviation: none
