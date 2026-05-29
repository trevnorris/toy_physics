---
unit_id: 173
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 173

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl:32-78`

**Issue:** The `.wl` is a line-by-line port of the `.py`: identical premise encoding (`d0A=d0+eps*lam*d01`, etc.), identical series-then-`(D[u2A,eps]/.eps->0)/lam` first-order-coefficient extraction (math:37-43 mirrors sympy:42-48), and identical sequential `Solve` choreography for `d21` and `d41` (math:67,71,74 mirrors sympy:84,88,91). A transliterated second engine cannot catch a SymPy-side bug because it reproduces the same algebra. The `.wl` must reach each result by a route distinct from the SymPy choreography while keeping the same final `expectZero` targets.

**Required change:**
Restructure the derivation path only (the six `expectZero` targets at math:49,60,61,65,77,78 must stay byte-for-byte the same — they encode the paper's closed forms). Concretely:

1. Replace the first-order-coefficient extraction at math:41-43

   ```
   u21 = FullSimplify[(D[u2A, eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];
   u41 = FullSimplify[(D[u4A, eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];
   p1  = FullSimplify[(D[p0A, eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];
   ```

   with direct series-coefficient extraction (a different route than differentiate-then-set-eps-zero), e.g.:

   ```
   u21 = FullSimplify[Coefficient[Series[-d2A/d0A, {eps, 0, 1}] // Normal, eps, 1]/lam, Assumptions -> $Assumptions];
   u41 = FullSimplify[Coefficient[Series[(d2A^2 - d0A*d4A)/d0A^2, {eps, 0, 1}] // Normal, eps, 1]/lam, Assumptions -> $Assumptions];
   p1  = FullSimplify[Coefficient[Series[n0A/d0A, {eps, 0, 1}] // Normal, eps, 1]/lam, Assumptions -> $Assumptions];
   ```

   This makes the `u2A,u4A,p0A` intermediates at math:37-39 unnecessary; you may remove math:37-39 and inline the expressions into the `Series[...]` calls above (do not change the underlying expressions `-d2A/d0A`, `(d2A^2-d0A*d4A)/d0A^2`, `n0A/d0A`).

2. Replace the sequential-`Solve` even-preserving choreography at math:67,71,74. Instead of `Solve[u41Can==8 u21Can/9, d41]` then `Solve[u21Can==0, d21]` then substituting, reach the even-preserving values by direct algebra distinct from the SymPy path:

   ```
   (* even-preserving condition u2^(1)=0 gives d21 directly *)
   u21ZeroD21 = FullSimplify[d21 /. First[Solve[u21Can == 0, d21]], Assumptions -> $Assumptions];
   ```
   keep this one (it is the natural Mathematica statement), but obtain `d41Even` by substituting the even-preserving `d21` into the *general* canonical `u41` relation and solving `u41Can == 8 u21Can/9` once with `d21` already fixed:
   ```
   d41Even = FullSimplify[
     d41 /. First[Solve[(u41Can == 8 u21Can/9) /. d21 -> u21ZeroD21, d41]],
     Assumptions -> $Assumptions];
   ```
   This removes the intermediate `d41Hidden` general-solve at math:67,74 that mirrors sympy:84,91. (If you prefer to keep a `d41Hidden` printout for parity with the transcript, that is acceptable, but the load-bearing `d41Even` must come from the route above, not from substituting into a previously stored `d41Hidden`.)

3. Leave the `expectZero` lines (math:49,60,61,65,77,78), the `xiLoad` definition (math:81), the lane-form prints (math:84-89), and the carry-forward prints (math:92-101) unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 173` and confirm: (a) the `.wl` no longer contains the `(D[u2A, eps] /. eps -> 0)/lam` extraction nor the general `d41Hidden = ... Solve[u41Can == 8*u21Can/9, d41]` step that mirrored sympy:42-48,84; (b) all six `expectZero` checks still PASS against the unchanged paper-form targets; (c) the script exits 0.

## Applied: F1

files_changed: mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl
summary: Replaced the differentiate-then-set-eps-zero first-order extraction with direct Series/Coefficient extraction (removing the u2A/u4A/p0A intermediates and inlining their expressions), and obtained d41Even by solving the canonical relation u41Can==8 u21Can/9 with d21 already fixed to u21ZeroD21 (removing the intermediate general d41Hidden Solve), leaving all six expectZero targets, xiLoad, lane-form, and carry-forward prints unchanged.
deviation: none

## F2 — stale_output (cosmetic banner mislabel)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl:26`

**Issue:** Both scripts open with `banner("STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT")`, but this unit is Stage 173. The wrong stage number propagates into both saved transcripts (output line 11). Cosmetic only; no math is affected.

**Required change:**
In both files, change the banner argument string from
`"STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT"`
to
`"STAGE 173 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT"`.
(sympy:30 uses `banner("...")`; math:26 uses `banner["..."]`. Edit only the string; keep the rest of the line.)

**Verification command:**
After Codex applies, the verifier re-runs both engines and confirms output line 11 reads "STAGE 173 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT" in both transcripts and both scripts exit 0.

## Applied: F2

files_changed: scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py, mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl
summary: Changed the banner string from "STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT" to "STAGE 173 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT" in both files.
deviation: none
