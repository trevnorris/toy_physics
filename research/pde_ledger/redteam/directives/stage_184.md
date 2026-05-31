---
unit_id: 184
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-30T01:29:37-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 184

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 and F2 share a single root cause and a single fix (route the Mathematica drift quantities through the already-defined branch composites). Applying F1 resolves F2 automatically; confirm F2's verification afterward.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch `paper/`, `notes/`, or any prose documents. The SymPy script is correct and must NOT be changed.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file `expectZero` checks passing. The orchestrator independently re-runs afterward.

## F1 — tautological_check (resolves F2 — insufficient_verification)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl:53,60,65,71`

**Issue:** The four load-bearing drift quantities are hardcoded to the expressions they are compared against, making the drift-law assertions tautological (`expr - expr == 0`), and leaving the branch composites `tTr/tTr0` (51-52), `nTr/nTr0` (58-59), `eComp/eComp0` (69-70), and `epsEtaVar` (46) as dead code. The Mathematica engine must instead derive each drift independently from its composite, exactly as the SymPy script does at lines 65-90 (which genuinely uses `sp.series(sp.log(composite/reference), small, 0, 2)`). The derived values must come out to the same closed forms already printed in the transcript, but they must be PRODUCED by the series operation, not restated.

**Required change:**

Use Mathematica's own series machinery (`SeriesCoefficient[Log[composite/reference], {small, 0, 1}]`) so each drift is genuinely a first-order log-drift of the composite. Make exactly these four substitutions (composite definitions `tTr/tTr0`, `nTr/nTr0`, `eComp/eComp0`, `epsEtaVar` are already defined above each site — do not re-define them):

1. Line 53 — tracking drift. Before:
   ```
   dlnTtr = FullSimplify[-cStar*theta1, Assumptions -> $Assumptions];
   ```
   After:
   ```
   dlnTtr = FullSimplify[SeriesCoefficient[Log[tTr/tTr0], {small, 0, 1}], Assumptions -> $Assumptions];
   ```

2. Line 60 — nontracking drift. Before:
   ```
   dlnNtr = FullSimplify[xi1 + bStar*theta1, Assumptions -> $Assumptions];
   ```
   After:
   ```
   dlnNtr = FullSimplify[SeriesCoefficient[Log[nTr/nTr0], {small, 0, 1}], Assumptions -> $Assumptions];
   ```

3. Line 65 — dressing drift. Before:
   ```
   dlnEpsEta = FullSimplify[sigmaEta, Assumptions -> $Assumptions];
   ```
   After:
   ```
   dlnEpsEta = FullSimplify[SeriesCoefficient[Log[epsEtaVar/epsEta], {small, 0, 1}], Assumptions -> $Assumptions];
   ```

4. Line 71 — selected-branch complement drift. Before:
   ```
   dlnEcomp = FullSimplify[-epsEta*sigmaEta/(1 - epsEta), Assumptions -> $Assumptions];
   ```
   After:
   ```
   dlnEcomp = FullSimplify[SeriesCoefficient[Log[eComp/eComp0], {small, 0, 1}], Assumptions -> $Assumptions];
   ```

Leave the `expectZero[...]` lines (55, 62, 67, 73, 76-78), the composite definitions (46, 51-52, 58-59, 69-70), the zero-map block (80-85), and everything else unchanged. The zero-map lines 80-85 will now substitute into genuinely-derived drifts, which is exactly the intent.

Optional (non-blocking) label fix: line 26 banner reads `"STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES"` but this is stage 184. If you touch it, change to `"STAGE 184 — EXACT BRANCH-INVARIANT COORDINATES"`. The SymPy script line 35 has the same stale banner; you may NOT edit the SymPy script under this directive, so leave it — flag in your Applied note if you want it tracked. Do not let this optional change gate the run.

**Expected post-fix derived values** (must match what the transcript already prints):
- `dlnTtr` -> `-(((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)*theta1)/(chi0*deltaU))`
- `dlnNtr` -> `(2*(1 + chi0 + deltaU)*theta1)/deltaU + xi1`
- `dlnEpsEta` -> `sigmaEta`
- `dlnEcomp` -> `-((epsEta*sigmaEta)/(1 - epsEta))` (equivalently `(epsEta*sigmaEta)/(-1 + epsEta)`)

All four `expectZero` drift-law checks (lines 55, 62, 67, 73), the three mirror checks (76-78), and the three zero-map checks (83-85) must still report `= 0` and `PASS`, and the script must `Exit[0]`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 184` and confirm: (a) `dlnTtr/dlnNtr/dlnEpsEta/dlnEcomp` are now assigned via `SeriesCoefficient[Log[.../...], {small, 0, 1}]` (genuinely consuming `tTr/tTr0`, `nTr/nTr0`, `epsEtaVar`, `eComp/eComp0`), so no composite remains write-only (F2); (b) the printed drift values are unchanged; (c) every `expectZero` PASSes and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`
- summary: Routed all four Mathematica drift quantities through `SeriesCoefficient[Log[...], {small, 0, 1}]` of the existing branch composites, which also resolves F2's shared insufficient-verification root cause.
- deviation: none
