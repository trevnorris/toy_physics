---
unit_id: 172
batch: V.1
created_at: 2026-05-28T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 172

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:42-77`

**Issue:** The `.wl` script reaches every slope (`deltaU2`, `deltaP0`, `deltaU2Star`, `deltaU4Star`) by the same `Normal[Series[...,{eps,0,1}]]` truncation divided by `eps*lam` that the SymPy script uses, with the same intermediate variable names. It is a line-by-line port, not an independent re-derivation, so it cannot catch an algebra/transcription error that was copied into both engines. The fix is to obtain the slopes by an independent route — implicit differentiation of the defining relations — while keeping the seven final residual assertions (and their printed forms) identical, so the verified claim does not change.

**Required change (independent re-derivation route):**

1. Replace the `deltaU2`/`deltaP0` extraction (currently wl L42-43):

   Before:
   ```
   deltaU2 = FullSimplify[(Normal[Series[u2A, {eps, 0, 1}]] - u2)/(eps*lam), Assumptions -> $Assumptions];
   deltaP0 = FullSimplify[(Normal[Series[p0A, {eps, 0, 1}]] - p0)/(eps*lam), Assumptions -> $Assumptions];
   ```
   After (implicit differentiation of `u2A*dA0 + dA2 == 0` and `p0A*dA0 - nA0 == 0` with respect to the perturbation parameter `t = eps*lam`):
   ```
   (* Independent route: differentiate the implicit defining relations, not the explicit quotient. *)
   du2sol = Solve[D[(u2 + t*deltaU2)*(d0 + t*dD0) + (d2 + t*dD2), t] == 0 /. t -> 0, deltaU2];
   deltaU2 = FullSimplify[deltaU2 /. First[du2sol], Assumptions -> $Assumptions];
   dp0sol = Solve[D[(p0 + t*deltaP0)*(d0 + t*dD0) - (n0 + t*dN0), t] == 0 /. t -> 0, deltaP0];
   deltaP0 = FullSimplify[deltaP0 /. First[dp0sol], Assumptions -> $Assumptions];
   ```
   (Here `d2 = -u2*d0` is already defined at L32; this uses the relation `u_2 = -D_2/D_0` <=> `u_2 D_0 + D_2 = 0`, and `P_0 = N_0/D_0` <=> `P_0 D_0 - N_0 = 0`, differentiated at the base point.)

2. Replace the `deltaU2Star`/`deltaU4Star` extraction (currently wl L76-77):

   Before:
   ```
   deltaU2Star = FullSimplify[(Normal[Series[u2AStar, {eps, 0, 1}]] - u2Star)/(eps*lam), Assumptions -> $Assumptions];
   deltaU4Star = FullSimplify[(Normal[Series[u4AStar, {eps, 0, 1}]] - u4Star)/(eps*lam), Assumptions -> $Assumptions];
   ```
   After (implicit differentiation of the `u_2` and `u_4` defining relations on the canonical branch, with `d2Star`, `d4Star` already defined at L67-68):
   ```
   du2StarSol = Solve[D[(u2Star + t*deltaU2Star)*(d0 + t*dD0) + (d2Star + t*dD2), t] == 0 /. t -> 0, deltaU2Star];
   deltaU2Star = FullSimplify[deltaU2Star /. First[du2StarSol], Assumptions -> $Assumptions];
   du4StarSol = Solve[
     D[(u4Star + t*deltaU4Star)*(d0 + t*dD0)^2 - ((d2Star + t*dD2)^2 - (d0 + t*dD0)*(d4Star + t*dD4)), t] == 0 /. t -> 0,
     deltaU4Star];
   deltaU4Star = FullSimplify[deltaU4Star /. First[du4StarSol], Assumptions -> $Assumptions];
   ```
   (Uses `u_4^(A) = (D_{A,2}^2 - D_{A,0} D_{A,4})/D_{A,0}^2` <=> `u_4 D_0^2 - (D_2^2 - D_0 D_4) = 0`, differentiated at the base point.)

3. Leave everything else unchanged: the `kA`/`gA` definitions (L48-49), all seven `expectZero[...]` checks (L52-101), the print statements, and `Exit[0]`. Do NOT change the SymPy script for F1.

**Self-test (auditor pre-checked):**
- `deltaU2` solved from `D[(u2 + t*deltaU2)(d0 + t*dD0) + (d2 + t*dD2), t]=0` at `t=0` gives `u2*dD0 + d0*deltaU2 + dD2 = 0` => `deltaU2 = -(dD2 + u2*dD0)/d0` — identical to the current series result.
- `deltaP0`: `D[(p0+t*deltaP0)(d0+t*dD0) - (n0+t*dN0),t]=0` at 0 => `p0*dD0 + d0*deltaP0 - dN0 = 0` => `deltaP0 = (dN0 - p0*dD0)/d0 = (d0*dN0 - n0*dD0)/d0^2` — identical.
- `deltaU4Star`: differentiating `u4 d0^2 - (d2^2 - d0 d4)` at the canonical base reproduces `-(5*dD0 + 18*dD2 + 81*dD4)/(81*d0)` — identical to the current series result (verified by hand in the report).
- The seven residual checks then evaluate to the same zeros. No `paper_misalignment` introduced (same constants, same final forms).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 172`, confirm the `Normal[Series[...]]` slope extractions are gone (replaced by `Solve`/`D` on the implicit relations), that the printed `deltaU2`/`deltaP0`/`deltaU4Star` forms are unchanged, that all seven `PASS:` lines remain, and that the script exits 0.

## Applied: F1

files_changed: mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl
summary: Replaced the `deltaU2`/`deltaP0` and `deltaU2Star`/`deltaU4Star` `Normal[Series[...]]` slope extractions with independent implicit-differentiation routes (`Solve`/`D` on the defining relations), matching the directive's "after" blocks verbatim; SymPy script left unchanged.
deviation: none

## F2 — stale_banner_label

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py:31` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:26`

**Issue:** Both scripts print a banner labeled `STAGE 155` although this is stage 172. The label is print-only and affects no assertion, but it mislabels the transcript.

**Required change:**
- In the SymPy script L31, change `banner("STAGE 155 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM")` to `banner("STAGE 172 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM")`.
- In the Mathematica script L26, change `banner["STAGE 155 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM"];` to `banner["STAGE 172 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM"];`.
- Change nothing else.

**Verification command:**
The verifier re-runs both engines; the first banner line in each `.txt` transcript must read `STAGE 172 — PHYSICAL SLOPE COLLAPSE...`.

## Applied: F2

files_changed: scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py, mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl
summary: Changed the stale `STAGE 155` banner label to `STAGE 172` in both the SymPy and Mathematica scripts.
deviation: none
