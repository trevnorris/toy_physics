---
unit_id: 043
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T07:55:56Z
findings_applied: 1
findings_blocked: 1
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 043

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage043_support_direction_extraction.md:119-121` quote:
  ```
  `D_phi := kappa_0 y_1 - kappa_1 y_0`

  `      = - kappa_0 kappa_1 g_B sigma_0 delta_U / (1 + delta_U).`
  ```

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-53` quote:
  ```
  dPhi = FullSimplify[Det[{{y0, y1}, {kappa0, kappa1}}], Assumptions -> $Assumptions];
  dPhiExpected = FullSimplify[kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];
  ```
  Note: `Det[{{y0, y1}, {kappa0, kappa1}}] = y0*kappa1 - y1*kappa0 = -(kappa_0 y_1 - kappa_1 y_0)`, the negative of the notes' definition. The companion SymPy script (`scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:77-78`) honors the notes' sign convention:
  ```
  Dphi = sp.factor(sp.expand(kappa0 * y1 - kappa1 * y0))
  Dphi_expected = sp.simplify(-kappa0 * kappa1 * gB * sigma0 * delta_U / (1 + delta_U))
  ```

## Resolve before fix_loop

The notes define `D_phi := kappa_0 y_1 - kappa_1 y_0` and give the explicit value with a leading minus sign. The SymPy script matches the notes' convention. The Mathematica script uses the OPPOSITE sign convention (`Det[{{y0,y1},{kappa0,kappa1}}]` = `-(notes' D_phi)`) and a sign-flipped expected. Both checks pass internally, but the printed `D_phi` in the Mathematica transcript has the opposite sign from the SymPy transcript and the notes. The vanishing condition `D_phi = 0 iff sigma_0 = 0 or delta_U = 0` is unaffected, but the named/printed quantity disagrees across engines and against the notes.

Which sign is canonical?

Possible directions (the user picks one):
- (a) Notes are canonical (`D_phi := kappa_0 y_1 - kappa_1 y_0`, leading minus sign in the explicit formula). Change Mathematica `wl:52` to `dPhi = FullSimplify[Det[{{kappa0, kappa1}, {y0, y1}}], Assumptions -> $Assumptions];` (equivalent to `kappa_0*y_1 - kappa_1*y_0`) AND change `wl:53` to `dPhiExpected = FullSimplify[-kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];`. No paper edit. Re-run mathematica.
- (b) Mathematica's convention is canonical (`D_phi := y_0 kappa_1 - y_1 kappa_0`, no leading minus). Update notes Section 3 (lines 119–121) AND SymPy `py:77-78` so all three sides agree. This requires a paper/notes-side edit; Codex will not apply either until directed.
- (c) The sign is conventional and a sign-flip is acceptable. Add a comment to the Mathematica script noting the convention difference. Discouraged: leaves the transcripts visibly disagreeing on a named quantity.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`

**Issue:**
The Mathematica script is structurally a line-for-line port of the SymPy script (`scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`). Both engines build `y = gB v + gU gS D_U v`, define `sigma_0 = gU gS/(KU gB)`, take the ratio `(y1/k1)/(y0/k0)`, and verify against the same closed-form `(1 + sigma0/(1+delta_U))/(1+sigma0)`. Both substitute `kappa_1^2 -> (2/11) sigma`, `kappa_0^2 -> (9/11) sigma` in the overlap. Both isolate `M_supp` via the same free-baseline `B = kappa_0^2` substitution and the same `B = 8/pi^2` evaluation. Variable names are case-only renames (`DU`↔`dU`, `Rphi`↔`rPhi`, `Dphi`↔`dPhi`, `Dzphi`↔`dPhiZ`, etc.). The existing Mathematica endpoint checks (deltaU=0, deltaU->Infinity) and the leading-order Series check on the mismatch are genuine value-adds but are decorations on top of the SymPy recipe, not a second derivation route. A transliterated second engine cannot catch a hidden bug in the shared recipe.

**Required change:**
Add TWO new `expectZero[...]` blocks to the Mathematica script. Each derives a load-bearing identity by a structurally distinct route from both the existing Mathematica code and the SymPy code. Do NOT remove or restructure any existing check.

### Insertion 1: Series-expansion derivation of R_phi (insert after line 59, before line 61)

Insert exactly the following block immediately after the existing `expectZero["D_phi - expected", dPhi - dPhiExpected];` (line 59), keeping a blank line of separation:

```
(* F2 cross-check 1: derive R_phi by Series expansion in sigma0, then resum.
   This route does NOT use the closed-form formula (1 + sigma0/(1+deltaU))/(1+sigma0)
   to construct rPhiExpected; it computes the Taylor coefficients of R_phi(sigma0)
   directly from the ratio (y1/kappa1)/(y0/kappa0), where y is parameterized by sigma0,
   then sums the resulting geometric series in closed form. *)
ClearAll[s];
yParam = {kappa0 gB (1 + s), kappa1 gB (1 + s/(1 + deltaU))};
rPhiOfS = FullSimplify[(yParam[[2]]/kappa1) / (yParam[[1]]/kappa0), Assumptions -> $Assumptions];
rPhiSeries = Normal@Series[rPhiOfS, {s, 0, 6}];
(* The closed-form expected, written as a single rational in s, should equal the
   resummation of the geometric series 1 + (sigma0/(1+deltaU) - sigma0)/(1+sigma0)
   = 1 - (deltaU/(1+deltaU)) sum_{k>=1} (-sigma0)^k. *)
rPhiResummed = FullSimplify[
  1 - (deltaU/(1 + deltaU)) Sum[(-s)^k, {k, 1, Infinity}],
  Assumptions -> $Assumptions && Abs[s] < 1];
expectZero["R_phi via Series resummation matches closed form (in s)",
  FullSimplify[rPhiOfS - rPhiResummed, Assumptions -> $Assumptions]];
(* Cross-check: Series truncation of rPhiOfS to order s^6 should match the corresponding
   truncation of the closed-form (1 + sigma0/(1+deltaU))/(1+sigma0). *)
rPhiClosedSeries = Normal@Series[(1 + s/(1 + deltaU))/(1 + s), {s, 0, 6}];
expectZero["R_phi Taylor series matches closed-form Taylor series to O(s^6)",
  FullSimplify[rPhiSeries - rPhiClosedSeries, Assumptions -> $Assumptions]];

```

Self-test for Insertion 1 (Codex should verify mentally before writing):
- `yParam` is just the notes' factored form of y with `sigma0 -> s`. Both the existing `y` and `yParam` reduce to the same expression after `sigma0 = gU gS/(KU gB)`. So `rPhiOfS` is the same R_phi but parameterized in `s`.
- `Sum[(-s)^k, {k, 1, Infinity}] = -s/(1+s)` for `|s| < 1`, so `rPhiResummed = 1 - (deltaU/(1+deltaU)) * (-s/(1+s)) = 1 + deltaU s / ((1+deltaU)(1+s))`. Verify this equals `(1 + s/(1+deltaU))/(1+s)`:
  - `(1 + s/(1+deltaU))/(1+s) = [(1+deltaU+s)/(1+deltaU)]/(1+s) = (1+deltaU+s)/[(1+deltaU)(1+s)]`.
  - `1 + deltaU s/[(1+deltaU)(1+s)] = [(1+deltaU)(1+s) + deltaU s]/[(1+deltaU)(1+s)] = [1+deltaU + s + deltaU s + deltaU s]/[...] = [1+deltaU + s + 2 deltaU s]/[...]`. Hmm, that's not equal. Let me recompute.
  - Actually `rPhiResummed = 1 - (deltaU/(1+deltaU)) * (-s/(1+s)) = 1 + (deltaU s)/[(1+deltaU)(1+s)]`. And `closed = (1+deltaU+s)/[(1+deltaU)(1+s)]`. Their difference: `closed - rPhiResummed = (1+deltaU+s)/[(1+deltaU)(1+s)] - 1 - (deltaU s)/[(1+deltaU)(1+s)] = [(1+deltaU+s) - (1+deltaU)(1+s) - deltaU s]/[(1+deltaU)(1+s)] = [1+deltaU+s - 1 - s - deltaU - deltaU s - deltaU s]/[...] = [-2 deltaU s]/[(1+deltaU)(1+s)]`. NOT zero.

The series-resummation route in Insertion 1 has an algebra error in the proposed `rPhiResummed` formula. The correct resummation of `R_phi = (1 + s/(1+deltaU))/(1+s)` as a power series in s is `(1+s/(1+deltaU)) * (1 - s + s^2 - s^3 + ...)` = `sum_{k>=0} (-s)^k * (1 + s/(1+deltaU))`, which is more involved. The simple "1 + deltaU s / (...)" form does not resum cleanly.

**Codex: skip Insertion 1.** The series-resummation route is more subtle than this directive can mechanically pin down. Append `## Blocked: F2-Insertion1` with the question:

> "The proposed series-resummation route for R_phi has an algebra error: `1 - (deltaU/(1+deltaU)) * sum_{k>=1}(-s)^k = 1 + deltaU s/[(1+deltaU)(1+s)]`, which does NOT equal `(1 + s/(1+deltaU))/(1+s)`. What is the correct independent-derivation route?"

Apply only Insertion 2 below.

### Insertion 2: Limit-based derivation of the mismatch sign (insert after line 160, before line 162)

Insert exactly the following block immediately after the existing `expectZero["mismatch formula", mismatch - mismatchExpected];` (line 160), keeping a blank line of separation:

```
(* F2 cross-check 2: derive the sign of (R_phi - R_U) at sigma0 < rho0 and sigma0 > rho0
   by direct numerical-symbolic limits, independent of the closed-form mismatch formula.
   This catches a hidden sign error in the closed-form expected. *)
ClearAll[sNum, rNum];
rPhiAt = FullSimplify[rPhiExpected /. sigma0 -> sNum, Assumptions -> $Assumptions];
rUAt = FullSimplify[rU /. rho0 -> rNum, Assumptions -> $Assumptions];
(* At deltaU = 1, sNum = 0, rNum = 1: rho0 > sigma0, mismatch should be positive. *)
diffCase1 = FullSimplify[(rPhiAt - rUAt) /. {deltaU -> 1, sNum -> 0, rNum -> 1}, Assumptions -> $Assumptions];
(* expected: deltaU*(rho0-sigma0)/((1+deltaU)(1+sigma0)(1+rho0)) at (1, 0, 1) = 1*1/(2*1*2) = 1/4 *)
expectZero["mismatch sign at deltaU=1, sigma0=0, rho0=1", diffCase1 - 1/4];
(* At deltaU = 1, sNum = 1, rNum = 0: sigma0 > rho0, mismatch should be negative. *)
diffCase2 = FullSimplify[(rPhiAt - rUAt) /. {deltaU -> 1, sNum -> 1, rNum -> 0}, Assumptions -> $Assumptions];
(* expected: 1*(-1)/(2*2*1) = -1/4 *)
expectZero["mismatch sign at deltaU=1, sigma0=1, rho0=0", diffCase2 - (-1/4)];
(* At sigma0 = rho0 (tracking): mismatch must vanish for ANY deltaU. *)
diffTracking = FullSimplify[(rPhiAt - rUAt) /. {sNum -> rNum}, Assumptions -> $Assumptions];
expectZero["mismatch vanishes at tracking sigma0=rho0", diffTracking];

```

Self-test for Insertion 2:
- At `deltaU=1, sigma0=0, rho0=1`: `R_phi = (1 + 0/2)/(1+0) = 1`. `R_U = (1 + 1/2)/(1+1) = (3/2)/2 = 3/4`. `R_phi - R_U = 1 - 3/4 = 1/4`. ✓ matches expected.
- At `deltaU=1, sigma0=1, rho0=0`: `R_phi = (1 + 1/2)/(1+1) = 3/4`. `R_U = (1 + 0)/(1+0) = 1`. `R_phi - R_U = 3/4 - 1 = -1/4`. ✓ matches expected.
- At `sigma0 = rho0` (any deltaU): both R_phi and R_U are the same function of the common ratio, so difference is 0. ✓

This insertion exercises specific numeric points where the existing closed-form `mismatchExpected` could (in principle) have a sign error and the rest of the script would not catch it (because the test is `mismatch - mismatchExpected == 0` symbolically — a sign-flipped `mismatchExpected` and a sign-flipped `mismatch` would both yield 0). The numeric-point check anchors the absolute sign.

**Codex: apply only Insertion 2.** Skip Insertion 1 and append a `## Blocked: F2-Insertion1` block as above.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 043` and confirm:
1. The script still exits 0.
2. The transcript shows three new `PASS:` lines (`mismatch sign at deltaU=1, sigma0=0, rho0=1`, `mismatch sign at deltaU=1, sigma0=1, rho0=0`, `mismatch vanishes at tracking sigma0=rho0`) immediately after the existing `mismatch formula` PASS line.
3. The `## Blocked: F2-Insertion1` block is appended to F2 in this directive file, with the question quoted above.

## Blocked: F2-Insertion1

- reason: The proposed series-resummation route for R_phi has an algebra error in the supplied formula.
- question: "The proposed series-resummation route for R_phi has an algebra error: `1 - (deltaU/(1+deltaU)) * sum_{k>=1}(-s)^k = 1 + deltaU s/[(1+deltaU)(1+s)]`, which does NOT equal `(1 + s/(1+deltaU))/(1+s)`. What is the correct independent-derivation route?"

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- summary: Added the limit-based mismatch sign checks after the existing mismatch formula assertion.
- deviation: none

## Applied: F2-Insertion2-iter2

- files_changed:
  - `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- summary: Rewrote the mismatch sign anchors to realize the sigma0/rho0 test triples through primitive coupling substitutions.
- deviation: Used `gW -> gU*gR/kU` for `rho0 = 1` because the current script defines `rho0 = gU*gR/(kU*gW)`.

## Applied: F2

files_changed: mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-54
summary: Swapped dPhi determinant rows to Det[{{kappa0, kappa1}, {y0, y1}}] (kappa-first convention matching SymPy + notes per stage 039 D_dir precedent). Updated dPhiExpected sign. Per user-approved Q1 (a) in batch III.1 v2.
deviation: none
