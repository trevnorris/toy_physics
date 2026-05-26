---
unit_id: 043
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 043 (batch III.1 v2)

Note on numbering: the auditor's report labels the findings F1 (`mathematica_transliteration`) and F2 (`paper_misalignment`, D_phi sign). The directive reverses the order (F1 = paper_misalignment, F2 = mathematica_transliteration) and the apply session at directive line 180 reuses an `## Applied: F2` header that actually documents the paper_misalignment fix. I use the auditor-report numbering below (F1 = transliteration, F2 = D_phi sign) and reconcile against the directive blocks.

This replaces an earlier batch III.1 verification of stage 043 that was tied to a different directive state (three findings: F1 transliteration, F2 / F3 tautological_check). The current report covers the active batch III.1 v2 directive.

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The Mathematica audit (`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:162-181`) gains a new block of three `expectZero[...]` checks immediately after the existing `mismatch formula` assertion (line 160), with no removals or restructuring of prior checks. The block corresponds to the directive's Insertion 2 (limit-based mismatch sign anchors). The earlier Insertion 1 (series-resummation derivation of `R_phi`) was correctly blocked by the auditor in the directive itself due to a known algebra error in the proposed `rPhiResummed` formula; Codex appended `## Blocked: F2-Insertion1` and skipped that insertion (`/var/projects/toy_physics/research/pde_ledger/redteam/directives/stage_043.md:161-164`). iter2 rewrote the test-triple substitutions to use primitive coupling substitutions (`gS -> 0, gW -> gU*gR/kU` for the (sigma0=0, rho0=1) triple; explicit `gU -> 1, kU -> 1, gB -> 1, gS -> 1, gR -> 0, gW -> 1` for the (sigma0=1, rho0=0) triple; `gS -> gB*gR/gW` for tracking) because `rPhi` and `rU` in the script are expressed in primitive couplings rather than in (sigma0, rho0) directly — the directive's iter1 substitutions `sigma0 -> sNum, rho0 -> rNum` would not flow back into `rPhi`/`rU`.

**Assessment:**
The substitutions are correct. Tracing them by hand:
- TestPoint1 `{deltaU -> 1, gS -> 0, gW -> gU*gR/kU}`: yields `sigma0 = gU*gS/(kU*gB) = 0` and `rho0 = gU*gR/(kU*gW) = gU*gR/(kU * gU*gR/kU) = 1`. With `deltaU=1`: `rPhi = (1+0/2)/1 = 1`, `rU = (1+1/2)/2 = 3/4`, so `rPhi-rU = 1/4`. Anchor expected `1/4`. PASS confirmed in exec log line 81.
- TestPoint2 `{deltaU -> 1, gU -> 1, kU -> 1, gB -> 1, gS -> 1, gR -> 0, gW -> 1}`: yields `sigma0 = 1`, `rho0 = 0`. `rPhi = (1+1/2)/2 = 3/4`, `rU = (1+0)/1 = 1`, `rPhi-rU = -1/4`. Anchor expected `-1/4`. PASS confirmed in exec log line 83.
- Tracking `{gS -> gB*gR/gW}`: enforces `sigma0 = gU*(gB*gR/gW)/(kU*gB) = gU*gR/(kU*gW) = rho0`. Anchor expects `(rPhi - rU) = 0` symbolically. PASS confirmed in exec log line 85.

The three new assertions are non-tautological. Crucially, they do NOT compare `mismatch` against `mismatchExpected`; they evaluate `rPhi - rU` (the raw rank-1 ratio difference built from the body of the script, not the closed-form expected) at fixed points and assert specific rational values (1/4, -1/4, 0). A sign error in `mismatchExpected` would not be hidden — these anchors catch any sign error in the body computation of the ratio difference, independent of the closed form. This is what makes them a structurally-distinct route from the symbolic `mismatch - mismatchExpected` check, which is the directive's intent.

Variable names (`mismatchAtTestPoint1`, `mismatchAtTestPoint2`, `mismatchAtTracking`) do not mirror SymPy names; no prior Mathematica assertion was removed or restructured. The deviation noted in `## Applied: F2-Insertion2-iter2` (using primitive coupling substitutions rather than the directive's `sigma0 -> sNum, rho0 -> rNum` placeholders) is forced by the script's structure (rho0 is *defined* as `gU*gR/(kU*gW)` and substituting into the symbol `rho0` after definition would not propagate into `rU`) and is a faithful realization of the directive's intent.

### F2 — paper_misalignment (D_phi sign)

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-53` changes (per the diff at `redteam/exec_logs/stage_043_diff.patch`):
- `dPhi = FullSimplify[Det[{{y0, y1}, {kappa0, kappa1}}], ...];` -> `dPhi = FullSimplify[Det[{{kappa0, kappa1}, {y0, y1}}], ...];`
- `dPhiExpected = FullSimplify[kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), ...];` -> `dPhiExpected = FullSimplify[-kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), ...];`

This is the directive's option (a): swap determinant rows so `dPhi = kappa0*y1 - kappa1*y0` (matching the notes' definition `D_phi := kappa_0 y_1 - kappa_1 y_0` and the SymPy script's `Dphi = sp.factor(sp.expand(kappa0 * y1 - kappa1 * y0))` at `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:77`), and add the leading minus sign to `dPhiExpected`. Applied per user-approved Q1 (a) in batch III.1 v2, recorded in the `## Applied: F2` block at directive line 180.

**Assessment:**
The change is correct. Substituting `y0 = kappa0*(gB + gS*gU/kU)`, `y1 = kappa1*(gB + gS*gU/(kU*(1+deltaU)))` into the new `dPhi = kappa0*y1 - kappa1*y0` gives `kappa0*kappa1*[(gB + gS*gU/(kU*(1+deltaU))) - (gB + gS*gU/kU)] = kappa0*kappa1*gS*gU/kU * [1/(1+deltaU) - 1] = -kappa0*kappa1*gS*gU*deltaU/(kU*(1+deltaU))`, which matches the new `dPhiExpected = -kappa0*kappa1*gB*sigma0*deltaU/(1+deltaU)` after substituting `sigma0 = gU*gS/(kU*gB)`. ✓

Mathematica exec log line 16 confirms the printed `D_phi = -((deltaU*gS*gU*kappa0*kappa1)/(kU + deltaU*kU))`, which now matches the SymPy log line 28-30 `-δ_U·g_S·g_U·κ₀·κ₁ / (K_U·(δ_U + 1))`. The two engines now print the same sign for the named `D_phi` quantity, and the notes' convention is honored on both sides. No collateral edits.

## Exec log assessment

**SymPy:** exit=0 (`redteam/exec_logs/stage_043_sympy.log:126`). Notable lines:
- `D_phi = -δ_U·g_S·g_U·κ₀·κ₁ / (K_U·(δ_U + 1))` (lines 28-30) — sign convention unchanged (SymPy was already correct per the notes).
- `D_phi - expected = 0` (line 31).
- `mismatch formula = 0` (line 108).
- All assertions PASS, exit 0.

The SymPy script was not edited; the SymPy log is unchanged from pre-fix and is provided for cross-engine sign comparison on `D_phi`.

**Mathematica:** exit=0 (`redteam/exec_logs/stage_043_mathematica.log:88`). Notable lines:
- `D_phi = -((deltaU*gS*gU*kappa0*kappa1)/(kU + deltaU*kU))` (line 16) — sign now matches SymPy and notes after F2 fix.
- `PASS: D_phi - expected` (line 20) — internal check still passes with the new sign-flipped `dPhiExpected`.
- `PASS: mismatch sign at deltaU=1, sigma0=0, rho0=1` (line 81) — new F1-Insertion2 anchor.
- `PASS: mismatch sign at deltaU=1, sigma0=1, rho0=0` (line 83) — new F1-Insertion2 anchor.
- `PASS: mismatch vanishes at tracking sigma0=rho0` (line 85) — new F1-Insertion2 anchor.
- Two `Limit::alimv` warnings (lines 28, 40) on pre-existing `Limit[...]` calls — unrelated to either finding.
- exit 0.

**Output freshness:** The exec logs are fresh (dated 2026-05-26T01:56:47). However, the saved `.txt` outputs under `mathematica/output/` and `scripts/output/` have mtimes of 2026-05-22 12:40, which are older than the current script mtime (2026-05-26 01:55 for the `.wl`). The exec logs themselves are fresh and correctly reflect the post-fix state — this freshness gap is a housekeeping issue for the saved `.txt` mirrors, not a verification blocker. The orchestrator should regenerate the saved `.txt` outputs to bring them in sync.

## Material-change assessment

`material_change`: false.

Reasoning: F2 is a sign-of-printed-quantity fix internal to the Mathematica script. The `D_phi = 0 iff sigma0 = 0 or delta_U = 0` vanishing condition is symmetric under sign and unchanged. The SymPy script (which downstream code might import or mirror) was untouched. F1 only adds three new `expectZero[...]` checks that test the existing computed `rPhi - rU` at three fixed points; no derived quantity used by downstream stages was modified. Specifically, `rPhi`, `rPhiExpected`, `dPhiZ`, `rho0`, `rU`, `mismatch`, and `mismatchExpected` are all unchanged. The `dPhi` printed name and `dPhiExpected` are sign-flipped relative to the previous (incorrect) Mathematica convention, but the SymPy convention (the canonical export) is unchanged. No downstream unit depends on the Mathematica printed sign of `D_phi`.

## Side observations (non-blocking)

- The `.txt` mirror outputs (`mathematica/output/...txt`, `scripts/output/...txt`) are stale relative to the current scripts. The exec logs are fresh; the `.txt` mirrors should be regenerated to keep them in sync. Not a verification concern.
- The directive's apply-block at line 180 reuses the `## Applied: F2` header that was already used at line 166 for the transliteration (Insertion 2) fix. The block at line 180 is actually the F1-in-directive / F2-in-report paper_misalignment D_phi sign fix from the Q1 apply session. The two blocks are textually distinct and unambiguously cover the two separate edits, so the labeling collision does not affect correctness, but it could confuse a future reader.
- F1-Insertion2 TestPoint2 uses seven explicit numeric substitutions (`gU -> 1, kU -> 1, gB -> 1, gS -> 1, gR -> 0, gW -> 1`) rather than the minimal substitutions that would force (sigma0=1, rho0=0). The result is unambiguous and matches the directive's self-test, so this is acceptable; a leaner formulation could use only the substitutions required to fix sigma0 and rho0 independently of the other couplings, but the present form is correct.
- Two `Limit::alimv` warnings remain in the Mathematica transcript on pre-existing `Limit[...]` calls (lines 28, 40 of the exec log). These are cosmetic; the limits evaluate correctly and residuals are 0.

## Verdict justification

Both findings are resolved. F2 (paper_misalignment, D_phi sign) was correctly fixed by the row-swap and sign-flip at `wl:52-53` per the user-approved option (a); the Mathematica printed `D_phi` now matches SymPy and the notes' convention. F1 (mathematica_transliteration) was addressed via Insertion 2 (three numeric-point sign anchors on `rPhi - rU`) after Insertion 1 was correctly blocked by the auditor due to an algebra error in the proposed series-resummation formula; iter2 corrected the iter1 realization to use primitive coupling substitutions that actually achieve the (sigma0=0, rho0=1), (sigma0=1, rho0=0), and (sigma0=rho0) test triples. The new anchors are non-tautological (they evaluate `rPhi - rU` from primitive couplings and assert specific rational values that catch a sign error independent of `mismatchExpected`), do not mirror SymPy variable names, and pass cleanly. Both engines exit 0. No collateral edits in the diff; no regressions in the exec logs.
