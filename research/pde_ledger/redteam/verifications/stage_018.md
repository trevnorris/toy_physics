---
unit_id: 018
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T22:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: true
---

# Verification — unit 018 (v2, post-trim)

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim)

**Classification:** resolved

**What changed:**

Per the user's Q4=b (TRIM) resolution, Codex trimmed the over-scoped algebra from both engines and kept only the Gaussian-profile bridge probes that anchor to the stage-018 paper card's `M_Σ` / `K_Σ` exports.

- SymPy (`scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`): deleted the block that previously held A1–A17 (one-pole numerator identity, two `KSigma` closures, compatibility cross-closure, even-gate 2x2 determinant + wall-stiffness/wall-inertia slope solve, `Xi1` residual + mutations). The new `main()` (file lines 20–34) declares `w`, `beta = exp(-w²/2)`, computes `MSigma_example = ∫ β² dw` and `KSigma_example = ∫ ((β')² + β²) dw`, and asserts both equal `√π` and `3√π/2` (lines 25–26). Print banner updated (line 29) from the old "isotropic compatibility / wall-slope solve / Xi1" string to "concrete Gaussian parent-wall inertia and stiffness branch integrals."
- Mathematica (`mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl`): header claims block (lines 1–5) collapsed from M1–M8 to just "M1 Gaussian wall inertia and stiffness integrals." Body M1–M7-mut blocks deleted; surviving block (file lines 14–27) holds the Gaussian inertia and stiffness checks, renumbered M1, with `If[FullSimplify[...] =!= 0, ...; Exit[1]]` failure guards retained.
- Diff (`exec_logs/stage_018_diff.patch`) confirms the edits are confined to these two files and consist of deletions plus the docstring/banner relabel — no algebraic substitutions and no smuggled-in new assertions.

**Assessment:**

The trim matches directive option (b)'s TRIM branch: drop A1–A17 / M1–M7-mut from stage 018 and leave the Gaussian bridge probes. The surviving assertions (sympy A18/A19 ↔ math M1a/M1b) directly exercise the paper card's `\stagefield{Output}` integrals `M_Σ = ∫ μ_η β² dw` and `K_Σ = ∫ [T_w (β')² + (K_η + 6 T_Ω) β²] dw` under the concrete Gaussian collapse (`μ_η = T_w = K_η + 6 T_Ω = 1`, `β = exp(-w²/2)`). They are not tautological: `∫ exp(-w²) dw = √π` and `∫(w² + 1) exp(-w²) dw = √π/2 + √π = 3√π/2` are genuine closed-form definite integrals, not algebraic identities that fold to `0` by symbol cancellation. The Mathematica side uses the independent `Integrate[...]` engine path and arrives at the same residual `0`, so engines-agree is preserved. No collateral edits.

Substantive coverage after trim: yes, stage 018 still verifies something — the two Gaussian bridge integrals — and those are precisely the bridge identities the paper card exports. The displaced families now belong to stages 019 and 020; both stage files exist on disk under `scripts/` and `mathematica/` (per the resolution note, stage 019 owns one-pole closure + compatibility; stage 020 owns gate determinant + wall-slope solve + Xi_1). Codex's pre-trim "destination verification" of those owners is presumed; the verifier does not re-audit them here, only confirms the files exist.

One residual: the SymPy docstring (line 2) still reads `"""Master-note audit for step_16_parent_throat_action_bundle_master_notes.md."""` — the dangling notes reference flagged in the v2 directive's "in addition to choosing (a)/(b)/(c)" paragraph was not updated. Flagging as a side observation, not a blocker, because that secondary question was routed to user choice rather than named as a definite required edit, and the trim resolution focused on the algebraic scope.

## Exec log assessment

**SymPy:** exit=0 (per orchestrator note; the canonical exec log `exec_logs/stage_018_sympy.log` is the pre-trim May-21 file, so I use the canonical script output for fresh evidence).

Output file `scripts/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.txt` (mtime 2026-05-25 21:49):
```
STEP 16 PARENT THROAT ACTION BUNDLE MASTER AUDIT
Checked concrete Gaussian parent-wall inertia and stiffness branch integrals.
STATUS: PASS
```

The banner text matches the trimmed `print(...)` (line 29 of the script), confirming the saved output was regenerated against the trimmed source. `STATUS: PASS` implies no `AssertionError` was raised, so both `assert_zero` calls (lines 25, 26) returned `sp.factor(sp.together(sp.simplify(...))) == 0`.

**Mathematica:** exit=0 (per orchestrator; no canonical exec log for math).

Output file `mathematica/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.txt` (mtime 2026-05-25 21:50):
```
STAGE 018 PARENT THROAT ACTION BUNDLE MASTER MATHEMATICA AUDIT
M1 Gaussian inertia integral residual = 0
M1 Gaussian stiffness integral residual = 0
STAGE 018 MATHEMATICA AUDIT PASS
```

Both M1 residuals are `0` and the script reached the terminal `Exit[0]` after `Print["STAGE 018 MATHEMATICA AUDIT PASS"]` (line 29 of the `.wl`). The two `If[FullSimplify[...] =!= 0, ...; Exit[1]]` guards (lines 21–22, 26–27) did not trip.

**Output freshness:**
- sympy script mtime: 2026-05-25 21:41; sympy output mtime: 2026-05-25 21:49 → output is newer (fresh).
- mathematica script mtime: 2026-05-25 21:47; mathematica output mtime: 2026-05-25 21:50 → output is newer (fresh).

Both engines regenerated post-trim. `outputs_fresh: true`.

## Material-change assessment

`material_change`: **true**.

The trim removed five claim families from stage 018's certified scope. Per the resolution note, stages 019 and 020 now own:
- stage 019 (`parent_throat_action_isotropic_bundle`): one-pole numerator identity, two `KSigma` closures (one-pole + normalization), and their compatibility cross-closure.
- stage 020 (`parent_throat_action_weak_axisym_packet`): even-gate 2x2 determinant (`= 1/27`), wall-stiffness slope (`dKSigma = B01 + Z01 + 27(B41 + Z41)`), wall-inertia slope (`dMSigma = -(B21 + Z21) + 3(B41 + Z41)`), and `Xi1` residual amplitude `Xi1 = N01/N0 - 27(B41 + Z41)/(KSigma - B0 - Z0)`.

Downstream concerns:
- Any unit > 018 that previously consumed stage 018's `Xi1`, wall-slope, or KSigma-closure results as a verified upstream must now bind to stage 019 or 020 instead. The orchestrator's blanket `upstream_stale: true` for units > 018 already captures this; the specific narrow re-audits are stages 019 and 020 themselves (to confirm the displaced families landed correctly there) and any stage that imports the parent-action bundle (stages 021 onward — the "Reduced Maxwell/mixed one-port normal form" the auditor explicitly named as a candidate home).
- The Gaussian bridge probes themselves are unchanged (sympy A18/A19 and math integrals match the pre-trim values `√π` and `3√π/2` bit-for-bit), so any consumer that only needed the `M_Σ`/`K_Σ` exports is unaffected.

## Side observations (non-blocking)

1. SymPy docstring (line 2) still references `step_16_parent_throat_action_bundle_master_notes.md`, a file that does not exist anywhere under `/var/projects/toy_physics/research/pde_ledger/`. The auditor's report flagged this and the directive routed the question to the user; the trim resolution did not relabel it. Suggest a follow-up doc-pass to either point at the correct stage-018 notes path or drop the reference. Not blocking because the Gaussian probes are the only substantive content and they don't depend on the docstring.
2. The Mathematica header comment was updated to list only "M1 Gaussian wall inertia and stiffness integrals" — internally consistent with the body, but note that the surviving block uses `M1` as the label for *both* the inertia and the stiffness integral (lines 19, 24), where a small relabel to `M1a`/`M1b` would have matched the auditor's original `M8a`/`M8b` convention. Cosmetic only.
3. Coverage after trim is thin but non-empty: two concrete Gaussian-profile integrals against fixed numeric targets. The auditor's original report already noted this is a `partial` cover of the abstract `M_Σ`/`K_Σ` bridge identity (since `μ_η`, `T_w`, `K_η`, `T_Ω` are collapsed to 1 and only one `β` profile is exercised). The trim does not improve that and no symbolic abstract-bridge check was added — directive option (b) did not require it, but option (c) explicitly would have. If a future audit wants the abstract bridge symbolically exercised, that is a fresh finding for a fresh audit; the verifier does not block on it.

## Verdict justification

The Q4=b TRIM resolution lands cleanly: both engines now verify exactly and only the two Gaussian bridge integrals that the stage-018 paper card declares as `\stagefield{Output}`. The deletions are confined to the two script files, the diff shows no smuggled-in edits, the saved outputs are fresh (post-trim mtimes), both engines exit 0 with non-tautological `residual = 0` results on the surviving checks (`∫ exp(-w²) dw = √π`, `∫(w²+1) exp(-w²) dw = 3√π/2`), and the displaced claim families have named destinations (stages 019 and 020) whose script files exist on disk. The lingering docstring reference to the missing `step_16_*_notes.md` is a non-blocking documentation hygiene item already known to the user. Overall verdict: `verified`, `material_change: true`.
