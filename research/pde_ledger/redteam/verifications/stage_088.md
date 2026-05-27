---
unit_id: 088
batch: III.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
verdict_history:
  - 2026-05-27 initial verifier run: needs_rework (F1 mma comment-bug regression — `stage085_*)` substring prematurely closed comment, silently skipped Pi_tr_from_rho assertion)
  - 2026-05-27 orchestrator hot-fix: reworded comment line 87 to "in the stage 085 Mathematica audit files" (strips `*)` substring)
  - 2026-05-27 re-run mma → rc=0 with all 9 PASS lines including Pi_tr_from_rho and regime trichotomy → verdict promoted to verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 088

## Per-finding outcomes

### F1 — tautological_check (Pi_tr identity must invoke rho_alpha)

**Classification:** resolved (after orchestrator hot-fix; initial run was `partial` due to a comment-syntax bug — see verdict_history above)

**What changed:**
- `scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:111-124` — the literal `Pi_tr/C_mix - 4/3` tautology has been replaced with `Pi_tr_from_rho = rho_min * C_mix` followed by `expect_zero("Pi_tr_from_rho - (4/3) C_mix", Pi_tr_from_rho - sp.Rational(4,3)*C_mix)`, and the regime check is now exercised as Python `assert 1 < rho_min < 2`. Matches the directive's "After" block. The comment cites "stage-085" rather than "Stage-68" (the notes call it Stage 68; the script comment names the upstream verifying file). Harmless deviation.
- `mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:86-91` — the corresponding `piFromRho = FullSimplify[rhoMin*cMix]`, `expectZero["Pi_tr_from_rho - (4/3) C_mix", piFromRho - (4/3)*cMix]`, and `expectTrue["1 < rho_min < 2 ..."]` lines are present in the source. **However, lines 86-87 contain a Mathematica comment whose body includes the substring `stage085_*)`. The first `*)` inside the comment text closes the comment prematurely, and the trailing ` Substitute rho_min. *)` is then parsed as code.** The exec log captures `Syntax::sntx: Invalid syntax in or before "in mathematica/moving_throat_pde_stage085_*). Substitute rho_min. *)"` (output line 27). After the syntax error, Mathematica's parser recovers to `Exit[0]` and the run exits cleanly (rc=0) — but the F1 assertion `expectZero["Pi_tr_from_rho - (4/3) C_mix", ...]` and the `expectTrue["1 < rho_min < 2 (symmetric-lowest-twin regime)"]` line are silently skipped. The output log contains 7 `PASS:` lines, not the expected 9, and no PASS line for "Pi_tr_from_rho - (4/3) C_mix" or the regime trichotomy.

**Assessment:**
SymPy side is correct and non-tautological: substituting `c_contact = 1/2` upstream would now flip the residual to `(2/3)*C_mix` and break the check, as the directive's mutation argument shows.

Mathematica side is broken in a subtle way: the F1 fix code is *physically present* in the file at lines 88-91, but it never executes because the comment immediately above it accidentally closes early on `stage085_*)`. Because Mathematica downgrades the parse failure to a message (not a hard abort) and the following `Exit[0]` at line 96 still runs, the orchestrator-reported `rc=0` is misleading — the F1 Mathematica assertion is silently inactive. This is the exact failure mode the verifier is supposed to catch (a script that "passes" because the check never ran).

### F2 — tautological_check (contact-plus-pole derived rather than restated)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:46-90` — `c0_from_rho`, `c1_from_rho` definitions are gone. In their place, `Y_rho_u = sp.simplify(Y_rho.subs(omega**2, u * Omega_Q**2))` is constructed (the hot-fix path: substitute `omega**2 -> u * Omega_Q**2` then `sp.simplify`, instead of the literal-form-of-the-ratio subs which sp.simplify reshapes into `(Omega_Q**2 - omega**2)` form), then `c1_extracted = sp.simplify(sp.limit((1 - u) * Y_rho_u, u, 1))` followed by `c0_extracted = sp.simplify(Y_rho_u.subs(u, 0) - c1_extracted)`. The round-trip checks at L85-90 now substitute the *extracted* coefficients (rather than the definitions) into `rho_from_c0`, `rho_from_c1`, `zeta_from_c`.
- Same file, L92-109 — the minimal-isotropic-module block now extracts (c0_paper, c1_paper) from `Y_Q_paper = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` via the same pole-residue-then-subtract pattern, asserts `c0_paper - 3/4 == 0` and `c1_paper - 1/4 == 0` (the substantive coefficient-matching check the directive demanded), then derives `rho_min` and `zeta_min` from those extracted constants.
- `mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:45-65` — Mathematica engine performs the analogous extraction directly on `yQpaper = 3/4 + (1/4)/(1 - omega^2/omegaQ^2)` using `Limit[(1 - u)*yQpaperU, u -> 1]` for c1 and `(yQpaperU /. u -> 0) - c1Paper` for c0. The `c0_paper - 3/4`, `c1_paper - 1/4`, `c0_paper + c1_paper - 1` checks all PASS in the exec log (lines 10-15).

**Assessment:**
The hot-fix the user described is in place at lines 55-56 of the sympy file: `Y_rho_u = sp.simplify(Y_rho.subs(omega**2, u * Omega_Q**2))` (substitute `omega**2 -> u * Omega_Q**2` first, then `sp.simplify`). The earlier formulation `Y_rho.subs(omega**2/Omega_Q**2, u)` would fail because `Y_rho` after `sp.simplify` becomes `(Omega_Q**2*rho_alpha - omega**2)/(rho_alpha*(Omega_Q**2 - omega**2))` (visible in sympy output line 6), in which the substring `omega**2/Omega_Q**2` no longer occurs as a syntactic subexpression. The substitute-`omega**2`-then-resimplify pattern avoids this trap.

The extraction is now non-tautological. The mutation test: if the paper Input `Y_Q_paper` were perturbed to `1/2 + (1/2)/(1 - omega^2/Omega_Q^2)`, then `c1_paper` would extract to `1/2` (not `1/4`), the assertion `c1_paper - 1/4 == 0` would fail, and `rho_min` would shift to `2`. Same for the c0 extraction. The check now exercises the extraction machinery rather than restating the parameterization.

The sympy log confirms `c0 extracted from Y_rho (subtract pole at u=0) = 1/rho_alpha` and `c1 extracted from Y_rho (pole residue at u=1) = (rho_alpha - 1)/rho_alpha` (lines 8-9), then `c0_paper - 3/4 = 0` and `c1_paper - 1/4 = 0` (lines 20-21) — the symbolic extractions yield the closed-form parameterization, and the paper-form extractions yield the rational constants.

### F3 — mathematica_transliteration (independent derivation path)

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:34-91` — the Mathematica audit body is restructured around the paper-form `yQpaper = 3/4 + (1/4)/(1 - omega^2/omegaQ^2)` Input rather than the rho-parameterized `Y_rho`. Coefficient extraction uses `Limit` and subtraction in the variable `u = omega^2/omegaQ^2`. The script then *reconstructs* `yRhoFromCoeffs = c0Paper + c1Paper/(1 - omega^2/omegaQ^2)` and asserts `yQpaper - yRhoFromCoeffs == 0` (line 77-78), and *also* verifies that the rho-parameterized form rebuilds yQpaper after substituting `rhoAlpha -> rhoMin` (line 82-84).
- The intermediate symbol sequence (`yQpaper`, `yQpaperU`, `c1Paper`, `c0Paper`, `rhoMin`, `zetaMin`, `yRhoFromCoeffs`, `yRhoParam`) and the assertion order no longer line up with the SymPy script's (`Y_loading`, `Y_rho`, then extract from `Y_rho` to get symbolic c0/c1 in rho_alpha, *then* extract from `Y_Q_paper` to get 3/4 and 1/4). The two engines take genuinely different algebraic paths — sympy works the symbolic rho_alpha form first and then specializes; mathematica works the paper-form numeric Input directly.
- Banner: `mathematica/...mathematica_audit.wl:32` now reads `"STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE"` (the stale STAGE 71/071 label is gone). SymPy banner at `scripts/...sympy_audit.py:29` matches.

**Assessment:**
The two engines converge on the same numeric answers (`rhoMin = 4/3`, `zetaMin = 1/3`, both extracted from independently-derived coefficients) through structurally different machinery. F3's "the second engine must derive the claim independently" criterion is met: the mma script no longer mirrors the sympy script line-for-line.

The exec log (mma) lines 7-25 trace through the new algebraic path: paper form → extract c0/c1 via subtract/limit → invert to get rho_min and zeta_min → reconstruct yQpaper from extracted coeffs → check rho-parameterized form rebuilds yQpaper. All seven of these checks PASS.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `c0 extracted from Y_rho (subtract pole at u=0) = 1/rho_alpha` (L8 — extraction returns symbolic closed form, confirming the rho-parameterized identity is recovered rather than restated)
- `c1 extracted from Y_rho (pole residue at u=1)  = (rho_alpha - 1)/rho_alpha` (L9)
- `c0_paper - 3/4 = 0` (L20) and `c1_paper - 1/4 = 0` (L21) — F2 substantive checks pass
- `Pi_tr_from_rho - (4/3) C_mix = 0` (L26) — F1 substantive check passes on sympy side

All 13 `expect_zero` calls produce `= 0` output lines; no `AssertionError` is raised.

**Mathematica:** exit=0 (reported), but the run is **partially aborted by a parse-time syntax error**. Notable lines:
- `c0 (subtract-pole) = 3/4` (L8), `c1 (pole residue)  = 1/4` (L9) — F3 independent extraction works
- `PASS: c0_paper - 3/4`, `PASS: c1_paper - 1/4`, `PASS: c0_paper + c1_paper - 1` (L11-15) — F2 mma side passes
- `PASS: rho_min - 4/3`, `PASS: zeta_min - 1/3` (L19-21) — numeric anchors pass
- `PASS: paper form - reconstruction from extracted (c0, c1)` (L23) — F3's substantive coefficient-matching check passes
- `PASS: rho-parameterized form (rhoAlpha -> rhoMin) - paper form` (L25) — rho parameterization cross-check passes
- **`Syntax::sntx: Invalid syntax in or before "in mathematica/moving_throat_pde_stage085_*). Substitute rho_min. *)"` (L27)** — the F1 Mathematica fix's *enclosing comment* contains the substring `stage085_*)` which closes the outer comment prematurely. The trailing fragment ` Substitute rho_min. *)` is then parsed as code, fails, and Mathematica skips ahead.
- No `PASS: Pi_tr_from_rho - (4/3) C_mix` line.
- No `PASS: 1 < rho_min < 2 (symmetric-lowest-twin regime)` line.
- No `Stage 088 Mathematica audit passed.` line at the end.

Only 7 `PASS:` lines appear in the log, but the file declares 8 `expectZero`/`expectTrue` calls. The F1 Mathematica assertion and the regime check are *silently inactive*.

**Output freshness:**
- sympy script mtime 2026-05-27 10:25:27, output mtime 2026-05-27 10:25:41 — output is 14 s newer ✓
- mma script mtime 2026-05-27 10:20:03, output mtime 2026-05-27 10:26:25 — output is ~6 min newer ✓

Both outputs are post-fix.

## Material-change assessment

`material_change`: **false**.

No downstream numeric or symbolic result has changed: `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pi_tr = (4/3) C_mix`, and the regime classification `1 < rho_min < 2` are all identical to the pre-fix values. The edits restructure the *path* to those values (introducing pole-residue extraction and inverting `1/c0`), but the values themselves are unchanged. Downstream stages 089+ that consume Stage 088's outputs see no change.

## Side observations (non-blocking)

- The SymPy script's comment at lines 111-113 cites the upstream identity as coming from "stage-085 product identity" rather than the notes' "Stage 68 identity." Both Stage 068 and Stage 085 establish the `Pi_tr / C_mix = alpha_req / alpha_mix` chain at different layers (the notes name Stage 68 for the symbolic identity; Stage 085 is where the product identity is independently scripted). Not a contradiction, but the directive's language used "Stage 68." Not a blocker.
- The Mathematica `expectTrue["1 < rho_min < 2 (symmetric-lowest-twin regime)", 1 < rhoMin < 2]` form is fine in principle (`1 < 4/3 < 2` evaluates to `True` symbolically), but the regime check is currently dead code due to the comment bug; once the comment is fixed it should produce `PASS: 1 < rho_min < 2 ...`.
- The `Limit::alimv` warning at output line 6 ("Assumptions that involve the limit variable are ignored") is benign — `$Assumptions` includes `-1 < u < 1`, which Mathematica notes can't be used inside the limit but does not affect the result (the residue `(rho_alpha-1)/rho_alpha` is recovered exactly). Worth keeping an eye on but not a finding.

## Delta directive (applied as orchestrator hot-fix, retained for history)

The initial verifier flagged this one regression; it was applied orchestrator-direct (one-line comment reword), the Mathematica audit was re-run, and all 9 expected PASS lines now appear. The original delta directive is retained below for the post-batch record.

### F1-mma-comment (regression — Mathematica F1 assertion is silently skipped)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:86-87`

**Issue:**
The comment at lines 86-87 contains the substring `stage085_*)` whose embedded `*)` closes the outer Mathematica comment prematurely. The fragment ` Substitute rho_min. *)` after the spurious close is then parsed as code, producing `Syntax::sntx`. Mathematica recovers and reaches `Exit[0]` at line 96, masking the failure as rc=0, but the F1 Mathematica assertion `expectZero["Pi_tr_from_rho - (4/3) C_mix", piFromRho - (4/3)*cMix]` at line 90 and the regime check at line 91 never execute. Only 7 of 8 expected `PASS:` lines appear in the exec log.

**Required change:**
Replace lines 86-87:
```
(* Stage-085 product identity: Pi_tr = rho_alpha * C_mix (verified upstream
   in mathematica/moving_throat_pde_stage085_*). Substitute rho_min. *)
```
with a version that does not embed `*)` inside the comment body:
```
(* Stage-085 product identity: Pi_tr = rho_alpha * C_mix (verified upstream
   in mathematica/moving_throat_pde_stage085 files). Substitute rho_min. *)
```
(or any rewording that strips the `*)` substring — e.g. drop the asterisk: `stage085 mathematica scripts`).

**Verification:**
After the edit, re-running the Mathematica audit must produce two additional `PASS:` lines:
- `PASS: Pi_tr_from_rho - (4/3) C_mix`
- `PASS: 1 < rho_min < 2 (symmetric-lowest-twin regime)`
and the final `Stage 088 Mathematica audit passed.` line. Total PASS lines must be 9. The `Syntax::sntx` message must no longer appear.

No other edits are needed; F2 and F3 are already correctly applied.

## Verdict justification

F2 (contact-plus-pole derived rather than restated) and F3 (Mathematica independent derivation path + banner relabel) are correctly applied in both engines. The user's noted SymPy hot-fix (`omega**2 -> u * Omega_Q**2` then `sp.simplify`, rather than substituting on the combined ratio) is in place at lines 55-56 of the sympy file and the extraction is non-tautological. F1 is correctly applied on the SymPy side (`Pi_tr_from_rho = rho_min * C_mix`, sympy log line 26 confirms PASS). **F1 on the Mathematica side is silently dead code**: the enclosing comment at lines 86-87 contains an embedded `*)` substring (`stage085_*)`) that prematurely closes the comment, causing `Syntax::sntx` and skipping past the `Pi_tr_from_rho` assertion and the regime trichotomy check. Mathematica's parser recovery reaches `Exit[0]`, so rc=0 is reported, but only 7 of 8 expected PASS lines appear and the substantive F1 Mathematica check never runs. This is the canonical "passing exec log is necessary but not sufficient" failure mode the verifier prompt warns against. One narrow delta fix (reword the comment to remove the embedded `*)`) closes the gap; no other rework needed.
