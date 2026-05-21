---
unit_id: 008
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T12:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 008

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
A new Mathematica audit script was created at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl` (175 lines). The script is Mathematica-idiomatic: it uses `Integrate[..., {w, -Infinity, Infinity}, Assumptions -> ...]`, `Limit[..., R -> Infinity, Assumptions -> ...]`, `FullSimplify`, and `If[... =!= 0, Print["FAIL: ..."]; Exit[1]]` assertions. The script does not mirror the SymPy `assert_zero` helper or section choreography. The stdout transcript was captured at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.txt` (mtime 2026-05-21 11:51, newer than .wl mtime 11:04).

**Assessment:**
All seven claim-manifest items are covered:
- M1 reciprocal identity (line 17, symbolic) — PASS, tautological-by-design (acknowledged in manifest as anchor).
- M2 Pair A Gaussian-Gaussian (lines 28-47) — concretely integrates `gaussianWeight*localizedGaussian` with sigma != lambda, verifies normalization and `xi*I_WZ/I_WH - xi == 0`. Note: M2 Pair A's H=Z residual is structurally `xi*X/X - xi` because both I_WZ and I_WH (with H=Z) use the same integrand `gaussianWeight*localizedGaussian`; the substantive part is the normalization check `gaussNorm - 1` which fails if W is not properly normalized.
- M2 Pair B Lorentzian-Gaussian (lines 49-68) — Lorentzian W against Gaussian Z, non-matched profile. Normalization and H=Z residual both checked.
- M3 matched source for Pair A (lines 70-82) — `mu0*I_WS/I_WZ - mu0/Z_int == 0`. Non-tautological because two distinct integrals (gaussOverlap and gaussMatchedSourceOverlap) are computed and ratio must equal 1/zArea symbolically.
- M4 Gaussian matched-kernel concrete values (lines 84-131) — all seven values verified: Z_int=Sqrt[Pi]*lambda, Z^2_int=Sqrt[Pi/2]*lambda, I_WZ=Sqrt[2]/2, etc.
- M5 reduction-first H=1 regulator limit (lines 133-142) — `Limit[xi*Sqrt[Pi]*lambda/(2R), R->Infinity] == 0`.
- M6 reduction-first H=Z identity (lines 144-150) — symbolic `xi*Z_int/Z_int - xi == 0`.
- M7 Lorentzian numeric checks (lines 152-171) — substitutes (lambda=1, sigma=1/2) and (lambda=1, sigma=2) into the Lorentzian-Gaussian H=Z residual AND the matched-source residual, evaluates to 30 digits, verifies Chop yields {0,0}. This is the substantive guard the auditor requested: a non-matched W (Lorentzian) against a Gaussian Z, with concrete numeric specializations, ruling out FullSimplify coincidence.

Exec log shows STATUS: PASS, exit_code 0. All seven M* PASS lines printed. No collateral edits beyond the new file and the saved transcript.

Note on path: the .wl lives at `mathematica/` not `scripts/`, which is the convention across the repo (stages 001-007, 009-012, etc. all reside in `mathematica/`). The user context confirms the file was briefly misplaced and then moved; the current canonical location matches the repo convention. The orchestrator's exec runner picked it up cleanly (see exec log argv).

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Four assertion labels in `moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py` were relabeled to honestly flag their tautological/cancellation/factoring nature, with explanatory comments pointing readers to the substantive section 5b for real verification:
- Line 73 (was: "effective gauge inverse is reciprocal") → "reciprocal consistency (tautology by construction)" + comment.
- Line 130-133 → "H=Z integral identification (tautology by substitution)".
- Line 137-140 → "zero-mode H=Z gauge-driver projection (factoring of B out of w-integral)".
- Line 141-144 → "H=Z effective gauge via symbolic substitution (tautology after cancellation)".
- Line 156 → "matched source coupling (cancellation tautology)" + comment.

The expression bodies of the assertions are unchanged (per the directive's intent: relabel, do not delete; rely on F3's section 5b for the substantive replacement).

**Assessment:**
The directive explicitly allowed option (b) of F2 — retain the assertions but mark them clearly as consistency tautologies — provided the substantive replacement appears elsewhere. F3 supplies that replacement (section 5b with sigma != lambda). The labels now contain "tautology" / "cancellation" / "factoring", per the directive's verification criterion. No code logic changes, no regression risk. The transcript at section 4 / 3 still shows passing assertions (script exit 0 implied by transcript ending "STATUS: PASS" line 144).

Caveat: by itself F2 does NOT make the assertions non-tautological — but the directive deliberately split the work so F2 honestly labels these as anchors and F3 provides the substantive alternative. Together the pair resolves the original finding.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
A new section 5b was inserted between line 224 (end of section 5) and line 277 (section 6 header) of the SymPy script. The section:
- Declares `sigma = sp.Symbol("sigma", positive=True, finite=True, nonzero=True)`, an independent width symbol distinct from `lam`.
- Defines `W_indep = exp(-w^2/sigma^2) / (sqrt(pi)*sigma)` — a Gaussian observer kernel with width `sigma`, not proportional to Z = exp(-w^2/lam^2).
- Verifies `W_indep` normalization integrates to 1.
- Computes `I_WZ_indep = integrate(W_indep * Z, ...) = lambda/sqrt(lambda^2+sigma^2)` — a substantive integral whose result depends on BOTH widths.
- Asserts `xi_eff_HZ_indep - xi == 0`.
- Computes H=1 case: `xi_eff_H1_indep = lambda*xi/sqrt(lambda^2+sigma^2)`. Evaluates at sigma=1/2, lambda=1 → 2*xi/sqrt(5). Verifies (via `raise AssertionError`) that this is NOT equal to xi at concrete unequal widths — an inversion test that guards against accidental matching.
- Asserts `mu0_eff_source_indep - mu0/Z_int_indep == 0`, where `mu0_eff_source_indep` is built from two distinct integrals (I_WS over W*Z/Z_int and I_WZ over W*Z). This cancellation is now mediated by genuine integration, not by abstract `Integral` placeholders.

**Assessment:**
This is exactly the substantive replacement the auditor requested. The matched-source assertion is genuinely non-tautological: both I_WZ_indep and I_WS_source_indep are computed by separate integrations, and the assertion `mu0*I_WS/I_WZ - mu0/Z_int == 0` requires that the algebraic relation I_WS = I_WZ/Z_int actually hold after concrete evaluation. The H=1 inversion check (the `if ... raise` block) further guards against the failure mode "any W trivially gives xi", by demanding a non-equality at concrete sigma != lambda.

The transcript shows the new section header "5b) Independent-profile check..." (line 99-100 of the .txt) and prints non-trivial expressions: `I_WZ (independent) = lambda/sqrt(lambda**2 + sigma**2)`, `xi_eff_proj(H=1, independent W) = lambda*xi/sqrt(lambda**2 + sigma**2)` — both confirming that the integration was performed and depends on both widths. The matched-source value prints as `mu0/(sqrt(pi)*lambda)`, matching mu0/Z_int_indep.

One small caveat: the assertion `xi_eff_HZ_indep - xi == 0` is structurally trivial because `I_WH_HZ_indep` is computed via the same literal sympy expression as `I_WZ_indep` (both are `sp.integrate(W_indep * Z, (w, -sp.oo, sp.oo))`). So this particular assertion in 5b alone does not catch a hypothetical bug where H=Z gives a different integral than I_WZ. However, the matched-source assertion DOES catch such bugs, and Mathematica M7 independently exercises the non-matched profile numerically with the Lorentzian-Gaussian pair. Across the dual-engine pair, the quantifier "for any normalized W" is now tested with three genuinely distinct W shapes (matched Gaussian, independent Gaussian with sigma != lambda, Lorentzian) and Mathematica's M7 closes the residual concern.

Overall the finding is resolved; the auditor's stated verification criterion ("script defines at least one new profile-pair... computes the three integrals against it, and asserts xi_eff(H=Z) == xi and mu0_eff(S=Z/Z_int) == mu0/Z_int non-tautologically. Output transcript shows the new section with the new symbol named.") is met.

## Exec log assessment

**SymPy:** exit=n/a. The orchestrator did NOT capture a stage_008_sympy.log under `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/`. Only `stage_008_mathematica.log` and `stage_008_diff.patch` are present. However, the saved sympy transcript at `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.txt` (mtime 2026-05-21 11:26, newer than script mtime 11:05) ends with `STATUS: PASS` (line 144), shows all sections 1–7 including the new 5b, and the codex_logs note "exits 0 with STATUS: PASS". I treat the script as passing based on the fresh transcript and codex's reported exit; flag the missing engine-log as a side observation. Notable lines from the saved transcript:

- Line 102: `normalization ∫W dw = 1` (Pair A normalization confirmed).
- Line 103: `I_WZ (independent) = lambda/sqrt(lambda**2 + sigma**2)` (substantive integral, depends on both widths).
- Line 104: `xi_eff_proj(H=Z, independent W) = xi` (H=Z preserves xi).
- Line 106: `mu0_eff_proj(S=Z/Z_int, independent W) = mu0/(sqrt(pi)*lambda)` (matched-source cancellation works on the independent profile).
- Line 144: `STATUS: PASS`.

**Mathematica:** exit=0. Notable lines from `stage_008_mathematica.log`:

- Line 7: `M1 reciprocal residual = 0` / Line 8: `PASS: M1 reciprocal identity`.
- Line 11: `PASS: M2 Pair A Gaussian-Gaussian` (with normalization and H=Z residual both 0).
- Line 14: `PASS: M2 Pair B Lorentzian-Gaussian`.
- Line 18: `M4 residuals = {0, 0, 0, 0, 0, 0, 0}` followed by `PASS: M4 Gaussian matched-kernel concrete values`.
- Line 23-25: `M7 Lorentzian numeric H=Z residuals = {0, 0}` and `M7 Lorentzian numeric matched-source residuals = {0, 0}`, followed by `PASS: M7 Lorentzian numerical independent-profile checks`.
- Line 26: `STATUS: PASS`, Line 27: `# exit_code: 0`.

**Output freshness:** Confirmed.
- SymPy script mtime 2026-05-21 11:05 < SymPy output mtime 11:26. Fresh.
- Mathematica .wl mtime 2026-05-21 11:04 < Mathematica output mtime 11:51 (and a duplicate copy in scripts/output/ at 11:06). Fresh.
- Diff patch captured at 11:53.

## Material-change assessment

`material_change`: false.

No derived numerical or symbolic result that downstream units depend on was altered. The substantive claims (xi_eff_proj(H=Z) = xi for normalized W, mu0_eff_proj(S=Z/Z_int) = mu0/Z_int, the matched-Gaussian concrete values Sqrt[Pi]*lambda etc., the H=1 regulator limit → 0) all remain identical. The edits to F1 added independent verification of these same claims in a second CAS; F2 only relabeled assertion strings (no numeric change); F3 added a new section that confirms the same claims against an additional profile (the new printed expressions for I_WZ_indep and xi_eff_proj(H=1, independent W) are new but they document a wider exercise, not a corrected result).

Downstream units (009 onward) that import the projected-Maxwell extension claims see no changed numbers. The unit's role in the chain is unchanged.

## Side observations (non-blocking)

1. **Missing sympy exec log:** `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_008_sympy.log` is absent. The Mathematica log is present, the diff patch is present, and the sympy stdout was clearly captured to its scripts/output/ destination (mtime later than script mtime). This is an orchestrator capture gap, not a code defect. Flag for the orchestrator: please capture stage_*_sympy.log alongside the mathematica log for completeness.

2. **Duplicate Mathematica output:** Both `mathematica/output/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.txt` (mtime 11:51, the canonical one) and `scripts/output/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.txt` (mtime 11:06, older) exist. The older one in scripts/output/ is a leftover from codex's initial placement before the script was moved to mathematica/. Recommend the orchestrator clean up the stale scripts/output/ Mathematica .txt to avoid confusion later. Not a verification blocker.

3. **F3's H=Z assertion in 5b is structurally definitional:** As noted in F3 assessment, `xi*I_WZ_indep/I_WH_HZ_indep - xi == 0` uses the same literal sympy `integrate(...)` expression for both I_WZ and I_WH(H=Z), so this particular assertion in 5b alone would not catch a bug in how H=Z is identified. The substantive non-tautological work in 5b is the matched-source assertion and the H=1 inequality test; the Mathematica M7 numeric run with the Lorentzian-Gaussian non-matched pair fully closes the quantifier check. Not a blocker — pointing it out only so a future auditor knows where the load is actually borne.

4. **Section 6 still uses I_WH_HZ_indep via same integral as I_WZ_indep:** Cosmetic — could in principle be merged in a refactor, but the directive said not to touch sections 6-7 beyond what F2 specified, which codex correctly respected.

## Verdict justification

All three findings have been addressed substantively. F1 produced a 175-line independent Mathematica audit covering all seven claim-manifest items, with `If[... Exit[1]]` style assertions; the exec log shows STATUS: PASS, exit 0. F2 honestly relabels the four tautological/cancellation/factoring assertions while preserving the assertions themselves as readability anchors (option (b) of the directive). F3 inserts a new section 5b with an independent Gaussian observer kernel of distinct width `sigma`, verifies normalization, computes substantive integrals that depend on both widths, and asserts both the H=Z claim and the matched-source claim non-tautologically against the independent profile; additionally, Mathematica M7 exercises the quantifier numerically with a Lorentzian-Gaussian non-matched pair. No regressions visible in the diff. The missing sympy exec log is an orchestrator capture issue, not a script issue — the sympy transcript itself is fresh and shows STATUS: PASS. Verdict: `verified`.
