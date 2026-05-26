---
unit_id: 014
batch: I.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T22:15:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: true
---

# Verification — unit 014

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim)

**Classification:** resolved

**What changed:**

User selected resolution direction (b) "Scripts over-cover; trim them." Codex (per directive follow-up) trimmed the Stage 014 scripts to retain only the even-gate mechanism sieve (K_1, H_even, compensation denominator). The `Xi_load`, `\delta P_2`, `\delta P_4` content and the printed-only `Compat` were removed from both engines, with downstream ownership reassigned: stage 010 owns `\delta P_2`/`\delta P_4` and stage 013 owns `Xi_load`.

Concrete edits visible in `exec_logs/stage_014_diff.patch`:

- `scripts/moving_throat_pde_stage014_..._sympy_audit.py`:
  - Module docstring (lines 1-15) rewritten: primitive packet trimmed from `{Q, S2, H_port, Delta, P, G_W}` to `{Q, S2, H_port, Delta}`; bottleneck list trimmed from `{Xi_load, K1, H_even}` plus `delta P2, delta P4` to just `{K1, H_even}`.
  - Symbol declarations trimmed (current lines 34-39): dropped `P, Gw, D2, D4, N0, Ptarget, p1, g1`. `D0` retained for the construction-level independence loop.
  - `n0, n2, n4` primitive definitions removed.
  - `subs_der` truncated to `{q1: mu1*Qx, s1: mu1*Sx, h1: mu1*Hx, d1: mu1*Dx}`.
  - `n0d, n2d, n4d` lifted-slip definitions removed; `Xi = ...` (old line 67), `Compat = ...` (old line 72), `deltaP2`, `deltaP4` (old lines 74-75), and the three `coef_*` derivative assignments (old lines 93-95) removed.
  - `assert_zero` calls for "d Xi_load / d Pprime", "d deltaP2 / d G_W prime", and "deltaP4 should depend on G_W prime" removed (old lines 126-128).
  - Linearity loop (current lines 80-82) iterates over `(z2d, z4d)` only, retaining the `(Px, Gx)` slips in the slip set so the construction-level independence check on K1/He still exercises `P'`/`G_W'` orthogonality.
  - Independence loop for `Sx, Hx, Gx` on `Xi` removed; the K1/H_even loop over `(Px, Gx)` retained.
  - Bundle pull-back `Xi_bundle` assertion removed; `K1_bundle` and `He_bundle` assertions retained (current lines 94-97).
  - Print sections renumbered: removed Section 4 "Isotropic compatibility and constant-prefactor transport" and the multi-channel readout language; current readout reads "repairs the conservative even gates K1 and H_even through the primitive slopes (Q', Delta', S2', H_port') when the compensation denominator is nonzero."
- `mathematica/moving_throat_pde_stage014_..._audit.wl`:
  - `Clear[...]` and `$Assumptions` lists trimmed to drop `P, Gw, D2, D4, Nbase`.
  - `N0[...], N2[...], N4[...]` primitive defs removed.
  - `n0d, n2d, n4d` lifted-slip defs removed.
  - `XiLoad`, `deltaP2`, `deltaP4` definitions removed.
  - `M1 XiLoad independence` loop, `M2 dXiLoad/dPx`, `M3 d(deltaP2)/dGx`, `M4 deltaP4 Gx dependence` assertions all removed.
  - Surviving guards: M1 K1/He independence (Px, Gx); M5 source/denom solve; M6 spectral solve; M7/M8 Jacobian determinants; M9 Hx and Sx compensation denominator factorizations; M10 sign-flip mutation.

**Assessment:**

The trim is exactly what direction (b) called for: every audit-flagged "extra" row in the original paper↔script cross-check table is gone. The surviving sympy assertion inventory (A1, A3, A5, A6, A10, A11, A12, A13, A15 from the original report) maps one-to-one onto the paper card's three deliverables (K_1 definition, H_even definition, compensation denominator proportional to `Q S_2 - \Delta H_{port}`) plus the mechanism-sieve negative-existence claim. No paper-side change was made, consistent with the no-prose-edit rule. The trim is downstream-relevant (see material-change section), but the per-finding outcome is "resolved": script scope now matches paper scope.

Minor wrinkle: the sympy file still declares `Px, Gx` as symbols (line 42) and exercises them in the linearity-and-independence loops (lines 80, 89-91). This is correct and intentional — those loops verify that `P'` and `G_W'` are absent from K1/H_even (which is a paper-relevant negative claim per the Section-3 print "P' and G_W' do not enter K1 or H_even"), so the symbols must exist to be substituted. No leftover load-bearing transport-ledger content.

### F2 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage014_..._sympy_audit.py`: the five-line block that was at the previous lines 139-143 (three-line comment plus two tautological `assert_zero("compensation K1", K1.subs(comp_surface))` / `assert_zero("compensation H_even", He.subs(comp_surface))` calls) is gone. Visible in the diff as the removal of the `# Note: ...` comment and the two `assert_zero("compensation K1/H_even", ...)` lines. Current file line 107 reads `raise AssertionError(f"Unexpected pure spectral solve: {sh_only}")` and is followed directly by line 108 `Z2_slot = (Q*S2 - Hport*Delta)/Delta**2`, exactly as the directive specified ("After deletion, line 138 should be followed directly by what is currently line 144"). Grep confirms no `assert_zero("compensation K1"` or `assert_zero("compensation H_even"` strings remain.

**Assessment:**

Matches the directive verbatim. The substantive compensation-surface verification is preserved by the surviving lines 100-102 (denominator factorization plus sign-flip mutation) and lines 108-110 (Z2-slot identification). `comp_surface` is still defined (line 63) and consumed by lines 73-74 (denominator extraction) and lines 169, 171 (print readout); removing the tautological pair has no side effects. `assert_zero` is silent on success so no print-transcript line was emitted by those calls. No collateral edits.

## Exec log assessment

**SymPy:** exit=0 (inferred — orchestrator-captured `stage_014_sympy.log` is stale from 21 May 14:58 and predates the fix; per the verification prompt the canonical artifact is the saved output transcript). Notable lines from `scripts/output/moving_throat_pde_stage014_..._sympy_audit.txt`:

- Lines 1-10: Section 1 "Mouth-local Taylor ansatz" reaches readout with the trimmed primitive list `{q1: mu1*Q'(0), s1: mu1*S2'(0), h1: mu1*H_port'(0), d1: mu1*Delta'(0)}` — no `p1, g1`.
- Lines 30-32: dependency pattern print confirms "K1 and H_even depend only on (Q', Delta', S2', H_port')" and "P' and G_W' do not enter K1 or H_even", which is the negative-claim flip side the surviving Px/Gx independence loops verify.
- Lines 40, 43: mechanism sieve `[{Deltax: 0, Qx: 0}]` and `[{Hx: 0, S2x: 0}]` — both trivial-only, anchoring the negative-existence claim.
- Lines 53, 55: comp_surface denominators visibly factor as `9*Delta**2*(Delta*H - Q*S2)` and `9*Delta*(Delta*H - Q*S2)`.
- Line 70: `STATUS: PASS`.

The script reaches the final `STATUS: PASS` print, so all surviving `assert_zero` / `assert_nonzero` / `raise AssertionError` checks held. The previously tautological "compensation K1/H_even" prints were never emitted (assert_zero is silent on success), so the transcript line count for those guards is unchanged — they simply have no representation in the output, which is the documented behavior.

**Mathematica:** exit=0 (inferred from saved transcript — orchestrator-captured `stage_014_mathematica.log` is not present in `redteam/exec_logs/`). Per `mathematica/output/moving_throat_pde_stage014_..._audit.txt`:

- Lines 2-9: M1 K1/He independence from Px and Gx — all four sub-checks residual 0, PASS.
- Lines 10-13: M5 source/denominator solve `{{Qx -> 0, Deltax -> 0}}`, M6 spectral solve `{{Sx -> 0, Hx -> 0}}` — both residual `True`, PASS.
- Lines 14-17: M7 source/denom Jacobian determinant nonzero (residual is a non-zero rational expression), M8 spectral Jacobian determinant residual `(-(Delta*Hport) + Q*S2)/Delta^4` — both PASS.
- Lines 18-21: M9 Hx and Sx compensation denominators residual 0, PASS.
- Lines 22-23: M10 sign-flip mutation residual `-18*Delta^2*Q*S2`, PASS.
- Line 24: `STATUS: PASS`.

The surviving Mathematica guard inventory M1/M5/M6/M7/M8/M9/M10 is exactly the mirror of the surviving SymPy inventory; M2/M3/M4 (XiLoad/deltaP2/deltaP4 transport-ledger guards) are gone, consistent with F1's trim.

**Output freshness:**
- SymPy script: mtime 2026-05-25 21:55. SymPy output: mtime 2026-05-25 21:56. Output is newer than script. Fresh.
- Mathematica script: mtime 2026-05-25 21:39. Mathematica output: mtime 2026-05-25 21:50. Output is newer than script. Fresh.

## Material-change assessment

`material_change`: true.

The trim of `Xi_load`, `\delta P_2`, `\delta P_4` from Stage 014 is a scope reassignment, not a derivation change — the actual coefficients (`d Xi_load / dP' = 2/P`, `d \delta P_2 / dG_W' = -2P/(D_0 \Delta^2)`) are unchanged where they live (downstream). However, *Stage 014 no longer verifies them*. Per the user note ("stage 010 owns δP_2/δP_4; stage 013 owns Xi_load"), the verification responsibility has moved. The orchestrator should confirm that stages 010 and 013 now do actually carry the relevant assertions and that no other stage was relying on Stage 014's transcript prints of `delta P2`, `delta P4`, or `d(Xi_load)/dP'` for cross-check anchoring. Stages numerically > 014 are unaffected by the *derived results* (K1, H_even, compensation surface) which are unchanged, but any stage that grepped Stage 014's output for transport-ledger numbers should be re-verified.

## Side observations (non-blocking)

1. The orchestrator-captured `stage_014_sympy.log` (mtime 21 May) and missing `stage_014_mathematica.log` are pre-fix; only the canonical saved-output transcripts under `scripts/output/` and `mathematica/output/` reflect the post-fix state. This is acknowledged in the verifier prompt and is not a defect.
2. The previous verification file at `redteam/verifications/stage_014.md` (verdict `verified`, 4-of-4 findings) corresponded to an earlier audit pass with a different finding set (F1-F4 all tautological / missing-mathematica). This new write overwrites that with the current 2-finding audit's outcome.
3. The trimmed sympy script still imports/declares `Px, Gx, D0` and the linearity loop iterates over `(Qx, Sx, Hx, Dx, Px, Gx) × (z2d, z4d)`. Since z2d and z4d are constructed from the four-primitive `subs_der` (no `p1, g1`), the `Px`/`Gx` second-derivative-zero checks are trivially zero by construction (z2d and z4d do not contain Px or Gx). These are weak checks but they are also non-load-bearing (the substantive linearity content is the four-slip cases); auditor will likely either accept or downgrade in a future pass. Not a blocker.
4. Verifier did not introduce new findings, per protocol.

## Verdict justification

Both findings are resolved with edits matching the user-chosen resolution direction (F1=trim) and the directive's verbatim deletion specification (F2). The surviving sympy and mathematica scripts both exit 0, produce fresh transcripts, and exercise exactly the three paper deliverables (K_1, H_even, compensation denominator proportional to `Q S_2 - \Delta H_{port}`) plus the mechanism-sieve negative-existence claim and the sign-anchor mutation. No collateral edits, no derived-result regressions. Material change is flagged true because Stage 014's *scope* contracted; downstream ownership for `Xi_load`/`\delta P_2`/`\delta P_4` lies with stages 013 and 010 per the user resolution and should be confirmed by the orchestrator's next audit sweep over those units. Verdict: `verified`.
