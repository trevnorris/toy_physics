---
unit_id: 090
batch: III.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 090

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl:34-41` replaces the hardcoded `rhoAlpha = 4/3; zetaReq = rhoAlpha - 1;` with `cContact = 3/4; cPole = 1/4; rhoAlpha = 1/cContact; zetaReq = cPole/cContact;`. The Print lines at .wl:46-49 now emit `c_contact`, `c_pole`, `rho_alpha`, `zeta_req`. The single tautology at the old line 46 has been replaced by three substantive `expectZero` calls at .wl:54-56: `rho_alpha - 4/3`, `zeta_req - 1/3`, and the identity `zeta_req - (rho_alpha - 1)`.

**Assessment:**
The edit matches the directive verbatim (including the docstring rewording). The identity `zeta_req - (rho_alpha - 1)` is no longer a definitional tautology because `rhoAlpha = 1/cContact = 4/3` and `zetaReq = cPole/cContact = 1/3` are derived independently from `(cContact, cPole) = (3/4, 1/4)`; the identity now depends on the substantive relation `cContact = cPole + 1/2` from the upstream module. The Mathematica exec log shows `PASS: rho_alpha - 4/3`, `PASS: zeta_req - 1/3`, and `PASS: zeta_req - (rho_alpha - 1)` with residual `0`. No collateral edits beyond what the directive asked for.

### F2 — script_missing_paper_claim

**Classification:** resolved

**What changed:**
Same edit as F1: the new `expectZero["rho_alpha - 4/3", ...]` and `expectZero["zeta_req - 1/3", ...]` at .wl:54-55 add the previously-missing paper-side coverage to the Mathematica engine. SymPy already covered these at .py:63-64.

**Assessment:**
Mathematica output now contains the two PASS lines for `rho_alpha - 4/3 = 0` and `zeta_req - 1/3 = 0`. Both engines now substantively certify the locked triple's first two components. The directive's note that F2 is resolved by F1's edit is correct.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py:101-105` adds the carry-forward anchor comment plus `Pe_req = sp.Integer(0)` and `expect_zero("Pe_req (carry-forward from Stage 075 transport map)", Pe_req)`. The matching edit appears at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl:61-65`: comment block plus `peReq = 0; expectZero["Pe_req (carry-forward from Stage 075 transport map)", peReq];`.

**Assessment:**
Both engines now define an explicit `Pe_req`/`peReq` symbol set to `0` with an inline comment linking the proxy inequality `zeta_req < A_F1` to the locked triple value, anchored to the upstream transport-map stage. Exec logs show the new PASS lines in both engines (`Pe_req (carry-forward from Stage 075 transport map) = 0`). Deviation from directive: the directive text says "Stage 062 transport map" but Codex wrote "Stage 075 transport map" in both engines; the sympy docstring at .py:19 still reads "Stage 62 transport map". This is a label-only difference (the stage number for the transport map carry-forward); since the user's hand-off note explicitly says "stage-075 transport map" this is an intentional renumber and not a regression. The substantive assertion (existence of `Pe_req = 0` symbol with anchored comment) is what the finding required.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `rho_alpha - 4/3 = 0`
- `zeta_req - 1/3 = 0`
- `zeta_req - (rho_alpha - 1) = 0`
- `Pe_req (carry-forward from Stage 075 transport map) = 0`

All three inequalities also evaluate `True`. The final ledger banner still concludes with `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pe_req = 0`.

**Mathematica:** exit=0. Notable lines:
- `PASS: rho_alpha - 4/3`
- `PASS: zeta_req - 1/3`
- `PASS: zeta_req - (rho_alpha - 1)`
- `PASS: Pe_req (carry-forward from Stage 075 transport map)`
- All three `expectTrue` inequality checks also PASS.

**Output freshness:** sympy `.py` mtime `2026-05-27 10:22:25`, output mtime `2026-05-27 10:24:55`; mathematica `.wl` mtime `2026-05-27 10:22:21`, output mtime `2026-05-27 10:26:37`. Both outputs post-date their scripts. Fresh.

## Material-change assessment

`material_change`: false.

No derived numeric values changed. `rho_alpha`, `zeta_req`, the three carried thresholds, and the resulting inequality verdicts are identical to the pre-fix state. The new `Pe_req = 0` symbol formalizes a carry-forward that was already implicit in the proxy inequality. Mathematica's `rho_alpha` is now derived rather than hardcoded but evaluates to the same `4/3`. Downstream units depending on Stage 090's locked triple `(rho_alpha = 4/3, zeta_req = 1/3, Pe_req = 0)` see no change in values; only the certification strength of the Mathematica engine has increased.

## Side observations (non-blocking)

- The directive's F3 text references "Stage 062 transport map" but Codex consistently wrote "Stage 075 transport map" in both engines and in the new comments. The sympy docstring at .py:19 still reads "Stage 62". Since the user's hand-off explicitly says "stage-075 transport map" this looks like an intentional renumber to match the active stage map; flagging for potential follow-up if a downstream tracker still cites stage 62.
- The notes file still references `scripts/moving_throat_pde_stage141_*` per the original auditor's mention; out of scope for this verification (verifier is scripts-only).
- Banner relabel from `STAGE 073` to `STAGE 090` is in place in both engines (sympy line 47, mathematica line 32 — the verifier confirms via grep that no `STAGE 073` string remains in either script).

## Verdict justification

All three findings are resolved with edits that match the directive (F1, F2 by structure; F3 by substance with an intentional stage-number relabel from 062 to 075). Both exec logs exit `0` with the new substantive PASS lines (`rho_alpha - 4/3`, `zeta_req - 1/3`, the now-non-tautological identity, and `Pe_req = 0`). The Mathematica engine independently derives `rho_alpha` and `zeta_req` from the contact-plus-pole coefficients `(3/4, 1/4)` rather than hardcoding `4/3`, satisfying the checkpoint-stage two-engine bar. Outputs are fresh. No regressions. No new findings introduced.
