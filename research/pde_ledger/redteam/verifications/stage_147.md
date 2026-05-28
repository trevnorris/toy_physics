---
unit_id: 147
batch: IV.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 147

## Per-finding outcomes

### F1 — missing_verification_script (sympy)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:52-75` adds the paper-literal anchors and chain-rule cross-check exactly as the directive prescribes (`AT_paper = sp.Float("-4.27263956256927", 30)`, `BT_paper`, `ratio_paper`, three numerical assertions, and the chain-rule reconstruction `AT_chain = sp.N(-dTm_dSigma * dSigma_dPi_at_star / gp_star, 30)`). The orchestrator's post-hoc fix is visible at line 72: `AT_30 = sp.N(AT, 30)` is now defined and compared against `AT_chain` at full 30-digit precision (line 73), replacing the prior `sp.N(AT)` default-15-digit truncation that produced the ~5e-16 residual. The centered-kernel structural assertion (`Wcenter_const == Wcenter_const_expected`) is present at lines 88-99.

**Assessment:**
The edits match the directive verbatim plus the documented orchestrator fix. The numerical anchors are non-tautological (independently computed left-hand sides compared against paper-quoted literals); the chain-rule check pairs an independently assembled `dTmDSigma` with `dSigma_dPi_at_star` and divides by `gp_star`, exercising a different evaluation path than the closed-form `AT` line; the centered-kernel constant-offset extraction is symbolic and not algebraically trivial. SymPy exec log shows all five `PASS:` lines for F1 contents (transcript lines 11-14 and 126). No collateral edits beyond what the directive named.

### F2 — tautological_check (mathematica)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:63-91` replaces the `R_q(g_minus) - 1/4` tautology with the three numerical anchors (`A_T vs paper`, `B_T vs paper`, `|A_T|/B_T vs paper`), the chain-rule cross-check (`aTChain` independently assembled from `dTmDSigma`, `dSigmaDPi`, `gPrimeStar`), and the centered-form constant-offset check. The Cluster-A F2(b) banner fix is in place at line 26 (`banner["STAGE 147 — FIRST-ORDER RIGIDITY KERNEL"]`). The F2(a) `Chop[wCenterConst - wCenterConstExpected, 10^-25]` wrap at line 91 is documented in the prior orchestrator notes and is the mechanism that lets `expectZero` see `0 === 0`.

**Assessment:**
The tautology is gone — there is no more `expectZero["R_q(g_minus)-1/4", ...]` line. The replacement assertions are anchored to externally-quoted literals and an independently assembled differential identity. The Mathematica exec log (lines 15-26) shows the corresponding `... = 0` / `PASS:` pairs for all six assertions (three paper anchors, chain-rule, centered form, moment stability for `g_*` and `S_*`). Banner at log line 3 correctly reads `STAGE 147 — FIRST-ORDER RIGIDITY KERNEL`.

### F3 — insufficient_verification (centering / kernel structure)

**Classification:** resolved

**What changed:**
SymPy script lines 101-113 add `g_star_resub = sp.N(gPi.subs(Pi, Pi_star), 40)`, `S_star_resub = sp.N(Sformula.subs(Pi, Pi_star), 40)`, and the two drift assertions (`< 1e-30`). Mathematica script lines 93-99 add the matching `gStarResub`, `sStarResub` with `expectZero[..., If[Abs[...] < 10^-30, 0, ...]]` wrappers. Both engines now exercise the moment-stability identity.

**Assessment:**
The resubstitution check is non-tautological: the cached `g_star`/`S_star` (from lines 25-26 in SymPy / 41-42 in Mathematica) are compared against fresh `.subs`/`/. p -> pStar` evaluations performed after the kernel-assembly block, so any precision drift between the family-1 anchor block and the kernel-assembly block would surface as a nonzero residual. Both exec logs show `PASS:` for both stability checks (sympy line 127; mathematica lines 27-30). The centered-kernel structural assertion under F1/F2 supplies the kernel-shape verification the auditor identified as the load-bearing half of the stage's claim.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `PASS: A_T matches paper-quoted -4.27263956256927 to 1e-12` (line 11)
- `PASS: A_T closed form agrees with chain-rule decomposition (residual < 1e-20)` (line 14)
- `PASS: Centered kernel W_*(x) has form A_T(c - g_*) + B_T(K_q - S_*)` (line 126)
- `PASS: g_*, S_* moment values stable across audit (drift < 1e-30)` (line 127)

**Mathematica:** exit=0. Notable lines:
- `A_T vs paper -4.27263956256927 = 0` / `PASS:` (lines 15-16)
- `A_T closed form vs chain-rule route = 0` / `PASS:` (lines 21-22)
- `W_*(x) centered form A_T(c-g_*) + B_T(K_q-S_*) = 0` / `PASS:` (lines 25-26)
- `g_* resubstitution drift = 0` / `PASS:` (lines 27-28)
- `Stage 147 Mathematica audit passed.` (line 32)

**Output freshness:** Confirmed. SymPy output (2026-05-27 20:09:55) is newer than the script (2026-05-27 20:08:13); Mathematica output (2026-05-27 19:59:45) is newer than the script (2026-05-27 19:57:07). Both transcripts contain the new `PASS:` lines named in the directive, including the orchestrator's chain-rule fix.

## Material-change assessment

`material_change`: false.

No numerical result changed. The edits added assertions and renamed a banner; no closed-form expression, no derived constant, no symbolic identity was altered. `A_T = -4.27263956256927...`, `B_T = 0.134875005736706...`, `|A_T|/B_T = 31.67851...` printed in the post-fix transcripts match the pre-fix values to all stated digits and match the paper-quoted literals from the appendix.

## Side observations (non-blocking)

- The directive's `## Applied: Fn` / `## Blocked: Fn` blocks were not appended to `redteam/directives/stage_147.md` by Codex; the edits are nevertheless verifiable directly against the file state. The orchestrator may want to backfill these blocks for audit-trail symmetry across stages.
- The Mathematica `Chop[..., 10^-25]` wrap on `wCenterConst - wCenterConstExpected` (line 91) is a precision-cushion against the residual `~10^-29` from finite-precision `N[...,40]` arithmetic; this is what makes `expectZero` see exact `0`. The pattern is documented from earlier batches and is not a tautology, since the unchecked residual would visibly exceed `10^-25` if either coefficient or the canonical-branch offset were perturbed.

## Verdict justification

All three findings are resolved with edits matching the directive (plus the documented `AT_30 = sp.N(AT, 30)` orchestrator fix and the F2(a) `Chop` wrap and F2(b) banner correction). Both engines exit 0; both transcripts contain every `PASS:` line named in the directive; no result downstream of stage 147 changes since no derived numerical or symbolic quantity was altered. Verdict: `verified`.
