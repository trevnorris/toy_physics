---
unit_id: 110
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 110

## Per-finding outcomes

The original auditor report `/var/projects/toy_physics/research/pde_ledger/redteam/reports/stage_110.md` carried `verdict: clean` with `findings_count: 0`. There is no directive file at `redteam/directives/stage_110.md` (consistent with no findings). Consequently there are no per-finding outcomes to classify; this verification only confirms that the auditor's flagged cosmetic note (banner / PASS-tag identifier drift to "STAGE 093") was addressed under the Cluster C banner sweep and that both engines still pass.

## Cluster C banner sweep

- SymPy `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.py:31` now reads `print('stage110: PASS')`. The stale `stage93: PASS` tag flagged in the auditor's cosmetic note is gone.
- Mathematica `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl:26` now reads `banner["STAGE 110 — EXPLICIT ISOTROPIC ROBIN OUTLET MODEL"];`. The stale `STAGE 093 — ...` banner flagged in the auditor's cosmetic note is gone.

Both edits match the user-supplied directive lines exactly. No other lines in either script were altered (the algebraic core — `Lambda_out`, `Lambda_R`, `Y_R`, the five closed-form assertions A1--A10 in the audit's assertion inventory — is unchanged).

## Exec log assessment

**SymPy:** exit=0 (terminal `stage110: PASS` printed; reaching that print requires all five `assert` statements at lines 26--30 to evaluate true). Notable lines:

```
c2 = -1/(3*rho - 9)
c4 = (4 - rho)/(9*(rho**2 - 6*rho + 9))
c5 = -1/(9*rho - 27)
chi_Q^R = -3/(rho - 3)
chi_Q^R linearized = rho**2/9 + rho/3 + 1
stage110: PASS
```

The trailing `stage110: PASS` tag is the renumbered tag — confirms the banner-sweep edit landed and that the script still exits cleanly.

**Mathematica:** exit=0 (script terminates with `Exit[0]` only after all five `expectZero` calls run; the transcript shows the explicit `PASS:` line for each). Notable lines:

```
STAGE 110 — EXPLICIT ISOTROPIC ROBIN OUTLET MODEL
PASS: c2 - 1/(9 - 3 rho)
PASS: c4 - (4 - rho)/(9 (3 - rho)^2)
PASS: c5 - 1/(27 - 9 rho)
PASS: chi_Q^R - 3/(3 - rho)
PASS: chi_Q^R linearized - (1 + rho/3 + rho^2/9)
Stage 110 Mathematica audit passed.
```

Banner line is the renumbered `STAGE 110 — ...` form, confirming the `.wl:26` edit landed.

**Output freshness:** both transcripts are newer than their scripts.
- SymPy script mtime 2026-05-27 15:08; SymPy output mtime 2026-05-27 15:18 (output newer).
- Mathematica script mtime 2026-05-27 15:08; Mathematica output mtime 2026-05-27 15:24 (output newer).

## Material-change assessment

`material_change`: false.

The only edits are a `print(...)` tag string and a `banner[...]` title string. Neither touches an algebraic symbol, a series order, an assumption, an assertion target, or any value that propagates downstream. Stages 111+ that reference Stage 110's `chi_Q^R = 3/(3 - rho_R)`, `c2 = 1/(9 - 3 rho)`, `c4 = (4 - rho)/(9 (3 - rho)^2)`, `c5 = 1/(27 - 9 rho)` see no change.

## Side observations (non-blocking)

None.

## Verdict justification

The auditor's report had zero findings; the only entry was a non-blocking cosmetic note pointing out residual STAGE-093 identifier drift in the SymPy `print` tag (`.py:31`) and the Mathematica `banner[...]` title (`.wl:26`). Both have been corrected to STAGE 110 / `stage110`. The five algebraic assertions in each engine (sympy `assert` block at lines 26--30; mathematica `expectZero` block at lines 49--53) are byte-for-byte unchanged. Both engines exit 0 with all five identities passing and outputs newer than scripts. No regressions, no collateral edits, no new findings raised. Verdict: `verified`.
