---
unit_id: 177
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 177

This unit was audited **clean** in pass-2 (0 findings, no directive written). Scripts were
NOT modified by Codex; only the committed `.txt` outputs were refreshed by the orchestrator's
reliability re-run. This is therefore a LIGHT confirmation of the clean verdict, not a
finding-by-finding rework check. The four confirmation gates (a)-(d) all hold.

## Per-finding outcomes

No findings were raised in the original report (`findings_count: 0`), so there is nothing to
classify per-finding. The report's positive claims were spot-confirmed below.

## (a) Audit verdict (clean) holds — independence re-confirmed

The `.wl` is a **genuinely independent route**, not a port of the `.py`. The decisive,
load-bearing independence is in the portwise-defect block:

- **`.py`** (sympy l.42, l.55, l.76): forms the load factor in **closed form**,
  `Lambda = (OU2*GW + R*GU)/(OU2*OW2 - R**2)`, then takes the slope of `Lambda_p**2/Kp` via
  `sp.series(..., eps, 0, 2).removeO()/eps`.
- **`.wl`** (l.66-68): forms the SAME amplitude from the **M/I/H factorization**,
  `lambdaSqOverKP = mCalP^2*(1+iCalP)^2/(1-hCalP)^2`, via `D[...,eps]/.eps->0`. The closed-form
  `lambdaP` IS defined (l.47) but is verifiably **never consumed** in the sigma block (grep
  confirms it appears only at its definition).

The two engines reach `sigma_r` through different symbolic objects, and the `.wl` makes their
equality an **explicit separate assertion** with no `.py` counterpart:
`expectZero["load-factor factorization lambda0^2/k = M^2(1+I)^2/(1-H)^2",
lambda0^2/k - mCal^2*(1+iCal)^2/(1-hCal)^2]` (wl l.62-63). The slope-extraction primitive
also differs (`series().removeO()/eps` vs `D[...,eps]/.eps->0`). Confirmed: not a transliteration.

## (b) Deliverable values still verify (residual->0 / PASS)

Every deliverable residual is `= 0` in the refreshed committed outputs, both engines:
d ln M/I/H = 0; `Sigma_{A,r} = lambda_A sigma_r = 0`; `grouped trace vanishes = 0`;
`a_Xi - eps Xi1/4 = 0`; `b_Xi - 3 eps Xi1/4 = 0`; `b_Xi - 3 a_Xi = 0`;
`P_A slope = P0*Xi_load = 0`; `(P1/P0) - Xi_load = 0`. The `.wl` additionally shows the bridge
`load-factor factorization ... = 0`. Every `.wl` line carries an explicit `PASS:`; closing line
`Stage 177 Mathematica audit passed.` No residual nonzero.

## (c) Banner / exit / FAIL gate

- Banner = `STAGE 177 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE` in BOTH committed outputs.
- Exec logs: sympy `# exit_code: 0`, mathematica `# exit_code: 0`.
- `grep -c FAIL` over both committed `.txt` outputs = 0 / 0. No FAIL anywhere.

## (d) Script diff empty

`redteam/pass2/exec_logs/stage_177_diff.patch` is 0 bytes (wc -l = 0). Scripts unchanged,
consistent with "audited clean, no Codex changes."

## Exec log assessment

**SymPy:** exit=0. `weak-axisymmetric d ln M/I/H = 0`; `Sigma_{A,r} = lambda_A sigma_r = 0`;
`(P1/P0) - Xi_load = 0`.

**Mathematica:** exit=0. `PASS: load-factor factorization lambda0^2/k = M^2 (1+I)^2/(1-H)^2`;
`PASS: Sigma_{A,r} = lambda_A sigma_r`; `Stage 177 Mathematica audit passed.`

**Output freshness:** confirmed. Both `.txt` mtimes = 2026-06-09 00:19:57; both script mtimes =
2026-05-30 01:09:47. Outputs are newer than scripts -> regenerated post-refresh.

## Material-change assessment

`material_change`: false. Scripts unchanged (empty diff); only output `.txt` regenerated with
identical content (same symbolic residuals, all zero). No derived result changed; no downstream
unit is affected.

## Side observations (non-blocking)

None.

## Verdict justification

`verified`. The clean audit verdict holds on the light re-confirmation: the `.wl` is a genuinely
independent route (M/I/H-factorized `lambdaSqOverK` plus an explicit bridge assertion the `.py`
lacks, distinct slope-extraction primitive, `lambdaP` defined-but-unused in the sigma block);
all deliverable residuals are 0 / every `.wl` check is `PASS` in the refreshed committed outputs;
both banners read STAGE 177; both engines exit 0; zero FAIL; and the script diff is empty.
