---
unit_id: 096
batch: IV.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T23:35:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 096

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
SECTION II in both engines was de-tautologized exactly as the directive specified.

- SymPy `scripts/...stage096..._sympy_audit.py:88-104`: the integer-literal hardcoding
  (`eps_2 = sp.Integer(0)`, `eps_4 = sp.Integer(0)`) is replaced by free symbols
  `eps_2, eps_4 = sp.symbols("eps_2 eps_4", real=True)` and a general formula
  `c_pole_gen = (1 + eps_4) / (4 * (1 + eps_2)**2)` (:88-89). Two CAN-FAIL probes
  are added (:90-97): `c_pole_gen.subs({eps_2:0, eps_4:1}) - 1/2` and
  `c_pole_gen.subs({eps_2:1, eps_4:0}) - 1/16`. The static limit is then taken
  FROM the general formula: `c_pole = sp.simplify(c_pole_gen.subs({eps_2:0, eps_4:0}))`
  (:104). All downstream definitions (`c_geom`, `Yhat_cons`, `Yhat_expected`,
  `rho_alpha`, `zeta_req`) and the five existing deliverable `expect_zero` assertions
  (:120-127) are byte-for-byte unchanged.
- Mathematica `mathematica/...stage096..._mathematica_audit.wl:57-61`: mirror change —
  `cPoleGen = (1 + eps4)/(4*(1 + eps2)^2)`, two `expectZero` probes for `1/2` and
  `1/16` (:58-59), and `cPole = FullSimplify[cPoleGen /. {eps2 -> 0, eps4 -> 0}]` (:61).
  The five deliverable `expectZero` checks (:78-82) are unchanged; the cosmetic
  `Print["eps2 = ", fmt[0]]` / `fmt[0]` preserves the prior eps printout (:76-77).

**Assessment:**
Correct and complete. The change matches the directive's required form in both engines,
symmetrically. The two new probes are genuinely non-tautological: I hand-checked the
literals against the helper definitions. `(1+1)/(4·(1+0)^2)=1/2` and
`(1+0)/(4·(1+1)^2)=1/16`, both distinct from `1/4`, so a wrong power on `(1+eps_2)`,
a wrong factor of 4, or an eps_2/eps_4 swap would drive a nonzero residual and trip
the gate. The gates themselves are can-fail: SymPy `expect_zero` (py:32-36) raises
`AssertionError` when `residual != 0`; Mathematica `expectZero` (wl:20-24) takes the
pass branch only on strict `res === 0` and otherwise calls `fail`→`Exit[1]` (wl:14-17),
so there is no nan/Piecewise/symbolic-leakage escape. The static limit is now derived
from the general formula rather than substituted as a literal, so the five deliverable
assertions are still arithmetic at eps=0 but are now fed by a structurally-verified
formula. SECTION I (orthogonality integrals + Laplace eigenvalue 6 + firewall
decoupling) is untouched in both engines (diff scope confirms; py:73-84, wl:46-55
unchanged). No collateral edits beyond the directive. The deliverable values
(1/4, 3/4, 4/3, 1/3, Yhat closed form) are preserved in both outputs.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `c_pole|eps_4=1,eps_2=0 - 1/2 = 0` (new probe PASS)
- `c_pole|eps_2=1,eps_4=0 - 1/16 = 0` (new probe PASS)
- `c_pole - 1/4 = 0`, `c_geom - 3/4 = 0`, `rho_alpha - 4/3 = 0`, `zeta_req - 1/3 = 0`,
  `Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)] = 0` (deliverables still PASS).
- All 15 SECTION I residuals = 0.

**Mathematica:** exit=0. Notable lines:
- `PASS: c_pole|eps4=1,eps2=0 - 1/2`
- `PASS: c_pole|eps2=1,eps4=0 - 1/16`
- `PASS: c_pole - 1/4`, `PASS: c_geom - 3/4`, `PASS: rho_alpha - 4/3`,
  `PASS: zeta_req - 1/3`, `PASS: Yhat_Q^cons - [...]`.
- All 15 SECTION I checks PASS; `Stage 096 Mathematica audit passed.`

**Output freshness:** confirmed. Both `.txt` outputs (mtime 23:22:01) are newer than
both edited scripts (.py 23:16:02, .wl 23:16:02), so the saved outputs were
regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit only de-tautologizes the verification path and adds
two can-fail probes; every deliverable value (eps=0, c_pole=1/4, c_geom=3/4,
rho_alpha=4/3, zeta_req=1/3, Yhat closed form) is bit-for-bit identical to pre-fix
output. No derived result that downstream units consume has changed.

## Side observations (non-blocking)

1. The cosmetic `Print["eps2 = ", fmt[0]]` / `Print["eps_2 =", eps_2_val]` retain the
   eps=0 display so the ledger printout is unchanged; this is the intended approach
   (eps_2/eps_4 are now free symbols, the printed 0 is the chosen static value). Fine.
2. The SymPy docstring (py:9) still attributes the formula to "the Stage 092 obstruction
   formula" vs the notes' "Stage-75" — this is a cross-reference provenance label, not a
   self-label, already deferred by the auditor to the dedicated numbering pass. Not a
   verification issue; noted only for continuity.

## Verdict justification

The single finding (F1, insufficient_verification at the checkpoint bar) is fully
resolved. Codex applied the directive verbatim and symmetrically in both engines:
the obstruction formula is now a free-symbol expression exercised by two distinct
can-fail literals (1/2 and 1/16, both ≠ 1/4) that would catch a wrong power, factor,
or symbol swap, and the static limit is derived from that general formula. SECTION I
is untouched, both engines exit 0, outputs are freshly regenerated, all card deliverable
values are preserved, and no regression appears in the diff. Verdict: verified.
