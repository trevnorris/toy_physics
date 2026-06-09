---
unit_id: 179
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 179

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Classification:** resolved

**What changed:**
The `.wl` slope route was re-authored. The old mirror — a hand-coded epsilon
perturbation `n0A = pA^2/deltaA^2` with `nuDirect = (D[Log[n0A], eps] /. eps -> 0)/lam`
(transcribing `.py`'s `series(log(N0A)).coeff(eps,1)/lam`) — is GONE. There is no `eps`,
`lam`, or rebuilt `n0A` anywhere in the new file. In its place (wl:28-31) a directional
log-derivative operator `logSlope[expr, pairs] = Total[ w_i * x_i * D[expr,x_i]/expr ]`
is applied directly to the closed-form shape `shape = (xWall + xCoupling*xMixed)/(1 - xCoupling^2)`
(wl:62): `tauFromShape = logSlope[shape, ...]` (wl:63-66) and
`nuFromFactor = logSlope[kk*shape^2, ...]` (wl:75-78).

**Assessment:**
Correct and genuinely independent. The `.py` still reaches the slope by perturbing the
six primitives with `(1 + eps*lam*.)` and Taylor-expanding `log(N0A)` (py:61-73, unchanged);
the `.wl` now reaches it by an analytic directional logarithmic derivative of the
closed-form shape — different objects (shape variables vs perturbed primitives), different
mechanism (attached-weight `x*D/expr` operator vs `eps`-series coefficient). This is not a
cosmetic relabel or a port: the construction the directive said to remove is absent, and a
can-fail independent construct is visibly present. D1 (`staticTransfer/k - shapeFromHats^2`,
wl:54), D2 (`nuFromFactor - (kappa1+2 tauFromShape)`, wl:83), and D3 (slippage equivalence,
wl:96-97) all rebuilt on the new shape objects. No collateral damage to the deliverable forms.

### F2 — insufficient_verification (both engines)

**Classification:** resolved

**What changed:**
Both engines' weighted-defect blocks are rebuilt. `.py` (py:107-137): the free symbols
`tau1,tau2,tau3` are gone; instead three port substitution dicts feed `tau_i = tau.subs(port_i)`
(the A2-validated closed-form `tau` from py:84) and `nu_i = nu_direct.subs(port_i)` (the
independently-extracted slope from py:73). Per-port `nu_i - (kappa1 + 2*tau_i) == 0`
(py:131-133) is asserted, then the collapse `Xi = sum rho_i(nu_i - kappa1)` vs
`2 sum rho_i tau_i` (py:135-136). `.wl` (wl:120-145): `tauPorts = tauFromShape /. portRules`,
`nuPorts = nuFromFactor /. portRules`, per-port identity checks (wl:133-139), then the same
collapse (wl:141-144).

**Assessment:**
Load-bearing in both engines now. The `Xi` and `Xi_expected` sides are no longer two copies
of the same `kappa1 + 2*tau_i` template over free symbols — `nu_i` is the actual port slope
and `tau_i` is the validated closed-form, distinct per port. The added per-port
`nu_i - (kappa1 + 2 tau_i)` assertions (out sympy:81-83, wl out:37-42) make the block fail
first if the D2 relation regresses, exactly the false-confidence gap the finding raised. Not
tautological. All four sub-checks `= 0` / PASS.

### F3 — stale_output

**Classification:** resolved

**What changed:** Both committed `.txt` regenerated post-fix.

**Assessment:** Output mtimes (Jun 9 00:21) are newer than script mtimes (Jun 8 22:18).
Line 8 of both transcripts reads `STAGE 179 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM`; the
stale `STAGE 162` banner and out-of-date `nu_direct` form are gone. `.py` was changed ONLY
in the weighted-defect block (diff py:107-137); D1/D2/D3 untouched.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `N0/K - T^2 = 0`
- `nu_direct - (kappa1 + 2 tau) = 0`
- `nu_1 - (kappa1 + 2 tau_1) = 0` ... `nu_3 - (kappa1 + 2 tau_3) = 0`
- `Xi_1 - 2 weighted tau = 0`

**Mathematica:** exit=0. Notable lines:
- `PASS: N0/K - T^2`
- `PASS: nu_factor - (kappa1 + 2 tau)` (slope reached via `logSlope`, printed `tau_shape`/`nu_from_factor` are closed-form shape rationals, not eps-series)
- `PASS: nu_1...`/`nu_2`/`nu_3 - (kappa1 + 2 tau_i)`, `PASS: Xi_1 - 2 weighted tau`
- `Stage 179 Mathematica audit passed.`

**Output freshness:** confirmed — both `.txt` mtimes (Jun 9 00:21) post-date the script mtimes
(Jun 8 22:18); banners read STAGE 179.

## Material-change assessment

`material_change`: false. Every deliverable form (`N0^(r)=K T_r^2`,
`T_r=(Ghat_W + Rhat Ghat_U)/(1-Rhat^2)`, `nu_r = kappa1 + 2 tau_r`,
`tau_r = m + I/(1+I) i + H/(1-H) h`, `Xi_1 = 2 sum_r rho_r tau_r`) is byte-identical in the
carry-forward block and unchanged in value; only the verification METHOD changed (independent
`.wl` route + strengthened F2 in both engines). No derived result moved, so no downstream unit
is affected.

## Side observations (non-blocking)

- The banner/docstring still carry the pre-renumber "Stage 176/160/161" label (py:10,93;
  wl:85) and notes use "176/177/178". This is the known script/output numbering-band item
  deferred to the dedicated pass (per the auditor's NOTE); not a verification blocker.
- `logSlope` uses `D[expr,x]/expr` (= `D[Log expr]`) symbolically; it never evaluates `Log`
  at a perturbation point, so the two engines remain mechanically distinct despite both being
  "log-slope" in spirit.

## Verdict justification

All three findings resolved. The F1 re-author is genuinely independent — the
`D[Log[n0A]]/.eps->0` mirror is removed and replaced by a directional log-derivative operator
on the closed-form shape, a can-fail construction the `.py` does not use; this clears the
173-lesson "insufficient re-author" risk. F2 is now load-bearing in both engines (per-port
`nu_i = kappa1 + 2 tau_i` fed from validated identities, not duplicated templates) and still
passes. All deliverables verify with residual 0, no value moved (material_change false),
banners read STAGE 179, and both engines exit 0. Verified.
