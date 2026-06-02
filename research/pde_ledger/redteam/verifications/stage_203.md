---
unit_id: 203
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 203

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
Section VI of both scripts was augmented to exercise the headline crossing theorem / slice scalar on the genuine graph lift instead of a hand-zeroed bare `\beta^5`.

- SymPy (`..._sympy_audit.py:282-339`): the demonstration path now moves TWO free coordinates — `gamma` via `beta_path = 2**(2*tau-1)` (line 282) and `cetaU` via `cetaU_path = cetaU_bar*exp(rho*tau)` (line 284); the other three (`lam_bar, KU_bar, KW_bar`) stay inert. The §IV graph maps `T_graph, Keta_graph, mu_graph` are substituted onto this path (`graph_path_subs`, lines 286-296), the three target monomials `Ctr_graph_lift, Cnt_graph_lift, epsEta_graph_lift` are rebuilt on the lift (lines 297-306), and three NEW assertions are added: `expect_zero("graph-lift target monomial q_tr/q_nt/q_eta", ...)` (lines 337-339) where `qtr_graph_lift = normalize(sp.log(Ctr_graph_lift / Ctr_tgt))` etc. The IVT scalar is now sourced from `beta_lift = y_graph_path[2]/gamma_bar` (line 315), i.e. the graph composition, not a fiat `beta_path`.
- Mathematica (`..._mathematica_audit.wl:265-308`): same two-coordinate path (`betaPath = 2^(2 tau-1)`, `cetaUPath = cetaUBar Exp[rho tau]`), graph maps lifted via `graphPathRules` (line 269), three NEW `expectZero["graph-lift target monomial q_tr/q_nt/q_eta", ...]` (lines 306-308), and `betaLift = yGraphPath[[3]]/gammaBar` (line 292) feeding the IVT scalar.

No collateral edits beyond §VI: the diff touches only the §VI block of each file (additions plus the `beta_path`→`beta_lift` swap, which is guarded by the retained `beta_path - gamma(tau)/gamma_bar = 0` equality check).

**Assessment:**
The fix matches the directive's claim manifest (M1 monomial-invariance on the lift, M2 crossing on the graph-sourced scalar). The new checks are FALSIFIABLE, not by-construction tautologies:

- `Ctr_graph_lift` is rebuilt as the FORWARD monomial `(gamma cetaU/KU)^(1+deltaUs) * (pi^2 T_graph_lift/(L^2 KU))^(1+chi0s)`, while `T_graph_lift` traces back to the INVERSE construction `deltaU_graph = (Ctr_tgt/(gamma cetaU/KU)^(1+deltaUs))^(1/(1+chi0s))` from §IV (sympy:194-198). The residual `log(Ctr_graph_lift/Ctr_tgt)` cancels to 0 only because the forward exponents `(1+deltaUs),(1+chi0s)` match the inverse graph-map exponents. A wrong dependent-coordinate exponent in `T_graph`/`deltaU_graph` (e.g. `chi0s` instead of `1+chi0s`) would leave a nonzero residual — exactly the defect the original report said §VI could not catch. Same logic holds for `q_nt` (depends on `Estar, Fstar` matching) and `q_eta` (depends on `Keta_graph`).
- The path moves ≥2 free coordinates and BOTH appear non-trivially in the residuals: `cetaU` enters `q_tr` and `q_eta` (and via `Keta_graph` into `q_nt`), `gamma` enters `q_tr` and `q_nt`. So the cancellation is a genuine multi-coordinate identity, not a single-coordinate path that could mask an untouched exponent (the self-test concern raised in the original report).
- The `Sigma0=Sigma5=0` reduction is now justified in-script by a comment tying it to the verified monomial invariance, and the load-bearing assertions are the three residual-bearing invariance checks (per directive item 1/2), not the bare reduction. The IVT sign/root/uniqueness checks remain and now apply to the graph-sourced `beta_lift`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `graph-lift target monomial q_tr = 0`, `q_nt = 0`, `q_eta = 0` (log lines 312-314) — the new invariance checks pass.
- `widehat Delta_Q(0) = -31/32`, `widehat Delta_Q(1) = 31`, `real crossing set = {1/2}` (lines 317-324); sign checks `< 0` / `> 0` and root `= 0` all confirm. No traceback; `STAGE 186 SYMPY AUDIT PASSED`.

**Mathematica:** exit=0. Notable lines:
- `PASS: graph-lift target monomial q_tr`, `q_nt`, `q_eta` (log lines 95, 97, 99) — independently derived in log-space.
- `PASS: real crossing set` (line 113), `widehat Delta_Q(0) = -31/32`, `widehat Delta_Q(1) = 31` (lines 104-105). 24 named `PASS:` lines, 0 `FAIL:`; `STAGE 186 MATHEMATICA AUDIT PASSED`.

Pass-line count: Mathematica emits 24 named `PASS:` assertions (the SymPy engine asserts via `expect_*`/`AssertionError`, 26 such checks, all passing). The orchestrator's reported "25 PASS" is its own aggregate roll-up; the substantive count of named Mathematica assertions is 24 and every one is a check that can fail (none is `x - x = 0`). No previously-passing check was weakened or removed — the diff is purely additive (3 new monomial-invariance checks) plus the `beta_path`→`beta_lift` swap, which is itself guarded by a retained equality assertion.

**Independence (load-bearing):** confirmed genuine, NOT transliteration. SymPy forms the monomials multiplicatively in power-space then takes `sp.log(Ctr_graph_lift/Ctr_tgt)` (sympy:297-310). Mathematica forms the residual directly additively in log-space: `qtrGraphLift = (1+deltaUs)(Log[gammaPath]+Log[cetaUPath]-Log[KUBar]) + (1+chi0s)(2 Log[Pi] + logTGraphLift - 2 Log[L] - Log[KUBar]) - Log[CtrTgt]`, with `logTGraphLift = logTGraph /. graphPathRules` built from the additive `logDeltaGraph` (wl:196-205, 271-287). This honors the directive's anti-transliteration guard (sidesteps `PowerExpand` branch traps by living in logs from the outset) and matches the original report's established §IV independence pattern.

**Output freshness:** confirmed. Both `.txt` outputs (mtime 1780425178) are newer than their scripts (sympy 1780424387, wl 1780424411), i.e. regenerated post-fix.

## Material-change assessment

`material_change`: false.

The fix only ADDED stronger composition checks around the same correct crossing result. The IVT demonstration path was swapped (`1 + rho(2tau-1)/(1+rho)` → `2^(2tau-1)`) and now moves two coordinates, which changed the illustrative endpoint values (`Delta_Q(0)=-31/32`, `Delta_Q(1)=31`). But `beta_path`/`hat_delta_graph`/`hat_chi_graph` are confined to §VI and consumed only by the §VI demonstration assertions (grep confirms no export, no downstream reference); they are within-script illustration artifacts, not carried/derived constants. The crossing theorem's structural content (endpoint sign change ⇒ interior root) is unchanged, and no quoted/derived value a downstream unit depends on was altered. No new constant was introduced (reuses `Ctr_tgt, Cnt_tgt, epsEta_tgt`).

## Side observations (non-blocking)

- Stale stage banners persist (`STAGE 186`, `Stage 192/197/202` references, `chiFromStage180`/`closureNumStage180` symbol names in the .wl). The original auditor recorded these as informational provenance noise from a renumbering; no number affects a verified identity. Not part of F1; flagged only.

## Verdict justification

`verified`. The single finding F1 is resolved in both engines. Section VI now genuinely composes `chi_Q` through the Stage-202 graph map over a two-coordinate free-quintuple path and verifies target-monomial invariance on the lift via three new residual-bearing assertions that would fail on any wrong graph exponent — closing the exact gap the auditor identified (the headline crossing theorem / slice scalar no longer runs on an inert hand-reduced `beta^5`). The Mathematica check is an independent log-space derivation, not a port of the SymPy power-space construction. Both engines exit 0 with all checks passing, outputs are fresh, no check was weakened, and the change is additive with no downstream-carried value altered (`material_change: false`). Checkpoint bar met.
