---
unit_id: 198
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 198

## Per-finding outcomes

The auditor reported `verdict: clean` with `findings_count: 0`. There were no
source-code findings to apply, so there are no per-finding `resolved/partial/…`
rows to roll up. The only non-source item flagged in the report was a stale
committed SymPy output banner (`STAGE 181` = 198−17, known pre-renumber drift),
which the auditor explicitly noted the orchestrator re-run would refresh. This
verification confirms that refresh landed and the clean disposition holds.

## Exec log assessment

**SymPy:** exit=0 (clean re-run; output regenerated). Canonical banner restored —
`STAGE 198 — EXACT FINITE ORBIT LAW FOR THE DEPENDENT TRIPLE` (out L3) and
`STAGE 198 LEDGER` (out L195); the stale `STAGE 181` banner is gone. Every
residual prints exactly `= 0`:
- `epsEta(Keta_orbit)/epsEta_star - 1 = 0`, `Ctr(TU_orbit)/Ctr_star - 1 = 0`,
  `Cnt(...)/Cnt_star - 1 = 0` (out L72–74) — deliverable (1).
- `Ctr ratio - m_T^(1+chi0_*) = 0`, `epsEta ratio - 1/m_K = 0`,
  `Cnt ratio - m_mu/(m_K m_T^F_*) = 0` (out L91–93) — deliverable (2).
- `M_* Delta x_mis - (q_tr,q_nt,q_eta) = [0,0,0]` (out L110–115) — deliverable (3).
- `T/K/mu_restore - orbit = 0` (out L157–159) — deliverable (4).
- `m_T/m_K/m_mu at q=0 - 1 = 0` (out L190–192) and direct/inverse chart
  consistency (out L174–176) — deliverable (5).

**Mathematica:** exit=0 (clean re-run; output regenerated). Banner
`STAGE 198 -- EXACT FINITE ORBIT LAW MATHEMATICA AUDIT` (out L3). Every check is
`PASS:` with `= 0` / `= {0,…}` / `= True`:
- `det(dependent log matrix) + (1+chiS) = 0` (out L12–13); the three
  `orbit agrees with SymPy target = 0` lines (out L14–19) — orbit produced by the
  independent `Coefficient`-matrix (out L10) + `LinearSolve` (out L11) route, then
  cross-checked against the SymPy closed forms.
- the three scaling-ratio `PASS` lines (out L30–35) — deliverable (2).
- `finite-ratio Jacobian agrees with SymPy M_* = {{0…},{0…},{0…}}` (out L44),
  computed as a `D[…]` Jacobian of `finiteLogPacket` (out L42) — deliverable (3).
- restoration column-`LinearSolve` `= {0,0,0}` and returns orbit (out L53–62) —
  deliverable (4).
- `finite orbit-lock equivalence in log coordinates = True` (out L78) —
  deliverable (5) biconditional, a real two-way `Equivalent`, not a tautology.

**Output freshness:** confirmed. Both `.txt` files re-generated at
2026-06-09 14:07, newer than the `.py`/`.wl` script sources (Jun 3). The stale
`STAGE 181` banner the auditor noted is gone from the refreshed SymPy output.

## Material-change assessment

`material_change`: false. No source-code edit was applied — this was a
clean-audit + output-refresh only. No derived result changed; downstream units
are unaffected.

## Side observations (non-blocking)

The paper-side card-text lag (`Mathematica audit: none yet` despite a passing,
genuinely-independent `.wl`) is out of red-team script scope and is
USER-DEFERRED to paper-cleanup. Non-blocking; does not affect this verdict.

## Verdict justification

`verified`. The audit had zero source findings; the only flagged item — a stale
committed SymPy output banner — is now resolved by the orchestrator re-run:
both engines carry the canonical `STAGE 198` banner with the stale `STAGE 181`
gone, both exit 0, and every check reads PASS / `= 0` / `= True` across all five
deliverables. The disposition holds: the `.wl` is genuinely independent (orbit
via log-linear `Coefficient`+`LinearSolve`, `M_*` via `D[]` Jacobian, restoration
via dependent-column `LinearSolve`), and the value reconciliation is complete with
0 of 10 deliverable values misaligned. The card-text lag is paper-side and
USER-DEFERRED, non-blocking.
