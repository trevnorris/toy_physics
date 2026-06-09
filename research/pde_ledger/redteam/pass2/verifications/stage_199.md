---
unit_id: 199
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 199

## Per-finding outcomes

### F1 — stale_output (SymPy `.txt` carried pre-renumber `STAGE 182` banner)

**Classification:** resolved

**What changed:**
No source-code edit (none was required). The orchestrator re-ran the SymPy script to refresh the committed transcript. The refreshed `scripts/output/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.txt` now has mtime `2026-06-09 14:09`, newer than the `.py` (`2026-06-03 15:59`).

**Assessment:**
The refresh fully closes the finding. The refreshed transcript's line 3 reads `STAGE 199 — EXACT PAIRWISE ORBIT-TRANSPORT LAW AND THE TWO-POINT ORBIT-LOCK THEOREM` and the ledger banner at line 420 reads `STAGE 199 LEDGER`. The stale `STAGE 182` / `STAGE 182 LEDGER` labels are gone (no occurrence of `182` anywhere in the file). Every residual is `= 0`: §I three pairwise ratios (lines 26–28), §II `ln Phi_mu - ln expanded = 0` plus three same-orbit ratios (67–70), §III mismatch collapse triple (105–107), §IV q-chart (138–140), the projector split `Q_quot/O_orb - expected`, `Q+O-Δx`, `M_* O_orb = 0`, `M_* Q_quot - qpair = 0` all 8-vector/3-vector zeros (309–368), §V restoration map (373–375), §VI cocycle Φ/m/Δx/q (380–407), §VII reduction-to-198 (412–417). No tautology introduced (no code touched). The Mathematica `.txt` was already fresh and is also re-confirmed clean: `STAGE 199` banner (line 3) and `STAGE 199 MATHEMATICA LEDGER` (line 135), every line `PASS` (no FAIL token anywhere), mtime `2026-06-09 14:09`.

### F2 — paper-side card-text lag (`\stagefield{Verification}` still says "Mathematica audit: none yet")

**Classification:** resolved (disposition holds — USER-DEFERRED, non-blocking)

**What changed:**
Nothing in scripts. This was explicitly an informational, paper-side item with no Codex action and no directive entry (the original report states "No directive entry is written for this item (paper-side, no math discrepancy)"). It is deferred to paper-cleanup.

**Assessment:**
Out of verifier (scripts-only) scope; correctly carried as a non-blocking prose-currency lag. The script and paper math agree; only the bookkeeping sentence lags. Disposition holds; not gated by this audit.

## Exec log assessment

**SymPy:** exit=0 (inferred from refreshed transcript; every check prints `= 0`, no exception/traceback). Notable lines:
- L3: `STAGE 199 — EXACT PAIRWISE ORBIT-TRANSPORT LAW AND THE TWO-POINT ORBIT-LOCK THEOREM` (canonical banner restored)
- L67: `ln Phi_mu - ln expanded monomial form = 0`
- L138–140: `q_tr - (1+chi0_*) ln m_T = 0`, `q_eta + ln m_K = 0`, `q_nt - [...] = 0`
- L357/363: `M_* O_orb Delta x = [0,0,0]`, `M_* Q_quot Delta x - qpair = [0,0,0]`
- L420: `STAGE 199 LEDGER`

**Mathematica:** exit=0. Notable lines:
- L10: `PASS: derived M_* rows - SymPy compiler rows` (M_* built from primitive weight vectors — independent route)
- L28: `PASS: ln Phi_T solve - ln SymPy transport` (Φ from native `Solve`)
- L73/75/77/79: `PASS: Q^2 - Q`, `PASS: O^2 - O`, `PASS: Q O`, `PASS: O Q` (augmented-system `LinearSolve` projectors + idempotency/orthogonality — independent route, no `.py` counterpart)
- L114: `PASS: det mismatch-to-q chart - (1+chi)` (extra determinant invariant, no `.py` counterpart)

**Output freshness:** confirmed. Both saved `.txt` outputs have mtime `2026-06-09 14:09`, newer than the SymPy `.py` (`2026-06-03 15:59`) and the Mathematica `.wl` (`2026-06-01 12:14`). Refresh landed; banners are canonical `STAGE 199` on both engines.

## Material-change assessment

`material_change`: false. No source-code change was made — only a transcript re-generation. No derived result changed; downstream units cannot be affected.

## Side observations (non-blocking)

None beyond the already-noted F2 paper-side card-text lag, which is USER-DEFERRED to paper-cleanup and out of scripts-only scope.

## Verdict justification

Both findings are resolved. F1's stale-output condition is fully cleared: the refreshed SymPy transcript carries the canonical `STAGE 199` banner (line 3) and `STAGE 199 LEDGER` (line 420), the stale `STAGE 182` labels are gone, and every residual is `= 0` with no code edited (so no tautology risk). The Mathematica transcript is re-confirmed fresh, `STAGE 199` throughout, every assertion `PASS`. The disposition holds: the `.wl` is genuinely independent — M_* built from primitive weight vectors, Φ derived via native `Solve`, projectors via augmented-system `LinearSolve` with idempotency/orthogonality and a determinant invariant absent from the `.py` — and value-reconciliation showed 0 misalignments. F2 is a non-blocking, USER-DEFERRED paper-side card-text lag, correctly outside scripts-only scope. Verdict: `verified`.
