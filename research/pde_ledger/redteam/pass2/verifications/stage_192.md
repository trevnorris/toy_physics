---
unit_id: 192
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 192

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No source-code change was required or made. The committed SymPy output
`scripts/output/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.txt`
was refreshed by the orchestrator's independent re-run. Its banner now reads the
canonical `STAGE 192` (line 3: `STAGE 192 — EXACT ORBIT/QUOTIENT PROJECTORS AND
THE MICROSCOPIC ORBIT-LOCK SPLIT`) and the ledger header reads `STAGE 192 LEDGER`
(line 626), matching `py:35` / `py:196`. The pre-renumber `STAGE 175` banner and
`STAGE 175 LEDGER` strings are gone.

**Assessment:**
Finding fully addressed. The fix is purely a banner refresh (renumber artifact);
the math content was already current. mtimes confirm the refresh landed: both
output `.txt` files are dated Jun 9 14:05, newer than the SymPy `.py` (Jun 3 15:59)
and the Mathematica `.wl` (Jun 1 11:09). No collateral edit; the `.py`/`.wl`
sources themselves were untouched, exactly as the directive specified.

## Exec log assessment

**SymPy:** exit=n/a (no exec log needed — finding was a re-run/refresh, not a code
change). The refreshed committed output shows every `expect_zero` residual as a
zero matrix/vector: `pivot block - expected = 0` (l.27–32), `det(P) - (1+chi0_*)
= 0` (l.35), `P_inv P - I = 0` and `P P_inv - I = 0` (l.52–63), `section -
expected = 0` (l.88–103), `M_* S - I_3 = 0` (l.104–109), all six projector
identities `Q²-Q / O²-O / QO / OQ / M_*O / M_*Q-M_*` = 0 (l.234–309), free-row
support (Q rows 0,1,2,3,5) = 0 (l.314–323), `Q Δx - expected dependent-triple
support = 0` (l.416–431), `M_* Δx_fail - q = 0` (l.432–437), `O Δx - expected
orbit law = 0` (l.502–517), `M_* Δx_orbit = 0` (l.518–523), `Δx - (orbit+fail) =
0` (l.524–539), and the representative/lock checks `M_*(Sq)-q / Q(Sq)-Sq / O(Sq) /
S*0 / orbit-lock packet` = 0 (l.564–623).

**Mathematica:** exit=n/a. The (already-fresh, also re-saved) Mathematica output
carries an explicit `PASS:` on every check — 27 PASS lines, zero failures.
Notable independent-route ties: `PASS: Q Delta x - constrained failure solve`
(l.52–53), `PASS: O Delta x - constrained kernel solve` (l.63–64), `PASS: Solve
laws for Q Delta x == 0 and M Delta x == 0 agree` (l.92–93). The section is built
by constrained `LinearSolve` (`S from constrained LinearSolve`, l.21) and only
then reconciled against the SymPy target (`PASS: section - SymPy target section`,
l.22–23), confirming the genuinely-independent route. Banner reads `STAGE 192 --
MATHEMATICA ORBIT/QUOTIENT PROJECTOR AUDIT` (l.3); ledger `STAGE 192 MATHEMATICA
LEDGER` (l.100).

**Output freshness:** confirmed. Both `.txt` outputs have mtime Jun 9 14:05,
newer than both scripts (SymPy `.py` Jun 3 15:59; Mathematica `.wl` Jun 1 11:09).
The stale-banner finding is closed.

## Material-change assessment

`material_change`: false.

No source-code change was made — only a committed-output banner refresh. No
derived result changed (the SymPy math was identical before and after; only the
stage label string differs). No downstream unit can depend on a banner string, so
no `upstream_stale` propagation is warranted on math grounds.

## Side observations (non-blocking)

- The audit's independence disposition holds on re-read of the refreshed
  Mathematica output: `S` is produced by a constrained `LinearSolve` on the
  augmented 8×8 system (output l.21) and the `… - SymPy target` lines are
  cross-engine tie-points, each paired with a from-construction check (`M S - I_3`
  l.26–27, `Q Delta x - constrained failure solve` l.52–53, the agreeing `Solve`
  laws l.92–93). This is independent derivation + cross-check, not transliteration.
  0 value-reconciliation misalignments confirmed.
- Paper-side card-text lag (`stage_192.tex:11` `\stagefield{Verification}` says
  "Mathematica audit: none yet" while a passing `.wl` exists) is the USER-DEFERRED
  paper-cleanup item — non-blocking, not routed to Codex, not a value/identity
  finding. Disposition holds.

## Verdict justification

The single low-severity `stale_output` finding is resolved: the orchestrator's
independent re-run refreshed the committed SymPy output, whose banner and ledger
header now read the canonical `STAGE 192`, matching the script source; output
mtimes are newer than both scripts. Both refreshed outputs show every check
passing (SymPy residuals all zero; Mathematica all 27 `PASS`), including the
independent-route cross-engine ties that establish the `.wl` as a genuine
re-derivation rather than a port. No source-code change was made, so
`material_change` is false. The card-text lag is the acknowledged USER-DEFERRED
paper-cleanup item and is non-blocking. Verdict: verified.
