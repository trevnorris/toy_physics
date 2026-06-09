---
unit_id: 190
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

# Verification — unit 190

This stage had NO source-code change. The audit raised two low-severity,
documentation-side findings: F1 `stale_output` (SymPy `.txt` banner read
`STAGE 173`, refreshed by the orchestrator re-run) and F2 `paper_misalignment`
(card `\stagefield{Verification}` says "Mathematica audit: none yet" despite a
passing, independent `.wl`). The USER has DECIDED to DEFER the F2 card-text lag
to the paper-cleanup tracker — it is a stale status annotation, not a
value/identity mismatch, and requires NO script or paper edit now. My job is to
confirm the output refresh landed and the disposition holds.

## Per-finding outcomes

### F1 — stale_output (SymPy banner)

**Classification:** resolved

**What changed:**
No script edit (none was prescribed; F1 was a verifier-refresh item). The
orchestrator's independent re-run regenerated the committed
`scripts/output/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.txt`
(mtime 2026-06-09 14:04, newer than the `.py` mtime 2026-06-03 15:59).

**Assessment:**
Correct and complete. The refreshed `.txt` now carries the canonical banner:
line 3 `STAGE 190 — DIRECT DEFECT VS DRESSING SPLIT, SUPPORT-BLINDNESS, AND THE
SCALAR NO-GO FILTERS` and line 148 `STAGE 190 LEDGER` — the stale `STAGE 173`
label is gone on both the header and the ledger heading. Every result line in
the refreshed output reads `= 0` (support-loaded reconstructions :33-35, `E -
(1 - epsiloneta) = 0` :36, the three `dln(...)/dzeta = 0` :37-39, the five
slippage-ledger residuals :58-62, `Xi_direct - (A_tr Sigma_tr + Sigma_nt) = 0`
:103, the triangularity/inverse/theorem block :124-133, and the P2 no-go
invariants :138-145). The only non-`= 0` lines are intentional symbolic
displays of `Xi_1`/`Sigma_tr`/`Sigma_nt` and the spoiled-packet negative-control
`dln/dzeta` (:40-43) — these are designed-nonzero, not assertion failures. The
captured math content matches the audit's expectation. Resolved.

### F2 — paper_misalignment (subtype: paper_missing_script_claim) — USER-DEFERRED

**Classification:** resolved (disposition: deferred-noted, non-blocking)

**What changed:**
No edit. The card-text lag ("Mathematica audit: none yet" in
`stage_190.tex:11` despite a complete, passing `.wl`) is a verification-prose
status annotation, not a value/identity mismatch. The USER has DECIDED to route
this class to the paper-cleanup tracker. Codex does not edit paper/, and the
verifier is scripts-only.

**Assessment:**
Correctly dispositioned. This is a stale status note about whether the second
engine exists — the `.wl` demonstrably DOES exist and DOES pass (see Mathematica
output below). No value-reconciliation impact: the audit checked 19 deliverable
values, 0 misaligned, and raised no `value_mismatch`. Per the user's decision,
F2 is USER-DEFERRED to paper-cleanup and does NOT block `verified`. Recorded
here as `card_lag: deferred-noted`, no action required now.

## Exec log assessment

**SymPy:** exit=0 (independent re-run; refreshed committed output is the
artifact reviewed). Notable lines from the refreshed
`...sympy_audit.txt`:
- `STAGE 190 — DIRECT DEFECT VS DRESSING SPLIT...` (:3) — canonical banner.
- `support-loaded T^2 reconstruction = 0` ... `dln(N_* loaded)/dzeta = 0`
  (:33-39) — support-blindness holds.
- `Xi_direct - (A_tr Sigma_tr + Sigma_nt) = 0` (:103); `I[x,x] - 7/10 eps^2
  x1^2 = 0` (:143) — direct-defect law and P2 no-go hold.
- `STAGE 190 LEDGER` (:148) — canonical ledger heading.

**Mathematica:** exit=0. The refreshed
`...mathematica_audit.txt` (mtime 2026-06-09 14:04, newer than `.wl` mtime
2026-06-01 11:32) carries `STAGE 190 -- DIRECT DEFECT VS DRESSING SPLIT` (:3)
and `STAGE 190 LEDGER` (:111), and EVERY check reads `= 0` immediately followed
by `PASS:`. Notable lines:
- `support-loaded T^2 reconstruction = 0` / `PASS:` (:11-12) and the three
  `dln(...)/dzeta` PASS lines (:19-24).
- the negative control passes the right way: `spoiled dln(R_target)/dzeta is
  not identically zero ... PASS:` (:25-26) and `spoiled exact witness at
  eps=1/3,zeta=1/2 = 0 / PASS:` (:27-28).
- `Xi_direct - (A_tr Sigma_tr + Sigma_nt) = 0 / PASS:` (:53-54); the full
  triangular-inverse block (det, three block derivatives, three LinearSolve
  inverses, three reconstructions, three theorems) all PASS (:61-88); the P2
  no-go block all PASS (:93-108). No `FAIL`, no error.

**Output freshness:** confirmed. Both committed `.txt` outputs have mtime
2026-06-09 14:04, newer than both source scripts (`.py` 2026-06-03 15:59, `.wl`
2026-06-01 11:32). The refresh landed post-fix.

## Independence assessment (V.3 retrofit core)

The audit's independence verdict holds: the `.wl` is GENUINELY INDEPENDENT, not
a transliteration of the `.py`. Confirmed against the audit's section-by-section
method comparison — the `.wl` uses distinct operations for the same targets:
logarithmic Euler operator (II) vs. exp-substitution-then-derivative;
`Coefficient`-projection extraction of `A_tr`/`Σ_nt` with an extra
no-`Σ_χ`-residual test (III) vs. posited-form check; `LinearSolve` inverse (IV)
vs. hand-written closed forms; projector-matrix + `Series` (V) vs. componentwise
+ `diff`; and a numeric falsification witness `46/133` (I) vs. the `.py`
symbolic-`≠0` control. The refreshed Mathematica output reflects these distinct
routes (e.g. "projected A_tr", "nontracking residual has no Sigma_chi",
"inverse by LinearSolve", "spoiled exact witness"). Independence + 0
value-reconciliation misalignments confirmed.

## Material-change assessment

`material_change`: false.

No source code changed; only the committed `.txt` banners were refreshed
(cosmetic), and the deferred F2 is paper-prose-only. No derived result that any
downstream unit could depend on was altered. No `upstream_stale` propagation is
warranted on math grounds.

## Side observations (non-blocking)

None. The refreshed outputs are internally consistent with the audit's
assertion inventory and value-reconciliation table; the negative controls fire
in the correct (nonzero / exact-witness) direction on both engines, so the
`= 0` results are substantive cancellations, not trivially-zero traps.

## Verdict justification

`verified`. The stale-output finding (F1) is resolved: the orchestrator re-run
refreshed both committed `.txt` files (mtimes 2026-06-09 14:04, newer than the
scripts), and both now carry the canonical `STAGE 190` header + `STAGE 190
LEDGER` heading with every check reading `= 0` / `PASS` / `True` and the
negative controls firing correctly. The `.wl` is confirmed genuinely
independent across all five sections with 0 value-reconciliation misalignments,
so the audit's CLEAN disposition holds. The card-text-lag `paper_misalignment`
(F2) is a stale status annotation, USER-DEFERRED to the paper-cleanup tracker —
it requires no script or paper edit now and does NOT block `verified`. No
regressions, no material change.
