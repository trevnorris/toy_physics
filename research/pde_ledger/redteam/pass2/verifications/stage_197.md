---
unit_id: 197
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T22:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 197

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No source-code edit (the directive correctly required none for Codex). The orchestrator's independent re-run refreshed the committed SymPy transcript `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.txt`. Both `.txt` outputs now carry mtime `Jun 9 14:07`.

**Assessment:**
The stale pre-renumber banner is gone. The refreshed SymPy output reads `STAGE 197 — EXACT CONDITIONAL PACKET-A CLOSURE THEOREM` (line 3) and `STAGE 197 LEDGER` (line 151); a `grep "STAGE 180"` on the file returns nothing (count 0). The Mathematica output banner reads `STAGE 197 - CONDITIONAL PACKET-A CLOSURE THEOREM` (line 3) and `STAGE 197 MATHEMATICA LEDGER` (line 105), also `STAGE 180`-free. The verification criterion from the report ("after re-run, grep STAGE 180 returns nothing") is met. This was a cosmetic/banner-only staleness — the symbolic content already agreed across engines and is unchanged.

## Exec log assessment

**SymPy:** exit=0 (orchestrator re-run). Refreshed `.txt` content, notable lines:
- `chi_Q extractor - deformation algebra formula = 0` (L77)
- `Delta_norm - P0^target(1/chi_Q - 1) = 0` (L78)
- `(3S - Sigma0)(chi_Q - 1) - closure numerator = 0` (L119)
- All Section I–VII identity checks emit `= 0`; banner reads `STAGE 197` (L3, L151).

**Mathematica:** exit=0 (orchestrator re-run). Refreshed `.txt` content, notable lines:
- `PASS: chi_Q extractor - SymPy deformation algebra formula` (L39)
- `Delta_norm == 0 iff chi_Q == 1 = True` → `PASS` (L62–63)
- `N_Q == 1 iff chi_Q == 1 = True` → `PASS` (L64–65)
- `Delta_norm at chi_Q = 6/5 is nonzero = True` → `PASS` (L71–72, negative control falsifiable)
- 24 `PASS:` lines, zero `FAIL`; banner reads `STAGE 197` (L3, L105).

**Output freshness:** confirmed. Both `.txt` files mtime `Jun 9 14:07`, newer than the source scripts; both banners now canonical `STAGE 197`; `grep "STAGE 180"` returns 0 on both.

## Material-change assessment

`material_change`: false. No source-code change occurred — only a committed-transcript refresh that corrected a cosmetic banner. The derived results (`χ_Q`, `Δ_norm`, closure gate, linearized map, higher-odd irrelevance) are byte-identical in substance to what downstream units already carried. No downstream unit is affected.

## Side observations (non-blocking)

The disposition holds as stated by the auditor. The `.wl` is genuinely independent, not a transliteration: it derives the DtN fingerprint natively via `FunctionExpand[z*D[SphericalHankelH1[2,z],z]/SphericalHankelH1[2,z]]` then `Series` to `z⁵`, `Solve`s the canonical-even matching for `{Σ₂,Σ₄}` with a uniqueness assertion, and re-derives `P₀^target` from the `Γ₅` normalization — whereas the `.py` hardcodes the expanded `Λ`-coefficients, the `Σ`-match closed forms, and the `P₀^target` literal. The iff is proved by `Reduce`/`Equivalent` over `Reals` on the `.wl` side vs. hand-built algebraic identities on the `.py` side. The output confirms this divergent choreography (e.g. mathematica L29–34 native canonical-even matching cross-checked against SymPy literals; L54–55 `Reduce[... ] = chiFree == 1`). 0 value-reconciliation misalignments (6 deliverables + the equivalence all MATCH).

The previously noted paper-side card-text lag (card still reads "Mathematica audit: none yet" despite the present, passing `.wl`) is a paper/prose file outside red-team scope; it is USER-DEFERRED to paper-cleanup and is non-blocking for this verification.

## Verdict justification

The single informational `stale_output` finding (F1) is resolved: the orchestrator re-run refreshed both committed transcripts, the canonical `STAGE 197` banner is now present and the stale `STAGE 180` banner is fully purged (grep count 0). Every check across both engines reports `= 0`/`PASS`/`True` with the negative control at `χ_Q=6/5` returning a nonzero defect (falsifiability intact), exit 0. The disposition holds — the `.wl` is genuinely independent (native `SphericalHankelH1` + `Solve` route vs. `.py` hardcoded closed forms) and value-reconciliation shows 0 misalignments. The card-text lag is correctly USER-DEFERRED and non-blocking. Verdict: `verified`.
