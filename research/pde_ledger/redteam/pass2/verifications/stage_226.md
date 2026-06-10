---
unit_id: 226
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 226

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim — stale verification pointer)

**Classification:** resolved (deferred, no script change required or possible)

**What changed:**
Nothing in scripts. The diff patch `redteam/pass2/exec_logs/stage_226_diff.patch` is empty (1 line); both `.py` and `.wl` are byte-unchanged from HEAD. The directive (`directives/stage_226.md`) correctly applies NOTHING: F1 is a paper-side card edit (`paper/stages/stage_226.tex:11` reads `... Mathematica audit: none yet.`) held for user resolution via `## Resolve before fix_loop`, and per the standing user decision it is deferred to PAPER_CLEANUP P4-51.

**Assessment:**
This is a stale STATUS annotation on the card, OUTSIDE scripts-only scope — not a script defect. The script side genuinely verifies MORE than the card admits: a complete second-engine `.wl` exists and passes all checks (`mathematica/output/...stage226...mathematica_audit.txt:93` "All Stage 226 Mathematica checks passed."). There is no math/value disagreement. Correctly classified as a deferred paper finding; it does NOT warrant needs_rework. No collateral edits exist (empty diff confirms this).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L46-rooted bridge proof: "Xi_load = Xi1 = N01/N0 - D01/D0" (genuine `simplify(Xi1 - Xi_load) == 0`, `Xi1` derived via `sp.diff(P0A,eps)`).
- "sigma_even = 2.67386816837172775", "sigma_transfer = 2.31561904386055442" (computed from `nullspace()` + QR projector, not hardcoded).
- "All Stage 226 checks passed."

**Mathematica:** exit=0. Notable lines:
- "OK  M1 Xi1 equals load-defect bridge" (via `Series[p0A,{eps,0,1}]`, a distinct expansion route from the `.py`'s `sp.diff`).
- "OK  M3 exact one-pole identity" (real residual `u4Value - 4 u2Value^2 == 0`, not assumed).
- "transfer null vector 1 constraint residuals = {0, 0, 0}" / "...vector 2... = {0, 0, 0}" (non-vacuous residual check on the actual null basis).
- "sigma_even = 2.67386816837172774819...", "sigma_transfer = 2.31561904386055441584..." (agree with SymPy to ~20 digits).
- "All Stage 226 Mathematica checks passed."

**Output freshness:** Orchestrator re-ran both engines directly (logs dated 2026-06-09T19:21:18, exit 0, deterministic; committed 226 `.txt` byte-identical to fresh). Outputs are fresh; no `stale_output`. Did NOT fail on committed `.txt` mtime.

## Independence spot-check (Mathematica)

Confirmed the `.wl` is a genuine independent recomputation, not an echo:
- M3 derives D0/D2/D4/u2/u4/P0 from the physical primitives (cWall, gU/gW, delta, qMix, ...) natively in its own kernel (L119-166), then `FullSimplify`s — it does not retype the `.py`'s numbers.
- M5/M6 build `evenMatrixExact` (2x5) and `transferMatrixExact` (3x5) via `coeffRow` (`CoefficientArrays`) of its own `FullSimplify`-derived `k1Mixed`/`hEvenMixed`/`d01Mixed`...`d41Mixed` (L214-228), and the nullspaces come from Mathematica's own `NullSpace` (L233, L244). The `.wl` pins FEWER literals than the `.py` — it never re-types the matrix entries or the null-basis vectors (the `.py` pins `expected_even`/`expected_w`/`expected_t`), so it is strictly less of an echo.
- M1 uses `Series` (L83) vs the `.py`'s `sp.diff` (L42) — different expansion routes.
- The projector swap (`Orthogonalize[N[basis,80]]` then `qRows^T.qRows`, L48-52, vs the `.py`'s `QRdecomposition`, L251) is mathematically cosmetic — the orthogonal projector onto the null-space column span is unique regardless of which orthonormalization produced the basis — but the overall script is a real independent recomputation. No `mathematica_transliteration` concern.

## Material-change assessment

`material_change`: false. No edits were applied; both scripts are unchanged from HEAD and both engines pass at byte-identical committed outputs. No derived result changed, so no downstream unit (>226) is affected by this verification.

## Side observations (non-blocking)

None that affect the verdict. (The honest caveat noted by the auditor — M2's `hOnePole` bakes in the `D4/D0 = u2^2 - u4` form rather than deriving it from `hEvenComp` — is independently corroborated by M3's `u4 == 4 u2^2` residual and the actual `D4/D0` value, so the one-pole closed form is still genuinely cross-checked. Not a finding.)

## Verdict justification

Scripts are clean: the diff is legitimately empty (a zero-correction batch, not a failure), both engines exit 0 with full PASS bodies, and the load-bearing assertions are non-tautological (the bridge identity, the one-pole residual, the `{0,0,0}` null-basis residual check, and the rank/nullity/sigma solves are all real computations, not assumptions). The `.wl` is a genuine independent recomputation that pins fewer literals than the `.py`. The sole finding (F1) is a stale card-side verification pointer — a paper_misalignment outside scripts-only scope, correctly deferred to PAPER_CLEANUP P4-51 by standing user decision, requiring no script change. Per the contract, only a genuine script-side defect rolls to needs_rework; a deferred paper finding alone does not. Verdict: `verified`.
