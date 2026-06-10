---
unit_id: 227
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 227

## Per-finding outcomes

The auditor report for unit 227 carries `verdict: clean` with `findings_count: 0`. There were no findings to resolve, no directive was written (confirmed: `redteam/pass2/directives/stage_227.md` does not exist — expected for a clean unit), and Codex applied no edits. The captured diff patch `stage_227_diff.patch` is empty (0 bytes) — this is legitimate, not a missing-capture failure: both the `.py` and `.wl` are unchanged from HEAD.

Because there are no `F1, F2, …` to classify, verification reduces to re-confirming the `clean` verdict still holds against the current scripts and the fresh exec logs.

### Re-confirmation of `clean` (no findings)

**Classification:** resolved (vacuous — nothing to fix; the clean audit is upheld).

**What I checked:**
- Re-read both current scripts in full: the SymPy `..._sympy_audit.py` and the Mathematica `..._mathematica_audit.wl`.
- Both exec logs end `# exit_code: 0` with every assertion line emitting `OK` / `passed`.
- The diff is empty; both scripts match HEAD.

**Assessment:** The load-bearing assertions are non-tautological and the `.wl` is a genuine independent recomputation, not a transliteration:

1. **Distinct first-variation machinery.** The `.py` dresses every primitive as `param*exp(eps*x)` and takes `sp.diff(EXPR_dressed, eps).subs(eps, 0)` (py:65–99). The `.wl` builds a hand-rolled logarithmic-derivative operator `deltaMixed[expr] := Total[(#[[1]] #[[2]] D[expr, #[[2]]]) & /@ logPairs]` over `logPairs` (wl:110–114) — `param·∂/∂param`, with no `eps` dressing variable anywhere. Same δln semantics, structurally different route.

2. **Λ factorization proved on fresh abstract symbols.** The `.wl` M2 (wl:175–182) proves `lambdaSym == (gW/omW^2)(1+I)/(1-H)` on fresh abstract symbols `gWSym, rSym, omUSym, omWSym` — not the sample-substituted primitives — with the pole guard `omUSym^2 omWSym^2 - rSym^2 != 0` in `$Assumptions` (wl:104). Strong independence; residual = 0.

3. **The `200+147π²` no-go factor is DERIVED by the `.wl`, not hard-pinned.** The `.py` hard-pins `expected_det` carrying `(200 + 147*sp.pi**2)` and asserts equality (py:215–216). The `.wl` does NOT carry the literal: it computes `detIH = FullSimplify[Det[reducedIH]]` from scratch (wl:250) and asserts only structural nonzero (wl:252) plus rank 2 (wl:253). The `.wl`'s emitted output independently produces `(... (200 + 147*Pi^2) ...)` (math.txt:66), and the SymPy output emits the same `(200 + 147*pi**2)` factor (sympy.txt:21). No surviving `251+215π²` anywhere in the current `.py`, `.wl`, or either committed `.txt` (repo grep over the unit-227 carriers returns empty). The pre-renumber paper_misalignment correction (`251+215π² → 200+147π²`) is re-confirmed holding across all carriers.

4. **Other load-bearing asserts are anchored.** The slope-law residual checks (A3/A5/A6 in `.py`; M1/M3 in `.wl`) compare three genuinely distinct constructions of `Xi_1` (direct `δN0/N0`, the `m,i,h` recombination, the specialized sample law), not an object against itself. Survivor directions are pinned eigen-directions of `red_i`/`red_h`/`red_m` nullspaces (py:225–235; wl:264–274), agreeing cross-engine to ~5e-9. The determinant nonzero check uses real positive literals (`200+147π²>0`, `98π²-25>0`), so the co-loading no-go is substantive, not vacuously passing.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Verified pure-transfer theorem:` with `I = 3/16`, `H = 25/(98*pi**2)` (sympy.txt:15–16).
- `det[(i,h)|_pure transfer] = -19*(-25 + 98*pi**2)*(200 + 147*pi**2)*(441*pi**2 + 4400)/(...)` (sympy.txt:21).
- `Stage 227 audit completed successfully.` (sympy.txt:34).

**Mathematica:** exit=0. Notable lines:
- `OK  M1 Xi1 direct N0 variation equals load-factor variation` (math.txt:18); `OK  M3 I sample equals 3/16` / `OK  M3 H sample equals 25/(98 Pi^2)` (math.txt:43,45).
- `det[(i,h)|pure transfer] = (-19*(-25 + 98*Pi^2)*(200 + 147*Pi^2)*(4400 + 441*Pi^2))/(...)` (math.txt:66), then `OK  M4 combined i=h determinant is structurally nonzero` (math.txt:68).
- `All Stage 227 Mathematica checks passed.` (math.txt:125).

**Output freshness:** The `.wl` output `.txt` mtime is 2026-06-09 19:20 — fresh from the orchestrator's just-completed close-out re-run, newer than the `.wl` source (2026-06-03). The `.py` output `.txt` mtime is 2026-06-02 17:25, which postdates the unchanged `.py` source (2026-05-11); its body matches the current source line-for-line (determinant, slope, survivor, ceiling literals all appear verbatim), so it is current. The `.wl` output refresh was label-only: its in-output cross-stage labels were previously stale (`Stage-242`/`Stage-243`) and the current `.wl` SOURCE already prints the canonical `Stage-225`/`Stage-226` (wl:154,243,244) — confirmed no surviving `Stage-242`/`Stage-243` in either committed `.txt`. No numeric change; not a defect.

## Material-change assessment

`material_change`: false. No script edits were applied; no derived result changed. The committed outputs are identical in substance to the prior committed state (only a label-only `.wl` output refresh). Downstream units are unaffected.

## Side observations (non-blocking)

- The card's `\stagefield{Verification}` reads "Mathematica audit: none yet" while a `.wl` now exists — a prose lag noted by the auditor, already tracked under the deferred numbering/output-band plan. Out of scripts-only scope; non-blocking.
- The SymPy output `.txt` (Jun 2) was not re-touched this close-out because the `.py` was unchanged; its content is verified current against the source, so no action needed.

## Verdict justification

`verified`. Unit 227 was audited `clean` with zero findings and zero directive, so there is nothing for Codex to have applied — the empty diff is correct, not a failure. I re-read both current scripts and confirmed the clean verdict holds: both exec logs exit 0 with all `OK`/passed lines; the load-bearing determinant `(200+147π²)` no-go factor is independently DERIVED by the `.wl` (`Det[reducedIH]` from scratch) rather than transliterated from the `.py` hard-pin; the `.wl` first-variation uses a distinct `param·∂/∂param` log-derivative operator with no `eps` dressing; the Λ factorization is proved on fresh abstract symbols; and no stale `251+215π²` or `Stage-242`/`Stage-243` survives anywhere. The only `.txt` refresh was label-only with no numeric change. No regressions, no material change.
