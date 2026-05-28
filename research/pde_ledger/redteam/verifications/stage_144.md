---
unit_id: 144
batch: IV.5
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 144

## Per-finding outcomes

### F1 — paper_misalignment (banner mislabelling)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:16` now reads `banner("STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH")`.
- Same file `:70` (post-insertion of F2's asserts; original line 46 in the directive) now reads `banner("STAGE 144 LEDGER")`.
- `mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:27` now reads `banner["STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH"];`.
- Same file `:73` (shifted from line 59 by F2 insertions) reads `banner["STAGE 144 LEDGER"];`.
- Both refreshed transcripts header `STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH` and `STAGE 144 LEDGER`. No `STAGE 127` substring remains in either script source or either output transcript.

**Assessment:**
Edit matches the directive verbatim. The renumber is purely cosmetic and does not affect derived results. Orchestrator note confirmed (Cluster A mass-renumber pass). Per orchestrator policy, treated as resolved after confirming STAGE 144 banner in both transcripts.

### F2 — insufficient_verification (numerical-target assertions)

**Classification:** resolved

**What changed:**
- SymPy script: lines 40-43 now compute and print `Sigma0_star` and `Sigma0_match`. Lines 49-68 add the requested numerical-target guards: `gplus > 1`, `2/pi < gminus < 1`, and five 1e-12-tolerance drift checks for `Pi_*`, `That_*`, `Sigma0_*`, `Pi_match`, `That_match`. A consolidated `PASS:` print follows on line 68.
- Mathematica script: lines 59-62 add `sigma0Star`/`sigma0Match` computation and prints. Lines 64-71 add the upper-branch, bracket, and five drift guards using the existing `pass`/`fail` helpers. Each emits a discrete `PASS: …` line.
- Both transcripts show the expected outputs: SymPy emits the consolidated `PASS: numerical-target assertions …` line and Mathematica emits all seven discrete `PASS:` lines (`upper branch`, `lower branch bracket`, `Pi_*`, `That(Pi_*)`, `Sigma0(Pi_*)`, `Pi_match`, `That(Pi_match)`).

**Assessment:**
The added assertions are non-tautological — they compare numerically derived quantities (from `nsolve`/`FindRoot`) against the notes-anchored decimal targets at tolerance `1e-12`, well looser than the 30-digit precision both engines use but tight enough to catch a sign flip, missing factor, or rational-coefficient error in `r`, `g_Π`, `s_Q`, `r_Q`, `Σ_0`, or `T̂`. The targets `1.50882951349316`, `0.901484054174205`, `1.80594111095636`, `1.90848600654854`, `1.01132972803599` exactly match the auditor's directive values. The SymPy `Sigma0(Pi_*) = 1.80594111095635380721796724713` and Mathematica `Sigma0(Pi_*) = 1.8059411109563538723736729091995054268193435351723899810506` both pass the `1e-12` tolerance versus the target `1.80594111095636`. No tolerance relaxation; no tautology; banner and assertion structure follow the directive exactly.

### F3 — insufficient_verification (paper-card checks (i) and (ii))

**Classification:** resolved (policy-accepted via paper-card downgrade)

**What changed:**
Per orchestrator note: paper-card edit at `paper/stages/stage_144.tex` lines 21-25 rewrote items (i) and (ii) as carry-forward citations of `\ref{stage:135}` and `\ref{stage:140}` (Cluster C downgrade). No script change is in scope or required.

**Assessment:**
Verifier policy disallows reading paper.tex, so the paper-side change is accepted on orchestrator authority. The script-side claim list is now consistent with the downgraded paper-card scope (numerical recording + F2's value targets), so no further script work is needed.

### F4 — mathematica_transliteration

**Classification:** resolved (policy-accepted; option a/c from directive)

**What changed:**
Per orchestrator note, F4 was accepted as written per user gate. No script change required. The Mathematica script remains a transliteration of the SymPy script in its symbolic forms for `r, gMinus, gPlus, gPi, sQ, rQ, sigma0, tHat`, but the project user has accepted that the second-engine value here is the independent numerical root-finder (Mathematica's 30-digit `FindRoot` vs. SymPy's `nsolve`).

**Assessment:**
Accepted on user gate. No verification action required.

## Exec log assessment

**SymPy:** exit=0 (no traceback; final ledger printed in full; `PASS: numerical-target assertions (upper/lower branches, Pi_*, Sigma0_*, That_*, Pi_match, That_match)` emitted at line 18 of the transcript). Notable lines:

```
STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH
Sigma0(Pi_*)  = 1.80594111095635380721796724713
PASS: numerical-target assertions (upper/lower branches, Pi_*, Sigma0_*, That_*, Pi_match, That_match)
STAGE 144 LEDGER
```

**Mathematica:** exit=0 (final `Stage 144 Mathematica audit passed.` line emitted; no `FAIL:` lines; banner says STAGE 144). Notable lines:

```
STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH
PASS: upper branch g_+^F1 > 1
PASS: lower branch bracket 2/pi < g_-^F1 < 1
PASS: Pi_* matches notes target
PASS: That(Pi_*) matches notes target
PASS: Sigma0(Pi_*) matches notes target
PASS: Pi_match matches notes target
PASS: That(Pi_match) matches notes target
Stage 144 Mathematica audit passed.
```

**Output freshness:** confirmed. mtimes:
- sympy script `19:50:02` → output `19:51:27` (output newer).
- mathematica script `19:50:13` → output `19:53:39` (output newer).

Both outputs were regenerated after the script edits.

## Material-change assessment

`material_change`: false.

The applied edits add cosmetic banner renames (F1) and additional non-tautological assertions (F2) but do not alter any numerical or symbolic derived quantity. The printed values for `g_∓^F1`, `Π_*`, `T̂(Π_*)`, `Σ_0(Π_*)`, `Π_match`, `T̂(Π_match)` are identical pre- and post-fix (they only changed in that `Σ_0(Π_*)` and `Σ_0(Π_match)` are now printed in addition). Downstream stages 146-153 that import `Π_*` and `T̂_m(Π_*)` are unaffected; in fact, F2's new asserts now actively guard the constants those stages inherit.

## Side observations (non-blocking)

- SymPy and Mathematica still diverge from each other at the 17th significant digit on `Π_*`, `T̂(Π_*)`, `Π_match`, `T̂(Π_match)` (e.g., SymPy `Π_* = 1.50882951349315552747…` vs. Mathematica `1.50882951349315558300…`). This is the same `nsolve`-vs-`FindRoot` precision mismatch the auditor already noted; both still pass the `1e-12` tolerance against the notes target, so it is not a verification failure. Per F4 policy-accept, no action needed.
- In the SymPy script, the `Sigma0_star`/`Sigma0_match` computations are inserted between the existing `Pi_*` print block and the legacy ordering assertion, so the transcript prints `Sigma0(Pi_*)` after the `That(Pi_*)` line rather than directly under the other `Pi_*` block. The directive permitted this placement; flagging only for surface awareness.

## Verdict justification

`verified`. All four findings closed: F1 (banner) confirmed in both scripts and both refreshed transcripts; F2 (numerical-target asserts) implemented exactly as the directive specified in both engines with non-tautological 1e-12 drift guards and a new `Σ_0` print, and both exec logs exit 0 with the expected PASS lines; F3 resolved by paper-side downgrade (accepted on orchestrator authority per verifier policy); F4 accepted as written per user gate. Exec logs are fresh (outputs newer than scripts). No regressions in the diff. `material_change: false` — derived results are bit-identical to pre-fix; only assertion coverage and cosmetic labelling changed.
