---
unit_id: 181
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T02:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
iteration: 2
---

# Verification — unit 181 (iteration 2)

Iteration 1 was `needs_rework` for a single residual defect: the F1 **secondary** M2 check
(`d/dzeta Xi_1 (support-loaded route)`) was vacuous on the Mathematica side because Codex had set
`t2LoadedPert = t2DirectPert` — a ζ-free object — so `D[xi1Loaded, zeta]` was trivially 0 regardless
of physics. F1's primary M1 check and all of F2/F3/F4 were already confirmed resolved in iteration 1.
The iteration-2 delta required rebuilding `t2LoadedPert` from a genuinely ζ-bearing route via
`rTargetLoaded`, with a `FreeQ[t2LoadedPert, zeta]` guard that fails if ζ is lost.

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed (iteration 2, Mathematica only):**
The Mathematica `d/dzeta Xi_1 (support-loaded route)` block (`.wl:100-109`) no longer sets
`t2LoadedPert = t2DirectPert`. It now reconstructs a genuinely ζ-bearing loaded shape:
- `mMixPert` (wl:100), `mSuppPert` (wl:101-102), `sSupportPert` (wl:103) — `mSuppPert` and
  `sSupportPert` both carry explicit `zeta`;
- `productLoadedPert = 8*lamNormPert*(1-epsPert)/Pi^2*sSupportPert` (wl:104);
- `rTargetLoadedPert = productLoadedPert/(mMixPert + mSuppPert)` (wl:105);
- `t2LoadedPert = lamNormPert/(omegaW2*(1+s*omegaW))*(1-epsEta)/rTargetLoadedPert` (wl:106).
A hard guard `If[TrueQ[FreeQ[t2LoadedPert, zeta]], fail[...]]` is inserted at wl:107 BEFORE the
`s`-log-derivative is taken, then `xi1Loaded = D[Log[t2LoadedPert], s] /. s -> 0` (wl:108) and
`expectZero["d/dzeta Xi_1 (support-loaded route)", D[xi1Loaded, zeta]]` (wl:109). The dead
`sPar = s;` line called out in the delta is absent (grep confirms no `sPar` in the file).

**Assessment:** Correct and non-tautological.
(a) `t2LoadedPert` is built from the perturbed ζ-bearing `rTargetLoaded` route, NOT from `t2DirectPert`.
(b) The `FreeQ[t2LoadedPert, zeta]` guard is present (wl:107) and did NOT fire — the Mathematica
exec exits 0 with no "t2LoadedPert lost zeta" failure printed, so ζ is genuinely retained in
`t2LoadedPert` prior to differentiation. (c) `D[xi1Loaded, zeta]` reduces to 0 (log line 37 `= 0`,
line 38 `PASS`); because the differentiated object provably still contained ζ, this is a real
algebraic cancellation, not the iteration-1 vacuous 0. (d) The primary M1 check
`Xi_1 derived from T^2 matches defect law = 0` (wl:99, log line 36) remains genuine — it perturbs
the closed-form `t2DirectPert` in all five primitives and matches the line-81 `xi1` closed form.
The SymPy M2 check was already genuine in iteration 1 (it perturbs the ζ-bearing `T2_loaded` via
`.subs`, sympy `.py:135-140`) and is untouched (SymPy script mtime unchanged at 01:15, iteration-2
delta was Mathematica-only). Both engines now exercise a genuinely ζ-bearing object for the
support-blindness-of-Ξ₁ headline. **resolved.**

### F2 — tautological_check

**Classification:** resolved (unchanged from iteration 1)

**What changed:** `R1` is independently derived via the perturbed closed-form `R_target`:
SymPy `.py:115-120` (`R1_derived = diff(log(R_target_pert), s).subs(s,0)`, `expect_zero(R1_derived - R1)`);
Mathematica `.wl:86-90` (`rTargetPert`, `r1Derived`, `expectZero["R_1 derived from R_target matches
closed form", r1Derived - r1]`). The original line-91/97 `selected-branch identity` check is retained
but is now backed by an anchored `R1`.

**Assessment:** Not regressed by iteration 2. `R_1 derived from R_target matches closed form = 0`
PASS in both transcripts (sympy log 22, math log 29-30); `selected-branch identity = 0` PASS
(sympy log 23, math log 31-32). The selected-branch identity is now a genuine cross-check between two
independently-anchored drifts, not vacuous. **resolved.**

### F3 — mathematica_transliteration

**Classification:** resolved (unchanged from iteration 1)

**What changed:** The Mathematica loaded route is reached through the `sSupport` cancellation
(`loadMassFromSupport = mMix*sSupport`, `rTargetLoaded = productLoaded/loadMassFromSupport`,
wl:42-43) rather than the SymPy `mMix + mSupp` term-by-term grouping; the F1/F2 perturbation
checks are inherently independent of the SymPy choreography; and the `spoiled` negative control was
ported (wl:53-64), producing a nonzero ζ-bearing drift.

**Assessment:** Not regressed by iteration 2. The iteration-2 loaded-perturbation block (wl:100-105)
adds further independent structure rather than re-mirroring SymPy. `spoiled d/dzeta ln R_target`
prints a nonzero rational in ζ in both transcripts (sympy log 16, math log 22) and the
`If[TrueQ[spoiledDrift === 0], fail[...]]` guard confirms the support-blindness check can fail.
The `.wl` no longer reproduces the SymPy `mMix/mSupp/mTr/productLoaded` grouping verbatim.
**resolved.**

### F4 — cosmetic banner mislabel

**Classification:** resolved (unchanged from iteration 1)

**What changed:** Both banners changed from `STAGE 164` to `STAGE 181` (sympy `.py:33`,
math `.wl:26`).

**Assessment:** Both transcripts now print `STAGE 181 — COHERENT TRACKING-BRANCH DEFECT LAW`
(sympy log 8, math log 8). **resolved.**

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Xi_1 derived from T^2 matches defect law = 0` (line 26) — primary M1, genuine.
- `d/dzeta Xi_1 (support-loaded route) = 0` (line 27) — M2 from ζ-bearing `T2_loaded`, genuine.
- `R_1 derived from R_target matches closed form = 0` (line 22), `selected-branch identity = 0` (line 23).

**Mathematica:** exit=0. Notable lines:
- `d/dzeta Xi_1 (support-loaded route) = 0` / `PASS` (lines 37-38) — now from the ζ-bearing
  `rTargetLoadedPert` route; the `FreeQ[t2LoadedPert, zeta]` guard did NOT fire (no
  "t2LoadedPert lost zeta" line anywhere in the transcript), so ζ was genuinely present before
  differentiation. This is the iteration-2 fix landing.
- `Xi_1 derived from T^2 matches defect law = 0` / `PASS` (lines 35-36).
- `R_1 derived from R_target matches closed form = 0` / `PASS` (lines 29-30);
  `selected-branch identity = 0` / `PASS` (lines 31-32).
- `spoiled d/dzeta ln R_target = (...)` nonzero in ζ (line 22) — negative control fires.

**Output freshness:** confirmed. The committed `.txt` outputs are newer than their scripts:
- math `.wl` 02:02:24 → math output `.txt` 02:07:03 (newer; output regenerated post-fix).
- sympy `.py` 01:15:36 → sympy output `.txt` 02:06:52 (newer). The SymPy script was not re-edited
  in iteration 2 (delta Mathematica-only), but the orchestrator re-ran both engines, refreshing the
  SymPy output too. The committed math `.txt` contents match the captured exec log (verified by grep
  of the F1/F2/F3/banner lines) — i.e. the saved output reflects the iter-2 FreeQ-guarded state with
  no guard failure.

## Material-change assessment

`material_change`: false. All edits add or de-tautologize verification checks; no closed-form target
(`Xi1`, `R1`, `Theta1`, `eps1`, the transfer-shape constant) was redefined. The printed `Xi_1`,
`R_1`, `Theta_1` forms are identical to iteration 1 / the original audit. No downstream unit can
depend on a changed derived result.

## Side observations (non-blocking)

- `epsPert` is defined twice in the Mathematica script (wl:87 and wl:96) and `eps_pert` twice in
  SymPy (`.py:116` and `.py:127`) — harmless redundant re-definitions to identical RHS, a minor
  tidiness artifact of inserting the F1 and F2 blocks separately. Not a correctness issue.
- The Mathematica `FreeQ` guard (wl:107) and the SymPy side give the support-blindness-of-Ξ₁ check
  genuinely different shapes (Mathematica reconstructs the loaded packet from `mMixPert`/`mSuppPert`;
  SymPy `.subs` into the precomputed `T2_loaded`), which incidentally strengthens the F3
  independence story for this specific check.

## Verdict justification

Iteration 2 closed the sole residual defect from iteration 1. The Mathematica F1 secondary check is
now genuine: `t2LoadedPert` is built from the ζ-bearing `rTargetLoadedPert` route (NOT
`t2DirectPert`), a `FreeQ[t2LoadedPert, zeta]` guard confirms ζ is present before differentiation and
did not fire, and `D[xi1Loaded, zeta] = 0` is therefore a real cancellation. F2/F3/F4 remain
resolved and were not regressed; SymPy was untouched and still passes. Both engines exit 0 and the
committed outputs are fresh and match the logs. All four findings `resolved`; verdict `verified`.
