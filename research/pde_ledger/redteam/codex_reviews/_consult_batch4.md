# Claude↔Codex math-direction consult — red-team remediation BATCH 4 (stages 130, 131)

- **Date:** 2026-05-29
- **codex_session_id:** `019e7594-ddea-7820-bc0d-8f48df46e523` (read-only, ephemeral)
- **Verdict:** Codex **CONCUR on all four** (Q1–Q4), no DISPUTE. Both orchestrator-review catches (131 R1 tautology, 131 F3 sign error) confirmed.
- Raw transcript (~90 KB) deleted after extracting this record (transcript-bloat policy; prompt forbade repo-wide grep, raw stayed small).
- 125 & 126 carried NO open question (clean directives, already Codex-fixed exit 0) — not part of this consult.

## Q1 — Stage 130 R1: global strict-monotonicity certificate → **CONCUR (use the covariance certificate)**

**CONCUR.** The FKG/Chebyshev symmetrized-covariance argument is a SOUND global proof of `dg_Π/dΠ > 0` for all Π>0 (standard covariance identity; the script already forms `dg/dPi + Cov/L` at `stage130_…sympy_audit.py:21`), strictly better than the 6-point sweep (`…py:38-44`, "not a proof"). Bracket uniqueness `2/π < g_- < 1` follows from the monotone range `2/π→1` (`notes/…stage130…:84-108`). Double-integral feasibility under 600s is "reasonable but not guaranteed" → the fallback ladder (`directive:352-369`) is correct, **especially the no-sampling HALT**. Resolution: F1 stays as drafted; Codex (fix-applier) RUNS it and uses the ladder (single-integral deviation form → guarded Reduce → HALT) if the double integral stalls; NEVER revert to sampling.

## Q2 — Stage 131 R1: proposed round-trip is STILL tautological → **CONCUR, DROP Anchor-3 (option ii)**

**CONCUR.** `threshold_residual` is defined at `stage131_…sympy_audit.py:34` and compared to the same expression at `:55-62` (Mathematica mirror `…wl:62-73`) — exact X−X. The notes define group/branch/threshold/split (`notes/…stage131…:17-65`) with NO independent upstream numeric values. The proposed round-trip hardcodes BOTH the value formula and its literal inverse (`directive:123-143`), so `(Π_*·Θ_σ/L)·(L/Θ_σ) ≡ Π_*` identically — it cannot catch a misquoted group def. **Option (ii) DROP (`directive:483-485`) is cleaner** than keeping a dressed-up tautology (and unlike 118's flagged `K_q`, there is no in-stage upstream anchor to keep it). → Rewrite F1 to delete the Anchor-3 block in both engines; record the threshold identity as purely definitional/out-of-scope; Anchors 1/2 (PASSED) + R2 carry the stage.

## Q3 — Stage 131 R3: sign error in the independent Mathematica Π_* residual → **CONCUR, corrected residual is right**

**CONCUR.** Clearing `gPi − gMinus = 0` (`stage131_…sympy_audit.py:23`) gives `40π·p(2p·eᵖ+π) − 20π·gMinus·(4p²+π²)(eᵖ−1) = 0`; the directive's `directive:409` had the sign flipped (`−(1−eᵖ)…` ⇒ `+(eᵖ−1)…`, ≈6366 at Π_* not 0). Corrected residual `40*Pi*p*(2*p*Exp[p]+Pi) - 20*Pi*gMinus*(4*p^2+Pi^2)*(Exp[p]-1)` has its root at Π_*. The cleared-denominator + bracketed-`FindRoot` route is **genuinely independent** of SymPy's rational-form single-seed `nsolve` (`…py:24`), matching the IV.5 139/143/144 precedent for transcendental numerical roots. → Fix the F3 residual (code + derivation text) before codex-invoke.

## Q4 — Stage 131 R2: branch-discrimination adequacy + upper-branch value → **CONCUR**

**CONCUR.** (4a) lower membership + (4b) separation from `g_nat=1` + (4c) separation from `g_+` together exercise paper Checks-item-3 (`directive:257-292`). The singular separation is notes-grounded by `Δg_-` (`notes/…stage122…:104`). The upper-branch value is correct: the exact `±37√3` pair is boxed at `notes/…stage122…:53-57`, with `g_+^{F1} ≈ 2.79795199200529` at `:60-66`.

## Actions taken (orchestrator, post-consult)
- 130: cleared `needs_user_resolution`; added RESOLVED block (certificate sound; Codex runs + uses no-sampling fallback ladder).
- 131 F1: rewritten to DROP Anchor-3 in both engines (definitional/out-of-scope).
- 131 F3: corrected the cleared-denominator residual sign in both the code block and the derivation prose.
- None conceptual → no user escalation. All remain `material_change: false` expected (verification-surface changes only).
