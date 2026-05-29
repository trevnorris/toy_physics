# Consult — batch 2 paper_misalignments (105 F2, 112 F1, 109 F3, 106 F1)

**Date:** 2026-05-29
**Mode:** Claude (orchestrator) + Explore evidence agents + Codex read-only consult
**codex_session_id:** 019e748e-9d0e-72c1-baff-ebdafb90f544
**Delegation basis:** user delegated `paper_misalignment` *direction* calls to "Claude + Codex agree" unless a fix changes the CONCEPTUAL nature. Codex verdict: **no item is conceptual** — two are label bookkeeping, two are coverage/citation.

## Verdicts

- **R1 — 105 F2 (label):** **CONCUR.** Stale internal labels "Stage 88" → canonical **"Stage 105"** (file path / `\label{stage:105}` / notes / the `.wl`'s own closing print). Direction (a). Script-only (`.py:3`, `.py:28`, `.wl:26`); no paper change. Rejecting (b) display-number "Stage 122" in scripts (contradicts the IV.4/IV.5 canonical-internal-number convention) and (c) keep-stale.
- **R2 — 112 F1 (label):** **CONCUR.** Same convention: "Stage 95" → canonical **"Stage 112"**; the *mathematics already matches* the `stage:112` card, only ID strings disagreed. Direction (a). Script-only (`.py:3`, `.py:54`, `.wl:26`; `.wl:70` already correct); no paper change.
- **R3 — 109 F3 (coverage):** **CONCUR.** Stage 109 (linearized branch-selection setup) should NOT absorb secondary checks owned elsewhere. **No 109 script change.** Paper-card cross-reference the three secondary Checks to their real owners: **(a) pure scale/argument → 108**, **(b) Robin → 110 + mixed-pole no-go → 111**, **(c) compensated even+odd → 112**. (Direction (c). The separate 109 F1/F2 *tautological_check* de-tautologizations still ride the fix loop.)
- **R4 — 106 F1 (coverage):** **DISPUTE → refined.** CONCUR that **no 106 script / Hankel re-derivation is needed** (chi_Q=1 is a legitimate carry-in, not a hidden assumption). BUT the citation must be corrected: **Stage 104 proves the outgoing l=2 DtN fingerprint but introduces no `chi_Q` symbol** — **chi_Q=1 is FIXED at Stage 105** ("EXACT FIXING OF chi_Q", inheriting 104's fingerprint). Therefore:
  - item (ii) [higher odd at O(ω⁷)] → **102** (`tauQ irrelevance at ω⁵`; `tauQ coeff at ω⁷ − 1/4`)
  - item (iii) [DtN fingerprint] → **104** (`Yhat_2^out` coeffs 1/9, 4/81, i/27); **chi_Q=1 → 105**
  - **Codex: "R4 is a citation/ownership correction, not a physics change."**

## Actions taken from this consult

- **105/112 labels:** direction (a) confirmed → rides each stage's fix loop (Codex updates the three label strings). `needs_user_resolution` cleared.
- **109 F3:** no 109 script change → `needs_user_resolution` cleared; card cross-ref (a→108, b→110/111, c→112) logged to PAPER_CLEANUP_TRACKER.
- **106 F1:** no NEW verification in 106 → `needs_user_resolution` cleared; card cross-ref (ii→102, iii→104+105) logged to PAPER_CLEANUP_TRACKER. **Plus a script-doc accuracy fix that now rides the 106 fix loop:** correct the SymPy docstring + `.wl` comment that attribute chi_Q=1 to "Stage 104" — chi_Q=1 is fixed at **Stage 105** (from Stage 104's fingerprint).

## Evidence anchors (read-only; quantities derived, not hardcoded)

- 106 item (ii): `scripts/...stage102...py:51-52` `expect_zero("D1: tauQ irrelevance at omega^5", tau5)`; `expect_zero("D2: tauQ coefficient at omega^7 - 1/4", tau7 - 1/4)`
- 106 item (iii) fingerprint: `scripts/...stage104...py:37-42` static DtN slot + `Y z^2=1/9`, `Y z^4=4/81`, `Y imag z^5=1/27` (104 introduces no `chi_Q` symbol)
- 106 chi_Q=1 fix: `scripts/...stage105...py` solve outgoing match + `chi_Q - 1` assertion; `.wl:69-72` Reduce-based `chi_Q - 1` (105 docstring: "coefficient inherited from the Stage 104 exact fingerprint")
- 109 (a) scale/argument: `...stage108...py:33` `pure scale invariance`; `...py:45-49` beta=±1, `chi_arg(beta=1) - 1`
- 109 (b) Robin: `...stage110...py:29` `chi_R - 3/(3-rho)`; mixed-pole: `...stage111...py:32-33` `kappa_match + 1/9`, `sigma_match`
- 109 (c) compensated: `...stage112...py:31-38` canonical-even solve branches; `...py:49` `chi_B(gamma=1/9) - 1`
