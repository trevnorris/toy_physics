---
unit_id: 166
batch: V.2
created_at: 2026-05-29T00:00:00-06:00
supersedes: 2026-05-28T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T22:46:14-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 166 (rewrite encoding Codex review of stage 166 + consult batch 7 Q5)

## What this rewrite does

The prior directive (now superseded) added two findings (F1 general-form asserts; F2
banner fix + an independent matrix-inverse cross-check block). A Codex read-only review
(`redteam/codex_reviews/stage_166.md`, findings_count = 1) confirmed the bulk of that
work is sound (see "Already applied" below), but found that the **matrix round-trip
assertion the prior directive added is tautological** (review R1, `tautological_check`).
This rewrite encodes exactly that one finding: replace the vacuous round-trip with a
genuine forward-transcription check. The resolution is the Claude+Codex consult
(`redteam/codex_reviews/_consult_batch7.md` Q5, DECIDED — REPLACE, do not keep).

Apply the finding below. After applying, append an `## Applied: F1` block under it with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none"). If the required
change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a
question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond the lines named.
Do NOT touch paper.tex, notes/, or any prose document. After editing, RUN the affected
script (`math -script <path>`) and iterate until it exits 0 with all in-file checks
passing, within `timeout 600` (a timeout is a failure — reformulate the math, never raise
the cap). Getting the script to run cleanly is your job; the orchestrator independently
re-runs afterward.

---

## Already applied (prior directive F1/F2 — confirmed PASS by review, no-op here)

The prior directive's F1 (the four `*** general` inversion asserts: `drho general`,
`da general`, `dcs general`, `dZ general` — live at wl:48-51, sympy:53-56) and F2 (the
`"STAGE 149" → "STAGE 166"` banner fix at wl:26 / py:33, plus the
`Independent matrix-inverse cross-check` block: `Mmat`, `inv`, `solVec`, and the four
`matrix drho/da/dcs/dZ` target asserts at wl:59-72) are **already applied and confirmed
PASS** by `redteam/codex_reviews/stage_166.md` (verdict table rows 1-12 all PASS). They
are NOT re-issued as findings here. The ONLY defect carried by that prior work is the
round-trip line the prior F2 introduced (wl:73-76), which review row 13 flags — that is
Finding 1 below.

---

## F1 — tautological_check (matrix round-trip is vacuous; replace with forward transcription)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl:73-76`

**Live anchor — current code (the lines to REPLACE):**
```wolfram
(* Round-trip: forward map of the matrix solution recovers the observables.   *)
(* Sum-of-squares scalarization (zero iff every component residual is zero;    *)
(* expectZero tests res === 0, which is False for a length-4 list).            *)
expectZero["matrix round-trip", Total[(Mmat . solVec - {dTheta, dKs, dKq, dP})^2]];
```
with `solVec = Inverse[Mmat] . {dTheta, dKs, dKq, dP}` defined at wl:67-68 (`inv = Inverse[Mmat]`; `solVec = inv . {dTheta, dKs, dKq, dP}`).

**Why it is tautological:** because `solVec = Inverse[Mmat] . v`, the round-trip computes
`Mmat . Inverse[Mmat] . v - v`, which is identically `0` for ANY invertible `Mmat`. A
wrong forward-map coefficient transcribed into `Mmat` (as long as the matrix stays
invertible) leaves this residual at zero — so the assertion canNOT fail on a transcription
error. The accompanying comment's claim that it "confirms `Mmat` was transcribed correctly
from eq1..eq4" is therefore FALSE: it only checks internal inverse self-consistency
(review R1, verdict-table row 13 = FINDING).

**Required change (consult batch 7 Q5, DECIDED — REPLACE, do not keep the misleading
line):** delete the round-trip assertion and its three-line comment, and in its place
assert a FORWARD-transcription check that does NOT reuse `Inverse[Mmat]` or `solVec`.
Insert exactly:
```wolfram
(* Forward-transcription check: confirm Mmat encodes eq1..eq4 by mapping the symbolic *)
(* drift vector forward and comparing to the HAND-TYPED boxed branch laws (notes      *)
(* section 1). fwdLaws is typed from the boxed laws, NOT built from Mmat/Solve/Inverse, *)
(* so a wrong Mmat coefficient yields a nonzero residual and FAILS. Total[(...)^2]      *)
(* scalarizes the length-4 residual (expectZero tests res === 0, False for a list).     *)
fwdLaws = {2*drho, drho + 2*da, dZ + 2*dcs - 2*da, 5*(dcs - da)};
expectZero["matrix forward-transcription", Total[(Mmat . {drho, da, dcs, dZ} - fwdLaws)^2]];
```

**Branch-law / Mmat confirmation (from consult Q5 — your `fwdLaws` RHS must match these
boxed laws exactly):** in column order `(drho, da, dcs, dZ)` the boxed laws of notes §1
(around lines 54-90) are
- eq1: `δlnΘ_w = 2 δlnρ_w` ⇒ `2*drho`
- eq2: `δlnK_s = δlnρ_w + 2 δln a` ⇒ `drho + 2*da`
- eq3: `δlnK_q = δlnZ_q + 2 δln c_s − 2 δln a` ⇒ `dZ + 2*dcs - 2*da`
- eq4: `δlnP_0 = 5(δln c_s − δln a)` ⇒ `5*(dcs - da)`

`Mmat` rows 3 and 4 (wl:64-65) are `{0, -2, 2, 1}` and `{0, -5, 5, 0}`; in column order
`(drho, da, dcs, dZ)` these correctly encode `dKq = dZ + 2 dcs − 2 da` and
`dP = 5(dcs − da)` — confirmed correct in the consult. (Rows 1,2 are `{2,0,0,0}` and
`{1,2,0,0}`, encoding eq1, eq2.) The forward-transcription residual is therefore
identically zero ⇒ PASS for the correct, live `Mmat`.

**Anti-tautology guard (state explicitly, must hold):**
- The check MUST use `Mmat . {drho, da, dcs, dZ}` — the symbolic *drift* vector, the same
  primitive symbols the matrix multiplies — NOT `Mmat . solVec` and NOT anything built
  from `Inverse[Mmat]`. Reusing the inverse re-creates the `X − X` identity-cancellation
  tautology that this finding is removing.
- The RHS `fwdLaws` must be **literally hand-typed** from the notes §1 boxed laws, NOT
  produced by any `Mmat`/`Solve`/`Inverse` operation. This is what gives the check teeth:
  a wrong `Mmat` coefficient (e.g. a row-3 sign flip or a row-4 `5 → 4`) makes the residual
  nonzero and the assertion FAILS.
- Do NOT "fix" a failing `matrix forward-transcription` check by editing `fwdLaws` to match
  a (possibly wrong) `Mmat`, or vice versa. If it fails, STOP and append `## Blocked: F1`
  with the residual — a failure here is a genuine transcription-error signal.

**No-fabricated-literal guard:** the `fwdLaws` coefficients (`2`; `1, 2`; `−2, 2, 1`;
`−5, 5`) are the documented branch-law coefficients from notes §1, not tuned numbers. Do
not introduce any numeric tolerance — this is an exact symbolic residual and must be
exactly `0`.

**SymPy is not implicated:** the SymPy script
(`scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py`) has NO
matrix route at all (only the `Solve`-based general/forward/bundle/frozen blocks). It is
NOT touched by this finding — leave it untouched.

**Verification command:**
Verifier runs `redteam exec-mathematica 166` (and `redteam exec-sympy 166` for exit 0).
Confirm:
- a new `matrix forward-transcription = 0` / `PASS` line appears in the Mathematica
  transcript;
- the old vacuous `matrix round-trip` line is GONE (replaced, not supplemented);
- the Mathematica script exits 0 (the prior run reported 19 PASS / 0 FAIL; the replaced
  line is still one assertion, so the PASS count is unchanged at 19);
- the SymPy script still exits 0 and is unchanged;
- the verifier should reason that a perturbed `Mmat` coefficient (a wrong row-3 or row-4
  entry) WOULD now make `matrix forward-transcription` nonzero and FAIL — i.e. the check
  can fail, unlike the round-trip it replaces.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl`
- summary: Replaced the vacuous matrix round-trip assertion with the hand-typed forward-transcription check from the boxed branch laws.
- deviation: none

---

## RESOLVED (consult batch 7)

`redteam/codex_reviews/_consult_batch7.md` Q5 (166-R1 matrix round-trip), Date 2026-05-29:
**CONCUR** — replace wl:73-76 (`Total[(Mmat . solVec − v)^2]`, vacuous since
`solVec = Inverse[Mmat] . v`) with a HAND-TYPED forward-transcription check
`Total[(Mmat . {drho, da, dcs, dZ} − fwdLaws)^2] == 0`, where
`fwdLaws = {2 drho, drho + 2 da, dZ + 2 dcs − 2 da, 5(dcs − da)}` from notes §1 boxed laws
(NOT built from `Mmat`/`Solve`/`Inverse`). Confirmed `Mmat` rows 3,4 = `{0,−2,2,1}`,
`{0,−5,5,0}` in column order `(drho, da, dcs, dZ)` (wl:61-66) are correct. **REPLACE** the
round-trip (do not keep the misleading line). This is a how-it's-checked fix; it does NOT
change Stage 166's published (paper-card) conceptual claim. Resolved Claude+Codex; NO user
escalation needed (`needs_user_resolution: false`).
