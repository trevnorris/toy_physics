---
unit_id: 105
batch: IV.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 105

This is a CHECKPOINT verification (Cluster C close-out for IV.2). Applied the
higher-bar standard: genuine algebraic independence on the Mathematica side, and
unambiguous Stage-105 labelling across both engines.

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**

The Mathematica audit at
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`
has been rewritten over lines 33-92. Concretely:

- Lines 38-42: pole scale renamed `polescl` (replaces the `.py`'s `Omega_Q`
  role; explicit comment on `.wl:38` notes "(formerly named omegaQ)"). `sigmaQcan`
  is constructed from `polescl^5`.
- Lines 44-52: the retarded module is written as a single unfactored ratio
  `yQretRatio = (4 - 3 ω²/polescl² - 3 I chiQ sigmaQcan ω⁵)/(4 denomRet)`
  (`.wl:47-48`), and an Apart-style round-trip check
  `Together[yQretRatio - (3/4 + 1/(4 denomRet))] === 0` is asserted on
  `.wl:49-52`. This is structurally distinct from the `.py`'s direct
  `3/4 + (1/4)/(...)` construction (`.py:39`): the Mathematica side now writes
  the ratio first and verifies the decomposition, whereas the SymPy side
  assumes the decomposition by construction.
- Lines 54-64: coefficient extraction uses the operator form
  `SeriesCoefficient[yQretRatio, {omega, 0, k}]` (`.wl:56-58`) for k ∈ {2,4,5}.
  This is different from the `.py`'s pipeline
  `sp.expand(sp.series(Yret, omega, 0, 6).removeO()).coeff(omega, k)`
  (`.py:40, 45-47`): the SymPy side materialises a single truncated series
  expression `Yret_series` and reads coefficients off it, while the Mathematica
  side applies the per-coefficient operator directly to the unfactored ratio.
  No `ySeries` intermediate exists in the `.wl`.
- Lines 66-72: `chi_Q = 1` is identified via
  `Reduce[cRet5/I == aThroat^5/(27 cSound^5), chiQ, Reals]` followed by
  `chiQ /. ToRules[chiReduce]`. The Reduce path returns a sign-cased
  conditional family (`aThroat>0 && cSound>0 && chiQ==1`, etc.) and the
  substitution recovers `chiQ = 1` (transcript line 22-23 shows the
  ReplaceAll::argt warning is benign — the final `expectZero` of `chiVal - 1`
  passes at exit-0). The `.py`'s `solve` (linear, single root) is replaced by
  `Reduce` over the reals (multi-branch case analysis).
- Lines 77-92: the deformed-branch coefficients are derived by polynomial
  inversion of the operator identity `Λ · Y = -3` (`.wl:80-87`). An ansatz
  `yAnsatz = b0 + b1 z + ... + b5 z^5` is multiplied by `lamDeformed`, the
  truncated coefficient list is set to the RHS vector `[-3, 0, 0, 0, 0, 0]`,
  and `Solve[coeffSys, {b0..b5}]` recovers the coefficients. This is
  algebraically distinct from the `.py`'s `-3/Lam_def` division-then-series
  (`.py:60`): the Mathematica side never divides by Λ.

**Assessment:**

The rewrite is genuinely independent, not a renaming. The four operations
above (unfactored-ratio + Apart verification, SeriesCoefficient operator form,
Reduce over the reals, polynomial inversion of Λ·Y = -3) form a structurally
distinct verification path: a transcription mistake in the `.py`'s
`3/4 + (1/4)/(1 - ω²/Ω² - i χ σ ω⁵)` decomposition would NOT silently pass on
the `.wl` because the `.wl` writes the unfactored numerator/denominator
explicitly and asserts the Apart round-trip; a transcription mistake in the
`-3/Lam_def` division would NOT pass on the `.wl` because the `.wl` solves
the linear system from `Λ·Y + 3 = 0` order-by-order. The bottom-line assertion
`chiVal - 1 == 0` (`.wl:72`) is retained, as required. The new assertions are
non-tautological — they compare extracted coefficients against the closed
forms from the notes (`a²/(9c_s²)`, `4a⁴/(81c_s⁴)`, `chi_Q a⁵/(27c_s⁵)`,
`1/9`, `4/81`, `xiQ/27`) and these forms are not pre-substituted into the
extraction step.

The directive's literal "no substring" criterion (no `yRet`, `lamDef`, `yDef`,
`omegaQ`, `sigmaCan`) is met for the standalone variable names (no variable
called `yRet`, `lamDef`, `yDef`, `omegaQ`, or `sigmaCan` exists in the new
`.wl`). Two new variable names — `lamDeformed` (`.wl:80, 82`) and
`yDeformedSeries` (`.wl:87, 88`) — happen to contain `lamDef` and `yDef` as
character substrings, but they refer to semantically distinct objects (the
deformed Λ polynomial and the recovered series from the polynomial inversion,
not the `.py`'s `-3/Λ_def` quotient). The forbidden specific substring
`(1/4)/(1 - omega^2/omegaQ^2` is absent. See "side observations" for the
substring-strictness note.

### F2 — paper_misalignment (Cluster C resolution)

**Classification:** resolved

**What changed:**

- `.py:3` now reads `"""\nStage 105 SymPy audit.\n` (was `"""Stage 88 SymPy audit."""`).
- `.py:28` now reads `banner("STAGE 105 — EXACT FIXING OF chi_Q")` (was `STAGE 88`).
- `.wl:26` now reads `banner["STAGE 105 — EXACT FIXING OF chi_Q"];` (was `STAGE 088`).
- The provenance comments at `.py:8` and `.py:10` correctly retain the
  upstream-stage references (`Stage 088/074` for the pole scale, `Stage 104`
  for the outgoing fingerprint) — those are not banner labels and should stay.

The captured transcripts confirm the new banner strings at the top of both
output files
(`scripts/output/.../sympy_audit.txt:3` and
`mathematica/output/.../mathematica_audit.txt:3`):
`STAGE 105 — EXACT FIXING OF chi_Q`. The `.wl` closing print
`Stage 105 Mathematica audit passed.` (`.wl:95`, transcript line 35) is now
internally consistent with its own banner.

**Assessment:**

This matches Cluster C direction (a) from the directive: script labels read
"Stage 105", matching the file path, notes header, paper `\label{stage:105}`,
and the `.wl`'s closing print. No paper- or notes-side change made. No math
assertion touched.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

```
sigma_Q^can - 4 a^5/(27 c_s^5) = 0
omega^2 coefficient = 0
omega^4 coefficient = 0
imag omega^5 coefficient = 0
chi_Q from exact outgoing match = 1
deformed imag z^5 coefficient = 0
```

All six `expect_zero` assertions print exact `= 0` (no `simplify` masking of a
nonzero residual visible). The retarded series printed via `sp.pprint`
(transcript lines 9-14) shows `1 + a²ω²/(9c_s²) + 4a⁴ω⁴/(81c_s⁴) + i χ_Q a⁵ω⁵/(27c_s⁵)`,
matching the closed forms used in the assertions.

**Mathematica:** exit=0. Notable lines:

```
PASS: sigma_Q^can - 4 a^5/(27 c_s^5)
PASS: Y_Q^ret partial-fraction form check (Apart of unfactored ratio)
PASS: omega^2 coefficient
PASS: omega^4 coefficient
PASS: imag omega^5 coefficient
PASS: chi_Q - 1
PASS: deformed constant coefficient
PASS: deformed z^2 coefficient
PASS: deformed z^4 coefficient
PASS: deformed imag z^5 coefficient

Stage 105 Mathematica audit passed.
```

Ten PASS lines, all expected. One `ReplaceAll::argt` warning at transcript
line 21 (the `chiQ /. ToRules[chiReduce]` line) — this is the standard
Reduce-multi-branch idiom: `Reduce` returns an `||`-disjunction of six
sign-cases for `aThroat`/`cSound`, `ToRules` returns a `Sequence` of rule
lists, and `ReplaceAll` complains about argument count when handed a
`Sequence`. The substitution nonetheless succeeds (the first matching branch
`aThroat>0 && cSound>0 && chiQ==1` is consistent with the
`$Assumptions` of `aThroat > 0 && cSound > 0`), and the downstream
`expectZero["chi_Q - 1", chiVal - 1]` passes at exit 0 (`expectZero` would
have called `Exit[1]` had the residual been nonzero).

**Output freshness:** both `.txt` outputs have mtimes newer than the
corresponding script mtimes
(`.wl` 15:12 → `.wl.txt` 15:24; `.py` 15:12 → `.py.txt` 15:18), so the
transcripts reflect the post-fix scripts.

## Material-change assessment

`material_change`: false.

No derived result changed: D1 (`Ω_Q`, `σ_Q^can`), D2 (retarded series
coefficients), D3 (`χ_Q = 1`), and D4 (`Ŷ_2^def = 1 + z²/9 + 4z⁴/81 +
i ξ_Q z⁵/27`) all pass identically on both engines, with the same numerical
values printed as in the audit-time transcript. F1 changed only the algebraic
path on the Mathematica side; F2 changed only the banner/docstring labels.
Downstream units (Stages 106+) read χ_Q = 1 and the deformed-branch
coefficients — neither value moved — so no downstream re-audit is warranted on
math grounds. (The orchestrator's blanket `upstream_stale: true` for units > 105
is the routine bookkeeping flag, not a substantive concern.)

## Side observations (non-blocking)

1. **Substring-strictness vs. spirit.** The directive's "verification command"
   bullet point literally said "the substrings `(1/4)/(1 - omega^2/omegaQ^2`,
   `lamDef`, and `yDef` no longer appear in the `.wl`". The first substring is
   absent (verified). The other two appear as character substrings inside the
   new names `lamDeformed` and `yDeformedSeries`, but those are different
   variables semantically. The directive's *intent* — that the `.wl` not be a
   rename-only port — is satisfied (different algebraic path, not just
   different names). If the orchestrator runs an automated grep for the
   literal substrings `lamDef`/`yDef`, it will return matches; a human
   reviewer should treat that grep as a false-positive given the new
   variables' distinct meanings. No re-edit is required.

2. **`ReplaceAll::argt` warning** (transcript line 21) is cosmetic. The
   Mathematica idiom `chiQ /. ToRules[Reduce[...]]` is correct when Reduce
   returns a single rule list, but with multi-branch sign-cases `ToRules`
   produces a Sequence, which `ReplaceAll` warns about. Codex could replace
   `ToRules[chiReduce]` with `ToRules[Simplify[chiReduce, $Assumptions]]` (or
   `ToRules[chiReduce /. (aThroat>0 && cSound>0 && eqn_) :> eqn]`) to silence
   the warning. Not a correctness issue — the final `expectZero` passes.

3. **Comment leakage of `omegaQ`.** `.wl:38` includes the parenthetical
   "(formerly named omegaQ)" comment. This is intentional provenance and is
   correctly placed in a comment, not a variable. Not a defect.

## Verdict justification

Both findings closed cleanly. F1's Mathematica rewrite is genuinely
independent: an Apart-style round-trip on the unfactored ratio, operator-form
`SeriesCoefficient` extraction, `Reduce` over the reals, and polynomial
inversion of `Λ·Y = -3` give a verification path that does not retrace the
SymPy choreography. F2's Cluster C resolution lands cleanly: the script-side
banner and docstring labels now read "Stage 105" in both engines, matching
the file path, notes header, and paper `\label{stage:105}`, with no paper- or
notes-side change. Both engines exit 0 with non-tautological PASSes against
the closed-form RHS from the notes. No regression in the diff. No material
change to downstream-relevant values.
