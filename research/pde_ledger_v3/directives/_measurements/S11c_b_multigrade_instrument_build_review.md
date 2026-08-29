# S11c-b multigrade instrument — BUILD review, two legs, consolidated (SOUND; one guard folded)

**Artifact:** `scripts/S11c_b_background_multigrade.py` (Codex-written, build log
`~/.s11_build/S11c_b_multigrade_build_codex.log`, 184,077 tok — real build). **Legs (Codex-written ⇒ fresh
Claude agent + Grok):** prompt `directives/_legs/S11c_b_multigrade_instrument_build_review.md`.
**Raw:** Grok `~/.s11_build/S11c_b_multigrade_buildreview_grok.txt`; Claude-agent scratch
`/tmp/.../scratchpad/s11cb_review/` (`independent_derivation.py/.out`, `mg_abl{1..4}_*.py`, `abl{1..4}.out`,
`abl1_micro.*`, `instrument_run.out`) — the agent's derivation/ablation artifacts; Grok's under
`~/.s11_build/multigrade_review/`. **Both legs: SOUND, NO BLOCKING FINDINGS.**

## CRUX — the emitted grades are the TRUE grades (both legs, independent routes)
Both legs re-extracted the `(eta_bg, sigma_W)` coefficients by a route the instrument does **not** use — Grok:
`sp.Poly.coeff_monomial` (polynomial leaves) + denom-clearing/geometric convolution of `(1+eta_bg·w1)^(-p)`
(rational coupling); the Claude agent: `sp.Poly` + an exact reciprocal/Neumann series of `1/den(eta_bg)`.
Neither is the instrument's `.diff().subs()/factorial`. Both diffed their coefficients against the emitted
`MULTIGRADE_*`:
- Grok: `TOTAL_DIFF_CELLS=0`, `ALL_CELLS_AGREE` for ADVECTIVE×4, ADMISSIBILITY-THETA×4, and a rational
  COUPLING leaf (`n_diff=0`); nonzero grades `[(0,1),(1,1),(2,1),(3,1),(4,1)]`, exact remainder leading
  `{(5,1)}` (content beyond `N=4` retained, not dropped).
- Agent: `ALL_INDEPENDENT_CHECKS_MATCH: True`, `coeff_mismatch=0` for ADVECTIVE, ADMISSIBILITY-THETA, COUPLING
  (rational), and all four KINETIC container leaves; `residual == layer A_minus_B: True` for every family,
  including the KINETIC per-label `leaf_A_minus_B`.
⇒ (i) the emitted grades are the true grades of the actual aligned operands; (ii) the graded residual is the
same object the layer routed as `A_minus_B`, keyed for KINETIC by the `_kinetic_pairs` semantic labels.
Full instrument run: exit 0, empty stderr, `EMITTED_CASES = Integer(20)` (Grok timed it at 470.24 s).

## Guards are genuine under ablation (both legs; FORM ablations mandatory)
- **FORM (swap `eta_bg`↔`sigma_W`; feed B into A; cyclically rotate kinetic leaves among the same labels):**
  grade structure moves and `WINDOW_CLEAN`/`GRADE_DIFFERENCE` go NONZERO — never byte-identical. The kinetic
  same-label `GRADE_DIFFERENCE` is non-vacuous (rotated `u_2_tt` shows up under the `[0]` label and is caught).
- **Remainder-hiding (`R:=leaf`, coeffs `0`):** `WINDOW_CLEAN` REJECTS it — content cannot hide in `R`.
- **Coefficient-formula (drop `1/(a!b!)`):** `WINDOW_CLEAN(a,b)` NONZERO on the coupling leaf (ratios
  `1,1,2,6,24`). ⚠ Agent's TEST-COVERAGE caveat: this ablation is a NO-OP on ADVECTIVE (all grades have
  `a,b∈{0,1}` ⇒ `a!b!=1`); it is only observable on a case carrying a grade with `a!b!>1`, i.e. a coupling
  case — which both legs exercised.
- Reuse/leak confirmed: operands only via `A.transform(bridge_a=True,bridge_d=False)` + `A._bridge_d` once;
  residual via `_arithmetic_residual` / `_kinetic_pairs`+`_kinetic_residual`; no `_classify_case`/Euler/
  divergence/`series`; no energy/operator/kernel rebuild; module+transcript header points at the committed
  layer/comparator/engine and the two real transcripts. No expected grade/engine verdict/coefficient target.

## The one finding — FOLDED
**`RECONSTRUCTION` was tautological (both legs independently).** `remainder := _normalise(leaf − window)`
(L193), so `reconstruction := _normalise(leaf − window − remainder) ≡ 0` for ANY coefficients — it stayed
exactly zero under all four corruptions and tests only that `_normalise` preserves value (a parser
round-trip). It is not a mis-measurement risk (the genuine coverage is `WINDOW_CLEAN` + `GRADE_DIFFERENCE`,
both shown live), but it is guard-theatre. FOLD (rule 12, build-skill "no tautological residual"): dropped the
`RECONSTRUCTION` field, its computation, `_reconstruction_association`, and its emission; left a comment naming
`WINDOW_CLEAN` as the genuine check. The coefficients and genuine guards are untouched, so the measurement is
unchanged — no re-leg (the removed guard was decorative). `py_compile` OK; committed instrument re-run to
capture the step-1 measurement (`~/.s11_build/S11c_b_multigrade_run.out`).

## Non-blocking (carry to step record)
- `_remainder_leading_grades` caps its scan at `2N` per axis (L164); a remainder whose first nonzero grade sat
  beyond `2N` in one axis would summarise as `EmptySet`. NOT triggered by any of the twenty cases (A/RESIDUAL
  leading `(5,1)`, B `EmptySet`) and it gates nothing — printed diagnostic only.
- `_coefficient` raises if a coefficient retains a bookkeeper (L141–144) — a fail-loud guard before emission;
  correct behaviour, hides nothing.

## Verdict
Instrument SOUND: it faithfully measures the `(eta_bg,sigma_W)` multigrade of the twenty aligned operands;
grades are the true grades (two independent re-derivations, zero mismatch), residuals equal the layer's routed
`A_minus_B`, kinetic leaves aligned by semantic label, genuine guards live. One tautological guard folded out.
NEXT (step 2): the per-case §2a/§3a/§3d adjudication reads the emitted grade spectra — which engine's
retention is spec-correct — off the committed instrument's stdout.
