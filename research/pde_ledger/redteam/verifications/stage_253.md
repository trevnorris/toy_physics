---
unit_id: 253
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T14:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 253

Stage 253 is the FINAL/CHECKPOINT capstone (physical calibration + material-threshold
companion). Higher bar applied to F2's dual-engine independence. All four findings resolved.

## Per-finding outcomes

### F1 — paper_misalignment (USER-RESOLVED, NOTES-ONLY)

**Classification:** resolved

**What changed:**
`notes/stages/...stage253..._sympy_audit.md:274` `\frac{187.23361317}{\zeta_{\rm ep}}.` →
`\frac{119.23361317}{\zeta_{\rm ep}}.`; and `:419` ` 10.95423247\,` → ` 10.95423248\,`.
The git working-tree diff of the notes file shows EXACTLY these two lines and nothing else.

**Assessment:**
Correct-to-script direction applied per the user-approved resolution. Verified scope is exactly
right:
- Card `paper/stages/stage_253.tex` is UNTOUCHED (`git status` clean; `git diff --stat` empty).
  It is 140 lines and contains NEITHER `187.23361317` NOR `119.23361317` NOR either a_int
  coefficient — confirming the original audit's `stage_253.tex:274/:419` citations were
  misattributions; both disputed values live only in the notes. (Card carries the symbolic
  product form Eq. app-part08-stage253-ep-product, not the benchmark decimal.)
- Script/.wl literals UNCHANGED: py:192 still asserts `119.23361317476524`; wl:280 still
  `119.23361317476524`; py:196 / wl:284 still `10.95423248`. Notes now agree with both engines.
- 119.2336… is corroborated independently: notes adjacent values give
  65.45193926/0.5489386551 = 119.2336…, matching both engines (`187` was a leading-digit typo
  sharing trailing `…23361317`). a_int `…248` = 4·K_turn = 4·2.73855812 (1-ulp notes typo fixed).

### F2 — mathematica_transliteration (CHECKPOINT higher bar)

**Classification:** resolved

**What changed:**
The `.wl` was re-authored from a 1:1 port to a genuinely non-mirrored route while preserving
every final `expectZero`/`expectNear` claim and all benchmark numbers:
- §III harmonic geometry: `V[rad_] := (1/2) kEff rad^2`, then `radialLogSlope = D[Log[V[r]], r]`
  (wl:157-161) derives `2/r` natively; `chiLambdaLattice = lambdaPhys (radialLogSlope /. r ->
  rTurnPhys)` (wl:162-165). The `.py` writes the pre-cancelled `2*lambda_phys/r_turn_phys`
  directly (py:100). True derivation, not transcription.
- `forceMatch = EStar (lambdaRef/lambdaPhys) VprimeTurnAbs == kEff rTurnPhys`; `kEffReq = kEff /.
  First[Solve[forceMatch, kEff]]` (wl:166-170). The `.py` writes the closed form directly (py:101).
- §I threshold via `Solve[gammaLatTurnPhys == zetaEp lambdaEpOmegaD, lambdaEpOmegaD]` (wl:96-101)
  vs the `.py`'s direct division `gamma_lat_turn_phys/zeta_ep` (py:53).
- §IV Korringa via `Solve[KCorr == TCeiling tCrossPhys, TCeiling]` (wl:211-214) vs `.py`'s
  `Kcorr/t_cross_phys` (py:126).
- Additional Solve-based substitutions `vPrimeFromKTurn` (wl:175-178) and `lambdaFromAInt`
  (wl:183-186) replace the `.py`'s literal `subs`. §V `PiEp` regrouped as
  `zetaEp lambdaEpOmegaD/gammaLatTurnPhys` (wl:226) vs the `.py`'s
  `Upsilon_lat*zeta_ep*lambda_ep_omega_D*t_star/gamma_lat_safe_eq` (py:142).

**Assessment:**
Checkpoint higher-bar independence is GENUINELY MET. The `.wl` now contains at least five native
`Solve[...]` steps and a `D[Log[V[r]],r]` log-derivative step that are ALL absent from the `.py`,
a different decomposition/grouping, and a different substitution mechanism (Solve vs subs). Crucially,
the construction is NOT re-asserting the same hardcoded closed forms as residuals: each `expectZero`
reaches the LHS via the independent route (Solve/D[Log]) and compares it against the paper's
*independently stated* closed form (e.g. `chiLambdaLattice - 2 lambdaRef/rTurn`,
`kEffReq - EStar lambdaRef^2 VprimeTurnAbs/(lambdaPhys^2 rTurn)`), so a shared-algebra error would now
be catchable. The `Solve` targets are non-vacuous (e.g. `kEff` appears multiplied by `rTurnPhys`, so
the solve actually divides through; the threshold solve divides by `zetaEp`). Engines still agree at
every check; benchmark numbers unchanged (119.23361317…, 10.95423248, etc.).

### F3 — tautological_check

**Classification:** resolved

**What changed:**
py:114 `assert sp.simplify(r_turn_phys - lambda_phys*r_turn/lambda_ref) == 0` →
`assert sp.simplify((lambda_ref/lambda_phys)*r_turn_phys - r_turn) == 0`.
wl:176 `expectZero["turning-point radius map", rTurnPhys - lambdaPhys rTurn/lambdaRef]` →
wl:199 `expectZero["turning-point radius round-trip", (lambdaRef/lambdaPhys) rTurnPhys - rTurn]`.

**Assessment:**
The `expr − expr ≡ 0` tautology is GONE in BOTH engines. The new assertion is a genuine
forward∘inverse round-trip against the inverse length map used in §3.1 (`r = (λ_ref/λ_phys) r_phys`):
it composes the forward map (`rTurnPhys = λ_phys r_turn/λ_ref`) through the inverse and must return
`r_turn`. It depends non-vacuously on λ_ref, λ_phys, r_turn and reduces to 0 only by genuine
cancellation — it would FAIL if the forward and inverse maps were inconsistent. The F2/F3 interaction
note was honored: the round-trip form is preserved and the literal self-comparison was not
reintroduced. Both engines PASS this check (sympy assert holds; wl:46-47 "turning-point radius
round-trip = 0 / PASS").

### F4 — stale label (informational)

**Classification:** resolved

**What changed:**
wl:65 `banner["STAGE 236 — …"]` → `banner["STAGE 253 — …"]`.

**Assessment:**
Correct. Mathematica transcript line 8 now reads "STAGE 253 — PHYSICAL CALIBRATION AND
MATERIAL-THRESHOLD COMPANION". String-only; no mathematical impact.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L26 `r_turn^phys = lambda_phys*r_turn/lambda_ref` (map intact)
- L53 `micro product thresh = 119.23361317477 / zeta_ep` (matches notes after F1 fix)
- L58 `a_int stiffness coeff = 10.95423248`; L65 "All symbolic and numerical checks passed."
All 23 SymPy asserts execute without AssertionError (clean exit through the final banner).

**Mathematica:** exit=0. 26 PASS lines counted (banner "STAGE 253"):
- 3 §I + 2 §II + 6 §III + 2 §IV + 3 §V = 16 symbolic `expectZero` PASS, plus 9 §VI `expectNear`
  PASS = 25 named checks; final L116 "All Stage 253 symbolic and benchmark-specialization checks
  passed."
- L46-47 round-trip PASS (F3); L89 `micro product threshold = 119.23361317476522` (max diff
  1.42e-14 ≤ 1e-9); L94 `a_int stiffness coefficient = 10.95423248`. Engines agree everywhere.

**Output freshness:** confirmed regenerated post-fix. Script mtimes 13:31:41 (.py) / 13:32:22 (.wl);
output `.txt` mtimes both 13:34:52 — newer than scripts, so outputs reflect the edited code.

## Material-change assessment

`material_change`: false. F2 re-authoring lands on identical final values and claims; F3 swaps a
vacuous self-comparison for a non-vacuous round-trip that still residual-equals 0; F4 is a print
string; F1 corrected two notes typos to match the already-asserted engine values. No carried/derived
constant changed. Stage 253 is the end-of-ladder companion and the disputed benchmark is a
calibration readback, not a forward constant — no downstream unit consumes it. No `upstream_stale`
propagation warranted on physics grounds (253 is the last stage regardless).

## Side observations (non-blocking)

- §V screening inequalities `Π_• ≥ 1` are correctly NOT asserted as theorems in either engine; only
  the ratio identities (`PiEp - lambdaEpOmegaD/thresholdLambda`, etc.) are tested — appropriate, as
  these are host-screening criteria with no host data in scope.
- §VI benchmark constants are declared benchmark inputs recombined arithmetically and asserted as
  readback; they are not hand-pasted derived results bypassing the symbolic path.
- F1 directive's `## Applied: F1` block lists the notes file as files_changed — consistent with the
  observed diff. Minor: the original report's F2 file path had a `toy_projects` typo in the body; the
  directive's F2 target path is correct and the right file was edited.

## Verdict justification

All four findings are `resolved`. F1's two notes typos are fixed to the engine-corroborated values
with the published card and script/.wl literals left untouched (verified by git diff and grep). F2's
checkpoint higher-bar is genuinely met: the `.wl` now reaches the same final claims by an independent
route (native `D[Log[V[r]]]`, multiple non-vacuous `Solve[...]` balances, regrouped Π_ep), comparing
independently-derived LHS against the paper's stated closed forms so a shared-algebra error is
catchable — not a transliteration and not a residual re-statement of hardcoded forms. F3's tautology
is eliminated in both engines and replaced by a real forward∘inverse round-trip. F4 is corrected.
Both engines exit 0 with all 23 SymPy asserts and 25 Mathematica checks passing; outputs are fresh.
No regressions in the diff; `material_change` is false. This completes the first end-to-end pass.
