---
unit_id: 196
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 196

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine gap)

**Classification:** resolved

**What changed:**
Codex created a new (untracked) independent Mathematica audit
`mathematica/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_mathematica_audit.wl`
plus its committed output `mathematica/output/...mathematica_audit.txt`. The SymPy
script `scripts/...sympy_audit.py` was NOT modified (mtime `May 11 11:58`, unchanged in
`git status`); only its `.txt` was refreshed by the orchestrator's reliability re-run.
The empty `stage_196_diff.patch` is expected — the `.wl` is a brand-new untracked file.

**Assessment — genuinely independent route (NOT a transliteration):**

The `.wl` is a real second-engine route, not a re-type of the `.py`. The decisive
evidence is Section III. SymPy *hardcodes* the Stage-194 deformation slots as literal
symbolic expressions (`L0_stage194 = -3*S + Sigma0`, `L2_stage194 = S*beta^2/3 + Sigma2`,
`L4_stage194 = S*beta^4/9 + Sigma4`, `L5_stage194 = S*beta^5/9 + Sigma5`, lines 123–126).
The Mathematica script instead *derives* those slots from first principles via the native
outgoing l=2 operator:
`lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2,x],x]/SphericalHankelH1[2,x]]`,
`lambdaWindow5 = Series[..., {x,0,5}]` → `-3 + x^2/3 + x^4/9 + (I/9)x^5` (output line 37),
then `deformedDtn = scale*(lambdaWindow5 /. x->stretch*z) + sig0 + sig2 z^2 + ...`.
I independently re-derived the spherical-Hankel l=2 log-derivative expansion and confirmed
it is exactly `-3 + x^2/3 + x^4/9 + (I/9)x^5`; substituting `x->stretch*z` and adding the
sig-slots reproduces SymPy's hardcoded `L0..L5` (`scale*(-3)+sig0`, `scale*stretch^2/3+sig2`,
etc.). So the matching slots are obtained by a genuinely different decomposition (native
DtN special-function expansion vs. SymPy's literal algebra). Likewise the canonical-even
matching is *solved* via `Solve[{coeff z^2 ==1/9, coeff z^4 ==4/81},{sig2,sig4},Reals]`
with a uniqueness guard (`If[Length[evenRules]=!=1, fail]`) rather than substituting the
notes' closed-form laws — a different route from SymPy, which hardcodes `Sigma2_match`/
`Sigma4_match` and only checks the consequence. Native primitives throughout: `Series`+
`Coefficient`, `D`, `Solve`/`Reduce`, `FullSimplify`, house idioms (`stripCE`, `$Assumptions`,
`math -script`). Not a line-by-line port; variable names and step order differ.

**Assessment — tail genuinely wired in and shown to cancel (the central check):**

Yes. The extra odd tail is a free symbol genuinely substituted into the denominators, and
the cancellation is demonstrated, not assumed:
- Section I: `retResponse[tail_] := 3/4 + 1/(4*(1 - frontCore - tail))` with
  `retWithTail = retResponse[I*tau7*w^7]` — `tau7` is a free real symbol wired into the
  pole. `retDifferenceLowCoeffs = Table[Coefficient[retDifferenceThrough7,w,k],{k,0,6}]`
  proven all-zero; `retThrough5 - retSympyLowTarget = 0` shows the *tailful* series through
  w^5 equals the tail-free one-pole form; and `retDifferenceThrough7 - I*tau7*w^7/4 = 0`
  confirms the tail re-enters at exactly w^7 (so it was genuinely present, not nulled).
- Section II: `dtnResponse[I*l7*z^7]` wires `l7` into the DtN denominator. Output line 24
  explicitly shows `- (I*l7*z^7)/l0` surviving in the z^7 series. `dtnDifferenceLowCoeffs`
  all-zero through z^6; `dtnThrough5 - dtnSympyLowTarget = 0`; and the l7-coefficient at z^7
  equals `-I/l0`. Tail present, cancels through z^5, re-enters at z^7.
- Section III (strongest): `deformedDtn` carries `+ I*l7*z^7`; `matchedDtn` (after the
  sig2/sig4 solve, which does not touch l7) still contains l7 at z^7. The check
  `expectZeroList["matched coefficients z^0..z^5 independent of L7",
  Table[D[Coefficient[matchedDtn,z,k],l7],{k,0,5}]]` = {0,0,0,0,0,0} differentiates each
  through-z^5 coefficient w.r.t. a tail symbol that is genuinely present — it would FAIL if
  any z^5-or-lower coefficient secretly depended on the tail. This is materially stronger
  than SymPy's A11/A12 (which the auditor flagged as trivially zero because their targets
  were already L7-free).
- χ_Q independently extracted: `z5Slot = Coefficient[matchedDtn,z,5]`,
  `chiFromSeries = z5Slot/(I/27)` (= `-27*I*z5Slot`, identical normalization to SymPy's
  `-27*I*z5_coeff`), yielding `3(scale*stretch^5 + 9 sig5)/(3 scale - sig0)`
  (output line 41) = χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀). Checked equal to the SymPy deformation
  formula and L7-independent (`D[chiFromSeries,l7]=0`).

No load-bearing check is vacuous. The only trivially-zero checks are the Section IV
`D[nNatural,l7]` / `D[deltaNatural,l7]` derivatives, which (exactly as the SymPy auditor
noted for A11/A12) are built from the already-l7-free `chiFromSeries`; they are redundant
confirmations, not the load-bearing evidence, so they do not weaken the verdict.

## Exec log assessment

**SymPy:** exit=0. The committed `.txt` ends with the Stage ledger and the log tail shows
`# exit_code: 0`. All `expect_zero` lines print `= 0` (e.g. `chi_Q extractor - deformation
algebra formula = 0`, `L7 enters normalized response first at z^7 = 0`). SymPy unchanged;
re-run only refreshed the output.

**Mathematica:** exit=0. Log tail shows `# exit_code: 0`. Every check prints a `PASS:` line
and the printed residuals are all `0` / all-zero lists, e.g.:
- `response-side tail absent through w^6 = {0, 0, 0, 0, 0, 0, 0}` → PASS
- `DtN tail absent through z^6 = {0, 0, 0, 0, 0, 0, 0}` → PASS
- `matched coefficients z^0..z^5 independent of L7 = {0, 0, 0, 0, 0, 0}` → PASS
- `chi_Q extractor - SymPy deformation formula = 0` → PASS
The `expectZero`/`expectZeroList`/`fail` helpers call `Exit[1]` on any nonzero residual, so
exit 0 confirms every assertion held; the helpers are not no-ops.

**Output freshness:** confirmed. `.wl` mtime `2026-06-01 12:00:35`; its `.txt` mtime
`2026-06-01 12:02:48` (newer than the script). SymPy `.txt` (`12:02:48`) also newer than the
`.py` (`May 11 11:58`). Both outputs regenerated post-fix.

## Cross-engine agreement

Yes — full agreement on every conclusion. Both engines produce identical closed forms:
σ_can = 9/(8 pole⁵); response series `1 + w²/(4 pole²) + w⁴/(4 pole⁴) + (9i χ/32) w⁵/pole⁵`;
DtN normalized response `1 - L2 z²/L0 + (L2²/L0² - L4/L0) z⁴ - i L5 z⁵/L0`; matching laws
Σ₂ = -(3Sβ²-3S+Σ₀)/9, Σ₄ = -(3Sβ⁴-3S+Σ₀)/27 (1/9 and 4/81 even slots); the L7-coefficient
at z⁷ = -i/L0; χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀); N_Q = (3S−Σ₀)/(3(Sβ⁵+9Σ₅)); and
Δ_norm = P0_target(1/χ_Q − 1) with P0_target = 54 G c_s⁵/(5 a⁵ c⁵). The `.wl` additionally
subtracts each of its derived forms against the explicit SymPy target and gets 0, so
cross-engine agreement is asserted in-script, not merely visually congruent.

## Material-change assessment

`material_change`: false. No derived result changed — the SymPy reference engine was not
touched and the new `.wl` independently reproduces every existing conclusion. Downstream
units are unaffected; no `upstream_stale` propagation is warranted on math grounds.

## Side observations (non-blocking)

- Stale internal labeling persists in the SymPy script (banner "STAGE 179", comments
  "Stage 194/195") vs. paper card "Stage 196". The new `.wl` correctly labels itself
  "STAGE 196". This is the same cosmetic mismatch the auditor recorded as non-blocking; the
  SymPy script was (correctly) not modified by this finding. Not a verification blocker.
- Section IV `D[nNatural,l7]` / `D[deltaNatural,l7]` are trivially zero (targets already
  L7-free), mirroring SymPy A11/A12 — redundant, not load-bearing. Noted only.

## Verdict justification

The single finding (missing second engine) is fully resolved: Codex added a genuinely
independent Mathematica route — deriving the Stage-194 slots from the native spherical-Hankel
l=2 DtN operator and solving the canonical-even matching, rather than re-typing SymPy's
hardcoded closed forms. The extra O(w⁷)/O(z⁷) odd tail is wired in as a free symbol in every
section, shown explicitly to survive into the z⁷ term, and proven to leave every coefficient
through O(z⁵)/O(w⁵) unchanged (the z⁰..z⁵ per-coefficient ∂/∂L7 = 0 check on the genuinely
tailful matched response is stronger than the SymPy analogue). χ_Q is independently extracted
and matches. Both engines exit 0 with all assertions zero, outputs are fresh, and conclusions
agree across engines. Verdict: verified.
