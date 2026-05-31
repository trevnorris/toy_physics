---
unit_id: 176
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T02:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 176

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy `scripts/...stage176...py:100-112`: the old vacuous check
  `Sigma_rigid = Sigma_factored.subs({dlnI:0, dlnH:0}); expect_zero(..., Sigma_rigid - 2*dlnM)`
  was replaced by a reduction of the *independent* exact drift:
  `rigid = {dGU: dOU + dGW - dR, dOW: 2*dR - dOU}`,
  `Sigma_exact_rigid = sp.simplify(Sigma_exact.subs(rigid))`,
  `dlnM_rigid = sp.simplify((2*dlnM).subs(rigid))`,
  `expect_zero("rigidity reduction of Sigma_exact to 2 d ln M", Sigma_exact_rigid - dlnM_rigid)`.
- Mathematica `mathematica/...stage176...wl:94-100`: analogous replacement using
  `rigid = {dGU -> dOU + dGW - dR, dOW -> 2*dR - dOU}`,
  `sigmaExactRigid = FullSimplify[sigmaExact /. rigid, ...]`,
  `dlnMRigid = FullSimplify[(2*dlnMExpr) /. rigid, ...]`,
  `expectZero["rigidity reduction of Sigma_exact to 2 d ln M", sigmaExactRigid - dlnMRigid]`.

**Assessment:**
Correct and non-tautological. The check now references `Sigma_exact`/`sigmaExact` (the
series/derivative of the actual `Lambda` definition), not the hand-built `Sigma_factored`/
`sigmaFactoredForm`. The rigidity substitutions are the correct primitive-variable solutions
of the two constraints stated in the script's own dln formulas: `dlnI = dR + dGU - dOU - dGW = 0`
gives `dGU = dOU + dGW - dR`, and `dlnH = 2*dR - dOU - dOW = 0` gives `dOW = 2*dR - dOU` — both
match the directive and the report's self-test. The comparison target `dlnM_rigid =
2*(dGW - (2*dR - dOU) - dK/2) = 2*dGW - 4*dR + 2*dOU - dK` is a nonzero literal in the primitive
drifts, so the residual is non-vacuous (it is NOT `a + b*0 + c*0 == a`). Because the new residual
flows from `Sigma_exact`, the auditor's falsifiability requirement is met: deliberately perturbing
a weight inside `Sigma_factored` would no longer leave this check passing — the check no longer
touches `Sigma_factored` at all. The Mathematica fix correctly uses `dlnMExpr` (the assigned
expression), not the unassigned bare symbol `dlnM`, per the directive caveat. Both transcripts show
`rigidity reduction of Sigma_exact to 2 d ln M = 0` (sympy log L30, math log L33-34 PASS).

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/...stage176...wl:60-61`: a two-line comment was added immediately above the
`sigmaExact = ...` block noting "Mathematica uses symbolic D[Log,eps] at eps->0; the SymPy audit
instead Taylor-expands and reads the eps^1 term." The `D[Log[...], eps] /. eps -> 0` extraction
(lines 62-65) is retained, distinct from the SymPy `sp.series(...).removeO()/eps` route.

**Assessment:**
Matches the directive exactly — the only mandated remediation for this low-severity flag was to
make the existing first-order-extraction divergence explicit via a comment, with no algebra rewrite
and no new physics. The distinct extraction mechanism is retained, the clarifying comment is
present, and the script exits 0. No collateral edits.

### F3 — stale_output

**Classification:** resolved

**What changed:**
- SymPy line 32: `banner("STAGE 159 — ...")` -> `banner("STAGE 176 — ...")`.
- Mathematica line 26: `banner["STAGE 159 — ..."]` -> `banner["STAGE 176 — ..."]`.
- Both saved transcripts regenerated; both now show `STAGE 176 — OUTGOING LOAD-FACTOR
  FACTORIZATION` in the top banner.

**Assessment:**
Correct. Both scripts and both committed `.txt` transcripts now read "STAGE 176". (Note: the banner
line is L8 of the exec logs and L3 of the saved `.txt` due to the leading blank line, not L11 as the
original report estimated — immaterial; the substantive correction is present.) Cosmetic-only change,
no numeric result affected. Both exit 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L8 `STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION` (F3 corrected)
- L10 `exact factorization of Lambda^2/K = 0`
- L18-19 `factored first-order defect formula = 0` / `expanded primitive-variable transport = 0`
- L30 `rigidity reduction of Sigma_exact to 2 d ln M = 0` (F1 new, non-tautological)
- L37 `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- L8 `STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION` (F3 corrected)
- L11/L15/L17 `PASS:` for factorization, factored defect, expanded transport
- L33-34 `rigidity reduction of Sigma_exact to 2 d ln M = 0` / `PASS:` (F1 new)
- L43 `# exit_code: 0`

**Output freshness:** confirmed. Script mtimes are both 2026-05-30 01:10:19; saved transcript
mtimes are 01:38:31 (sympy .txt) and 01:38:40 (math .txt) — both newer than their scripts. The
committed `.txt` files match the exec logs line-for-line (corrected banner + new check name + PASS).

## Material-change assessment

`material_change`: false.

No derived result changed. F1 strengthened a verification (replaced a vacuous check with a
falsifiable reduction of the same already-established quantity) but produced the identical residual
`= 0`; the underlying factorization, first-order decomposition, and rigidity reduction are unchanged.
F2 added only a comment. F3 corrected a banner string. No downstream unit can depend on any of these
edits, so no units > 176 need be marked upstream_stale on account of this stage.

## Side observations (non-blocking)

- The Mathematica `dlnM` symbol remains unassigned and is used inside `sigmaFactoredForm` (line 72)
  before being substituted by `dlnMExpr` (line 73) for the A2 check; the F1 fix correctly side-steps
  this by using `dlnMExpr` directly. This pre-existing pattern is intentional and not a defect.
- Both engines exercise the single one-port case (subscript `r` dropped), consistent with the
  original report's scope; not a regression.

## Verdict justification

All three findings are resolved. F1 (the one substantive finding) was fixed exactly as the directive
prescribed: the rigidity-corollary check now derives its residual from the independently computed
`Sigma_exact` reduced under the correct primitive rigidity constraints, compared against the same-
surface reduction of `2*dlnM`, yielding a non-vacuous residual that genuinely exercises deliverable
D4 and would fail under a broken drift relationship. F2 and F3 were the minimal documentation/cosmetic
remediations mandated. Both engines exit 0, both transcripts are fresh and match the exec logs, and the
diff shows no collateral edits or regressions. Verdict: verified.
