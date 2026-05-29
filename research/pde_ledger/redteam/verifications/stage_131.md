---
unit_id: 131
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 131 (red-team v2: R1/R2/R3 review findings)

This is the v2 verification covering the three findings raised by the read-only Codex review
(`redteam/codex_reviews/stage_131.md`): R1 (F1) tautological Anchor-3, R2 (F2) weak off-star
branch check, R3 (F3) Mathematica Π_* transliteration. Resolutions were agreed Claude+Codex in
`redteam/codex_reviews/_consult_batch4.md` (Q2 drop Anchor-3; Q3 corrected sign; Q4 branch
checks + g_+ value). The orchestrator independently re-ran both engines; I verified against the
committed `.txt` transcripts, the post-fix scripts, and `redteam/exec_logs/stage_131_diff.patch`.
(This file supersedes the 2026-05-27 v1 verification, which covered the original F1–F4 directive.)

## Per-finding outcomes

### F1 — tautological_check (R1)

**Classification:** resolved

**What changed:**
The tautological "Anchor 3 / Anchor (3)" parent-threshold-identity block was DELETED in both
engines with NO replacement (consult Q2, option (ii) DROP — even the proposed round-trip was
itself tautological).
- SymPy `scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py`: the old
  `:55–62` block (`threshold_at_star = threshold_residual.subs(...)`, `expected_form = ...`,
  `assert sp.simplify(threshold_at_star - expected_form) == 0`, and
  `print("PASS: parent threshold identity at Pi = Pi_* matches notes Sec. 2")`) is gone. The
  script now runs straight from Anchor 2's PASS (`:53`) into the F2 "Anchor (4)" block (`:55`).
- Mathematica `mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`:
  the old `:62–73` block (`thresholdAtStar = FullSimplify[...]`, `expectedForm`,
  `identityResidual = Chop[Simplify[thresholdAtStar - expectedForm], 10^-30]`, and the
  `If[TrueQ[identityResidual === 0], pass["...(notes Sec. 2)"], ...]`) is gone. The script runs
  from Anchor 2 (`:66`) straight into the F2 "Anchor 4" block (`:68`).
The diff confirms both deletions, including removal of the old PASS / `pass` lines mentioning
"parent threshold identity" / "(notes Sec. 2)".

**Assessment:**
Clean drop, no dressed-up tautology substituted. Neither committed transcript contains the
strings "parent threshold identity", "Sec. 2", or any new round-trip PASS line. The
`threshold_residual` definition (SymPy `:34`, Mathematica `:52`) is retained but feeds only an
informational `print` ("Parent bias mismatch formula = ..."), not any assertion — so no X−X
check survives. Matches consult Q2 exactly. Correct.

### F2 — insufficient_verification (R2)

**Classification:** resolved

**What changed:**
The weak single off-star check (`gPi(2*Pi_*) - g_-` asserted `> 1e-3`) was replaced by THREE
substantive branch-discrimination checks in both engines.
- SymPy `:55–88`: (4a) lower-branch membership `|g_Pi(Pi_*) - g_-^{F1}| < 1e-30`; (4b) singular
  branch excluded via `g_nat - g_Pi(Pi_*) ≈ Δg_- = 0.241964921055337` (tol 1e-12) AND `> 1e-3`;
  (4c) upper branch `|g_Pi(Pi_*) - g_+^{F1}| > 1` with
  `g_plus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) + 37*sp.sqrt(3))/(20*sp.pi)`.
- Mathematica `:68–95`: 4a `lowerResidual < 10^-30`; 4b `Abs[singSep - deltaGMinusNotes] < 10^-12
  && singSep > 10^-3` with `deltaGMinusNotes = 0.241964921055337`; 4c `upperSep > 1` with the
  same `gPlusExact` closed form.

**Assessment:**
All three checks present in both engines; old single off-star PASS line gone from both
transcripts (no "gPi(2*Pi_*)", no "lower-branch discrimination (paper Checks item 3)"). The Δg_-
literal is notes-grounded, not fabricated: `1 − g_-^{F1} = 1 − 0.758035078944663 =
0.241964921055337`, exactly the asserted constant (notes stage122:104). Substantive — each
assertion can fail: SymPy output shows `|g_Pi(Pi_*) - g_-^F1| = 5.51e-42` (membership; fails if
Π_* is a wrong root), `g_nat - g_Pi(Pi_*) = 0.241964921055337173...` (matches Δg_- to printed
precision), and `|g_Pi(Pi_*) - g_+^F1| = 2.0399...` (>1). The g_+^{F1} closed form ≈ 2.79795 is
the consult-Q4-confirmed boxed value. Together they exercise paper Checks item 3
(lower vs singular vs upper). Correct.

### F3 — transliteration (R3)

**Classification:** resolved

**What changed:**
Mathematica `:41` (old `piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, ...]`) was replaced
by an independent cleared-denominator route (`:45–47`):
`gThresholdResidual[p_] := 40*Pi*p*(2*p*Exp[p] + Pi) - 20*Pi*gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1);`
solved via `FindRoot[gThresholdResidual[piM] == 0, {piM, 1.4, 1.6}, ...]` (bracketing seed pair).
SymPy is unchanged (`sp.nsolve(gPi - g_minus, 1.5, ...)` at `:24`).

**Assessment:**
Sign is correct: the gMinus term is `- 20*Pi*gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1)`, i.e.
`(Exp[p] - 1)`, NOT the sign-flipped `-(1 - Exp[p])` draft that evaluated to ≈6366 at Π_*
(consult Q3 confirmed this corrected form). The route is genuinely independent of SymPy: a
different equation FORM (polynomial-in-(piM, Exp[piM]) vs the rational `gPi - g_minus`) AND a
different solver path (bracketing `{1.4, 1.6}` FindRoot vs single-seed Newton `nsolve`). `gMinus`
flows in symbolically from F4 rather than as a re-typed literal — not X−X. Mathematica output
confirms `Pi_* = 1.5088295134931555830...` still prints and Anchor 1 PASSES (residual ~4.4e-15 <
1e-14); it agrees with SymPy's `1.50882951349315558300...` to all printed digits. Anchors 2 and
4 (which reuse `piStar`/`gPi`) still PASS. Matches the IV.5 139/143/144 precedent for
transcendental numerical roots. Correct.

## Reconcile (unchanged surface)

- **Anchor 1** — Π_* vs `1.50882951349316`: SymPy output line 6 PASS; Mathematica output line 12
  PASS (residual ~4.4e-15).
- **Anchor 2** — g'(Π_*) vs `0.0714453558083195`: SymPy output line 7 PASS; Mathematica output
  line 14 PASS.
- **F4** — `g_-^{F1}` closed form `(2*sqrt(4107 - 100*pi**2) - 37*sqrt(3))/(20*pi)` vs
  `0.758035078944663`: SymPy output line 1 PASS; Mathematica output line 6 PASS. Untouched (not
  in diff).
- **Banner** — Mathematica `:26` reads "STAGE 131 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH
  COMPENSATION"; output line 3 confirms. Untouched (not in diff).

All four are present, correct, and outside the diff.

## Exec log / PASS-line assessment

**SymPy:** exit=0 (no exception → all asserts passed). 6 PASS lines:
1. `PASS: g_-^F1 closed form matches literal 0.758035078944663` (F4)
2. `PASS: Pi_* matches notes Sec. 1 value 1.50882951349316` (Anchor 1)
3. `PASS: g'(Pi_*) matches notes Sec. 3 value 0.0714453558083195` (Anchor 2)
4. `PASS: Pi_* on lower branch — |g_Pi(Pi_*) - g_-^F1| = 5.51...E-42` (F2 4a)
5. `PASS: singular equal-normalized branch excluded — g_nat - g_Pi(Pi_*) = 0.241964921055337...` (F2 4b)
6. `PASS: upper branch excluded — |g_Pi(Pi_*) - g_+^F1| = 2.0399...` (F2 4c)

**Mathematica:** exit=0 (`Exit[0]`; "Stage 131 Mathematica audit passed." printed). 6 PASS lines:
1. `PASS: g_-^F1 closed form vs literal` (F4)
2. `PASS: piStar notes Sec. 1 value` (Anchor 1)
3. `PASS: slope at piStar notes Sec. 3 value` (Anchor 2)
4. `PASS: Pi_* on lower branch (membership)` (F2 4a)
5. `PASS: singular equal-normalized branch excluded (notes Delta g_-)` (F2 4b)
6. `PASS: upper branch excluded` (F2 4c)

Both engines: 6 PASS lines each. The three new F2 lines appear in both; the removed lines (F1
"parent threshold identity ... Sec. 2", old F2 single off-star "lower-branch discrimination
(paper Checks item 3)") are absent from both. A passing log is necessary but not sufficient — I
confirmed each surviving assertion is falsifiable and tied to an independent notes/closed-form
anchor, not a restatement of itself.

**Output freshness:** Both committed `.txt` transcripts reflect the post-fix scripts — they print
the three F2 PASS lines and the F3-route Π_* (matching), contain no threshold-identity / off-star
line, and carry the corrected STAGE 131 banner. Consistent with the post-fix source and the diff.
(Per instructions I did not re-run; the orchestrator's independent re-run produced these.)

## Material-change assessment

`material_change`: false.

Verification-surface change only. F1 removed a definitional tautology (no value). F2 added
discrimination checks against imported/derived constants (Δg_- and g_+^{F1}, used only for branch
separation — not new ledger outputs). F3 changed HOW Mathematica derives Π_* but it still equals
1.50882951349316 (verified: Mathematica prints 1.5088295134931555830..., agreeing with SymPy and
the notes literal; Anchor 1 PASS). The load-bearing numerics (Π_*, g'(Π_*), g_-^{F1}) are
unchanged. No downstream unit depends on any newly-introduced quantity; no re-audit warranted.

## Side observations (non-blocking)

- `threshold_residual` (SymPy `:34`) / `thresholdResidual` (Mathematica `:52`) and the `V1`/`v1`
  expressions remain defined and printed but are now used only in informational `print`
  statements after the Anchor-3 deletion. The F1 directive explicitly permitted leaving them
  ("you MAY leave it — do not chase unused-variable cleanup"). Not a defect.
- The SymPy comment header still reads "Anchor (4)" even though Anchor (3) is gone; purely
  cosmetic numbering, no effect on checks.

## Verdict justification

All three review findings are resolved. F1: the X−X parent-threshold-identity block is cleanly
deleted in both engines with no replacement and no dressed-up tautology; both transcripts confirm
the old PASS line is gone. F2: the weak off-star check is replaced by three substantive,
falsifiable branch checks (lower membership, singular exclusion against the notes-grounded Δg_- =
0.241964921055337, upper exclusion via the g_+^{F1} closed form), all passing with the right
numbers. F3: the Mathematica Π_* derivation now uses the correctly-signed cleared-denominator
residual solved over a bracketing seed, structurally independent of SymPy's rational single-seed
nsolve, and Π_* still equals 1.50882951349316. Reconciled Anchors 1/2, F4, and the STAGE 131
banner are present, correct, and untouched. Both engines exit 0. No regressions in the diff;
material_change false. Verdict: verified.
