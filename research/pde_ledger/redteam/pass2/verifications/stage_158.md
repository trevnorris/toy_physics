---
unit_id: 158
batch: IV.6
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T12:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 158

## Per-finding outcomes

### F1 — mathematica_transliteration (orchestrator-elevated; user-authorized re-author)

**Classification:** resolved

**What changed:**
The `.wl` derivation route was re-authored to use analytic differentiation at the
base point instead of Series-expanding the same shifted closed form the SymPy
script uses. Concretely (current `.wl`):

- `delta R` law (wl:33-35): `rBase = rFun /. g->gStar`; `rSlope = D[rFun, g] /. g->gStar`;
  `rLin = Expand[rBase + rSlope*dg]`. Replaces HEAD's
  `rLin = Normal[Series[rFun /. g -> gStar + dg, {dg,0,1}]]`.
- `delta Mq` law (wl:42-48): `mQLin` built from `D[mQBase, sigmaVar]` and
  `D[mQBase, rVar]` partials evaluated at the base point. Replaces HEAD's
  `Expand[mQ] /. dSigma0*dR -> 0`.
- `delta Pi` law (wl:54-61): `piLin` built from the three first partials of `piBase`.
  Replaces HEAD's `Expand[piExpr] /. {dSigma0*dR->0, dSigma0*dS->0, dR*dS->0}`.
- composed-law slope `dRFromDg` (wl:67-70): now derived as
  `D[(gSym-rSym)^2/(1+rSym^2), gSym] /. gSym -> rSym - Sqrt[1+rSym^2]/2 * dgSym`
  instead of the hardcoded `-dgSym/Sqrt[1+rSym^2]`.
- `Delta_Q (chi)` law (wl:88-90): `chiBase = chi/.eps->0`; `chiSlope = D[chi,eps]/.eps->0`;
  `chiLin = Expand[chiBase + chiSlope*eps]`. Replaces HEAD's
  `chiLin = Expand[Normal[Series[chi, {eps,0,1}]]]`.
- Numerical-coefficient block (wl:101-131): every coefficient is now produced by
  differentiating the underlying closed form (`dRdGCan = D[rNumFun,gNum]/.gNum->gNumStar`,
  `D[mQNum,...]`, `D[piNum,...]`, `D[sigmaThatNum,tNum]`) and substituting the canonical
  values, instead of re-typing the SymPy closed forms (`-1/sqrt1`, `sigma0Can/sqrt1`, ...).

**Assessment:**
Correct and complete against the directive's acceptance criteria:

1. **Genuinely independent of the `.py`.** The two load-bearing checks (`delta R`,
   `linear Delta_Q`) no longer Series-expand the shifted closed form — they take
   analytic derivatives at the base point. The `.py` still uses `sp.series(R_shift, dg, 0, 2)`
   (py:40-41) and `sp.series(chi, eps, 0, 2)` (py:90), so the algorithmic pathways now
   genuinely differ. The gain/slope checks likewise switched from expand-and-drop-cross-term
   (the `.py` mechanism, py:50/59) to analytic partials. The numeric block now differentiates
   rather than re-types the `.py` closed forms. The second engine re-derives rather than echoes.

2. **Assertion set UNCHANGED.** All six `expectZero[...]` calls (names and expected-zero
   targets) are byte-identical to HEAD: `linear delta R law` (rLin-rExpected),
   `delta Mq law`, `delta Pi law`, `composed delta Mq law`, `composed delta Pi law`,
   `linear Delta_Q law`. The diff touches only derivation lines, never the assertion calls.

3. **Numeric values UNCHANGED.** The fresh Mathematica exec log reproduces the committed
   transcript to full precision: dR/dg = -0.49021604438762603754..., dMq/dSigma0 = -1/4,
   dMq/dg = 2.28001126927792351405..., dPi/dSigma0 = 0.83240947108163457213...,
   dPi/dS = -1.16275838754221894078..., dPi/dg = 1.52843317823248362127...,
   dSigma0/dThat = 6.42981496203005499347..., dPi/dThat = 5.35223887169621835652... —
   each matching the directive's required value exactly.

4. **No new tautology.** Each re-authored law is compared to a SEPARATELY hand-typed closed
   form: `rLin` (from D[rFun,g]) vs `rExpected = 1/4 - dg/Sqrt[1+r^2]`; `chiLin`
   (from D[chi,eps]) vs `chiExpected = 1 + eps*(5b+a0/3+9a5)`; the partial-built `mQLin`/`piLin`
   vs the hand-typed `-rStar*dSigma0 - sigma0*dR` / `dPiExpected`. These are distinct
   expressions, so a wrong derivative would make the residual nonzero — the checks are can-fail,
   not X-X. The composed checks still glue `dRFromDg` (now itself a derivative) into the boxed
   forms, also can-fail.

No collateral edits beyond the named derivation route: the diff is confined to the `.wl`;
the `.py`, paper, and notes are untouched (directive scope respected). The only structural
additions are scratch symbols (`sigmaVar`, `rVar`, `piSigmaVar`, `gNum`, etc.) cleared with
`Clear[...]`, used only to host the symbolic differentiation.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`linear delta R law = 0` / `delta Mq law = 0` / ... / `linear Delta_Q law = 0`;
`dR/dg = -0.49021604438762605982`; `# exit_code: 0`. (SymPy script unchanged this iteration;
log confirms it still passes and supplies the engine-agreement reference values.)

**Mathematica:** exit=0. Notable lines:
`PASS: linear delta R law` ... `PASS: linear Delta_Q law` (all six PASS);
`dR/dg = -0.49021604438762603754...`; `Stage 158 Mathematica audit passed.`; `# exit_code: 0`.
Engine agreement with SymPy preserved (coefficients match to ~17 digits, as in HEAD).

**Output freshness:** confirmed. Committed `.txt` mtimes (mathematica + sympy outputs:
2026-06-08 11:56:32) are newer than the `.wl` script mtime (11:49:45), so the saved transcript
was regenerated post-fix. `git diff HEAD` on the committed `.txt` is empty — the transcript is
byte-identical to HEAD, exactly as required (route changed, every value preserved).

## Material-change assessment

`material_change`: false.

This was a derivation-route change only. The asserted identities, their expected-zero targets,
and every printed numeric coefficient are unchanged in value (committed transcript byte-identical
to HEAD). No derived result that a downstream unit could depend on changed, so no downstream unit
is invalidated by this edit.

## Side observations (non-blocking)

- The new numeric block re-derives `dRdGCan` via `D[rNumFun,gNum]/.gNum->gNumStar` and reuses
  it for both `coefDRDG` and `coefDMQDg`/`coefDPiDG`; this is the genuinely-independent numeric
  route the directive offered as option (1)-style differentiation, and it additionally hardens
  the numerics (they are now computed, not re-typed). No issue.
- `coefDSigmaDT` is now `D[(20/9)*tNum^2, tNum] /. tNum->tCan = (40/9)*tCan`, matching the
  `.py` `(40/9)*T_can` and the notes' traction relation; the `40/9` is recovered by
  differentiation rather than hardcoded. No issue.

## Verdict justification

All four acceptance criteria hold: the script exits 0 with all six `expectZero` printing
`PASS:` and the final pass line; the asserted-identity set and expected-zero targets are
byte-identical to HEAD; every printed numeric coefficient is unchanged in value (committed
transcript byte-identical to HEAD, full precision matching the directive); and the two
load-bearing laws (`delta R`, `Delta_Q`) plus the gain/slope/numeric blocks no longer mirror
the SymPy Series/expand-and-drop mechanism — they use analytic differentiation at the base
point, a genuinely independent route, with no tautology introduced. F1 is resolved.
