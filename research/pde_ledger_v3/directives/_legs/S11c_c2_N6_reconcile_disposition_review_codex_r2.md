# Codex-sol (gpt-5.6-sol xhigh) review — S11c-c2 N6 reconcile DISPOSITION (round 2, post-fold)

⚠ Raw codex transcript (~668 KB, echoed source-file contents) trimmed for tree hygiene (regenerable). Review
content verbatim below. Verdict: **NOT-SOUND** on ONE header-scoping line (now folded); everything else confirmed
clean. [Round-2 folded that line → cleared.]

---

The disposition is accurate except for one overbroad sentence that must be narrowed.

## Must-fix finding — universal matched-agreement claim
**Artifact:** disposition line 35: "Cross-engine, every MATCHED comparison agreed"
**Sources:** the committed tally records matched symbolic disagreements — `N6RC_CARRIER_EULERIAN`/`CARRIER_MATERIAL`
40 NONZERO each (tally:19); `SOURCE_ACTUAL`/`BASELINE`/`PREDICTED` 76 NONZERO each (tally:12); `FROZEN_PHI` 18
NONZERO (tally:6). The question doc calls these "SURFACED computed residuals — matched keys, some NONZERO"
(question:41).
**Required change:** replace with e.g. "Cross-engine, every matched comparison in the following
vanishing/control/premise subset agreed…". Without that restriction, line 35 wrongly claims agreement for the
surfaced-nonzero operands, even though §2 later reports them correctly.
[FOLDED round-2: header restricted to "…IN THE VANISHING / CONTROL / PREMISE SUBSET…" + an explicit "⛔ this is
NOT 'every matched key agreed' — carrier/source/Φ have matched keys, 40/76/18 of THOSE are NONZERO (§2)".]

## Other checks (all confirmed clean)
- All numerical counts match the committed tally: 160/0 covariance families + control delta, 320/0 carrier bridge,
  400/0 support, 40/76/18 operand gaps, 4+4 census gaps.
- The `(0)−(0)` treatment is correct. Both engines form `R_cov` + the carrier bridge as internal residuals
  (covariance_sympy.py:194, reconcile_sympy.py:250, .wl:874–886); their per-engine vanishings make the cross-engine
  zeros residuals-of-zeros. The disposition correctly separates the four nonzero `FROZEN_RELATIONS` premise
  agreements + the 400 support agreements from these trivial zeros.
- The EL calculation reproduces (independent re-derivation):
  `L = a θ'e' + b eW'θ'`, `T(a)=kR, T(b)=−kR/W, R=W_0/W`; `EL_θ(L) = −a e'' − b(e'W' + eW'')`; since `R' = −RW'/W`,
  `T(L) = kθ'(Re'−(R/W)eW') = kθ'(Re)'`; therefore
  `EL_θ(TL) − T(EL_θ L) = −k(Re)'' + kR e'' − (kR/W)(e'W' + eW'') = (kW_0/W²)W'e' − (2kW_0/W³)e(W')²`.
  First term = one background derivative, survives at σ_W^1; second is quadratic. Matches the counterexample
  (PATH_B:14–22).
- The grading obstruction is supported: WL applies profile rules before `gradePart` (.wl:127–141), SymPy replaces
  profiles before coefficient extraction (diagnostic:245–248); the first-order identity needs `(RW)_1 = R_0W_1 +
  R_1W_0`, excluded by the frozen no-cross-grade contract. The disposition properly qualifies that an upstream
  replay validates the replay, not the old streams, absent a separately-established relation.
- The density table is not upgraded to source agreement; "just thickness" is expressly forbidden; the debt is not
  dismissed; all required caveats + earlier carries appear (disposition:146–155). Reading B and c1 stand.

**NOT-SOUND** — line 35 must restrict "every MATCHED comparison agreed" to the enumerated vanishing/control/premise
subset. [Folded → cleared.]
