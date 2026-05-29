---
unit_id: 108
batch: IV.2
verifier_model: claude-opus-4-8
verify_date: 2026-05-29T00:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: true
---

# Red-team verification — unit 108

Scope: confirm Codex's edits closed the in-scope findings F2, F3, F4. F1 is a
`paper_misalignment` that is DELIBERATELY GATED (held for user resolution as a
paper-card cleanup; checks #2/#3 covered at sibling stages 110/111/112). F1 is
classified `blocked_legitimate` and EXCLUDED from the pass/fail rollup. Overall
verdict depends only on F2/F3/F4.

Note: this file supersedes a stale v1 verification (dated 2026-05-27) whose
F1/F2/F3 referred to a different finding set than the current v2 directive.

## F1 — paper_misalignment (GATED / EXCLUDED)

**Classification:** `blocked_legitimate`

Per the directive (`## F1` + `## Resolve before fix_loop`), Codex was instructed
to do nothing on F1; it is a paper-side scoping question (do the Robin /
standalone mixed-pole / compensated-Robin–mixed checks belong to stage 108 or to
a sibling stage ≈109–113?). It is being resolved separately as a paper-card
cleanup, NOT a script change. The diff confirms Codex touched ONLY the F2/F3/F4
lines and made no attempt to add Robin/mixed/compensated material to either
script. This is the correct, intended behavior. Excluded from the rollup as
instructed.

## F2 — mathematica_transliteration

**Classification:** `resolved`

**What changed:** An additive independent-route cross-check was inserted in the
`.wl` Class D block (`mathematica/...stage108..._mathematica_audit.wl:102-110`),
AFTER the existing series-route `chiGen` (line 100, untouched) and printed
`chi_gen(beta)` (line 101), and BEFORE the `sigma5PresGen` solve (line 111):

```wolfram
L0raw = -3*sNorm + sigma0;
L2raw = sNorm*beta^2/3 + sigma2;
L4raw = sNorm*beta^4/9 + sigma4;
L5raw = sNorm*beta^5/9 + sigma5;
solAlt = First[Solve[{-L2raw/L0raw == 1/9, L2raw^2/L0raw^2 - L4raw/L0raw == 4/81},
                     {sigma2, sigma4}, Reals]];
chiGenAlt = FullSimplify[(27*(-L5raw/L0raw)) /. solAlt, Assumptions -> $Assumptions];
expectZero["chiGen independent-route agreement", chiGenAlt - chiGen];
```

**Assessment (directive checks (a)-(d)):**
- (a) The new route EXISTS and does NOT use `Series` of the rational function.
  The raw L-coefficients are written directly as symbolic polynomials in
  `{sNorm,beta,sigma0,sigma2,sigma4,sigma5}`; the even-fingerprint constraints are
  plain algebraic equations on those symbols (`-L2raw/L0raw == 1/9`,
  `L2raw^2/L0raw^2 - L4raw/L0raw == 4/81`); `chiGenAlt = 27*(-L5raw/L0raw)`.
  No `Series`, no `Coefficient`, no `D[...]`/derivative. Confirmed.
- (b) It shares only the physical premise — the literal `lambdaOut` coefficients
  embedded in `-3*sNorm`, `sNorm*beta^2/3`, `sNorm*beta^4/9`, `sNorm*beta^5/9` —
  not the algorithm. The series route reaches `chiGen` via
  `Series → Coefficient → Solve → assemble` (lines 85-100); the raw route hardcodes
  the coefficients and solves plain equations. Different choreography. Confirmed.
- (c) The existing series-route `chiGen` (line 100) and its print (line 101) are
  untouched. The diff shows the F2 block is purely additive (inserted after line
  101); no series-route line was deleted or altered. Confirmed.
- (d) The exec log (`stage_108_mathematica.log:40-41`) shows
  `chiGen independent-route agreement = 0` followed by
  `PASS: chiGen independent-route agreement`. Confirmed.

**Genuine second-engine cross-check, not a tautology:** The two routes encode the
same physics (same `lambdaOut` premise) by structurally different algebra. The
residual `chiGenAlt - chiGen` is 0 ONLY because both correctly reproduce
`3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0)`. A wrong coefficient in either
route (e.g. a mistyped truncation order in the series route, or a wrong `L4raw`
in the raw route) would make the residual nonzero. This is a real
algorithm-independent confirmation, satisfying the second-engine policy.

## F3 — tautological_check (Class A)

**Classification:** `resolved`

**What changed:** The bare S-cancellation comparison was removed in both engines
and replaced with an anchor to the literal canonical fingerprint
(appendix eq:app-part04-Yout-dtn):
- SymPy (`scripts/...stage108..._sympy_audit.py:31-33`): `Y_can` (the series of
  `-3/Lambda_out`) is GONE; replaced by
  `Y_can_literal = 1 + z**2/9 + 4*z**4/81 + I*z**5/27`, asserted via
  `expect_zero('pure scale invariance', sp.expand(Y_scale) - Y_can_literal)`.
- Mathematica (`...stage108..._mathematica_audit.wl:35-38`): `yCan` (the series of
  `-3/lambdaOut`) is GONE; replaced by
  `yCanLiteral = 1 + z^2/9 + (4*z^4)/81 + (I/27)*z^5`, asserted via
  `expectZero["pure scale invariance", yScale - yCanLiteral]`.

The diff (`stage_108_diff.patch:8-13`, `59-66`) confirms both `yCan`/`Y_can`
self-construction lines were deleted, not merely supplemented.

**Assessment:** Falsifiable. `Y_scale` is still computed from `Lambda_out`
(`series((-3*S)/(S*Lambda_out))`), so a wrong `Lambda_out` coefficient now shifts
`Y_scale` while the literal target is fixed, breaking the check. The scale `S`
still drops out (it cancels inside `Y_scale`), so Class A's claim is still
exercised — but now against an external paper anchor rather than against itself.
Exec logs: `stage_108_sympy.log:10` `pure scale invariance = 0`;
`stage_108_mathematica.log:10-11` `pure scale invariance = 0` / `PASS`.

## F4 — tautological_check (Class C and Class D loci)

**Classification:** `resolved`

**What changed (directive checks (a)-(c)):**
- (a) Class C locus is now anchored to the VALUE, not a round-trip:
  - SymPy (`...sympy_audit.py:69-73`): the old
    `expect_zero('preservation locus check', chi_add.subs(Sigma5, chi_pres) - 1)`
    round-trip is GONE; replaced by
    `expect_zero('Sigma5 locus (Class C) = -Sigma0/27', sp.simplify(chi_pres) - (-Sigma0/27))`.
  - Mathematica (`...mathematica_audit.wl:79-81`): the old round-trip
    `expectZero["preservation locus check", (chiAdd /. sigma5 -> sigma5Pres) - 1]`
    is GONE; replaced by
    `expectZero["Sigma5 locus (Class C) = -sigma0/27", sigma5Pres - (-sigma0/27)]`.
    (Codex took the directive's accepted alternative: the prior anchored line
    `Sigma5 preservation locus + sigma0/27` was equivalent; it consolidated to the
    `-sigma0/27` form and deleted the round-trip. Either path was explicitly
    permitted by the directive's F4 note.)
- (b) The general/Class D round-trip is DEMOTED to a print in both engines:
  - SymPy (`...sympy_audit.py:105-106`): now
    `print('general preservation locus check =', sp.simplify(chi_gen.subs(Sigma5, chi_pres_gen) - 1))`
    — no longer an `expect_zero`.
  - Mathematica (`...mathematica_audit.wl:117`): now
    `Print["general preservation locus check = ", fmt[FullSimplify[(chiGen /. sigma5 -> sigma5PresGen) - 1, ...]]]`
    — no longer an `expectZero`.
- (c) The submanifold anchor (A7) still asserts and passes in both engines:
  - SymPy (`...sympy_audit.py:101-104`):
    `expect_zero('general preservation submanifold = S(1 - beta^5)/9 - Sigma0/27', ...)`.
  - Mathematica (`...mathematica_audit.wl:113-116`):
    `expectZero["general preservation submanifold = S(1 - beta^5)/9 - sigma0/27", ...]`.

**Assessment:** Falsifiable. The Class C anchor `chi_pres - (-Sigma0/27)` compares
the SOLVED locus value against the independently-stated notes value
`-Sigma0/27` (the beta=1 reduction of the submanifold). If `chi_add`'s closed form
were wrong, `solve(chi_add==1)` would yield a different `chi_pres` and the residual
would be nonzero — unlike the old round-trip, which was guaranteed by `solve`.
Demoting the Class D round-trip removes no real coverage: the submanifold anchor
(A7) is the genuine, load-bearing Class D test and remains an assertion. Exec
logs confirm: SymPy `stage_108_sympy.log:22` `Sigma5 locus (Class C) = -Sigma0/27 = 0`,
line 27 submanifold `= 0`, line 28 demoted print `= 0`; Mathematica
`stage_108_mathematica.log:35-36` Class C locus `= 0` / `PASS`, lines 43-44
submanifold `= 0` / `PASS`, line 45 demoted print `general preservation locus check = 0`
(a bare print, no PASS line, confirming it is no longer an assertion).

## Exec-log assessment

Both logs are complete (not truncated) and exit 0.

SymPy (`stage_108_sympy.log`, exit_code: 0):
- `pure scale invariance = 0` (F3)
- `Sigma5 locus (Class C) = -Sigma0/27 = 0` (F4 Class C)
- `general preservation submanifold = S(1 - beta^5)/9 - Sigma0/27 = 0` (F4 anchor intact)
- `general preservation locus check = 0` (F4 demoted print — informational, no assertion banner)

Mathematica (`stage_108_mathematica.log`, exit_code: 0):
- `pure scale invariance = 0` / `PASS: pure scale invariance` (F3)
- `chiGen independent-route agreement = 0` / `PASS: chiGen independent-route agreement` (F2)
- `Sigma5 locus (Class C) = -sigma0/27 = 0` / `PASS` (F4 Class C)
- `general preservation submanifold = S(1 - beta^5)/9 - sigma0/27 = 0` / `PASS` (F4 anchor intact)
- `general preservation locus check = 0` (F4 demoted print — no PASS banner, confirming demotion)
- Final: `Stage 108 Mathematica audit passed.`

The demoted Class D round-trip printing `0` with NO `PASS` line in the Mathematica
log, and with no `expect_zero` banner formatting in the SymPy log, is the
signature that it is correctly informational rather than load-bearing.

## Material-change assessment

`material_change: true`. The edits are substantive, not cosmetic: F2 adds a real
algorithm-independent second-engine cross-check (`chiGenAlt`); F3 replaces a
self-referential cancellation with an external paper anchor (the literal
fingerprint); F4 replaces two solve-then-resubstitute round-trips that could never
fail with one falsifiable value anchor (Class C) plus a demotion that preserves
the genuine submanifold anchor (Class D). The previously-asserted physics
quantities (`chi_add`, `chi_gen`, the loci) are unchanged in value; the
verification surface is strictly strengthened. No downstream unit's result
changes — this is a verification-hardening pass, not a correction.

## Side observations (non-blocking)

- The original auditor noted a paper-prose label artifact: the stage-108 card
  display title reads `Stage~125` while the `\label` and filename are 108. This is
  outside script scope and not a script finding; flagged here only for
  completeness, no action in this verification.
- Mathematica's `Sigma4(beta=1)` prints in an un-cancelled form (log line 26),
  but the subsequent `Sigma4(beta=1) + sigma0/27 = 0` PASS confirms it reduces to
  `-sigma0/27`; cosmetic display artifact only.

## Verdict justification

All three in-scope findings are `resolved`. F2's new route is genuinely
algorithm-independent — it builds the L-coefficients by hand and solves plain
algebraic equations with no `Series`/`Coefficient`/derivative, sharing only the
`lambdaOut` physical premise with the series route, so the agreement assertion
catches a route-specific coefficient error rather than rubber-stamping the same
algebra; the existing series route is untouched (additive). F3 removes the bare
S-cancellation and anchors `Y_scale` to the literal appendix fingerprint in both
engines, making it falsifiable against a wrong `Lambda_out`. F4 replaces both
solve-then-resubstitute round-trips with a falsifiable Class C value anchor and a
demotion-to-print, leaving the load-bearing submanifold anchor intact. Both
scripts exit 0 with every relevant assertion printing `0`/`PASS` and the demoted
round-trips printing as bare informational lines. F1 is correctly gated and
excluded. Verdict: `verified`.
