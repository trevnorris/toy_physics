---
unit_id: 170
batch: V.1
verifier_model: claude-opus-4-8
verify_date: 2026-05-29T00:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 170

(This round verifies the single `tautological_check` finding raised in the current
`reports/stage_170.md` against the Section-5 weak-axisymmetric block. An earlier verification of a
different finding set — paper_misalignment / transliteration / banner — has been superseded by this
file.)

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Codex deleted the re-typed helper defs in both engines and routed the Section-5 lane inputs through
the derived Section-2 map objects.

- SymPy (`scripts/...stage170...sympy_audit.py`): deleted `kappa_map`/`gamma_map` (formerly
  L144-148). The lane loop (now L151-158) computes
  `dkW[A] = dkappa_from_du2.subs({dD2: eps_l*lam*D2_1, dD0: eps_l*lam*D0_1})` and
  `dgW[A] = sp.simplify(dgamma_from_dP0.subs({dN0: eps_l*lam*N0_1, dD0: eps_l*lam*D0_1}).subs(P0, N0/D0))`.
  The gamma amplitude assertion (L157-158) now reconciles the target via `gamma1.subs(P0, N0/D0)`.
  `kappa1`/`gamma1` (L144-145) are unchanged paper targets. Signature-ratio checks (L160-163)
  unchanged.
- Mathematica (`mathematica/...stage170...mathematica_audit.wl`): deleted `kappaMap`/`gammaMap`
  (formerly L138-139). The six lane assignments (now L140-145) substitute lane-scaled inputs into
  `dkappaFromdu2`/`dgammaFromdP0`, with `P0 -> N0/D0` and `FullSimplify` on the gamma side. The
  three gamma amplitude assertions (L151-153) now use `gamma1 /. P0 -> N0/D0`. Kappa amplitude
  assertions (L148-150) and signature ratios (L154-157) unchanged.

**Assessment:**
The edit matches the directive exactly and addresses the finding.

(a) The re-typed helpers no longer exist. `grep -E "kappa_map|gamma_map|kappaMap|gammaMap"` across
both files returns no matches (exit 1).

(b) The lane checks now reference the derived map objects: SymPy uses
`dkappa_from_du2`/`dgamma_from_dP0`; Mathematica uses `dkappaFromdu2`/`dgammaFromdP0`.

(c) The checks still pass. Both exec logs exit 0 with every Section-5 line present (sympy log
L49-58, mathematica log L60-79 all PASS / `= 0`).

(d) The check is now genuine, not tautological. `dkappa_from_du2` (sympy L69) is a *derived* object:
`sp.solve(sp.Eq(du2sym, du2_hyb), dkappa)[0].subs(du2sym, du2)`, where `du2` (L52) is the
first-order series coefficient of the genuine definition `u2_full = -(D2+eps*dD2)/(D0+eps*dD0)`
(L48). `kappa1` (L144) is the separately hand-written paper closed form. The amplitude assertion
`dkW[A] - eps_l*lam*kappa1` therefore now flows the LHS through the Section-2 derivation while the
RHS is the independent target — a sign or coefficient error injected into the series expansion or
the hybrid inversion would make the residual nonzero. The Mathematica side is the same structure
with an independent first-order primitive (`D[u2Full, eps] /. eps -> 0`, wl L54) and native `Solve`
inversion (wl L67-68). The gamma side correctly applies `P0 -> N0/D0` to both the derived-map output
and the `gamma1` target, mirroring the Section-2 reconciliation (sympy L79 / wl L83), so the
comparison is apples-to-apples rather than spuriously passing.

No collateral edits. The diff (`stage_170_diff.patch`) is confined to the helper deletions, the
lane-assignment replacements, and the gamma-target reconciliations the directive named. Kappa
amplitude assertions and all four signature ratios are byte-for-byte unchanged, exactly as
instructed. The signature-ratio checks now test linearity of the *derived* map, which the directive
explicitly accepts as an acceptable secondary check.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L49 `delta kappa_W^(20) - eps lambda kappa1 = 0`
- L50 `delta gamma_W^(20) - eps lambda gamma1 = 0`
- L55-58 signature ratios all `= 0`
- L67 `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- L60-61 `delta kappa_W^(20) - eps kappa1 = 0` / `PASS`
- L66-67 `delta gamma_W^(20) - eps gamma1 = 0` / `PASS`
- L72-79 four signature-ratio `PASS` lines
- L88 `# exit_code: 0`

**Output freshness:** confirmed. The saved `.txt` outputs (sympy and mathematica) carry mtime
2026-05-29 00:05:18, after the exec runs timestamped 00:03:18 / 00:03:22 in the logs. The sympy
`.txt` Section-5 block (L44-53) matches the exec-log content.

## Material-change assessment

`material_change`: false. The fix only strengthens the Section-5 verification path (routing through
the derived map instead of a re-typed copy). No derived constant, transport law, outlet-map formula,
or carry-forward result changed value — every printed residual remains `0`/PASS and the
carry-forward formulas at the end of both outputs are identical to pre-fix. Downstream units
depending on Stage 170's outlet map (K_A, G_A) see no numeric or symbolic change.

## Side observations (non-blocking)

- The SymPy Section-5 log labels (lines 49-54) print the literal token "lambda" for all three
  lanes, so the 20/21/22 amplitude-check lines are visually indistinguishable in the transcript even
  though the assertions differ (correct `lam` substituted each iteration). Purely cosmetic; the
  Mathematica side labels them distinctly. Not raised as a finding.
- The signature-ratio checks now verify linearity of the derived map — a legitimate but weaker
  secondary check; the load-bearing amplitude checks are the genuine ones.

## Verdict justification

The single finding F1 is fully resolved: the re-typed `kappa_map`/`gamma_map`/`kappaMap`/`gammaMap`
helpers are gone from both engines, the Section-5 lane checks now substitute lane-scaled inputs into
the derived Section-2 maps (`dkappa_from_du2`/`dgamma_from_dP0`,
`dkappaFromdu2`/`dgammaFromdP0`) while comparing against the independently-written `kappa1`/`gamma1`
targets with correct `P0 -> N0/D0` reconciliation, both exec logs exit 0 with every Section-5 PASS
line intact, the saved outputs were regenerated post-fix, and the diff contains no collateral edits
or regressions. The check is now genuine rather than tautological. Verdict: verified.
