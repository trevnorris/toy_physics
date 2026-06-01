---
unit_id: 193
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 193

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy §IV (`scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py:127-162`)
was rewritten exactly as the directive specified. The hand-written `K0`/`D0`-based
quadratic `Deff = D0*I3 - chi**2*Cvec.T*Cvec/K0` is gone. Instead:
- Two distinct block symbols are introduced: `D0scalar` (the eliminated scalar/geometry
  block, line 127) and `D2blk` (the leading l=2 block, line 128) — no longer conflated
  with the §II conservative `D0`.
- An explicit block operator with a **linear** chi coupling is built (lines 137-139):
  `Dblock = [[D0scalar, chi*Cvec],[chi*Cvec.T, D2blk*I3]]`. The printed matrix
  (output lines 70-76) shows off-diagonal entries `c₂₀⋅χ, c₂₁⋅χ, c₂₂⋅χ` — linear in chi,
  no `**2`.
- `Deff` is produced by the **genuine Schur-complement expression** (line 142):
  `D2blk*I3 - (chi*Cvec.T)*(Matrix([[D0scalar]]).inv())*(chi*Cvec)`. The chi² factor is
  *generated* by the Schur algebra contracting two linear-in-chi off-diagonals, not typed
  by hand.
- The load-bearing assertion (lines 153-156) subtracts the Schur output from the chi²
  target `D2blk*I3 - chi**2*Cvec.T*Cvec/D0scalar`; output (lines 107-112) shows the 3×3
  zero matrix.
- The ledger text (line 175) was updated to `D_eff,l=2 = D2 I - chi^2 C^T C / D0scalar`.

The git diff (`stage_193_diff.patch`) confirms the ONLY edits are §IV (lines 124-162) and
the one matching ledger print line. Sections §I, §II, §III are byte-for-byte untouched.

**Assessment:**
The fix is GENUINE. The first new `expect_zero` is non-tautological: it can only pass
because the off-diagonal was linear in chi and the Schur contraction
`(chi C^T)(D0scalar)^{-1}(chi C)` genuinely produces `chi^2 C^T C / D0scalar`. If the
result of the Schur algebra were chi^0 or chi^1, the subtraction against the chi^2 target
would be nonzero and the check would FAIL — confirmed by the explicit non-diagonal cross
terms `-c₂₀⋅c₂₁⋅χ²/D0scalar` etc. in the printed Schur output (output lines 77-91), which
match the rank-1 outer product `C^T C` exactly, not a diagonal-only hand form. The
`d/dchi Deff|_{chi=0}` check is now load-bearing because `Deff` derives its chi-dependence
from the Schur algebra, not from a hand-typed `chi**2`; it vanishes precisely because the
Schur complement is quadratic (the firewall content). No new tautology was introduced; the
three §IV checks now exercise the linear→quadratic Schur theorem rather than restating a
hand-written quadratic. Symbol roles are explicit and disambiguated from §II.

### F2 — missing_mathematica

**Classification:** resolved

**What changed:**
A new independent-engine script was created:
`mathematica/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_mathematica_audit.wl`
(output `.txt` all PASS, exit 0). It re-verifies every load-bearing SymPy assertion via
Mathematica-native primitives through a genuinely different derivation route.

**Assessment:**
This is a GENUINE INDEPENDENT ROUTE, not a transliteration. Concrete evidence of
independence (different decomposition / native primitives, not a line-by-line port):

- §I: SymPy uses scalar formula functions `grouped_trace_anomaly`/`grouped_inverse`. The
  `.wl` instead encodes the map as a 3×3 **matrix** `projector = {{1/5,2/5,2/5},
  {1/5,-1/10,-1/10},{0,1/2,-1/2}}` and obtains the inverse with native `Inverse[projector]`,
  then checks it against the inverse-map target matrix (lines 81-98). Matrix route vs.
  closed-form route — different decomposition.
- §II: This is the strongest independence. SymPy *asserts* the closed forms
  `nu2_common=-D2/D0`, `nu4_common=(D2²-D0·D4)/D0²`. The `.wl` *derives* them from the
  response function `responseShape[w]=d0/(d0+d2 w²+d4 w⁴)` via Taylor derivatives
  `D[..,{w,2}]/2` and `D[..,{w,4}]/24` at `w->0` (lines 107-118), and solves the one-pole
  surface with native `Solve[deltaPole==0, d4]` (line 121) instead of back-substituting a
  hand-chosen `D4_onepole`. Genuinely different generation of the same moments.
- §III: SymPy `sp.series`; the `.wl` uses native `Series`+`Normal` and additionally
  isolates the w² and w⁴ coefficients individually with `Coefficient` (lines 129-138) —
  finer-grained, native, and feeds the §II-derived (not hand-typed) moments.
- §IV: Mirrors the CORRECTED firewall. Block built with native `ArrayFlatten`
  (lines 148-151), Schur complement via native `Inverse[scalarBlock]` (line 152). Adds
  checks SymPy lacks: `D[blockOperator,chi]` equals C and its 2nd derivative vanishes
  (confirming linear coupling, lines 157-158), and the Schur **determinant identity**
  `Det[block] = Det[scalar]·Det[Schur]` (line 160). These strengthen rather than re-type.

No check is vacuous/X−X. Every `expectZero` subtracts an independently computed quantity
from a target (e.g. `nu4FromD - (d2^2-d0 d4)/d0^2` compares a Taylor-derived moment against
the closed form; `carrierCoeff2 - nu2FromD` compares a series coefficient against a
derivative-derived moment). The §IV checks mirror the genuine Schur complement (via
`Inverse`), never a hand-written quadratic. House idioms (`expectZero`, `stripCE`,
`$Assumptions` positivity, `Exit[0]`) follow the established pattern; the math route is the
verifier's own.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Schur complement - (D2 I - chi^2 C^T C / D0scalar) =` → 3×3 zero matrix (output 107-112)
- `d/dchi D_eff at chi=0 (linear-order firewall) =` → 3×3 zero matrix (output 113-118)
- Printed `D_block(chi)` shows linear off-diagonals `c₂₀⋅χ` (output 70-76); printed `Deff`
  shows the chi²/D0scalar rank-1 structure (output 77-91), confirming the Schur algebra ran.

**Mathematica:** exit=0. Notable lines:
- `PASS: Schur complement - (D2 I - chi^2 C^T C / D0scalar)` (output 58-59)
- `PASS: det block - det scalar det Schur` (output 60-61)
- `PASS: d/dchi D_eff at chi=0 (linear-order firewall)` (output 62-63)
- `Delta_pole = -((3*d2^2 + d0*d4)/d0^2)` derived from response derivatives (output 25),
  agreeing with SymPy's `-(3 D2^2 + D0 D4)/D0^2`.
All PASS lines present; no FAIL.

**Output freshness:** Confirmed. SymPy script mtime 11:39, its output `.txt` 11:43.
Mathematica script mtime 11:41, its output `.txt` 11:43. Both outputs newer than their
scripts — regenerated post-fix. Both exec logs end `# exit_code: 0`.

## Material-change assessment

`material_change`: false.

§IV proves the same firewall result the paper already states (`D_eff,l=2 = D2 I -
chi^2 C^T C / D0scalar`, linear-order vanishes); no derived constant or downstream value
changed — only the *strength* of the in-script proof was raised from assertion to
derivation. The new symbols `D0scalar`/`D2blk` are local placeholders disambiguating the
block roles and do not propagate. §I–§III are untouched. No downstream unit depends on a
changed result.

## Side observations (non-blocking)

- The SymPy banner/ledger still read "STAGE 176" while the file is stage 193 (pre-existing,
  noted in the auditor report; cosmetic, not part of any finding). Not blocking.
- Cross-engine variable naming differs (`D0scalar`/`D2blk` in SymPy vs `d0geom`/`d2iso` in
  the `.wl`), which is appropriate for an independent route and reinforces non-transliteration.

## Verdict justification

Both findings are resolved. F1's §IV now builds an explicit block operator with a
demonstrably linear chi coupling (`chi*Cvec`, no `**2`) and produces `Deff` via the genuine
Schur-complement expression `D2blk*I3 - (chi*Cvec.T)*inv(D0scalar)*(chi*Cvec)`; the firewall
assertion passes only because the Schur algebra turns the linear coupling into a chi²
correction and would fail for any chi⁰/chi¹ result — no new tautology, §I–§III untouched.
F2's new `.wl` is a genuinely independent Mathematica-native route (matrix `Inverse` for the
projector, Taylor-derivative derivation of the moments, native `Solve`/`Series`/`Coefficient`,
`ArrayFlatten` block + `Inverse` Schur + determinant identity), not a port, with no vacuous
checks, mirroring the corrected firewall. Cross-engine conclusions agree and both engines
exit 0 with fresh outputs. Verdict: verified.
