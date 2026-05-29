---
unit_id: 112
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T11:50:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 112

Note: this supersedes the earlier 2-finding verification (claude-opus-4-7, 2026-05-27).
The directive was subsequently amended to fold in F3 (the sigma_W != 0 qualifier from
Codex review R1) and the F2 `.wl` block was revised to add `bDef`, the `presCond`
factorization, the `Reduce[... && sigma != 0]` branch isolation, and the degenerate
sigma=0 case. This verifies the current 3-finding state (F1 label / F2 independent route /
F3 sigma!=0 qualifier).

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script; resolved direction (a), label only)

**Classification:** resolved

**What changed:**
The stale "Stage 95 / STAGE 095" identifier strings were re-anchored to the canonical
internal "Stage 112":
- `scripts/...sympy_audit.py:3` → `Stage 112 SymPy audit.`
- `scripts/...sympy_audit.py:66` → `print('stage112: PASS')` (was `.py:54`; the F3
  insertion shifted line numbers down)
- `mathematica/...mathematica_audit.wl:26` → `banner["STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW"];`
- `mathematica/...mathematica_audit.wl:90` already read `Stage 112 Mathematica audit passed.` and was left untouched per direction (a).

**Assessment:**
Correct and complete. Direct inspection of the current files (grep) shows no residual
`Stage 95`/`STAGE 095`/`Stage 129` strings; all four identifier strings read "Stage 112".
The exec-log banner prints `STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW` and the SymPy
log ends `stage112: PASS`, so both runtime transcripts are now self-consistent. No
paper/notes edit — direction (a) was an APPROVED override of the paper_misalignment hold
(Claude+Codex consult 019e748e). These label edits do not appear in `stage_112_diff.patch`
because that diff was captured against a baseline in which F1 had already been applied (the
`## Applied: F1` block predates the F2/F3 iteration); their current presence is verified
directly from the working tree.

### F2 — mathematica_transliteration (independent route)

**Classification:** resolved

**What changed:**
`.wl:55-71` inserts — ALONGSIDE the still-present chi_Q block (`.wl:75-87`) — an
independent Stage-92 linearized cross-check deriving `gamma_W = 1/9` from the deformation
data `(b, a0, a5)`:
- `bDef = 9*(-L2/L0) - 1` on solB (canonical-even consistency via the z^2/z^0 ratio),
- `a0Def = (z^0 of lambdaHyb on solB) - (z^0 of lambdaOut)` (the z^0 deformation difference),
- `a5Def = (z^5/I of lambdaHyb on solB) - (z^5/I of lambdaOut)` (the z^5 deformation difference),
then `presCond = a0Def/3 + 9 a5Def`, checked to factor as `sigma*(1 - 9 gamma)`, and
`Reduce[presCond == 0 && sigma != 0, gamma, Reals]` → `gamma = 1/9`.

**Assessment — GENUINELY independent (complementary path, NOT a silent restatement):**

1. **Different intermediates.** The chi route uses only the z^0 and z^5 coefficients and
   forms the single rational ratio `chi = (-l5/l0)/(1/27)`, comparing it to
   `(1-9 sigma gamma)/(1-sigma)`. The linearized route instead extracts three separate,
   individually-checked deformation data — `a0` (a z^0 DIFFERENCE `lambdaHyb - lambdaOut`,
   never formed by the chi route), `a5` (a z^5 difference), and `b` (a z^2/z^0 ratio, an
   intermediate the chi route never touches) — and combines them through a DIFFERENT
   functional `a0/3 + 9 a5`. The only shared primitive is series/coefficient extraction
   from the same Lambda_hyb, which is unavoidable; the path from coefficients to
   `gamma_W = 1/9` is structurally distinct.

2. **The algebraic equivalence is exactly what makes it a valid cross-check, not a
   tautology.** Yes, `a0/3 + 9 a5 = sigma(1-9 gamma)` is equivalent to
   `chi_B - 1 = 0 ⟺ sigma(1-9 gamma) = 0`. But equivalence of the FINAL RESULT is the
   point of cross-route verification — two independent paths should land the same answer.
   The tautology test is whether they share the load-bearing computation, and they do not:
   the linear route reaches `sigma(1-9 gamma)` by contracting the deformation vector
   (`a0/3 + 9 a5`), whereas the chi route reaches it via the rational `-l5/l0` scaled by 27
   and minus 1. A transliteration bug in the shared chi choreography (sign slip in the `/I`
   normalization, miscopied `1/27` scale, wrong branch substitution) would NOT propagate
   into `a0Def`/`a5Def`, because those are direct coefficient differences against
   `lambdaOut` with their own normalizations — so the linear route would diverge and one of
   its `expectZero` checks would fail. That is real complementary content.

3. **a0/a5/b extraction checks are NON-TAUTOLOGICAL.** Each is extracted-vs-normalization,
   not identity-by-construction. On `solB = {rho->4 sigma, kappa->1/3}` the z^0 coeff of
   lambdaHyb is `-3 + 3 sigma`, the z^2 coeff `(1-sigma)/3`, the z^5/I coeff `1/9 - gamma sigma`:
   - `a0Def - 3*sigma`: a0Def = `(-3 + 3 sigma) - (-3) = 3 sigma` → check = 0; would FAIL if
     solB carried a wrong rho or the lambdaOut constant were mis-stated.
   - `a5Def + sigma*gamma`: a5Def = `(1/9 - gamma sigma) - 1/9 = -gamma sigma` → check = 0;
     would FAIL on a sign error or wrong gamma coupling in z^5.
   - `bDef = 9*(-L2/L0) - 1` = `9*(1/9) - 1 = 0` only because the canonical-even z^2/z^0
     ratio holds on solB; a wrong kappa breaks it.
   These match the notes' boxed `(b,a0,a5)=(0, 3 sigma_W, -sigma_W gamma_W)` and condition
   `a0/3 + 9 a5 = sigma_W(1-9 gamma_W) = 0` (notes lines 96-106), so the prescribed targets
   are correct, not invented to pass.

4. **Existing chi_Q branch checks REMAIN.** Directive said "alongside, not replacing" —
   confirmed: `.wl:75-87` still computes `chiA`/`chiB` and runs all four original chi checks
   (branch A, branch B formula, branch B at gamma=1/9, scaled identity), unchanged.

### F3 — symbol_assumption_error (sigma_W != 0 qualifier)

**Classification:** resolved

**What changed:**
- `.wl` side (inside the F2 block): `presCond` factorizes to `sigma*(1-9 gamma)`
  (`.wl:67-68`, COMPUTED and checked), the iff is closed only on the nontrivial branch via
  `Reduce[presCond == 0 && sigma != 0, gamma, Reals]` (`.wl:69-71`), and the degenerate
  case is `(chiBgen /. sigma -> 0) - 1` (`.wl:72-73`) where chiBgen is symbolic in gamma.
- `.py` side (`.py:55-61`): adds `chi_B_general = (1-9*sigma*gamma)/(1-sigma)`, a comment
  stating preservation holds iff `sigma*(1-9 gamma)=0`, the factorization check
  `together(chi_B_general - 1).as_numer_denom()[0] - sigma*(1-9 gamma)` (numerator of
  chi_B-1 = `sigma(1-9 gamma)` up to SymPy sign normalization, landing at 0 per the log),
  and the degenerate check `chi_B_general.subs(sigma, 0) - 1`.

**Assessment:**
Both real, not asserted by construction. The factorization genuinely shows
`chi_B - 1 ⟺ sigma(1-9 gamma) = 0` by extracting and comparing the numerator (a wrong
factorization leaves a nonzero residual and fails). The degenerate assertion is real: at
sigma=0, `chi_B_general` is `1/1 = 1` INDEPENDENTLY of the still-symbolic gamma, genuinely
exercising gamma-independence rather than substituting gamma to a fixed value. This
directly closes Codex R1, which complained that the original
`Solve[a0/3 + 9 a5 == 0, gamma]` proved only the generic branch while the declared
assumptions still permitted sigma=0; the `sigma != 0` restriction is now explicit in the
Reduce and the sigma=0 case is handled separately on both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `chi_B - 1 factorizes as sigma(1-9 gamma) = 0`
- `degenerate sigma=0: chi_B = 1 for any gamma = 0`
- `scaled identity on branch B = 0`
- `stage112: PASS`

**Mathematica:** exit=0, 14 PASS / 0 FAIL. All 14 are content-bearing (none is an
identity-by-construction `0`):
1. branch A rho - sigma; 2. branch A kappa; 3. branch B rho - 4 sigma;
4. branch B kappa - 1/3; 5. independent: b = 0 on solB;
6. independent: a0 - 3 sigma on solB; 7. independent: a5 + sigma*gamma on solB;
8. independent: preservation cond = sigma (1 - 9 gamma);
9. independent: gamma_W = 1/9 on nontrivial branch (sigma != 0);
10. degenerate sigma=0: chi_B = 1 for any gamma;
11. chi_Q branch A - (1 - 9 sigma gamma); 12. chi_Q branch B - (1 - 9 sigma gamma)/(1 - sigma);
13. chi_Q branch B at gamma=1/9; 14. scaled identity on branch B.
Banner correctly reads `STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW`.

**Output freshness:** confirmed. Script mtimes (sympy 11:28:29, wl 11:31:27) precede the
saved `.txt` outputs (sympy out 11:35:04, wl out 11:35:37), so the committed outputs were
re-generated post-fix.

## Material-change assessment

`material_change`: false. No derived result changed. All edits are (a) cosmetic identifier
strings (F1), (b) an additive independent cross-check re-deriving the already-established
`gamma_W = 1/9` (F2), and (c) additive honesty qualifiers around the already-true
`chi_B = 1` claim plus its degenerate exception (F3). The five original paper deliverables
(two-branch solve, branch-B chi_Q, preservation at gamma=1/9, collapse identity) are
computed exactly as before; no downstream unit's carried value is altered.

## Side observations (non-blocking)

- The new `.wl` comment calls the route "Stage-92 linearized" (dropping the older
  "(= stage 109)" parenthetical); this matches the notes box wording (line 96, "Stage-92
  deformation data"). Cosmetic.
- `a0Def` peels the Lambda_out constant via `Coefficient[lambdaOut, z, 0]` rather than the
  literal `-3` used in the earlier draft — slightly more robust to a future lambdaOut edit.
  Non-blocking.

## Verdict justification

All three findings are resolved. F1's label re-anchoring to "Stage 112" is present in all
four identifier strings across both engines and confirmed in the runtime transcripts. F2's
linearized `(b, a0, a5)` route is a genuine complementary computational path, not a silent
restatement: it uses distinct intermediates (z^0 and z^5 deformation differences against
lambdaOut, plus a z^2/z^0 ratio) the chi route never forms, its a0/a5/b extraction checks
are non-tautological (each would fail under a wrong branch, sign, or normalization), and the
existing chi_Q block is retained alongside it. F3 closes Codex R1 by making the `sigma != 0`
restriction explicit (Reduce) and handling the degenerate sigma=0 case separately on both
engines, with the factorization computed and checked rather than asserted. Both exec logs
exit 0 (SymPy clean; Mathematica 14/14 PASS, all content-bearing), outputs are fresh, and no
regressions appear in the diff. Verdict: verified.
