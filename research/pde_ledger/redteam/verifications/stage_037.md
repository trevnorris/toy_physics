---
unit_id: 037
batch: III.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 037

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex restructured the `.wl` file at the locations the directive named.

- Section 3 (`mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:95-127`): the hand-supplied `deltaUW`, `xiTerm`, `alphaTerm`, and `sigmaExpected` literals were removed. The Schur block is now derived by (a) computing `sigmaWall = FullSimplify[cMat . LinearSolve[bMat, Transpose[cMat]]]` (line 114, unchanged), then (b) solving for `alphaSolved` from the off-diagonal `sigmaWall[[1,2]] == alpha*kappa0*kappa1` (lines 116-119) and for `xiSolved` from the on-diagonal `sigmaWall[[1,1]] == xi + alphaSolved*kappa0^2` (lines 120-123). The substantive cross-check is `expectZero["Sigma_wall (2,2) consistency with ansatz", sigmaWall[[2,2]] - xiSolved - alphaSolved*kappa1^2]` at line 125. `xiSolved` and `alphaSolved` are printed at lines 126-127.
- Section 4 (`mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:129-195`): the prior `aExpected` assignment is replaced with `aDerived = FullSimplify[Together[k0 - gUCont^2/omegaU2], ...]` (line 150), and `deltaExpected` is replaced with `deltaDerived = FullSimplify[Together[deltaKAx/aDerived], ...]` (line 165). Two new substantive identities are asserted: `"A numerator matches Schur form"` (lines 168-174) tests `Numerator[Together[aDerived]]*(muEta*kU) - Denominator[Together[aDerived]]*(kU*kEtaEff - cEtaU^2) == 0`, and `"delta numerator matches closed form"` (lines 180-187) tests the analogous identity for `delta`. The sanity step `"A reduces to closed form" (a - aDerived)` is kept at line 167. The print block at lines 189-195 uses the derived/expected variables.
- Sections 1, 2, 5 are unchanged (consistent with directive item 3). `Chi`, `beta0`, `alphaMix`, `mMix` retain their `*Expected` checks (consistent with directive item 4).

**Assessment:**
The edit matches the directive exactly. The forbidden literals (`aExpected`, `deltaExpected`, `xiTerm`, `alphaTerm`, `sigmaExpected`, `deltaUW`) are absent from the file (grep returns empty). The new Schur check is non-tautological: `xiSolved` and `alphaSolved` are recovered from only two entries of `sigmaWall`, and the `(2,2)` entry must then independently agree with `xiSolved + alphaSolved*kappa1^2` for the ansatz `xi I + alpha v v^T` to be the correct factorization — that is a true cross-check, not a re-statement. The recovered closed form printed in the output (`Xi (recovered) = gU^2/aU`) matches the SymPy `Xi = gU^2/A_U`, and the `alpha (recovered)` numerator/denominator expansion is the same rational form as the SymPy `alpha` (after `Together`). The new Section-4 identities are also non-tautological: they test that `aDerived` and `deltaDerived` reduce algebraically to the SymPy closed-form `(kU*kEtaEff - cEtaU^2)/(muEta*kU)` and `(kU*Pi^2*tw)/(ell^2*(kU*kEtaEff - cEtaU^2))` by extracting numerator/denominator from the `Together` reduction and asserting cross-multiplied equality, rather than by setting that closed form as a hand-supplied target. No collateral edits beyond the F1 scope appear in the diff.

## Exec log assessment

**SymPy:** exit=0 (inferred from the post-fix transcript at `scripts/output/moving_throat_pde_stage037_continuum_kernel_sympy_audit.txt`; the fix loop halts on FAIL/Traceback). Notable lines:

```
kappa0 = 2*sqrt(2)/pi
sigma  = 88/(9*pi**2)
Sigma_wall - [Xi I + alpha v v^T] =
⎡0  0⎤
⎣0  0⎦
A continuum formula = 0
delta continuum formula = 0
```

SymPy is unchanged by this directive (F1 only targets the `.wl` mirror); the script still runs to completion and emits the same closed forms as before.

**Mathematica:** exit=0 (inferred from `mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt:93`, which ends with `Stage 037 Mathematica audit passed.`). Notable lines:

```
Sigma_wall (2,2) consistency with ansatz = 0
PASS: Sigma_wall (2,2) consistency with ansatz
Xi (recovered) = gU^2/aU
alpha (recovered) = (-88*aU*gB^2*gR^2 + 9*(aPhi*gR^2*gU^2 + aU*(aU*aW*gB^2 + 2*aPhi*gR*gU*gW + aPhi*aU*gW^2))*Pi^2)/(aPhi*aU*(-88*gR^2 + 9*aU*aW*Pi^2))
A numerator matches Schur form = 0
PASS: A numerator matches Schur form
delta numerator matches closed form = 0
PASS: delta numerator matches closed form
```

Every new and retained assertion prints `= 0` and `PASS`. `Xi (recovered) = gU^2/aU` is identical to the SymPy `Xi = gU^2/A_U`. The `alpha (recovered)` rational form expands the SymPy `alpha = gB^2/A_phi + (A_U gW + gR gU)^2/(A_U Delta_UW)` over a common denominator with `Delta_UW = A_U A_W - gR^2 sigma = (9*aU*aW*Pi^2 - 88*gR^2)/(9*Pi^2)` after substituting `sigma = 88/(9*Pi^2)`. Engines agree.

**Output freshness:** `mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl` mtime is 2026-05-22 12:16; `mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt` mtime is 2026-05-22 12:18 (newer, fresh). `scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py` mtime is 2026-04-01 12:39; `scripts/output/moving_throat_pde_stage037_continuum_kernel_sympy_audit.txt` mtime is 2026-05-22 12:17 (newer, fresh).

## Material-change assessment

`material_change`: false.

The change is restricted to how the Mathematica script *derives* the same closed forms; the headline quantities (`kappa0`, `kappa1`, `sigma`, `Xi`, `alpha`, `A`, `Delta0`, `Chi`, `beta0`, `alphaMix`, `mMix`, `delta`) print identical values before and after, and the SymPy engine was not touched. No downstream-visible result has changed; downstream units that consume these closed forms see the same numbers.

## Side observations (non-blocking)

- The directive instructed Codex to keep the existing `Chi`, `beta0`, `alphaMix`, `mMix` `*Expected` checks. Codex did keep them, but these four targets remain hand-supplied closed forms (lines 152-164), so the "transliteration" pattern is now *partially* broken (Section 3 + the `A`/`delta` half of Section 4 are independently derived) rather than fully broken. This is exactly what the directive asked for — the directive explicitly stated "once `A` and `delta` are derived independently, the rest follow" — so it does not block verification. It is noted only for context.
- The variable `v = {{kappa0}, {kappa1}}` is defined at line 115 but never used after that point (the verifier check at line 125 indexes `sigmaWall` directly). Harmless dead code.

## Verdict justification

The single finding F1 is fully resolved. The hand-supplied Schur targets (`xiTerm`, `alphaTerm`, `sigmaExpected`) and the hand-supplied `aExpected`/`deltaExpected` literals are gone; in their place the Mathematica script (i) recovers `xi` and `alpha` by solving against two entries of `sigmaWall` and cross-checks the third, and (ii) derives the `A` and `delta` closed forms via `Together` reductions and proves equality with the Schur closed form by a numerator/denominator cross-multiplication identity. Both substantive new assertions print `= 0`, the script exits 0, and the saved output is newer than the script source. No regressions in the diff; SymPy untouched.
