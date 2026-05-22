---
unit_id: 044
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T13:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 5
findings_total: 5
material_change: false
---

# Verification — unit 044

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:63-77`, Codex added the directive's block verbatim: `xiRoots = xi /. Solve[branchEq == 0, xi]`, then `SelectFirst` with a `FullSimplify === 0` test for the zero-load branch, the documented `MissingQ` fallback using `Select` + `Simplify`, a `Print["xi_phys (from Solve) = ", ...]`, and `expectZero["xiPhys solve match", xiPhysSolve - xiPhys]`. Inserted after the existing `zero-load root` check (line 61) and before `subbanner["2. Exact ..."]` (now line 79), as instructed. The pre-existing manual block (lines 42-61) was preserved (augmented, not replaced).

**Assessment:**
Edit matches the directive exactly. The new check is non-tautological: `Solve` enumerates the algebraic roots of `branchEq`, and the directive selects whichever root satisfies the zero-load criterion `(# /. {mMix->0, mSupp->0}) === 0`. This is an algebraically distinct path from the manual `(-bCont + Sqrt[deltaDisc])/2` construction (`Solve` would use companion-matrix / closed-form coefficient extraction internally; both must agree on the same physical root). The Mathematica output line 16 shows `xi_phys (from Solve)` printed with a numerically different rearrangement of the same root: the manual form factors a common `9` out of the sqrt, while the Solve form distributes it (`Sqrt[81*delta^2 + ...]/18`). Output line 17-18: `xiPhys solve match = 0 / PASS`. No collateral edits.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:148-154`: The previous `mismatch = sp.symbols("Delta_R", real=True)` line plus `C_rewritten = C_cont.subs(Rphi, RU - mismatch)` and `expect_zero("quadratic mismatch penalty", ...)` are gone. Replaced with `C_from_branch_eq = sp.simplify(sp.Poly(branch_eq, xi).nth(0) / 9)` against `C_expected = -delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU-Rphi)**2`, asserted as `expect_zero("mismatch penalty in C coefficient", ...)`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:140-148`: Parallel change using `CoefficientList[branchEq, xi][[1]]/9`. The previous `cRewritten = cCont /. rPhi -> rU - deltaR` block removed.
- WL line 35: `Clear[xi, delta, mMix, mSupp, rU, rPhi];` — `deltaR` removed.
- WL line 37: `Element[{xi, delta, mMix, mSupp, rU, rPhi}, Reals]` — `deltaR` removed.

**Assessment:**
Edit matches the directive. The new assertion is non-tautological: `branch_eq` was independently constructed at line 65 as `Numerator[Together[n_req - Msupp]]`, so extracting its `xi`-constant coefficient and comparing to the manually-written `C_expected` actually tests whether `C_cont` has the claimed mismatch-penalty form. If `C_cont` had been written with a sign error in `(RU-Rphi)^2`, the residual would be nonzero. Sympy output line 137 shows `mismatch penalty in C coefficient = 0`; Mathematica output line 48-49: `mismatch penalty in C coefficient = 0 / PASS`. No collateral edits.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`: lines 127-129 (the `# Setting actual support baseline ...` comment, `tracking_equation = sp.simplify(...)`, and `expect_zero("tracking total-loading law", tracking_equation)`) are deleted. The remaining `27.4` block keeps `n_track`, `G_q`, `tracking collapse of n_req`, then `F_track`, etc.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`: the `trackingEquation = FullSimplify[...]` line is deleted and the `expectZero["tracking total-loading law", trackingEquation]` assertion is deleted. Surrounding `nTrack`, `gQ`, `fTrack`, `fTrackExpected`, and the `tracking collapse of n_req` / `tracking F collapse` asserts remain.

**Assessment:**
Matches the directive. The redundant check (algebraically equal to the preceding `tracking collapse` assertion) is removed. Output for sympy section 27.4 shows just `tracking collapse of n_req = 0` and `tracking F collapse = 0`; the Mathematica output similarly shows only those two PASS lines under section 4. No phantom of `tracking total-loading law` remains in either transcript.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
- Sympy lines 97-109: new `Rphi_lit = sp.Integer(2)` block inserted after `print("F_cont built successfully.")` and before `subbanner("27.3 …")`. `F_lit = sp.simplify(F_cont.subs(Rphi, Rphi_lit))` compared against the manually-transcribed `F_lit_expected` template at `Rphi_lit = 2`, then `expect_zero("third-slice F at Rphi=2", F_lit - F_lit_expected)`.
- WL lines 94-103: parallel block with `fLit = FullSimplify[fCont /. rPhi -> 2, ...]` and `fLitExpected` written out with literal 2 substituted into the same template. `expectZero["third-slice F at rPhi=2", fLit - fLitExpected]`.

**Assessment:**
Matches the directive exactly. The original auditor's report notes (self-test paragraph) that this is a sanity check — `F_lit_expected` is the same template evaluated at `Rphi=2`, so the residual is identically 0 unless `F_cont`'s definition itself contains an `Rphi`-dependent typo that fails to survive substitution. The non-trivial content is that `F_lit_expected` is an independent textual transcription of the template with the literal `2` plugged in — a coefficient mistype in `F_cont` would break this. Sympy output line 90: `third-slice F at Rphi=2 = 0`. Mathematica output line 24-25: `third-slice F at rPhi=2 = 0 / PASS`. Both engines.

### F5 — symbol_assumption_error

**Classification:** resolved

**What changed:**
- Sympy: line 50 (`sigma0 = sp.symbols("sigma_0", real=True)`) deleted entirely.
- Sympy line 13 (docstring): now reads `4. The minimal-kernel limit R_phi=1 gives the exact source-tied closure.` — directive's required change applied.
- WL line 35: `Clear[xi, delta, mMix, mSupp, rU, rPhi];` (both `sigma0` and `deltaR` removed; the `deltaR` removal coordinates with F2).
- WL line 37: `Element[{xi, delta, mMix, mSupp, rU, rPhi}, Reals]` (same).

**Assessment:**
Matches the directive. `grep -i 'sigma\|deltaR\|Delta_R'` against both files returns nothing. The docstring now accurately describes the substitution that the script actually performs (`Rphi -> 1`). The directive permitted touching the SymPy docstring (this is the script's own header, not a prose document under paper.tex/notes/), so this is in scope. Both engines still exit 0 since `sigma0` was unused.

## Exec log assessment

**SymPy:** exit=0 (inferred from clean transcript with no traceback / AssertionError). Notable lines from `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.txt`:

- L21: `quadratic branch equation = 0`
- L52: `zero-load root = 0`
- L90: `third-slice F at Rphi=2 = 0`
- L125-126: `source-tied n = 0` / `source-tied F = 0`
- L131-132: `tracking collapse of n_req = 0` / `tracking F collapse = 0`
- L137: `mismatch penalty in C coefficient = 0`
- No occurrence of `tracking total-loading law` or `quadratic mismatch penalty` (the deleted asserts).

**Mathematica:** exit=0 (transcript ends with `Stage 044 Mathematica audit passed.`; the script's `Exit[0]` is reached). Notable lines from `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.txt`:

- L12: `PASS: quadratic branch equation`
- L15: `PASS: zero-load root`
- L17-18: `xi_phys (from Solve) = ...` / `xiPhys solve match = 0 / PASS: xiPhys solve match`
- L24-25: `third-slice F at rPhi=2 = 0 / PASS`
- L33-35: `PASS: source-tied n` / `PASS: source-tied F`
- L41, L43: `PASS: tracking collapse of n_req` / `PASS: tracking F collapse`
- L48-49: `mismatch penalty in C coefficient = 0 / PASS`
- L51: `Stage 044 Mathematica audit passed.`
- No occurrence of `tracking total-loading law` or `quadratic mismatch penalty`.

**Output freshness:** Confirmed.
- SymPy: script mtime `2026-05-22 12:41:49`, output mtime `2026-05-22 12:45:51` (4 minutes newer → fresh).
- Mathematica: script mtime `2026-05-22 12:43:12`, output mtime `2026-05-22 12:46:08` (3 minutes newer → fresh).

## Material-change assessment

`material_change`: false.

Rationale: none of the five edits altered a derived numerical or symbolic result that downstream units consume. F1 added a redundant cross-check on `xiPhys` (the value is unchanged; an additional independent derivation matches it). F2 replaced a tautological assertion with a substantive one that still passes against the unchanged `C_cont`. F3 deleted a redundant assertion. F4 added a sanity-check assertion at `Rphi=2`. F5 removed a dead symbol declaration and corrected a docstring. The exported algebraic outputs (`n_req^(cont)`, `xi_phys`, `F_cont`, `n_source`, `F_source`, `F_track`, `B_cont`, `C_cont`) are unchanged in form. No constant value, sign, or coefficient was mutated. Downstream units 045+ that depend on these forms see no algebraic change.

## Side observations (non-blocking)

- The new `Print["xi_phys (from Solve) = ", fmt[xiPhysSolve]]` line uses a different printed form than the manual `xi_phys` (factor-of-9 redistribution between the leading rational and the sqrt). Both forms are algebraically equal, as evidenced by `xiPhys solve match = 0`. Non-issue.
- The sympy `t = sp.symbols("t", real=True)` declaration at line 49 is still in place but `t` is now unused in the script body (the docstring references `q = t R_U`, `r = t R_phi`, but the script uses `q = sqrt(lambda0) * RU` and `r = sqrt(lambda0) * Rphi` at lines 56-57; the local `t` symbol declaration is dead). This is identical to the `sigma0` issue F5 raised, but for a different symbol. NOT a new finding from this verifier — flagging only as a non-blocking observation.

## Verdict justification

All five findings are resolved with edits that match the directive precisely (no deviations, no collateral edits). Both engines re-ran post-fix and produced fresh outputs (mtimes verified newer than scripts) with every assertion reporting `= 0` / `PASS` and no `FAIL`/Traceback/`$Failed` lines. The new F1 `xiPhys solve match` provides genuine algebraic-route independence between SymPy and Mathematica; the new F2 `mismatch penalty in C coefficient` non-tautologically anchors docstring claim 6 by extracting `C` from the independently-built `branch_eq`; F3's deletion removes the redundant tautology; F4's literal third slice at `Rphi=2` adds bivariate constraint coverage; F5 cleans up dead scaffolding. No regressions visible in the diff. Verdict: verified.
