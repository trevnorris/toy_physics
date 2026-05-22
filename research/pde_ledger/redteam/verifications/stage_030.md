---
unit_id: 030
batch: II.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 030

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:113-117` — the assertion `expect_zero("Gamma5_sel - G5*P0_sel", Gamma5_sel - G5 * P0_sel)` is removed and replaced with the directive's exact explanatory comment block. The remaining `expect_zero("P0_sel + beta0*d(log lambda_-)/d alpha", ...)` at lines 118-121 is preserved.
- `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:102-107` — the assertion `expectZero["Gamma5_sel - G5*P0_sel", gamma5Sel - g5*p0Sel]` is removed and replaced with the directive's exact `(* Note: ... *)` comment block. The downstream `P0_sel + beta0*d(log lambda_-)/d alpha` check at line 107 is preserved.

**Assessment:**
Edit matches the directive verbatim. No collateral edits in the surrounding lines (definitions of `C5_sel`, `Gamma5_sel`, `P0_sel` are unchanged, so the surviving Part-I generic identity check `Gamma5 - C5/D0` at sympy:55 / wl:41 retains its original substantive content — that check is non-tautological because it exercises the series expansion of `Y_-(omega)` against a hand-written closed form). The saved sympy output no longer contains `Gamma5_sel - G5*P0_sel = 0`, and the saved Mathematica output no longer contains `PASS: Gamma5_sel - G5*P0_sel`. Removal is correct; no residual claim is made.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:147-157` — the single `expect_zero("spectral condition rewrite", ...)` call is removed and replaced with the directive's exact multi-line comment block. The `lambda_req = sp.simplify(...)` binding (line 147) and the `print("lambda_req =")` / `sp.pprint(lambda_req)` (lines 156-157) are preserved.
- `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:126-131` — the `expectZero["spectral condition rewrite", ...]` block is removed and replaced with the directive's exact `(* Note: ... *)` comment. The `lambdaReq = FullSimplify[...]` binding (line 126) and `Print["lambda_req = ", fmt[lambdaReq]]` (line 131) are preserved.

**Assessment:**
Edit matches the directive verbatim. `lambda_req` is still computed and printed, just no longer asserted against a tautological residual. Saved sympy output no longer contains `spectral condition rewrite = 0`; saved Mathematica output no longer contains `PASS: spectral condition rewrite`. The `lambda_req` value still prints in both outputs (sympy lines 107-122; Mathematica line 40) at the expected canonical form `5*beta0*c^5*mhat^2*a^5*(s_-)/108/G/cs^5` (up to symbol renaming `a -> radius`, `c -> cSpeed` in the WL). No false-confidence claim remains.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:48-73` — Part II is restructured: a 2x2 matrix `mMat = {{a + dK - alpha*x1, -alpha*Sqrt[x0*x1]}, {-alpha*Sqrt[x0*x1], a - alpha*x0}}` is defined (lines 55-56), `eigVals = Eigenvalues[mMat]` is computed (line 63), and `lamMinus` / `lamPlus` are bound to `First[Select[eigVals, ...]]` results (lines 64-73). The literal closed-form expression `(2*a + dK - alpha*sigma - r)/2` no longer appears as an *assignment* to `lamMinus` / `lamPlus`; it appears only inside `Select`'s predicate as a label/disambiguator to pick which of Mathematica's two returned eigenvalues is the "minus" branch.
- `sMinusHF`, `sMinusClosed` (lines 78-82) and the existing overlap / weak-loading checks (lines 83-85) are preserved verbatim, as the directive instructed.

**Assessment:**
Edit matches the directive verbatim. The matrix trace and determinant — checked independently by me against the directive's algebra — give:
- `Tr[mMat] = (a + dK - alpha*x1) + (a - alpha*x0) = 2a + dK - alpha*(x0+x1)`.
- `Det[mMat] = (a + dK - alpha*x1)(a - alpha*x0) - alpha^2*x0*x1`. Expand: `a^2 + a*dK - alpha*(a+dK)*x0 - alpha*a*x1 + alpha^2*x0*x1 - alpha^2*x0*x1 = a*(a+dK) - alpha*((a+dK)*x0 + a*x1)`, matching `t0` used in the det identity at wl:110.

So the matrix uniquely determines `lamMinus`/`lamPlus` up to ordering, and `Eigenvalues[]` is genuinely independent of the SymPy `(2A+DK-alpha*sigma ∓ R)/2` expression — the SymPy file types those eigenvalues by hand; the Mathematica file now derives them by solving `det(M - lambda I) = 0`. The post-fix Mathematica output still prints `lambda_- = (2*a + dK - alpha*(x0 + x1) - Sqrt[dK^2 + 2*alpha*dK*(x0 - x1) + alpha^2*(x0 + x1)^2])/2`, identical canonical form to before — as the directive predicted — confirming the spectrum is the same but the derivation route is different. All four downstream Part-II/III assertions (`selected overlap: HF - closed form`, `weak-loading overlap limit`, `P0_sel + beta0*d(log lambda_-)/d alpha`, `det identity`) still PASS in the Mathematica output (lines 18-19, 21-22, 30-31, 32-33).

One nuance: `Select` uses the closed-form expression as a disambiguating predicate, so the closed form is still referenced in the file. The directive explicitly anticipated this — it instructs `Select` + `First` to pick the right eigenvalue. The value bound to `lamMinus` comes from `eigVals[[k]]` (Mathematica's solver output), not from the typed expression. The disambiguator is a label, not a source. This satisfies the directive's stated criterion: "the literal expression `(2*a + dK - alpha*sigma - r)/2` must no longer appear as an *assignment* to `lamMinus`". A stronger remedy would be e.g. ordering by `Re[#]` or by free-parameter limits, but the directive did not require that and the second-engine policy is satisfied by routing the value through `Eigenvalues[]`.

## Exec log assessment

**SymPy:** exit=n/a (no `stage_030_sympy.log` exists in `redteam/exec_logs/`). I cannot verify a passing exit code from a captured log. I confirmed instead that the saved output transcript at `scripts/output/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.txt` is post-fix (mtime 17:17, vs script mtime 17:15) and contains the expected passing residuals:

- `u2 coefficient check = 0`, `u4 coefficient check = 0`, `Gamma5 coefficient check = 0` (lines 11-13).
- `selected overlap: HF - closed form = 0`, `weak-loading overlap limit = 0` (lines 30, 47).
- `P0_sel + beta0*d(log lambda_-)/d alpha = 0`, `det identity = 0` (lines 100-101).
- `normalization equivalence = 0` (line 106).
- The removed assertions (`Gamma5_sel - G5*P0_sel`, `spectral condition rewrite`) do NOT appear in the output, as required.

**Mathematica:** exit=n/a (no `stage_030_mathematica.log` exists). I confirmed the saved output at `mathematica/output/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.txt` is post-fix (mtime 17:17, vs script mtime 17:15) and contains all expected `PASS:` lines:

- `PASS: u2 coefficient check`, `PASS: u4 coefficient check`, `PASS: Gamma5 coefficient check` (lines 7, 9, 11).
- `PASS: selected overlap: HF - closed form`, `PASS: weak-loading overlap limit` (lines 19, 22).
- `PASS: P0_sel + beta0*d(log lambda_-)/d alpha`, `PASS: det identity` (lines 31, 33).
- `PASS: normalization equivalence` (line 39).
- `lambda_- = (2*a + dK - alpha*(x0 + x1) - Sqrt[...])/2` printed on line 16 — confirming `Eigenvalues[mMat]` returned the expected canonical form.
- The removed assertion lines (`PASS: Gamma5_sel - G5*P0_sel`, `PASS: spectral condition rewrite`) do NOT appear.

**Output freshness:** confirmed.
- sympy script mtime: 2026-05-21 17:15. sympy output mtime: 2026-05-21 17:17.
- Mathematica script mtime: 2026-05-21 17:15. Mathematica output mtime: 2026-05-21 17:17.
- diff patch mtime: 2026-05-21 17:17.
All outputs are newer than their corresponding scripts, so the saved transcripts reflect the post-fix state.

## Material-change assessment

`material_change`: false.

No derived quantity changes. The only edits remove two tautological assertions (with no effect on any printed value) and re-source the Mathematica eigenvalues from `Eigenvalues[mMat]` instead of typing the closed form (the value of `lamMinus`/`lamPlus` is identical, confirmed by the unchanged printed canonical form in the post-fix Mathematica output). Downstream Part-III/IV bindings (`C5_sel`, `Gamma5_sel`, `P0_sel`, `lambda_req`) print at identical canonical forms. No downstream unit should be affected.

## Side observations (non-blocking)

- In the Mathematica Part IV (line 122), `cond1 = FullSimplify[mhat^2*g5Phys*p0Sel - 2*gConst/(5*cSpeed^5), ...]`, whereas in the SymPy Part IV (line 139), `cond1 = mhat**2 * Gamma5_phys - 2 * G / (5 * c**5)` where `Gamma5_phys = G5_phys * P0_sel`. Both reduce to the same expression after `cond1 - G5_phys*cond2` simplification (both engines report `normalization equivalence = 0`), but the WL builds `cond1` slightly differently. Pre-existing and not part of any finding.
- The Mathematica file uses `radius`, `cSpeed`, `gConst` where SymPy uses `a`, `c`, `G` (line 119-120). The bindings are consistent within each engine. Pre-existing.
- The `Select` predicate in the new Part II calls `FullSimplify[# - ((2*a + dK - alpha*sigma) - r)/2, ...] === 0 &`. This is robust as long as Mathematica's `Eigenvalues` returns expressions that simplify under `$Assumptions` to the typed root. The post-fix output confirms it does. If Mathematica ever returns the eigenvalues in a form `FullSimplify` cannot reduce to that exact root, `Select` would return `{}` and `First[]` would error — the directive notes this is the correct failure mode. No issue observed in the current output.

## Verdict justification

All three findings are resolved. The two tautological assertions (F1, F2) are removed and replaced with explanatory comments verbatim per the directive; the saved output transcripts no longer carry their PASS lines. The Mathematica transliteration (F3) is broken by re-sourcing `lamMinus`/`lamPlus` from `Eigenvalues[mMat]` on an explicit matrix whose trace and determinant I verified independently match the wall-block identity. All surviving substantive checks (Part-I expansion identities, HF overlap, weak-loading limit, log-derivative identity, det identity, normalization equivalence) continue to pass in both engines, with identical canonical forms for all derived quantities. No material change downstream. Saved outputs are fresh (mtime 17:17 vs script mtime 17:15). Exec logs are not present at `redteam/exec_logs/stage_030_{sympy,mathematica}.log`, but the orchestrator's post-fix re-run is evidenced by the fresh output transcripts, and substance verification is satisfied.
