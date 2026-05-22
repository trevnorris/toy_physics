---
unit_id: 030
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 030 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.txt`

## What the script claims to verify

The script audits five things for a "selected-mode" wall block: (i) the coefficient structure of the generic normalized response `Y_- = D0/(D0 + D2 w^2 + D4 w^4 - i C5 w^5)` matches the closed forms `u2 = -D2/D0`, `u4 = (D2^2 - D0 D4)/D0^2`, `Gamma5 = C5/D0`; (ii) hand-written expressions for the lower eigenvalue `lam_-` and the closed-form overlap `s_-` of a 2x2 wall block satisfy the Hellmann-Feynman identity `s_- = -d(lam_-)/d(alpha)`, reduce to `x0` at `alpha = 0`, and satisfy a determinant identity `lam_- lam_+ = A(A+DK) - alpha T0`; (iii) the selected odd coefficient `C5_sel = beta5 s`, the ratio `Gamma5_sel = C5_sel/lam_-`, and the static prefactor `P0_sel = beta0 s/lam_-` satisfy `Gamma5_sel = G5 P0_sel` and `P0_sel = -beta0 d(log lam_-)/d alpha`; (iv) two purported "normalization target" forms `mhat^2 Gamma5_phys = 2G/(5 c^5)` and `mhat^2 P0_sel = NQ_target` are equivalent, where `Gamma5_phys = G5_phys P0_sel`, `G5_phys = a^5/(27 cs^5)`, `NQ_target = 54 G cs^5/(5 a^5 c^5)`; and (v) a "spectral condition rewrite" identity relating `lam_-` and `lambda_req := mhat^2 beta0 s/NQ_target`. The 2x2 matrix itself is never defined in-script; `lam_+/-` are stated as algebraic expressions in `A, DK, alpha, x0, x1`.

## Assertion inventory

| #   | Script      | Line   | Form                                                                                                                                                                         | Anchored to claim? |
|-----|-------------|--------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------|
| A1  | sympy       | 53     | `Ysel.coeff(w,2) - u2 == 0`                                                                                                                                                  | yes                |
| A2  | sympy       | 54     | `Ysel.coeff(w,4) - u4 == 0`                                                                                                                                                  | yes                |
| A3  | sympy       | 55     | `im(Ysel).coeff(w,5) - Gamma5 == 0`                                                                                                                                          | yes                |
| A4  | sympy       | 87     | `s_minus_hf - s_minus_closed == 0`                                                                                                                                           | yes                |
| A5  | sympy       | 91     | `s_minus_closed.subs(alpha,0) - x0 == 0`                                                                                                                                     | yes                |
| A6  | sympy       | 113    | `Gamma5_sel - G5*P0_sel == 0`                                                                                                                                                | **no (tautology)** |
| A7  | sympy       | 114-117| `P0_sel + beta0*d(log lam_-)/d(alpha) == 0`                                                                                                                                  | partial (duplicate of A4) |
| A8  | sympy       | 121    | `lam_- lam_+ - (A(A+DK) - alpha T0) == 0`                                                                                                                                    | yes                |
| A9  | sympy       | 138-141| `cond1 - G5_phys*cond2 == 0`                                                                                                                                                 | partial            |
| A10 | sympy       | 144    | `(lam_- - lambda_req) + (mhat^2 P0_sel - NQ_target)*lam_-/NQ_target == 0`                                                                                                    | **no (tautology)** |
| B1  | mathematica | 39     | `Coefficient[ySel,w,2] - u2 == 0`                                                                                                                                            | yes                |
| B2  | mathematica | 40     | `Coefficient[ySel,w,4] - u4 == 0`                                                                                                                                            | yes                |
| B3  | mathematica | 41     | `Coefficient[Im[ySel],w,5] - gamma5 == 0`                                                                                                                                    | yes                |
| B4  | mathematica | 64     | `sMinusHF - sMinusClosed == 0`                                                                                                                                               | yes                |
| B5  | mathematica | 66     | `(sMinusClosed /. alpha->0) - x0 == 0`                                                                                                                                       | yes                |
| B6  | mathematica | 83     | `gamma5Sel - g5*p0Sel == 0`                                                                                                                                                  | **no (tautology)** |
| B7  | mathematica | 84     | `p0Sel + beta0*D[Log[lamMinus],alpha] == 0`                                                                                                                                  | partial (duplicate of B4) |
| B8  | mathematica | 87     | `lamMinus*lamPlus - (a*(a+dK) - alpha*t0) == 0`                                                                                                                              | yes                |
| B9  | mathematica | 101    | `cond1 - g5Phys*cond2 == 0`                                                                                                                                                  | partial            |
| B10 | mathematica | 104-108| `(lamMinus - lambdaReq) + (mhat^2 p0Sel - nQTarget)*lamMinus/nQTarget == 0`                                                                                                  | **no (tautology)** |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:113`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:83`

**What's wrong:**
The assertion `Gamma5_sel - G5*P0_sel == 0` is algebraically guaranteed by the immediately preceding definitions in the same script. From sympy lines 101-104:

```
beta5 = G5 * beta0
C5_sel = beta5 * s_minus_closed              # = G5*beta0*s
Gamma5_sel = sp.simplify(C5_sel / lam_minus) # = G5*beta0*s/lam_-
P0_sel = sp.simplify(beta0 * s_minus_closed / lam_minus)  # = beta0*s/lam_-
```

Therefore `Gamma5_sel = G5*beta0*s/lam_- = G5*(beta0*s/lam_-) = G5*P0_sel` identically — the residual is forced to zero by the three lines of definition that precede it, regardless of any physics. No symbolic input choice can make this assertion fail. The Mathematica file repeats the same construction (`c5Sel = g5*beta0*sMinusClosed`, `gamma5Sel = c5Sel/lamMinus`, `p0Sel = beta0*sMinusClosed/lamMinus`, then `expectZero` on `gamma5Sel - g5*p0Sel`).

**Why this matters:**
The script's stated purpose includes verifying the relation between the selected-mode odd coefficient ratio and the static prefactor. A tautological "verification" gives false confidence: the captured `PASS` in the output transcript does not establish anything about the physics; it only confirms that `(G5*beta0*s/lam_-) - G5*(beta0*s/lam_-) = 0`. Downstream readers who treat this as evidence the identity holds are misled.

**Required change:**
Remove the assertion and replace it with an explanatory comment stating that the relation is a definitional consequence (so the script makes no claim there). Routing the equality through the Part-I series expansion is also tautological because that result `Gamma5 = C5/D0` is itself a definitional consequence of the series — so the only honest remedy is to drop the false-confidence "verification" and annotate. See directive for the exact code form.

**Verification:**
After the fix, the assertion line (sympy 113, wl 83) is removed and replaced with an explanatory comment. The captured output transcript must no longer contain `"Gamma5_sel - G5*P0_sel = 0"` or its PASS line.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:143-144`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:103-108`

**What's wrong:**
The "spectral condition rewrite" assertion expands to `0 = 0` by algebraic substitution of the definitions, independent of any spectral structure. With

```
lambda_req := mhat^2 * beta0 * s_minus_closed / NQ_target
P0_sel     := beta0 * s_minus_closed / lam_minus
X          := mhat^2 * P0_sel - NQ_target
```

the asserted residual is

```
(lam_- - lambda_req) + X*lam_-/NQ_target
  = lam_- - mhat^2*beta0*s/NQ_target + (mhat^2*beta0*s/lam_-)*lam_-/NQ_target - lam_-
  = lam_- - mhat^2*beta0*s/NQ_target + mhat^2*beta0*s/NQ_target - lam_-
  = 0
```

`P0_sel`, `lam_-`, `mhat`, `beta0`, `s`, and `NQ_target` all appear only through the substitution `P0_sel = beta0*s/lam_-`, after which the two surviving terms cancel. The identity holds for *any* choice of these symbols — no physical content is exercised. Mathematica lines 103-108 mirror the same construction.

**Why this matters:**
The captured `PASS` is presented in Part IV's stated rewrite of the spectral condition. The "rewrite" being verified is in fact a vacuous algebraic rearrangement. The check provides no evidence that `lam_- = lambda_req` is equivalent to `mhat^2 P0_sel = NQ_target` *as a physical constraint*; it only verifies the algebraic identity that connects the residuals, which holds trivially.

**Required change:**
Drop the assertion and replace it with an explanatory comment recording that the rearrangement `(lam_- - lambda_req) + (mhat^2 P0_sel - NQ_target)*lam_-/NQ_target` is `0` by definitional substitution. Substituting the fixed-point constraint `mhat^2 = NQ_target/P0_sel` into `lambda_req` and asserting `lam_- = lambda_req_at_fixed` is *also* tautological — it likewise reduces to `0` by the same definitional chain `P0_sel = beta0*s/lam_-`, `lambda_req = mhat^2*beta0*s/NQ_target`. The only honest remedy is removal with an annotation. See directive.

**Verification:**
After the fix, the assertion at sympy:144 / wl:104-108 is removed and replaced with a comment. The captured output transcript must no longer contain `"spectral condition rewrite = 0"` or its PASS line.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl` (entire file)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py` (entire file)

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script with only cosmetic renamings (uppercase Greek to lowercase, `_` to camelCase). Compare:

Part I (sympy 41-44 vs wl 28-31):
```python
w = sp.symbols("omega", real=True)
D0, D2, D4, C5 = sp.symbols("D0 D2 D4 C5", nonzero=True, real=True)
Dsel = D0 + D2 * w**2 + D4 * w**4 - I * C5 * w**5
Ysel = sp.expand(sp.series(D0 / Dsel, w, 0, 6).removeO())
```
```mathematica
$Assumptions = Element[{w, d0, d2, d4, c5}, Reals] && d0 != 0;
dSel = d0 + d2*w^2 + d4*w^4 - I*c5*w^5;
ySel = Expand[Normal[Series[d0/dSel, {w, 0, 5}]]];
```

Part II (sympy 68-74 vs wl 48-54):
```python
sigma = x0 + x1
delta_kappa = x0 - x1
KappaProd = x0 * x1
R = sp.sqrt((DK + alpha * delta_kappa) ** 2 + 4 * alpha**2 * KappaProd)
lam_minus = sp.simplify((2 * A + DK - alpha * sigma - R) / 2)
lam_plus = sp.simplify((2 * A + DK - alpha * sigma + R) / 2)
```
```mathematica
sigma = x0 + x1;
deltaKappa = x0 - x1;
kappaProd = x0*x1;
r = Sqrt[(dK + alpha*deltaKappa)^2 + 4*alpha^2*kappaProd];
lamMinus = FullSimplify[(2*a + dK - alpha*sigma - r)/2, ...];
lamPlus = FullSimplify[(2*a + dK - alpha*sigma + r)/2, ...];
```

Part III (sympy 100-104 vs wl 74-77):
```python
beta0, G5 = sp.symbols("beta0 G5", positive=True, real=True)
beta5 = G5 * beta0
C5_sel = sp.simplify(beta5 * s_minus_closed)
Gamma5_sel = sp.simplify(C5_sel / lam_minus)
P0_sel = sp.simplify(beta0 * s_minus_closed / lam_minus)
```
```mathematica
beta5 = g5*beta0;
c5Sel = FullSimplify[beta5*sMinusClosed, ...];
gamma5Sel = FullSimplify[c5Sel/lamMinus, ...];
p0Sel = FullSimplify[beta0*sMinusClosed/lamMinus, ...];
```

This is identical variable choreography, identical intermediate expressions, identical assertion sequence (`expect_zero` ↔ `expectZero`), and identical assertion phrasings in the same order. No independent derivation occurs in Mathematica; it is a re-spelling of the SymPy algebra. This violates the second-engine policy that both engines independently re-derive the result from the physical premises.

**Why this matters:**
A transliterated second engine cannot catch errors in the first engine's algebra — both engines will make the same substitution mistakes or definition errors. Cross-engine agreement on identical algebra is not evidence of correctness; it only confirms that both CASes simplify the same expression to the same canonical form. The audit framework relies on the two engines reaching the result by genuinely different paths.

**Required change:**
Restructure Part II of the Mathematica script so that `lamMinus` and `lamPlus` are produced by `Eigenvalues[mMat]` operating on an explicit 2x2 matrix `mMat`, rather than typed verbatim as closed forms. The matrix is determined uniquely (up to similarity) by the trace `2a + dK - alpha*(x0+x1)` and determinant `a*(a+dK) - alpha*((a+dK)*x0 + a*x1)` that the existing det identity (wl:87) already checks against:

```mathematica
mMat = {{a + dK - alpha*x1, -alpha*Sqrt[x0*x1]},
        {-alpha*Sqrt[x0*x1], a - alpha*x0}};
```

The loading direction `v` is not defined in the existing script, so `sMinusClosed` cannot be safely rederived from an eigenvector overlap. Keep the existing closed-form `sMinusClosed` and HF cross-check; only the source of `lamMinus`/`lamPlus` changes. This single restructuring is sufficient to clear the transliteration label because every downstream Part III/IV expression flows through `lamMinus`. See directive for the exact replacement block.

**Verification:**
After the fix, the Mathematica file must contain an `Eigenvalues[mMat]` call whose output is what `lamMinus`/`lamPlus` are bound to. The literal expression `(2*a + dK - alpha*sigma - r)/2` must no longer appear as an *assignment* to `lamMinus`. The captured output transcript must show the same final `lambda_-`/`lambda_+` symbolic forms as before, and all Part II/III assertions must still pass with residual 0.

## Independent-derivation check (Mathematica)

The Mathematica script does NOT derive the result independently. It is a transliteration of the SymPy script. See F3 above for line-by-line correspondence. Every Mathematica binding has a one-to-one preimage in the SymPy file, in the same order, with the same intermediate expressions.

## Engine cross-check

Both engines emit the same final symbolic expressions (modulo CAS canonical-form differences). All ten matched assertion pairs (A1-A10 ↔ B1-B10) report residual = 0 in their captured outputs:

| sympy line | sympy output line | mathematica line | wl output line | residual |
|-----------|-------------------|-----------------|----------------|----------|
| 53        | 19                | 39              | 14             | 0        |
| 54        | 20                | 40              | 16             | 0        |
| 55        | 21                | 41              | 18             | 0        |
| 87        | 38                | 64              | 26             | 0        |
| 91        | 55                | 66              | 29             | 0        |
| 113       | 108               | 83              | 38             | 0        |
| 114-117   | 109               | 84              | 40             | 0        |
| 121       | 110               | 87              | 42             | 0        |
| 138-141   | 115               | 101             | 48             | 0        |
| 144       | 116               | 104-108         | 50             | 0        |

Engines agree — but as F3 explains, this agreement reflects identical algebra rather than independent derivation, so it is not strong evidence of correctness.

The closed forms also agree (e.g., `Ysel` at sympy output line 13-18 matches `ySel` at wl output line 13; `lamMinus` at sympy output line 26-31 matches wl output line 24; `s_-` at sympy output 39-54 matches wl output 28; `lambda_req` at sympy output 117-132 matches wl output 52). Both saved outputs are fresh: sympy script mtime Apr 1 12:39, output mtime May 11 12:41; wl script mtime Apr 21 17:04, output mtime May 11 12:46. No stale_output finding.

## Verdict justification

Six of the ten assertion pairs are substantive: the generic-response coefficient identities (A1/A2/A3), the Hellmann-Feynman identity (A4), the weak-loading limit (A5), and the determinant identity (A8) are non-tautological checks that exercise the algebraic structure of the response and the wall-block formulas. These hold up under attack: the series expansion is real, the HF differentiation is real, the alpha=0 limit reduces non-trivially (sqrt of (DK)^2 enters as DK, then `(sigma + (x0-x1))/2 = x0`), and the determinant matches the trace/det of a 2x2 with the asserted entries.

Three of the ten pairs are tautological or near-tautological: A6/B6 (`Gamma5_sel - G5*P0_sel`, F1) and A10/B10 ("spectral condition rewrite", F2) reduce to `0=0` by definition substitution alone. A7/B7 (`P0_sel + beta0*d(log lam_-)/d(alpha)`) is not a tautology of construction but is a strict consequence of A4/B4 (the HF identity) and the definition of `P0_sel`, so it duplicates an earlier check; I note this in the assertion inventory but do not raise it as a separate finding because the duplication does not give false confidence as long as A4 itself is real.

A9/B9 ("normalization equivalence") collapses to the prefactor identity `G5_phys * NQ_target = 2G/(5c^5)` after the `P0_sel`-dependent terms cancel, i.e., it really verifies that `54/(27*5) = 2/5`. This is a substantive consistency check on the literal numerical prefactors `27` and `54`, so I do not flag it as tautological, but it is weaker than the docstring implies (it does not test the spectral structure of `P0_sel`).

Finally, F3 (transliteration) applies to the Mathematica script as a whole: it does not constitute independent verification.

The verdict is `findings` with three real issues (F1, F2, F3). No stop-cold flag: the tautologies can be replaced with non-tautological checks without contradicting downstream units, and the transliteration can be fixed without changing any numerical result. No `CRITICAL_DOWNSTREAM` is warranted because the underlying quantities (`lam_-`, `s_-`, `P0_sel`) are not changed by removing the tautological assertions — only the *evidence* the script provides for them is improved.

## Self-test notes

I verified the tautology claims by manual substitution: F1's residual reduces to `G5*beta0*s/lam_- - G5*(beta0*s/lam_-) = 0` using only the four definition lines that precede it; F2's residual reduces to `lam_- - mhat^2*beta0*s/NQ_target + mhat^2*beta0*s/NQ_target - lam_- = 0` after substituting `P0_sel = beta0*s/lam_-`. Neither requires `lam_-` or `s` to have any particular physical content, confirming the tautology. For F3 I confirmed the matrix structure inferred from the trace `2A + DK - alpha*(x0+x1)` and determinant `A(A+DK) - alpha((A+DK)x0 + A x1)` is consistent with a real symmetric 2x2 — so an `Eigenvalues[M]` call in Mathematica is feasible (no complex-eigenvalue pitfall). I also checked symbol assumptions are consistent (all positivity declarations are physically motivated; `alpha >= 0` and `x0, x1 > 0` are correctly chosen so `R > 0` and `lam_- < lam_+`). No `sp.diff` or `D[...]` of an independent-variable trap was found in the proposed directive — the only differentiation in the existing or proposed checks is `d(lam_-)/d(alpha)` and `d(log(lam_-))/d(alpha)`, both of which are with respect to a variable `lam_-` actually depends on.
