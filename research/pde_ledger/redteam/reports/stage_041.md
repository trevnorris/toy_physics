---
unit_id: 041
batch: III.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 041 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage041_rank2_support_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage041_rank2_support_mathematica_audit.txt`

## What the script claims to verify

The audit asserts five things for a 2x2 reduced loaded matrix `M(m,n,q,r)` with eigen-parameter `lambda = 1 - xi`:
(1) the determinant `det(M - (1-xi)I)` expands to the rank-2 decomposition `xi(delta+xi) - m(delta+(1+q^2)xi) - n(delta+(1+r^2)xi) + mn(q-r)^2`;
(2) the solution `n_req` of that determinant equation matches the hand-written rearrangement;
(3) `d n_req/dm` equals `-(delta+(1+qr)xi)^2 / (delta+(1+r^2)xi - m(q-r)^2)^2`, establishing monotonicity;
(4) collapsing to `r=q` yields the Stage-23 form `n_req = G_q - m`;
(5) substituting the source-tied specialization `q = t R_U, r = t, t^2 = lambda0 = 2/9` reproduces an explicit closed form `n_req^(src)` and its derivative.

## Assertion inventory

| #  | Script       | Line | Form                                            | Anchored to claim? |
|----|--------------|------|-------------------------------------------------|--------------------|
| A1 | sympy        | 72   | `simplify(Det - Det_expected) == 0`             | yes                |
| A2 | sympy        | 82   | `simplify(n_req - n_expected) == 0`             | partial (solve of linear eq) |
| A3 | sympy        | 94   | `simplify(dn_dm - monotone_expected) == 0`      | yes                |
| A4 | sympy        | 102  | `simplify(n_track - (G_q - m)) == 0`            | yes                |
| A5 | sympy        | 124  | `simplify(n_src - n_src_expected) == 0`         | **no — tautological** |
| A6 | sympy        | 143  | `simplify(dn_dm_src - dn_dm_src_expected) == 0` | yes                |
| B1 | mathematica  | 56   | `expectZero["determinant decomposition", ...]`  | yes                |
| B2 | mathematica  | 65   | `expectZero["n_req - expected", ...]`           | partial            |
| B3 | mathematica  | 76   | `expectZero["dn/dm - expected", ...]`           | yes                |
| B4 | mathematica  | 84   | `expectZero["tracking collapse", ...]`          | yes                |
| B5 | mathematica  | 112  | `expectZero["source-tied formula", ...]`        | **no — tautological** |
| B6 | mathematica  | 116  | `expectZero["source-tied dn/dm", ...]`          | yes                |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py:108-124`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl:88-112`

**What's wrong:**
Section 24.4's source-tied check writes the same algebraic expression twice and asserts their difference is zero. In sympy:

```python
q2_src = lam0 * R_U**2
r2_src = lam0
qr_src = lam0 * R_U            # defined but never used
qr_diff_sq = lam0 * (R_U - 1)**2

n_src = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + q2_src) * xi))
    / (delta + (1 + r2_src) * xi - m * qr_diff_sq)
)
n_src_expected = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + lam0 * R_U**2) * xi))
    / (delta + (1 + lam0) * xi - m * lam0 * (R_U - 1)**2)
)
expect_zero("source-tied formula", n_src - n_src_expected)
```

After expanding the `*_src` placeholders, `n_src` and `n_src_expected` are character-for-character the same rational function. The check is algebraically guaranteed by construction — it does not verify that the source-tied formula follows from the general `n_expected` derived in section 24.1 under the physical substitution `q = t R_U, r = t, t^2 = lambda0`. The author's own comment on lines 110–111 states that intent ("Because only q^2, r^2, q*r, and (q-r)^2 appear, we can substitute directly") but the substitution is never applied to `n_expected`; a new expression is constructed independently and compared to itself.

The Mathematica script repeats the same pattern at lines 88–101: `nSrc` and `nSrcExpected` are structurally identical, and `expectZero["source-tied formula", nSrc - nSrcExpected]` at line 112 is likewise vacuous.

**Why this matters:**
Theorem 5 of the unit (source-tied feasibility window) rests on `n_req^(src)` being the actual specialization of the general `n_req` formula. As written, neither engine ever performs that specialization; both just confirm a hand-written expression equals itself. If the rearranged formula contained a sign error or wrong coefficient, A5/B5 would still pass — there's no algebraic anchor to the section 24.1 result.

**Required change:**
Replace the second-definition pattern with an actual substitution of `n_expected` from section 24.1. Concretely:

- In `sympy_audit.py` around line 113, derive `n_src` by substituting the source-tied invariants into `n_expected`. The cleanest form keeps the script's own promise (lines 110–111: "only q^2, r^2, q*r, (q-r)^2 appear"):

```python
n_src = sp.simplify(
    n_expected.subs({q**2: q2_src, r**2: r2_src, q*r: qr_src, (q - r)**2: qr_diff_sq})
)
```

Then `expect_zero("source-tied formula", n_src - n_src_expected)` becomes a genuine cross-check: sympy must verify that substituting `q^2=lambda0 R_U^2, r^2=lambda0, q*r=lambda0 R_U, (q-r)^2=lambda0(R_U-1)^2` into the general formula reproduces the hand-written `n_src_expected`. (Sympy's `subs` on `(q-r)**2` will only fire if `n_expected` carries that factor literally; if it has been expanded, prefer instead the safer route of introducing a real symbol `t`, substituting `q -> t*R_U, r -> t`, expanding, then substituting `t**2 -> lam0`.)

A robust alternative that does not rely on `subs` matching the expanded form:

```python
t = sp.symbols("t", real=True)
n_src = sp.simplify(
    n_expected.subs({q: t * R_U, r: t}).rewrite(sp.Pow)
)
n_src = sp.simplify(n_src.subs(t**2, lam0))
expect_zero("source-tied formula", n_src - n_src_expected)
```

Apply the same fix in the Mathematica script around line 93. Replace the standalone `nSrc` construction with a substitution of `nExpected`:

```
nSrc = FullSimplify[
  nExpected /. {q -> t rU, r -> t},
  Assumptions -> $Assumptions
];
nSrc = FullSimplify[nSrc /. t^2 -> lambda0, Assumptions -> $Assumptions];
expectZero["source-tied formula", nSrc - nSrcExpected];
```

(or use the direct-substitution form `nExpected /. {q^2 -> lambda0 rU^2, r^2 -> lambda0, q r -> lambda0 rU, (q - r)^2 -> lambda0 (rU - 1)^2}` if `nExpected` is left in factored form). Add `t` to the `Clear[...]` list at line 35 and to `$Assumptions` at lines 37–38.

The `dn_dm_src` derivative test in section 24.4 should be left intact (it differentiates `n_src_expected`, so its claim — that the derivative simplifies to the stated form — is non-tautological).

**Verification:**
After the fix, the sympy script will compute `n_src` by substituting into `n_expected` (defined at line 78). The captured output at line 58 of `moving_throat_pde_stage041_rank2_support_sympy_audit.txt` (currently `n_req^(src) = ...`) must still print and the line `source-tied formula = 0` must still appear. The script must still `Exit 0`. Verifier should also inspect the diff to confirm `n_src` is now constructed from `n_expected` rather than as a fresh literal. Analogous changes to the Mathematica output file (line 41).

## Independent-derivation check (Mathematica)

The `.wl` script shares the same overall structure as the `.py` script (same matrix `mMat`, same `Det[...]` vs hand-written `detExpected`, same `Solve[...]` for `n`, same `D[nExpected, m]`, and the same standalone construction of `nSrc`). However, each engine's symbolic algebra (`FullSimplify` vs sympy's `simplify`) does independent expansion/cancellation work on the determinant, on the solve, on the derivative quotient, and on the `r=q` collapse. The hand-written `*_expected` forms are necessary because they are the *claims under test*. I do not classify this as a `mathematica_transliteration` finding: both engines independently expand and cancel the same algebraic content from the same matrix premise.

Quoting comparable sections:

- SymPy `Det = sp.expand((M - lam * sp.eye(2)).det())` (line 62) vs Mathematica `detExpr = Expand[Det[mMat - lam IdentityMatrix[2]]]` (line 50). Both engines independently compute the 2x2 determinant.
- SymPy `dn_dm = sp.simplify(sp.diff(n_expected, m))` (line 86) vs Mathematica `dndm = FullSimplify[D[nExpected, m], ...]` (line 69). Both independently apply the quotient rule and simplify.

The shared design flaw in section 24.4 (`n_src` vs `n_src_expected` being structurally identical) is captured under F1 above; that is a methodology bug in both scripts, not a transliteration of one into the other.

## Engine cross-check

Both engines produce the same final symbolic results, and both report `PASS`/`exit 0`:

| Quantity                  | SymPy output (line)              | Mathematica output (line)        |
|---------------------------|----------------------------------|----------------------------------|
| `det decomposition`       | `= 0` (txt line 23)              | `= 0`, PASS (txt line 18–19)    |
| `n_req - expected`        | `= 0` (txt line 30)              | `= 0`, PASS (txt line 21–22)    |
| `dn/dm - expected`        | `= 0` (txt line 42)              | `= 0`, PASS (txt line 28–29)    |
| `tracking collapse`       | `= 0` (txt line 53)              | `= 0`, PASS (txt line 35–36)    |
| `source-tied formula`     | `= 0` (txt line 64)              | `= 0`, PASS (txt line 42–43)    |
| `source-tied dn/dm`       | `= 0` (txt line 82)              | `= 0`, PASS (txt line 47–48)    |

Engines agree. The agreement is real for A1/A3/A4/A6/B1/B3/B4/B6 (genuine symbolic computations) and trivial for A5/B5 (identical expressions on both sides — F1).

Side-by-side numerator/denominator forms for `n_req`:

- SymPy (txt 24–29): `n_req = (delta*m - delta*xi + m*q^2*xi + m*xi - xi^2) / (-delta + m*q^2 - 2*m*q*r + m*r^2 - r^2*xi - xi)` (overall sign on both numerator and denominator differs from the hand-written form, but the ratio is identical).
- Mathematica (txt 20): `n_req = (delta*(-m + xi) + xi*(-(m*(1 + q^2)) + xi)) / (delta - m*(q - r)^2 + xi + r^2*xi)`.

Multiplying SymPy's by `-1/-1` yields Mathematica's form. Agreement confirmed.

For `dn/dm`: SymPy (txt 35–41) and Mathematica (txt 27) both report `-(delta + (1+qr)xi)^2 / (delta - m(q-r)^2 + xi + r^2 xi)^2`. Agree.

For source-tied derivative: SymPy `-(2 R_U xi + 9 delta + 9 xi)^2 / (2 R_U^2 m - 4 R_U m - 9 delta + 2 m - 11 xi)^2` (txt 75–81) and Mathematica `-((9 delta + (9+2 rU) xi)^2 / (9 delta - 2 m (-1 + rU)^2 + 11 xi)^2)` (txt 46). Both are `-(9 delta + (9 + 2 rU) xi)^2 / (9 delta + 11 xi - 2 m (rU - 1)^2)^2` after factoring/expanding the denominator: `9 delta + 11 xi - 2 m (R_U^2 - 2 R_U + 1) = 9 delta + 11 xi - 2 m R_U^2 + 4 m R_U - 2 m`, which is `-(2 R_U^2 m - 4 R_U m - 9 delta + 2 m - 11 xi)`. Agree.

## Verdict justification

The script verifies four substantive theorems (rank-2 determinant decomposition, the support-loading formula's rearrangement, the monotonicity derivative, and the tracking-collapse identity), and both engines agree on those results. One section — 24.4, the source-tied specialization — is a tautological self-comparison in both engines; the script never algebraically connects `n_req^(src)` to the general `n_req` derived in section 24.1. That is a real finding: the source-tied formula is announced as a *consequence* of the general theorem under the substitution `q = t R_U, r = t, t^2 = lambda0`, but neither engine performs that substitution. Verdict: `findings`, count 1.

Attacks attempted that failed:
- Hand-expanded `det(M - (1-xi)I)` against `Det_expected` — matches.
- Hand-derived quotient-rule numerator for `dn/dm` and matched against `-(delta + (1+qr)xi)^2`. Verified `-(delta+(1+q^2)xi)(delta+(1+r^2)xi) + xi(delta+xi)(q-r)^2 = -delta^2 - 2 delta xi (1+qr) - xi^2 (1+qr)^2`.
- Substituted `r=q` into `n_expected` by hand and recovered `xi(delta+xi)/(delta+(1+q^2)xi) - m = G_q - m`. Matches.
- Substituted `q -> t R_U, r -> t, t^2 -> lambda0` into the general `n_expected`. Recovered exactly `n_src_expected`. So Theorem 5 is *mathematically true*; the finding is solely that the script does not verify it — only restates it.
- Variable-domain check: `delta > 0, xi > 0, R_U > 0` are used only as labels; no `simplify` step requires them. `m, n, q, r` are real with no sign assumption, consistent with the formulas.
- Stale-output check: sympy output 2026-05-11 12:42 (script 2026-04-01 12:39); Mathematica output 2026-05-11 12:49 (script 2026-05-11 11:56). Both outputs newer than scripts — fresh.

## Self-test notes

- **Variable independence:** verified `n_expected` depends on `m, n` (implicit via solve), `q, r, xi, delta`; `sp.diff(n_expected, m)` is non-degenerate (numerator after quotient rule is `-(delta+(1+qr)xi)^2`, a non-zero symbolic expression). Similarly, the proposed fix uses `n_expected.subs(...)` which preserves dependence on `m, xi, delta, R_U` after substitution; `n_src_expected` depends on the same variables.
- **Parity/symmetry:** not applicable — no integrals over symmetric domains in this audit unit. The algebraic identities are checked structurally.
- **Trivial-case substitution:** at `R_U = 1`, `n_src_expected` numerator becomes `xi(delta+xi) - m(delta + (1 + 2/9)xi) = xi(delta+xi) - m(delta + 11 xi/9)`, denominator `delta + (1+2/9)xi - 0 = delta + 11 xi/9`, giving `n_src = xi(delta+xi)/(delta+11xi/9) - m`. The proposed substitution-based construction must reproduce this same value — and at `R_U=1` we have `q = r = t`, so `n_expected` collapses to its `r=q` form `G_q - m` with `q^2 = lambda0`, i.e. `xi(delta+xi)/(delta + (1+lambda0)xi) - m = xi(delta+xi)/(delta+11xi/9) - m`. Matches. So the proposed fix would yield 0 in `expect_zero("source-tied formula", ...)`.
- **Path specifications:** no `missing_verification_script` findings; no new file paths needed in the directive. Edits target existing `scripts/.../stage041_..._sympy_audit.py` and `mathematica/.../stage041_..._mathematica_audit.wl`.
