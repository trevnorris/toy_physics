---
unit_id: 013
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 013 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

Per the docstring, this script audits the "Master-note" notes file `step_11_projected_maxwell_mouth_taylor_master_notes.md`. The assertions split into three groups: (1) a one-sided Taylor projection check that uses the moments of `exp(-u)` and `u*exp(-u)` on `[0, oo)` and confirms the first-moment recovery to order `ell^1`; (2) a hardcoded family of chain-rule expansions `z0, z2, z4, n0, n2, n4` for a bottleneck quantity and its `P`-shift, with derived combinations `Xi, K1, H_even, deltaP2, deltaP4`, which the script then probes for `S2/Hport/Gw/P`-prime dependence and tests a single substantive coefficient (`dXi/dPprime = 2/P` and `d(deltaP2)/dGprime = -2*P/(D0*Delta^2)`); and (3) a mechanism "sieve": a 2x2 linear-solve argument that `K1 = H_even = 0` admits only the trivial source/denominator pair `(Qx, Dx) = 0` and only the trivial spectral pair `(Sx, Hx) = 0`. The output transcript is a single "STATUS: PASS" line.

## Assertion inventory

| #  | Script | Line  | Form                                                                                          | Anchored to claim? |
|----|--------|-------|-----------------------------------------------------------------------------------------------|--------------------|
| A1 | sympy  | 26    | `series(Xproj, ell, 0, 2).removeO() - (X0 + ell*X1) == 0`                                     | yes                |
| A2 | sympy  | 27-30 | mutant of A1 must be nonzero                                                                  | yes                |
| A3 | sympy  | 34    | `mu1_W2 - 2 == 0` (closed-form moment of `u^2*exp(-u)`)                                       | yes                |
| A4 | sympy  | 35    | `series(Xproj_W2, ell, 0, 2).removeO() - (X0 + ell*mu1_W2*X1) == 0`                           | yes                |
| A5 | sympy  | 36-39 | mutant of A4 must be nonzero                                                                  | yes                |
| A6 | sympy  | 86-87 | `diff(Xi, Sx) == 0`, `diff(Xi, Hx) == 0`, `diff(Xi, Gx) == 0`                                 | **no (tautological)** |
| A7 | sympy  | 88-90 | `diff(K1, Px) == 0`, `diff(K1, Gx) == 0`, `diff(H_even, Px) == 0`, `diff(H_even, Gx) == 0`    | **no (tautological)** |
| A8 | sympy  | 91    | `diff(Xi, Px) - 2/P == 0`                                                                     | yes                |
| A9 | sympy  | 109   | `diff(deltaP2_der, Gx) + 2*P/(D0*Delta**2) == 0`                                              | yes                |
| A10| sympy  | 110-111| `diff(deltaP4_der, Gx) != 0`                                                                  | yes                |
| A11| sympy  | 131   | `diff(Xi, Px) != 0` (duplicate of A8 minus the specific coefficient)                          | partial            |
| A12| sympy  | 132   | `diff(deltaP4_der, Gx) != 0` (duplicate of A10)                                               | partial            |
| A13| sympy  | 133-134| `qd_matrix.det() != 0`, `sh_matrix.det() != 0`                                               | yes                |
| A14| sympy  | 135-138| linear-solve `qd_only == [{Qx:0, Dx:0}]`, `sh_only == [{Sx:0, Hx:0}]`                        | partial (redundant with A13) |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl` (does not exist)

**What's wrong:**
Unit 013 has `is_status_only_candidate: False` and is not a checkpoint, so the policy requires both engines. No `.wl` file exists for this stage. Adjacent stages 010, 011, 012 each carry a paired `.wl` audit; stage 013 breaks the pattern. There is therefore no second-engine cross-check on any of A1-A14.

**Why this matters:**
Without an independent Mathematica derivation, the one-sided Taylor first-moment relation `int_0^inf u*exp(-u) du = 1`, the closed-form second moment `int_0^inf u^2*exp(-u) du = 2`, the hardcoded chain-rule structures of `z0, z2, z4, n0, n2, n4`, the substantive coefficient checks `dXi/dPprime = 2/P` and `d(deltaP2)/dGprime = -2*P/(D0*Delta^2)`, and the 2x2 sieve determinants are confirmed by only one CAS. The second-engine policy exists precisely because a single CAS's `simplify` path can silently agree with the formula it was handed.

**Required change:**
Add a new file `mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl` that independently rederives and asserts each of M1-M6 listed in the directive. The script must derive the chain-rule expansions from a master `Q/(Delta - S2 - ...)` style primitive (i.e., write the source expression and let Mathematica do the Taylor expansion), not transliterate the literal `z0, z2, z4` polynomials from the `.py` file.

**Verification:**
Verifier runs `redteam exec-mathematica 013`; the new script must exist, exit 0, and its saved transcript must show `Print` lines confirming each of M1-M6.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:86-90`

**What's wrong:**
The script declares Xi, K1, H_even by writing literal symbolic formulas:

  Xi = `(2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)).subs(subs_der) / mu1`
  K1 = `(-(z2 + z0/9)).subs(subs_der) / mu1`
  H_even = `((-z4 + Rational(2,3)*z2 - z0/27)).subs(subs_der) / mu1`

By inspection of the source code, none of these formulas contains the symbols `s1, h1, g1` (Xi), nor `p1, g1` (K1, H_even). After `subs_der`, that translates to: Xi contains no `{Sx, Hx, Gx}`, and K1/H_even contain no `{Px, Gx}`. Consequently the assertions

  `assert_zero(f"Xi independence from {sym}", sp.diff(Xi, sym))` for `sym in (Sx, Hx, Gx)` — lines 86-87
  `assert_zero(f"even gates independence from {sym}", sp.diff(K1, sym))` for `sym in (Px, Gx)` — line 89
  `assert_zero(f"even gates independence from {sym}", sp.diff(H_even, sym))` for `sym in (Px, Gx)` — line 90

each compute `sp.diff(EXPR, VAR)` where `VAR` literally does not appear in `EXPR`. The derivative is identically zero by SymPy's variable-presence rule, and `assert_zero` cannot fail. These checks confirm only that Python's substring `Sx`/`Hx`/`Gx`/`Px` is absent from the formula the developer wrote two lines earlier; they exercise zero physics.

**Why this matters:**
The script's "STATUS: PASS" transcript is partly purchased by checks that cannot fail. A reader scanning the assertion count will mistake five tautologies for five physical confirmations. If a future edit accidentally introduces a `Gw`/`Sx`/etc. contamination into the master formulas, these particular assertions will catch it only because the symbol token would then *appear* in the formula; that is a much weaker guarantee than the one the comment "Xi independence from Sx" implies.

**Required change:**
Delete the five tautological independence assertions at lines 86-90. They check `sp.diff(EXPR, VAR)` where `VAR` is not a symbol of `EXPR` and so are guaranteed to pass regardless of the physics. The surviving substantive coefficient checks at line 91 (`dXi/dPprime - 2/P == 0`) and line 109 (`d(deltaP2_der)/dGprime + 2*P/(D0*Delta**2) == 0`), together with the determinant and linear-solve assertions at lines 133-138, continue to exercise the K1/H_even/Xi forms non-trivially via routes that *can* fail. Replacing the five lines with additional substantive coefficient assertions is a scope extension and is therefore out of scope for this correction; deletion is the minimal fix.

**Verification:**
Verifier checks (a) lines 86-90 (the two `for sym in ...:` loops on independence) are gone, (b) the substantive assertion at line 91 (`dXi/dPprime`) is preserved, (c) the script still exits 0 and the determinant / solve assertions at lines 133-138 still execute.

### F3 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:49-80`

**What's wrong:**
The chain-rule expansion expressions

  `z0 = (Delta*q1 - Q*d1)/Delta**2` (line 49)
  `z2 = (...)/Delta**3` (line 50)
  `z4 = (...)/Delta**4` (lines 51-60)
  `n0, n2, n4` (lines 62-80)

are written as literal symbolic polynomials. There is no in-script derivation of where they come from (no Taylor-expansion of a primitive like `Q/(Delta - S2 - Hport*ell^2 - ...)`), no comment naming the source notes section beyond the docstring's generic mention of `step_11_projected_maxwell_mouth_taylor_master_notes.md`, and no comment naming the upstream stage script that verified them. The downstream "physics" assertions A6, A7, A8, A9, A10 then test properties of these hardcoded expressions; if a sign or factor in (say) `z4`'s `-3*Q*S2**2*d1` term is wrong, this script will happily certify the wrong answer because it never derives the term independently.

**Why this matters:**
The script is positioned as a "master-note audit" — the bottom line is that someone reading the transcript should believe the master-note's chain-rule formulas have been independently verified. But the chain-rule formulas are exactly what the script accepts on faith. The downstream sieve check, the determinant check, and the coefficient checks all sit on top of these hardcoded forms; they cannot detect an error inside `z0, z2, z4, n0, n2, n4`. A "PASS" here misrepresents the audit's coverage.

**Required change:**
Add an in-script derivation for `z0, z2, z4, n0, n2, n4` *or* add a precise upstream anchor. Concretely:

Option A (derivation): Above line 49, introduce a primitive `Zsource = Q / (Delta - S2*ell**2 - Hport*ell**4 + sp.O(ell**6))`-style symbolic source (or whatever the master-note actually defines; absent an explicit master form the auditor cannot specify it). Use `sp.series(Zsource.subs(...), ell, 0, 5)` plus an explicit `t` substitution to extract `z0, z2, z4` and assert each one matches the literal value the script currently uses. Likewise for `n0, n2, n4` via the `P/(Delta - S2*ell^2 - Gw*ell^4 + ...)` analog.

Option B (anchor): Insert a comment block immediately above line 49 of the form

    # z0, z2, z4 are the order-ell^0, ell^2, ell^4 coefficients of the master
    # bottleneck Q/(Delta - S2*ell^2 - Hport*ell^4 + ...) chain rule, derived
    # and verified in stage 0XX (file: scripts/moving_throat_pde_stage0XX_...py
    # assertion `...`). Likewise n0, n2, n4 come from the P-analog in stage 0YY.

filling in the actual upstream stage numbers and assertion labels. The auditor cannot confirm those numbers from the script alone, so the directive asks Codex to verify the upstream anchor exists; if it does not, Codex is blocked.

**Verification:**
Verifier confirms either (a) new derivation code appears before line 49 and asserts each `z*` and `n*` against its literal form, or (b) a citation comment block names a concrete upstream stage + assertion label, and that upstream script/assertion actually exists.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` file exists. See F1.

## Engine cross-check

Not applicable — only one engine present.

## Verdict justification

The one-sided Taylor projection block (A1-A5) is real algebra and holds up: the moment integrals `int_0^inf exp(-u) du = 1`, `int_0^inf u*exp(-u) du = 1`, `int_0^inf u^2*exp(-u) du = 2` are non-trivial closed forms, and the assertions confirm specific coefficient relations rather than tautologies. The single substantive coefficient checks at line 91 (`dXi/dPprime = 2/P`) and line 109 (`d(deltaP2)/dGprime = -2*P/(D0*Delta^2)`) are also real, and I verified by hand that the expected factor of 2 comes from the `2*Delta^2*P*g1` term in n2's numerator. The 2x2 sieve linear-algebra checks (A13, A14) are real conditional on the K1/H_even forms being correct. What does not hold up is (F2) the five "independence from variable that the formula doesn't contain" tautological assertions, (F3) the lack of any in-script derivation or upstream anchor for the literal `z0..n4` polynomials, and (F1) the absence of a Mathematica counterpart. Verdict: `findings`, three findings, no stop-cold. Output file is fresher than the script (script May 4, output May 11), so no `stale_output` finding.

## Self-test notes

I checked the variable-independence trap (audit prompt step 1): five of the existing assertions take `sp.diff(EXPR, VAR)` where `VAR` is not among `EXPR`'s symbols — that becomes F2. I checked parity (step 2): the only unbounded integrals are over `[0, oo)`, not symmetric, so the parity trap does not apply here. I checked the trivial-case substitution (step 3) for the surviving substantive checks: with `Q = D0 = Delta = 1`, `dXi/dPprime = 2/P` matches; with `Q = D0 = Delta = P = 1`, `d(deltaP2_der)/dGprime = -2*P/(D0*Delta^2) = -2`, and the g1 term in n2 contributes `-2*P*g1/Delta^2` which after the `D0^2/D0^3 = 1/D0` factor and `/mu1` gives `-2*P*Gx/(D0*Delta^2)`, derivative `-2*P/(D0*Delta^2)`. These match. I checked path specifications (step 4): the missing-mathematica directive targets `mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl`, matching the naming pattern of adjacent stages 010-012.
