---
unit_id: 003
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 003 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage003_bdg_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt`

mtime check: SymPy script Apr 21, output May 11 (fresh). Mathematica script Apr 21, output May 11 (fresh).

## What the script claims to verify

Per the module docstring, the unit claims to verify (a) the axisymmetric (qa, qL) wall + scalar BdG Lagrangian and its Euler-Lagrange equations, (b) exact elimination of the matter modes to produce the effective wall kernel `D_eff(omega) = K - omega^2 M - C (Omega_m^2 - omega^2 I)^{-1} C^T`, (c) the low-frequency renormalizations `K_eff = K - C Omega_m^{-2} C^T` and `M_eff = M + C Omega_m^{-4} C^T`, (d) the exact two-pole spectrum for one wall + one BdG mode and its perturbative shift `delta Omega^2 = -g^2/[M(varpi^2 - Omega^2)] + O(g^4)`, (e) channelwise grouped-P2 self-energies with trace/a2/b2 anisotropy invariants and isotropy preservation, and (f) the l=0 / l=2 selection rule from spherical harmonic orthogonality on an isotropic background. This is a checkpoint stage; both engines are required.

## Assertion inventory

| #   | Script      | Line     | Form                                                                                  | Anchored to claim? |
|-----|-------------|----------|---------------------------------------------------------------------------------------|--------------------|
| A1  | sympy       | 113-116  | `EL_qa.lhs + (expected EL terms) == 0`                                                | yes                |
| A2  | sympy       | 117-120  | `EL_qL.lhs + (expected EL terms) == 0`                                                | yes                |
| A3  | sympy       | 121-124  | `EL_xa.lhs + (expected EL terms) == 0`                                                | yes                |
| A4  | sympy       | 125-128  | `EL_xb.lhs + (expected EL terms) == 0`                                                | yes                |
| A5  | sympy       | 151      | `Deff - manual == 0` (matrix-algebra vs hand-expansion of the same formula)           | partial            |
| A6  | sympy       | 177      | `Deff_series - target_series == 0`                                                    | partial            |
| A7  | sympy       | 209-210  | `roots[i] - closed_form == 0` (sp.solve vs hand quadratic formula)                    | yes                |
| A8  | sympy       | 224-225  | wall-like / matter-like O(g^2) shift matches hand formula                             | yes                |
| A9  | sympy       | 299-302  | isotropic a2 / b2 / D20-D21 / D21-D22 vanish under iso substitution                   | partial            |
| A10 | sympy       | 328-331  | Y00-Y20, Y20-Y21c, Y20-Y22c cross integrals zero; Y20 norm = 1                        | partial            |
| B1  | mathematica | 70-72    | kinetic coefficient extraction (vqa^2, vqL^2, vqa vqL) matches Lagrangian             | no                 |
| B2  | mathematica | 91       | `dEff - manual == 0` (matrix-algebra vs hand-expansion of same formula)               | partial            |
| B3  | mathematica | 106      | `dEffSeries - (kEff - omega^2 mEff - omega^4 nEff) == 0`                              | partial            |
| B4  | mathematica | 121-122  | `dispersion[k -> m om2, w2 -> root] == 0` (root substitution into dispersion)         | yes                |
| B5  | mathematica | 123-124  | Vieta sum/product of roots                                                            | yes                |
| B6  | mathematica | 133-134  | wall-like / matter-like O(g^2) shifts                                                 | yes                |
| B7  | mathematica | 166-169  | isotropic a2 / b2 / D20-D21 / D21-D22 vanish under iso substitution                   | partial            |
| B8  | mathematica | 187-190  | Y00-Y20, Y20-Y21c, Y20-Y22c cross integrals; Y20 norm                                 | partial            |

"Partial" rows feed `insufficient_verification` (the assertion is true but does not exercise the full claim made by the docstring/comments) or `mathematica_transliteration` findings. "No" rows feed `tautological_check`.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:130-151`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:189-210` (one-mode dispersion analogue)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:78-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:115-122` (one-mode analogue)

**What's wrong:**

The script's docstring (lines 7-18 of the .py) says it verifies "exact elimination of the matter modes to obtain the effective wall kernels." The subbanner at line 130 says "Frequency-space elimination". But the elimination is never carried out in either engine. Both scripts simply *define*

```
Deff = Kmat - omega**2 * Mmat - Cmat * (Omat - omega**2 * I)**(-1) * Cmat.T
```

and then check that this matches a hand-expanded form (`manual` in both scripts). Both `Deff` and `manual` are constructed from the same hand-input matrices — the assertion `expect_zero("D0_eff - manual form", Deff - manual)` only verifies that SymPy's matrix inverse and multiplication agree with the textbook two-by-two formula. It does NOT verify that the formula `K - omega^2 M - C (Omega^2 - omega^2 I)^{-1} C^T` is what you actually get from algebraically eliminating xa, xb out of the Euler-Lagrange equations derived in I.1 (lines 113-128 of the .py).

The same pattern repeats in Section II: the dispersion `(K - M w2)(varpi2 - w2) - g^2 = 0` is asserted at line 197 (`.py`) / line 115 (`.wl`) without being derived from a single-mode Lagrangian. The Section I.1 EL machinery exists in the .py but is never plugged into Section I.2's elimination; Mathematica's Section I never even constructs the EL equations.

In short: I.1 verifies the EL equations exist with the right form; I.2 verifies the named formula is what SymPy/Mathematica say it is. There is no step that joins them — the elimination assertion is the link, and that link is missing.

**Why this matters:**

A sign error or off-by-one in the elimination — for example, dropping a +/- on the inverse, putting `+ c^2/(w^2 - omega^2)` instead of `- c^2/(w^2 - omega^2)` in the manual form (the Mathematica output at line 23 already shows `+c1a^2/(omega^2 - wa^2)` which is the same as `-c1a^2/(wa^2 - omega^2)` after sign-flipping the denominator) — would silently pass because both Deff and manual are derived from the same hand-input. The verification cannot detect a wrong sign in the published claim; it only detects sympy/mathematica computing matrix inverses incorrectly. The same hazard applies to the Section II dispersion: a wrong relative sign of `g^2` would not be caught.

**Required change:**

Add an explicit elimination step in both engines:

In SymPy, after building the EL equations at lines 107-110, take the frequency-space ansatz `qa, qL, xa, xb -> Qa, QL, Xa, Xb * exp(-I*omega*t)`, derive linear-in-amplitudes equations, solve the (Xa, Xb) block in terms of (Qa, QL), substitute back into the (Qa, QL) block, and compare the resulting effective 2x2 matrix to `Deff` (or to `manual`). Use the existing `EL_xa`, `EL_xb` `EL_qa`, `EL_qL` symbols, not hand-typed expressions. The new assertion should be `expect_zero("elimination derives D0_eff", Deff_from_elimination - Deff)`.

For Section II, similarly: construct a one-mode-each Lagrangian `L = M*qd^2/2 - K*q^2/2 + xd^2/2 - varpi^2*x^2/2 + g*q*x`, derive both EL equations, take the freq-space ansatz, solve for x in terms of q, substitute, and obtain `(K - M omega^2)(varpi^2 - omega^2) - g^2 == 0` as the derived condition, not a hand-input.

Mirror the new derivation in the `.wl` script (a structurally independent path is preferable — see F4 — but at minimum reproduce the elimination from the EL equations in both engines, not by transliterating one another).

**Verification:**

After Codex applies, the .py file should contain a new section between I.1 and I.2 (and a one-mode analogue in Section II) that uses the existing EL_qa / EL_xa expressions to algebraically eliminate the matter amplitudes and produces an expression that simplifies to `Deff` symbolically. The output transcript should include a new `PASS:`-level line such as `derived D0_eff - Deff = 0`. Same in the .wl output.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:64-72`

**What's wrong:**

The Mathematica script declares the Lagrangian at lines 54-59 with explicit `1/2 maa D[qa, t]^2 + maL D[qa, t] D[qL, t] + 1/2 mLL D[qL, t]^2 + ...`, then at lines 68 renames the time derivatives `D[qa,t] -> vqa` etc., and then asserts

```
expectZero["qa kinetic coefficient", Coefficient[lVel, vqa^2] - maa/2];
expectZero["qL kinetic coefficient", Coefficient[lVel, vqL^2] - mLL/2];
expectZero["qa-qL mixed kinetic coefficient", Coefficient[lVel, vqa vqL] - maL];
```

Each assertion extracts the coefficient of a term that was just literally written into the Lagrangian with that coefficient. `Coefficient[1/2 maa vqa^2 + ..., vqa^2]` is `maa/2` by construction. These three checks cannot fail.

**Why this matters:**

Three of the six "PASS" lines in the Mathematica transcript convey no information. They give the false impression that the Lagrangian structure has been independently verified, when nothing about that structure has been exercised. The Mathematica script (unlike the SymPy script) never derives Euler-Lagrange equations from `lRed`, so the only Lagrangian-level checks it makes are these three tautological ones.

**Required change:**

Replace the three tautological coefficient checks with actual Euler-Lagrange derivations, analogous to the SymPy script's lines 107-128. Use Mathematica's `D[L, q] - D[D[L, D[q, t]], t]` (or equivalent VariationalMethods package), then assert that each EL equation equals the expected `Maa qa'' + MaL qL'' + Kaa qa + KaL qL - c1a xa - c1b xb` (and analogues for qL, xa, xb). The new check should non-trivially fail under a wrong sign in `lRed`.

**Verification:**

After Codex applies, lines 64-72 of the .wl should compute EL equations from `lRed` (not extract literal coefficients), and the output transcript should report `qa equation = 0`, `qL equation = 0`, `xa equation = 0`, `xb equation = 0` — matching SymPy's Section I.1.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:312-334`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:171-190`

**What's wrong:**

Section IV's claim per the print statement (line 333 .py / line 79 of the .wl summary) is "the angular overlap is diagonal in (l,m)" — i.e. every cross-integral between distinct spherical harmonics vanishes and every diagonal one is unity. The script defines four harmonics (Y00, Y20, Y21c, Y22c) and checks only:

- `<Y00, Y20> = 0` (line 328 .py / line 187 .wl)
- `<Y20, Y21c> = 0`
- `<Y20, Y22c> = 0`
- `<Y20, Y20> = 1`

Three orthogonality combinations among the declared harmonics are not exercised:

- `<Y00, Y21c>` and `<Y00, Y22c>` — these are the relevant integrals for the *l=0 wall couples only to l=0 matter* selection rule between Y00 and the two grouped P2 channels c1c, c2c. If they happened to be nonzero, the claim "l=0 wall motion couples only to l=0 matter modes" would fail.
- `<Y21c, Y22c>` — m-orthogonality inside l=2, relevant to "grouped real l=2 wall motion couples channelwise inside l=2".
- `<Y00, Y00>`, `<Y21c, Y21c>`, `<Y22c, Y22c>` norms — only Y20's norm is verified.

The selection rule conclusion is broader than what the four checks actually establish.

**Why this matters:**

A typo in the normalization of Y21c or Y22c would not be detected; a wrong `cos(2 ph)` vs `sin(2 ph)` would not be detected unless it accidentally made the Y20-Y22c cross-integral nonzero. The claim "couples channelwise inside l=2" is exactly the statement that `<Y21c, Y22c> = 0`, but that integral isn't done.

**Required change:**

In both the .py and the .wl, add the missing cross-integrals and norms. Specifically:

In `.py` between lines 326 and 328, compute `I00_21c, I00_22c, I21c_22c, N00, N21c, N22c` analogously to the existing four, and add `expect_zero(...)` for each.

In `.wl` between lines 185 and 190, add analogous `i0021c, i0022c, i21c22c, norm00, norm21c, norm22c` integrals and `expectZero` calls.

**Verification:**

After Codex applies, the SymPy output should contain at least seven `... = 0` lines (or `norm - 1 = 0`) in Section IV (currently has four), and the Mathematica output should likewise have at least seven PASS lines in section IV.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:136-169` (Section III)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:171-190` (Section IV)

**What's wrong:**

The Mathematica script's Sections III and IV are line-by-line transliterations of the SymPy script's corresponding sections, with no independent derivation strategy.

Comparison of Section III, SymPy lines 244-302 vs Mathematica lines 143-169:

| SymPy                                                            | Mathematica                                                       |
|------------------------------------------------------------------|-------------------------------------------------------------------|
| `D20 = sp.simplify(K20 - M20*omega**2 - g20**2/(w20**2 - omega**2))` | `d20 = FullSimplify[k20 - m20 omega^2 - g20^2/(w20^2 - omega^2), ...]` |
| `D20s = sp.expand(sp.series(D20, omega, 0, 5).removeO())`         | `d20s = Expand[Normal[Series[d20, {omega, 0, 4}]]]`                |
| `d220 = sp.simplify(D20s.coeff(omega, 2))`                       | `d220 = FullSimplify[Coefficient[d20s, omega, 2], ...]`            |
| `d2bar = sp.simplify((d220 + 2*d221 + 2*d222) / 5)`               | `d2Bar = FullSimplify[(d220 + 2 d221 + 2 d222)/5, ...]`            |
| `a2 = sp.simplify((2*d220 - d221 - d222) / 10)`                   | `a2 = FullSimplify[(2 d220 - d221 - d222)/10, ...]`                |
| `b2 = sp.simplify((d221 - d222) / 2)`                             | `b2 = FullSimplify[(d221 - d222)/2, ...]`                          |
| `iso_subs = {K20: K2, K21: K2, K22: K2, M20: M2, ...}`            | `isoSubs = {k20 -> k2, k21 -> k2, k22 -> k2, m20 -> m2, ...}`      |
| `expect_zero("isotropic a2", a2.subs(iso_subs))`                  | `expectZero["isotropic a2", a2 /. isoSubs]`                        |

This is the same algebraic choreography, the same coefficient extractions, the same a2/b2/d2bar definitions, and the same isotropy substitution — re-written in Mathematica syntax. The Mathematica script makes no independent algebraic choice. If the SymPy author wrote a wrong combination (e.g. `(d220 + 2 d221 + 2 d222)/5` instead of the correct trace weight), Mathematica reproduces it verbatim and cannot detect the error.

Section IV is similarly a literal mirror: same four spherical harmonics with the same normalization conventions, same four integrals in the same order.

This violates the second-engine policy: both engines are supposed to derive the result independently from the physical premises.

**Why this matters:**

A second engine that mirrors the first is not a cross-check; it is duplicate algebra. The transliteration is most concerning in Section III, where the *definition* of the anisotropy invariants (`d2bar`, `a2`, `b2` as specific linear combinations of the channelwise d2 coefficients) is a methodological choice that should be re-derived from the underlying group-theoretic decomposition in the second engine, not copied from the first.

**Required change:**

Replace Section III in the .wl with an independent derivation. Two acceptable paths:

- (Preferred) Define the channelwise self-energy `D20, D21, D22` from a freshly written one-channel Lagrangian (one wall mode + one BdG mode each), derive each via the same elimination route F1 fixes, then construct `d2bar, a2, b2` via projection onto a representation-theoretic basis (e.g. by spelling out the Cartesian components of the real P2 channels and forming `Tr`, scalar anisotropy `a2`, and bi-axial anisotropy `b2` from explicit Y_20, Y_21c, Y_22c overlaps). The resulting symbolic forms should match SymPy's at the end, but the construction path should not be the same.
- (Acceptable) At minimum, change the order of operations: extract `d220, d221, d222` first, project to `d2bar, a2, b2` second, verify isotropy substitution third — but use a different basis or different coefficient labeling so that an off-by-one in SymPy's basis cannot be propagated by copy.

For Section IV, instead of computing the four hand-picked integrals, integrate the full pairwise overlap matrix over the four declared harmonics and assert that it equals the 4x4 identity. This requires no per-pair handpicking and structurally cannot be a transliteration of the SymPy version (which picks pairs).

**Verification:**

After Codex applies, the .wl Section III should contain at least one definition or assertion that does not appear (in any direct translation) in the .py — for instance, a representation-theoretic projection step. The Mathematica Section IV should compute a single overlap matrix rather than four separate integrals.

## Independent-derivation check (Mathematica)

The .wl script is **mixed**: Section I uses a different (and weaker — see F2) verification strategy than the .py (kinetic coefficient extraction vs Euler-Lagrange equations); Section II is mostly independent (Vieta sum/product plus direct dispersion-substitution, where the .py uses `sp.solve`); Sections III and IV are line-by-line transliterations of the .py (see F4 for the quoted side-by-side).

So the Mathematica script does *some* independent derivation work in Sections I and II, but Sections III and IV are not independent. This is sufficient to file F4 against III and IV, but the early sections do partially honor the second-engine policy.

## Engine cross-check

Both engines report exit 0 and all assertions PASS. The corresponding outputs:

- Section I, D_eff: SymPy prints the (Kaa - Maa omega^2 + c1a^2/(omega^2 - wa^2) + ...) form (after sympy's denominator-sign convention); Mathematica prints the same form (line 23 of the .wl output). Algebraically identical modulo trivial sign convention `1/(w^2-omega^2) = -1/(omega^2 - w^2)`.
- Section I, series expansion: both engines produce the same Keff/Meff/Neff structure, both report `series match = 0`.
- Section II, roots: SymPy expresses as `Om2/2 + varpi2/2 ± sqrt(M*Om2^2 - 2 M Om2 varpi2 + M varpi2^2 + 4 g^2)/(2 sqrt(M))`; Mathematica as `(om2 + varpi2 ± sqrt(4 g^2/m + (om2 - varpi2)^2))/2`. These differ in cosmetic form but are equal once the sqrt is factored — both pass their respective `expectZero` against the documented closed form `(Om2 + varpi2 ± sqrt((Om2 - varpi2)^2 + 4 g^2/M))/2`.
- Section II, shifts: identical numerical/symbolic forms `Om2 - eps^2 g^2/(m delta)` and `delta + om2 + eps^2 g^2/(m delta)` in both engines.
- Section III, channelwise pieces: identical d2bar, a2, b2 expressions (modulo notation), identical isotropy results.
- Section IV: identical zero results for the four integrals.

Engines agree. No `engine_disagreement`.

## Verdict justification

Four findings: one medium-severity insufficient_verification (F1 — elimination from EL not derived in either engine; the heart of the claim is asserted, not produced), one medium-severity mathematica_transliteration (F4 — Sections III and IV are line-by-line copies), one low tautological_check (F2 — Mathematica's only Lagrangian-level checks are coefficient extractions of just-written terms), and one low insufficient_verification (F3 — Section IV omits several relevant orthogonality integrals).

Attacks that the script holds up against: the EL equations in SymPy Section I.1 are real and would catch a sign error in `lRed`. The Vieta and dispersion-substitution checks in Mathematica Section II are real and would catch a sign error in the closed-form roots. The perturbative O(g^2) wall/matter shifts in both engines are real Taylor-expansion checks against the predicted shift formula and would catch a factor-of-2 or sign error in the published shift. The freshly typed isotropy substitution in Section III is structurally constrained to vanish — not an attack survivor, but tagged informationally rather than as a separate finding.

No `stop_cold` flag. The findings are all script-level (verification scope and second-engine independence), not internal mathematical inconsistencies, so they're fixable in place; nothing about the underlying physical claim is contradicted by the script itself. Engines agree; outputs are fresh.
