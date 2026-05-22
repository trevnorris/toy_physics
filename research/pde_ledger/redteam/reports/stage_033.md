---
unit_id: 033
batch: II.1
auditor_model: claude-opus-4-7
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

# Audit unit 033 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt`

## What the script claims to verify

The pair of scripts assemble a "selected-branch" normalization product `N_-(alpha0) = beta0 * s_-^2 / (kappa0^2 * lambda_-)` built from finite-throat overlap constants `kappa0^2 = 8/pi^2`, `kappa1^2 = 16/(9 pi^2)` and the auxiliary radical `R = sqrt((DeltaK + alpha0*delta_kappa)^2 + 4 alpha0^2 Kprod)`. They claim to verify five physics results: (i) an exact monotonicity identity `dN/dalpha = beta0(2 s_- s_-' lambda_- + s_-^3)/(kappa0^2 lambda_-^2)`, which is equivalent to the identity `dlambda_-/dalpha = -s_-`; (ii) the closed-form finite-throat critical loading `alpha_crit = 9 pi^2 A (A+DeltaK)/(8(11A+9DeltaK))`; (iii) `N_-(0) = beta0 kappa0^2 / A`; (iv) the weak-loading slope at alpha0 = 0 in both a generic-coupling and a closed-form expression; and (v) the microscopic K0-onset value obtained by solving `N_-(0)|_mic = NQ`. Stage 16.6 then claims the fully substituted stability gate `alpha_crit(mic) - alpha_0(mic)` can be written as `gate_num/gate_den` with `gate_den = 8 varpi^2 Omega_U^2 Delta0 (11 A + 9 DeltaK)`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 54 | `simplify(dN - dN_formula) == 0` | yes (identity dlambda/dalpha = -s) |
| A2 | sympy | 60 | `simplify(alpha_crit - alpha_crit_target) == 0` | partial (algebraic rearrangement of hand-stated forms; original `alpha_crit` is not derived in-script) |
| A3 | sympy | 65 | `simplify(N0 - beta0 kappa0^2/A) == 0` | yes |
| A4 | sympy | 73 | `simplify(coef1 - coef1_target) == 0` | yes |
| A5 | sympy | 74 | `simplify(coef1 - coef1_target_closed) == 0` | yes |
| A6 | sympy | 97-100 | `simplify(K0_onset - target) == 0` (with K0_onset from sp.solve) | yes |
| A7 | sympy | 112-115 | `simplify(alpha_crit_mic - alpha0_mic - gate_num/gate_den) == 0` | **no — tautological** (gate_num is defined as `(alpha_crit_mic - alpha0_mic)*gate_den`) |
| M1 | mathematica | 60 | `FullSimplify[dN - dNFormula] == 0` | yes (mirror of A1) |
| M2 | mathematica | 65 | `FullSimplify[alphaCrit - alphaCritClosed] == 0` | partial (mirror of A2) |
| M3 | mathematica | 69 | `FullSimplify[n0 - beta0*kappa0Sq/A] == 0` | yes |
| M4 | mathematica | 81 | `FullSimplify[coef1 - coef1Target] == 0` | yes |
| M5 | mathematica | 82 | `FullSimplify[coef1 - coef1Closed] == 0` | yes |
| M6 | mathematica | 105 | `FullSimplify[(n0Mic /. K0 -> k0Onset) - NQ] == 0` | yes (substitution-back consistency check) |
| M7 | mathematica | 106-109 | `FullSimplify[k0Onset - (gU^2/OmegaU^2 + kappa0Sq*chi^2/(NQ*delta0^2))] == 0` | **no — tautological** (k0Onset is defined as exactly the target expression at line 95) |
| M8 | mathematica | 116-119 | `FullSimplify[alphaCritMic - alpha0Mic - gateNum/gateDen] == 0` | **no — tautological** (mirror of A7) |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:104-115`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:112-119`

**What's wrong:**
The Stage 16.6 "fully substituted microscopic stability gate" assertion is algebraically guaranteed by construction. The SymPy script defines

```python
gate_num = sp.simplify(sp.expand((alpha_crit_mic - alpha0_mic) * gate_den))
```

at line 105, then at lines 112-115 asserts

```python
expect_zero(
    "alpha_crit(mic) - alpha_0(mic) - gate_num/gate_den",
    alpha_crit_mic - alpha0_mic - gate_num / gate_den,
)
```

By substitution `gate_num/gate_den == (alpha_crit_mic - alpha0_mic) * gate_den / gate_den`, so the residual is identically zero in *any* commutative algebra — it would pass even if `alpha_crit_mic` or `alpha0_mic` were arbitrary symbols with no physical content. The Mathematica script (lines 113, 116-119) reproduces the same tautological structure with `gateNum = FullSimplify[Expand[(alphaCritMic - alpha0Mic)*gateDen], ...]` followed by `expectZero["alpha_crit(mic) - alpha_0(mic) - gate_num/gate_den", alphaCritMic - alpha0Mic - gateNum/gateDen]`. The transcript banner promises a verified closed-form numerator-over-denominator decomposition of the gate, but the assertion only confirms that division undoes multiplication.

**Why this matters:**
Stage 16.6 is the unit's final and most substantive claim: that the fully substituted stability gate has the specific factored form `8 varpi^2 Omega_U^2 Delta0 (11 A_mic + 9 DeltaK)` in its denominator (and the long polynomial in its numerator shown in the saved output). A genuine error in either `alpha_crit_mic` or `alpha0_mic` would not be caught — the script would still report PASS. There is therefore no verification at all of the gate structure beyond the printed transcript.

**Required change:**
Replace the constructed `gate_num`/`gateNum` with an *independent* check. The substantive claim is that `together(alpha_crit_mic - alpha0_mic)` has denominator equal to `gate_den` (after clearing common factors and overall sign). Verify this by extracting the denominator of the simplified difference and checking it equals `gate_den` up to a multiplicative rational constant; the numerator can then be defined as `simplify((alpha_crit_mic - alpha0_mic) * gate_den)` *for printing only*, with a separate assertion that this numerator is a polynomial in `{K0, gU, gB, gR, gW, OmegaU, OmegaW, varpi, DeltaK}` (i.e. has trivial denominator after `together`). See directive F1 for exact instructions.

**Verification:**
After the fix, the SymPy assertion at lines 112-115 must reference a target denominator (e.g. `Poly(together(alpha_crit_mic - alpha0_mic).as_numer_denom()[1])`) compared against `gate_den` directly, *not* a self-constructed `gate_num/gate_den`. The Mathematica assertion at lines 116-119 must do the analogous denominator extraction with `Denominator[Together[...]]`. The new check must be capable of failing if `alpha_crit_mic` or `alpha0_mic` is perturbed.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:95,106-109`

**What's wrong:**
In Mathematica the K0_onset value is *defined* directly at line 95:

```
k0Onset = FullSimplify[gU^2/OmegaU^2 + kappa0Sq*chi^2/(NQ*delta0^2), Assumptions -> $Assumptions];
```

Then at lines 106-109 the script asserts:

```
expectZero[
  "K0_onset - [gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)]",
  k0Onset - (gU^2/OmegaU^2 + kappa0Sq*chi^2/(NQ*delta0^2))
];
```

The asserted residual is `expr - expr` where `expr` is the line-95 right-hand side — identically zero. (The companion SymPy script at line 87 does the real work: `K0_onset = sp.simplify(sp.solve(sp.Eq(N0_mic, NQ), K0)[0])`, then compares the *solved* value against the same target at line 99.) The substantive Mathematica check is at line 105 (`N_-(0) at K0_onset - NQ`), which does verify consistency by back-substitution; but the assertion at lines 106-109 adds nothing and masquerades as an independent verification of the K0_onset formula.

**Why this matters:**
Reading the Mathematica transcript, one would conclude that K0_onset was derived twice (once by formula, once by check). In fact, only line 105 is substantive; the line 106-109 check is dead weight. If a future maintainer edits the line-95 `k0Onset` expression, the line 106-109 check will silently track the edit and never alert. The second-engine policy expects both engines to independently *derive* the result; here Mathematica hardcodes it while SymPy solves.

**Required change:**
Replace the hardcoded `k0Onset` at line 95 with a `Solve` call analogous to SymPy's, then keep the comparison at lines 106-109 as a non-trivial check. Specifically:

```
k0Onset = K0 /. First@Solve[n0Mic == NQ, K0];
k0Onset = FullSimplify[k0Onset, Assumptions -> $Assumptions];
```

Then the line 106-109 assertion checks that the solved value equals the closed form, which is the intended physics claim.

**Verification:**
After the fix, line 95 must invoke `Solve[n0Mic == NQ, K0]` (or equivalent inversion). The line 106-109 assertion stays unchanged; its content is then non-tautological. The line 105 assertion remains an additional cross-check.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:33-119`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:33-115`

**What's wrong:**
The Mathematica script is a line-by-line algebraic mirror of the SymPy script. Compare:

Constants block:
```python
# sympy lines 33-37
kappa0_sq = sp.Rational(8) / sp.pi**2
kappa1_sq = sp.Rational(16) / (9 * sp.pi**2)
sigma = sp.simplify(kappa0_sq + kappa1_sq)
delta_kappa = sp.simplify(kappa0_sq - kappa1_sq)
Kprod = sp.simplify(kappa0_sq * kappa1_sq)
```
```
(* mathematica lines 33-37 *)
kappa0Sq = 8/Pi^2;
kappa1Sq = 16/(9*Pi^2);
sigma = FullSimplify[kappa0Sq + kappa1Sq, ...];
deltaKappa = FullSimplify[kappa0Sq - kappa1Sq, ...];
kProd = FullSimplify[kappa0Sq*kappa1Sq, ...];
```

Branch definitions:
```python
# sympy lines 41-44
R = sp.sqrt((DeltaK + alpha0 * delta_kappa)**2 + 4 * alpha0**2 * Kprod)
lambda_minus = sp.simplify((2 * A + DeltaK - alpha0 * sigma - R) / 2)
s_minus = sp.simplify(sp.Rational(1, 2) * (sigma + ((DeltaK + alpha0 * delta_kappa) * delta_kappa + 4 * alpha0 * Kprod) / R))
Nminus = sp.simplify(beta0 * s_minus**2 / (kappa0_sq * lambda_minus))
```
```
(* mathematica lines 39-48 *)
r = FullSimplify[Sqrt[(DeltaK + alpha0*deltaKappa)^2 + 4*alpha0^2*kProd], ...];
lambdaMinus = FullSimplify[(2*A + DeltaK - alpha0*sigma - r)/2, ...];
sMinus = FullSimplify[1/2*(sigma + ((DeltaK + alpha0*deltaKappa)*deltaKappa + 4*alpha0*kProd)/r), ...];
nMinus = FullSimplify[beta0*sMinus^2/(kappa0Sq*lambdaMinus), ...];
```

dN/dalpha identity:
```python
# sympy lines 51-54
ds = sp.simplify(sp.diff(s_minus, alpha0))
dN = sp.simplify(sp.diff(Nminus, alpha0))
dN_formula = sp.simplify(beta0 * (2 * s_minus * ds * lambda_minus + s_minus**3) / (kappa0_sq * lambda_minus**2))
expect_zero("dN/dalpha - monotonicity formula", dN - dN_formula)
```
```
(* mathematica lines 54-60 *)
ds = FullSimplify[D[sMinus, alpha0], ...];
dN = FullSimplify[D[nMinus, alpha0], ...];
dNFormula = FullSimplify[beta0*(2*sMinus*ds*lambdaMinus + sMinus^3)/(kappa0Sq*lambdaMinus^2), ...];
expectZero["dN/dalpha - monotonicity formula", dN - dNFormula];
```

The same correspondence holds for stages 16.2, 16.4, the microscopic-coupling block (`aMic`, `delta0`, `chi`, `beta0Mic`, `alpha0Mic` ↔ `A_mic`, `Delta0`, `Chi`, `beta0_mic`, `alpha0_mic`), and the gate block (`gateNum = FullSimplify[Expand[...]] ↔ gate_num = sp.simplify(sp.expand(...))`). The Mathematica engine never derives any quantity from a physically distinct starting point — variable names are mechanical renames, intermediate steps line up one-for-one. The only structural divergence (`k0Onset` hardcoded vs `K0_onset` solved) is the wrong direction: SymPy does the harder thing, Mathematica does the easier one (F2).

**Why this matters:**
The second-engine policy exists so that algebra mistakes in one engine's transcription get caught by an independently-derived counterpart. A pure transliteration cannot catch any SymPy algebra error, because the same algebra is being re-typed. The gate-identity tautology (F1) is the perfect example: the same wrong check appears in both engines because they share the same construction logic. A genuine cross-check would have caught the tautology.

**Required change:**
Insert at least one *independent* derivation path into the Mathematica script — typically a numerical-substitution cross-check at randomized parameter values that bypasses the analytic identity entirely. See directive F3 for exact instructions.

**Verification:**
After the fix, the Mathematica script must contain at least one assertion that exercises the Stage 16.1 monotonicity identity and the Stage 16.6 gate identity at concrete numerical values (e.g. `{A->2, DeltaK->1, alpha0->1/3, beta0->1, varpi->1, OmegaU->1, OmegaW->1, gB->1/2, gU->1/3, gW->1/4, gR->1/5, K0->3, NQ->1}` plus a second distinct numeric assignment), checking the residual is zero to machine precision in floating-point. This is structurally distinct from the SymPy approach (closed-form simplification) and would catch algebra errors.

## Independent-derivation check (Mathematica)

The `.wl` file is **not** an independent derivation. The constants, the radical `r`, `lambdaMinus`, `sMinus`, `nMinus`, the monotonicity formula `dNFormula`, the microscopic substitutions, and the gate construction all mirror SymPy's algebraic choreography one-for-one. The only structural divergence is at line 95 (`k0Onset` hardcoded vs `sp.solve` in SymPy) — and that divergence makes Mathematica's check at lines 106-109 tautological (F2). Three quoted blocks above show the parallel structure; this satisfies the `mathematica_transliteration` finding category.

## Engine cross-check

Both engines report PASS on every assertion they share. The numerical outputs printed (e.g. `N_-(0) = 8*beta0/(pi^2*A)`, `alpha_crit = 9*pi^2*A*(A+DeltaK)/(8*(11*A + 9*DeltaK))`, `K0_onset = gU^2/OmegaU^2 + 648*pi^2*(...)^2/(NQ*(9*pi^2*OmegaU^2*OmegaW^2 - 88*gR^2)^2)`) match between SymPy (lines 40, 46, 65 of output) and Mathematica (lines 18, 21, 35 of output) up to formatting differences (`8*beta0/(pi^2*A)` vs `(8*beta0)/(A*Pi^2)`). No numeric or sign disagreement detected. Outputs are fresh (both `.txt` files mtime newer than corresponding script mtimes: sympy script Apr 3 12:05 / output May 11 12:41; mathematica script May 11 11:56 / output May 11 12:47).

## Verdict justification

The unit's first five claims (monotonicity, alpha_crit closed form, N_-(0), weak-loading slope generic + closed, K0_onset by solve) hold up under attack in SymPy: I tried to construct a hidden tautology, a missing branch, or a symbol-domain error and found none — the derivative-formula check at line 54 is a genuine identity (encoding `dlambda/dalpha = -s_minus`), the K0_onset solve at line 87 is a real inversion, and the alpha=0 substitution at lines 63-65 is a direct evaluation. Stage 16.2 has the modest concern that the `alpha_crit` formula is not derived in-script (both forms are hand-stated and then equated), but the algebraic identity between them is itself non-trivial after `kappa` substitution, so I do not file a separate hardcoded_result finding. Stage 16.6 is the failure: the gate-identity assertion is tautological by construction in both engines (F1). Mathematica's K0_onset block has a separate tautology (F2) that the SymPy script avoids. The whole Mathematica script is a transliteration (F3) which is *why* F1 reproduces verbatim in both engines. Verdict: `findings`, no stop-cold (all three findings are fixable mechanically).

## Self-test notes

I checked the four required traps. (1) Variable-independence: every `sp.diff(EXPR, alpha0)` and `D[expr, alpha0]` in the scripts has `alpha0` genuinely present in `EXPR` (`s_minus` and `Nminus` both contain `alpha0` via `R`); no zero-derivative-by-construction artifacts. (2) Symmetry/parity: no integrals in this unit; trap N/A. (3) Trivial-case pre-check: at alpha0=0 we get `R = DeltaK`, `lambda_- = A`, `s_- = (sigma+delta_kappa)/2 = kappa0^2`, so `N_-(0) = beta0*kappa0^2/A` matches the script's claim; I also verified the F1 tautology by hand-substituting `gate_num/gate_den` (gate_den cancels exactly, leaving `alpha_crit_mic - alpha0_mic`), confirming the residual is structurally zero. (4) Path specifications: both directive targets use the canonical paths (`scripts/...sympy_audit.py`, `mathematica/...mathematica_audit.wl`) — no new files are introduced.
