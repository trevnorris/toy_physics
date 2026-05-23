---
unit_id: 064
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 064 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.txt`

## What the script claims to verify

Per the docstring, the script claims to verify five things: (1) a closure law `chi_sigma(y) = g_phi * chi_phi(y) / H(y)` from a local static linear-response argument; (2) that the overlap invariants reduce to `O = g_phi*I1`, `N_ss = g_phi^2*I2` with `I1 = ∫chi_phi^2/H` and `I2 = ∫chi_phi^2/H^2`; (3) that the coherence factor is `C^2 = I1^2/(N_pp*I2)`; (4) that in the constant-compressibility limit `H = H_w`, the coherence becomes exactly 1 and the gain reproduces the Stage-45/46 formula `G_eq = g_phi^2 N_pp/(K_X H_w)`; (5) the discrete two-point Cauchy gap identity; and (6) an eliminated-source softening identity `F_eff = (1/2)(K_X - Lambda^2/Theta) phi^2` with `Lambda^2/Theta = g_phi^2 N_pp/H_w` under the matched branch.

## Assertion inventory

| #  | Script       | Line | Form                                                                                                  | Anchored to claim? |
|----|--------------|------|-------------------------------------------------------------------------------------------------------|--------------------|
| A1 | sympy        | 77   | `expect_zero("matched-layer coherence", C2_const - 1)`                                                | no (tautological under `const_subs`) |
| A2 | sympy        | 78   | `expect_zero("matched-layer gain ...", Geq_const - g_phi**2 * Npp / (KX * Hw))`                       | no (tautological under `const_subs`) |
| A3 | sympy        | 95   | `expect_zero("two-point Cauchy gap identity", gap_disc - gap_expected)`                               | yes (non-trivial algebraic identity) |
| A4 | sympy        | 107  | `expect_zero("effective support softening", F_eff - (1/2)(K_X - Lambda^2/Theta) phi^2)`               | yes (real derivation from F) |
| A5 | sympy        | 118  | `expect_zero("equilibrium softening equals g_phi^2 I1", soft_eq - g_phi**2 * (Npp/Hw))`               | no (tautological under `const_subs`) |
| A6 | mathematica  | 49   | `expectZero["matched-layer coherence", c2Const - 1]`                                                  | no (mirror of A1) |
| A7 | mathematica  | 50   | `expectZero["matched-layer gain ...", gEqConst - gPhi^2*npp/(kX*hw)]`                                 | no (mirror of A2) |
| A8 | mathematica  | 64   | `expectZero["two-point Cauchy gap identity", gapDisc - gapExpected]`                                  | yes (mirror of A3, but holds) |
| A9 | mathematica  | 77   | `expectZero["effective support softening", fEff - 1/2*(kX - lambda^2/theta)*phi^2]`                   | yes (mirror of A4) |
| A10| mathematica  | 84   | `expectZero["equilibrium softening equals g_phi^2 I1", softEq - gPhi^2*(npp/hw)]`                     | no (mirror of A5) |

Of the 10 assertions, six (A1, A2, A5, A6, A7, A10) are guaranteed-true by direct substitution; only A3/A8 and A4/A9 carry real algebraic content.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:71-78`
- `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:43-50`

**What's wrong:**
The "matched-layer coherence = 1" and "matched-layer gain" assertions (A1/A2 and their Mathematica mirrors A6/A7) are forced true by the substitution dictionary `const_subs = {I1: Npp/Hw, I2: Npp/Hw**2}`. Tracing the algebra:

```
C2 = (g_phi*I1)^2 / (Npp * g_phi^2 * I2) = I1^2 / (Npp * I2)
C2.subs(I1=Npp/Hw, I2=Npp/Hw^2) = (Npp/Hw)^2 / (Npp * Npp/Hw^2)
                                = (Npp^2/Hw^2) / (Npp^2/Hw^2) = 1
```

The result `C2_const = 1` is an arithmetic consequence of choosing `I2 = I1^2/Npp` — and the substitution dictionary *is* exactly that choice in disguise (`(Npp/Hw)^2/Npp = Npp/Hw^2`). The physical content the docstring claims to verify — that "the matched-layer limit `H(y) = H_w` forces `I1 = Npp/Hw` and `I2 = Npp/Hw^2`" — is asserted by the substitution rather than derived. No `chi_phi(y)` profile is ever defined, no `H(y)` function is ever defined, and `Npp = ∫chi_phi^2` is never tied to the discrete `Npp` symbol.

The gain assertion `Geq_const - g_phi**2 * Npp / (KX * Hw) == 0` is the same trick one step later: `Geq = g_phi^2 * I1 / KX` and `I1` is substituted to `Npp/Hw`, so the RHS is just `g_phi^2 * (Npp/Hw) / KX`, identical to the subtracted constant. The assertion cannot fail.

**Why this matters:**
Two of the five "main checks" the docstring advertises (coherence reaching 1 in the matched-layer limit, gain reproducing Stage-45/46) are verified only by a tautology. If `Hw` were typoed as `Hw^2` in the substitution, the assertion would still pass with a corresponding redefinition. The script provides no defense against a wrong derivation upstream.

**Required change:**
Add an explicit, parameterised concrete-profile derivation. Introduce a smooth `chi_phi(y)` (e.g., a Gaussian `chi_phi = exp(-y^2/(2*L^2))`) and a constant `H(y) = Hw`, define `Npp_int = sp.integrate(chi_phi**2, (y, -sp.oo, sp.oo))`, `I1_int = sp.integrate(chi_phi**2 / Hw, (y, -sp.oo, sp.oo))`, `I2_int = sp.integrate(chi_phi**2 / Hw**2, (y, -sp.oo, sp.oo))`, then `expect_zero` on `I1_int - Npp_int/Hw` and `I2_int - Npp_int/Hw**2`. This verifies that the matched-layer reductions follow from the integral definitions — not from a hand-substitution. Only after that step should the existing `C2_const - 1` and `Geq_const - ...` checks run (and at that point they would consume `Npp_int, I1_int, I2_int` rather than abstract symbols). Mirror the same derivation in the `.wl` script with `Integrate[Exp[-y^2/L^2], {y, -Infinity, Infinity}]`-style calls.

**Verification:**
New checks `matched-layer I1 reduction = 0` and `matched-layer I2 reduction = 0` should appear in the sympy output ahead of "matched-layer coherence", with a printed concrete value of `Npp_int = sqrt(pi)*L` (for a unit-amplitude Gaussian) or similar.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:114-118`
- `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:79-84`

**What's wrong:**
The "equilibrium softening equals g_phi^2 I1" assertion (A5/A10) is tautological by the same `const_subs` mechanism. Trace:

```
Theta_eq = Hw * Nss.subs(const_subs) = Hw * g_phi^2 * (Npp/Hw^2) = g_phi^2 * Npp/Hw
Lambda_eq = g_phi * Osp.subs(const_subs) = g_phi * g_phi * (Npp/Hw) = g_phi^2 * Npp/Hw
soft_eq = Lambda_eq^2 / Theta_eq = (g_phi^2*Npp/Hw)^2 / (g_phi^2*Npp/Hw) = g_phi^2 * Npp/Hw
```

`soft_eq - g_phi**2 * (Npp/Hw) == 0` is therefore an algebraic identity once `const_subs` is applied; no physics is being checked. The docstring says this verifies that the "equilibrium softening equals `g_phi^2 * I1`" — but again `I1` has been hand-substituted to `Npp/Hw`, so the equality `soft_eq = g_phi^2 * I1` is just `g_phi^2 * Npp/Hw = g_phi^2 * (Npp/Hw)`.

**Why this matters:**
The script presents this as the "equilibrium softening" identity tying `Lambda^2/Theta` back to the `g_phi^2 I1` overlap. A genuine check would compare `Lambda^2/Theta` derived from the abstract `Osp, Nss` (before substitution) against `g_phi^2 I1` (also before substitution) using only the closure relations.

**Required change:**
Replace lines 114-118 of the sympy script (and the mirroring block in the .wl) with a check at the abstract level:
```
Theta_abs = Hw * Nss   # = Hw * g_phi^2 * I2
Lambda_abs = g_phi * Osp   # = g_phi^2 * I1
soft_abs = sp.simplify(Lambda_abs**2 / Theta_abs)
# In the matched layer, the closure gives I2 = I1^2/Npp; only THEN should it reduce to g_phi^2 * I1.
soft_matched = sp.simplify(soft_abs.subs(I2, I1**2 / Npp))
expect_zero("equilibrium softening (closure form)", soft_matched - g_phi**2 * I1)
```
This separates two facts the current code conflates: (i) `Lambda^2/Theta = g_phi^2 * I1^2 / (Npp * I2)` always, and (ii) the matched-layer closure forces `I2 = I1^2/Npp`. The assertion `soft_matched - g_phi^2 * I1 == 0` then non-trivially exercises the closure.

**Verification:**
A new printed line `soft_abs = I1**2*g_phi**2/(I2*Npp)` (or equivalent) must appear before the assertion, demonstrating that the pre-substitution form is non-trivial; the assertion then fires after a single `subs(I2, I1**2/Npp)`.

### F3 — insufficient_verification

**Severity:** high
**Files:**
- `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:30-50`
- `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:26-41`

**What's wrong:**
The docstring's first and second main checks are: (1) the closure law `chi_sigma(y) = g_phi * chi_phi(y) / H(y)` arising from a "local static linear-response closure"; (2) the overlap invariants `O = g_phi * I1`, `N_ss = g_phi^2 * I2` with `I1 = ∫chi_phi^2/H, I2 = ∫chi_phi^2/H^2`. The script declares `I1, I2, Npp, g_phi, KX, Hw` as abstract positive symbols and writes:

```python
Osp = sp.simplify(g_phi * I1)
Nss = sp.simplify(g_phi**2 * I2)
```

These are *definitional* — `sp.simplify(g_phi*I1)` does no algebraic work since the product is already simplified. The "closure law" itself (chi_sigma = g_phi chi_phi/H) is never written down as a sympy expression, never verified by minimising a free-energy functional in sigma, never tied to any local linear response argument. Neither integral expression for `I1` or `I2` ever appears. The audit therefore does not exercise main checks 1-2 at all; the assertions only fire on later algebraic manipulations of the abstract symbols.

**Why this matters:**
The script presents itself as verifying the physics closure that *leads to* the overlap forms. In fact it only verifies arithmetic on the assumed forms. A reader of the saved output would conclude the closure was checked — it wasn't.

**Required change:**
Add an explicit derivation block before line 49 of the sympy script:
```
# Local static closure: minimize (1/2) H(y) sigma^2 - g_phi chi_phi(y) sigma over sigma
# at fixed y.  H(y) > 0, real.
y = sp.symbols('y', real=True)
chi_phi_sym = sp.Function('chi_phi')(y)
H_sym = sp.Function('H')(y)
sigma_sym = sp.symbols('sigma_local', real=True)
local_F = sp.Rational(1,2) * H_sym * sigma_sym**2 - g_phi * chi_phi_sym * sigma_sym
chi_sigma_pred = sp.solve(sp.diff(local_F, sigma_sym), sigma_sym)[0]
expect_zero("closure law chi_sigma = g_phi chi_phi/H", chi_sigma_pred - g_phi*chi_phi_sym/H_sym)
```
Then, with a concrete `chi_phi = exp(-y^2/(2 L^2))` and `H = Hw` constant, define `Npp_int, I1_int, I2_int` as integrals and `Osp_int = sp.integrate(chi_phi*chi_sigma_pred.subs(...), y)` and confirm `Osp_int - g_phi*I1_int == 0` etc. Mirror the same construction in the `.wl` (using `D[]` and `Integrate[]`).

**Verification:**
A line `closure law chi_sigma = g_phi chi_phi/H = 0` must appear in the new sympy output, followed by `Osp_int - g_phi*I1_int = 0` and `Nss_int - g_phi^2*I2_int = 0` with concrete integral values printed.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:26-84`
- (against) `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:49-118`

**What's wrong:**
The Mathematica script is a line-by-line port of the SymPy script rather than an independent derivation. Quoting corresponding blocks:

SymPy lines 57-65:
```
Osp = sp.simplify(g_phi * I1)
Nss = sp.simplify(g_phi**2 * I2)
C2 = sp.simplify(Osp**2 / (Npp * Nss))
Geq = sp.simplify(g_phi**2 * I1 / KX)
```
Mathematica lines 33-36 (same order, same lowercased names):
```
osp = FullSimplify[gPhi*i1, ...];
nss = FullSimplify[gPhi^2*i2, ...];
c2 = FullSimplify[osp^2/(npp*nss), ...];
gEq = FullSimplify[gPhi^2*i1/kX, ...];
```

SymPy lines 71-72:
```
const_subs = {I1: Npp / Hw, I2: Npp / Hw**2}
C2_const = sp.simplify(C2.subs(const_subs))
```
Mathematica lines 43-44:
```
constSubs = {i1 -> npp/hw, i2 -> npp/hw^2};
c2Const = FullSimplify[c2 /. constSubs, ...];
```

SymPy lines 82-91 vs. Mathematica lines 52-60: identical `w1, w2, H1, H2` choreography, identical `Npp_disc/i1Disc/i2Disc` definitions, identical `expected_gap` form.

SymPy lines 100-110 vs. Mathematica lines 71-77: identical free-energy `F = (Theta/2)sigma^2 - Lambda sigma phi + (KX/2) phi^2`, identical Solve step, identical target form.

This mirror-image structure violates the second-engine policy: when both engines make the same definitional choices (e.g., the `const_subs` dictionary appears verbatim), they cannot independently catch a derivation error upstream.

**Why this matters:**
The point of a second engine is *adversarial cross-check* on the algebra. If the .wl makes every choice the .py makes, the second engine's `PASS` adds no information beyond confirming that Mathematica's algebra agrees with sympy's on the same inputs — which we already know.

**Required change:**
Refactor the `.wl` so its derivation proceeds independently. Concretely:
- In the matched-layer block, do not use a `constSubs` dictionary. Instead define a concrete `chiPhi[y_] := Exp[-y^2/(2 L^2)]` and `H[y_] := hw`, then `nppInt = Integrate[chiPhi[y]^2, {y, -Infinity, Infinity}]`, `i1Int = Integrate[chiPhi[y]^2/H[y], ...]`, `i2Int = Integrate[chiPhi[y]^2/H[y]^2, ...]`. Then `c2Const` is computed from these integrals, not from a substitution dict.
- In the softening block, derive `Lambda^2/Theta` from a `Reduce` of `D[f, sigma] == 0` followed by `Eliminate[..., sigma]` rather than reproducing the sympy `Solve`/`subs` pattern.

**Verification:**
The refactored `.wl` should print concrete integral values (e.g., `nppInt = Sqrt[Pi]*L`) and the matched-layer `c2Const` and `gEqConst` checks should consume those integrals — not symbolic `i1/i2`. The directive for F1 and F3 dovetails with this requirement.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** an independent re-derivation; it is a transliteration of the SymPy script (see F4 for paired excerpts). Same variable order, same `constSubs == const_subs`, same five-assertion structure, same labels on `expectZero`. Both engines arrive at the same `PASS` lines via the same algebraic path; the cross-check therefore degenerates to "two CASes agree on substitution arithmetic."

## Engine cross-check

Both engines print the same numerical and symbolic results:
- `O_(sigma phi) = I1*g_phi` (py) / `gPhi*i1` (wl) — agree
- `C_(sigma phi)^2 = I1**2/(I2*N_pp)` (py) / `i1^2/(i2*npp)` (wl) — agree
- `G_eq | H=const = N_pp*g_phi**2/(H_w*K_X)` (py) / `(gPhi^2*npp)/(hw*kX)` (wl) — agree
- `N_pp I2 - I1^2 = w1*w2*(H1-H2)^2/(H1^2*H2^2)` — agree
- `F_eff = phi^2*(K_X*Theta - Lambda^2)/(2*Theta)` (py) / `(phi^2*(-lambda^2 + kX*theta))/(2*theta)` (wl) — agree
- `(Lambda^2/Theta)_eq = N_pp*g_phi**2/H_w` — agree

The engines agree, but because the .wl is a transliteration of the .py, the agreement is mechanical, not adversarial.

## Verdict justification

The script passes its `expect_zero` calls because six of ten assertions (A1, A2, A5, A6, A7, A10) are tautologies under the `const_subs` dictionary that hand-encodes the answers; the only non-trivial checks are the discrete Cauchy gap identity (A3/A8) and the free-energy softening algebra (A4/A9). The docstring's claimed main checks 1 and 2 (the closure law and the integral definitions of `I1`, `I2`) are never actually exercised — `chi_phi(y)` and `H(y)` are never instantiated. Add to this the fact that the Mathematica script mirrors the SymPy script structurally rather than re-deriving from first principles, and the cross-engine cross-check provides no additional safety. Verdict: **findings**; not stop-cold because the corrections are local (replace abstract-symbol assertions with concrete-profile integral derivations) and do not propagate to downstream units in a way that would invert a constant or sign.

## Self-test notes

- Variable independence: F3's proposed `sp.diff(local_F, sigma_sym)` — `local_F` depends on `sigma_sym` via the `(1/2) H sigma^2 - g_phi chi_phi sigma` form, so the derivative is `H*sigma - g_phi*chi_phi`, non-trivially zero at `sigma = g_phi chi_phi/H`. Confirmed.
- Symmetry/parity: F1's proposed Gaussian `chi_phi = exp(-y^2/(2L^2))` is even; `chi_phi^2/H` is even, integral over `(-∞, ∞)` is nonzero (`Hw>0`, gives `sqrt(pi)*L/Hw`). No parity-driven zero traps.
- Trivial-case pre-check: With `L=1, Hw=1`, `Npp_int = sqrt(pi)`, `I1_int = sqrt(pi)`, `I2_int = sqrt(pi)`, so `I1_int - Npp_int/Hw = 0` and `I2_int - Npp_int/Hw^2 = 0` — proposed assertions reduce correctly.
- Path specifications: F-block targets reference `scripts/...sympy_audit.py` and `mathematica/...mathematica_audit.wl`, consistent with the rest of the project layout.
