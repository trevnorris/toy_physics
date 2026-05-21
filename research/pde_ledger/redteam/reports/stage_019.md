---
unit_id: 019
batch: I.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 019 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The script audits the bridge from the parent wall block `(KSigma, MSigma)` into the isotropic full-bundle moment expansion `(P0, P2, P4)`. Concretely it asserts: (i) the one-pole defect `u4 - 4*u2^2` reduces to the compact numerator `D0*(B4+Z4) - 3*(MSigma+B2+Z2)^2` over `D0^2`; (ii) solving the one-pole defect for `KSigma` reproduces the closed form `B0 + Z0 + 3*(MSigma+B2+Z2)^2/(B4+Z4)` and matches the normalization branch `K_from_norm = B0 + Z0 + N0/P0_target` (with `P0_target = 54*G*cs^5/(5*a^5*c^5*mhat0^2)`); (iii) the `P2=P4=0` constant-prefactor branch closes algebraically with `N2_const, N4_const` whose Jacobian determinant equals `D0^3`, plus mutation guards that detect perturbed closures; (iv) the one-pole quadratic in `MSigma` factors as `-3*(MSigma - M_+)*(MSigma - M_-)` with `M_± = ±sqrt(D0*(B4+Z4)/3) - (B2+Z2)`, and the positive root carries `u2 = +gap/D0` while the negative root carries `u2 = -gap/D0`; (v) three numerical samples of the response signs and a concrete Gaussian profile `beta=exp(-w^2/2)` giving `MSigma=sqrt(pi)`, `KSigma=3*sqrt(pi)/2`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 39 | `assert_zero one_pole_defect - one_pole_numerator/D0^2` | yes |
| A2 | sympy | 51 | `assert_zero K_from_one_pole - (B0+Z0+3*(MSigma+B2+Z2)^2/(B4+Z4))` | yes |
| A3 | sympy | 52 | `assert_zero K_from_norm - (B0+Z0+N0/P0_target)` | yes |
| A4 | sympy | 56 | `assert_zero compatibility - (3*(MSigma+B2+Z2)^2/(B4+Z4) - N0/P0_target)` | yes |
| A5 | sympy | 85 | `assert_zero N2_const - N2_const_closed` | yes |
| A6 | sympy | 86 | `assert_zero N4_const - N4_const_closed` | yes |
| A7 | sympy | 87 | `assert_zero const_prefactor_matrix.det() - D0^3` | yes |
| A8 | sympy | 88 | `assert_zero P2 - (N2 - N2_const_closed)/D0` | yes |
| A9 | sympy | 89 | `assert_zero P4.subs(N2,N2_const_closed) - (N4 - N4_const_closed)/D0` | yes |
| A10 | sympy | 90 | `assert_zero const_prefactor_solution[N2] - N2_const_closed` | yes |
| A11 | sympy | 91 | `assert_zero const_prefactor_solution[N4] - N4_const_closed` | yes |
| A12 | sympy | 92-95 | `assert_nonzero mutated-N2 closure gives eps/D0` | yes |
| A13 | sympy | 96-99 | `assert_nonzero mutated-N4 closure gives eps/D0` | yes |
| A14 | sympy | 101 | `assert_zero N4_const(K_from_one_pole) - N4_md_one_pole(K_from_one_pole)` | yes |
| A15 | sympy | 107-110 | `assert_zero one_pole_numerator + 3*(MSigma-M+)*(MSigma-M-)` | yes |
| A16 | sympy | 111 | `assert_zero M+ + M- + 2*(B2+Z2)` (Vieta sum) | yes |
| A17 | sympy | 112-115 | `assert_zero M+ * M- - ((B2+Z2)^2 - gap^2)` (Vieta product) | yes |
| A18 | sympy | 118 | `assert_zero u2(M+) - gap/D0` | yes |
| A19 | sympy | 119 | `assert_zero u2(M-) + gap/D0` | yes |
| A20 | sympy | 120 | `assert_nonzero u2(M+) - u2(M-)` | yes |
| A21 | sympy | 121-122 | exactly two M roots | yes |
| A22 | sympy | 165-172 | three numeric stable-pole sign guards | yes |
| A23 | sympy | 184 | `assert_zero MSigma_example - sqrt(pi)` | yes |
| A24 | sympy | 185 | `assert_zero KSigma_example - 3*sqrt(pi)/2` | yes |

Every SymPy assertion is non-tautological and exercises the stated claim. None are trivially-true (no `x := expr; assert x == expr` patterns; assignments and assertions go through algebraic transformations).

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl` (does not exist)

**What's wrong:**
The unit ships a SymPy auditor but no Mathematica counterpart. Per project policy a non-status, non-checkpoint unit requires two independent algebraic engines. The single-engine result for the entire `(P0, P2, P4)` bridge — including the closed-form `K_from_one_pole`, `N2_const`, `N4_const`, the `det = D0^3` Jacobian, the M-root factorization, and the concrete Gaussian-profile integrals — has no second-engine confirmation.

**Why this matters:**
A single-engine algebraic identity can pass for the wrong reason (mis-applied `simplify`, an undeclared assumption silently picked up by SymPy's solver, or a one-engine convention on root ordering). Without a Mathematica re-derivation the closed forms `K_from_one_pole`, `N2_const_closed`, `N4_const_closed`, and the M-root closed forms are uncorroborated.

**Required change:**
Add `moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl` under `/var/projects/toy_physics/research/pde_ledger/mathematica/` that independently derives (not transliterates) the manifest items M1-M12 in the directive and aborts via `Print["FAIL: ..."]; Exit[1]` on any mismatch.

**Verification:**
After Codex writes the file, the verifier runs `redteam exec-mathematica 019`; the saved output must contain status "PASS", exit code 0, and explicit confirmations for M1-M12.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be checked. The directive demands the new Mathematica script derive each claim from the bundle definition (`P_k = expansion of N/(D0+D2*x+D4*x^2)` at small x) rather than copying the SymPy variable choreography (`D2/D0`, `D2^2 - D0*D4`, etc.). The directive lists what must be checked, not how.

## Engine cross-check

n/a — only one engine present.

## Verdict justification

The SymPy script holds up under attack: every assertion is anchored to a non-tautological algebraic identity. I attempted to break it by (a) checking whether `K_from_one_pole = solve(... , KSigma)[0]` could pick a wrong root — but the equation is linear in `KSigma` (a single root); (b) by checking whether the symbol assumption `nonzero=True` on `Z4` is violated by the `Z4=0` numeric sample — sympy treats `subs` as opaque to declared assumptions, and the downstream guard checks `B4+Z4>0` (which holds with `B4=1`), so the algebra of the assertions is undisturbed; (c) by checking whether `M_roots` ordering matters — the script's `M_root_positive/negative` are written by hand symbolically, not pulled by index from `M_roots`; (d) by checking the concrete Gaussian integrals by hand (`integral w^2 exp(-w^2) = sqrt(pi)/2`, `integral exp(-w^2) = sqrt(pi)`, total `3*sqrt(pi)/2` — matches A24); (e) by checking the mutation guards in A12/A13 algebraically (`P2 - (N2 - (N2_const + eps))/D0 = eps/D0`, nonzero when `eps` nonzero — sound); (f) by recomputing the Jacobian determinant from the expanded P2/P4 forms (`partial P2_zero_eq / partial N2 = D0`, `partial P2_zero_eq / partial N4 = 0`, `partial P4_zero_eq / partial N4 = D0^2`, det = `D0 * D0^2 = D0^3`). The single sustained finding is the absence of a second engine.

## Self-test notes

- Variable independence: every `sp.diff(P2_zero_eq, N2)` and `sp.diff(P4_zero_eq, N4)` is differentiated w.r.t. a symbol that actually appears in the expanded polynomial — confirmed by manual expansion of P4 (the N4 coefficient is `(B0-KSigma+Z0)^2 = D0^2`). No identically-zero derivative trap.
- Symmetry/parity of the Gaussian integrals: `(d beta/dw)^2 = w^2 exp(-w^2)` is even, integrand is even, integral is nonzero (`sqrt(pi)/2`). `beta^2 = exp(-w^2)` even, integral `sqrt(pi)`. Both nonzero as the script claims.
- Trivial-case substitution for the manifest items in the directive: with `B0=1, Z0=2, KSigma=13, B2=3, MSigma=1, Z2=4, B4=5, Z4=7` (script's baseline) the one-pole numerator becomes `10*12 - 3*(1+3+4)^2 = 120-192 = -72`; the defect `(120 - 192)/100 = -0.72`, matching `(u4 - 4 u2^2)` computed from `u2 = -8/10 = -0.8`, `u2^2 = 0.64`, `4*u2^2 = 2.56`, `u4 = (64 - 10*-12)/100 = 184/100 = 1.84`, `u4 - 4 u2^2 = -0.72`. M1 holds. All manifest items produce nonzero residuals on trivial substitutions.
- Path specification: the new `.wl` must live under `mathematica/` with the canonical name; the directive Target line names it explicitly.
