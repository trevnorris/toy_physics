---
unit_id: 008
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 008 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The docstring says the script (1) writes the exact projected brane equation with a weighted gauge-fixing profile H(w), (2) derives the zero-mode projection-first effective coupling mu0_eff_proj = mu0 * I_WS/I_WZ and effective gauge parameter xi_eff_proj = xi * I_WZ/I_WH, and (3) compares the cases H=1 and H=Z, including Gaussian matched-kernel examples. The bottom-line claims the assertions are meant to test are: (a) for H=Z, xi_eff_proj = xi for ANY normalized projection kernel W; (b) for S = Z/Z_int, mu0_eff_proj = mu0/Z_int for any normalized W; (c) explicit Gaussian-kernel numerical instantiation; (d) reduction-first H=1 gauge parameter vanishes in the R→∞ regulator limit.

## Assertion inventory

| #  | Script | Line    | Form                                                                                                | Anchored to claim? |
|----|--------|---------|-----------------------------------------------------------------------------------------------------|--------------------|
| A1 | sympy  | 73      | `simplify(inv_xi_eff_proj * xi_eff_proj - 1) == 0`                                                  | no (tautology)     |
| A2 | sympy  | 125-128 | `simplify(I_WH_generic.subs(H,Z) - I_WZ_generic) == 0`                                              | no (substitution)  |
| A3 | sympy  | 132-135 | `simplify(∫W*H*B/xi).subs(H,Z) - B*∫W*Z/xi) == 0`                                                   | no (factoring)     |
| A4 | sympy  | 136-139 | `simplify((xi*I_WZ_generic/I_WH_generic).subs(H,Z)) - xi == 0`                                      | no (cancellation)  |
| A5 | sympy  | 149     | `simplify(mu0*(I_WZ/Zint)/I_WZ - mu0/Zint) == 0`                                                    | no (cancellation)  |
| A6 | sympy  | 191     | `Z_int_gauss - sqrt(pi)*lambda == 0`                                                                | yes                |
| A7 | sympy  | 192     | `Z2_int_gauss - sqrt(2*pi)*lambda/2 == 0`                                                           | yes                |
| A8 | sympy  | 193     | `I_WZ_match - sqrt(2)/2 == 0`                                                                       | yes                |
| A9 | sympy  | 194     | `xi_eff_HZ_match - xi == 0`                                                                         | partial            |
| A10| sympy  | 195     | `xi_eff_H1_match - xi/sqrt(2) == 0`                                                                 | yes                |
| A11| sympy  | 196     | `mu0_eff_source_match - mu0/Z_int_gauss == 0`                                                       | yes                |
| A12| sympy  | 197     | `mu0_eff_delta_match/(mu0/Z_int_gauss) - sqrt(2) == 0`                                              | yes                |
| A13| sympy  | 227     | `limit(sqrt(pi)*lambda*xi/(2*R), R, oo) == 0`                                                       | yes                |

A1–A5 are algebraic identities guaranteed by how `xi_eff_proj`, `inv_xi_eff_proj`, `I_WH`, `mu0_eff_match_source` were just defined. A9 is partial — for a normalized matched-Gaussian kernel it is a real numerical check, but the docstring's general "for ANY normalized W" claim is not exercised; the actual symbolic proof of that claim is A4, which is a tautology.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/` (no `*_stage008_*_mathematica_audit.wl`)

**What's wrong:**
The unit's manifest entry has `is_status_only_candidate: False` and the audit prompt confirms the mathematica script path is `(missing)`. There is no `.wl` script for stage 008. The instructions are explicit: a non-status-only, non-checkpoint unit requires both engines for independent verification. Currently the projected-Maxwell extension is verified by a single engine (SymPy) only.

Verified by directory listing: only `moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py` exists; no Mathematica counterpart of any name.

**Why this matters:**
The whole point of the dual-engine policy is to catch transcription errors and sympy-specific quirks (e.g., `simplify` succeeding under hidden assumptions, Integral-objects being collapsed by `simplify` heuristics) by re-deriving from physical premises in an independent CAS. With only one engine, every check below — including the substantive Gaussian-kernel evaluations on lines 191–197 — has no cross-engine confirmation.

**Required change:**
Add a Mathematica audit script `moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl` that independently verifies the unit's claims listed in the Claim manifest of the directive. It must use `FullSimplify`-backed `If[... Exit[1]]` style assertions (i.e., must actually exit non-zero on failure), not just `Print`.

**Verification:**
After Codex applies, the file exists at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl`, and `redteam exec-mathematica 008` exits 0.

---

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:73`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:125-128`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:136-139`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:149`

**What's wrong:**
Each of these assertions is algebraically guaranteed by the immediately preceding definition. They cannot fail for any choice of integrand profile or any choice of physics.

- Line 70–73: `mu0_eff_proj` and `xi_eff_proj` are defined as `mu0*I_WS/I_WZ` and `xi*I_WZ/I_WH`; `inv_xi_eff_proj = I_WH/(xi*I_WZ)`. The assertion `inv_xi_eff_proj * xi_eff_proj - 1 == 0` is `(I_WH/(xi*I_WZ))*(xi*I_WZ/I_WH) - 1`, which is the trivial reciprocal identity. The author wrote down two expressions that are reciprocals by construction; this assertion only confirms that fact.

- Lines 119–128: `I_WH_generic = Integral(W*H, w)`; the assertion `I_WH_generic.subs(H_generic, Z_generic) - I_WZ_generic == 0`. After substitution, the integrand becomes `W*Z`, identical to `I_WZ_generic`. This is sympy's substitution machinery, not a physical claim about the H=Z limit.

- Lines 136–139: `(xi * I_WZ_generic / I_WH_generic).subs(H_generic, Z_generic) - xi == 0`. After substituting H→Z, both integrals are identical, so the ratio is 1 and the residue is 0. This is xi/xi = xi, not a verification that the H=Z choice is the correct gauge-aligning weight.

- Line 148–149: `mu0_eff_match_source = simplify(mu0 * (I_WZ / Zint) / I_WZ)`, then `mu0_eff_match_source - mu0/Zint == 0`. The `I_WZ` cancels by elementary algebra. This asserts mu0/Zint = mu0/Zint.

**Why this matters:**
The comments above these assertions and the section headings claim they verify "H=Z effective gauge", "matched source coupling", "effective gauge inverse is reciprocal" — claims that are presented to the verifier (and downstream paper-reader) as substantive. They are not. The only substantive verifications of the H=Z and S=Z/Z_int claims are the concrete Gaussian numerics on lines 194–196. The symbolic "for any normalized W" assertion is decorative.

**Required change:**
Either (a) replace these tautological checks with substantive ones, or (b) remove them and rely on the Gaussian numerical checks. A substantive replacement for A1 (line 73) would be to start from the derived equation `I_WZ ∂_μ f^{μν} + (I_WH/ξ) ∂^ν B = μ0 I_WS j^ν` and verify that dividing through by I_WZ produces an equation whose gauge-driver coefficient is `1/(xi*I_WZ/I_WH)` — i.e., derive xi_eff_proj rather than write it down and then divide back out. A substantive replacement for A2–A4 is to define a concrete two-parameter family of profiles `W(w)` and `Z(w)` (e.g., two different Gaussians with different widths, or a Gaussian and an exponential) and check the H=Z claim numerically across that family, not by symbolic substitution into a placeholder Integral.

If the author insists on keeping the symbolic-substitution checks, mark them clearly in comments as "consistency tautologies, not verification" so future readers do not mistake them for evidence.

**Verification:**
After the rewrite, each of the four checks should: (1) start from a derivation, not from the answer; or (2) exercise a concrete second profile distinct from the Gaussian already used on lines 173–197. Output should mention at least one new profile pair name.

---

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:119-139`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:171-197`

**What's wrong:**
The script's stated claims are framed "for ANY normalized projection kernel W" (lines 115, 158, 75, 95-115 commentary, and the section-7 summary on lines 257-260). However, the actual exercises are:
- the abstract "ANY W" claims (lines 125–139) are tautological substitutions, as documented in F2;
- the only concrete W actually integrated against is the matched Gaussian kernel `W = Z/Z_int` with `Z = exp(-w²/λ²)` (lines 173–197) — a single, highly symmetric, self-matched choice.

This is a single-branch check standing in for a quantifier claim. The matched-Gaussian case has W∝Z, which makes I_WZ, I_WH (H=Z), and I_WS (S=Z/Z_int) all proportional to ∫Z²; the algebraic structure of "any normalized W" is not exercised — only the structure where W and Z share the same functional form is.

Furthermore, the normalization condition `∫W=1` is invoked in the prose (lines 51, 100, 158) but never enforced by an assertion or check; the Gaussian `W_match = Z/Z_int_gauss` does satisfy it (by construction of Z_int), but no profile with ∫W=1 and W not proportional to Z is tested.

**Why this matters:**
"For ANY normalized W" is a strong claim. It is the form of the result that downstream units will quote. Verifying it on a single self-aligned profile leaves open the possibility that an unnoticed sympy assumption or an algebraic error specific to W∝Z would pass these checks but fail on a generic W. The H=1 case on line 195 also uses the matched-Gaussian W, so the only H=1 numerical evidence is again specific to W∝Z.

**Required change:**
Add at least one additional concrete profile pair where W is NOT proportional to Z. Two practical choices:

1. `W(w) = (1/(sqrt(pi)*sigma)) * exp(-w²/sigma²)` with `sigma != lambda` (independent Gaussian kernel). Verify ∫W=1, then compute I_WZ, I_WH (for H=Z and H=1), I_WS (for S=Z/Z_int and S=delta), and assert:
   - xi_eff_proj(H=Z) - xi == 0
   - mu0_eff_proj(S=Z/Z_int) - mu0/Z_int == 0
   - xi_eff_proj(H=1) is sigma- and lambda-dependent (not equal to xi)

2. Optionally: a Lorentzian W (e.g., `W = (sigma/pi)/(w²+sigma²)`) paired with the Gaussian Z, for a second non-matched check.

Place these as a new section between current sections 5 and 6 (around line 218). The existing section-5 Gaussian matched-kernel checks remain.

**Verification:**
After fix, the script defines at least one new profile-pair (e.g., `W_indep = exp(-w²/sigma²)/(sqrt(pi)*sigma)` with a fresh symbol `sigma` declared positive), computes the three integrals against it, and asserts `xi_eff(H=Z) == xi` and `mu0_eff(S=Z/Z_int) == mu0/Z_int` non-tautologically. Output transcript shows the new section with the new symbol named.

## Independent-derivation check (Mathematica)

No `.wl` script exists for this unit (see F1). Independent-derivation comparison is impossible.

## Engine cross-check

No mathematica engine present. Cross-check is `n/a`.

## Verdict justification

The script's Gaussian-matched-kernel section (lines 173–197) is substantive: it computes real integrals and checks concrete numbers (sqrt(pi)*lambda for Z_int, sqrt(2)/2 for I_WZ, sqrt(2) factor for the delta-source mismatch, xi/sqrt(2) for H=1, mu0/Z_int for the matched source). The regulator-limit check on line 227 is also substantive. The reciprocal-identity check at line 73, the H→Z substitution checks at lines 125–139, and the matched-source-coupling check at line 149 are tautological — they verify their own definitions. The general "for ANY normalized W" claim is only exercised on a single W∝Z profile, leaving the quantifier underdetermined. And there is no independent Mathematica verification at all, which is the dominant finding. Hence verdict is `findings`, with the missing-mathematica as the high-severity item and tautology/insufficiency as medium-severity items.

Attacks tried that did NOT produce findings: (1) symbol-domain check — `xi`, `lam`, `mu0`, `R`, `Z_int`, `I_WZ`, `I_WH`, `I_WS` all declared positive/nonzero, which is consistent with their integral-of-positive-function interpretation; no `simplify` step relies on hidden complex-vs-real assumptions; (2) sign-and-factor check on Gaussian integrals — `Z_int = sqrt(pi)*lam`, `Z2_int = sqrt(pi/2)*lam = sqrt(2*pi)*lam/2`, both correct; ratio `Z2_int/Z_int = sqrt(2)/2` matches `I_WZ_match`; (3) regulator-limit direction — `lim_{R→∞} sqrt(pi)*lambda*xi/(2*R) = 0` is correct, and the interpretation "unsafe noncompact gauge fixing" follows; (4) stale-output — script mtime May 4, output mtime May 11, output is fresher than script, not stale.
