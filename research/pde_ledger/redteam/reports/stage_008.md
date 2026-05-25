---
unit_id: 008
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 008 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_008.tex`
- notes: (none — per prompt; em_projected per-stage notes not committed)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.txt`

## What the paper claims

Stage 008 ("Projected Maxwell extension and clean matching") exports two matching conditions
between the projection-first weighted Maxwell extension and the reduction-first localized
gauge-fixed Maxwell form. Quoting verbatim, the card states `H(w)=Z(w), \xi_{\rm eff}^{\rm proj}=\xi`
(eq:stage005-hz) and `S(w)=Z(w)/Z_{\rm int}, \mu_{\rm eff}^{\rm proj}=\mu_0/Z_{\rm int}`
(eq:stage005-source-match), and the `\paragraph{Output.}` line says "Stage~008 exports the
matching conditions \eqref{eq:stage005-hz}--\eqref{eq:stage005-source-match}. These conditions
are the firewall between exact projection and reduced brane coupling." The card also says
"The audit also checks the Gaussian case and the regulator limit, including mutation guards
distinguishing H=1 from H=Z." Part-appendix row 008 summarizes the deliverable as "Slot-level
distinction between source-flux projection and matched reduction."

## What the script claims to verify

The SymPy script asserts (1) the zero-mode projected effective parameters `mu0_eff_proj = mu0 I_WS/I_WZ`
and `xi_eff_proj = xi I_WZ/I_WH`, (2) that H=Z preserves xi exactly (both as a symbolic substitution
identity and concretely on a matched-Gaussian kernel and on an independent-Gaussian kernel with
sigma != lambda), (3) that S=Z/Z_int reproduces `mu0/Z_int` on those same two profiles, (4) a
mutation-guard contrast that H=1 yields `xi/sqrt(2)` (matched Gaussian) or `xi*lambda/sqrt(lambda^2+sigma^2)`
(independent Gaussian) -- i.e., not xi -- and (5) that the reduction-first H=1 case has zero
regulator limit. The Mathematica script asserts the same family of identities using a different
second-profile choice (Lorentzian observer kernel) plus concrete Gaussian numerics and a numerical
H=Z/matched-source check at two distinct (lambda,sigma) parameter pairs.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `H=Z => xi_eff_proj = xi` (eq:stage005-hz), for any normalized W | sympy s5 Gaussian, s5b independent Gaussian; .wl M2 Pair B Lorentzian (independent profile), M4 matched Gaussian | match (substantive structurally, but F1 flags that the concrete-profile H=Z ratios collapse to X/X by construction in both engines) |
| `S=Z/Z_int => mu_eff_proj = mu0/Z_int` (eq:stage005-source-match), for any normalized W | sympy s4 + s5 matched + s5b independent; .wl M3, M4 distributed, M7 Lorentzian numeric | match |
| "Gaussian case" | sympy s5; .wl M4 | match |
| "regulator limit" | sympy s6 (xi*sqrt(pi)*lambda/(2R) -> 0 as R -> inf); .wl M5 | match |
| "mutation guards distinguishing H=1 from H=Z" | sympy s5 (xi/sqrt(2) != xi) + s5b (xi*lambda/sqrt(lambda^2+sigma^2) != xi at sigma=1/2,lambda=1); .wl M4 H=1 (xi/sqrt(2)) | match |
| Slot-level distinction source-flux vs matched reduction (appendix row) | sympy s5 case A vs case B (sqrt(2) ratio), s7 summary point 4; .wl M4 delta-source ratio sqrt(2) | match |

`paper_alignment: aligned`. Every paper-side deliverable has a script-side check that exercises it,
on at least one profile per engine. The defect is local: two of the four "H=Z" checks on the
Mathematica side (M2 Pair A, M2 Pair B) and one inside M4 (line 119) are computed as
`Integrate[X]/Integrate[X]` where the two `Integrate` calls have identical integrands; those
particular ratios are 1 regardless of the underlying physics. The substantive matched-source
and H=1 checks on the Mathematica side, and the SymPy mutation-guard / matched-source checks,
do exercise the paper claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74 | `inv_xi_eff_proj * xi_eff_proj - 1 == 0` | none (reciprocal scaffolding) | no (labeled "tautology by construction") |
| A2 | sympy | 130-133 | `I_WH.subs(H,Z) - I_WZ == 0` | H=Z claim (symbolic substitution) | partial (labeled "tautology by substitution") |
| A3 | sympy | 137-140 | gauge-driver factoring with H->Z | H=Z claim | partial (labeled "tautology by substitution") |
| A4 | sympy | 141-144 | `xi*I_WZ/I_WH .subs(H,Z) - xi == 0` | H=Z claim (substitution) | partial (labeled tautology) |
| A5 | sympy | 156 | matched-source coupling cancellation | S=Z/Z_int claim (substitution) | partial (labeled tautology) |
| A6 | sympy | 198 | `Z_int_gauss - sqrt(pi)*lambda == 0` | Gaussian-case concrete value | yes |
| A7 | sympy | 199 | `Z2_int_gauss - sqrt(2pi)*lambda/2 == 0` | Gaussian-case concrete value | yes |
| A8 | sympy | 200 | `I_WZ_match - sqrt(2)/2 == 0` | Gaussian-case concrete value | yes |
| A9 | sympy | 201 | `xi_eff_HZ_match - xi == 0` (W=Z/Z_int, H=Z) | H=Z claim, matched Gaussian | partial -- I_WH_HZ_match and I_WZ_match are integrated by literally identical expressions, so this ratio is 1 by construction; see F1 |
| A10 | sympy | 202 | `xi_eff_H1_match - xi/sqrt(2) == 0` | H=1 mutation guard (Gaussian) | yes |
| A11 | sympy | 203 | `mu0_eff_source_match - mu0/Z_int_gauss == 0` | S=Z/Z_int claim, matched Gaussian | yes |
| A12 | sympy | 204 | `mu0_eff_delta_match / (mu0/Z_int_gauss) - sqrt(2) == 0` | source-flux vs matched-reduction slot distinction | yes |
| A13 | sympy | 235 | independent W normalization | scaffolding for s5b | yes |
| A14 | sympy | 241 | `xi_eff_HZ_indep - xi == 0` (W=Gaussian sigma, H=Z) | H=Z claim, independent profile | partial -- I_WZ_indep and I_WH_HZ_indep also computed by identical Integrate calls; ratio is 1 by construction; see F1 |
| A15 | sympy | 250-253 | `xi_eff_H1_indep - xi != 0` at concrete (sigma=1/2,lambda=1) | H=1 mutation guard, independent profile | yes |
| A16 | sympy | 259-262 | `mu0_eff_source_indep - mu0/Z_int_indep == 0` | S=Z/Z_int claim, independent profile | yes |
| A17 | sympy | 282 | reduction-first H=1 regulator limit = 0 | regulator-limit deliverable | yes |
| M1 | wl | 17-22 | reciprocal residual | none (scaffolding) | no (algebraic tautology) |
| M2a-norm | wl | 37,41 | `gaussNorm - 1 == 0` | scaffolding | yes |
| M2a-HZ | wl | 38,44 | `xi*gaussOverlap/gaussGaugeWeight - xi` (Gaussian W, Gaussian Z) | H=Z claim, Gaussian-Gaussian | no -- `gaussOverlap` and `gaussGaugeWeight` are computed by literally identical `Integrate` calls (see F1); ratio is X/X = 1 |
| M2b-norm | wl | 58,62 | `lorentzNorm - 1 == 0` | scaffolding | yes |
| M2b-HZ | wl | 59,65 | `xi*lorentzOverlap/lorentzGaugeWeight - xi` (Lorentzian W) | H=Z claim, Lorentzian profile | no -- `lorentzOverlap` and `lorentzGaugeWeight` are also identical `Integrate` calls (see F1) |
| M3 | wl | 76,79 | `mu0*gaussMatchedSourceOverlap/gaussOverlap - mu0/zArea` | S=Z/Z_int claim, Gaussian profile | yes |
| M4a | wl | 100,110 | `zArea - sqrt(pi)*lambda` | Gaussian-case value | yes |
| M4b | wl | 100,113 | `zSelfArea - sqrt(pi/2)*lambda` | Gaussian-case value | yes |
| M4c | wl | 100,116 | `matchedOverlap - sqrt(2)/2` | Gaussian-case value | yes |
| M4d | wl | 100,119 | `xi*matchedOverlap/matchedOverlap - xi` | H=Z (matched Gaussian) | no -- `X/X = 1` algebraic tautology with no Integrate diversity (see F1) |
| M4e | wl | 100,122 | `xi*matchedOverlap/matchedNorm - xi/sqrt(2)` | H=1 mutation guard | yes |
| M4f | wl | 100,125 | `(mu0*matchedDeltaSource/matchedOverlap)/(mu0/zArea) - sqrt(2)` | slot-distinction (delta vs matched) | yes |
| M4g | wl | 100,128 | `mu0*matchedDistributedSource/matchedOverlap - mu0/zArea` | S=Z/Z_int (matched Gaussian) | yes |
| M5 | wl | 134,138 | `Limit[xi4UnweightedReg[R], R->Infinity] - 0` | regulator-limit deliverable | yes |
| M6 | wl | 145,147 | `xi*zIntSym/zIntSym - xi` | H=Z (algebraic, reduction-first) | partial (substitution identity, by section a reduction-first algebra check) |
| M7-gauge | wl | 155-159,165 | numeric `xi*lorentzOverlap/lorentzGaugeWeight - xi` at two (lambda,sigma) | H=Z, Lorentzian, numeric | no -- same X/X structure as M2b-HZ, numerically (see F1) |
| M7-src | wl | 160-162,168 | numeric `mu0*lorentzMatchedSourceOverlap/lorentzOverlap - mu0/zArea` | S=Z/Z_int, Lorentzian, numeric | yes |

## Findings

### F1 -- tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl:31-46` (M2 Pair A)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl:49-67` (M2 Pair B)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl:88-124` (M4d, line 119)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl:152-167` (M7 gauge)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:184-201` (matched Gaussian H=Z, lines 184/186/201)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py:237-241` (independent W H=Z, lines 237/239/241)

**What's wrong:**

The paper claim under audit is `H = Z => xi_eff_proj = xi for any normalized W` (eq:stage005-hz).
The general derivation `xi_eff_proj = xi * I_WZ / I_WH` is fine; the script then attempts to
verify the H=Z specialization on concrete profiles. The defect is that on every concrete
profile, the script computes both I_WZ and I_WH_HZ by calling `Integrate[W*Z, ...]` twice
with literally identical arguments and storing the result in two distinct symbols, then
checks `xi*I_WZ/I_WH_HZ - xi == 0`. Because the two `Integrate` calls produce the same closed
form (they are the same integral), the ratio is 1 by simple algebra of identical symbolic
values, regardless of whether `Integrate` even returned the correct integral.

Quoting Mathematica M2 Pair A, lines 31-38:
```
gaussOverlap =
  Integrate[gaussianWeight*localizedGaussian, {w, -Infinity, Infinity}, ...];
gaussGaugeWeight =
  Integrate[gaussianWeight*localizedGaussian, {w, -Infinity, Infinity}, ...];
m2aGaugeResidual = FullSimplify[xi*gaussOverlap/gaussGaugeWeight - xi];
```
The `gaussOverlap` and `gaussGaugeWeight` integrals are textually identical. The H=Z check
`xi*X/X - xi` evaluates to 0 by algebra alone; it does not test that integration of W*Z equals
integration of W*H with H=Z, because the script never assembles a separate `H` and never
integrates `W*H` against a distinct (even if equal) expression. This pattern repeats at
M2 Pair B (lines 52-59), M4 line 119 (`xi*matchedOverlap/matchedOverlap - xi`), and M7 gauge
residual (lines 155-159, `xi*lorentzOverlap/lorentzGaugeWeight - xi` at two numeric parameter
pairs).

The SymPy script has the same structural defect on its two H=Z gauge-parameter checks:
- Line 184: `I_WZ_match = simplify(integrate(W_match * Z, ...))`
- Line 186: `I_WH_HZ_match = simplify(integrate(W_match * Z, ...))`  (the two integrands are
  identical strings)
- Line 201: `assert_zero("Gaussian H=Z gauge parameter", xi_eff_HZ_match - xi)` where
  `xi_eff_HZ_match = simplify(xi * I_WZ_match / I_WH_HZ_match)` collapses to xi by `X/X`.

Same in s5b:
- Line 237: `I_WZ_indep = simplify(integrate(W_indep * Z, ...))`
- Line 239: `I_WH_HZ_indep = simplify(integrate(W_indep * Z, ...))`  (identical)
- Line 241: `assert_zero("independent-profile H=Z gauge parameter", xi_eff_HZ_indep - xi)`
  collapses to `X/X` algebra.

Note that the SymPy script's s3 already acknowledges the symbolic-level H=Z check is "tautology
by substitution"; the s5 / s5b concrete-Gaussian checks were intended to be the substantive
verification but inherit the same structural defect because the `H` weight is never carried as
a separate, computed-against expression.

**Why this matters:**

The paper card explicitly says "mutation guards distinguishing H=1 from H=Z". The H=1 mutation
guards ARE genuine: line 202 (sympy) and line 122 (.wl) compute `xi/sqrt(2)` and assert it
differs from xi. The H=Z side of the mutation guard, however, is currently inert on every
profile -- if the matched-Gaussian or Lorentzian `Integrate` returned the wrong value (or even
the literal symbol `Indeterminate`), `X/X` would still equal 1 and the check would still pass.
That means the only thing the script actually verifies for H=Z on concrete profiles is that
SymPy/Mathematica can perform the algebraic cancellation of identical expressions. The exact
physical claim -- that integration of `W*H` with `H=Z` yields the same value as integration of
`W*Z` -- is asserted by construction rather than computed.

This is the v2 question-1 risk inverted: the script's assertion does not exercise the paper's
claim at the concrete-profile level. The claim is true by direct substitution, but the
"verification on concrete kernels" the paper card promises is currently absent. A second-engine
audit is supposed to catch precisely this kind of trivial circle.

**Required change:**

Make the H=Z integral textually distinct from the I_WZ integral on each concrete profile, so
that the assertion forces the integrator to independently produce the same closed form. The
cleanest pattern is: keep the original I_WZ computation, add a parallel I_WH_HZ computation
whose integrand uses an alias for Z (or a re-expressed copy), then add a substantive equality
assertion `I_WH_HZ - I_WZ == 0` BEFORE the ratio assertion.

For the SymPy script:
- After line 180 (`Z = sp.exp(-w**2 / lam**2)`), introduce a separate `H_Z_expr = sp.exp(-w**2 / lam**2)`
  (textually re-typed, not the Python `Z` object) and compute `I_WH_HZ_match = simplify(integrate(W_match * H_Z_expr, ...))`
  to replace line 186. Add an assertion just before line 201:
  `assert_zero("Gaussian H=Z integrals match (independent computation)", I_WH_HZ_match - I_WZ_match)`.
  Same for `I_WH_HZ_indep` at line 239.

For the Mathematica script:
- M2 Pair A (lines 31-46): replace lines 34-36 with
  `gaussGaugeWeight = Integrate[gaussianWeight*Exp[-w^2/lambda^2], {w, -Infinity, Infinity}, Assumptions -> lambda > 0 && sigma > 0];`
  and add a check `If[FullSimplify[gaussGaugeWeight - gaussOverlap] =!= 0, Print["FAIL: M2 Pair A H=Z integrals match"]; Exit[1]];`
  before the existing H=Z gauge check.
- M2 Pair B (lines 49-67): same treatment with `lorentzGaugeWeight = Integrate[lorentzWeight*Exp[-w^2/lambda^2], ...]`.
- M4 line 119: introduce `matchedGaugeOverlap = Integrate[matchedWeight*Exp[-w^2/lambda^2], {w, -Infinity, Infinity}, Assumptions -> lambda > 0];`
  computed earlier in the M4 block, then replace `xi*matchedOverlap/matchedOverlap - xi` with
  `xi*matchedOverlap/matchedGaugeOverlap - xi` (and add an equality check
  `matchedGaugeOverlap - matchedOverlap == 0` to the residual list before line 119).
- M7 (lines 152-167): recompute a fresh `lorentzGaugeOverlap = Integrate[lorentzWeight*Exp[-w^2/lambda^2], ...]`
  and use it in the M7 numeric gauge residual instead of `lorentzGaugeWeight`.

After the changes, each concrete-profile assertion order is: (a) the H integral and the Z
integral agree as closed forms, (b) the gauge ratio is xi. Both must pass.

**Verification:**

After Codex applies, the verifier runs both engines fresh. Each affected H=Z assertion should
still pass (because the integrals really are equal). The check now exercises the integrator on
two textually-distinct integrands and confirms they collapse, which is the substantive content
of the paper claim. New "M2 Pair A H=Z integrals match" print line should appear at the
re-computed location; similarly for Pair B, M4 H=Z, M7. SymPy output should gain analogous
"Gaussian H=Z integrals match" and "independent-profile H=Z integrals match" lines.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT a transliteration of the SymPy script. SymPy uses two profiles for
its concrete checks: (a) a matched Gaussian W=Z/Z_int, (b) an independent Gaussian W with width
sigma != lambda. Mathematica uses (a) a generic Gaussian W with width sigma, (b) a Lorentzian W
`(sigma/Pi)/(w^2+sigma^2)`. The Lorentzian profile is a structurally different choice of
observer kernel that the SymPy script does not test. The numerical M7 cross-check at two
parameter pairs (`{lambda->1, sigma->1/2}` and `{lambda->1, sigma->2}`) is also not present in the
SymPy script. The Mathematica script independently derives the matched-source coupling using
its own profile machinery rather than echoing the SymPy intermediate `mu0_eff_source_match`
symbol-by-symbol. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines exit cleanly with `STATUS: PASS`. Where they verify overlapping claims (Gaussian
matched source coupling, regulator limit, mutation-guard H=1 = xi/sqrt(2), Gaussian Z_int and
Z2_int values), the numerical/symbolic outputs agree:
- `Z_int = sqrt(pi)*lambda` (sympy line 198, .wl M4 line 110)
- `Z2_int = sqrt(pi/2)*lambda` (sympy line 199 (sqrt(2pi)*lambda/2 ≡ sqrt(pi/2)*lambda), .wl line 113)
- `I_WZ_match = sqrt(2)/2` (sympy line 200, .wl line 116)
- H=1 matched: `xi/sqrt(2)` (sympy line 202, .wl line 122)
- matched source `mu0/Z_int` (sympy line 203, .wl line 128 and M3 line 79)
- delta-source ratio `sqrt(2)` (sympy line 204, .wl line 125)
- regulator limit 0 (sympy line 282, .wl line 138)

No `engine_disagreement` finding.

## Verdict justification

The paper claims are fully covered by the assertion set in both engines: the matching conditions
H=Z and S=Z/Z_int are present, the Gaussian explicit case is computed, the regulator limit is
derived symbolically and confirmed via Limit, the mutation-guard contrast against H=1 produces
a quantitatively distinct result (xi/sqrt(2) for matched Gaussian, lambda*xi/sqrt(lambda^2+sigma^2)
for independent Gaussian), and the source-flux vs matched-reduction slot distinction (the sqrt(2)
delta-vs-matched ratio) is checked on both engines. paper_alignment is `aligned`.

The single finding is local: the H=Z gauge-parameter check on every concrete profile (matched
Gaussian and independent Gaussian on the SymPy side; Gaussian-Gaussian, Lorentzian-Gaussian,
M4 matched, M7 Lorentzian numeric on the Mathematica side) is currently structured as `X/X - 1`
because both I_WZ and I_WH (with H=Z) are computed from textually identical Integrate calls.
The check passes by algebraic cancellation regardless of what the integrator returns. The fix
is mechanical and does not alter any paper claim -- it just makes the H=Z verification on
concrete profiles substantive rather than algebraically forced. The mutation-guard H=1 side
is genuine, so the symmetry the paper card promises is half-real; the fix completes it.

Verdict: `findings`, `stop_cold: null`. One script-side finding to be applied by Codex.

## Self-test notes

I checked the trap categories: (1) `sp.diff` variable-independence -- no derivatives in the
proposed fix, only integrations of textually-distinct integrands. (2) parity/symmetry -- the
proposed re-computed integrals are even-times-even Gaussian/Lorentzian × Gaussian, the integrals
are nonzero finite, no spurious vanishing. (3) trivial-case pre-check -- for the matched
Gaussian H=Z fix, `Integrate[(Z/Z_int)*Exp[-w^2/lambda^2], ...] = Integrate[(Z/Z_int)*Z, ...]`
must give the same `sqrt(2)/2`; both forms reduce to the same closed form via standard
Gaussian-product integration, so the new assertion passes substantively. (4) path specifications
-- both scripts already exist; no `missing_verification_script` paths to construct. (5) paper
round-trip -- the fix does not introduce any new constant or change any existing one; it only
restructures how H enters the H=Z integral, which the paper specifies as `H(w) = Z(w)` directly
(eq:stage005-hz). No new paper_misalignment is created.
