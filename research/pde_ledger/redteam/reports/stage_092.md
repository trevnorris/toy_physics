---
unit_id: 092
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage092_dynamic_geometry_obstruction.md
  paper_appendix: present
---

# Audit unit 092 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_092.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage092_dynamic_geometry_obstruction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1218 references the stage via `\input{stages/stage_092}`; no extra appendix prose specific to this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.txt`

## What the paper claims

Stage 092 derives an exact obstruction formula for the grouped-`P2` plus geometry split when the geometry lane carries dynamic `O(omega^2)` and `O(omega^4)` moments. The card's load-bearing block-quote states: "If geometry carries dynamic even moments, the pole fraction becomes `(1+eps_4)/[4(1+eps_2)^2]`." (`stage_092.tex:16`). The branch identity `K_0 K_4 = 4 K_2^2` is the imported constraint (`stage_092.tex:9`). The notes elaborate four deliverables: (i) the exact closed-form for `K_(g,0)` on the branch (Section 2), (ii) the dimensionless rewriting `K_0 = 4 K_pole (1+eps_2)^2/(1+eps_4)` giving the `c_pole = (1+eps_4)/[4(1+eps_2)^2]` identity (Section 3), (iii) the static-limit collapse to `c_pole = 1/4`, `c_geom = 3/4` when `eps_2 = eps_4 = 0` (Section 3), and (iv) the small-contamination first-order expansion `c_pole = (1/4)[1 + eps_4 - 2 eps_2 + O(eps^2)]` (Section 4). The stage card's Checks list also mentions an `l=0 / l=2` orthogonality precondition; that is an imported hypothesis from the upstream minimal-isotropic module, not an in-stage algebraic check.

## What the script claims to verify

Both scripts construct `K(omega) = K_(g,0) + K_(g,2) omega^2 + K_(g,4) omega^4 + K_pole/(1-omega^2/Omega_Q^2)`, series-expand to order `omega^4`, read off `K_0`, `K_2`, `K_4`, form the branch obstruction `K_0 K_4 - 4 K_2^2`, solve `branch == 0` for `K_(g,0)`, assert that the static-geometry limit (`K_(g,2)=K_(g,4)=0`) gives `K_(g,0) = 3 K_pole` (the equivalent of `c_pole = 1/4`), build `c_pole`, substitute the dimensionless variables `eps_2 = Omega_Q^2 K_(g,2)/K_pole` and `eps_4 = Omega_Q^4 K_(g,4)/K_pole`, and assert the closed-form `c_pole - (1+eps_4)/(4(1+eps_2)^2) == 0`. They additionally compute and display `c_geom = 1 - c_pole`, the small-`(eps_2, eps_4)` expansion, and the residual difference from the predicted linear part — but the residual is only printed, not asserted as second-order.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Pole fraction formula `c_pole = (1+eps_4)/[4(1+eps_2)^2]` (stage_092.tex:16; notes §3) | `expect_zero("c_pole - (1+eps4)/(4(1+eps2)^2)", ...)` at sympy:67; `expectZero[...]` at .wl:59 | match |
| Branch identity `K_0 K_4 = 4 K_2^2` (stage_092.tex:9; notes §2) | `branch = K0*K4 - 4*K2^2`, solved for `K_(g,0)` (sympy:47–50; .wl:45–48) | match |
| Static-limit collapse `eps_2=eps_4=0 ⇒ c_pole = 1/4` (Checks; notes §3) | `expect_zero("static-geometry limit gives 3 K_pole", ...)` (sympy:53–56; .wl:51) — equivalent because `K_(g,0)=3 K_p ⇒ K_0=4 K_p ⇒ c_pole=1/4` | match (indirect) |
| Exact `K_(g,0)` closed form (notes §2) | `Kg0_sol` printed (sympy:50; .wl:48–50) — no explicit assertion against a hardcoded form, but the form is derived in-script from `branch == 0`, so the closed form is exhibited not assumed | match |
| Small-contamination first-order expansion `c_pole = (1/4)[1+eps_4-2 eps_2 + O(eps^2)]` (notes §4) | `small_series`, `linear_part`, `remainder` are computed and printed (sympy:72–78; .wl:64–68); no `expect_zero` on the order of the remainder | partial |
| `l=0`, `l=2` orthogonality (Checks bullet, stage_092.tex:23) | Not exercised — imported upstream hypothesis carried forward from the minimal isotropic module | not-a-check (upstream import) |

Set `paper_alignment: aligned` — the load-bearing identity is exact in both engines; the gap is a missing assertion on a side deliverable, not a mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53–56 | `expect_zero(Kg0_sol(eps=0) - 3*Kp)` | static-limit `c_pole=1/4` (equivalent via `K_(g,0)=3 K_p`) | yes |
| A2 | sympy | 67 | `expect_zero(cpole_dimless - (1+eps4)/(4(1+eps2)^2))` | main `c_pole` formula | yes |
| A3 | mathematica | 51 | `expectZero["static-geometry limit gives 3 K_pole", ...]` | static-limit `c_pole=1/4` | yes |
| A4 | mathematica | 59 | `expectZero["c_pole - (1+eps4)/(4(1+eps2)^2)", cPoleDimless - cPoleExpected]` | main `c_pole` formula | yes |

Both load-bearing assertions are non-tautological: `cpole_dimless` is built from `Kg0_sol` (derived by solving the branch identity), with the substitution `Kg2 → eps2*Kp/Omega_Q^2`, `Kg4 → eps4*Kp/Omega_Q^4` and then compared against an independently constructed candidate `(1+eps4)/(4*(1+eps2)^2)`. The static-limit assertion is non-tautological in the same way: `Kg0_sol` is the algebraic solution of `K_0 K_4 = 4 K_2^2`, and we then assert that substituting the static values yields `3 K_pole`.

The small-`(eps_2, eps_4)` expansion (sympy:72–78; .wl:64–68) is computed but only printed; the difference between `small_series` and the predicted linear part is computed and printed, but no assertion checks that `remainder` is genuinely second-order or higher. The displayed output happens to show `-eps_2 eps_4/2`, which is second-order, but a regression that broke this would only be visible to a human reading the transcript, not via exit code.

## Findings

### F1 — mathematica_transliteration

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl:1-73`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py:1-83`

**What's wrong:**
The `.wl` is structurally a line-by-line port of the `.py`. Same choreography: build `kCons`, expand `Series[..., {omega, 0, 4}]`, read off `k0`/`k2`/`k4`, build `branch = k0*k4 - 4*k2^2`, `Solve[branch == 0, kg0]`, assert static-limit gives `3 kPole`, build `cPole`, substitute `eps2 = omegaQ^2 kg2/kPole` / `eps4 = omegaQ^4 kg4/kPole`, assert against `(1 + eps4)/(4*(1+eps2)^2)`, then the same small-series tail computation. The print labels are pairwise identical: e.g. `"Branch obstruction = "` (sympy:48; .wl:46), `"K_g0 on branch = "` (sympy:51; .wl:50), `"static-geometry limit gives 3 K_pole"` (sympy:54; .wl:51), `"c_pole - (1+eps4)/(4(1+eps2)^2)"` (sympy:67; .wl:59), `"First-order expansion of c_pole = "` (sympy:76; .wl:66), `"Linear part = "`, `"Dropped higher-order tail = "`. Variable names are the only place they diverge (`Kg0` vs `kg0`, `OmegaQ` vs `omegaQ`).

Quoted pair (algorithmic core, sympy:35–50 vs .wl:33–48):

```
# sympy
K = sp.simplify(Kg0 + Kg2 * omega**2 + Kg4 * omega**4 + Kp / (1 - omega**2 / OmegaQ**2))
series = sp.expand(sp.series(K, omega, 0, 6).removeO())
K0 = sp.simplify(series.coeff(omega, 0))
K2 = sp.simplify(series.coeff(omega, 2))
K4 = sp.simplify(series.coeff(omega, 4))
...
branch = sp.simplify(K0 * K4 - 4 * K2**2)
Kg0_sol = sp.solve(sp.Eq(branch, 0), Kg0)[0]
```

```
(* mathematica *)
kCons = FullSimplify[kg0 + kg2*omega^2 + kg4*omega^4 + kPole/(1 - omega^2/omegaQ^2), ...];
series = Expand[Normal[Series[kCons, {omega, 0, 4}]]];
k0 = FullSimplify[Coefficient[series, omega, 0], ...];
k2 = FullSimplify[Coefficient[series, omega, 2], ...];
k4 = FullSimplify[Coefficient[series, omega, 4], ...];
...
branch = FullSimplify[k0*k4 - 4*k2^2, ...];
kg0Sol = FullSimplify[First[Solve[branch == 0, kg0]], ...];
```

**Why this matters:**
The second-engine policy requires the Mathematica script to derive the result independently. Here both engines perform the same series→coefficient→solve→substitute pipeline with identical intermediate variables, so a logic bug in the chosen derivation path would be reproduced in both engines and the cross-check would silently pass. This is a low-severity finding because the identity is algebraically simple (one rational equation in five symbols), so an independent re-derivation would likely arrive at the same closed form, but the Mathematica script should at minimum start from a different algebraic vantage point — e.g., directly substituting the dimensionless `eps_2`, `eps_4` into `K_0`, `K_2`, `K_4` and verifying the branch identity in those variables, or computing `c_pole = K_pole/K_0` after using the branch-solved `K_0 = 4 K_pole (1+eps_2)^2/(1+eps_4)` from notes Section 3 directly without going through `Solve[branch == 0, kg0]`.

**Required change:**
Restructure the Mathematica derivation so it does not mirror the sympy step-by-step. A minimal independent path:
1. Define `eps2 = omegaQ^2 * kg2/kPole`, `eps4 = omegaQ^4 * kg4/kPole` as inputs.
2. Compute `k0 = kPole*(1+eps_2)` style expressions in the dimensionless variables directly.
3. Form `c_pole` analytically and verify the closed form `(1+eps4)/(4(1+eps2)^2)` follows from the branch identity in eps-form.
4. Independently verify the static limit by substitution on the dimensionless expression rather than through the `kg0`-elimination route.

See the directive for the suggested replacement (`mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl:33–68`).

**Verification:**
After Codex applies, the rewritten `.wl` should: (a) not contain `Solve[branch == 0, kg0]` (the new path eliminates this intermediate); (b) still print `c_pole - (1+eps4)/(4(1+eps2)^2) = 0` and pass; (c) still cover the static limit `eps_2 = eps_4 = 0 ⇒ c_pole = 1/4`; (d) exit 0.

### F2 — insufficient_verification

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py:72-78`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl:64-68`

**What's wrong:**
The notes Section 4 explicitly lists the first-order expansion `c_pole = (1/4)[1 + eps_4 - 2 eps_2 + O(eps^2)]` as a deliverable. The scripts compute `small_series`, predict `linear_part = (1/4)(1 + eps_4 - 2 eps_2)`, and compute `remainder = small_series - linear_part`, but only PRINT the residual:

```
# sympy:78
print("Dropped higher-order tail  =", remainder)
```

```
(* .wl:68 *)
Print["Dropped higher-order tail = ", fmt[Expand[smallSeries - linearPart]]];
```

There is no `expect_zero` / `expectZero` (or any check that `remainder` lies in the ideal generated by `eps_2^2, eps_2*eps_4, eps_4^2`). A regression that, say, changed `linear_part` to `(1/4)(1 + eps_4 - eps_2)` would still execute and print, and would only be caught by a human reading the transcript.

**Why this matters:**
The first-order expansion is a stated deliverable. The script claims to verify the small-contamination behavior but no assertion enforces it. Either the assertion should be added, or the prints should be honest about being illustrative-only — but since the notes call out this expansion as a substantive part of the result, the assertion is appropriate.

**Required change:**
Add an assertion in both scripts that the lowest-order content of `small_series` agrees with `linear_part`. The simplest formulation: assert that the coefficient of `eps_2^0 eps_4^0`, the coefficient of `eps_2`, and the coefficient of `eps_4` in `small_series` match `1/4`, `-1/2`, and `1/4` respectively. Concretely in SymPy (insert after line 78):

```python
expect_zero("first-order eps^0 coefficient", small_series.coeff(eps2, 0).coeff(eps4, 0) - sp.Rational(1, 4))
expect_zero("first-order eps2 coefficient", small_series.coeff(eps2, 1).coeff(eps4, 0) - sp.Rational(-1, 2))
expect_zero("first-order eps4 coefficient", small_series.coeff(eps4, 1).coeff(eps2, 0) - sp.Rational(1, 4))
```

And mirror those checks in the Mathematica file via `Coefficient[smallSeries, eps2^i eps4^j]` or equivalent.

**Verification:**
After Codex applies, the sympy transcript should print three new `... = 0` lines and pass; the Mathematica transcript should print three new `PASS:` lines. If any first-order coefficient mismatches, the script exits 1.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration; see F1. It uses the same `kg0 + kg2 omega^2 + ...` template, the same `Series[..., {omega, 0, 4}]` expansion, the same `Solve[branch == 0, kg0]`, and the same `cPoleDimless - cPoleExpected` residual. The only distinguishing feature is the use of `FullSimplify` with explicit `$Assumptions`, but the algebraic path is identical.

## Engine cross-check

Both engines independently report:

- Static-limit residual: `0` (sympy:20 of output; .wl output:20)
- `c_pole - (1+eps4)/(4(1+eps2)^2)` residual: `0` (sympy:23; .wl:24)
- Dropped higher-order tail: `-eps_2*eps_4/2` in both (sympy:27; .wl:29) — they agree on the residual content too.

Engines agree.

## Verdict justification

The mathematical claim — `c_pole = (1+eps_4)/[4(1+eps_2)^2]` derived from `K_0 K_4 = 4 K_2^2` plus the conservative-carrier pole — holds up. Both engines verify it via non-tautological assertions, and the static-geometry limit collapsing to `c_pole = 1/4` is asserted in both. The notes-side derivation chain (Section 2 closed form for `K_(g,0)`, Section 3 dimensionless rewrite, Section 4 first-order expansion) is faithfully reproduced; the first three deliverables have load-bearing checks, and the fourth is computed but not asserted. The Mathematica file is a structural transliteration of the SymPy file. Verdict: `findings`, but only low-severity — the load-bearing identity is solid, and both findings are about hardening the verification surface rather than disputing the math.

## Self-test notes

- Checked variable independence: no spurious `sp.diff` chains — the script only does series expansion and `Solve`, both used correctly.
- Checked symbol assumptions: `K_p > 0`, `Omega_Q > 0` are physically motivated and used consistently in both engines.
- Checked the small-series remainder: substituting `eps_2 = eps_4 = 0` into `small_series` gives `1/4` (matches `linear_part`); the proposed coefficient checks pass mentally; substituting concrete `eps_2 = 0.01, eps_4 = 0` into `(1+eps_4)/[4(1+eps_2)^2]` and expanding gives `1/4 - eps_2/2 + ...`, consistent.
- Checked engine-disagreement: residuals match exactly across engines.
