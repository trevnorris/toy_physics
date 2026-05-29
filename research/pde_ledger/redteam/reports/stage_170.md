---
unit_id: 170
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T16:30:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage170_linear_grouped_outlet_map.md]
  paper_appendix: present
---

# Audit unit 170 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_170.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 47, 71, 417-432, 554)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.txt`

## What the paper claims

The stage reduces the linear grouped-`P2` outlet problem to two exact microscopic combinations. The card `\stagefield{Output}` states verbatim: "Reduces the direct outlet deformations to \(\mathcal K_A=\delta D_{A,2}+\delta D_{A,0}/9\) and \(\mathcal G_A=\delta N_{A,0}-P_0\delta D_{A,0}\)." The notes elaborate the deliverables: (i) the first-order grouped transport laws for `δu_2`, `δu_4`, `δP_0` around the canonical compensated branch (`u_2=1/9`, `u_4=4/81`); (ii) the exact outlet map `δκ_W^(A)=3(1-σ*)/(σ*D_0)·K_A` and `δγ_W^(A)=-(1-σ*)/(9σ*N_0)·G_A`; (iii) the one-parameter even-consistency relation `δD_{A,4}=(2/3)δD_{A,2}+(1/27)δD_{A,0}`, equivalent to `δu_4=(8/9)δu_2`; (iv) the grouped trace/anomaly (a/b) projector forms; and (v) the weak-axisymmetric `Y_20` branch collapse to two scalar amplitudes `κ_1`, `γ_1` carrying the grouped signature `(λ_20,λ_21,λ_22)=(1,1/2,-1)`. Card Checks item 2 explicitly requires verifying that signature "before reducing grouped defects to a scalar."

## What the script claims to verify

The docstring lists checks 1-4 (transport laws, outlet map, even-consistency, trace/anomaly); Section 5 (added per the inline comment "paper Sec. 5 / card Checks item 2") claims to verify the weak-axisymmetric signature `(1,1/2,-1)` and the scalar-amplitude collapse to `κ_1`, `γ_1`. Sections 1-4 derive `δu_2`, `δu_4`, `δP_0` by Taylor-expanding the genuine definitions `u_2=-D_2/D_0`, `u_4=(D_2²-D_0D_4)/D_0²`, `P_0=N_0/D_0` and comparing the first-order coefficients to the paper's transport laws; they then solve the carried-forward hybrid relations for `δκ_W`/`δγ_W` and check against the paper's closed forms, and verify the `δu_4=(8/9)δu_2` even-consistency relation and the a/b projector forms. Section 5 feeds lane-scaled inputs into a re-typed copy of the outlet formula and compares to the same re-typed formula.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (i) transport laws `δu_2`, `δu_4`, `δP_0` | sympy L57-59 / wl L59-61, derived from definitions | match |
| (ii) outlet map `δκ_W`, `δγ_W` | sympy L73-80 / wl L77-84, solved from hybrid relations | match |
| (iii) even-consistency `δD_{A,4}=(2/3)δD_{A,2}+(1/27)δD_{A,0}` | sympy L89-91 / wl L90-95 | match |
| (iv) trace/anomaly a/b projectors | sympy L120-129 / wl L118-127 | match |
| (v) weak-axisymmetric signature + scalar collapse `κ_1`, `γ_1` | sympy L160-166 / wl L150-159 | mismatch (tautological — see F1) |

`paper_alignment: aligned` — every paper deliverable has a corresponding check and all formulas match the paper verbatim. The single finding is a verification-quality defect (the Section-5 check is tautological), not a formula/claim disagreement, so it is `tautological_check`, not `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 57 | `du2 + (dD2+dD0/9)/D0 == 0` | (i) δu_2 | yes |
| A2 | sympy | 58 | `du4 + (dD4+2dD2/9+5dD0/81)/D0 == 0` | (i) δu_4 | yes |
| A3 | sympy | 59 | `dP0 - (dN0-P0 dD0)/D0 == 0` | (i) δP_0 | yes |
| A4 | sympy | 73-76 | `dkappa_from_du2 - 3(1-σ)(dD2+dD0/9)/(σD0) == 0` | (ii) δκ_W | yes |
| A5 | sympy | 77-80 | `dgamma_from_dP0 + (1-σ)(dN0-P0 dD0)/(9σN0) == 0` | (ii) δγ_W | yes |
| A6 | sympy | 89 | `du4_from_kappa - (8/9) du2 == 0` | (iii) δE_4/δE_2=8/9 | yes |
| A7 | sympy | 91 | `relation - ((2/3)dD2 + dD0/27) == 0` | (iii) even-consistency | yes |
| A8 | sympy | 120-128 | a/b κ,γ `from_map - target == 0` | (iv) projectors | yes |
| A9 | sympy | 160-161 | `dkW[A] - eps lam kappa1 == 0` (etc.) | (v) scalar collapse | **no — tautological** |
| A10 | sympy | 163-166 | signature ratios `dkW[21]-(1/2)dkW[20]==0` (etc.) | (v) signature | **no — tautological** |
| B1-B8 | mathematica | 59-95, 118-127 | mirror of A1-A8 (independent route) | (i)-(iv) | yes |
| B9 | mathematica | 150-155 | `dkW20 - epsL kappa1 == 0` (etc.) | (v) scalar collapse | **no — tautological** |
| B10 | mathematica | 156-159 | signature ratios | (v) signature | **no — tautological** |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:144-166`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:138-159`

**What's wrong:**
Section 5 claims (inline comment, sympy L132-140 / wl L130-136) to verify card Checks item 2 — the weak-axisymmetric signature `(1,1/2,-1)` collapse to scalar amplitudes `κ_1`, `γ_1`. But the check is tautological: the amplitude target and the "mapped" output are built from the *same hand-retyped closed-form expression*, not from the separately-derived map object.

SymPy:
```python
def kappa_map(dD2_, dD0_):
    return 3*(1 - sigma)*(dD2_ + dD0_/9)/(sigma*D0)        # L144-145
kappa1 = 3*(1 - sigma)*(D2_1 + D0_1/9)/(sigma*D0)          # L150
dkW[A] = kappa_map(eps_l*lam*D2_1, eps_l*lam*D0_1)          # L158
expect_zero(f"...", dkW[A] - eps_l*lam*kappa1)             # L160
```
Because `kappa_map` and `kappa1` are the *identical* linear closed form, `kappa_map(eps_l·lam·D2_1, eps_l·lam·D0_1) = eps_l·lam·kappa1` holds *by construction* for any value of the physics — the assertion cannot fail no matter whether the outlet map derived in Section 2 is right or wrong. Likewise the signature checks at L163-166 (`dkW[21]-(1/2)dkW[20]`, `dkW[22]+dkW[20]`, gamma analogues) merely test that a linear expression is linear. The `.wl` repeats the same pattern verbatim: `kappaMap` (L138), `kappa1` (L140), `dkW20=kappaMap[epsL*1*D21,epsL*1*D01]` (L142), `expectZero[..., dkW20-epsL*kappa1]` (L150), and signature ratios L156-159.

The Section-2 maps `dkappa_from_du2` / `dgamma_from_dP0` (sympy L69-70) are *derived* objects (solved from `du2=du2_hyb` where `du2` itself came from series-expanding the definition `u_2=-D_2/D_0`). Section 5 ignores those derived objects and re-types a fresh copy of the answer, so it adds zero verification value beyond confirming the formula is linear.

**Why this matters:**
Card Checks item 2 is a named deliverable, and the saved transcripts print PASS for all six amplitude checks and four signature ratios (sympy output L44-53; mathematica output L55-74) — giving false confidence that the signature collapse was machine-verified. If a future edit broke the Section-2 derivation (e.g. a sign error in `du2_hyb`), Section 5 would still pass, because it never touches the derived map. A genuine check must flow the lane-scaled inputs through the derived map and compare to an independently written `κ_1`/`γ_1` target.

**Required change:**
In both scripts, route the lane-scaled inputs through the *derived* Section-2 map objects instead of the re-typed `kappa_map`/`gamma_map` (and `kappaMap`/`gammaMap`), keeping `kappa1`/`gamma1` as the independently-written paper target. Concretely, replace the `kappa_map`/`gamma_map` calls with substitution into `dkappa_from_du2` / `dgamma_from_dP0`:

SymPy (replace L144-148 helper defs and L158-159 calls):
```python
# delete kappa_map / gamma_map helper defs; use the derived Section-2 maps
for A, lam in lanes.items():
    dkW[A] = dkappa_from_du2.subs({dD2: eps_l*lam*D2_1, dD0: eps_l*lam*D0_1})
    dgW[A] = sp.simplify(
        dgamma_from_dP0.subs({dN0: eps_l*lam*N0_1, dD0: eps_l*lam*D0_1}).subs(P0, N0/D0)
    )
    expect_zero(f"delta kappa_W^({A}) - eps lambda kappa1", dkW[A] - eps_l*lam*kappa1)
    expect_zero(f"delta gamma_W^({A}) - eps lambda gamma1",
                sp.simplify(dgW[A] - eps_l*lam*gamma1.subs(P0, N0/D0)))
```
Mathematica (replace L138-139 helper defs and L142-147 assignments analogously): drop `kappaMap`/`gammaMap`, and set each `dkW2x` / `dgW2x` by substituting the lane-scaled inputs into the derived `dkappaFromdu2` / `dgammaFromdP0` (applying `P0 -> N0/D0` on the gamma side and on `gamma1`, mirroring the Section-2 gamma check at wl L83), then compare to `epsL*lam*kappa1` / `epsL*lam*gamma1`.

The signature-ratio checks (sympy L163-166 / wl L156-159) can stay as-is once `dkW`/`dgW` come from the derived map — they then test linearity of the *derived* map, which is acceptable as a secondary check, but the load-bearing amplitude checks above are now genuine.

**Verification:**
After the fix, the six amplitude assertions still print `= 0` / PASS (the derived map equals the paper formula, so they pass), but they now depend on `dkappa_from_du2` / `dgamma_from_dP0`. Confirm the new lines reference `dkappa_from_du2`/`dgamma_from_dP0` (sympy) and `dkappaFromdu2`/`dgammaFromdP0` (wl) rather than `kappa_map`/`gamma_map`, and that `redteam exec-sympy 170` / `exec-mathematica 170` exit 0 with all PASS lines intact.

## Independent-derivation check (Mathematica)

The `.wl` is not a pure transliteration of the `.py`. The load-bearing first-order extraction uses a different primitive: SymPy does `series(...).coeff(eps,1)` (sympy L52-54) whereas Mathematica uses the analytic derivative `D[u2Full, eps] /. eps -> 0` (wl L54-56), with the divergence explicitly noted in comments (wl L52-53, L66). The map inversion likewise drops SymPy's `du2sym`/`dP0sym` placeholder idiom in favor of native `dkappa /. First[Solve[...]]` (wl L66-74). The section structure is parallel and the final residuals match, but the core mechanism differs, so this does not rise to a `mathematica_transliteration` finding. Note that the Section-5 tautology (F1) is present identically in both engines, so the fix must be applied to both.

## Engine cross-check

Both engines emit identical PASS/`= 0` results for every check (sympy output L9-53 vs mathematica output L9-74). The non-trivial simplified forms also agree, e.g. `a_kappa = -(aD0+9*aD2)*(sigma-1)/(3*D0*sigma)` (sympy L28) vs `-1/3*((aD0+9*aD2)*(-1+sigma))/(D0*sigma)` (mathematica L35) — algebraically identical. No `engine_disagreement`.

## Verdict justification

Sections 1-4 hold up under attack: I tried to break them by checking whether `du2`/`du4`/`dP0` are derived from genuine definitions (they are — series-expansion of `u_2=-D_2/D_0` etc.), whether the outlet-map and even-consistency checks are tautological (they are not — they solve carried-forward relations against independent targets), and whether symbol assumptions are unsound (they are fine; `nonzero` on perturbations is harmless for these linear identities). All four sections faithfully exercise paper deliverables (i)-(iv) and the formulas match the card/notes/appendix verbatim. Section 5, however, is tautological in both engines: it feeds lane-scaled inputs through a re-typed copy of the answer formula and compares to the same formula, so it cannot fail regardless of the physics. The formula it claims is correct (matches notes Sec 5), so this is a verification-quality defect (`tautological_check`), not a `paper_misalignment`. Verdict: `findings` (one medium finding), no stop-cold — fixing it only strengthens an already-correct claim and propagates to no downstream constant.

## Self-test notes

Trap 1 (variable independence): the prescribed fix uses `.subs` / `/.` only — no `diff`/`D[]` on a non-dependent variable — so the identically-zero-derivative trap cannot fire; the existing Section-1 `series`/`D[u2Full,eps]` is on `eps`, which `u2Full` genuinely depends on. Trap 2 (parity): no unbounded integrals, N/A. Trap 3 (trivial-case): substituting `D2_1=1, D0_1=9` into the fixed check gives derived-map `eps·lam·6(1-σ)/(σD0)` equal to `eps·lam·kappa1`, residual 0 (passes), and a hypothetical sign bug in `dkappa_from_du2` would now make it nonzero (so the check is genuinely live); the `P0 -> N0/D0` reconciliation on the gamma side mirrors Section 2's wl L83 / sympy L79. Trap 5 (paper round-trip): the fix keeps `kappa1`/`gamma1` as the notes-Sec-5 boxed formulas and introduces no new constant, so no new `paper_misalignment`.
