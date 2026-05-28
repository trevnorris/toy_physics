---
unit_id: 172
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage172_physical_slope_collapse.md]
  paper_appendix: present
---

# Audit unit 172 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_172.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage172_physical_slope_collapse.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 75, 461-487)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.txt`

## What the paper claims

Stage 172 ("Physical slope collapse", anchor MTDC-T9.1) rewrites the prior-stage linear grouped outlet obstruction pair in physical grouped-response variables. The card's `\stagefield{Output}` states verbatim: "Identifies \(\mathfrak K_1=-D_0u_2^{(1)}\) and \(\mathfrak G_1=D_0P_1\), reducing the linear grouped outlet problem to physical slopes." The notes (and appendix subsection "Physical slope collapse", eqs app-part05-u2-u4-P0-defs, KG-physical, dk-dg-physical) expand this into a chain of exact first-order identities: (1) the slope formulas `δu_2^(A) = -(δD_2 + u_2 δD_0)/D_0` and `δP_0^(A) = (δN_0 - P_0 δD_0)/D_0`; (2) the obstruction rewrite `K_A = -D_0 δu_2^(A) + (1/9 - u_2)δD_0`, `G_A = D_0 δP_0^(A)`, collapsing to `K_A = -D_0 δu_2^(A)` on the canonical branch `u_2 = 1/9`; (3) the hidden-even consistency relation `δu_4^(A) = (8/9) δu_2^(A)` under `δD_4 = (2/3)δD_2 + (1/27)δD_0`; (4) the direct outlet maps `δκ_W = -3(1-σ_*)/σ_* · δu_2` and `δγ_W = -(1-σ_*)/(9σ_*) · δP_0/P_0`; and (5) the even-preserving residual `Δ_Q = δP_0/P_0`. These are the distinct deliverables.

## What the script claims to verify

Both scripts (docstring "Checks 1-6") build the perturbed grouped quantities `D_{A,0}=D_0+ελ δD_0`, `D_{A,2}=D_2+ελ δD_2` (with `D_2=-u_2 D_0`), `N_{A,0}=N_0+ελ δN_0`, extract the first-order slopes `δu_2` and `δP_0` via a `Series`/`Normal` truncation divided by `ελ`, define the carried-forward obstruction pair `K_A = δD_2 + (1/9)δD_0` and `G_A = δN_0 - P_0 δD_0`, and then assert seven residual-zero identities: the two obstruction rewrites, the canonical collapse at `u_2=1/9`, the `δu_4 = (8/9)δu_2` even-consistency reduction (using independently extracted `δu_4` on the canonical `(u_2,u_4)=(1/9,4/81)` branch), the two direct outlet maps, and `Δ_Q = δP_0/P_0`. This is exactly the deliverable chain above.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `δu_2^(A) = -(δD_2 + u_2 δD_0)/D_0` | `delta_u2` extracted (sympy L52, wl L42); printed `(-dD0*u2-dD2)/D0` | match |
| `δP_0^(A) = (δN_0 - P_0 δD_0)/D_0` | `delta_P0` extracted (sympy L53, wl L43) | match |
| `K_A = -D_0 δu_2 + (1/9-u_2)δD_0` | expect_zero `K_A + D0*delta_u2 - (1/9-u2)*dD0` (sympy L62, wl L52) | match |
| `G_A = D_0 δP_0` | expect_zero `G_A - D0*delta_P0` (sympy L66, wl L56) | match |
| canonical collapse `K_A=-D_0 δu_2` at u2=1/9 | expect_zero subs u2->1/9 (sympy L72, wl L59) | match |
| `δu_4 = (8/9)δu_2` under even relation | expect_zero (sympy L101, wl L83) | match |
| `δκ_W = -3(1-σ)/σ · δu_2` | expect_zero (sympy L115, wl L92) | match |
| `δγ_W = -(1-σ)/(9σ) · δP_0/P_0` | expect_zero (sympy L120, wl L96) | match |
| `Δ_Q = δP_0/P_0` | expect_zero (sympy L126, wl L101) | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological script-side residual check, with matching constants (1/9, 4/81, 8/9, 2/3, 1/27, factors 3 and 1/9 in the outlet maps).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `simplify(K_A + D0*delta_u2 - (1/9-u2)*dD0) == 0` | obstruction rewrite (deliv 2) | yes |
| A2 | sympy | 66 | `simplify(G_A - D0*delta_P0) == 0` | obstruction rewrite (deliv 2) | yes |
| A3 | sympy | 72 | `simplify((K_A+D0*delta_u2).subs(u2,1/9)) == 0` | canonical collapse (deliv 2) | yes |
| A4 | sympy | 101 | `simplify((delta_u4 - 8/9 delta_u2).subs(dD4,even)) == 0` | hidden-even (deliv 3) | yes |
| A5 | sympy | 115 | `simplify((delta_kappa + 3(1-σ)/σ delta_u2).subs(u2,1/9)) == 0` | outlet κ (deliv 4) | yes |
| A6 | sympy | 120 | `simplify(delta_gamma + (1-σ)/(9σ P0) delta_P0) == 0` | outlet γ (deliv 4) | yes |
| A7 | sympy | 126 | `simplify(DeltaQ - delta_P0/P0) == 0` | even-preserving residual (deliv 5) | yes |
| B1-B7 | mathematica | 52,56,59,83,92,96,101 | `expectZero[...]` mirror of A1-A7 | same deliverables | yes |

All A-rows verified by hand below; none is tautological (each combines a derived slope with an independently-defined obstruction/coefficient, so a sign or factor error would produce a nonzero residual).

## Findings

### F1 — mathematica_transliteration

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:28-102`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py:37-129`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script rather than an independent re-derivation. Corresponding sections:

- Symbol setup. SymPy L37-43:
  `D2 = -u2*D0; P0 = N0/D0; DA0 = D0+eps*lam*dD0; DA2 = D2+eps*lam*dD2; NA0 = N0+eps*lam*dN0; u2A = -DA2/DA0; P0A = NA0/DA0`.
  Mathematica L32-40: `d2 = -u2*d0; p0 = n0/d0; dA0 = d0+eps*lam*dD0; dA2 = d2+eps*lam*dD2; nA0 = n0+eps*lam*dN0; u2A = -dA2/dA0; p0A = nA0/dA0`. Same variable choreography, same intermediate names, same order.
- Slope extraction. SymPy L52: `delta_u2 = simplify((sp.series(u2A,eps,0,2).removeO() - u2)/(eps*lam))`.
  Mathematica L42: `deltaU2 = FullSimplify[(Normal[Series[u2A,{eps,0,1}]] - u2)/(eps*lam), ...]`. Identical extraction recipe (truncate the series, subtract the base value, divide by `eps*lam`).
- Obstruction definitions and the seven checks. SymPy L58-59: `KA = dD2 + Rational(1,9)*dD0; GA = dN0 - P0*dD0`. Mathematica L48-49: `kA = dD2 + dD0/9; gA = dN0 - p0*dD0`. The seven `expectZero[...]` residuals (wl L52-101) are character-for-character the same residual expressions as the seven `expect_zero(...)` calls (py L62-129).

The only differences are CAS syntax (`FullSimplify[Together[Expand[...]]]` vs `simplify(expand(...))`) and casing of symbol names. No quantity is re-derived from an independent route (e.g., the `.wl` does not obtain `δu_2` by an alternative method such as implicit differentiation of `u_2 D_0 + D_2 = 0`, nor does it re-derive `K_A` from the microscopic bundle definitions). Both engines therefore exercise the *same algebra*, which defeats the purpose of a second engine.

**Why this matters:**
The second-engine policy requires both engines to reach the claim from the physical premises independently, so that a transcription or algebra error in one is caught by the other. A transliteration cannot catch such an error: any mistake copied into both scripts would pass in both. The agreement of the two outputs here is therefore weak evidence.

**Required change:**
Re-derive at least the load-bearing slope identities in the `.wl` by a route structurally distinct from the SymPy script, then assert the same final residuals. Minimal independent route that does not echo the `.py` choreography:
- Obtain `δu_2` by implicit differentiation of the defining relation `u_2^(A) D_{A,0} + D_{A,2} = 0` (i.e. `Solve`/`D` on `u2A*dA0 + dA2 == 0` rather than series-truncating `-dA2/dA0`), and `δP_0` from `P_0^(A) D_{A,0} - N_{A,0} = 0`.
- Obtain `δu_4` by implicit differentiation of `u_4^(A) D_{A,0}^2 - (D_{A,2}^2 - D_{A,0} D_{A,4}) = 0` rather than series-truncating the explicit quotient.
- Keep the seven final `expectZero` residual checks (wl L52-101) unchanged so the verified claim is identical; only the derivation path of `deltaU2`, `deltaP0`, `deltaU4Star` changes.

**Verification:**
After the fix, the `.wl` must no longer extract `deltaU2/deltaP0/deltaU4Star` via `Normal[Series[...]]/(eps*lam)`; instead via `Solve`/`D` on the implicit defining relations. The seven `PASS:` lines in the output must remain, and `Exit[0]` must hold. The printed slope forms (`deltaU2 = -((dD2+dD0*u2)/d0)`, etc.) must be unchanged.

### F2 — stale_banner_label (cosmetic; reported under tautological/insufficient-adjacent housekeeping)

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py:31`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:26`

**What's wrong:**
Both scripts print the banner `"STAGE 155 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM"` while this is unit/stage 172. The docstring and the seven assertions are unambiguously stage-172 content (the KA/GA/u4 algebra matches the stage-172 notes exactly), so this is a copy-paste label, not a wrong-stage script. It affects no assertion and no math result; it is purely a print string. (For reference, the notes file header itself reads "Stage 240" — the notes/in-paper numbering offset is a separate prose matter, not a script defect.)

**Why this matters:**
A misleading banner can mask which stage an output transcript belongs to during later auditing or when transcripts are diffed across units. No mathematical impact.

**Required change:**
In both scripts change the banner literal from `STAGE 155` to `STAGE 172` (sympy L31, wl L26). Do not change anything else.

**Verification:**
The first banner line in both `.txt` transcripts should read `STAGE 172 — PHYSICAL SLOPE COLLAPSE...` after re-run.

## Independent-derivation check (Mathematica)

Not independent — see F1. The `.wl` mirrors the `.py` symbol-for-symbol and step-for-step: same `D2=-u2*D0` construction, same `Series`-truncation slope extraction divided by `eps*lam`, same `KA=dD2+dD0/9`, and the same seven residual expressions in the same order. The two scripts differ only in CAS syntax and identifier casing. This is the textbook signature of transliteration.

## Engine cross-check

Both engines produce identical results and both `Exit 0`:
- `δu_2`: SymPy `(-dD0*u2 - dD2)/D0`; Mathematica `-((dD2 + dD0*u2)/d0)` — same expression.
- `δP_0`: SymPy `(D0*dN0 - N0*dD0)/D0**2`; Mathematica `(d0*dN0 - dD0*n0)/d0^2` — same.
- `δu_4` (canonical): SymPy `(-5*dD0 - 18*dD2 - 81*dD4)/(81*D0)`; Mathematica `-1/81*(5*dD0 + 18*dD2 + 81*dD4)/d0` — same.
- All seven residual checks print `0` / `PASS` in both.
No `engine_disagreement`. (But per F1 the agreement is weakened by transliteration.)

## Verdict justification

The math is fully sound and exactly aligned with the paper. I re-derived all seven identities by hand: the first-order slope formulas (`δu_2 = -(δD_2 + u_2 δD_0)/D_0`, `δP_0 = (δN_0 - P_0 δD_0)/D_0`), the obstruction rewrites (`K_A + D_0 δu_2 - (1/9-u_2)δD_0 ≡ 0`, `G_A - D_0 δP_0 ≡ 0`), the canonical collapse at `u_2=1/9`, the `δu_4 = (8/9)δu_2` reduction (I confirmed `δu_4 = -(5δD_0+18δD_2+81δD_4)/(81D_0)` and the `δD_4=(2/3)δD_2+(1/27)δD_0` substitution gives exactly zero), and the two outlet maps plus `Δ_Q = δP_0/P_0`. Constants (1/9, 4/81 with `D_4=-D_0/27` verified consistent, 8/9, factors 3 and 1/9) all match the notes and appendix. Attacks that failed: I tried to find a hidden tautology (none — each check couples a derived slope to an independently defined obstruction, so sign/factor errors would survive `simplify`); I checked the `Series`/`eps*lam` extraction for an off-by-one in series order (correct: the perturbation enters as `ε·λ`, so the `O(ε)` coefficient over `ε·λ` is the slope, and `λ` is nonzero symbolic); I checked the `D_4` value against the `u_4` definition (consistent); I checked symbol domains (`D0,N0,u2,sigma` nonzero — required for the divisions, and physically justified). The verdict is `findings` solely because (F1) the Mathematica script is a transliteration of the SymPy script rather than an independent derivation, weakening the second-engine guarantee, and (F2) both scripts carry a stale `STAGE 155` banner label. Neither is a math error and neither is stop-cold. F1's severity/applicability is partly a project-policy question (a deliberate mirror policy would downgrade it), so the directive frames it as an independence-restoring re-derivation that leaves the seven verified residuals unchanged.

## Self-test notes

I checked the variable-independence trap (the `Series` extractions depend genuinely on `eps`; no identically-zero derivatives), the trivial-case substitution for each `assert_zero` (all reduce to 0 only via real cancellation, not by construction), the symbol-domain trap (`nonzero` assumptions are necessary and physically justified), and the paper round-trip (F1's prescribed implicit-differentiation route reaches the same slope expressions and uses the same constants as the paper, introducing no new `paper_misalignment`; F2 touches only a print literal). No self-test failures.
