---
unit_id: 011
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
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

# Audit unit 011 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_011.tex`
- notes: `(none)` (per-stage source notes under `notes/em_projected/step_011_*.md` were never committed; manifest reflects this)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_011.tex`) presents the projected Maxwell `P_2` bridge as the application of "bundle shifts to the one-pole/normalization compatibility surface". It defines `S = M + B_2 + Z_2`, `T = B_4 + Z_4` (eq:stage011-compat-base), states the isotropic compatibility surface after eliminating the static wall stiffness as `C = N_0/P_{0,target} - 3 S^2/T` (eq:stage011-compat), and gives the projected first variation `δC = n_0/P_{0,target} - 6 S z_2/T + 3 S^2 z_4/T^2` (eq:stage011-dcompat). The card explicitly states that the static conservative shift `z_0` cancels from δC (changing `P_0` but not the eliminated compatibility surface). The `\stagefield{Output}` line declares verbatim: "Stage~011 exports the `z_0`-cancellation and the projected `P_2` compatibility bridge \eqref{eq:stage011-dcompat}." The Part I appendix row supplements: "`P_2`-lane bridge from projected electromagnetic slots into the grouped normalization language."

## What the script claims to verify

The SymPy script perturbs the isotropic full-bundle backbone (`D_n, N_n`) by `D_n → D_n - eps z_n`, `N_n → N_n + eps n_n`, derives first-order shifts in `u_2, u_4, P_0, P_2, P_4` and the one-pole defect, then verifies the K-eliminated compatibility surface and its `z_0`-independent first variation in both fixed-target and transported-target forms. Beyond that it also verifies a constant-prefactor branch factorization (P2=P4=0 transport), the real-Y20 lane lambdas `(1, 1/2, -1)`, the weak-axisymmetric grouped trace identity `xbar = x_0` with `b = 3a`, the static Xi1 prefactor slope `na/N_0 + za/D_0`, and the per-lane `u_2` slope. The Mathematica `.wl` mirrors M1, M2, M3, M5 (and z0-independence), M6, M7, M8 (lane lambdas), M9 (grouped trace + b=3a), M10 (Xi1), M11 (u2 slope) — but does not duplicate the constant-prefactor branch factorizations or the dK shifts.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check(s) | Status |
|---|---|---|
| `S = M+B_2+Z_2, T = B_4+Z_4` (eq:stage011-compat-base) | implicit via symbols `S, T` and explicit `(M+B_2+Z_2)`/`(B_4+Z_4)` substitutions (sympy L110, L117) | match |
| `C = N_0/P_target - 3 S^2/T` (eq:stage011-compat) | sympy A4 (L141) + Mathematica M4 (L87) verify K-eliminated surface equals `compat_direct_p = (N_0+eps n_0)/P_target - 3 (S+eps z_2)^2/(T+eps z_4)`; eps→0 limit matches | match |
| `δC = n_0/P_target - 6 S z_2/T + 3 S^2 z_4/T^2` (eq:stage011-dcompat) | sympy A9 (L146) and Mathematica M5 (L93) | match |
| z_0 cancels from δC | sympy A13, A14 (L156-157), A15 (negative control for `dK_norm`); Mathematica M5/M7 z0-independence (L95, L116) | match |
| P_0 still moved by z_0 | sympy A2 (L103) records dP0 in closed form carrying z_0; A15 confirms dK_norm carries z_0; no explicit dP0/dz0 ≠ 0 assertion, but the recorded closed form makes the dependence manifest | match |
| (script-only) constant-prefactor branch P2/P4 transport (sympy L169-176) | sympy A18, A19 | extra |
| (script-only) real-Y20 lane lambdas `(1, 1/2, -1)` and grouped trace + b=3a | sympy A20-A24, Mathematica M8-M9 | extra |
| (script-only) static Xi1 prefactor slope `na/N_0 + za/D_0` | sympy A25, Mathematica M10 | extra |
| (script-only) per-lane u2 slope | sympy A26, Mathematica M11 | extra |
| (script-only) transported-target normalization K surface + compatibility | sympy A10, A11, A12, A14 and Mathematica M6, M7 | extra |

The paper's `Output` line is fully covered by sympy A1-A17 and Mathematica M1-M7 (the core projected-Maxwell P2 bridge with `z_0`-cancellation). Everything else in the script (constant-prefactor branch, weak-axisymmetric Y20 anomalies, Xi1 slope, u2 lane slope, transported-target variant) is downstream scaffolding the stage card does not enumerate as a deliverable. The appendix row's phrase "into the grouped normalization language" is loose enough to cover the dP0/dP2/dP4 closed forms, but the Y20/Xi1/u2 sections are clearly grouped-bundle anisotropy machinery that the stage card does not promise.

Set `paper_alignment: partial` — the script faithfully covers the paper's enumerated Output, but the script delivers substantially more than the paper card declares.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 102 | `assert_zero("delta u2", du2 - (D0 z2 - D2 z0)/D0^2)` | scaffolding for grouped normalization language (appendix row) | yes |
| A2 | sympy | 103 | `assert_zero("delta P0", dP0 - (D0 n0 + N0 z0)/D0^2)` | appendix row; supports "z0 changes P0 but not C" | yes |
| A3 | sympy | 113 | `assert_zero("one-pole first variation", d_pole - (D0 z4 - z0(B4+Z4) - 6 z2(M+B2+Z2)))` | substrate for eq:stage011-compat-base | yes |
| A4 | sympy | 141 | `assert_zero("compatibility surface after eliminating K", (K_norm_p - K_one_pole_p) - compat_direct_p)` | eq:stage011-compat (after eps→0) | yes |
| A5 | sympy | 142 | `assert_zero("one-pole K shift", dK_one_pole - (z0 + 6 S z2/T - 3 S^2 z4/T^2))` | substrate for eq:stage011-dcompat | yes |
| A6 | sympy | 143 | `assert_zero("normalization K shift", dK_norm - (z0 + n0/Ptarget))` | substrate for eq:stage011-dcompat | yes |
| A7 | sympy | 144 | `assert_zero("compatibility shift from competing K surfaces", d_compat - (dK_norm - dK_one_pole))` | substrate for eq:stage011-dcompat | yes |
| A8 | sympy | 145 | `assert_zero("compatibility first variation from eliminated surface", d_compat - d_compat_direct)` | substrate for eq:stage011-dcompat | yes |
| A9 | sympy | 146 | `assert_zero("compatibility first variation", d_compat_direct - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2))` | eq:stage011-dcompat (direct) | yes |
| A10 | sympy | 147 | transported-target K surface | none (extra) | yes |
| A11 | sympy | 148-151 | transported-target compat surface | none (extra) | yes |
| A12 | sympy | 152-155 | transported-target first variation | none (extra) | yes |
| A13 | sympy | 156 | `sp.diff(d_compat_direct, z0) == 0` | "z0 cancels from δC" | yes |
| A14 | sympy | 157 | `sp.diff(d_compat_transport, z0) == 0` | none (extra, transported variant) | yes |
| A15 | sympy | 158 | `assert_nonzero(sp.diff(dK_norm, z0))` | negative control protecting A13 | yes |
| A16 | sympy | 159-162 | sign-mutation negative control on fixed-target | negative control protecting A9 | yes |
| A17 | sympy | 163-166 | sign-mutation negative control on transported-target | negative control protecting A12 | yes |
| A18 | sympy | 175 | constant-prefactor P2 factorization | none (extra) | yes |
| A19 | sympy | 176 | constant-prefactor P4 factorization | none (extra) | yes |
| A20 | sympy | 183 | `lam20 - 1` | none (extra: Y20 lane) | yes |
| A21 | sympy | 184 | `lam21 - 1/2` | none (extra) | yes |
| A22 | sympy | 185 | `lam22 + 1` | none (extra) | yes |
| A23 | sympy | 195 | weak-axisymmetric grouped trace | none (extra) | yes |
| A24 | sympy | 196 | `bx - 3 ax` | none (extra) | yes |
| A25 | sympy | 204 | static Xi1 slope | none (extra) | yes |
| A26 | sympy | 211 | u2 lane slope | none (extra) | yes |
| M1 | mathematica | 59 | `deltaU2 - (D0 z2 - D2 z0)/D0^2` | appendix row | yes |
| M2 | mathematica | 62 | `deltaP0 - (D0 n0 + N0 z0)/D0^2` | appendix row | yes |
| M3 | mathematica | 68-71 | one-pole first variation | substrate for eq:stage011-compat-base | yes |
| M4 | mathematica | 86-88 | `compatFixed - compatFixedClosed` | eq:stage011-compat | yes |
| M5 | mathematica | 91-95 | fixed-target first variation + z0 independence | eq:stage011-dcompat + z0 cancellation | yes |
| M6 | mathematica | 105-108 | transported-target K surface | none (extra) | yes |
| M7 | mathematica | 110-116 | transported-target first variation + z0 independence | none (extra) | yes |
| M8 | mathematica | 129-133 | lane lambdas + same-sign selection | none (extra) | yes |
| M9 | mathematica | 141-142 | grouped trace + b=3a | none (extra) | yes |
| M10 | mathematica | 146 | static Xi1 slope | none (extra) | yes |
| M11 | mathematica | 150-153 | u2 lane slope | none (extra) | yes |

A13 and A14 differentiate closed-form expressions that already do not contain z0 by construction — the non-trivial K-elimination work is done by A4/A7/A8 (and Mathematica M4). Given that chain, A13 and A14 document the consequence rather than carry independent algebraic content. Acceptable because A4/A7/A8 form a substantive chain.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** `paper_missing_script_claim`
**Files:**
- paper card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_011.tex:37-39`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:117-211` (sections beyond the K-eliminated compatibility surface)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl:97-153` (transported-target variant, Y20 lanes, Xi1, u2 lane slope)

**What's wrong:**

Paper side, stage_011.tex L37-39:
> Stage~011 exports the `z_0`-cancellation and the projected `P_2` compatibility bridge \eqref{eq:stage011-dcompat}.

This is the only deliverable the card enumerates. The appendix row (L44) adds the phrase "into the grouped normalization language" but does not enumerate specific identities.

Script side, the sympy script delivers (beyond the card's declared output):
- transported-target K-surface and compatibility variant (L132-155);
- constant-prefactor branch factorizations imposing `P_2 = 0`, `P_4 = 0` (L168-176);
- real-Y20 weak-axisymmetric lane lambdas `(1, 1/2, -1)` and the grouped trace `xbar = x_0` with `b = 3a` (L179-196);
- static Xi1 prefactor slope `na/N_0 + za/D_0` (L199-204);
- per-lane `u_2` slope `(D_0 z_{2a} - D_2 z_a)/D_0^2` (L207-211).

The Mathematica `.wl` mirrors most of these (M6-M11). None of these appear in the paper card or the appendix row.

**Why this matters:**

This is the inverse direction of the usual paper_misalignment: the script does more than the paper claims, not less. Three concrete risks:
1. If the paper later cites stage 011 as the source of the Y20 lane signature `(1, 1/2, -1)` or of `Xi1^{(proj,static)} = na/N_0 + za/D_0`, that citation will fail because the paper card does not document the result there. Downstream stages or main-text chapters may already do this.
2. If the script was written against a stale paper revision that did enumerate these items and the card was later trimmed (or vice-versa), there is no notes/ trail to disambiguate.
3. The constant-prefactor branch factorization (A18, A19) and the transported-target variant (A10-A12, A14) carry algebraic content that may be load-bearing for stages 012-021 ("Projected Maxwell primitive bridge", "Reduced Maxwell/mixed one-port", etc.) — if so, those results deserve a paper anchor here rather than being hidden in this script.

The script's internal math is sound; nothing in the extra sections is wrong. The misalignment is purely an enumeration gap between paper and script.

**Required change:**

Goes to user resolution via the directive's `## Resolve before fix_loop` block. Codex must not silently edit either side.

**Verification:**

Once the user has chosen a direction, the verifier confirms either (a) the stage_011.tex card lists the additional deliverables in its `Output` paragraph (and any new equation labels), or (b) the extra script sections are removed and any downstream consumer is re-pointed to the stage that does carry them.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally close to the SymPy script but not a line-by-line transliteration:
- SymPy uses `sp.series(expr, eps, 0, 2).removeO()` then divides by `eps`; Mathematica uses `Coefficient[Normal[Series[expr - (expr /. var → 0), {var, 0, 1}]], var, 1]`. Different idiom, same operation.
- SymPy imports `sympy.physics.wigner.gaunt` directly; Mathematica builds the Gaunt coefficient by hand as `Sqrt[5/(4 Pi)] * ThreeJSymbol[{2,0},{2,0},{2,0}] * ThreeJSymbol[{2,m1},{2,m2},{2,m3}]` (L118-124). This is a genuine independent derivation — the 3j-product is constructed from elementary objects rather than borrowed from a packaged Gaunt routine, so engine bugs in `gaunt` are not shared.
- SymPy verifies the constant-prefactor branch factorizations (A18, A19); Mathematica omits them entirely. Mathematica also omits the `dK_one_pole` / `dK_norm` intermediate steps and goes straight to `compatFixed = FullSimplify[kNorm - kPole]` then asserts the closed form (M4). This is independent algebra reaching the same answer.

Verdict on transliteration: **not a transliteration**. Same target identities, parallel topology, but different machinery for Gaunt and a different intermediate path through the K-elimination. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce `STATUS: PASS` on the overlapping claims:

| Claim | SymPy residual (output) | Mathematica residual (output) |
|---|---|---|
| delta u2 | implicitly 0 (assertion passes) | `M1 delta u2 residual = 0` |
| delta P0 | implicitly 0 | `M2 delta P0 residual = 0` |
| one-pole first variation | implicitly 0 | `M3 one-pole first variation residual = 0` |
| fixed-target K-eliminated surface | implicitly 0 | `M4 fixed-target K-eliminated surface residual = 0` |
| fixed-target compatibility first variation | implicitly 0 | `M5 fixed-target compatibility first variation residual = 0` |
| z0 independence (fixed) | implicitly 0 | `M5 fixed-target z0 independence residual = 0` |
| transported-target K surface | implicitly 0 | `M6 transported-target normalization K surface residual = 0` |
| transported-target first variation | implicitly 0 | `M7 transported-target compatibility first variation residual = 0` |
| z0 independence (transported) | implicitly 0 | `M7 transported-target z0 independence residual = 0` |
| Y20 lane lambdas (1, 1/2, -1) | implicitly 0 | `M8 real-Y20 lambda(0,1,2) residual = 0` |
| same-sign m=1,2 selection | implicitly 0 | `M8 same-sign selection m=1,2 residual = 0` |
| grouped trace, b=3a | implicitly 0 | `M9 weak-axisymmetric grouped trace, b=3a residual = 0` |
| static Xi1 slope na/N0 + za/D0 | implicitly 0 | `M10 static Xi1 slope residual = 0` |
| u2 lane slope (D0 z2a - D2 za)/D0^2 | implicitly 0 | `M11 u2 projected-Maxwell lane slope residual = 0` |

Engines agree. `engines_agree: true`.

## Verdict justification

The script's core mathematics holds up against the paper. All eight identities the stage card enumerates (compatibility-base symbols, K-eliminated compatibility surface, projected first variation, z0 cancellation, and the implicit "P0 still moved by z0") are independently verified by both engines, with parallel but non-transliterated derivations. Negative-control assertions (A15-A17) protect the core claims against trivial passes. Both saved outputs are fresh (May 21) vs script mtimes (May 4 sympy, May 21 .wl).

The single finding is `paper_misalignment / paper_missing_script_claim`: the script delivers substantially more than the stage card declares as Output (constant-prefactor branch, Y20 weak-axisymmetric grouped signature, Xi1 slope, u2 lane slope, transported-target variant). The script's extras are mathematically sound but are not anchored anywhere in the paper card or appendix row. The direction of fix (extend paper card vs. trim script) is the user's call, so the directive routes a `## Resolve before fix_loop` block rather than instructing Codex to edit either side.

Verdict: `findings`. No stop-cold flag — the math is correct, only the paper-side enumeration is incomplete.

## Self-test notes

I checked: (1) The compatibility surface closed form is bit-for-bit identical between script (`compat_direct_p = (N0+eps n0)/Ptarget - 3 (S+eps z2)^2/(T+eps z4)`) and paper card (eq:stage011-compat / dcompat with eps→0 / first variation). (2) `S = M+B2+Z2, T = B4+Z4` are abstract symbols in the script but the one-pole defect at L110 spells `D0*(B4+Z4) - 3*(M+B2+Z2)^2` explicitly, so the symbol identification is unambiguous. (3) The Y20 lane lambdas (1, 1/2, -1) are computed via independent gaunt/3j routes in the two engines; both report 0 residual against the targets. (4) The `sp.diff(d_compat_direct, z0) == 0` check is structural (since `compat_direct_p` is built without z0) — the load-bearing work is in A4/A7/A8 showing that the K-eliminated `K_norm_p - K_one_pole_p` actually equals `compat_direct_p`, which is non-trivial. (5) No `simplify(positive=True)` traps; perturbation symbols `z0, z2, z4, n0, n2, n4` are unconstrained, which is correct because their signs are not fixed by the physical setup.
