---
unit_id: 121
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage121_geometric_r_selection.md]
  paper_appendix: present
---

# Audit unit 121 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_121.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage121_geometric_r_selection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1276 is an `\input{stages/stage_121}` include only — no separate status row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.txt`

## What the paper claims

The card body states the stage's bottom-line result verbatim: "Identifying the mixed D/N tube with the throat span gives \(\mathfrak r_{\rm geom}(L/a)=\sqrt{12(L/a)^2/\pi^2-1}\)." The notes elaborate the full derivation chain: (1) the identification \(L_W=L\) of the auxiliary mixed D/N tube with the actual throat axial span; (2) substituting the upstream (Stage 116) tube-length formula \(L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}\) and solving \(L_W=L\) gives the exact geometric branch law \(\mathfrak r_{\rm geom}(L/a)=\sqrt{\frac{12}{\pi^2}(L/a)^2-1}\); (3) the existence condition \(L/a\ge \frac{\pi}{2\sqrt3}\approx 0.9069\); (4) the preferred-aspect value at \(L/a=37/20\): \(\mathfrak r_{F1}=\frac{\sqrt{4107-100\pi^2}}{10\pi}\approx 1.778\) and \(r_c^{F1}=\mathfrak r_{F1}^2\approx 3.161\); (5) the mixed-tube pole \(\Omega_W=\frac{\pi c_s}{2L}\). The card is terse and reports only deliverable (2) in body form; the numeric deliverables (3)–(5) live in the notes.

## What the script claims to verify

Both scripts derive `r_geom` by solving the upstream tube-length relation `L_W = (pi*a/2)*sqrt((1+r^2)/3)` for the radius `r` and substituting `L_W -> L`, then verify six identities: (A1) the solved root squares to the closed form `12 L^2/(pi^2 a^2) - 1`; (A2) re-substituting that root into the original tube-length formula returns `L` (round-trip on the correct branch); (A3) at `L/a = 37/20` the value equals `sqrt(4107 - 100 pi^2)/(10 pi)`; (A4) its square equals `4107/(100 pi^2) - 1`; (A5) the half-wave pole at `L_W = L` is `pi c_s/(2 L)`; (A6) `r_geom` vanishes at the existence threshold `L/a = pi/(2 sqrt 3)`. Outputs print the numeric values `r_F1 ≈ 1.77799353547…` and `r_c ≈ 3.16126101219…` and the threshold `≈ 0.90689968…`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `r_geom(L/a) = sqrt(12(L/a)^2/pi^2 - 1)` (card body + notes box) | A1 (closed form) + A2 (round-trip) | match |
| identification `L_W = L` | implicit via `.subs(LW, L)` + A2 round-trip | match |
| existence threshold `L/a >= pi/(2 sqrt 3) ≈ 0.9069` (notes) | A6 (`r_geom` vanishes there) + printed threshold | match |
| preferred aspect `L/a = 37/20` → `r_F1 = sqrt(4107-100pi^2)/(10pi) ≈ 1.778` (notes) | A3 + printed `r_F1` | match |
| `r_c^F1 = r_F1^2 ≈ 3.161` (notes) | A4 + printed `r_c` | match |
| `Omega_W = pi c_s/(2 L)` (notes box) | A5 | match (trivial-substitution check; see Self-test) |

Every paper-side deliverable has a corresponding script-side check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 25-28 | `simplify(r_geom^2 - (12L^2/(pi^2 a^2) - 1)) == 0` | `r_geom` closed form | yes |
| A1 | math | 62-65 | `expectZero[derivedRadius^2 - (12 ell^2/(Pi^2 throatScale^2) - 1)]` | `r_geom` closed form | yes |
| A2 | sympy | 30-33 | `LW_formula.subs(r, r_geom) - L == 0` | `L_W=L` / round-trip | yes |
| A2 | math | 67-70 | `expectZero[(stage99TubeLength/.radius->derivedRadius) - ell]` | `L_W=L` / round-trip | yes |
| A3 | sympy | 40-44 | `simplify(r_F1 - sqrt(4107-100pi^2)/(10pi)) == 0` | `r_F1` value | yes |
| A3 | math | 81-84 | `expectZero[rF1Derived - rF1Target]` | `r_F1` value | yes |
| A4 | sympy | 50-54 | `simplify(rc_F1 - (4107/(100 pi^2) - 1)) == 0` | `r_c^F1` value | yes |
| A4 | math | 91-94 | `expectZero[rcF1Derived - rcF1Target]` | `r_c^F1` value | yes |
| A5 | sympy | 61-64 | `OmegaW.subs(LW,L) - pi c_s/(2L) == 0` | `Omega_W` pole | partial (trivial substitution) |
| A5 | math | 103-106 | `expectZero[omegaAtEqualTubeLength - Pi soundSpeed/(2 ell)]` | `Omega_W` pole | partial (trivial substitution) |
| A6 | sympy | 69-72 | `simplify(r_geom.subs(L, a*threshold)) == 0` | existence threshold | yes |
| A6 | math | 112-116 | `expectZero[derivedRadius/.ell->thresholdAspect*throatScale]` | existence threshold | yes |

## Findings

None. The trivial-substitution character of A5 (both engines) and the stale `stage99` variable-name label in the `.wl` are discussed below as sub-threshold observations; neither rises to a finding.

### Sub-threshold observation 1 — A5 (`Omega_W`) is a near-tautological substitution

In both engines A5 computes `pi*c_s/(2*L_W)`, substitutes `L_W -> L`, and then asserts the result equals `pi*c_s/(2*L)` — i.e. it checks `pi c_s/(2L) - pi c_s/(2L) == 0` after the substitution. The check cannot fail and merely documents the boxed notes identification `Omega_W = pi c_s/(2L)`. This is the *only* genuinely trivial check in the unit. It does not rise to an `insufficient_verification` finding because the paper's own claim for this deliverable is itself just the substitution `L_W=L` into a pole formula carried unchanged from upstream — there is no nontrivial algebra the script is failing to exercise. Noted, not filed.

### Sub-threshold observation 2 — stale `stage99TubeLength` label in `.wl`

The `.wl` names the upstream tube-length formula `stage99TubeLength` (lines 45/46/69), but the notes attribute that formula to **Stage 116** (notes lines 24, 33). "99" is the pre-renumber (−17) image of "116", i.e. the known numbering-drift in the script/output band. The label is functionally inert — the formula `(Pi*throatScale/2)*Sqrt[(1+radius^2)/3]` is mathematically correct and matches the notes verbatim — and per the deferred-dedicated-pass policy for the script/output numbering band this is not filed as a blocking finding. Recorded here for that pass.

## Independent-derivation check (Mathematica)

The `.wl` was added in the 2026-06-01 retro-sweep (commit `251639c`, "red-team retro-sweep {121,122,123}: dual-engine .wl backfill") to an originally SymPy-only stage, and it is a **genuine independent route**, not a transliteration. The load-bearing root-selection step is materially different in the two engines:

- SymPy (line 21): `r_geom = sp.solve(sp.Eq(LW, LW_formula), r)[0]` — blindly takes the first root, no branch logic, no existence guard.
- Mathematica (lines 46-54): `Solve[..., radius, Reals]`, then `Select`s the root satisfying `# > 0` *under explicit branch assumptions*, and `fail`s unless exactly one positive root survives (`If[Length[positiveRadiusRoots] =!= 1, fail[...]]`).

The Mathematica script additionally constructs its own `branchAssumptions` block (lines 36-43) encoding the existence-threshold inequalities `ell/throatScale > Pi/(2 Sqrt[3])`, `-3 Pi^2 throatScale^2 + 36 ell^2 > 0`, etc., that have no counterpart in the SymPy script. Variable naming is independent (`ell/throatScale/tubeLength/radius/soundSpeed` vs `L/a/L_W/r/c_s`), and the `expectZero` harness uses `FullSimplify[Together[Expand[stripCE[...]]]]` with a `ConditionalExpression` strip, structurally unlike SymPy's `simplify(expand(...))`. The downstream checks (A1–A6) verify the same *targets* — which is correct and required: both engines must confirm the same physics — but they are constructed from independently-derived expressions. No `mathematica_transliteration` finding.

## Engine cross-check

The two saved outputs agree on every emitted value:

| value | sympy `.txt` | mathematica `.txt` |
|---|---|---|
| `r_geom(L/a)` | `sqrt(12*L**2 - pi**2*a**2)/(pi*a)` | `Sqrt[12*ell^2 - Pi^2*throatScale^2]/(Pi*throatScale)` |
| `r_F1` (symbolic) | `sqrt(4107 - 100*pi**2)/(10*pi)` | `Sqrt[4107 - 100*Pi^2]/(10*Pi)` |
| `r_F1` (numeric) | `1.7779935354749781185` | `1.77799353547497811851…` |
| `r_c(F1)` | `-1 + 4107/(100*pi**2)` | `-1 + 4107/(100*Pi^2)` |
| `r_c(F1)` (numeric) | `3.1612610121908122732` | `3.16126101219081227320…` |
| `Omega_W(L_W=L)` | `pi*c_s/(2*L)` | `(Pi*soundSpeed)/(2*ell)` |
| threshold | `0.90689968211710892530` | `0.90689968211710892529703912…` |

All six `expectZero`/`expect_zero` residuals are `0` in both engines; the `.wl` prints six `PASS:` lines plus the final passed banner. `engines_agree: true`.

## Verdict justification

Clean. I read the card, the notes, and the appendix include, built the model of the geometric branch law, then attacked the scripts. Attacks tried and failed: (1) checked whether A1/A2 are tautological — they are not; `r_geom` is produced by `solve` and the two checks (closed-form square + round-trip) would both expose a wrong-root selection, which is exactly the failure mode the Mathematica positive-branch `Select`/`fail` guards against. (2) Verified the `4107` constant: `12·(37/20)² = 4107/100`, so the F1 value and `r_c` checks are anchored to the correct upstream aspect ratio `37/20` (= the notes' preferred `L/a`). (3) Checked symbol domains — all symbols `positive=True, real=True` (SymPy) / `> 0, Reals` (Mathematica), consistent with physical lengths/radius/sound-speed and with the notes' setup; positivity is what makes the single-positive-root selection valid. (4) Checked the threshold check A6 substitutes the exact zero-locus `L = a·pi/(2 sqrt 3)` and correctly returns 0. (5) Confirmed the retro-sweep `.wl` is an independent route, not a port. The only soft spots — the trivial A5 substitution and the stale `stage99` label — are sub-threshold and do not affect any verified deliverable. Every emitted deliverable value reconciles to the card and/or notes (see below).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_geom(L/a) = sqrt(12(L/a)^2/pi^2 - 1)` | py L23/27 + wl L60/64; sympy.txt L5, math.txt L5 | `.tex:16` (card body), `.md:39-44` (boxed) | MATCH |
| `L_W = L` (identification) | py L22 `.subs(LW,L)`; A2 round-trip | `.md:27` (boxed `L_W=L`) | MATCH |
| existence threshold `pi/(2 sqrt 3) ≈ 0.90689968…` | py L66-67 + wl L108-109; sympy.txt L16, math.txt L21 | `.md:48-54` (boxed `>= pi/(2 sqrt 3) ≈ 0.9069`) | MATCH |
| `L/a = 37/20` (preferred aspect) | py L35 + wl L72 | `.md:56-61` (`L/a = 37/20`) | MATCH |
| `r_F1 = sqrt(4107-100 pi^2)/(10 pi) ≈ 1.77799…` | py L40/37-38 + wl L77-79; sympy.txt L8-9, math.txt L10-11 | `.md:62-72` (boxed, `≈ 1.77799353547498`) | MATCH |
| `r_c^F1 = 4107/(100 pi^2) - 1 ≈ 3.16126…` | py L50/46-48 + wl L87-89; sympy.txt L11-12, math.txt L14-15 | `.md:74-79` (boxed, `≈ 3.16126101219081`) | MATCH |
| `Omega_W = pi c_s/(2 L)` | py L58-63 + wl L96-105; sympy.txt L14, math.txt L18 | `.md:84-91` (boxed `Omega_W = pi c_s/(2L)`) | MATCH |

All seven deliverable values appear in the notes (and the headline `r_geom` form also in the card body) and agree to full printed precision in both engines. The terse `.tex` card legitimately reports only the `r_geom` closed form; the numeric deliverables living in the `.md` notes count as MATCH per the augmentation guards (a value correctly in `.md` is a MATCH even when the card omits it).

INTERNAL (scaffolding, no finding): `LW_formula`/`stage99TubeLength` (upstream-carried tube-length formula, not a Stage-121 deliverable), `radiusRoots`/`positiveRadiusRoots`/`derivedRadius` (branch-selection intermediates), all `expectZero`/`expect_zero` residuals (`= 0`), `PASS:` flags, and the final passed-banner lines.

## Self-test notes

Variable-independence: no `sp.diff`/`D[]` derivatives in either script, so the identically-zero-derivative trap does not apply. Symmetry/parity: no unbounded integrals, so the parity trap does not apply. Trivial-case pre-check: A1/A2/A3/A4/A6 each substitute concrete values (`37/20`, the threshold) and the residuals reduce to literal `0`, confirmed against both saved outputs; A5 is a pure substitution and is the one near-tautological check (noted, sub-threshold). No directive is written — `findings_count: 0`, verdict clean.
