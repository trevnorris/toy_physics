---
unit_id: 170
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md]
  paper_appendix: present
---

# Audit unit 170 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_170.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md` (only one matching file)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 71, 417-433, 554)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.txt`

## What the paper claims

Stage 170 (`Linear grouped outlet map`, \StatusExactClosure) solves the linear grouped-`P2`
outlet problem left open by Stage 169. `\stagefield{Output}`: *"Reduces the direct outlet
deformations to \(\mathcal K_A=\delta D_{A,2}+\delta D_{A,0}/9\) and \(\mathcal G_A=\delta
N_{A,0}-P_0\delta D_{A,0}\)."* The notes carry the full derivation (more detail than the
terse card) and enumerate these distinct deliverables:
(1) the linear grouped transport laws `δu2 = -(δD_{A,2}+δD_{A,0}/9)/D0`,
`δu4 = -(δD_{A,4}+(2/9)δD_{A,2}+(5/81)δD_{A,0})/D0`, `δP0 = (δN_{A,0}-P0·δD_{A,0})/D0`;
(2) the exact direct-outlet map `δκ_W^{(A)} = 3(1-σ*)(δD_{A,2}+δD_{A,0}/9)/(σ*·D0)` and
`δγ_W^{(A)} = -(1-σ*)(δN_{A,0}-P0·δD_{A,0})/(9σ*·N0)`;
(3) the one-parameter even-consistency relation `δD_{A,4} = (2/3)δD_{A,2} + (1/27)δD_{A,0}`
(equiv. `δu4 = (8/9)δu2`);
(4) the grouped trace/anomaly (a-/b-) projector versions of the κ/γ maps;
(5) the weak-axisymmetric `(λ20,λ21,λ22)=(1,1/2,-1)` signature collapsing to two scalar
amplitudes `κ1`, `γ1`. The appendix Sec. "Direct outlet map" (eqs.
`app-part05-KG-defs-section`, `app-part05-dk-dg-outlet`) reproduces deliverables 1+2 verbatim.

## What the script claims to verify

Both scripts (docstring, banners, `expectZero`/`expect_zero` labels) verify exactly those
five deliverables. They build the canonical compensated branch `u2=1/9, u4=4/81, D2=-D0/9,
D4=-D0/27`, take the first-order coefficient of `u2,u4,P0` in the perturbation `eps`,
assert each equals the boxed transport law (Sec.1), invert the Stage-159 hybrid relations to
recover the κ_W/γ_W maps and assert the closed forms (Sec.2), solve `δu4=(8/9)δu2` for `dD4`
and assert the even-consistency relation (Sec.3), re-derive the a-/b- projector maps two ways
and assert agreement (Sec.4), and finally feed the lane-scaled `(1,1/2,-1)` defects through
the same maps and assert the κ1/γ1 amplitudes plus the signature ratios (Sec.5). Every check
is a non-tautological symbolic identity (residual must FullSimplify/simplify to 0).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| δu2,δu4,δP0 transport (notes Sec.1) | py 57-59 / wl 59-61 | match |
| δκ_W, δγ_W outlet map (Output; notes Sec.2; app eq dk-dg-outlet) | py 73-80 / wl 77-84 | match |
| even-consistency δD4=(2/3)δD2+(1/27)δD0 (notes Sec.3) | py 89-91 / wl 90-95 | match |
| grouped trace/anomaly a_κ,b_κ,a_γ,b_γ (notes Sec.4) | py 120-129 / wl 118-127 | match |
| weak-axisym signature + κ1,γ1 (notes Sec.5; app Xiload lanes) | py 156-163 / wl 148-157 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 57 | `expect_zero(du2 + (dD2+dD0/9)/D0)` | transport δu2 | yes |
| A2 | sympy | 58 | `expect_zero(du4 + ...)` | transport δu4 | yes |
| A3 | sympy | 59 | `expect_zero(dP0 - (dN0-P0 dD0)/D0)` | transport δP0 | yes |
| A4 | sympy | 73 | `expect_zero(dkappa_from_du2 - 3(1-σ)(dD2+dD0/9)/(σD0))` | κ_W map | yes |
| A5 | sympy | 77 | `expect_zero(dgamma_from_dP0 + ...)` | γ_W map | yes |
| A6 | sympy | 89 | `expect_zero(du4_from_kappa - (8/9)du2)` | even consistency | yes |
| A7 | sympy | 91 | `expect_zero(relation - ((2/3)dD2+dD0/27))` | δD4 relation | yes |
| A8-A11 | sympy | 120-128 | a/b κ,γ map-vs-projector agreement | trace/anomaly | yes |
| A12-A21 | sympy | 156-163 | lane κ_W/γ_W = εκ1/εγ1 + signature ratios | weak-axisym | yes |
| B1-B21 | mathematica | 59-61,77-84,90-95,118-127,148-157 | `expectZero[...]` mirrors of A1-A21 | all 5 deliverables | yes |

Every row traces to a specific paper-side deliverable; no orphaned assertions, none tautological.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent re-derivation**, not a transliteration. The two
load-bearing computational steps use genuinely different mechanisms from the `.py`, and the
script's inline comments document that this divergence was deliberate (this file was
de-transliterated in the BATCH 1 remediation commit `bda2107`):

1. **First-order coefficient extraction.**
   - SymPy (py 52-54): `sp.expand(sp.series(u2_full, eps, 0, 2).removeO()).coeff(eps, 1)` —
     Taylor-series expansion, then pluck the `eps^1` coefficient.
   - Mathematica (wl 54-56): `FullSimplify[(D[u2Full, eps] /. eps -> 0), ...]` — symbolic
     derivative evaluated at `eps=0`. Comment wl 52-53: *"Independent route: first-order
     coefficient via D[..., eps] /. eps -> 0, a different mechanism than the SymPy
     series().coeff(eps,1) this once mirrored."*

2. **Inversion of the hybrid relations.**
   - SymPy (py 69-70): placeholder-symbol idiom — `sp.solve(sp.Eq(sp.Symbol('du2sym'),
     du2_hyb), dkappa)[0].subs(sp.Symbol('du2sym'), du2)`.
   - Mathematica (wl 67-70): direct solve — `dkappa /. First[Solve[du2 == du2Hyb, dkappa]]`.
     Comment wl 66: *"Invert directly (no du2sym/dP0sym placeholder idiom — that was a SymPy
     tell)."*

The independence shows through in the saved outputs: the trace/anomaly forms are simplified
to different surface representations by the two engines (e.g. SymPy `a_kappa = -(aD0 +
9*aD2)*(sigma - 1)/(3*D0*sigma)` vs. Mathematica `a_kappa = -1/3*((aD0 + 9*aD2)*(-1 +
sigma))/(D0*sigma)`), and the engines independently reduce every residual to 0. The shared
banner titles, `expectZero` labels, and printed carry-forward block are cosmetic
presentation, not algebraic choreography. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit identical residual = 0 for all 21 checks (SymPy `... = 0`; Mathematica
`... = 0` + `PASS:` for each). Both transcripts end `# exit_code: 0`. The labeled
trace/anomaly forms agree up to engine-specific normalization (same rational function). No
`engine_disagreement`.

## Verdict justification

`clean`. I read the card, the (single) notes file, and the part-05 appendix rows before the
scripts, then attacked the algebra. Independently verified: the four transport laws
(δu2,δu4,δP0 first-order on the canonical branch D2=-D0/9, D4=-D0/27), the κ_W/γ_W maps
(inverting the Stage-159 hybrid relations), and the even-consistency relation
(`δu4=(8/9)δu2 ⟹ δD4=(2/3)δD2+(1/27)δD0`) — all reproduce the boxed notes/card/appendix
forms exactly. Every assertion is non-tautological (each is an algebraic identity the engine
must reduce, e.g. `dkappa_from_du2` is *solved* from the hybrid relation and then compared
against the independently-written closed form, so a wrong solve would fail). Symbol domains
are correct (`D0,N0,σ≠0`, `σ≠1` so the `(1-σ)` and `1/D0` denominators are well-posed; all
real). The weak-axisym section reuses the SAME verified maps with lane-scaled inputs — not a
re-postulation. The `.wl` is an independent re-derivation. Outputs are fresh (both `.txt`
mtimes 2026-05-29 00:05, after the 2026-05-28 23:58 script mtimes). Attacks that failed:
looked for a tautological `solve-then-assert-same-symbol` (no — solve target vs. asserted
closed form are written independently); looked for the lane section silently asserting
`0==0` because it reused already-derived `dkW[A]` on both sides (no — RHS is the
independently-written `eps·λ·κ1`); searched for stale `168`/`4107`/`100π²` constants (none
present in any stage-170 file). Paper claim matches script claim.

## Value Reconciliation (pass-2 augmentation)

All emitted RESULT values are closed-form symbolic identities (this stage pins no numeric
benchmark constants beyond the canonical rationals `u2=1/9, u4=4/81`, which are inputs/branch
definitions, not deliverables). Reconciliation against `.tex` card, notes `.md`, and the
part-05 appendix:

| value | source (py / wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δu2 = -(δD2+δD0/9)/D0` | py 57 / wl 59 (out L9) | notes Sec.1 box L99-104 | MATCH |
| `δu4 = -(δD4+(2/9)δD2+(5/81)δD0)/D0` | py 58 / wl 60 (out L10) | notes Sec.1 box L114-120 | MATCH |
| `δP0 = (δN0-P0·δD0)/D0` | py 59 / wl 61 (out L11) | notes Sec.1 box L106-111 | MATCH |
| `δκ_W^{(A)} = 3(1-σ*)(δD2+δD0/9)/(σ*·D0)` | py 73-75 / wl 77-79 (out L16) | Output (tex L15); notes Sec.2 L217-225; app eq dk-dg-outlet (tex L427-429) | MATCH |
| `δγ_W^{(A)} = -(1-σ*)(δN0-P0·δD0)/(9σ*·N0)` | py 77-79 / wl 81-83 (out L17) | Output (tex L15); notes Sec.2 L228-237; app eq dk-dg-outlet (tex L431-432) | MATCH |
| `δD4 = (2/3)δD2 + (1/27)δD0` (even consistency) | py 91 / wl 95 (out L23) | notes Sec.3 box L277-282; app eq hidden-even-D41-general (tex L530) | MATCH |
| `a_κ,b_κ = 3(1-σ*)(a/bD2 + a/bD0/9)/(σ*·D0)` | py 120-121 / wl 118-119 (out L28-29,36-37) | notes Sec.4 box L313-327 | MATCH |
| `a_γ,b_γ = -(1-σ*)(a/bN0 - P0·a/bD0)/(9σ*·N0)` | py 122-128 / wl 120-126 (out L30-31,38-39) | notes Sec.4 box L333-348 | MATCH |
| `κ1 = 3(1-σ*)(D2^(1)+D0^(1)/9)/(σ*·D0)` | py 144,156 / wl 138,148 (out L44) | notes Sec.5 box L395-403 | MATCH |
| `γ1 = -(1-σ*)(N0^(1)-P0·D0^(1))/(9σ*·N0)` | py 145,157 / wl 139,151 (out L45) | notes Sec.5 box L416-424 | MATCH |
| weak-axisym signature `(1,1/2,-1)` on κ_W/γ_W | py 159-163 / wl 154-157 (out L50-53) | card Checks item 2 (tex L21); notes Sec.5 L386-411; app Xiload-lanes (tex L547-551) | MATCH |

INTERNAL (scaffolding, no prose expectation): `u2=1/9`, `u4=4/81`, `D2=-u2·D0`,
`D4=-D0/27` (canonical-branch definitions — these DO appear in notes Sec.1 L88-95 anyway),
`eps`/`eps_l`/`epsL` perturbation symbols, `du2_hyb`/`dP0_over_P0_hyb`/`du4_from_hyb`
(Stage-159 hybrid relations carried in as inputs), `P0_ref=N0/D0`, the placeholder symbols
`du2sym`/`dP0sym`, and all `expectZero` residual-equals-0 check values.

reconciliation: complete; 11 deliverable values checked, 0 misaligned

## Self-test notes

Checked: (1) variable-independence — all first-order coefficients are taken w.r.t. `eps`,
on which `u2Full/u4Full/P0Full` genuinely depend, so no identically-zero-derivative trap;
the `du4=(8/9)du2` solve for `dD4` is well-posed since `du4` depends on `dD4`. (2) No
integrals, so parity is n/a. (3) Trivial-case: the solved κ_W/γ_W are compared against
independently written closed forms, and the lane section's RHS (`eps·λ·κ1`) is written from
the standalone `kappa1`/`gamma1` definitions rather than reusing the LHS, so the signature
checks are substantive, not `0==0`. Conclusion: no directive needed; verdict clean.
