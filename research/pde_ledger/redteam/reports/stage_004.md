---
unit_id: 004
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

# Audit unit 004 red-team report (v2, paper-grounded)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_004.tex`
- notes: (none) — `notes/em_projected/step_01_*.md` is not committed; the script docstring names `step_01_projected_maxwell_readme.md` as the source note but no such file exists in the repo. Audit proceeds from paper card + appendix row + scripts.
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row for stage 004 at lines 30-31)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.txt`

## What the paper claims

The paper card for stage 004 is a brief bundle-index card. The Output line states verbatim: "Stage~004 fixes the source ordering for the projection-first EM bundle and anchors the projection-by-parts identity used in the next stage." The body of the card adds two pieces: (i) a "Bundle role" paragraph enumerating three rules — projection by parts creates a real transverse leakage/source term, homogeneous Maxwell projects in measured fields while inhomogeneous equations carry source-coupled flux fields, and reduction is a matched channel (not the primary parent derivation); and (ii) the projection identity, written algebraically as `∂_w(W Q) − W ∂_w Q − (∂_w W) Q = 0` (eq:stage004-projection-ibp), with the comment that after integration this is the bookkeeping identity converting the parent w-flux into a boundary term plus a kernel-gradient leakage term. The appendix row summarizes the stage as "Bundle ordering and projection-by-parts identity for the projection-first electromagnetic sector," status `\StatusExact{}`. The stage card explicitly defers the homogeneous/inhomogeneous Maxwell split, the matched zero-mode reduction, and the projection/reduction coupling-mismatch quantification to later stages (005, 006, 007 respectively per the appendix table). So the load-bearing deliverables here are: (a) the projection-by-parts identity, and (b) the bundle source ordering (i.e., stages 005-007 exist and are sequenced).

## What the script claims to verify

The sympy script frames itself as a "bundle-level audit for step_01_projected_maxwell_readme.md" — i.e., the README that ties stages 004-007 together — and performs six substantive checks plus one inventory check: (1) the three downstream scripts (stage 005 covariant, stage 006 vector, stage 007 projection/reduction comparison) exist on disk; (2) the integrated projection-by-parts identity `∫W·Q_w dw + ∫W_w·Q dw = 0` on concrete decaying Gaussian profiles; (3) cyclic Bianchi reduces to vector Faraday signs for three components, with `F_{23}=B_1, F_{30}=E_3, F_{02}=-E_2` and cyclic; (4) `∫exp(-w²/λ²) dw = √π·λ`; (5) `∫exp(-w²/λ²)² dw = (√(2π)/2)·λ`; (6) matched-kernel overlap `∫(Z/∫Z)·Z dw = √2/2`; (7) delta-source projection/reduction ratio `μ_{proj,δ}/μ_{red} = √2`. The Mathematica script (M1–M6) mirrors checks (2)–(7) (no inventory check) with different variable names but identical mathematical content.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| eq:stage004-projection-ibp `∂_w(WQ)-W∂_wQ-(∂_wW)Q=0` (algebraic Leibniz, integrated as boundary + leakage) | sympy IBP block lines 38-44 / mathematica M1 lines 16-26: integrated decay version `∫W Q_w + ∫W_w Q = 0` | match (script verifies the integrated decay case, which is the substantive bookkeeping content the paper says the identity supports) |
| Bundle source ordering: stages 005, 006, 007 sequenced as the projection-first EM bundle | sympy lines 23-30: inventory FileNotFoundError on missing 005/006/007 scripts | match |
| Bundle role rule (i): projection by parts → transverse leakage/source term | covered by IBP check above | match |
| Bundle role rule (ii): homogeneous Maxwell in measured fields, inhomogeneous in source-coupled flux | sympy Faraday block lines 54-67 / mathematica M2 lines 28-69 — only writes the homogeneous side in two notations and checks they match; does NOT exercise any projection/source-flux content | mismatch (tautological — see F1) |
| Bundle role rule (iii): reduction is a matched channel, not primary | sympy lines 69-82 / mathematica M3-M6 — Gaussian closed forms + matched-kernel √2 coupling ratio | extra (this is forward-touching content properly belonging to stage 007/008/009 per the appendix table; not a paper-card claim for stage 004 itself, but consistent with the script's "bundle README" framing in its docstring) |

Dominant pattern: the stated `Output` (IBP identity + bundle ordering) is faithfully covered. The Faraday block is structurally aligned with rule (ii) but is tautological (F1 below). The Gaussian/√2 block is bundle-scope scaffolding that the paper card defers to later stages; the script's docstring honestly frames itself as a bundle-README audit, and there is no direct contradiction with the paper, so `paper_alignment: aligned` (extras are scope-honest, not misaligned).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A0 | sympy | 28-30 | `FileNotFoundError` if 005/006/007 scripts missing | bundle ordering | yes (meta inventory) |
| A1 | sympy | 41-44 | `∫W·Q_w + ∫W_w·Q == 0` on Gaussian W, odd Gaussian Q | IBP identity (integrated decay case) | yes |
| A2 | sympy | 58-60, 66-67 | Faraday component 1: F-form expression `−` vector-form expression `== 0` | bundle role (ii) — claimed | no (tautological substitution; see F1) |
| A3 | sympy | 61-62, 66-67 | Faraday component 2 (same shape) | bundle role (ii) | no (tautological; see F1) |
| A4 | sympy | 63-64, 66-67 | Faraday component 3 (same shape) | bundle role (ii) | no (tautological; see F1) |
| A5 | sympy | 79 | `Z_int − √π·λ == 0` | extra (stage 007 scaffold) | yes |
| A6 | sympy | 80 | `Z2_int − √(2π)·λ/2 == 0` | extra (stage 007 scaffold) | yes |
| A7 | sympy | 81 | `I_WZ − √2/2 == 0` | extra (stage 007 scaffold) | yes |
| A8 | sympy | 82 | `μ_proj_δ/μ_red − √2 == 0` | extra (stage 007 scaffold) | yes |
| M1 | mathematica | 18-26 | integrated IBP on decaying Gaussian/odd-Gaussian | IBP identity | yes |
| M2a | mathematica | 32-43 | Faraday comp 1: cyclic sum on F-symbols `−` vector form `== 0` | bundle role (ii) | no (tautological; see F1) |
| M2b | mathematica | 45-56 | Faraday comp 2 (same shape) | bundle role (ii) | no (tautological; see F1) |
| M2c | mathematica | 58-69 | Faraday comp 3 (same shape) | bundle role (ii) | no (tautological; see F1) |
| M3 | mathematica | 72-82 | profileArea = √π·λ | extra (stage 007 scaffold) | yes |
| M4 | mathematica | 85-94 | profileSelfMass = √(2π)·λ/2 | extra | yes |
| M5 | mathematica | 97-107 | overlapValue = √2/2 | extra | yes |
| M6 | mathematica | 110-121 | pointSourceCoupling/volumeReducedCoupling = √2 | extra | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py:54-67`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl:32-69`

**What's wrong:**

Both engines' "Faraday / Bianchi" blocks (A2-A4 in sympy, M2a-M2c in mathematica) substitute the antisymmetric-tensor components by hand and then compare the substituted expression against the vector-form expression. The sympy block defines

```
F23, F31, F12 = B1, B2, B3
F10, F20, F30 = E1, E2, E3
F01, F02, F03 = -E1, -E2, -E3

faraday[0] = ( sp.diff(F23, t) + sp.diff(F30, y) + sp.diff(F02, z) )
           - ( sp.diff(B1, t)  + sp.diff(E3, y)  - sp.diff(E2, z) )
```

After substitution, the LHS is literally `sp.diff(B1,t) + sp.diff(E3,y) + sp.diff(-E2,z) = sp.diff(B1,t) + sp.diff(E3,y) - sp.diff(E2,z)`, which is exactly the RHS character-for-character. The residual is identically zero by symbol substitution, independent of any physical content. The mathematica M2 block has the identical structure (line 32-43): `twoForm23 = B1[...]; twoForm30 = E3[...]; twoForm02 = -E2[...]; m2LeftOne = D[twoForm23,t] + D[twoForm30,y] + D[twoForm02,z]; m2RightOne = D[B1[...],t] + D[E3[...],y] - D[E2[...],z]; assert m2LeftOne - m2RightOne =!= 0`. The two sides differ only in notation, not algebra.

The paper-side claim these blocks purport to support is the appendix row's "homogeneous Maxwell projects in measured fields, while the inhomogeneous equations carry source-coupled flux fields" (bundle role rule ii). A meaningful test of that distinction would either (a) derive the homogeneous block from `dF = 0` applied to `F = dA` for some potential A (then the Bianchi cyclic identity follows from `d² = 0` / Schwarz, which is non-trivial because `sp.diff` does NOT a priori commute on user-declared `Function` objects unless declared real and smooth — this is testable content), or (b) explicitly verify the cyclic-sum identity using an antisymmetric F-tensor parameterized by an independent set of components (not by substituting the components directly into the vector form on both sides).

The check as written cannot fail. It does not exercise the Bianchi identity, the cyclic-sum structure, or the homogeneous/inhomogeneous distinction. The only way it could fail is a literal typo in the manual transcription on the RHS — and that would catch a notational bug, not a physics bug.

**Why this matters:**

The paper card lists three "rules" the bundle is supposed to fix; the Faraday block is the only piece of the script that addresses rule (ii). With this block being tautological, rule (ii) is in fact *unverified* by any check in the stage 004 scripts. Either stage 004 should exercise rule (ii) for real, or the bundle-role enumeration should acknowledge that rule (ii) is checked by stage 005/006 scripts and the Faraday block here is decorative (and its assert should be removed or relabeled). The current state credits coverage that does not exist.

A weaker but still useful fix: replace the trivial substitution with an actual cyclic-Bianchi test using an antisymmetric `F[μ,ν]` symbolic array (e.g., a `MutableDenseNDimArray` in sympy / `Array` in Mathematica) whose components are *independent* function symbols, and verify that imposing `F[i,j] = -F[j,i]` plus assigning the E,B components by the standard map yields the Maxwell-Faraday equation as a cyclic-sum identity over derivatives of those independent components. Then the check has content: a wrong sign in the E,B↔F map would produce a nonzero residual.

**Required change:**

Replace the tautological substitute-and-compare blocks with a check that exercises the cyclic Bianchi `∂_[α F_{βγ]} = 0` for `F = dA` from a vector potential, then specializes to the E,B map and reads off the three vector-form Maxwell-Faraday components. Concretely (sympy):

```python
A0, A1, A2, A3 = (sp.Function(f"A{i}")(t, x, y, z) for i in range(4))
coords = (t, x, y, z)
def F(mu, nu):
    return sp.diff([A0, A1, A2, A3][nu], coords[mu]) - sp.diff([A0, A1, A2, A3][mu], coords[nu])
# Cyclic Bianchi: dF = 0  =>  ∂_α F_{βγ} + ∂_β F_{γα} + ∂_γ F_{αβ} = 0.
for (a, b, c) in [(0, 2, 3), (0, 3, 1), (0, 1, 2)]:
    cyc = (sp.diff(F(b, c), coords[a])
         + sp.diff(F(c, a), coords[b])
         + sp.diff(F(a, b), coords[c]))
    assert_zero(f"cyclic Bianchi ({a},{b},{c})", cyc)
# Specialize to E,B via F_{0i} = -E_i and F_{ij} = epsilon_{ijk} B_k; verify
# the three cyclic identities reduce to ∂_t B + curl(E) = 0 component-wise.
```

This is non-tautological (Schwarz symmetry of mixed partials on smooth `Function` arguments must hold — sympy will simplify it because the coords are declared real). The Mathematica block should be rewritten analogously, using `D[a_, b_] := D[a_, b]` on `A0[t,x,y,z]` etc., NOT by transliterating the sympy block line for line — pick different intermediate symbols (e.g., build an antisymmetric `4×4` array `Fmat` with `Fmat[[mu, nu]] := D[A[nu], coords[[mu]]] - D[A[mu], coords[[nu]]]` and loop over cyclic index triples).

If the maintainer prefers to keep the stage 004 script minimal and consider Maxwell-Faraday a stage-006 responsibility, the alternative is to delete the A2-A4 / M2 block entirely and remove the comment claim that this script tests "vector Bianchi signs." Either fix removes the tautology; do not leave the current block in place with the current comment.

**Verification:**

After fix, the verifier confirms:
- the new check at the named line range no longer contains the pattern `F23, F31, F12 = B1, B2, B3` followed by `sp.diff(F23, t) - sp.diff(B1, t)`-style direct substitution-comparison;
- if the replacement uses `F = dA`, the script still exits 0 (sympy/Mathematica should simplify the cyclic Bianchi to zero via Schwarz);
- if the replacement is deletion, the script docstring and the print line on (sympy 85, mathematica 6) no longer claim "vector Bianchi signs" coverage.

## Independent-derivation check (Mathematica)

The Mathematica script has different variable names from the sympy script (`decayingWindow`/`oddProbe` vs `W_ex`/`Q_ex`; `localizedProfile`/`profileArea`/`profileSelfMass`/`matchedWeight`/`overlapValue` vs `Z`/`Z_int`/`Z2_int`/`W_match`/`I_WZ`; `pointSourceCoupling`/`volumeReducedCoupling` vs `mu_proj_delta`/`mu_red`). The intermediate algebra is the same, which is expected for closed-form Gaussian integrals where the manipulation space is narrow. The check ordering M1-M6 exactly mirrors the sympy script's ordering (IBP, Bianchi-substitution, Z_int, Z²_int, I_WZ, μ ratio), and the M2 Faraday block has the same tautological structure as the sympy A2-A4 block (`twoForm23 = B1[...]` then `m2LeftOne - m2RightOne` with both sides expanding to the same expression). The naming is deliberately decorrelated but the structural choreography is one-to-one. This is a borderline `mathematica_transliteration` case; I do not file it as a separate finding because (a) the structural ordering was prescribed by the v1 directive's claim manifest M1-M6, so an independent script following the same spec would land in the same shape, and (b) the Gaussian-closed-form portion has limited room for genuine algorithmic divergence. The structural parallel of the Faraday block is best treated as a symptom of the same tautology issue (F1) rather than a separate transliteration finding.

## Engine cross-check

Both engines produce final residuals of exactly `0` for every check (sympy: implicit via `assert_zero` raising on nonzero residue and the script printing "STATUS: PASS"; mathematica output transcript shows `M1 residual = 0` through `M6 residual = 0` and `STATUS: PASS`). The Gaussian closed forms agree numerically: `Z_int = √π·λ`, `Z²_int = √(2π)·λ/2`, `I_WZ = √2/2`, `μ_proj_δ/μ_red = √2` in both engines. Engines agree at the level they claim to.

## Verdict justification

The substantive paper-side deliverables — the projection-by-parts identity and bundle source ordering — are correctly verified in both engines. The integrated IBP check on Gaussian decay profiles non-tautologically exercises the boundary-vanishing case the paper says drives the downstream bookkeeping. The Gaussian / matched-kernel block (A5-A8 / M3-M6) is mathematically sound and the √2 ratio is the correct value (`μ_proj_δ/μ_red = (W_match(0)/I_WZ)/(1/Z_int) = ((1/(√π·λ))/(√2/2))·√π·λ = 2/√2 = √2`).

The single new finding is the Faraday/Bianchi block: v1's audit waved this through as "partial" (a sign-convention spot-check); reading it against the paper's bundle role rule (ii) and against the substitution structure of the code reveals that it is purely a renaming exercise with no math content. Both engines have the same defect. Verdict is `findings` (not `stop_cold`): F1 is mechanically fixable and does not propagate to a downstream-quoted constant.

Paper alignment is `aligned`: every paper-side deliverable in the stage 004 card has a script-side check (the IBP and inventory checks cover the stated Output; the Faraday block notionally covers rule ii but the cover is empty per F1; the Gaussian/√2 extras are bundle-README scaffolding consistent with the script's docstring framing). The notes file the script docstring names (`step_01_projected_maxwell_readme.md`) is not committed; if it ever lands and disagrees with the script's scope, this verdict should be revisited.

## Self-test notes

I checked: (1) the IBP integrand `W·Q_w + W_w·Q` for `W = exp(-w²/λ²)`, `Q = w·exp(-w²/λ²)` equals `d/dw(WQ) = d/dw(w·exp(-2w²/λ²))`, whose integral over (-∞,∞) is zero by the boundary vanishing of `w·exp(-2w²/λ²)`; this confirms A1/M1 are non-tautological and pass for the right reason; (2) the F1 fix proposal (cyclic Bianchi from `F = dA`) genuinely exercises Schwarz symmetry of mixed partials, which sympy/Mathematica enforce only when the function arguments are real coords — both scripts already declare `t, x, y, z` real, so the proposed check will simplify to zero non-trivially; (3) the √2 ratio is consistent with the paper's matched-channel rule (iii) but is properly anchored in stages 007-009 by the appendix table, so leaving A5-A8/M3-M6 as bundle scaffolding does not introduce a new paper_misalignment.
