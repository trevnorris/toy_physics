---
unit_id: 012
batch: I.1
auditor_model: claude-opus-4-7-1m
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

# Audit unit 012 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_012.tex`
- notes: `(none)` (per prompt; `notes/em_projected/step_NN_*.md` not committed for EM-projected stages 004-020)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 46 + `\input{stages/stage_012}` line 103)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.txt`

## What the paper claims

The stage card (Stage 012, "Projected Maxwell primitive bridge", Part I anchor MTDC-T4, status `StatusExactClosure`) identifies the primitive packet as built from the denominator `Delta = A W - R^2`, the conservative numerator `Q`, the mixed transfer numerator `P`, and the omega^2 denominator slope `S2`. Per the card, "at the mouth, projected perturbations of these primitives determine the `z_n` and `n_n` slots used in Stage 010." The card explicitly states the audit "checks both fixed-target and transported-target compatibility shifts." The Output line reads: "Stage 012 exports the primitive-to-bundle bridge used by the mouth-Taylor stages." The Part I appendix row 46 reinforces: "Primitive finite-throat bridge data induced by the projected electromagnetic sector." No source notes file is available for this EM-projected stage.

## What the script claims to verify

Both scripts perturb the primitive one-port data `(Q, S2, H, Delta, P, Gw)` by small slips `(q1, s1, h1, d1, p1, g1)` and derive the first-order induced bundle corrections `(z0, z2, z4, n0, n2, n4)` as both a series-coefficient extraction and a Frechet/partial-derivative sum, asserting the two routes agree with independently typed closed forms. They then build the K-surface algebra and verify (a) the static `Xi1 = n0/N0 + z0/D0` closed form, (b) the fixed-target linear compatibility shift `n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2`, (c) the transported-target linear compatibility shift `-6 S z2/T + 3 S^2 z4/T^2`, (d) that `z0` cancels from the `q1, d1` partials of `K_norm - K_one` (both fixed-target and transported-target), (e) that the normalization K surface alone retains `z0` contributions in `q1, d1`, and (f) two mutation negative controls in which `z4`'s sign is flipped and the residuals must be nonzero.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Primitive packet identification (Q, S2, P, Delta from `A W - R^2`) | sympy 51, 56-62; mathematica 41-48 | match (Delta entered as opaque symbol; the `A W - R^2` decomposition is upstream, not load-bearing at this bridge level) |
| Projected perturbations determine z_n, n_n slots | sympy 76-146 (series + Frechet against closed forms); mathematica 77-177 (Series-Coefficient + partial-derivative routes against closed forms) | match |
| Fixed-target compatibility shift checked | sympy 185-188 ("primitive compatibility shift from competing K surfaces"); mathematica M5 (222-225) | match |
| Transported-target compatibility shift checked | sympy 194-197 ("primitive transported-target compatibility shift"); mathematica M7 (254-261) | match |
| z_n/n_n exportable as Stage 010 inputs | derived forms exposed in transcript section 2 (sympy output lines 28-39); independent verification belongs to Stage 010 | match (consistent with bridge role) |
| Static Xi1 primitive prefactor (supporting context) | sympy 149-153; mathematica M3 (186-191) | extra-but-consistent (no paper disagreement; the card discusses bridge role generally) |

Dominant pattern: all paper-side deliverables are exercised; no script-side check contradicts or materially exceeds the paper card. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 83 | `assert_zero` z0 closed form vs series | bridge (z0) | yes |
| A2 | sympy | 84-87 | `assert_zero` z2 closed form vs series | bridge (z2) | yes |
| A3 | sympy | 88-102 | `assert_zero` z4 closed form vs series | bridge (z4) | yes |
| A4 | sympy | 103 | `assert_zero` n0 closed form vs series | bridge (n0) | yes |
| A5 | sympy | 104-115 | `assert_zero` n2 closed form vs series | bridge (n2) | yes |
| A6 | sympy | 116-134 | `assert_zero` n4 closed form vs series | bridge (n4) | yes |
| A7 | sympy | 141-146 | Frechet route reproduces each z/n | independent re-derivation of bridge | yes |
| A8 | sympy | 150-153 | static Xi1 closed form | Xi1 primitive prefactor support | yes |
| A9 | sympy | 177-184 | K-surface solve round-trips | K-surface algebra used by compatibility | yes |
| A10 | sympy | 185-188 | fixed-target compatibility shift closed form | fixed-target compatibility | yes |
| A11 | sympy | 189-193 | transported normalization K surface + compatibility surface | transported-target setup | yes |
| A12 | sympy | 194-197 | transported compatibility shift closed form | transported-target compatibility | yes |
| A13 | sympy | 203-211 | q1, d1 partial cancellation of z0 in both compat differences | z0-channel cancellation | yes |
| A14 | sympy | 214-221 | `assert_nonzero` normalization-K retains q1, d1 channel from z0 | dual to A13 | yes |
| A15 | sympy | 222-229 | mutation tests (z4 sign flip → nonzero) | negative-control | yes |
| M1 | math | 50-66 | `expectZero[Z0form - Q/Delta]` etc. | (none — algebraic identity by definition) | no (tautology) |
| M2 | math | 166-177 | series + partial routes vs closed forms (z0..n4) | bridge identification (both routes) | yes |
| M3 | math | 186-191 | static Xi1 | Xi1 prefactor support | yes |
| M4 | math | 201-209 | K solve round-trips | K-surface algebra | yes |
| M5 | math | 222-225 | fixed-target linear compatibility shift | fixed-target compatibility | yes |
| M6 | math | 232-240 | transported normalization K surface + round-trip | transported-target setup | yes |
| M7 | math | 254-261 | transported compatibility surface + linear shift | transported-target compatibility | yes |
| M8 | math | 263-286 | q1, d1 cancellation in compat differences; retention in normalization K | z0-channel cancellation | yes |
| M9 | math | 288-295 | mutation tests (z4 sign flip → nonzero) | negative-control | yes |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl:50-66`

**What's wrong:**
The M1 block defines each primitive one-port form on lines 41-48 and then on lines 50-66 asserts that each form equals its own definition.

- Line 41: `Z0form = Q/Delta;`
- Line 50: `expectZero["M1 Z0 primitive one-port", Z0form - Q/Delta];`

The residual `Z0form - Q/Delta` is identically zero by construction; the same pattern repeats for `Z2form`, `Z4form`, `N0form`, `N2form`, `N4form` on lines 51-66. The Mathematica output transcript confirms each residual is literally `0` (`mathematica/output/...txt` lines 2-13). These six PASS lines cannot fail no matter what the physics is.

**Why this matters:**
Six PASS lines in the Mathematica transcript convey false signal — they look like substantive verifications of the primitive one-port forms but verify nothing beyond `FullSimplify[0]==0`. The *correctness* of `Z0..N4` as the right "primitive one-port" representations is delegated entirely to the M2 series and partial routes (which take the M1 expressions as ground truth). The paper card is silent on the M1 forms' provenance; the SymPy docstring (line 55) calls them the "Stage-4 / Stage-5 primitive one-port formulas." A reader inspecting the PASS count is misled.

**Required change:**
Remove the six M1 `expectZero` calls at lines 50-66 (they verify nothing). Replace with a single `Print` group that emits the assumed forms and labels them as upstream carry-forwards:

```mathematica
Print["M1 primitive one-port forms (carried from Stage 4 / Stage 5):"];
Print["  Z0 = ", fmt[Z0form]];
Print["  Z2 = ", fmt[Z2form]];
Print["  Z4 = ", fmt[Z4form]];
Print["  N0 = ", fmt[N0form]];
Print["  N2 = ", fmt[N2form]];
Print["  N4 = ", fmt[N4form]];
```

The substantive M2..M9 checks all stand on their own and are unaffected by this trimming.

**Verification:**
Re-run mathematica audit; the output should no longer contain the six tautological `PASS: M1 ...` lines but should still PASS M2..M9 and exit 0. The final `STATUS: PASS` line is unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration of the `.py`. Evidence:

- Different helpers (`expectZero`/`expectNonZero` with `Module`+`FullSimplify[Together[...]]` vs SymPy's `assert_zero` using `sp.factor(sp.together(sp.simplify(...)))`).
- Different linearization mechanic: `Coefficient[Normal[Series[..., {ell, 0, 1}]], ell, 1]` (mathematica lines 79, 219, 251) extracts the first-order coefficient directly; SymPy's `(sp.series(expr.subs(subs), ell, 0, 2).removeO() - expr) / ell` (sympy line 74) subtracts the zeroth-order term and divides — different algebraic paths to the same coefficient.
- Different K-surface symbol naming: mathematica uses `kVar, bBare, zSlot, nSeed, dTarget, pGoal, shapeS, shapeT` (line 30, 38); SymPy uses `K, B0, Z0slot, N0base, D0target, Ptarget, S, T` (line 155).
- Different solve idiom: `kVar /. First[Solve[..., kVar]]` (mathematica 196, 199, 230) vs `sp.solve(sp.Eq(...), K)[0]` (sympy 156, 160).
- The Mathematica script provides both a series route (M2 "series closed form" lines 77-106 → 166-171) AND an explicit partial-derivative route (M2 "partial route" lines 108-143 → 172-177) for each z/n; both are written out independently rather than reusing a helper, consistent with fresh authorship rather than mechanical translation.

The choreography parallels the SymPy structure because the paper claim is the same, but the algebraic derivations are independent. The single exception is M1, which is its own block on the Mathematica side with no SymPy counterpart (and is flagged as F1).

## Engine cross-check

Both engines pass and produce equivalent symbolic results:

- z0 closed form: both `(Delta*q1 - Q*d1)/Delta^2` (sympy output line 29; mathematica M2 residual 0 at line 14).
- z2, z4, n0, n2, n4 closed forms: residuals 0 in both engines (sympy run completes through `STATUS: PASS`; mathematica M2 lines 14-37).
- Static Xi1: identical form via M3 (mathematica output line 38-39) and sympy section 4 (output lines 53-55).
- Fixed-target compatibility shift: both match `n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2`.
- Transported-target compatibility shift: both match `-6 S z2/T + 3 S^2 z4/T^2`.
- z0-channel cancellation in compat differences: both pass (sympy `for slip in (q1, d1)` block; mathematica M8 lines 54-57).
- Normalization K retains q1/d1 via z0: mathematica residual `ell/Delta` for q1 and `-(ell*(2*P^2 + Delta*pGoal*Q))/(Delta^3*pGoal)` for d1 (nonzero, PASS); SymPy `assert_nonzero` passes.
- Mutation tests (z4 sign flipped): both engines report nonzero residual (mathematica M9 residuals on lines 66-68; sympy `assert_nonzero` passes).

Engines agree.

## Verdict justification

The scripts faithfully verify what the paper card claims: they identify the same primitive packet `(Q, S2, P, Delta` with auxiliaries `H, Gw)`, derive the induced bundle corrections `z0..n4` from primitive slippages by two independent routes (series + Frechet/partials), and check both fixed-target and transported-target compatibility shifts — the two specific deliverables the card calls out explicitly. The substantive assertion ladder (sympy 83-229 / mathematica M2-M9) is non-tautological: the closed forms targeted by each `assert_zero`/`expectZero` are stated independently of how the LHS is computed, and the z0-channel cancellation logic plus the mutation negative controls give real attack-survival evidence.

The one finding is the M1 Mathematica block (six self-equating identities). It is `low` severity because the substantive verification is intact downstream; it is a real `tautological_check` and inflates the apparent PASS count. No `paper_misalignment`: the paper's deliverables map cleanly onto script-side assertions; the static Xi1 block is supplementary context (not flagged as a paper-side gap because the card itself describes the bridge role generally without enumerating every exported quantity).

Outputs are fresh (sympy output 11:52 > script 11:37; mathematica output 11:51 > script 11:39).

Attacks tried that failed: (i) checked whether `Delta`'s lack of decomposition into `A W - R^2` is load-bearing — it is not at this bridge level; the script treats `Delta` as opaque and the primitive perturbation `d1` covers any composite variation. (ii) Checked whether `assert_zero`'s use of `sp.factor(sp.together(sp.simplify(...)))` could hide branch errors — the closed forms are rational in symbolic constants and require no domain assumptions beyond `nonzero=True` (set on sympy line 51); no `sqrt` or domain-sensitive structure. (iii) Checked whether the `q1, d1` cancellation in `K_norm - K_one` is tautological given `z0` enters both K-surfaces with the same coefficient — verified the cancellation is non-trivial because `z0 = (Delta*q1 - Q*d1)/Delta^2` itself carries `q1, d1`, so the cancellation depends on `z0` entering symmetrically, which is a genuine algebraic property of the K-surface construction (not a typed-in identity). (iv) Checked whether the `assert_nonzero` mutation tests could pass for spurious reasons — they are sign-flipped versions of the actual identity, so a nonzero residual cleanly confirms the sign was substantive (mathematica's `M9` residual `(6 * (...) * shapeS^2)/(Delta^4 * shapeT^2)` exhibits the expected `6 S^2 z4/T^2` doubled magnitude). (v) Reread the paper card looking for a deliverable the scripts miss — the card mentions Stage 010 and Stage 021 by reference but does not require this stage to verify the downstream wiring (that belongs to Stage 010's audit), so no `script_missing_paper_claim`.

## Self-test notes

Walked through the proposed F1 fix mentally: replacing six tautological `expectZero` calls with `Print` statements does not introduce any new check, so no paper round-trip concern arises (the assumed forms are unchanged in semantics). Verified the M1 expressions are pure rational functions in symbolic constants, so `Z0form - Q/Delta` and analogues are guaranteed zero — F1 is real. Parity/branch traps not applicable (no integrals or square roots). Path specification names the existing `.wl` file in the `mathematica/` directory; Codex edits only that file, no new file creation needed.
