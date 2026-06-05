---
unit_id: 012
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 012 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_012.tex`
- notes: `(none)` (no `notes/stages/moving_throat_pde_stage012_*.md` exists)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row at line 46; chapter context lines 1-12)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_012.tex`) is terse and qualitative. It states the stage "rewrites the projected bundle shifts in the primitive variables used by the Stage~021 reduced normal form" (lines 11-12). The primitive packet is built from the denominator `\Delta = A W - R^2`, the conservative numerator `Q`, the mixed transfer numerator `P`, and the `\omega^2` denominator slope `S_2` (lines 13-21); at the mouth, projected perturbations of these primitives determine the `z_n` and `n_n` slots used in Stage~010. The "Transport role" paragraph states "The audit checks both fixed-target and transported-target compatibility shifts" connecting projection-first mouth data to the same primitive variables that appear in the one-port Schur complement (lines 23-26). The `\stagefield{Output}` (lines 28-30) reads verbatim: "Stage~012 exports the primitive-to-bundle bridge used by the mouth-Taylor stages." The appendix row (line 46) describes "Primitive finite-throat bridge data induced by the projected electromagnetic sector," status `\StatusExactClosure{}`. The deliverables are therefore: (D1) closed forms for the induced first-order bundle corrections `z0,z2,z4,n0,n2,n4` from primitive slippages `q1,s1,h1,d1,p1,g1`; (D2) the fixed-target compatibility shift; (D3) the transported-target compatibility shift; (D4) the static `Xi1` (P1/P0) primitive contribution. The card states NO numeric constants — it is a purely algebraic-identity stage.

## What the script claims to verify

The SymPy docstring (lines 2-23) states the script perturbs the primitive one-port data `Q,S2,H,Delta,P,Gw` by slippages `q1,s1,h1,d1,p1,g1` and derives the induced bundle corrections `z0,z2,z4,n0,n2,n4` exactly at first order, which then feed the full-bundle formulas (D0..P4, Xi1, isotropic compatibility surface). The assertions: (a) each of the six bundle corrections matches a pinned closed form (lines 83-134); (b) each correction agrees between the series route `dlin` and the analytic Frechet/partial-derivative route (lines 141-146); (c) the static Xi1 closed form (lines 149-153); (d) the fixed-target compatibility shift `dCompat_direct` (lines 185-188); (e) the transported-target normalization surface, compatibility surface, and shift (lines 189-197); (f) z0 cancels from both compatibility differences in the q1/d1 channels while the bare normalization K-surface retains those channels (lines 203-221, including two `assert_nonzero` positive controls); (g) two `assert_nonzero` mutation/negative-control checks that flip the z4 sign (lines 222-229). The Mathematica `.wl` independently reproduces (a)-(g) under M1-M9 labels.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| D1 bundle corrections z0..n4 closed forms | py 83-134 (closed form) + 141-146 (Frechet cross-route); wl M2 156-167 | match |
| D2 fixed-target compatibility shift | py 185-188 (`dCompat_direct`); wl M5 212-215 | match |
| D3 transported-target compatibility shift | py 189-197; wl M6/M7 191-251 | match |
| D4 static Xi1 (P1/P0) primitive contribution | py 149-153; wl M3 176-181 | match |
| "primitive-to-bundle bridge exported" (Output) | the closed forms in §2 of output + z0 cancellation structure (py 203-229; wl M8/M9) | match |

No script-side assertion lacks a paper-side home; no paper deliverable lacks a script-side check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1-A6 | sympy | 83-134 | `factor(together(simplify(z*−closed)))==0` | D1 | yes |
| A7-A12 | sympy | 141-146 | `dlin route − Frechet route ==0` | D1 (cross-route) | yes |
| A13 | sympy | 150-153 | `Xi1_static − closed ==0` | D4 | yes |
| A14-A15 | sympy | 177-184 | Solve round-trips for K_one, K_norm | D2/D3 scaffolding | yes |
| A16 | sympy | 185-188 | `dCompat_direct − closed ==0` | D2 | yes |
| A17 | sympy | 189 | transported K surface closed form | D3 | yes |
| A18 | sympy | 190-193 | transported compatibility surface | D3 | yes |
| A19 | sympy | 194-197 | transported compatibility shift | D3 | yes |
| A20-A23 | sympy | 203-211 | no-z0-channel in q1/d1 (4 checks) | Output (cancellation) | yes |
| A24-A25 | sympy | 214-221 | K_norm retains q1/d1 channel (nonzero) | Output (positive ctrl) | yes |
| A26-A27 | sympy | 222-229 | z4-sign mutation must NOT vanish (negative ctrl) | D2/D3 (anti-tautology) | yes |
| M1-M9 | mathematica | 156-285 | mirror of A1-A27 via Series+partial routes, plus expectNonZero controls | D1-D4 + Output | yes |

## Findings

None.

### Attacks attempted and survived

1. **Tautology on closed-form checks (A1-A6).** The closed forms are pinned literals, so in isolation A1-A6 could be "write the answer, assert the answer." But the corrections `z0..n4` are computed two independent ways (`dlin` series, lines 73-74; `frechet` analytic partials, lines 138-139) and A7-A12 (lines 141-146) assert the two routes agree. A pinned literal that were wrong would have to be wrong identically in both the hand-typed closed form AND match both computed routes — not a tautology. Verified z0 by hand: `Z0=Q/Δ`, `dZ0 = q1/Δ − Q d1/Δ² = (Δq1−Qd1)/Δ²`, matches line 83 and output line 29. Verified Xi1: `n0/N0 = 2(Δp1−Pd1)/(ΔP) = 2p1/P − 2d1/Δ`, plus `z0/D0 = (Δq1−Qd1)/(D0Δ²)`, matches line 152 and output line 55. Holds.

2. **Anti-tautology controls present.** A26-A27 (lines 222-229) deliberately flip the z4 sign and assert the residual is nonzero; the Mathematica M9 residual (output lines 61-63) is `(6·z4_numerator·shapeS²)/(Δ⁴shapeT²)` — a genuine nonzero quantity equal to `6S²z4/T²`, confirming the compatibility shift truly carries the `+3S²z4/T²` term and the assertion is not vacuously satisfied. This is the correct sign-control design.

3. **z0-cancellation claim (the substantive Output claim).** A20-A23 assert the q1/d1 partials of `K_norm − K_one` equal those of `compat_direct` (which has no z0 dependence), proving z0 cancels from the compatibility difference; A24-A25 (positive controls) assert `K_norm` alone keeps q1/d1 channels (Mathematica residuals `ell/Delta` and `−ell(2P²+ΔpGoal Q)/(Δ³pGoal)`, output lines 57-59 — nonzero). The pair (cancels in difference, survives in bare surface) is a real, non-circular structural test. Holds.

4. **Variable-independence trap (self-test #1).** The `assert_nonzero` partials (`sp.diff(K_norm_p, q1)`, line 217; `diff(K_norm_p, d1)`, line 220) — checked that `K_norm_p` genuinely depends on q1 and d1 through `z0` and `n0`. `K_norm_p` solves `(N0base+ell n0)/(K−B0−Z0slot−ell z0) = Ptarget` ⇒ `K = B0+Z0slot+ell z0 + (N0base+ell n0)/Ptarget`, which depends on q1 (via z0) and d1 (via z0 and n0). Nonzero is correct; Mathematica residuals confirm. No identically-zero-derivative trap.

5. **Symbol-domain check (self-test, symbol_assumption_error).** SymPy declares `Q,S2,H,Delta,P,Gw` and the K-surface symbols `nonzero=True` (lines 51, 148, 155); slippages `q1..g1` and `ell` unrestricted (lines 52-53). This is correct: the denominators are `Delta`, `Delta²..Delta⁵`, `D0`, `Ptarget`, `T`, so `Delta≠0`, `D0≠0`, `Ptarget≠0`, `T≠0` are exactly the non-vanishing assumptions needed; slippages legitimately may be any sign. Mathematica declares everything `Reals` plus `Delta!=0 && P!=0 && D0!=0 && pGoal!=0 && dTarget!=0 && shapeT!=0` (lines 33-39) — consistent set, no over-strong assumption masking a branch. No `assume(positive=...)` anywhere, so no positivity over-assumption risk.

6. **Parity/symmetry (self-test #2).** No integrals in this stage; all checks are rational-function algebraic identities. Trap N/A.

7. **Trivial-substitution pre-check (self-test #3).** Mentally set `q1=s1=h1=d1=p1=g1=0`: all `assert_zero` residuals (which are linear in the slippages) reduce to `0` trivially, while the `assert_nonzero` controls would then ALSO vanish — but those controls are nonzero for generic slippages, which is what `simplify==0?` tests (it asks whether the symbolic expression is identically zero, not zero at a point). The negative-control residuals are nonzero generic symbolic expressions (confirmed by Mathematica output lines 57-63). Correct.

8. **stale_output.** sympy.py mtime May 21 11:37 < sympy.txt May 25 17:15 (fresh); mathematica.wl May 25 17:14 < mathematica.txt May 25 17:15 (fresh). Both outputs newer than their scripts. No `stale_output` finding.

## Independent-derivation check (Mathematica)

The `.wl` is an independent second engine, not a transliteration. Evidence:
- Route divergence: SymPy derives corrections via a Taylor-series helper `dlin` (`series(...,ell,0,2).removeO()`, line 74) for the load-bearing values, cross-checked against `frechet` analytic partials (line 139). Mathematica derives via `Coefficient[Normal[Series[... , {ell,0,1}]], ell, 1]` (`z0Series`, lines 67-71) AND `D[form,var]*slip` partials (`z0Partial`, lines 98-103), asserts the two agree (M2 series-route vs partial-route, lines 156-167), then uses the partial route downstream (lines 169-174). SymPy uses the series route downstream. The two engines thus carry the result forward through *different* computed objects.
- Distinct symbol vocabulary: WL uses `kVar,bBare,zSlot,nSeed,dTarget,pGoal,shapeS,shapeT` (lines 28-31) where py uses `K,B0,Z0slot,N0base,D0target,Ptarget,S,T` (line 155) — not a mechanical rename of the same identifiers.
- The pinned closed forms (`z0Expected` etc., lines 135-154) match py's pinned forms, which is expected and correct: both engines independently confirm the same physics target; that is agreement, not echo. The choreography overlaps because the underlying derivation is the same first-order perturbation — acceptable for a second engine.

Conclusion: independent. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass all checks and produce identical closed forms. SymPy output §2 (lines 28-39) gives the six corrections; Mathematica M2 (output lines 9-32) confirms `residual = 0` for both the series and partial routes against the same expected forms. The Xi1 form (py output line 55) matches the WL M3 target (lines 179-181, residual 0, output line 33). The z0-channel positive controls match in spirit: py `assert_nonzero` and WL `expectNonZero` with explicit residuals `ell/Delta`, `−(ell(2P²+ΔpGoal Q))/(Δ³pGoal)` (output lines 57-59). The z4-mutation negative controls give the same `6S²z4/T²` residual in WL (output lines 61-63). Engines agree.

## Verdict justification

`clean`. I read the paper card, the appendix row and chapter context, and confirmed there are no source notes for this stage (terse algebraic-identity stage, status `\StatusExactClosure{}`). I built the deliverable model (D1 bundle corrections, D2 fixed-target shift, D3 transported-target shift, D4 static Xi1, plus the z0-cancellation Output claim) and mapped every one to a non-tautological, anchored script-side check present in BOTH engines. I attacked the closed-form pins (defeated by the dual-route cross-checks A7-A12), the cancellation claim (real positive+negative controls), the symbol domains (correct non-vanishing set, no over-strong positivity), the variable-independence trap (K_norm genuinely depends on q1,d1), and the negative controls (genuinely nonzero generic residuals). All survived. The paper states no numeric constants, so there is nothing to value-mismatch. Outputs are fresh. Both engines pass independently.

## Self-test notes

Checked: (#1) variable independence — `diff(K_norm_p, q1/d1)` are genuinely nonzero because K_norm carries z0(q1,d1) and n0(d1); not a zero-derivative trap. (#2) parity — N/A, no integrals, all rational-function identities. (#3) trivial substitution — `assert_zero` residuals vanish at zero slippage as expected; `assert_nonzero` controls are generic symbolic nonzeros (confirmed by WL residuals in output), so the simplify-is-zero test correctly fails for them. (#5 paper round-trip) — no fix prescribed, no new paper_misalignment introduced. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

The scripts emit purely **symbolic** closed-form results; there are no named numeric constants, benchmarks, or figures of merit. The stage card is deliberately terse (an algebraic-identity stage that "exports the primitive-to-bundle bridge") and there are no per-stage notes. Per the augmentation Guards, a terse `.tex` legitimately omits intermediate quantities, and there are no notes to carry them. None of the symbolic closed forms is reproduced verbatim in the `.tex` (the card names the primitive packet `\Delta, Q, P, S_2` and the slots `z_n, n_n` symbolically but does not display the derived closed forms). The card's `\stagefield{Output}` does not box any specific formula — it states the bridge is "exported," which the script's closed forms collectively constitute. Because every emitted value is a symbolic derivation rather than a stated deliverable constant, and the card's stated deliverable is the *existence/export* of the bridge (which the scripts supply), there are no MISMATCH or MISSING-DELIVERABLE items.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| z0 = (Δq1−Qd1)/Δ² | py:83 / py.txt:29; wl:135 / wl.txt:9-10 | card lines 13-21 (symbolic packet, form not displayed) | INTERNAL (symbolic derivation; card terse, no notes) |
| z2 closed form | py:84-87 / py.txt:31; wl:136-138 | not displayed in card | INTERNAL |
| z4 closed form | py:88-102 / py.txt:33; wl:139-142 | not displayed in card | INTERNAL |
| n0 = 2P(Δp1−Pd1)/Δ³ | py:103 / py.txt:35; wl:143 | not displayed in card | INTERNAL |
| n2 closed form | py:104-115 / py.txt:37; wl:144-147 | not displayed in card | INTERNAL |
| n4 closed form | py:116-134 / py.txt:39; wl:148-154 | not displayed in card | INTERNAL |
| Xi1^(primitive,static) = 2p1/P − 2d1/Δ + (Δq1−Qd1)/(D0Δ²) | py:149-153 / py.txt:55; wl:176-181 / wl.txt:33 | not displayed in card | INTERNAL |
| fixed-target compat shift = n0/Ptarget − 6S z2/T + 3S²z4/T² | py:185-188; wl:212-215 / wl.txt:39 | "checks both fixed-target ... compatibility shifts" card line 24 (named, not displayed) | INTERNAL |
| transported-target compat shift = −6S z2/T + 3S²z4/T² | py:194-197; wl:248-251 / wl.txt:47 | "and transported-target compatibility shifts" card line 24 (named, not displayed) | INTERNAL |
| transported K surface = B0+Z0slot+ell z0+D0target | py:189; wl:221-225 / wl.txt:41 | not displayed in card | INTERNAL |

Internal scaffolding (no prose expectation): pass/fail flags (`STATUS: PASS`), residual-zero check values, the `dlin`/`frechet` cross-route intermediates, the K_one/K_norm Solve round-trips, the z0-channel positive controls (`ell/Delta`, `−ell(2P²+ΔpGoal Q)/(Δ³pGoal)`), and the z4-sign mutation negative controls.

reconciliation: complete; 10 deliverable values checked, 0 misaligned. (All emitted results are symbolic derivations the terse card names but does not display; no numeric constants exist to mismatch, and the card's stated deliverable — exporting the bridge — is satisfied by the scripts.)
