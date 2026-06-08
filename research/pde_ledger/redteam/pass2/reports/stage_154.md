---
unit_id: 154
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T17:33:08Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage154_coevolving_core_mouth_map.md]
  paper_appendix: present
---

# Audit unit 154 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_154.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage154_coevolving_core_mouth_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 32 "Stages 154--163: co-evolving core--mouth branch, linear defect transport, D/N similarity, and the final off-family normal coordinate"; `\input{stages/stage_154}` at line 1342)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.txt`

## What the paper claims

Stage 154 replaces the "fixed compensated core + corrected mouth source" closure by a single self-consistent co-evolving core–mouth map. The card's designated audit target (stage_154.tex:7 "Its audit target is the verification output quoted below") is the single quoted body equation (stage_154.tex:16): the nonlinear fixed-point map `Σ ∝ e^{−Φ_Σ0[Σ]}` with `R[Σ]=(g[Σ]−r_F1)²/(1+r_F1²)`. The notes expand this into the deliverables: (D1) the exact Family-1 core overlap law `R[Σ] = (g−r_F1)²/(1+r_F1²)`, `r_F1 ≈ 1.77799353547498` (notes:36-42); (D2) the canonical compensation equivalence `g = g_* ⟺ R = 1/4` on the lower branch, where `g_* = r_F1 − ½√(1+r_F1²) ≈ 0.758035…` (notes:90-107); (D3) the exact first-order defect transport `R = 1/4 − δg/√(1+r²) + δg²/(1+r²)`, hence `δR = −δg/√(1+r²) + O(δg²)`, with `√(1+r_F1²) ≈ 2.039916913…`, `δR ≈ −0.490215 δg` (notes:114-153); and (D4) the local slope/bias identity `Π = Φ'_Σ0[Σ](0) = Σ0[1 − R·S]`, with the imported closure `Σ0 = (20/9) T̂_m²` from Stage 140 (notes:158-188). The map definitions — the potential `Φ_Σ0[Σ] = Σ0[T_s − R·T_q]` and the fixed-point equation `Σ = e^{−Φ}/∫e^{−Φ}` (notes:59-81) — are postulates/definitions, not identities. The card `Checks` list (deviations about the renormalized point; even-preservation before reading the odd defect; tangent motion → `δ⊥=0`) reads as a block-level (154–163) narrative template; the `δ⊥=0` item in particular corresponds to the block's "final off-family normal coordinate" (a later stage, per appendix:32), not to stage 154 itself.

## What the script claims to verify

The SymPy docstring (py:7-12) enumerates four checks: (1) the exact Family-1 core ratio law `R(g)`; (2) the compensation equivalence `g=g_* ⟺ R=1/4`; (3) the exact defect-transport expansion of `R` around the lower compensated branch; (4) the linearized slope identity `dPi` in terms of `dSigma0, dS, dR`. Three load-bearing `expect_zero`/`expectZero` assertions implement these: `R(g_*) − 1/4 == 0` (py:36 / wl:35), the shifted-`R` polynomial equals `1/4 − dg/√(1+r²) + dg²/(1+r²)` (py:38-40 / wl:37-39), and the linearized `dPi` equals `(1−R*S*)dΣ0 − Σ0(R*dS + S*dR)` (py:45-57 / wl:46-53). The `R(g)` form itself is printed for display, not asserted (the `R(g_*)=1/4` and shift checks exercise it). All checks are fully symbolic in a free real symbol `r`; the script never pins the numeric `r_F1`, so it verifies the identities for ALL `r`, which is strictly stronger than checking the one numeric value.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: core ratio law `R[Σ] = (g−r_F1)²/(1+r_F1²)` (notes:36-42; card:16) | `R = (g-r)**2/(1+r**2)` printed py:32-33 / wl:31-32, and exercised by the `R(g_*)` and shift checks | match |
| D2: `g=g_* ⟺ R=1/4` on lower branch, `g_*=r−½√(1+r²)` (notes:90-107) | py:36 / wl:35 `R(g_*) − 1/4 == 0` with `g_star = r − √(1+r²)/2` | match |
| D3: defect transport `R = 1/4 − δg/√(1+r²) + δg²/(1+r²)` (notes:121-141) | py:38-40 / wl:37-39 shifted-`R` minus the boxed RHS `== 0` | match |
| D4 (base): slope identity `Π = Φ'(0) = Σ0[1 − R·S]` (notes:161-167) | only its LINEARIZED form is tested (py:45-57); base `Φ'(0)=Σ0(1−R·S)` not derived from the actual `Φ`/operators | partial — see Verdict; card-designated target is D1, script self-scopes to the linearized `dPi`, notes box is honored at first order |
| D4 (linearized): `δΠ = (1−R*S*)δΣ0 − Σ0(R*δS + S*δR)` (notes:170-181) | py:57 / wl:53 `dPi − dPi_expected == 0` | match |
| Map definitions: `Φ = Σ0[T_s − R·T_q]`, fixed-point `Σ=e^{−Φ}/∫e^{−Φ}` (notes:59-81) | none (these are definitions/postulates, not identities) | n/a — not checkable assertions |
| `Σ0 = (20/9)T̂_m²` (notes:185) | not in script (imported from Stage 140) | n/a — upstream carry-forward |
| numeric `r_F1`, `g_*`, `√(1+r²)`, `δR≈−0.490215` (notes) | not emitted; script keeps `r` free symbolic | n/a — not script outputs (symbolic identities subsume them) |

The card's explicit audit target (the `R`-law) and the script's four self-declared checks (D1, D2, D3, D4-linearized) are all faithfully exercised. The base slope identity `Φ'(0)=Σ0(1−R·S)` is the one notes deliverable the script abstracts (it verifies the linearized product-rule consequence rather than re-deriving `Φ'(0)` from the integral operators); this is discussed in Verdict and is not raised as a hard finding (the card-designated deliverable is met and the script honestly scopes itself to the linearized form).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `expect_zero(R(g_*) − 1/4)` | D2 (compensation ⟺ R=1/4) | yes — `g_*` independently constructed, residual non-trivially → 0 |
| A2 | sympy | 40 | `expect_zero(R(g_*+dg) − [1/4 − dg/√(1+r²) + dg²/(1+r²)])` | D3 (defect transport) | yes — RHS written independently; falsifiable on any coefficient error |
| A3 | sympy | 57 | `expect_zero(dPi − dPi_expected)` | D4-linearized (slope) | yes — `dPi_expected` written independently; falsifiable on a sign |
| B1 | mathematica | 35 | `expectZero[(R /. g→g_*) − 1/4]` | D2 | yes |
| B2 | mathematica | 39 | `expectZero[Series(R)/.g→g_*+dg − boxedRHS]` | D3 | yes |
| B3 | mathematica | 53 | `expectZero[dPi − dPiExpected]` | D4-linearized | yes |

All three load-bearing assertions (per engine) trace to a stage-154 deliverable and are non-tautological: each compares an engine-computed expression (`R(g_*)`, the polynomial shift, the syntactic/series linearization) against an independently authored target. None is an `x==x` round-trip; each fails on a wrong sign/coefficient.

## Findings

None.

### Examined and dismissed (NOT findings)

- **Dead code, not a defect.** py:46 computes `Pi_lin = sp.expand(sp.series(Pi, dSigma0, 0, 2).removeO())` but py:48-53 immediately OVERWRITES `Pi_lin` with a `subs`-based linearization, discarding the series result. Harmless noise (the live linearization is correct, see below); not a finding category.
- **`subs`-based linearization is correct.** Expanding `Π=(Σ0+δΣ0)(1−(R*+δR)(S*+δS))` yields `Σ0(1−R*S*) − Σ0R*δS − Σ0S*δR − Σ0δRδS + δΣ0(1−R*S*) − δΣ0R*δS − δΣ0S*δR − δΣ0δRδS`. The four `subs` keys (`dSigma0*dR`, `dSigma0*dS`, `dR*dS`, `dSigma0*dR*dS`) zero out all and only the second-order/cross terms (SymPy `subs` matches `dR*dS` inside `Sigma0*dR*dS`), leaving exactly `dPi_expected = (1−R*S*)δΣ0 − Σ0(R*δS + S*δR)`. The `.wl` uses an `epsLin` order-tracking `Series[..., {epsLin,0,1}]` (wl:47-49) — a genuinely different linearization mechanism — and lands on the same target.
- **Single-branch is correct, not `missing_branch`.** `R(g)=1/4` has two roots `g = r ± ½√(1+r²)`; the script tests only the lower `g_* = r − ½√(1+r²)`. The notes explicitly scope D2/D3 to the "lower compensated branch" / "lower compensation moment" (notes:91, 110, 153), so the upper branch is out of scope by design.
- **No hardcoded numerics.** The script keeps `r` a free symbol and verifies the identities for all `r`; it never pins `r_F1 ≈ 1.778`. So there is no literal-against-itself check and nothing to mis-reconcile numerically.

## Independent-derivation check (Mathematica)

Borderline-port, recorded as PARTIAL, NOT raised as a finding — consistent with the standing pass-2 disposition for stages whose only operations are substitution/series/simplify on given closed forms (cf. stage 150).

Three corresponding sections:
1. The encoded laws — py:32,35 `R = (g - r)**2 / (1 + r**2)` … `g_star = r - sp.sqrt(1 + r**2) / 2` vs wl:31,34 `rFun = (g - r)^2/(1 + r^2)` … `gStar = r - Sqrt[1 + r^2]/2`. These are the paper's boxed premises (notes:39, 100); BOTH engines must encode them — re-typing a physical premise is not transliteration of a derivation.
2. The shifted-`R` check — py:38-40 uses exact polynomial `R.subs(g, g_star+dg)` then `expand`; wl:37-39 uses `Normal[Series[rFun /. g→gStar+dg, {dg,0,2}]]`. Different mechanism (exact expansion vs truncated series; they coincide only because the object is degree-2).
3. The linearization — py:48-53 uses syntactic `subs` to drop cross terms; wl:46-49 uses `epsLin` order-flagging + `Series[..., {epsLin,0,1}]`. Genuinely different algorithmic routes to the same `dPi`.

So the one non-premise algorithmic step (the linearization) is performed by two distinct methods, and the shifted-`R` step by two distinct methods. The shared structure is confined to encoding the paper's premise laws. This is not a `mathematica_transliteration` defect.

## Engine cross-check

Both engines pass all three checks and agree:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| `R(g)` printed form | `(g - r)**2/(r**2 + 1)` (line 5) | `(g - r)^2/(1 + r^2)` (line 5) |
| `R(g_*) − 1/4` | `0` (line 6) | `0` + PASS (lines 6-7) |
| shifted-`R` residual | `0` (line 7) | `0` + PASS (lines 8-9) |
| `dPi` identity | `0` (line 12) | `0` + PASS (lines 14-15) |

All three residuals are `0` in both engines; no `engine_disagreement`. Carry-forward formula blocks are byte-identical across the two transcripts (sympy lines 14-18, math lines 17-21).

## Verdict justification

Clean. I read the card, the notes, and the appendix row before the scripts. The card's designated audit target — the core overlap law `R[Σ]=(g−r_F1)²/(1+r_F1²)` and its compensation consequence `R=1/4` — is faithfully and non-tautologically verified, as are the defect-transport expansion and the linearized slope identity. Attacks tried and failed: (a) the `R(g_*)=1/4` check is non-tautological — `g_*` is an independently written closed form and `(g_*−r)² = (1+r²)/4` only because of the specific `½√(1+r²)` offset; a wrong factor would not vanish; (b) the shifted-`R` and `dPi` targets are authored independently of the engine-computed LHS and fail on any sign/coefficient error; (c) transliteration is borderline but the two engines use distinct mechanisms for both non-premise steps (exact-expand vs series, syntactic-subs vs epsLin-series), so it is PARTIAL not a defect; (d) no stale `168π²`/`168%` label and no `4107`/`100π²` Family-1 radius constant appears anywhere in the stage-154 sources (grep clean); (e) outputs are fresh (sympy out 2026-05-28 11:30 > script 2026-05-27 23:11; math out 11:34 > script 11:33). The one soft gap — the script verifies the *linearized* slope identity rather than re-deriving the base `Π=Φ'(0)=Σ0(1−R·S)` from the integral operators — is not raised as a finding: the card's explicit audit target is the `R`-law (covered), the script's own docstring scopes item 4 to the linearized `dPi` (covered), and re-deriving `Φ'(0)` would require importing the full `T_s`/`T_q`/`S` operator machinery the stage deliberately abstracts (scope expansion, which the auditor must not prescribe). The notes' Section-4 box is honored at first order.

## Self-test notes

Variable-independence: the only differentiation-like operation is the `Series`/`subs` linearization in `dSigma0`/`dg`/`epsLin`, and `Pi`/`rFun` genuinely depend on those variables, so no expansion is identically the input — the `dPi` and shifted-`R` checks are substantive. Trivial-case/parity: substituting the closed form `g_* = r − ½√(1+r²)` into `R` gives `((1+r²)/4)/(1+r²)=1/4` exactly (verified by hand), and the degree-2 shift coefficients `−1/√(1+r²)` (linear) and `+1/(1+r²)` (quadratic) match the boxed RHS term-by-term. Paper round-trip: no fix is prescribed, so no new misalignment is introduced; the `r`-free symbolic checks subsume the notes' numeric `r_F1`/`g_*`/`√(1+r²)`/`δR` values without contradicting them. No directive written — zero findings.

## Value Reconciliation (pass-2 augmentation)

The stage-154 scripts are purely SYMBOLIC. They emit NO named numeric constants — every check is performed in a free real symbol `r`, which subsumes (and is strictly stronger than) the notes' numeric `r_F1`, `g_*`, `√(1+r²)`, and `δR` figures. The deliverable values are therefore the symbolic identities themselves; each reconciles to a boxed equation in the notes (the terse `.tex` card legitimately quotes only the `R`-law).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `R(g) = (g−r)²/(1+r²)` (core overlap law) | py:32 / wl:31; out (sympy) line 5, (math) line 5 | notes:36-42 (boxed); card stage_154.tex:16 | MATCH (notes + card) |
| `g = g_* ⟺ R = 1/4`, `g_* = r − ½√(1+r²)` | py:35-36 / wl:34-35; out (sympy) line 6, (math) lines 6-7 | notes:100,104-107 (boxed `R=1/4`) | MATCH (notes) |
| defect transport `R = 1/4 − δg/√(1+r²) + δg²/(1+r²)` | py:38-40 / wl:37-39; out (sympy) line 7, (math) lines 8-9 | notes:121-141 (two boxed forms) | MATCH (notes) |
| linearized slope `δΠ = (1−R*S*)δΣ0 − Σ0(R*δS + S*δR)` | py:45-57 / wl:46-53; out (sympy) line 12, (math) lines 14-15 | notes:170-181 (coupled conditions); base box notes:161-167 | MATCH (notes, at first order) |

The notes' base slope identity `Π = Φ'(0) = Σ0[1 − R·S]` (notes:161-167) is a STATED deliverable but lives correctly in the notes (and the script verifies its linearized consequence), so it is a MATCH in the docs, not a MISSING-DELIVERABLE — the gap is script-coverage (discussed in Verdict, intentionally not raised), not a doc-reflection defect. `Σ0 = (20/9)T̂_m²` (notes:185) is an upstream (Stage 140) carry-forward, not a stage-154 script output.

INTERNAL scaffolding (no finding): `R`/`rFun`, `g_star`/`gStar`, `R_shift`/`rShiftSeries`, `R_shift_expected`/`rShiftExpected`, `Pi`/`piExpr`, `Pi_lin`/`piLin`, `Pi0`/`pi0`, `dPi_expected`/`dPiExpected`, `epsLin` order-flag; pass/zero residual flags.

reconciliation: complete; 4 deliverable values checked, 0 misaligned
