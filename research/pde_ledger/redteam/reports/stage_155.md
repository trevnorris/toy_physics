---
unit_id: 155
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage155_frozen_traction_fixedpoint.md]
  paper_appendix: present
---

# Audit unit 155 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_155.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage155_frozen_traction_fixedpoint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 1344 — bare `\input{stages/stage_155}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.txt`

## What the paper claims

The stage card states (quoted block, line 16):
> "At old traction, fixed point survives but moves off exact compensation: \(\mathcal R_{\rm fp}\approx0.282714\)."

The notes elaborate the full deliverable set: solve the exact nonlinear fixed-point map
\(\Sigma\mapsto e^{-\Phi_{\Sigma_0^*}[\Sigma]}/\int e^{-\Phi}\) at frozen canonical traction
\(\Sigma_0^*\approx1.80594111095636\), and report the fixed-point moments
\(\mathfrak g_{\rm fp}\approx0.693352419668063\),
\(\mathcal S_{\rm fp}\approx0.6216013167514007\),
\(\mathcal R_{\rm fp}\approx0.2827139049082381\),
\(\Pi_{\rm fp}\approx1.4885734438300713\), and verify that
\(\delta\mathcal R_{\rm fp}=\mathcal R_{\rm fp}-1/4\) matches the exact first-order
transport-law prediction \(\delta\mathcal R_{\rm pred}=
-\delta\mathfrak g/\sqrt{1+\mathfrak r_{F1}^2}+(\delta\mathfrak g)^2/(1+\mathfrak r_{F1}^2)\)
with \(\delta\mathfrak g_{\rm fp}\approx-0.0646826592766000\). The paper card's bottom-line
claim is the R value; the notes elevate the moment quartet and the transport-law identity
to load-bearing deliverables.

## What the script claims to verify

The SymPy script (docstring lines 5–11) claims three checks: (1) solve the exact nonlinear
fixed-point map, (2) report g_fp, S_fp, R_fp, Pi_fp, (3) verify the exact shifted-R transport
law from "Stage 154". In practice only the transport-law consistency is enforced by an
assertion (`raise AssertionError` at line 101 if `|pred − direct| > 1e-10`); the four moments
are printed but never asserted against any paper-stated value. The Mathematica script
performs the same fixed-point solve and additionally asserts five numeric targets via
`expectApprox` (lines 101–105): g_fp, S_fp, R_fp, Pi_fp against paper-stated values, plus the
transport-law consistency.

## Paper ↔ script cross-check

| Paper-side deliverable | SymPy coverage | Mathematica coverage |
|---|---|---|
| \(\mathcal R_{\rm fp}\approx0.2827139049082381\) (card + notes §2) | partial — printed only, no assertion | match — `expectApprox` line 103 |
| \(\mathfrak g_{\rm fp}\approx0.693352419668063\) (notes §1) | partial — printed only | match — `expectApprox` line 101 |
| \(\mathcal S_{\rm fp}\approx0.6216013167514007\) (notes §1) | partial — printed only | match — `expectApprox` line 102 |
| \(\Pi_{\rm fp}\approx1.4885734438300713\) (notes §2) | partial — printed only | match — `expectApprox` line 104 |
| Transport-law identity \(\delta\mathcal R_{\rm pred}=\delta\mathcal R_{\rm direct}\) (notes §3) | match — `AssertionError` line 100–101 | match — `expectApprox` line 105 |

Dominant pattern: Mathematica fully anchors all five deliverables; SymPy anchors only one of
five. Setting `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 100–101 | `if abs(float(pred - direct)) > 1e-10: raise AssertionError(...)` | transport-law identity | partial — uses same-run g_fp/R_fp so it tests internal self-consistency of the algebraic law, not the absolute moment values |
| A2 | mathematica | 101 | `expectApprox["g_fp numeric check", gFp, 0.693352419668063, 10^-12]` | g_fp | yes |
| A3 | mathematica | 102 | `expectApprox["S_fp numeric check", sFp, 0.6216013167514007, 10^-12]` | S_fp | yes |
| A4 | mathematica | 103 | `expectApprox["R_fp numeric check", rFp, 0.2827139049082381, 10^-12]` | R_fp (card bottom line) | yes |
| A5 | mathematica | 104 | `expectApprox["Pi_fp numeric check", piFp, 1.4885734438300713, 10^-12]` | Pi_fp | yes |
| A6 | mathematica | 105 | `expectApprox["transport law consistency", deltaRPred, deltaRDirect, 10^-10]` | transport-law identity | yes |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:87-99`

**What's wrong:**
The SymPy script's docstring claims it "Report[s] g_fp, S_fp, R_fp, Pi_fp" (line 10) and the
paper notes elevate these four moments plus the transport-law identity to the stage's
load-bearing deliverables. However, the SymPy script only enforces one assertion: the
transport-law consistency check at lines 100–101. The four moments g_fp, S_fp, R_fp, Pi_fp
are printed (lines 89–92) but never compared to the paper-stated values. The single
`AssertionError` cannot detect any error that perturbs g_fp and R_fp in a way that preserves
the analytical relation between them — both quantities are read off the same fixed-point
solve, so a bug in the numerical integrator or the fixed-point iteration that shifted both
g_fp and R_fp consistently would leave the check passing. The Mathematica twin (lines
101–104) does anchor all four moments with `expectApprox`; the SymPy side is asymmetric.

Quoted comment claim (sympy line 10): "Report g_fp, S_fp, R_fp, Pi_fp."

Sympy bottom-line assertion (line 100–101):
```
if abs(float(pred - direct)) > 1e-10:
    raise AssertionError("Exact transport law failed numerically.")
```

Mathematica anchor (line 103): `expectApprox["R_fp numeric check", rFp, 0.2827139049082381, 10^-12];`

**Why this matters:**
The paper card's load-bearing number is \(\mathcal R_{\rm fp}\approx0.282714\). A numerical
audit that prints this number but does not assert it against the paper value leaves the
sympy side as a status-print rather than a verification. Any silent change to the operators
`Ts`/`Tq`, the grid resolution `N`, the initial seed, or the convergence tolerance that
shifted the converged moments would not trip the existing assertion as long as the algebraic
relation \(\delta\mathcal R_{\rm pred}=\delta\mathcal R_{\rm direct}\) survived.

**Required change:**
After the existing prints (sympy line 92, just before the `dg = sp.Float(...)` line), add
four explicit assertions anchoring the moment quartet to the paper-stated values, matching
the Mathematica twin:

```python
assert abs(g_fp - 0.693352419668063) < 1e-12, ("g_fp mismatch", g_fp)
assert abs(S_fp - 0.6216013167514007) < 1e-12, ("S_fp mismatch", S_fp)
assert abs(R_fp - 0.2827139049082381) < 1e-12, ("R_fp mismatch", R_fp)
assert abs(Pi_fp - 1.4885734438300713) < 1e-12, ("Pi_fp mismatch", Pi_fp)
```

These constants are the exact values both the notes (§1, §2) and the Mathematica
`expectApprox` calls already use, so no new value is introduced.

**Verification:**
Re-running the SymPy script must still exit 0 (the existing Mathematica output confirms
these constants match the current fixed-point solve to ≤2.22e-16). The new assertions appear
between line 92 and line 94, and a fresh transcript should show the four moments still
printed identically to the current `.txt` output.

### F2 — paper_misalignment

**Subtype:** target_mismatch (cosmetic / label-only)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl:26`

**What's wrong:**
Both scripts print a banner titled `"STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"`,
but the paper stage card (`stage_155.tex` line 1) declares this is `Stage 155`. The output
transcripts both also embed the wrong stage tag (`scripts/output/...txt` line 11 and
`mathematica/output/...txt` line 11).

Sympy line 80: `banner("STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT")`
Mathematica line 26: `banner["STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"];`

This is a non-load-bearing label, but it is a script-side cross-reference to the wrong unit
and suggests the file was copy-pasted from an earlier stage (138). The math itself is
unaffected, but a reviewer reading the transcript would be unable to confirm by header alone
that it belongs to stage 155.

Additionally, the SymPy docstring at line 11 cites "Stage 154" as the source of the
transport law, but the notes (`moving_throat_pde_stage155_frozen_traction_fixedpoint.md`
§3) attribute it to "Stage 239". This is a separate textual inconsistency between the
script's docstring and the notes; it does not affect the math.

**Why this matters:**
Tracker hygiene and provenance. The verifier and downstream auditors index on these banners
when scanning transcripts; mislabeled output complicates correlation between stage-155
artifacts and stage-138 artifacts.

**Required change:**
Sympy line 80: change `"STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"` to
`"STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"`.
Mathematica line 26: same change.

The docstring "Stage 154" vs. notes "Stage 239" cross-reference is a paper-side / notes-side
inconsistency about which earlier unit hosts the transport law; resolving it is the user's
call (the math is unchanged either way). Do not auto-edit; surface as a Resolve-before
question.

**Verification:**
After the patch, fresh sympy and mathematica transcripts both begin with
`STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT`.

## Independent-derivation check (Mathematica)

The Mathematica `.wl` is structurally a direct re-implementation of the SymPy
fixed-point iteration. Compare:

- Sympy line 32–38 (grid + cos/cosh weights) vs Mathematica line 33–42: identical grid
  construction, identical kappa, identical `cGrid`/`kqGrid` definitions.
- Sympy `Ts`/`Tq` (lines 43–51) vs Mathematica `tsOperator`/`tqOperator` (lines 46–56):
  identical `Accumulate`/`cumsum` choreography.
- Sympy `solve_fixed_point` (lines 70–78) seeds with `normalize(Pi_star * np.exp(-Pi_star * x))`;
  Mathematica `solveFixedPoint` (line 70) seeds with `normalize[piStar*Exp[-piStar*xGrid]]`.
  Same seed, same tolerance (1e-13), same maxIt (400).

The Mathematica version differs only by the `phShift = ph - Min[ph]` safeguard
(line 64) before exponentiation, a numerical-stability tweak that cancels under `normalize`.

This is a numerical fixed-point solve — independent symbolic re-derivation is not really
applicable here (both engines must run the same iteration to land on the same fixed point).
The cross-check value is in the agreement at ~1e-16 level (engine disagreement test below).
This is not flagged as `mathematica_transliteration` because no alternative independent
derivation exists for "iterate this nonlinear operator to convergence"; both engines have to
implement the iteration. The Mathematica side is, however, the only place where the
load-bearing moments are pinned to paper values — F1 still applies on the sympy side.

## Engine cross-check

Both engines converge to identical iteration counts (16) and identical residuals
(6.084e-14). Reported moments:

| Quantity | SymPy | Mathematica | |Δ| |
|---|---|---|---|
| g_fp | 0.693352419668063 | 0.6933524196680632 | ≲1.1e-16 |
| S_fp | 0.6216013167514007 | 0.621601316751401 | ≲2.2e-16 |
| R_fp | 0.2827139049082381 | 0.2827139049082381 | 0 |
| Pi_fp | 1.4885734438300713 | 1.488573443830071 | ≲2.2e-16 |
| ΔR_direct | 0.0327139049082381117017348515219 | 0.03271390490823811 | ≲1.7e-17 |
| ΔR_transport | 0.0327139049082376621515011616188 | 0.032713904908237605 | ≲5.7e-17 |
| Engine disagreement | — | — | n/a (agreement at FP precision) |

Engines agree to floating-point precision on every moment. No engine_disagreement finding.

## Verdict justification

The math holds: both engines independently converge to the same fixed-point moments at
1e-16 agreement, and both confirm the transport-law identity at the ~5e-13 residual level
asserted by the script. R_fp = 0.2827139… matches the paper card's bottom-line claim
\(\mathcal R_{\rm fp}\approx 0.282714\) and the notes' more precise value
0.2827139049082381 to all printed digits. The notes' moment quartet and the
\(\delta\mathcal R\)-transport identity are all numerically verified by the Mathematica
script. The SymPy script, however, only enforces the transport-law consistency and leaves
the four moments unanchored to paper values (F1 — `insufficient_verification`). The banner
in both scripts is mislabeled as "STAGE 138" (F2 — cosmetic `paper_misalignment`).
Verdict: `findings`, count 2, with no stop-cold flag — both findings are mechanical
script-side fixes (add four `assert` lines and correct the banner string) and do not
require user resolution beyond confirming the labeling convention.

## Self-test notes

Variable-independence trap: F1 adds plain numerical `assert abs(value - constant) < tol`
checks against values both produced by the same fixed-point solve and stated in the notes;
no derivatives involved, no `simplify` under assumptions, no zero-trivially traps. The
constants 0.693352419668063, 0.6216013167514007, 0.2827139049082381, 1.4885734438300713
are taken verbatim from the Mathematica script's existing `expectApprox` calls and from
notes §1/§2, so paper round-trip is preserved. Trivial-case pre-check: each `assert` would
pass at the script's current numerical output (the Mathematica twin's residuals are
≤2.22e-16 against these same constants), so adding them does not break the current `EXIT_CODE: 0`.
Path specifications: F2 names script lines explicitly (sympy line 80, mathematica line 26).
