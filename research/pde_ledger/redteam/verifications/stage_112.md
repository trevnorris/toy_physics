---
unit_id: 112
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 112

Note on finding ordering: the auditor report numbered F1 = mathematica_transliteration and F2 = paper_misalignment. The directive swapped them (F1 = paper_misalignment / Cluster C, F2 = mathematica_transliteration) because Cluster C resolution was held first pending user direction. This verification follows the directive's ordering, which also matches the close-out summary from the orchestrator.

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script, Cluster C)

**Classification:** resolved

**What changed:**

Cluster C resolution along direction (a) — "internal unit number 'Stage 112' is canonical for scripts". The three identifying strings flagged in the directive are now updated:

- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:3` reads `Stage 112 SymPy audit.` (previously `Stage 95 SymPy audit.`).
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:54` reads `print('stage112: PASS')` (previously `print('stage95: PASS')`).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:26` reads `banner["STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW"];` (previously `STAGE 095`).
- The trailing Mathematica print at `.wl:83` already said `Stage 112 Mathematica audit passed.` and is unchanged, consistent with direction (a) in the directive.

In addition the SymPy provenance comment block at `.py:5-13` was refreshed; line 7 now reads `` `Lambda_out` is the exact canonical outgoing DtN branch from Stages 104/105. `` (was `Stage 087-088` per the orchestrator summary). Stage 104 is the canonical Lambda_out source after the renumber chain, so this aligns the provenance comment with the rest of the script's stage references and with the directive's "Stage 104" anchor at the top of the file. This is an additive cosmetic update beyond the three lines named in the directive, but it is internally consistent with direction (a) and does not touch any derivation, assertion, or assumption.

**Assessment:**

The three required label edits are in place and align across both engines and across all banner/intermediate/final strings:

- SymPy transcript shows `stage112: PASS` (`scripts/output/...sympy_audit.txt:6`).
- Mathematica transcript shows `STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW` in the banner (`mathematica/output/...mathematica_audit.txt:3`) and `Stage 112 Mathematica audit passed.` at the tail (`...:34`).

Both engines now consistently identify the unit as Stage 112. No algebraic content was touched by the label sweep. The provenance comment update is editorial, references stages that the rest of the script already cites, and is not load-bearing for any assertion. Direction (a) was the path the directive itemized most explicitly, so the closure matches the directive's offered resolutions.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**

A new independent Stage-92 linearized cross-check block was inserted in `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:55-66`, sitting between the four branch-identity `expectZero` calls (lines 50–53) and the `chiA = ...` definition (line 68). The block:

1. Extracts `a0Def` as the deviation of the constant-piece (`z^0` coefficient) of `lambdaHyb /. solB` from the Lambda_out constant `-3` (line 60), and `a5Def` as the imaginary-part deviation of the `z^5` coefficient from `1/9` (line 61).
2. Asserts `a0Def - 3*sigma == 0` and `a5Def + sigma*gamma == 0` (lines 62–63), matching the notes' branch-selection data `(a_0, a_5) = (3 sigma_W, -sigma_W gamma_W)`.
3. Solves the preservation condition `a0Def/3 + 9*a5Def == 0` for `gamma` (line 64), prints the result (line 65), and asserts `gammaFromLinear - 1/9 == 0` (line 66).

The existing chi_Q-based derivation (`chiA`/`chiB`, the chi_Q assertions at lines 74–76, and the scaled-identity check at lines 78–80) is left in place unchanged, so the two derivation paths now run as complementary engine paths to `gamma_W = 1/9` rather than as a transliterated echo.

**Assessment:**

The block uses a genuinely different algebraic route to land `gamma_W = 1/9`:

- The chi_Q path uses `(-L5/L0)/(1/27)` evaluated on `solB`, asserts the closed form `(1 - 9 sigma gamma)/(1 - sigma)`, substitutes `gamma -> 1/9`, and checks the residual is 1.
- The new block reads off `(a_0, a_5)` directly from `lambdaHyb /. solB - lambdaOut` (no chi_Q construction, no division-and-substitution choreography), asserts the notes' branch-selection identifications, and `Solve`s the linearized preservation constraint for `gamma`.

These are algebraically independent: a bug in the chi_Q normalization (`(-L5/L0)/(1/27)`) would not also break the direct coefficient read-off, and conversely a sign error in the `(a_0, a_5)` extraction would not also break the chi_Q closed form. The Mathematica transcript shows the three new lines `independent: a_0 - 3 sigma on solB = 0`, `independent: a_5 + sigma gamma on solB = 0`, and `independent: gamma_W from a_0/3 + 9 a_5 = 0 = 0` all PASS, with `gamma_W from linearized preservation = 1/9` printed in between (transcript lines 15–21). The assertions are non-tautological: `a0Def` and `a5Def` are computed before any substitution that would force them to 3*sigma or -sigma*gamma, and `gammaFromLinear` is the symbolic output of `Solve`, not a constant the script writes in.

The SymPy script is unchanged (the auditor's F1 in the original report only required an *independent* path in the Mathematica twin; the SymPy script's chi_Q-based derivation remains the primary engine A path). Engine independence is now restored.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

```
Lambda_hyb(z) = -I*gamma*sigma*z**5 - kappa**2*sigma*z**4 - kappa*sigma*z**2 + rho - sigma + I*z**5/9 + z**4/9 + z**2/3 - 3
canonical-even solutions = [{kappa: 0, rho: sigma}, {kappa: 1/3, rho: 4*sigma}]
chi_Q branch B = (9*gamma*sigma - 1)/(sigma - 1)
scaled identity on branch B = 0
stage112: PASS
```

**Mathematica:** exit=0. 11 PASS lines:

```
PASS: branch A rho - sigma
PASS: branch A kappa
PASS: branch B rho - 4 sigma
PASS: branch B kappa - 1/3
PASS: independent: a_0 - 3 sigma on solB
PASS: independent: a_5 + sigma gamma on solB
PASS: independent: gamma_W from a_0/3 + 9 a_5 = 0
PASS: chi_Q branch A - (1 - 9 sigma gamma)
PASS: chi_Q branch B - (1 - 9 sigma gamma)/(1 - sigma)
PASS: chi_Q branch B at gamma=1/9
PASS: scaled identity on branch B
```

The three new independent-route lines (a_0, a_5, gamma_W) appear exactly where the directive placed them — between the branch-identity block and the chi_Q closed-form block — and `gamma_W from linearized preservation = 1/9` is printed verbatim (transcript line 19).

**Output freshness:** both transcripts are newer than their scripts.

- SymPy script mtime 2026-05-27 15:15; SymPy output mtime 2026-05-27 15:18 (output newer).
- Mathematica script mtime 2026-05-27 15:16; Mathematica output mtime 2026-05-27 15:24 (output newer).

## Material-change assessment

`material_change`: false.

F1 is a string-only label sweep plus a provenance-comment refresh; no algebra changed. F2 adds an independent verification path that lands the same `gamma_W = 1/9` already verified by the chi_Q route — it expands the verification surface area but does not change any derived result. Downstream units that consume `gamma_W = 1/9`, the branch B identifications `(rho_R = 4 sigma_W, kappa_W = 1/3)`, the chi_Q^hyb closed form, or the collapse identity `Lambda_hyb = (1 - sigma_W) Lambda_out` will see exactly the same values as before. No downstream re-audit is required on derivation-content grounds.

## Side observations (non-blocking)

- The Mathematica block's `a0Def` line uses `(Coefficient[lambdaHyb /. solB, z, 0]) - (-3)` to peel off the Lambda_out constant. This is correct (Lambda_out's `z^0` coefficient is `-3`) but is implicit — a future reader without the notes open might miss why `-3` appears as a literal. Not a finding; just a readability note.
- The directive offered direction (a), (b), or (c) for F1 and the close-out chose (a). Direction (a) is reflected consistently across all four flagged strings (sympy:3, sympy:54, math:26, math:83), so the choice is internally consistent.

## Verdict justification

Both findings are resolved with no regressions. F1: the three flagged identifier strings (sympy:3 docstring, sympy:54 final print, math:26 banner) all now reference Stage 112, and the trailing math:83 print that already read "Stage 112" is unchanged, exactly per direction (a) of the directive; a small additive provenance-comment refresh on the SymPy side is internally consistent and non-load-bearing. F2: the Mathematica script now derives `gamma_W = 1/9` via the Stage-92 linearized `(a_0, a_5)` cross-check using a genuinely independent algebraic path (direct coefficient read-off + `Solve` on the linear preservation condition), and the three new `expectZero` lines (a_0, a_5, gamma_W) PASS alongside the seven pre-existing PASS lines for a total of 11 PASS lines in the Mathematica transcript. SymPy exits 0 with `stage112: PASS`. Outputs are fresher than scripts in both engines. No material change to any derived result downstream. Verdict: `verified`.
