---
unit_id: 101
batch: IV.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage101_natural_source_map_reduction.md"]
  paper_appendix: present
---

# Audit unit 101 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_101.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage101_natural_source_map_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows around line 1236 and the upstream "Retarded factorization" subsection, lines 260-343)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.txt`

## What the paper claims

The stage card (stage_101.tex line 13) says the computation "isolates the reduced product `\widehat m_0^{\,2}\chi_Q N_Q=1` and the canonical condition `\chi_Q=1`." The bottom-line quote (lines 15-17) states: "On the natural point-particle source map, `\widehat m_0\to1`, so the last obstruction is purely `\chi_Q`." The notes elaborate: at `mhat0 -> 1` the Stage-83 factorization becomes `N_Q = 1/chi_Q`; on the canonical branch `chi_Q = 1` gives `N_Q = 1`; with `Delta_Q := chi_Q - 1`, the small-`Delta_Q` expansion gives `N_Q - 1 = -Delta_Q + O(Delta_Q^2)`. The stage card also lists three `Checks`: (1) `mhat0^2 * chi_Q * N_Q` keeps source/conservative/outgoing factors separate; (2) higher odd terms begin beyond the point-particle 2.5PN coefficient; (3) outgoing `l=2` DtN fingerprint against the normalized `z = omega a / c_s` expansion.

## What the script claims to verify

Both scripts solve `mhat0^2 * chiQ * NQ = 1` for `NQ`, then substitute `mhat0 -> 1` (yielding `1/chiQ`), then `mhat0 -> 1, chiQ -> 1` (yielding `1`), and print a series of `1/(1+DeltaQ) - 1` around `DeltaQ = 0`. The Mathematica script wraps the substitutions in `expectZero` calls (which print and `Exit[1]` on nonzero). The SymPy script issues only `print` statements and contains no `assert` calls at all. Neither script tests checks (2) or (3) from the stage card.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| `mhat0^2 chiQ NQ = 1` is the starting factorized constraint (input from Stage 83) | both scripts encode and solve it (sympy line 8; wl line 33) | match (input only) |
| `mhat0 -> 1` reduces to `NQ = 1/chiQ` | sympy line 10 (print); wl line 40 (expectZero) | partial (mathematica check is tautological; see F1) |
| Canonical compact outgoing branch `chiQ = 1` gives `NQ = 1` | sympy line 11; wl line 41 | partial (mathematica check is tautological; see F1) |
| `Delta_Q := chiQ - 1`, then `N_Q - 1 = -Delta_Q/(1+Delta_Q)` exactly | sympy lines 13-15; wl lines 43-45 | partial (series printed but never compared to `-Delta_Q + Delta_Q^2`; see F3) |
| `N_Q - 1 = -Delta_Q + O(Delta_Q^2)` (small-`Delta_Q` linearization) | series_delta computed and printed; never asserted | partial / no assertion |
| Check (1) — product keeps source/conservative/outgoing factors separate | `mhat0`, `chiQ`, `NQ` declared as independent symbols | implicit match |
| Check (2) — higher odd terms begin beyond point-particle 2.5PN coefficient (appendix line 295: "`O(omega^7)` invisible to point-particle 2.5PN") | not exercised in either script | missing (see F4) |
| Check (3) — outgoing `l=2` DtN fingerprint vs `z = omega a/c_s` expansion (appendix lines 305-326) | not exercised in either script | missing (see F4) |
| SymPy script has at least one falsifiable assertion | none present — only `print` calls | missing (see F2) |

`paper_alignment: partial` — the load-bearing reduction (`mhat0->1 ⇒ NQ = 1/chiQ`; `chi_Q = 1 ⇒ N_Q = 1`) is present in both scripts at the symbolic level, but in a tautological form on the Mathematica side and with no assertion on the SymPy side. Two listed `Checks` are absent.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 8 | `sp.solve(sp.Eq(mhat0**2 * chiQ * NQ, 1), NQ)[0]` (no assert) | input factorization | n/a (no assertion) |
| A2 | sympy | 9-11 | `print(...)` of solved `NQ` and substitutions (no assert) | reduction to `1/chiQ` and `1` | no (no assertion) |
| A3 | sympy | 13-16 | `print(...)` of series and exact replacement (no assert) | small-`Delta_Q` linearization | no (no assertion) |
| A4 | sympy | 17 | `print(...)` of `sp.simplify(... - 1/(1+DeltaQ))` (no assert) | exact `1/(1+Delta_Q)` form | no (no assertion) |
| A5 | mathematica | 33-34 | `nQSol = First[Solve[mHat0^2*chiQ*nQ == 1, nQ]]; nQExact = ...` | input factorization | n/a (definition) |
| A6 | mathematica | 40 | `expectZero["point-particle natural branch reduction", (nQExact /. mHat0 -> 1) - 1/chiQ]` | reduction `mhat0->1 ⇒ NQ = 1/chiQ` | **no — tautological by construction** |
| A7 | mathematica | 41 | `expectZero["canonical compact outgoing branch gives NQ=1", (nQExact /. {mHat0 -> 1, chiQ -> 1}) - 1]` | `chi_Q=1 ⇒ N_Q=1` | **no — tautological by construction** |
| A8 | mathematica | 47 | `expectZero["exact replacement chiQ=1+DeltaQ", FullSimplify[(nQExact /. {mHat0 -> 1, chiQ -> 1 + deltaQ}) - 1/(1 + deltaQ), ...]]` | exact `N_Q = 1/(1+Delta_Q)` | **no — tautological by construction** |

The Mathematica `expectZero` calls all reduce to checking `f(x) - f(x) = 0` after the substitution: since `nQExact = 1/(chiQ * mHat0^2)`, applying `mHat0 -> 1` yields literally `1/chiQ`, so `(1/chiQ) - 1/chiQ = 0` is guaranteed by construction, not by the physics. Same for the other two checks.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:40`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:41`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:47`

**What's wrong:**
All three Mathematica `expectZero` checks reduce to identities by construction. After
```
nQExact = nQ /. First[Solve[mHat0^2*chiQ*nQ == 1, nQ]]
```
`nQExact === 1/(chiQ * mHat0^2)`. Substituting `mHat0 -> 1` then literally yields `1/chiQ`, so line 40's residual `(nQExact /. mHat0 -> 1) - 1/chiQ` is `1/chiQ - 1/chiQ ≡ 0` regardless of what physics the symbols are supposed to represent. Lines 41 and 47 follow the same pattern (`1 - 1 ≡ 0` and `1/(1+deltaQ) - 1/(1+deltaQ) ≡ 0`). None of these checks can fail even if Stage 83's factorization were the wrong equation.

**Why this matters:**
The Mathematica `expectZero` calls give the impression of an engine-independent verification of the reduction `NQ = 1/chiQ` on the natural source-map branch. They do not. They verify substitutions of a solved expression back into itself — i.e., that `Solve` and `/.` are mutual inverses, which is a built-in property of the CAS, not a check of the stage's claim. The non-tautological content of the stage (the LINK between `mhat0^2 chiQ NQ = 1` and the reduction to `1/chiQ`, the small-`Delta_Q` linearization to the right order) is not anchored anywhere.

**Required change:**
Replace each `expectZero` so that the residual being simplified is anchored to an independent statement of the same identity, not a re-substitution of `nQExact`:
- Line 40: `expectZero["point-particle natural branch reduction", (1/(chiQ*mHat0^2) /. mHat0 -> 1) - 1/chiQ]` is still tautological; instead test that the original factorization equation `mHat0^2 * chiQ * nQ == 1` with `mHat0 -> 1` and `nQ -> 1/chiQ` simplifies to `0`, anchoring to the input equation rather than the solved form. Concretely: `expectZero["point-particle natural branch reduction", (mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, nQ -> 1/chiQ}]`.
- Line 41: anchor to the input equation: `expectZero["canonical compact outgoing branch gives NQ=1", (mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, chiQ -> 1, nQ -> 1}]`.
- Line 47: anchor to the input equation: `expectZero["exact replacement chiQ=1+DeltaQ", (mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, chiQ -> 1 + deltaQ, nQ -> 1/(1 + deltaQ)}]`.

These remain symbolic identities (the equation `mhat0^2 chiQ NQ = 1` is the stage input from Stage 83), but they test the consistency `input ↔ proposed solution` rather than `solved expression ↔ substituted solved expression`. A typo in either side of the proposed reduction (e.g., `nQ -> chiQ` instead of `1/chiQ`) now fails, where before it could not.

**Verification:**
After the edit, lines 40/41/47 should still pass (the proposed solutions are correct), but if the SymPy script's proposed reduction at line 10 were mutated to `1/chiQ**2`, the corresponding Mathematica `expectZero` would fail. Verifier confirms the new `expectZero` arguments name `mHat0^2*chiQ*nQ - 1` (or equivalent) as the residual rather than `nQExact` substitutions.

### F2 — missing_verification_script (subtype: script_doesnt_cover_claim)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py:1-18`

**What's wrong:**
The SymPy script contains zero falsifiable assertions. Every line that could exercise a claim is wrapped in `print(...)`:
- Line 9 prints `sol_NQ` (no `assert sp.simplify(...) == 0`).
- Lines 10-11 print substitutions of `sol_NQ` (no `assert`).
- Lines 13-16 print the series and exact expression (no `assert`).
- Line 17 prints `sp.simplify(sol_NQ.subs({mhat0:1, chiQ:1+DeltaQ}) - 1/(1+DeltaQ))`. This is the only candidate for a substantive check, but it is not an `assert`. The script exits 0 regardless of the printed values.

A typo (e.g., `sol_NQ.subs({mhat0:2, ...})`) would print a wrong value and still pass the audit. The script's contract with the framework — "scripts must verify the stage's claim, not narrate it" — is unmet.

**Why this matters:**
The unit's manifest says `is_status_only_candidate: False`, so a script is required, not just a transcript. Without assertions, the "PASS" status of `output/.../stage101_..._sympy_audit.txt` is unrelated to whether the stage's reduction is correct. Engines-agree status is also bogus: both scripts trivially "agree" because neither side can fail.

**Required change:**
Add explicit `assert` statements anchored to non-tautological residuals. Mirror the structure recommended in F1 (anchor to the input equation `mhat0**2 * chiQ * NQ = 1` rather than to `sol_NQ`):
- After line 11, add:
  ```python
  assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, NQ: 1/chiQ})) == 0
  assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, chiQ: 1, NQ: 1})) == 0
  ```
- After line 16, add:
  ```python
  assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, chiQ: 1 + DeltaQ, NQ: 1/(1 + DeltaQ)})) == 0
  ```
- Also assert the small-`Delta_Q` linearization to the right order (this is the stage's stated form `N_Q - 1 = -Delta_Q + O(Delta_Q^2)`):
  ```python
  assert sp.expand(series_delta - (-DeltaQ + DeltaQ**2)) == 0
  ```
  The expected series is `-DeltaQ + DeltaQ**2` because `1/(1+DeltaQ) - 1 = -DeltaQ + DeltaQ^2 - DeltaQ^3 + ...`, and the script asks for order 3 (which `removeO()` truncates to degree 2). Matches the printed line in the transcript.

**Verification:**
After the edit, `redteam exec-sympy 101` shows the script still exits 0, but mutating any of the proposed solutions on the RHS (e.g., changing `1/chiQ` to `1/chiQ**2`) now makes the script exit 1. Verifier confirms `assert` keyword appears at least 4 times in the file.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py:14-16`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:44-46`

**What's wrong:**
The small-`Delta_Q` expansion `series_delta` (sympy line 14) and `seriesDelta` (wl line 44) are computed and printed, but never compared to the paper's stated form `-Delta_Q + O(Delta_Q^2)` (notes line 77; appendix line 1061 cites `N_Q - 1` in this form). The transcripts show `DeltaQ**2 - DeltaQ` and `-deltaQ + deltaQ^2` respectively — a human reading the transcript can verify, but the script cannot fail if the series turns out to be `Delta_Q + Delta_Q**2` (a sign flip).

This is distinct from F2 (which is about the SymPy script having NO assertions at all) and from F1 (which is about the Mathematica `expectZero`s being tautological): F3 is the substantive content of the small-`Delta_Q` linearization, which neither script currently anchors.

**Why this matters:**
The stage's downstream-use note (stage_101.tex line 27) says this card feeds Stages 107-113 and the 2.5PN/4PN outgoing bridge. The sign and order of the linearization in `Delta_Q` is the bridge between `chi_Q = 1` (exact closure) and the off-canonical case. A sign flip here would propagate to every downstream defect-transport step. The current scripts would not catch such a flip.

**Required change:**
Add the series-comparison assertions described in F2's final bullet (for SymPy) and the matching Mathematica check:
- SymPy: `assert sp.expand(series_delta - (-DeltaQ + DeltaQ**2)) == 0` after line 14.
- Mathematica: insert after line 44, `expectZero["small-DeltaQ series matches paper", Expand[seriesDelta - (-deltaQ + deltaQ^2)]];`.

**Verification:**
After the edit, the SymPy script and the Mathematica script both contain an explicit comparison of the truncated series to the literal expression `-DeltaQ + DeltaQ**2` (resp. `-deltaQ + deltaQ^2`). Verifier confirms both literals appear at the matching lines.

### F4 — paper_misalignment (subtype: script_missing_paper_claim)

**Severity:** medium
**Files:**
- Paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_101.tex:21-25` (the `\stagefield{Checks}` block)
- Paper: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:295` (Check 2 source) and `:305-326` (Check 3 source)
- Script: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py` (entire file)
- Script: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl` (entire file)

**What's wrong:**
The stage card lists three `Checks`:
> "(1) Check the product `\widehat m_0^{\,2}\chi_QN_Q` keeps source, conservative, and outgoing factors separate. (2) Check that higher odd terms begin beyond the point-particle 2.5PN coefficient. (3) Check the outgoing `l=2` DtN fingerprint against the normalized `z=\omega a/c_s` expansion."

Of these:
- Check (1) is implicitly honored by the script (three independent symbols), though no explicit assertion enforces separability.
- Check (2) is asserted in the appendix (line 295: "Higher odd denominator terms beginning at `O(\omega^7)` are invisible to the point-particle 2.5PN coefficient") but **never exercised in the script**. There is no `omega` symbol or `O(\omega^7)` expansion anywhere in either audit file.
- Check (3) is derived in the appendix (lines 305-326: `\Lambda_2^{\rm out}(z) = -3 + z^2/3 + z^4/9 + i z^5/9 + O(z^6)` and `\widehat Y_2^{\rm out}(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)`) but **never exercised in the script**. There is no `z` or Hankel-function expansion anywhere in either audit file.

These are checks the stage card itself promises but its scripts do not perform.

**Why this matters:**
The verification status of Stage 101 — including the load-bearing identification `chi_Q = 1` for the canonical compact branch — depends on Check (3) (the DtN fingerprint match). Without Check (3) being scripted in this unit, the chain of "I've verified `chi_Q = 1`" rests on Stage 80 (which the notes file cites: "Stage 80 already fixed the unique minimal passive/outgoing grouped-`P2` one-pole completion"). The user must decide whether Check 2 and 3 belong in stage 101 (so the script should be extended) or whether they have already been verified upstream (so the stage card should be trimmed to point to those upstream verifications instead of listing them as its own checks).

**Required change:**
This is a paper-vs-script disagreement that the auditor cannot mechanically resolve. See `## Resolve before fix_loop` in the directive.

**Verification:**
Pending user resolution. After the user decides the direction (extend script or trim paper), a follow-up directive will specify the concrete edit and verifier path.

## Independent-derivation check (Mathematica)

The Mathematica `.wl` is structurally a near line-for-line port of the SymPy `.py`: both solve `mhat0^2 * chiQ * NQ = 1` (sympy line 8 / wl line 33), substitute `mhat0 -> 1` (sympy line 10 / wl lines 37, 40), substitute `mhat0 -> 1, chiQ -> 1` (sympy line 11 / wl lines 38, 41), then form `1/(1 + DeltaQ) - 1` and series-expand to order 2/3 (sympy lines 13-14 / wl lines 43-44), then compare exact substitution to `1/(1 + DeltaQ)` (sympy line 17 / wl line 47). Variable names differ only in case (`mhat0`/`mHat0`, `DeltaQ`/`deltaQ`). The Mathematica script does add `expectZero` wrappers that the SymPy script lacks, but they wrap the same residuals.

Given the stage's content is a single algebraic substitution (solve a linear-in-`NQ` equation, then substitute), there is little room for genuinely independent derivation, so I do not raise `mathematica_transliteration` as a separate finding. The defect is captured by F1+F2: the structure of the verification is the same in both, and on both engines it is non-substantive.

## Engine cross-check

Both engines produce identical symbolic results:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `NQ` from factorization | `1/(chiQ*mhat0**2)` | `1/(chiQ*mHat0^2)` |
| `mhat0 -> 1` | `1/chiQ` | `chiQ^(-1)` |
| `mhat0 -> 1, chiQ -> 1` | `1` | `1` |
| `NQ - 1` in `DeltaQ` | `-DeltaQ/(DeltaQ + 1)` | `-1 + (1 + deltaQ)^(-1)` (same after simplification) |
| small-`DeltaQ` series | `DeltaQ**2 - DeltaQ` | `-deltaQ + deltaQ^2` (same after reordering) |
| exact replacement check | `0` | `0` |

`engines_agree: true`. Both outputs are fresh (sympy.txt mtime 2026-05-11 12:45 > sympy.py mtime 2026-04-01 12:39; wl.txt mtime 2026-05-11 13:05 > wl mtime 2026-05-11 11:56). `outputs_fresh: true`. (Minor cosmetic note, not filed as a finding: the Mathematica banner at wl line 26 says `STAGE 084 — NATURAL SOURCE-MAP REDUCTION` while the stage is unit 101; the terminal print at line 50 correctly says "Stage 101 Mathematica audit passed." This is a copy-paste leftover that doesn't affect math.)

## Verdict justification

Both scripts encode the stage's load-bearing algebra (`NQ = 1/(chiQ*mhat0^2)` and its limits) and agree numerically. They are, however, non-substantive: the Mathematica `expectZero` calls are tautological re-substitutions of the solved expression (F1), the SymPy script has no `assert` statements at all (F2), and the small-`Delta_Q` linearization that downstream stages depend on is never compared to the paper's claimed form `-Delta_Q + O(Delta_Q^2)` (F3). Additionally, two of the three `Checks` listed on the stage card (higher odd terms past 2.5PN; outgoing `l=2` DtN fingerprint vs `z` expansion) are not exercised by either script (F4). F1-F3 are mechanical fixes Codex can apply; F4 is a `paper_misalignment` and routes to the user.

Attacks I tried that did NOT yield findings: (a) sign of `Delta_Q` linearization — the printed series `DeltaQ**2 - DeltaQ` correctly equals `-DeltaQ + DeltaQ**2` after reordering, so no sign flip; (b) positivity assumptions on `chiQ`, `mhat0`, `NQ`, `DeltaQ` (sympy line 5; wl lines 29-31) — the paper has `Delta_Q = chiQ - 1` which can be negative if `chi_Q < 1`, BUT the script declares `deltaQ` positive (sympy line 5: `positive=True, real=True`; wl line 30 only Reals not positive); the SymPy declaration is over-restrictive but does not cause a wrong answer here because `1/(1+DeltaQ)` is regular for `DeltaQ > -1`, the relevant neighborhood — I considered raising `symbol_assumption_error` for the SymPy `positive=True` on `DeltaQ` but the printed series and exact form are unaffected, so this is at most cosmetic and I do not file it; (c) engine disagreement — none, outputs match; (d) DtN expansion sign convention vs appendix line 312 (`i z^5 / 9` vs `i z^5 / 27` in `\widehat Y_2^{\rm out}`) — not in scope because neither script exercises Check (3).

## Self-test notes

I checked (1) Variable independence — no `sp.diff` or `D[...]` calls; no derivative-zero traps. (2) Symmetry/parity — no integrals; no parity traps. (3) Trivial-case pre-check — mentally substituted the proposed F1/F2 anchored residuals: `(mhat0**2*chiQ*NQ - 1).subs({mhat0:1, NQ:1/chiQ})` simplifies to `chiQ * (1/chiQ) - 1 = 0`, correct; `(mhat0**2*chiQ*NQ - 1).subs({mhat0:1, chiQ:1+DeltaQ, NQ:1/(1+DeltaQ)})` simplifies to `(1+DeltaQ)/(1+DeltaQ) - 1 = 0`, correct; if a typo replaced `1/chiQ` with `1/chiQ**2`, the residual becomes `chiQ/chiQ**2 - 1 = 1/chiQ - 1`, which is non-zero — so the proposed assertion is genuinely falsifiable. (4) Path specifications — directives keep all script paths verbatim, no missing-script invention. (5) Paper round-trip — proposed fixes anchor to the equation `mhat0^2 chiQ NQ = 1` already stated in stage_101.tex line 13 and the notes lines 7, 65-73, so they do not introduce a fresh paper-side constant or claim.
