---
unit_id: 122
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage122_mouth_source_compensation_test.md]
  paper_appendix: present
---

# Audit unit 122 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_122.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1278 `\input{stages/stage_122}`; no other narrative mentions this unit specifically)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py`
- mathematica: (missing — paper card explicitly says "Mathematica audit: none yet")
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The paper card (stage_122.tex line 15–17) gives the bottom-line claim in a quote block: "Equal-normalized mouth source misses the lower compensated branch but only by a modest traction renormalization." The card itself is terse and provides no closed-form expressions; the notes (which the audit prompt names as the authoritative source when the .tex is terse) make this concrete with five boxed identities:

1. The natural equal-normalized local mouth source defined by `g_s = g_m √K_s`, `g_q = g_m √K_q` gives `𝔤_nat = 1`.
2. The exact compensation family from Stage 221 is `𝔤_comp^±(𝔯) = 𝔯 ± (1/2)√(1+𝔯²)`, equivalent to `1+𝔯² = 4(𝔤−𝔯)²`.
3. Inserting the Stage 223 geometric value `𝔯_F1 = √(4107−168π²)/(10π)` gives the two exact compensated couplings `𝔤_±^F1 = (2√(4107−168π²) ± 37√3)/(20π)`, with numerical values `𝔤_−^F1 ≈ 0.758035078944663` and `𝔤_+^F1 ≈ 2.79795199200529`.
4. The compensation defect on the natural branch is `𝒞_nat = 1 + 𝔯_F1² − 4(1−𝔯_F1)² = (−12321 + 80π√(4107−168π²))/(168π²) ≈ 1.74016524722739`.
5. The traction renormalization ratios are `𝒯_m^(±)/𝒯_m^nat = 1/𝔤_±^F1`, numerically `≈ 1.31920016339112` and `≈ 0.357404273860789`.

The notes' "Result" block enumerates four deliverables: `𝔯_F1 ≈ 1.778`, `𝔤_nat = 1`, "does not lie exactly on the compensation surface", and "nearest compensated branch is only a modest traction renormalization away."

## What the script claims to verify

The script (`stage122_mouth_source_compensation_test_sympy_audit.py`) computes:
- `R = 37/20`, `rF = sqrt(12 R²/π² − 1)` (this simplifies to `sqrt(4107 − 100π²)/(10π)`, not `sqrt(4107 − 168π²)/(10π)` as the notes claim — see F1 below).
- `gminus = simplify(rF − sqrt(1+rF²)/2)`, `gplus = simplify(rF + sqrt(1+rF²)/2)`.
- `g_nat = 1` (set as integer, not asserted).
- Compensation defect `1 + rF² − 4(g_nat − rF)²` (printed only, no assertion).
- Distances `g_nat − gminus`, `gplus − g_nat` (printed only).
- Traction ratios `1/gminus`, `1/gplus` (printed only).
- Two `expect_zero` assertions: `gminus − (2√(4107−100π²) − 37√3)/(20π) == 0` and `gplus − (2√(4107−100π²) + 37√3)/(20π) == 0`.

Only the two algebraic-identity assertions are load-bearing. Everything else is `print` output, so the script's effective verified claim is "the compensated branch values derived from the geometric `rF` match the explicit closed form `(2√(4107−100π²) ± 37√3)/(20π)`."

## Paper ↔ script cross-check

| Paper-side deliverable (notes) | Script check | Status |
|---|---|---|
| `𝔯_F1 = √(4107−168π²)/(10π)` (notes §2 box) | Script defines `rF = sqrt(12·R²/π² − 1)` which equals `sqrt(4107 − 100π²)/(10π)` (algebraically different from notes) | mismatch — see F1 |
| `𝔤_nat = 1` from `g_s = g_m √K_s, g_q = g_m √K_q` (notes §1 box) | Script sets `g_nat = sp.Integer(1)` by fiat; no derivation, no assertion | partial — see F2 (the value is hardcoded; the notes' derivation from channel weights is not exercised) |
| `𝔤_±^F1 = (2√(4107−168π²) ± 37√3)/(20π)` (notes §2 box) | A1/A2 assert `gminus − (2√(4107−100π²) − 37√3)/(20π) == 0` and same for `gplus` | mismatch in the algebraic form (100 vs 168) — but the numerical value of the script matches the notes' numerical claim (`0.758035...`, `2.79795...`); the script's form is the correct one and the notes have a typo. See F1. |
| `𝒞_nat = (−12321 + 80π√(4107−168π²))/(168π²) ≈ 1.74016524722739` (notes §3 box) | Script prints `comp_def = (−12321 + 80π√(4107−100π²))/(100π²)`, no assertion | mismatch in algebraic form + missing assertion; see F1 + F2 |
| `𝒯_m^(±)/𝒯_m^nat = 1/𝔤_±^F1` (notes §5 box) | Script prints `T_ratio_minus = 1/gminus`, `T_ratio_plus = 1/gplus`; no assertion | missing — see F2 |
| Compensation family identity `1+𝔯² = 4(𝔤−𝔯)²` defining `𝔤_±` (notes §2) | Not asserted; never checked that `gminus`/`gplus` actually solve this quadratic | missing — see F2 |
| Numerical claim `𝔯_F1 ≈ 1.778` (notes Result block) | Not printed; numeric `rF` value not surfaced in the output | partial — informational only |

The script's two algebraic assertions cover only one of the notes' five boxed deliverables, and even that one is exercised with a different literal closed-form than the notes state. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 47 | `expect_zero("gminus exact form", gminus - gminus_exact)` with `gminus_exact = (2*sqrt(4107 - 100*pi**2) - 37*sqrt(3))/(20*pi)` | Closed form of `𝔤_−^F1` (notes §2 box) — but with `100π²` where notes have `168π²` | partial — the identity is non-tautological (SymPy must rationalize the nested `sqrt(1+rF²)` against `37·sqrt(3)/(10π)`), but the right-hand side disagrees algebraically with the notes |
| A2 | sympy | 48 | `expect_zero("gplus exact form", gplus - gplus_exact)` with `gplus_exact = (2*sqrt(4107 - 100*pi**2) + 37*sqrt(3))/(20*pi)` | Closed form of `𝔤_+^F1` (notes §2 box) — same `100` vs `168` mismatch | partial — same reasoning as A1 |

Neither assertion exercises (a) the definition of `𝔤_nat = 1` from equal-normalized channel weights, (b) the compensation-family quadratic `1+𝔯² = 4(𝔤−𝔯)²` itself, (c) the closed form of the defect `𝒞_nat`, (d) the traction-ratio identities `𝒯_m^(±)/𝒯_m^nat = 1/𝔤_±^F1`. These are the bulk of the notes' deliverables.

## Findings

(`paper_misalignment` first per the prompt ordering.)

### F1 — paper_misalignment

**Subtype:** `notes_contradicts_script` (the notes' algebraic form contradicts the script's algebraic form; the script's form is internally consistent and matches the numerics, so the notes appear to have a typo)

**Severity:** high

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:50` (boxed `𝔯_F1 = √(4107−168π²)/(10π)`)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:56` (boxed `𝔤_±^F1 = (2√(4107−168π²) ± 37√3)/(20π)`)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:88` (boxed `𝒞_nat = (−12321+80π√(4107−168π²))/(168π²)`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:18–22` (`R = 37/20`, `rF = sqrt(12*R**2/pi**2 - 1)`, which equals `sqrt(4107 − 100π²)/(10π)`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:45–46` (`gminus_exact = (2*sqrt(4107 - 100*pi**2) - 37*sqrt(3))/(20*pi)`)

**What's wrong:**

The notes consistently write the throat geometric ratio and derived quantities with `4107 − 168π²` and `168π²` as denominator. The script consistently uses `4107 − 100π²` and `100π²`. Numerically:

- Notes claim: `𝔯_F1 ≈ 1.778` (Result block). `sqrt(4107 − 168π²)/(10π) = sqrt(2448.90)/31.416 ≈ 49.487/31.416 ≈ 1.575`, which does NOT match `1.778`. `sqrt(4107 − 100π²)/(10π) = sqrt(3120.04)/31.416 ≈ 55.857/31.416 ≈ 1.778`, which matches.
- Notes claim: `𝔤_−^F1 ≈ 0.758035078944663`. The script's `gminus_exact = (2*sqrt(4107 − 100π²) − 37√3)/(20π) ≈ (111.71 − 64.09)/62.83 ≈ 0.758` matches. The notes' literal formula `(2*sqrt(4107 − 168π²) − 37√3)/(20π) ≈ (98.97 − 64.09)/62.83 ≈ 0.555` does not match its own quoted numeric.
- Notes claim: `𝒞_nat ≈ 1.74016524722739`. Script's printed defect `(−12321 + 80π√(4107−100π²))/(100π²) ≈ 1.74016` matches. The notes' literal formula `(−12321 + 80π√(4107−168π²))/(168π²) ≈ (−12321 + 12442.9)/1658.1 ≈ 0.074` does NOT match `1.74016`.

So the notes' algebraic forms in §2 and §3 contradict their own numerical values; the script's algebraic forms agree with all the numerical values in the notes. The script appears correct; the notes appear to have a typographic substitution of `168` for `100`.

Independent algebraic check: with `R = 37/20`, `12·R² = 12·1369/400 = 16428/400 = 4107/100`, so `12·R²/π² − 1 = 4107/(100π²) − 1 = (4107 − 100π²)/(100π²)`. Thus `rF² = (4107 − 100π²)/(100π²)` and `rF = √(4107 − 100π²)/(10π)`. The factor `100` (= `400/4`, from `(20)²/(20·20)·10` in the script's setup) is correct; `168` has no derivation path from `R = 37/20`.

**Why this matters:**

This is a documentation-vs-script algebraic disagreement on the load-bearing closed forms quoted by the stage. Downstream stages 125–139 are stated to depend on these values. If a downstream author copy-pastes the notes' formulas verbatim, they will compute numbers that do not match this stage's actual output.

**Required change:**

Routed to user (no Codex action). See `## Resolve before fix_loop` in the directive.

**Verification:**

After user resolves direction, either (a) the notes are corrected to `4107 − 100π²` and `100π²` and the numerical claims are unchanged, or (b) the script and its outputs are updated to use `168π²` and the numerical values in the notes are also updated to match (this would invalidate the present `𝔯_F1 ≈ 1.778`, `𝔤_−^F1 ≈ 0.758`, `𝒞_nat ≈ 1.740` numbers). The natural reading is (a): the script's algebra is self-consistent with the numerics; the notes have a typo.

### F2 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:28–42` (computed but not asserted: `g_nat`, `comp_def`, `delta_g_minus`, `delta_g_plus`, `T_ratio_minus`, `T_ratio_plus`)

**What's wrong:**

The script's two `expect_zero` assertions only cover the closed-form equivalence of `gminus`/`gplus` to a rewritten expression of the same two quantities. The notes enumerate five distinct boxed identities (`𝔤_nat = 1`, the compensation-family quadratic `1+𝔯² = 4(𝔤−𝔯)²` evaluated by `𝔤_±`, the defect `𝒞_nat`, the closed form of `𝔤_±^F1`, and the traction ratios), plus the result-block claim `1 ≠ 𝔤_±^F1`. Of these, only the closed form of `𝔤_±^F1` is checked. In particular:

- The compensation-family identity itself is never asserted: the script never checks `simplify(1 + rF**2 - 4*(gminus - rF)**2) == 0` (and the same for `gplus`). This is the *definition* of the compensation surface; without this check, the script does not actually verify that `gminus, gplus` lie on the surface.
- The defect closed-form `(−12321 + 80π√(4107−100π²))/(100π²)` is printed but not asserted equal to `1 + rF² − 4(1 − rF)²` algebraically (the print just shows whatever `simplify` chooses to render).
- The traction-ratio identities `T_m^(±)/T_m^nat = 1/𝔤_±^F1` are printed but never asserted.
- The notes' core result statement "`1 ≠ 𝔤_±^F1`" (i.e., the natural branch is not on the compensation surface) is never expressed as `expect_nonzero(comp_def)` or equivalent.

**Why this matters:**

Three of the four boxed claims and the headline "Result" block claim of the notes are unexercised. The script's PASS only confirms one algebraic rewrite; if a future edit silently broke the compensation-family quadratic (e.g., a sign flip on `4(g−r)²`), nothing in the script would catch it. The paper's Output line ("Equal-normalized mouth source misses the lower compensated branch") is not exercised as an inequality.

**Required change:**

Add the following assertions to the SymPy script, after the existing `expect_zero` calls. (Use the script's existing 100-vs-168 form, since F1 indicates the script is the correct side; the user resolution of F1 may later require these to be re-expressed, but the algebraic identities below are framing-independent.)

1. Compensation-family quadratic identities for `gminus` and `gplus`:
   - `expect_zero("compensation quadratic at gminus", 1 + rF**2 - 4*(gminus - rF)**2)`
   - `expect_zero("compensation quadratic at gplus",  1 + rF**2 - 4*(gplus  - rF)**2)`

2. Defect closed form (asserting the script's printed expression matches the derivation):
   - `defect_exact = (-12321 + 80*sp.pi*sp.sqrt(4107 - 100*sp.pi**2)) / (100*sp.pi**2)`
   - `expect_zero("defect closed form", comp_def - defect_exact)`

3. Non-trivial defect (the headline "Result" statement that natural ≠ compensated):
   - Add `expect_nonzero` helper that raises if `simplify(expr) == 0`, then `expect_nonzero("natural off compensation", comp_def)`.

4. Traction-ratio identities (one assertion per branch, asserting `gminus * T_ratio_minus == 1`):
   - `expect_zero("traction ratio (-) identity", gminus * T_ratio_minus - 1)`
   - `expect_zero("traction ratio (+) identity", gplus  * T_ratio_plus  - 1)`

**Verification:**

After Codex applies, the output `.txt` should show six new lines: `compensation quadratic at gminus = 0`, `compensation quadratic at gplus = 0`, `defect closed form = 0`, `natural off compensation = <nonzero>` (e.g., the explicit `(-12321 + 80*pi*sqrt(...))/(100*pi**2)`), `traction ratio (-) identity = 0`, `traction ratio (+) identity = 0`. Script must exit 0. The verifier runs `redteam exec-sympy 122` and confirms.

### F3 — missing_verification_script

**Subtype:** `missing_mathematica`

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_122.tex:11` — "SymPy audit: ... Mathematica audit: none yet."

**What's wrong:**

No `.wl` exists for this unit. The unit's manifest entry has `is_status_only_candidate: False`, which per the prompt would normally require both engines. However:

- The user's audit message states: "no Mathematica `.wl` is present for this stage; only a SymPy script. Do not file `missing_mathematica` unless an upstream-referenced result truly requires it."
- The paper card explicitly acknowledges the absence: "Mathematica audit: none yet."
- The upstream-referenced results (`𝔯_F1` from Stage 223, the compensation family quadratic from Stage 221) are themselves checked or claimed in those upstream units. This stage is essentially a numerical/algebraic evaluation step plugging the Stage 223 geometric value into the Stage 221 family.

Therefore I am filing this as a low-severity informational finding only — the absence is acknowledged by the paper card and no upstream Mathematica-only result is contradicted. Per the user's instruction I am **not** marking this as blocking. (The verifier and orchestrator should treat this as informational; the principal red-team issues are F1 and F2.)

**Why this matters:**

If the project's policy in a later sweep promotes "Mathematica audit: none yet" to "Mathematica audit required," this stage will need a second-engine derivation. For now, the absence is documented.

**Required change:**

None (per user instruction). If a Mathematica script is later added, it must independently derive `gminus`/`gplus` from the compensation quadratic and check (1)–(4) of F2's claim manifest below, not transliterate the SymPy algebra.

**Verification:**

n/a (informational).

## Independent-derivation check (Mathematica)

n/a — no Mathematica script. The paper card line 11 acknowledges "Mathematica audit: none yet." See F3.

## Engine cross-check

n/a — single-engine unit.

Output mtime: sympy script `Apr 18 16:42`, sympy output `May 11 12:45`. Output is newer than script; `outputs_fresh: true`. No `stale_output` finding.

## Verdict justification

The script's two assertions are mathematically correct and non-tautological (SymPy must rationalize `sqrt(1+rF²)` against `37·sqrt(3)/(10π)` to confirm them), and the script's algebraic forms match its own numerical values and match the notes' numerical values. However, two issues prevent a clean verdict:

1. The script's `rF` and `gminus_exact`/`gplus_exact` use `4107 − 100π²`, while the notes' boxed formulas use `4107 − 168π²`. The script is self-consistent and the script's forms reproduce the notes' numerical claims; the notes' formulas do not reproduce their own quoted numbers. This is `paper_misalignment` (subtype `notes_contradicts_script`) and requires user resolution before any Codex edit.

2. Even granting the script's algebraic frame, only one of the notes' five boxed deliverables is checked. The compensation-family identity itself (which defines `gminus, gplus`), the defect closed form, the traction ratios, and the headline "natural ≠ compensated" inequality are all computed and printed but never asserted. This is `insufficient_verification` and can be remedied without touching the paper or notes.

Attacks tried: (a) tried to break the existing assertions by checking that `gminus_exact - gminus = 0` is not algebraically a tautology — it isn't, because the construction `rF - sqrt(1+rF**2)/2` does not syntactically equal `(2*sqrt(4107-100*pi**2) - 37*sqrt(3))/(20*pi)`; SymPy must simplify the nested square root and recognize `sqrt(4107) = 37*sqrt(3)`. (b) Checked whether `rF` could be negative on some branch — `sqrt(...)` is non-negative in SymPy by default, and the notes state `𝔯_F1 ≈ 1.778 > 0`, consistent. (c) Checked the `R = 37/20` derivation — `R = sqrt(4107)/(20) = 37*sqrt(3)/20` from `sqrt(4107) = 37*sqrt(3)`; the script's `R = 37/20` is the rational normalization with `sqrt(3)` absorbed implicitly into `sqrt(1+rF²) = R/π·... ` (this is fine algebraically). (d) Checked the banner "STAGE 105" — it is stale text; not load-bearing.

## Self-test notes

(1) Variable independence: no derivatives in any proposed assertion, so the early-units derivative-zero trap does not apply. (2) Symmetry/parity: no integrals; not applicable. (3) Trivial-case pre-check for F2's proposed assertions: `compensation quadratic at gminus = 1 + rF² − 4(gminus − rF)² = 1 + rF² − 4·((rF − sqrt(1+rF²)/2) − rF)² = 1 + rF² − 4·(sqrt(1+rF²)/2)² = 1 + rF² − (1+rF²) = 0` ✓; symmetric for gplus. `defect_exact - comp_def`: both have identical numerator and denominator, simplify to 0 ✓. `gminus * (1/gminus) − 1 = 0` ✓. `comp_def != 0` because numerically ≈ 1.74 ≠ 0 ✓. (4) Path specs: F3 has no required action; F2 edits only `scripts/` per project policy. (5) Paper round-trip: F2's additions use the script's existing `4107 − 100π²` frame and do not introduce any new constants the paper/notes contradict (within the script-correct interpretation of F1).
