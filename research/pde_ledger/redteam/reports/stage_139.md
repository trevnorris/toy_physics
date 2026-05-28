---
unit_id: 139
batch: IV.5
auditor_model: claude-opus-4-7[1m]
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
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md
  paper_appendix: present
---

# Audit unit 139 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_139.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only contains `\input{stages/stage_139}` at line 1312; no inline prose)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.txt`

## What the paper claims

Stage 139 evaluates the explicit Family-1 mouth gains at the canonical compensation point `Pi_*`. The paper card's `Purpose` field says this is a "coupled mouth fixed point and gain selection ledger step" and the body quote says "Natural and exact-compensated gain pairs are evaluated on the Family--1 branch." The notes give the operational deliverables explicitly: for the natural equal-normalized branch (`g_c = 1`), `R_q^nat ≈ 0.145454452260421`, `M_s^nat,* ≈ 1.66854252965624`, `M_q^nat,* ≈ -0.242696939724365`; for the exact-compensated branch (`g_c = r_F1 - sqrt(1+r_F1^2)/2`), `R_q^comp = 1/4`, `M_s^comp,* ≈ 1.80594111095636`, `M_q^comp,* ≈ -0.451485277739090`. The card's `Checks` field enumerates three requirements: (a) check the gain pair `(M_s,M_q)` against outlet consistency `Pi_* = M_s + M_q * S_q(Pi_*)`; (b) check the self-matched susceptibility closure before using the one-scalar branch law; (c) record numerical fixed points as numerically located. Inputs are stated to include `Pi_*` and `S_q(Pi_*)` from upstream (Stage 236) and `r_F1` from Stage 223.

## What the script claims to verify

The SymPy script computes `r_F1` from the Stage 223 closed form, takes `Pi_* = 1.50882951349316` and `S_q(Pi_*) = 0.658075937605429` as hardcoded literals (matching the notes' Stage 236 values), then constructs `R_q^nat = (1-r_F1)^2/(1+r_F1^2)`, `M_s^nat = Pi_*/(1 - R_q^nat * S_q)`, `M_q^nat = -R_q^nat * M_s^nat` for the natural branch and the analogous trio (with `g_-^F1 = r_F1 - sqrt(1+r_F1^2)/2`, `R_q^comp = (g_-^F1 - r_F1)^2 / (1 + r_F1^2)`) for the compensated branch. It prints all values and the two comparison ratios. The *only* numerical assertion is `assert abs(R_q^comp - 1/4) < 1e-25` on line 30. The Mathematica script is a line-by-line transliteration with the same single check (`expectApprox["R_q^comp - 1/4", rQComp, 1/4, 10^-25]`). Nothing else is asserted; the rest of the output is informational printing.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side coverage |
|---|---|
| Natural-branch numerical values `M_s^nat,*`, `M_q^nat,*`, `R_q^nat` | partial — computed and printed but not asserted against any target |
| Compensated-branch numerical values `M_s^comp,*`, `M_q^comp,*` | partial — computed and printed but not asserted against any target |
| `R_q^comp = 1/4` | match — asserted in both engines |
| Check (a) outlet consistency `Pi_* = M_s + M_q * S_q(Pi_*)` | missing — never checked; by construction `M_q = -R_q M_s` and `M_s = Pi_*/(1 - R_q S_q)` makes this algebraically built in, so even if "checked" it would be tautological |
| Check (b) self-matched susceptibility closure | missing — no test exists |
| Check (c) numerical fixed points recorded | match — they appear in the printed transcript |
| `r_F1` derived from Stage 223 form | match (implicit) — computed and printed; not asserted against the notes' numeric 1.77799353547498 |

Dominant pattern: partial. Set `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 30 | `assert abs(R_q^comp - 1/4) < 1e-25` | R_q^comp = 1/4 | partial — assertion is algebraically true for any `rF` (see F1) |
| A2 | mathematica | 53 | `expectApprox["R_q^comp - 1/4", rQComp, 1/4, 10^-25]` | R_q^comp = 1/4 | partial — same algebraic-identity issue as A1 |

No assertions test the natural-branch values, the compensated-branch numerical M values, outlet consistency at Pi_*, or self-matched susceptibility closure.

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:53`

**What's wrong:**
The only check in either engine is `R_q^comp - 1/4 == 0`. With `g_minus = rF - sqrt(1+rF^2)/2`, we have
`g_minus - rF = -sqrt(1+rF^2)/2`, so `(g_minus - rF)^2 / (1 + rF^2) = (1+rF^2)/4 / (1+rF^2) = 1/4`,
which holds for *any* value of `rF`. The assertion therefore does not exercise anything Family-1-specific; substituting any other geometry would still pass. Furthermore, none of the paper's headline numerical deliverables (`M_s^nat,* ≈ 1.66854252965624`, `M_q^nat,* ≈ -0.242696939724365`, `R_q^nat ≈ 0.145454452260421`, `M_s^comp,* ≈ 1.80594111095636`, `M_q^comp,* ≈ -0.451485277739090`, `r_F1 ≈ 1.77799353547498`) are asserted; they are only printed. The card's three `Checks` items (outlet consistency, self-matched susceptibility closure, numerical fixed points) are not exercised by a single comparison.

**Why this matters:**
A regression that mistyped `Pi_*` or `S_q(Pi_*)` (both are bare-literal imports without provenance assertion), or that swapped the natural/compensated `R_q` formulas, would not be caught. The script's "PASS" line is supported by a single algebraic identity that is independent of all imported physics.

**Required change:**
Add `expectApprox`-style closeness checks (and SymPy analogs) for each numerical deliverable named in the notes:
1. `r_F1` against the notes' numeric 1.77799353547498 (or, better, against the exact `sqrt(12/pi^2 * (37/20)^2 - 1)` closed form so a deviation in either side fails).
2. `R_q^nat` against `(1 - r_F1)^2 / (1 + r_F1^2)` AND against the notes' literal `0.145454452260421` (the second check pins the imported r_F1 value).
3. `M_s^nat,*` against `1.66854252965624` (literal from notes).
4. `M_q^nat,*` against `-0.242696939724365` (literal from notes).
5. `M_s^comp,*` against `1.80594111095636` (literal from notes).
6. `M_q^comp,*` against `-0.451485277739090` (literal from notes).
7. Outlet consistency: assert `|Pi_star - (M_s^nat + M_q^nat * S_q)| < 10^-25` AND the same for the compensated pair. NOTE: by construction `M_s = Pi_*/(1 - R_q S_q)` and `M_q = -R_q M_s` makes this algebraically trivial, but it still verifies that the *imported* Pi_* and S_q literals satisfy the form (this is the operational meaning of paper check (a)).
8. Cross-check `R_q^nat = (1 - r_F1)^2 / (1 + r_F1^2)` against the numeric `0.145454452260421`; if the bare-literal `Pi_*` were perturbed by a regression, an absolute-value check here is the only thing that flags it.

Choose tolerances `10^-12` on the literals (notes only give ~15 digits) and `10^-25` on algebraic identities.

**Verification:**
The verifier confirms by counting assertions: the sympy script must contain at least 6 `assert` statements (one per numerical deliverable above plus the existing R_q^comp check) and the Mathematica must contain matching `expectApprox` calls. Output transcript must show at least 6 PASS lines.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:28-39`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:5-16`

**What's wrong:**
The `.wl` is a syntactic transliteration of the `.py`. Comparing the two construction blocks:

SymPy (lines 5-16):
```
rF = sp.N(sp.sqrt(sp.Rational(12,1)/sp.pi**2 * (sp.Rational(37,20))**2 - 1), 30)
Pi_star = sp.N('1.50882951349316', 30)
Sq_star = sp.N('0.658075937605429', 30)
Rq_nat = sp.N((1 - rF)**2 / (1 + rF**2), 30)
Ms_nat = sp.N(Pi_star / (1 - Rq_nat * Sq_star), 30)
Mq_nat = sp.N(-Rq_nat * Ms_nat, 30)
g_minus = sp.N(rF - sp.sqrt(1 + rF**2)/2, 30)
Rq_comp = sp.N((g_minus - rF)**2 / (1 + rF**2), 30)
Ms_comp = sp.N(Pi_star / (1 - Rq_comp * Sq_star), 30)
Mq_comp = sp.N(-Rq_comp * Ms_comp, 30)
```

Mathematica (lines 28-39):
```
rF = N[Sqrt[(12/Pi^2)*(37/20)^2 - 1], 30];
piStar = SetPrecision[1.50882951349316, 30];
sQStar = SetPrecision[0.658075937605429, 30];
rQNat = N[(1 - rF)^2/(1 + rF^2), 30];
mSNat = N[piStar/(1 - rQNat*sQStar), 30];
mQNat = N[-rQNat*mSNat, 30];
gMinus = N[rF - Sqrt[1 + rF^2]/2, 30];
rQComp = N[(gMinus - rF)^2/(1 + rF^2), 30];
mSComp = N[piStar/(1 - rQComp*sQStar), 30];
mQComp = N[-rQComp*mSComp, 30];
```

Line-by-line correspondence is exact: same variable order, same intermediate names (rF↔rF, Pi_star↔piStar, etc.), same single assertion target. The Mathematica script does not re-derive `g_-` or `R_q` from first principles (e.g., it does not solve the compensated-branch condition to recover `g_c = r - sqrt(1+r^2)/2`); it just rewrites the SymPy algebra in Mathematica syntax. This violates the second-engine policy: the second engine must derive, not echo.

**Why this matters:**
If the SymPy formula were wrong (e.g., wrong sign in `g_minus`), the Mathematica would carry the identical error and the "engines agree" property would be vacuous.

**Required change:**
In the Mathematica script, replace the bare assignment of `gMinus` with a derivation: solve `(gc - rF)^2 / (1 + rF^2) == 1/4` symbolically via `Solve[(gc - rF)^2 == (1 + rF^2)/4, gc]`, select the branch with `gc < rF` (matching the notes' lower-compensated branch), and assign that solution to `gMinus`. Then assert that the resulting `gMinus` equals `rF - Sqrt[1+rF^2]/2` numerically. This makes the Mathematica side an independent re-derivation rather than an echo.

**Verification:**
After the fix, `gMinus` definition in the .wl file must contain a `Solve[...]` or `Reduce[...]` call rather than a direct `rF - Sqrt[...]/2` assignment, and a new `expectApprox` line must check the agreement against the closed form.

### F3 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:6-7`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:29-30`

**What's wrong:**
`Pi_star = 1.50882951349316` and `Sq_star = 0.658075937605429` are introduced as bare numeric literals with no comment citing their upstream provenance. The paper card's `Inputs` field says these come from the mouth-source / outlet-consistency block (and the notes attribute them to Stage 236), but the script provides no in-file anchor (no comment string `# from stage 236` and no cross-check assertion). If an upstream value drifted, this script would silently use the stale literal forever.

**Why this matters:**
This stage's entire downstream numerical content rests on these two literals. The paper card's check (b) ("self-matched susceptibility closure before using the one-scalar branch law") implicitly requires that `S_q(Pi_*)` be verified against the susceptibility law at the chosen `Pi_*`. Neither happens in this script.

**Required change:**
Add comments at lines 6-7 (sympy) and 29-30 (mathematica) explicitly citing the upstream source: `# from Stage 236 (Pi_*, S_q(Pi_*))` and `# from Stage 223 closed form for r_F1`. This is a minimum; ideally the literals would be loaded from a shared upstream-results module, but that is a refactor and out of scope for this directive.

**Verification:**
The .py and .wl files contain provenance comments adjacent to the `Pi_star` and `Sq_star` definitions naming the source stage.

### F4 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`

**What's wrong:**
Script mtime: 2025-04-01 21:36. Output mtime: 2025-05-11 12:45 (so output is newer than script — fine). For the Mathematica side, script mtime: 2026-05-11 11:56, output mtime: 2026-05-11 13:13 (fine). Also, the Mathematica banner at line 26 reads `"STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS"`, but this file is the Stage 139 audit. The misleading banner is reflected in the saved output transcript line 11 (`STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS`). Outputs are not stale by mtime, but the transcripts contain a stage-number label that is wrong. This is informational rather than mtime-stale, but flagging here because the auditor expects the output to reflect the current script faithfully and the banner is misleading enough that a reviewer might think they are looking at Stage 122's output.

**Why this matters:**
Output transcripts are an artifact the verifier reads; a wrong stage label increases the chance of a downstream cross-reference error.

**Required change:**
Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:26` from
`banner["STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS"];`
to
`banner["STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS"];`
The verifier will re-run mathematica and the transcript banner will then read `STAGE 139`.

**Verification:**
After the fix, the Mathematica output transcript's banner line reads `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script (see F2 for quoted side-by-side excerpts). The same `g_minus = rF - sqrt(1+rF^2)/2` assignment appears in both with no independent derivation in either; the Mathematica script does not solve the compensated-branch defining equation. This is the basis for the F2 `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce numerically identical results to ≥28 digits across all printed quantities. The single PASS line (`R_q^comp - 1/4 ~ 0`) appears in both transcripts with residual `0` (Mathematica) and assertion passing (SymPy). Engines agree; but as F2 explains, the agreement is vacuous because they share the same construction.

## Verdict justification

The scripts run, the engines numerically agree, the printed numbers reproduce the notes' boxed values to ~15 digits, and the one assertion (`R_q^comp = 1/4`) passes. However, that assertion is an algebraic identity holding for any `rF` (not Family-1 specific), and none of the paper's headline numerical deliverables (six numbers boxed in the notes) or the three `Checks` items on the card are asserted. The Mathematica script is a transliteration of the SymPy script rather than an independent derivation, and the imported literals `Pi_*` and `S_q(Pi_*)` carry no in-file provenance. The Mathematica banner labels the unit as Stage 122. Verdict: `findings`, four findings, no stop-cold (no `paper_misalignment` because the paper and script claims are mutually compatible — the script just doesn't actually verify most of what the paper says it does).

## Self-test notes

I checked: (1) the `R_q^comp = 1/4` identity by direct algebra and confirmed it holds for any `rF` (not Family-1 specific), justifying the F1 insufficiency claim. (2) The outlet-consistency identity `Pi_* = M_s + M_q * S_q` is built in by construction `M_s = Pi_*/(1 - R_q S_q)` and `M_q = -R_q M_s`, so a literal coding of "check this" would be tautological — F1's required change adds it anyway because it serves to pin the imported `Pi_*` and `S_q` literals. (3) I did not invent new derivations; I only ask Codex to add `expectApprox`/`assert` calls against literal targets stated in the notes, and one `Solve` call for the Mathematica `g_minus` definition. (4) Mathematica banner typo is a string edit only; no math impact.
