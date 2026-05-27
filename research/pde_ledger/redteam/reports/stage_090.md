---
unit_id: 090
batch: III.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage090_updated_reduced_status.md]
  paper_appendix: present
---

# Audit unit 090 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_090.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage090_updated_reduced_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.txt`

## What the paper claims

Stage 090 is a checkpoint-status update for Part III's reduced moving-throat program. The card's `\stagefield{Output}` states verbatim: "The support/source side is finished under the minimal module. The next theorem gate is derivation of the module itself." The body enumerates six settled facts inside the reduced Family-1 branch: (i) explicit Family-1 support/source has a hard constructive ceiling; (ii) outgoing normalization factors cancel from the support theorem; (iii) the theorem depends only on `rho_alpha = alpha_req/alpha_mix`; (iv) the natural contact-plus-pole reading gives `rho_alpha = 4/3`; (v) this lies strictly inside the Family-1 success region; (vi) the explicit branch succeeds at zero transport bias. The notes additionally lock the triple `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pe_req = 0`. The Part III appendix row (line 158) summarizes the status as "Support/source side finished under the minimal module; remaining gap is deriving that module."

## What the script claims to verify

The SymPy script (docstring lines 2-20) declares the audit a checkpoint-consistency check that reconfirms `rho_alpha = 4/3` and `zeta_req = 1/3` and places the carried Family-1 thresholds (`rho_suff^(chi)`, `zeta_max^(F1)`, `A_F1`) on the correct side of the branch values. It does this by computing `rho_alpha = 1/c_contact` and `zeta_req = c_pole/c_contact` from `c_contact = 3/4`, `c_pole = 1/4`, then asserting equality to `4/3` and `1/3`, the identity `zeta_req = rho_alpha - 1`, and three inequality margins. The Mathematica script asserts only the identity `zeta_req = rho_alpha - 1` (with both quantities hardcoded) and the three inequalities; it does not derive `rho_alpha` from the contact-plus-pole coefficients and does not assert `rho_alpha = 4/3` or `zeta_req = 1/3` numerically.

## Paper ↔ script cross-check

| Paper deliverable | SymPy coverage | Mathematica coverage |
|---|---|---|
| Body (iii): theorem depends on `rho_alpha` | implicit (uses `rho_alpha`) | implicit |
| Body (iv) / notes: `rho_alpha = 4/3` from contact-plus-pole | match (derived from `c_contact = 3/4`, asserted = 4/3 at line 63) | missing — value is hardcoded at .wl:34, never derived from `(c_contact, c_pole)` |
| Notes: `zeta_req = 1/3` | match (derived `c_pole/c_contact` and asserted = 1/3 at line 64) | missing — `zetaReq = rhoAlpha - 1` is a definition, not a derivation; never asserted = 1/3 |
| Notes identity: `zeta_req = rho_alpha - 1` | match (line 65) | tautology — `expectZero[zetaReq - (rhoAlpha - 1)]` where `zetaReq := rhoAlpha - 1` (lines 35, 46) |
| Body (v): `rho_alpha < rho_suff^(chi)` | match (line 88-91) | match (.wl:47) |
| Body (i)/(vi): support ceiling and zero-bias success | partial — represented via `zeta_req < zeta_max_f1` and `zeta_req < A_F1` (lines 92-99) | partial (.wl:48-49) |
| Notes: `Pe_req = 0` | missing — no `Pe_req` symbol; the `zeta_req < A_F1` check is treated as a proxy in the print only | missing |

Dominant pattern: partial paper alignment. SymPy faithfully exercises every deliverable except the explicit `Pe_req = 0` symbol. Mathematica is materially weaker: it skips the contact-plus-pole derivation and replaces the substantive `rho_alpha = 4/3` and `zeta_req = 1/3` checks with a tautology.

## Assertion inventory

| #  | Script        | Line     | Form                                                | Exercises which paper claim?                          | Anchored to claim? |
|----|---------------|----------|-----------------------------------------------------|-------------------------------------------------------|--------------------|
| A1 | sympy         | 63       | `expect_zero("rho_alpha - 4/3", rho_alpha - 4/3)`   | body (iv); notes triple                               | yes                |
| A2 | sympy         | 64       | `expect_zero("zeta_req - 1/3", zeta_req - 1/3)`     | notes triple                                          | yes                |
| A3 | sympy         | 65       | `expect_zero("zeta_req - (rho_alpha - 1)", ...)`    | identity used in body (iii)/(iv)                      | yes                |
| A4 | sympy         | 88-91    | `expect_true(rho_alpha < rho_suff_chi)`             | body (v)                                              | yes                |
| A5 | sympy         | 92-95    | `expect_true(zeta_req < zeta_max_f1)`               | body (i) ceiling                                      | yes                |
| A6 | sympy         | 96-99    | `expect_true(zeta_req < A_F1)`                      | body (vi) zero-bias (proxy for `Pe_req = 0`)          | partial            |
| M1 | mathematica   | 46       | `expectZero["zeta_req - (rho_alpha - 1)", ...]`     | identity                                              | no — tautological  |
| M2 | mathematica   | 47       | `expectTrue[rhoAlpha < rhoSuffChi]`                 | body (v)                                              | yes                |
| M3 | mathematica   | 48       | `expectTrue[zetaReq < zetaMaxF1]`                   | body (i) ceiling                                      | yes                |
| M4 | mathematica   | 49       | `expectTrue[zetaReq < aF1]`                         | body (vi) (proxy for `Pe_req = 0`)                    | partial            |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl:34-46`

**What's wrong:**
The Mathematica script defines

```
rhoAlpha = 4/3;
zetaReq = rhoAlpha - 1;
```

at lines 34-35 and then asserts

```
expectZero["zeta_req - (rho_alpha - 1)", zetaReq - (rhoAlpha - 1)];
```

at line 46. By construction `zetaReq - (rhoAlpha - 1)` is identically `0` for any value assigned to `rhoAlpha`. The assertion cannot fail regardless of whether `rho_alpha = 4/3` or `zeta_req = 1/3` actually follow from the contact-plus-pole coefficients `(c_contact, c_pole) = (3/4, 1/4)`. The SymPy script (lines 52-65) does the substantive thing: it computes `rho_alpha = 1/c_contact` and `zeta_req = c_pole/c_contact` from the upstream module and then asserts both equal `4/3` and `1/3` respectively.

**Why this matters:**
This is a checkpoint stage. The paper's body item (iv) and the notes' locked triple `(rho_alpha = 4/3, zeta_req = 1/3, Pe_req = 0)` are exactly the contact-plus-pole arithmetic that this stage is supposed to certify. The Mathematica engine does not certify any of it — it accepts `rho_alpha = 4/3` as input. The two-engine policy demands independent re-derivation of the load-bearing identity, not an echo of the SymPy result.

**Required change:**
Introduce `cContact = 3/4` and `cPole = 1/4` as inputs at `.wl:34`. Derive `rhoAlpha = 1/cContact` and `zetaReq = cPole/cContact`. Replace the tautology at line 46 with two substantive `expectZero` calls anchored to `4/3` and `1/3`, then keep the identity check as a third assertion that now exercises something (since `zetaReq` and `rhoAlpha` are independently derived).

**Verification:**
Re-run Mathematica audit. New lines should compute `rhoAlpha` and `zetaReq` from `(cContact, cPole)` and assert `rhoAlpha - 4/3 == 0`, `zetaReq - 1/3 == 0`, and `zetaReq - (rhoAlpha - 1) == 0` non-tautologically. Script must exit 0.

### F2 — script_missing_paper_claim

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl:34-49`

**What's wrong:**
The notes (lines 17-21) state the locked triple `rho_alpha = 4/3, zeta_req = 1/3, Pe_req = 0`, and the paper body item (iv) requires `rho_alpha = 4/3` from the contact-plus-pole reading. The Mathematica script never asserts `rhoAlpha - 4/3 == 0` or `zetaReq - 1/3 == 0` — both values are merely defined as constants. The SymPy script does assert both (lines 63-64). Mathematica is therefore missing the paper-side deliverables that SymPy covers; this fails the "both engines required" bar for a checkpoint stage.

**Why this matters:**
For checkpoint stages, paper alignment must be exact and both engines must substantively verify the same paper-side claims. As written, the Mathematica script's load-bearing content reduces to inequality checks on numerically hardcoded thresholds. A Mathematica-only run cannot certify the `(rho_alpha, zeta_req)` triple at all.

**Required change:**
Same edit as F1 — once `rhoAlpha` and `zetaReq` are derived from `(cContact, cPole)`, add the explicit `expectZero["rho_alpha - 4/3", rhoAlpha - 4/3]` and `expectZero["zeta_req - 1/3", zetaReq - 1/3]` assertions at `.wl:46` so the Mathematica script mirrors SymPy lines 63-64.

**Verification:**
Mathematica output must show `PASS: rho_alpha - 4/3` and `PASS: zeta_req - 1/3` lines. Engines now cross-check on the same paper-side deliverables.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py:96-106`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl:49`

**What's wrong:**
The notes (line 19) and the script's printed conclusion (sympy lines 104-106) state the deliverable `Pe_req = 0`. The scripts represent this only as the inequality `zeta_req < A_F1` (sympy line 96-99, mathematica line 49). No `Pe_req` symbol is defined; no equation links the inequality to `Pe_req = 0`. The implication "zero-bias baseline lying above `zeta_req` ⇒ `Pe_req = 0`" is left as a print-only assertion.

**Why this matters:**
The paper card body item (vi) and the notes' locked triple include `Pe_req = 0` as an explicit deliverable. The proxy inequality is necessary but not sufficient to certify the value: the reader cannot trace from the assertion to the symbolic `Pe_req = 0` without an external lemma. For a checkpoint stage this is below the "substantive assertion" bar, though it is the lightest of the three findings because the proxy is a defensible carry-forward if Stage 062's transport map already certifies the equivalence.

**Required change:**
Either (a) add a literal `Pe_req = sp.Integer(0)` symbol with the carried derivation comment (e.g. "from Stage 062 transport map: `zeta_req < A_F1` ⇒ `Pe_req = 0`") and `expect_zero("Pe_req", Pe_req)`, mirroring in Mathematica; or (b) keep the proxy but add an inline comment at sympy:96 and mathematica:49 stating "Stage 062 transport map: `zeta_req < A_F1` ⇒ `Pe_req = 0`; this inequality is the carry-forward proxy for that conclusion." Option (b) is the lighter-touch fix and adequate for a status-checkpoint stage.

**Verification:**
If option (a): new output lines `Pe_req = 0` and `PASS: Pe_req` appear in both engines. If option (b): the comment appears at the cited lines and the output is unchanged.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** a transliteration in the strict sense — it has fewer assertions than the SymPy script, not the same ones in Mathematica syntax. However, it is also **not** an independent derivation: the load-bearing identity `rho_alpha = 4/3` is hardcoded at `.wl:34` rather than derived from `(cContact, cPole)`. SymPy does the derivation at lines 52-56 (`rho_alpha = sp.simplify(1 / c_contact)`); Mathematica skips that step entirely. The correct fix is to make Mathematica do the same derivation independently — this would convert M1 from a tautology into a substantive check and add the missing M2/M3 assertions for `rho_alpha - 4/3` and `zeta_req - 1/3`. Marking this as a `tautological_check` (F1) and `script_missing_paper_claim` (F2) rather than `mathematica_transliteration` since the .wl is shorter than the .py rather than mirroring it.

## Engine cross-check

Both engines exit 0 and agree on every inequality:

- `rho_alpha < rho_suff^(chi)`: both `True` (3.466... > 1.333...)
- `zeta_req < zeta_max^(F1)`: both `True` (2.467... > 0.333...)
- `zeta_req < A_F1`: both `True` (1.000... > 0.333...)
- `zeta_req = rho_alpha - 1`: both `0`

The disagreement is in coverage, not in value: SymPy additionally asserts `rho_alpha = 4/3` and `zeta_req = 1/3` (residuals `0`); Mathematica does not. There is no numerical disagreement on the assertions that both engines run. The numerical thresholds match to full quoted precision (Stage 69 / 63-64 / 62 carry-forward values are identical between scripts). Outputs are fresh: scripts mtime is `May 11 11:56`; sympy output mtime is `May 11 12:44`; mathematica output mtime is `May 11 13:02`. No stale_output finding.

## Verdict justification

The SymPy script holds up: every paper-side deliverable except the explicit `Pe_req = 0` symbol has a non-tautological, well-anchored assertion that exercises the contact-plus-pole arithmetic, the carry-forward thresholds, and the inequality margins. The Mathematica script does not hold up to the checkpoint-stage bar: the identity check is a definitional tautology, and the load-bearing `rho_alpha = 4/3` and `zeta_req = 1/3` are nowhere asserted. The `Pe_req = 0` deliverable is represented in both engines only as a proxy inequality with no symbol. None of the findings are mathematically wrong — the SymPy script and the paper agree on every value and every inequality — so neither `UNFIXABLE` nor `CRITICAL_DOWNSTREAM` is warranted. Two cosmetic notes (not findings): both scripts label their banner "STAGE 073" rather than "STAGE 090", and the notes file references `scripts/moving_throat_pde_stage141_*` rather than `stage090_*`; these are mislabels, not substantive errors. Verdict: `findings` (three), `paper_alignment: partial`, fixable by tightening the Mathematica engine and optionally adding a `Pe_req` symbol in both engines.

## Self-test notes

Walked through the proposed F1/F2 fix: replacing `rhoAlpha = 4/3` with `cContact = 3/4; cPole = 1/4; rhoAlpha = 1/cContact; zetaReq = cPole/cContact` makes `zetaReq - (rhoAlpha - 1) = 1/3 - (4/3 - 1) = 0` (no longer tautological — it now depends on `cContact = cPole + 1/2`, which is what the upstream module fixes), and `rhoAlpha - 4/3 = 1/(3/4) - 4/3 = 0` and `zetaReq - 1/3 = (1/4)/(3/4) - 1/3 = 0`. All three reduce to `0` non-trivially. No derivative or symmetry traps in this stage (no differentiation, no integration). Path specifications are correct: the .wl already exists at `mathematica/`, so no missing-script reconstruction is needed. Paper round-trip: the fix introduces no new constants beyond `(3/4, 1/4)` which are already on the paper side (notes line 17-19, body item iv). No new `paper_misalignment` introduced.
