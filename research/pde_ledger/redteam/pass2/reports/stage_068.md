---
unit_id: 068
batch: III.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage068_resonance_thresholds.md]
  paper_appendix: present
---

# Audit unit 068 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_068.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage068_resonance_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 114)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 068 takes the Stage-066 matched-branch wall threshold window and shows that inserting the Stage-067 sech–Gaussian benchmark coherence penalty changes the matched thresholds only mildly. The card's `\stagefield{Output}` is: "Threshold correction factor \(P_{\rm res}\)", with the body equation \(P_{\rm res}\simeq1.005612487760576\) (\eqref{eq:app-stage068-Pres}), i.e. about a 0.56% raise of the required wall figure of merit — "not by an order-one factor". The notes are authoritative on intent and enumerate four deliverables: (1) the resonance-corrected wall figure `W_res(r) = C^2(r) W_wall`; (2) exact threshold translation `W_fail^(res) = Pe_req/[C^2 Delta_inf]`, `W_suff^(res) = Pe_req/[C^2 Delta_0]`; (3) the resonance-point penalty `P_res = 1/C_res^2 = 1.005612487760576...` derived from `C_res^2 = 0.994418836451529...` carried from Stage 067; (4) the profile-sensitive band widths, each ≈ 0.56% on both the success side (`Pe_req/Delta_0 <= W_wall < P_res Pe_req/Delta_0`) and failure side. The appendix row (line 114) summarizes: "Matched thresholds modified by \(P_{\rm res}\simeq1.005612\)."

## What the script claims to verify

Both scripts (a) build `W_res` and `W_wall` independently from a microscopic gain decomposition `G_match = rho_star g_phi^2 N_phiphi/(m c_s^2 K_X)` with `kappa`, and assert `W_res - C2*W_wall == 0` plus the matched limit `C2->1`; (b) derive `P_res = 1/C_res^2` from the ratio of required wall figures and assert it equals `1/Cres^2`, then numerically cross-check that `1/0.994418836451529` matches the paper's `1.005612487760576` to < 1e-12; (c) solve the matched and profile-family Peclet balances for the thresholds and assert the profile/matched ratio equals `1/C2` (and scales by `Pres` at `C2->1/Pres`); (d) compute the band widths two ways (C-form via `Cres^2`, P-form via `(Pres-1)`) and assert they agree under `Pres = 1/Cres^2`. This is what the verdict applies to.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `W_res(r) = C^2(r) W_wall` | py:62-63 / wl:43-44 `W_res - C2*W_wall == 0` (built from independent gain components), plus matched limit py:65-66 / wl:45-46 | match |
| `W_fail^(res)=Pe_req/[C^2 Delta_inf]`, `W_suff^(res)=Pe_req/[C^2 Delta_0]` | py:117-127 / wl:70-91 (Solve/Reduce balances), printed at py:129-132 / wl:93-96; ratio asserts py:136-137 / wl:99-100 | match |
| `P_res = 1/C_res^2 = 1.005612487760576` | symbolic py:79-89 / wl:49-55; numeric anchor py:91-101 / wl:57-66 | match |
| Resonance-point scaling by `P_res` | py:140-143 / wl:103-104 | match |
| Profile-sensitive band widths ≈0.56% (success+failure) | py:157-174 / wl:107-125 two-form cross-check | match |

`paper_alignment: aligned` — every paper deliverable has a faithful, non-tautological script-side check, and the constants match.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 62-63 | `simplify(Wres - C2*Wwall) == 0` | deliverable 1 | yes |
| A2 | sympy | 65-66 | matched limit `C2->1` == 0 | deliverable 1 (sanity) | yes |
| A3 | sympy | 88-89 | `Pres_from_ratio - 1/Cres^2 == 0` | deliverable 3 (symbolic) | yes |
| A4 | sympy | 98-101 | `|1/Cres^2 - 1.005612...| < 1e-12` | deliverable 3 (numeric) | yes |
| A5 | sympy | 136-137 | `Wfail_res*C2 - Wfail_match == 0` | deliverable 2 | yes |
| A6 | sympy | 140-143 | scale by `Pres` at `C2->1/Pres` | deliverable 4 (resonance pt) | yes |
| A7 | sympy | 171-174 | band C-form vs P-form under `Pres=1/Cres^2` | deliverable 4 | yes |
| A8 | math | 43-46 | `Wres - C2*Wwall == 0` + matched limit | deliverable 1 | yes |
| A9 | math | 54-55 | `PresFromRatio - 1/Cres^2 == 0` | deliverable 3 (symbolic) | partial (see Independent-derivation note) |
| A10 | math | 61-65 | numeric anchor < 1e-12 | deliverable 3 (numeric) | yes |
| A11 | math | 99-104 | threshold ratio + resonance scaling | deliverables 2,4 | yes |
| A12 | math | 122-125 | band C-form vs P-form | deliverable 4 | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.txt` (mtime 2026-05-22 20:00:11)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.txt` (mtime 2026-05-22 20:00:17)

**What's wrong:**
Both saved outputs predate the current scripts (`.py`/`.wl` mtime 2026-06-03 15:59:11), and their CONTENT disagrees with what the current scripts would emit:
- SymPy output line 3 banner reads `STAGE 51 — ...`; current script line 37 prints `STAGE 68 — ...`. Mathematica output lines 3 banner reads `STAGE 051 — ...`; current script line 26 prints `STAGE 068 — ...`.
- SymPy output line 7 contains `P_res*C_res^2 - 1 = 0`; the current script does not emit any such line (it prints `P_res numeric residual = ...` at py:97, which is absent from the saved output).
- Band-width labels in both outputs read `(A)`/`(B)` (e.g. SymPy line 16 `Success-side band width (A) = ...`); current scripts print `(C-form)`/`(P-form)` (py:165-168, wl:117-120).
These confirm the captured transcripts are stale.

**Why this matters:**
The committed transcript is the auditable record. It currently shows a different stage number and a different set of printed lines than the present scripts produce, so a reader cannot trust it as the script's output.

**Required change:**
Re-run both scripts and overwrite the two `.txt` files with the fresh transcripts (orchestrator exec-refresh / sed-refresh per the pass-2 infra note). No script-math change required for this finding.

**Verification:**
After refresh, SymPy `.txt` line 3 reads `STAGE 68 — RESONANCE-CORRECTED THRESHOLDS`, contains a `P_res numeric residual = ...` line, and band labels read `(C-form)`/`(P-form)`; Mathematica `.txt` line 3 reads `STAGE 068 — ...` with `(C-form)`/`(P-form)` labels and all `PASS:` lines present.

### F2 — stale_output (self-label in SymPy docstring)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:4`

**What's wrong:**
The module docstring's self-label line reads `moving_throat_pde_stage51_resonance_thresholds_sympy_audit.py` — a stale "stage51" self-label inside a Stage-068 script (filename and Mathematica banner are correctly `068`; the SymPy banner at py:37 is also a stale `STAGE 68`, but the docstring filename label is the unambiguous self-label). Per the in-loop Reading-2 policy, a verdict:findings stage gets its unambiguous self-labels fixed.

**Why this matters:**
The mislabeled docstring filename misrepresents which stage this script audits and is the root cause of the stale `STAGE 51`/`STAGE 068` banner mismatch propagating into the saved outputs.

**Required change:**
At `py:4`, change the docstring filename from `moving_throat_pde_stage51_resonance_thresholds_sympy_audit.py` to `moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`. At `py:37`, change the banner string `"STAGE 68 — RESONANCE-CORRECTED THRESHOLDS"` to `"STAGE 068 — RESONANCE-CORRECTED THRESHOLDS"` to match the Mathematica banner (wl:26) and the canonical 3-digit stage label. Re-run so the refreshed SymPy `.txt` banner reads `STAGE 068`.

**Verification:**
`py:4` filename ends in `stage068`; `py:37` banner reads `STAGE 068`; refreshed SymPy `.txt` line 3 reads `STAGE 068 — RESONANCE-CORRECTED THRESHOLDS`.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a line-by-line transliteration — it uses genuinely different machinery in the threshold section: SymPy uses `sp.solve(Eq(W*Delta, Pe_req), W)` (py:117,125) whereas Mathematica uses `Reduce[W*Delta == PeReq && W>0, W, Reals]` with `Cases`/`HoldPattern` extraction (wl:70-91), an independent route deliberately chosen "to keep the positivity premise explicit". The gain-decomposition build (wl:36-46) and band-width two-form check (wl:107-125) mirror the SymPy structure conceptually but are re-expressed in native Mathematica. One soft spot: the `PresFromRatio` step (wl:49-53) hand-writes the ratio `(PeReq/(C2*Delta1))/(PeReq/Delta1)` instead of solving the two balances as SymPy does (py:79-86); it still collapses to `1/C2` and the assertion is non-tautological against `1/Cres^2`, so this is a weaker-but-valid derivation, not a transliteration and not a finding. Overall verdict: independent.

## Engine cross-check

Both engines reach identical final forms (modulo Mathematica's `(-1+Pres)` vs SymPy's `(Pres-1)` ordering and the stale stage labels): `W_res - C2*W_wall = 0`, `P_res - 1/C_res^2 = 0`, matched thresholds `Pe_req/Delta_inf`, `Pe_req/Delta_0`, profile thresholds `Pe_req/(C2*Delta_inf)`, `Pe_req/(C2*Delta_0)`, and band widths `Pe_req*(P_res-1)/Delta_{0,inf}`. All cross-relation and band-equality checks pass in both transcripts. Engines agree.

## Verdict justification

The scripts hold up against every paper-side deliverable: the four claims (W_res factorization, exact threshold translation, P_res = 1/C_res^2 with the exact numeric anchor 1.005612487760576 ↔ 0.994418836451529, and the two-route band widths) are each exercised by a non-tautological, well-anchored assertion in both engines, and the constants match the card/notes/appendix exactly. Attacks tried that failed: (i) checked the numeric anchor for self-reference — it independently links the two Stage-067-carried numbers via the reciprocal, not a number-against-itself; (ii) checked the threshold ratio for tautology — the matched and profile balances are solved independently and their ratio is a genuine `1/C2`, not constructed; (iii) checked symbol domains — all symbols positive/real, consistent with physical figures of merit, no hidden assumption that forces a pass; (iv) checked the band two-form cross-check — perturbing `Pres ≠ 1/Cres^2` would break the equality, so it is sensitive. The only defects are non-math: both committed `.txt` transcripts are stale (older mtime AND divergent content — old stage banners, an obsolete printed line, old `(A)/(B)` band labels) and the SymPy docstring carries a stale `stage51` self-label (with a matching `STAGE 68` non-3-digit banner). Hence `verdict: findings`, all low-severity, `material_change` false; no `paper_misalignment`, so no user gate is required by the standard audit.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `P_res = 1.005612487760576` | py:95 `Pres_paper`, wl:60; numeric anchor py:97 / wl:62 | `.tex:16` (`P_{\rm res}\simeq1.005612487760576`), `.tex:18`/`.tex:20` Output, appendix `:114` (`\simeq1.005612`); `.md:16`,`.md:88`,`.md:119`,`.md:161` | MATCH |
| `C_res^2 = 0.994418836451529` | py:93 `Cres_sq_numeric`, wl:59 | `.md:12` (`C_res^2 = 0.994418836451529...`); card states "carried from stage 067" | MATCH (lives in notes) |
| `W_res(r) = C^2(r) W_wall` (symbolic) | py:178 / wl:128, asserted py:62-63 / wl:43-44 | `.md:51-55` (`W_res(r) = C^2(r) W_wall`) | MATCH |
| `W_fail^(res) = Pe_req/[C^2 Delta_inf]` | py:179 / output; asserted py:126,136 | `.md:77` (`W_wall <= Pe_req/[C^2(r) Delta_inf]`) | MATCH |
| `W_suff^(res) = Pe_req/[C^2 Delta_0]` | py:180 / output; asserted py:127,137 | `.md:78` | MATCH |
| Resonance-point: `W_{fail,suff}^(res,*) = P_res Pe_req/Delta_{inf,0}` | py:182-183 / output; asserted py:140-143 | `.md:82-84` | MATCH |
| Band widths `(P_res-1) Pe_req/Delta_{0,inf}` (≈0.56%) | py:162-168 / output 16-19; wl:114-120 | `.md:128-141` (band width ≈ `0.56%`) | MATCH |

INTERNAL scaffolding (accounted for, no finding): gain-decomposition component symbols (`rho_star,g_phi,N_phiphi,m,c_s,K_X,kappa` and `Gmatch_expr/Wwall_expr/Gres_expr/Wres_expr`), the solve-helper dummy symbols (`Wm,Wp,W_match,W_prof,Delta,Delta1,Ain,Atrans`), residual/PASS flags, the `1e-12` tolerance, and the matched-limit `C2->1` sanity check.

reconciliation: complete; 7 deliverable values checked, 0 misaligned.

## Self-test notes

Checked variable independence (no `diff`/`D` in this stage — all checks are algebraic identities, so the zero-derivative trap is N/A). Checked symbol domains: all positive/real, matching the physical figure-of-merit setup; no over-strong assumption forces a pass — perturbing `Pres ≠ 1/Cres^2` breaks the band cross-check (py:171-174 / wl:122-125), confirming sensitivity. Trivial-case pre-check on the numeric anchor: `1/0.994418836451529 = 1.005612487760...`, matching the paper literal to < 1e-12, so the assertion is real and would fail if either Stage-067-carried number were altered. The two findings are output-freshness/self-label only and require no math edit, so no new `paper_misalignment` is introduced.
