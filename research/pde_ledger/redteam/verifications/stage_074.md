---
unit_id: 074
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification -- unit 074 (v2)

This is the v2 verifier pass. The v1 audit (resolved in the prior pass)
hardened the SymPy substitution chain (F1 `tautological_check`) and added
a provenance comment (F2 `insufficient_verification`); both were verified
clean. The v2 audit raised one new finding -- a paper-side internal
inconsistency in the boxed `alpha` value -- which is the subject of this
verification.

## Per-finding outcomes

### F1 -- paper_misalignment (alpha = 128/sqrt(5) vs 179 vs 111)

**Classification:** resolved

**What changed:**

Per the directive front-matter and the captured diff
(`redteam/exec_logs/stage_074_diff.patch`), the orchestrator applied the
user's chosen direction (a) directly (no Codex iteration needed):

- Paper-side: `paper/stages/stage_074.tex:31`, `\frac{128}{\sqrt5}` ->
  `\frac{111}{\sqrt5}` (recorded in directive; verifier does not re-read
  prose per scope rules).
- Notes-side: `notes/stages/.../stage074..._family1_healing_lock.md:117`,
  `179/sqrt(5)` -> `111/sqrt(5)`; and the collateral typo at
  `notes/stages/.../stage075..._family1_threshold_window.md:63`,
  `179/sqrt(5)` -> `111/sqrt(5)`.
- Script-side (SymPy):
  `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:70`
  now reads
  `expect_zero("alpha_ref - 111/sqrt(5)", alpha_ref - sp.Rational(111) / sp.sqrt(5))`,
  added immediately after the `kappa_ref - 12321/5` assertion.
- Script-side (Mathematica):
  `mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl:59`
  now reads
  `expectZero["alpha_ref - 111/sqrt(5)", alphaRef - 111/Sqrt[5]]`,
  added immediately after the `kappaRef - 12321/5` assertion.

The diff against both script files matches what the directive
prescribed exactly: one new assertion line per engine, no collateral
script edits.

**Assessment:**

The edit is correct and addresses the finding. Critically, the new
assertion is non-tautological: `alpha_ref` is derived inside each
script as `sqrt(kappa_ref)` where `kappa_ref` itself comes from
specializing the upstream derivation chain
`kappa_lock = 4*chi_lock**2 + (4/5)*Lambda_ell**2` (SymPy line 48 /
Mathematica `kappaLock` line 38) at `Lambda_ref = 37`. The literal
`111/sqrt(5)` in the assertion call is *not* used anywhere upstream of
`alpha_ref`, so the comparison `alpha_ref - 111/sqrt(5) == 0` is a real
test of the upstream chain -- if any step were wrong, the assertion
would fail.

Both engines' regenerated outputs confirm the assertion passes:

- SymPy
  (`scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt:19`):
  `alpha_ref - 111/sqrt(5) = 0`.
- Mathematica
  (`mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt:21-22`):
  `alpha_ref - 111/sqrt(5) = 0` followed by `PASS: alpha_ref - 111/sqrt(5)`.

All previously passing assertions (`chi_s - Lambda_ell/2`,
`kappa - (9/5) Lambda_ell^2`, `chi_ref - 37/2`, `kappa_ref - 12321/5`)
still pass on both engines, so no regression.

The paper/notes edits themselves cannot be re-verified from inside the
scripts-only scope, but the directive front-matter records that they
were applied; the script-side hard assertion now locks the engine
output to the canonical `111/sqrt(5)`, so any future drift on the
paper/notes prose away from that literal will be caught the next time
this stage is audited.

## Exec log assessment

**SymPy:** exit=0. No dedicated
`redteam/exec_logs/stage_074_sympy.log` was captured because the
orchestrator did not invoke the Codex `$RT exec-sympy` wrapper for this
resolution path (the fix was applied directly). The regenerated
`scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
serves as the run record. Notable lines:

> `chi_s - Lambda_ell/2 = 0`
> `kappa - (9/5) Lambda_ell^2 = 0`
> `chi_ref - 37/2 = 0`
> `kappa_ref - 12321/5 = 0`
> `alpha_ref - 111/sqrt(5) = 0`

The script reaches its trailing `Final ledger:` block, which the
`expect_zero` helper would have prevented via `AssertionError` if any
check had failed, so exit status is effectively 0.

**Mathematica:** exit=0. No dedicated
`redteam/exec_logs/stage_074_mathematica.log` was captured for the
same reason; the regenerated
`mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt`
is the run record. Notable lines:

> `PASS: chi_s - Lambda_ell/2`
> `PASS: kappa - (9/5) Lambda_ell^2`
> `PASS: chi_ref - 37/2`
> `PASS: kappa_ref - 12321/5`
> `PASS: alpha_ref - 111/sqrt(5)`
> `Stage 074 Mathematica audit passed.`

`Exit[0]` is reached on the final script line.

**Output freshness:** confirmed.

- SymPy `.py` mtime = 1779869150, `.txt` mtime = 1779869734 (txt newer
  by 584 s).
- Mathematica `.wl` mtime = 1779869153, `.txt` mtime = 1779869887 (txt
  newer by 734 s).

Both outputs are post-fix.

## Material-change assessment

`material_change`: false.

The only script-side change is an additional assertion that locks the
already-computed `alpha_ref` to its already-correct symbolic form
`111/sqrt(5)`. No constant, symbolic intermediate, or carry-forward
quantity changes value. Stage 075 already consumes the script-computed
`alpha`, not the (now-fixed) paper literal, so the paper/notes typo
fix likewise does not change any downstream numeric. No need to mark
units > 074 as upstream_stale.

## Side observations (non-blocking)

- The orchestrator's "user-resolved paper_misalignment" path skips
  capture of `redteam/exec_logs/stage_074_{sympy,mathematica}.log`
  because no Codex iteration runs. The verifier had to fall back on
  the regenerated `.txt` outputs to confirm exit status. Future
  resolutions on this path might either still capture a one-shot log
  or document explicitly in the directive that the `.txt` outputs are
  the substitute log record. This is not a verification blocker.
- The cosmetic "STAGE 57"/"STAGE 057" banner strings still appear in
  both script outputs (inherited from the pre-rename Stage-57 file).
  Already noted in the v1 verifier report; still not part of any
  finding.

## Verdict justification

The single F1 paper_misalignment was resolved as direction (a) and the
orchestrator's applied edits match the directive: paper line 31 and
both notes files updated to `111/sqrt(5)`; both engines gained a
non-tautological `alpha_ref - 111/sqrt(5) == 0` assertion that locks
the engine-computed `sqrt(kappa_ref)` to the canonical literal. Both
engine outputs are fresh and show every assertion passing, including
the new one. No regressions, no material change to any
downstream-relevant quantity. Verdict: `verified`.

stage 074: verified
