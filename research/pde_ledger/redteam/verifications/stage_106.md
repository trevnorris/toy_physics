---
unit_id: 106
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 106

## Per-finding outcomes

Note on directive bookkeeping: the directive file `redteam/directives/stage_106.md`
does not carry the standard `## Applied: F<n>` annotation blocks Codex usually
appends. The verifier nonetheless cross-checks the current state of the two
script files against the directive's "Required change" specifications, and the
human-provided close-out note in the verifier prompt enumerates exactly which
finding each edit addresses. All four findings can be reconciled from script
content alone.

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
- SymPy module docstring at
  `scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:1-13`
  was extended with a carry-forward annotation naming Stage 102
  (`higher_odd_irrelevance`) as the upstream verifier of Check (ii) and Stage 104
  (`outgoing_dtn_fingerprint`) as the upstream verifier of Check (iii) and the
  source of `chi_Q = 1`. The docstring explicitly limits the present stage's
  responsibility to Check (i) (factorization separability and closure to
  `N_Q = 1`).
- The corresponding annotation appears as a comment block in
  `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:29-32`
  ("Paper-card Checks (ii) and (iii) are exercised upstream at stages 102 and
  104; this engine uses chi_Q = 1 as carry-in …").
- No new SymPy or Mathematica `expect_zero` / `expectZero` calls were
  introduced specifically for F1 itself; the resolution is the user's chosen
  direction (b) — delegate Checks (ii) and (iii) to upstream stages 102 and 104.

**Assessment:**
The auditor's F1 was a `paper_misalignment` blocked on user resolution. The
prompt confirms the user chose option (b) (carry-forward) rather than option
(a) (extend scripts to own Checks (ii)/(iii)). Per the directive's "Resolve
before fix_loop" block, option (b) requires the paper card to delegate, and the
red-team is instructed not to edit paper.tex. The scripts have been annotated
with an explicit carry-forward note pointing at stages 102 and 104, which is
the appropriate script-side companion edit for the chosen direction.
Independently, the Mathematica script (lines 47-60) does in fact series-expand
`Y_Q^ret` to omega^7 and exhibits the next odd term at omega^7
(`(8 I/243) a^7/cS^7`) — this is a bonus script-side trace of Check (ii) that
does not contradict the carry-forward annotation, since the load-bearing
assertion ownership still lives in stage 102. No new fragile assertion was
added to either script for F1. The cross-stage chain is internally
consistent.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy
  `scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:64-71`
  now asserts
  `expect_zero("target identity K0_target K4_target - 4 K2_target^2", K0_target*K4_target - 4*K2_target**2)`
  and
  `expect_zero("target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3)", Gamma5_target - 9*sqrt(K2_target**5/K0_target**3))`,
  replacing the previous tautological `K4 - 4 K2^2/K0` and
  `Gamma5 - 9 sqrt(K2^5/K0^3)` checks.
- Mathematica
  `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:74-81`
  has the structurally equivalent `expectZero` calls on `k0Target*k4Target -
  4*k2Target^2` and `gamma5Target - 9*Sqrt[k2Target^5/k0Target^3]`.

**Assessment:**
Matches the directive's "After" patches verbatim. The new residuals depend on
the four hardcoded literal coefficients (`54/5, 6/5, 8/15, 2/5`) rather than on
the script-defined `K2 = K0/(4 OmegaQ^2)` chain, so substituting any one literal
with a wrong value would now cause the assertion to fail. The exec logs show
both new lines printing `0` and (in the Mathematica case) `PASS:`. Non-
tautological.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
was rewritten body-down. Compared to the SymPy script the choreography is now:

1. lines 39-41: `OmegaQ = 3 cS/(2 a)`, `sigmaQcan = (9/8)/OmegaQ^5`, with an
   independent sanity check `sigmaQcan - 4 a^5/(27 cS^5) == 0`.
2. line 45: defines `Yret[omSym_, chiSym_] := 3/4 + (1/4)/(1 - omSym^2/OmegaQ^2
   - I*chiSym*sigmaQcan*omSym^5)` — the appendix's retarded one-pole form.
3. lines 49-60: series-expands `Yret` to omega^7, extracts the omega^5/6/7
   coefficients, and asserts `omega^5 coefficient form = (I chiQ sigmaQcan)/4`.
4. lines 65-81: defines the four target literals and tests their mutual
   identity (the F2 fix).
5. lines 86-95: closes `N_Q = 1` via the source-map relation rather than via
   `Solve[constraint, NQ]`. Variable names `nqNatural`, `gamma5OnNatural`,
   `gammaEffCanonical` replace the SymPy `NQ_general`, `NQ_canonical`, `K0`,
   `K2`, `K4`, `Gamma5` chain.
6. lines 101-112: F4 sensitivity addition.

No `nqGeneral`, `k0`, `k2`, `k4`, or `gamma5` intermediates exist in the
rewritten file (the only `k0Target/k2Target/k4Target/gamma5Target` symbols are
hardcoded literal constants, not intermediate quantities derived from the
constraint). The bottom-line claim remains `N_Q = 1` and
`gamma_quad^eff = 2 G/(5 c^5)` as required.

**Assessment:**
The rewritten path uses a structurally different starting point (retarded one-
pole form + series expansion) from the SymPy path (algebraic constraint solve
+ chi_Q substitution). A bug in `Solve` or in the SymPy branch selection would
not be reproduced in the Mathematica path because the Mathematica path never
calls `Solve` and never picks a branch. The directive's "Required change"
items 1-7 are all present in the rewritten file. No transliteration concern
remains.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy
  `scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:79-96`
  defines `Delta_Q = sp.symbols("Delta_Q", real=True)`, substitutes
  `m0hat=1, chi_Q=1+Delta_Q` into `m0hat**2 * Gamma5`, series-expands to order
  1, and asserts `linear_coeff - (-2 G/(5 c^5)) == 0` and
  `zeroth_coeff - Gamma5_target == 0`. Code matches the directive's "After"
  block.
- Mathematica
  `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:101-112`
  has the structurally equivalent block on `gammaEffOff` / `gammaEffSeries`
  with the same two `expectZero` calls.

**Assessment:**
Both engines print the new pair of residuals (=0) per the exec logs
(`Delta_Q first-order sensitivity coefficient = 0` and `Delta_Q zeroth-order
coefficient equals Gamma5_target = 0`). The assertions are non-tautological:
the linear coefficient depends on the sign and form of `N_Q = 1/(chi_Q
m0hat^2)` (a sign flip or factor error in N_Q would change the slope), and the
zeroth coefficient depends on `Gamma5_target = 2G/(5 c^5)`. Spot-check:
`gamma_eff = m0hat^2 * (1/(chi_Q m0hat^2)) * Gamma5_target` at `m0hat=1` and
`chi_Q = 1 + Delta_Q` is `Gamma5_target/(1+Delta_Q) = Gamma5_target -
Gamma5_target * Delta_Q + O(Delta_Q^2)`, giving zeroth = `Gamma5_target`,
linear = `-Gamma5_target = -2 G/(5 c^5)`. Matches.

## Exec log assessment

**SymPy:** exit=0 (inferred; no traceback, no `AssertionError`, all
`expect_zero` lines print `= 0`, RESULT block emitted in full). Notable lines:

```
target identity K0_target K4_target - 4 K2_target^2 = 0
target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3) = 0
canonical gamma_eff - target = 0
Delta_Q first-order sensitivity coefficient = 0
Delta_Q zeroth-order coefficient equals Gamma5_target = 0
```

**Mathematica:** exit=0 (script ends with `Exit[0]`; no `FAIL:` lines; every
`expectZero` was followed by a `PASS:`). Notable lines:

```
PASS: omega^5 coefficient form
PASS: target identity k0Target * k4Target - 4 k2Target^2
PASS: target identity gamma5Target - 9 Sqrt[k2Target^5/k0Target^3]
PASS: Delta_Q first-order sensitivity coefficient
PASS: Delta_Q zeroth-order coefficient equals Gamma5_target
```

Engine cross-check: SymPy reports `K0 = 54 G c_s^5 /(5 a^5 c^5 chi_Q m0hat^2)`,
…, `Gamma5 = 2 G/(5 c^5 chi_Q m0hat^2)`. Mathematica derives `gamma_eff =
2 G/(5 c^5)` on the canonical branch and the same Delta_Q slope of
`-2 G/(5 c^5)`. The two engines reach the same bottom line via structurally
different paths (constraint-solve vs. retarded one-pole series expansion), so
this engine agreement is now substantive evidence.

**Output freshness:** confirmed — SymPy log mtime `2026-05-27 15:18:01` is
later than its script mtime `15:13:22`; Mathematica log mtime
`2026-05-27 15:24:11` is later than its script mtime `15:14:04`. Both outputs
were regenerated after the post-fix script edits.

## Material-change assessment

`material_change`: false.

The fixes consist of: (i) docstring annotations and one structural rewrite of
the Mathematica audit whose bottom-line printed values are unchanged
(`N_Q = 1`, `gamma_quad^eff = 2 G/(5 c^5)`); (ii) two replaced assertions in
each engine that now test the four `*_target` literals' mutual algebraic
consistency rather than a definitional identity (residual is still
identically zero in both forms); (iii) two new sensitivity assertions whose
correctness is implied by the existing closure and which do not export any
new symbol or change any downstream-visible constant. No
constant, no derived expression, no carry-forward symbol was modified.
Downstream units that import `chi_Q = 1`, `N_Q = 1`, or
`gamma_quad^eff = 2 G/(5 c^5)` see no change.

## Side observations (non-blocking)

- The SymPy script's `K0/K2/K4/Gamma5` intermediate variables (lines 54-62) are
  retained from the pre-fix structure and are no longer load-bearing for any
  assertion now that F2 replaced the branch-identity checks with target-only
  checks. They are still printed for narrative purposes; this is not a defect,
  but they could be removed in a future cleanup. Not blocking.
- The Mathematica script's omega^6 coefficient (`(16 a^6)/(729 cS^6)`) is
  printed but not asserted on. It is mentioned for completeness; the load-
  bearing odd-coefficient assertion is on omega^5. Acceptable.
- The auditor's secondary observation about stale banner labels in the pre-fix
  scripts is now moot — both scripts banner `STAGE 106 — REDUCED 2.5PN
  CLOSURE ON CANONICAL OUTGOING DtN BRANCH` correctly.

## Verdict justification

All four auditor findings have been addressed correctly. F1 was resolved by
the user's chosen carry-forward direction (option (b) in the directive) and is
documented in both engines' header comments naming stages 102 and 104 as the
upstream owners of Checks (ii) and (iii). F2's replacement assertion now tests
the mutual consistency of the four hardcoded target literals and would fail
under perturbation of any one of them. F3's Mathematica rewrite uses a
structurally independent retarded one-pole + series-expansion path with no
shared intermediate names against the SymPy file, so the engine cross-check is
no longer hostage to transliteration. F4's Delta_Q sensitivity adds two non-
tautological assertions in each engine whose residuals depend on the sign and
form of `N_Q` and on `Gamma5_target`. Both exec logs run clean. No regression
is visible and no downstream-visible constant changed; `material_change` is
false.
