---
unit_id: 074
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 074

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:29-41` (per
the diff at `redteam/exec_logs/stage_074_diff.patch`). The previous block

```
chi_s = sp.symbols("chi_s", positive=True)
kappa = sp.symbols("kappa", positive=True)
chi_lock = sp.simplify(Lambda_ell / 2)
```

was replaced with the directive's substitution chain: declare positive symbols
`hbar, m_psi, c_s, ell, L`; define `chi_def = m_psi * c_s * L / hbar` (line 33);
apply `subs(c_s, hbar/(2*m_psi*ell))` to get `chi_in_ell = L/(2*ell)` (line 37);
print the intermediate; then apply `subs(L, Lambda_ell*ell)` to produce
`chi_lock = Lambda_ell/2` (line 41). The unused literal `chi_s`/`kappa` symbol
declarations were removed (they were only definitional and not referenced
downstream). The redundant first `chi_def` line shown in the directive's
illustrative snippet was correctly omitted — only the single canonical
assignment exists at line 33.

**Assessment:**
Edit is correct and matches the directive exactly. The assertion at line 53
(`chi_s - Lambda_ell/2`) is now non-tautological: `chi_lock` arrives via the
substitution chain `m_psi*c_s*L/hbar -> L/(2*ell) -> Lambda_ell/2`. If any link
in the chain were wrong, the assertion would fail. The output transcript
confirms the intermediate print `chi (after healing-length substitution) =
L/(2*ell)` and the final `chi_s (locked) = Lambda_ell/2`, both expected.
Reference branch lines 56-58 still derive `chi_ref` and `kappa_ref` from the
derived `chi_lock`/`kappa_lock` via subs(Lambda_ell, 37), so they inherit the
non-tautological derivation as the directive intended. No collateral edits
beyond the symbol-declaration cleanup that was required for the new chain.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:43-47`
now contains the required four-line provenance comment immediately above the
`kappa_lock = sp.simplify(4 * chi_lock**2 + sp.Rational(4, 5) * Lambda_ell**2)`
definition (line 48). The comment text matches the directive verbatim.

**Assessment:**
Correct. The comment anchors the `4` and `4/5` coefficients to the Family-1
Euler-Lagrange branch from earlier stages and clarifies the scope of this
stage's check. Combined with the F1 chain, the kappa assertion at line 54 is
now a meaningful test that the locked `chi_s = Lambda_ell/2` substituted into
the Family-1 branch formula yields `(9/5) Lambda_ell^2`. The assertion still
passes in the refreshed output.

## Exec log assessment

**SymPy:** exit=0 (inferred). The captured exec log file
`redteam/exec_logs/stage_074_sympy.log` is absent from the exec_logs
directory; only `stage_074_diff.patch` is present. However, the refreshed
output transcript at
`scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
shows all six expected lines and four `= 0` assertion outputs without any
`AssertionError` trace, which the `expect_zero` helper would have raised on
failure. Output mtime is May 22 23:12, newer than the script mtime May 22
23:11, confirming the post-fix run. Notable lines:

```
chi (after healing-length substitution) = L/(2*ell)
chi_s (locked) = Lambda_ell/2
kappa(Lambda_ell) = 9*Lambda_ell**2/5
chi_s - Lambda_ell/2 = 0
kappa - (9/5) Lambda_ell^2 = 0
chi_ref - 37/2 = 0
kappa_ref - 12321/5 = 0
alpha (numeric) = 49.640709100495331260
```

**Mathematica:** exit=0 (inferred). Log file
`redteam/exec_logs/stage_074_mathematica.log` is also absent, but the
Mathematica script was untouched by this directive (F1 and F2 both target the
SymPy file only), so a re-run was not strictly required. The Mathematica
output at
`mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt`
shows four `PASS:` lines and the trailing `Stage 074 Mathematica audit
passed.` banner; mtime May 22 23:12. Notable lines:

```
PASS: chi_s - Lambda_ell/2
PASS: kappa - (9/5) Lambda_ell^2
PASS: chi_ref - 37/2
PASS: kappa_ref - 12321/5
Stage 074 Mathematica audit passed.
```

**Output freshness:** SymPy output mtime (May 22 23:12) > SymPy script mtime
(May 22 23:11). Mathematica output mtime (May 22 23:12) > Mathematica script
mtime (May 11 11:56). Both outputs are post-fix.

## Material-change assessment

`material_change`: false.

The edit only changes the derivation route for `chi_lock`; the final symbolic
and numeric values are identical to the pre-fix audit (compare row-by-row in
the original report's engine-cross-check table — every printed quantity
matches: `chi_s = Lambda_ell/2`, `kappa = 9*Lambda_ell**2/5`, `chi_ref = 37/2`,
`kappa_ref = 12321/5`, `alpha = 111*sqrt(5)/5`, `alpha (numeric) =
49.640709100495331260`). No downstream-quoted constant or symbolic result has
changed. This is a provenance-only fix: the SymPy script now derives the value
the Mathematica script already derived, rather than declaring it.

## Side observations (non-blocking)

- The exec log files `stage_074_sympy.log` and `stage_074_mathematica.log` are
  not present in `redteam/exec_logs/`. I inferred exit=0 from the refreshed
  output transcripts (the `expect_zero` helper raises `AssertionError` on
  failure, which would have prevented the trailing `Final ledger:` print on
  the SymPy side and the trailing `Stage 074 Mathematica audit passed.` line
  on the Mathematica side; both are present). Orchestrator may want to confirm
  log-capture behaviour, but this does not impede verification.
- The banner string on SymPy line 27 still reads "STAGE 57" (and Mathematica
  output line 3 reads "STAGE 057"), inherited from the original Stage-57 file
  the script was renamed from. The docstring header on line 3 of the SymPy
  file also still names `stage57_family1_healing_lock_sympy_audit.py`. Cosmetic
  only; not part of either finding.

## Verdict justification

Both findings were applied exactly as the directive specified, with no
deviation and no collateral edits beyond the symbol-declaration cleanup
necessitated by F1's substitution chain. The refreshed outputs confirm all
four assertions still pass on both engines, and all printed quantities are
identical to the pre-fix values — the fix changes provenance, not results.
Verdict is `verified` with `material_change: false`.

stage 074: verified
