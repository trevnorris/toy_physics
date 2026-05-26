---
unit_id: 033
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 033

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:125-128`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:128-131,146-160`

**Issue:**
The Stage 16.6 final identity assertions in both engines are tautological by construction. `gate_num_target` (sympy line 124) and `gateNumTarget` (mma line 126) are defined as `gate_num_actual / den_ratio`, where `den_ratio = gate_den_actual / gate_den_claim`. Therefore `gate_num_target / gate_den_claim = gate_num_actual / gate_den_actual = alpha_crit_mic - alpha0_mic` identically. The labeled final `expect_zero` / `expectZero` is decorative — the substantive verification is in the `den_ratio.is_number` guard (sympy line 120-122) / `NumericQ[denRatio]` guard (mma line 122-124), which DO catch a structurally wrong claimed denominator. The Mathematica numerical cross-check for the gate identity at lines 153-158 inherits the same tautology: substituting concrete rationals into a symbolically-zero expression still yields zero, so it adds no coverage. (The monotonicity numerical check at lines 147-152 IS substantive and must be kept.)

The fix re-labels the decorative assertion to make its tautological status explicit, and drops the redundant numerical gate-identity check (keeping the monotonicity numeric check).

**Required change:**

### SymPy file edit

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`, edit the `expect_zero` call at lines 125-128. Change only the label string; do not change the residual expression.

Before (lines 125-128):
```python
expect_zero(
    "alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim",
    sp.simplify(alpha_crit_mic - alpha0_mic - gate_num_target / gate_den_claim),
)
```

After:
```python
expect_zero(
    "gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) "
    "(tautological by reconstruction; substantive check is den_ratio.is_number above)",
    sp.simplify(alpha_crit_mic - alpha0_mic - gate_num_target / gate_den_claim),
)
```

### Mathematica file edit (label change)

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl`, edit the `expectZero` call at lines 128-131. Change only the label string.

Before (lines 128-131):
```
expectZero[
  "alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim",
  alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim
];
```

After:
```
expectZero[
  "gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) (tautological by reconstruction; substantive check is NumericQ[denRatio] above)",
  alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim
];
```

### Mathematica file edit (drop redundant gate numeric check)

Also in the same `.wl` file, in the `Do[...]` block at lines 146-160, remove the gate-identity numeric residual block (lines 153-158 inclusive). Keep the monotonicity-numeric block and the loop iterator unchanged.

Before (lines 146-160):
```
Do[
  monotonicityNumeric = N[(dN - dNFormula) /. rule, 30];
  Print["monotonicity numeric residual = ", fmt[monotonicityNumeric]];
  If[Abs[monotonicityNumeric] > 10^-20,
    fail["monotonicity numeric residual nonzero", monotonicityNumeric],
    pass["monotonicity numeric residual zero at rule"]
  ];
  gateNumeric = N[(alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim) /. rule, 30];
  Print["gate-identity numeric residual = ", fmt[gateNumeric]];
  If[Abs[gateNumeric] > 10^-20,
    fail["gate-identity numeric residual nonzero", gateNumeric],
    pass["gate-identity numeric residual zero at rule"]
  ],
  {rule, {numericRule1, numericRule2}}
];
```

After:
```
Do[
  monotonicityNumeric = N[(dN - dNFormula) /. rule, 30];
  Print["monotonicity numeric residual = ", fmt[monotonicityNumeric]];
  If[Abs[monotonicityNumeric] > 10^-20,
    fail["monotonicity numeric residual nonzero", monotonicityNumeric],
    pass["monotonicity numeric residual zero at rule"]
  ],
  {rule, {numericRule1, numericRule2}}
];
```

**Self-test:**

1. Variable independence: the edits do not introduce new `sp.diff` / `D[...]` calls. The retained `monotonicityNumeric = N[(dN - dNFormula) /. rule, 30]` already has `alpha0` in the differentiated expressions (`dN = D[nMinus, alpha0]`); unchanged.
2. Symmetry/parity: no integrals; N/A.
3. Trivial-case pre-check: the monotonicity numeric check evaluates `(dN - dNFormula)` at concrete rules. The SymPy/Mathematica analytic checks have already passed (residual 0 in transcripts), so at any rule the numeric residual stays 0 to machine precision. Unchanged.
4. Path specifications: edits target the existing `.py` in `scripts/` and `.wl` in `mathematica/`; no new files.
5. Paper round-trip: the label changes are cosmetic-script-side only; no paper-side claim is touched, no new `paper_misalignment` risk.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 033` and `redteam exec-mathematica 033` and confirm:
(a) the sympy `expect_zero` at line 125-128 now uses the longer label including the substring `"tautological by reconstruction"`;
(b) the mma `expectZero` at line 128-131 now uses the analogous longer label;
(c) the mma `Do[...]` block no longer contains a `gateNumeric` block;
(d) both scripts still exit 0;
(e) the SymPy transcript still shows `"PASS"`-equivalent (zero residual) for the (now-relabeled) gate identity, and the Mathematica transcript still shows `PASS: monotonicity numeric residual zero at rule` twice but no longer shows `gate-identity numeric residual`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl`
- summary: Relabeled the reconstructed gate identity as tautological and removed the redundant Mathematica gate numeric residual check.
- deviation: none
