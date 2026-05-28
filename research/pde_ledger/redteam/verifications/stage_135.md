---
unit_id: 135
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 135

This batch was orchestrator-direct (Codex bypassed); the directive has no
`## Applied:` blocks, so verification proceeds by reading the current state of
the edited files. Findings F1 and F2 in the directive are tightly coupled and
were resolved by a single SymPy edit batch; the Mathematica script was not
touched (it already contained all six checks).

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
In `scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py`,
the original conditional `raise` block

```python
if abs(residual) > 1e-12:
    raise AssertionError("Outlet-consistent threshold did not close.")
```

was removed. The residual is still computed and printed for the transcript
(lines 85-86) as `Pi_star - Sigma_star*(4 - S_star) = ...`. The Mathematica
script was not edited (per directive); its `expectApprox["closure residual",
residual, 0, 10^-14]` at line 78 remains as a shadowed check (substantive
content lives elsewhere in the .wl per the original report).

**Assessment:**
Correct. The tautological gate is gone, and the residual print remains as a
numerical sanity probe (`-7.22e-17` in the sympy output, exactly as the
original audit reported). The load-bearing assertions are now the anchored
ones added under F2.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
The SymPy script now contains five substantive assertions plus the mixed-lane
anchor:

1. Line 40: `assert residual_sub == 0` — symbolic substitution check
   (`M_s + M_q*S_q - Sigma_m*(4 - S_q) == 0` after substituting
   `M_s -> 4*Sigma_m, M_q -> -Sigma_m`). The output prints
   `M_s + M_q*S_q - Sigma_m*(4 - S_q) = 0` confirming the assertion fires
   on actual symbolic content (not a no-op).
2. Line 58: `assert s_in_range` — `0 < S_q(Pi_*) < 1` (output: `True`).
3. Line 69: `assert abs(Sigma_star - Sigma_target) < 1e-12` — anchors
   `Sigma_m^* ≈ 0.451485277739090` (output shows
   `0.451485277739089696513730132210`, diff ~1.1e-16).
4. Line 71: `assert abs(M_s_star - M_s_target) < 1e-11` — anchors
   `M_s^* ≈ 1.80594111095636`.
5. Line 73: `assert abs(M_q_star - M_q_target) < 1e-12` — anchors
   `M_q^* ≈ -0.451485277739090`.
6. Line 81: `assert abs(mixed_correction - mixed_target) < 1e-11` — anchors
   mixed-lane correction `M_q^* * S_q(Pi_*) ≈ -0.297111597463199` (output:
   `-0.297111597463198745446779533971`).

**Assessment:**
All six claim-manifest items (M1-M6) are now asserted. The targets come from
the notes file (lines 86-101, 124-127 per original audit); both engines agree
to ~14 sig figs on every quantity. The substitution check in step 1 is
non-tautological because it uses independent symbols `M_s_sym, M_q_sym`
substituted into the generic two-channel law and then compared against the
reduced form — exactly what the original report and directive prescribed
(self-test note 3 of the audit confirmed this would pass under any `S_q`).

Cluster A rename verified independently: SymPy banner at line 27 reads
`STAGE 135 — OUTLET-CONSISTENT ONE-PARAMETER CLOSURE`; Mathematica banner at
line 38 matches; notes H1 reads `Stage 135` (not 237). The 118/237 cosmetic
labels flagged in the original audit's closing note are now consistent across
all three artifacts.

## Exec log assessment

**SymPy:** exit=0 (script ran cleanly; no AssertionError raised; output ends
after the final residual print, which is the last statement in the script).
Notable lines:
- `M_s + M_q*S_q - Sigma_m*(4 - S_q) = 0` (M1 confirmed)
- `0 < S_q(Pi_star) < 1 -> True` (M2 confirmed)
- `Sigma_m^* = 0.451485277739089696513730132210` (M3)
- `M_s^* = 1.80594111095635878605492052884` (M4)
- `M_q^* = -0.451485277739089696513730132210` (M5)
- `Sigma_m^*, M_s^*, M_q^* anchored to notes values within tolerance.`
- `M_q^* * S_q(Pi_*) = -0.297111597463198745446779533971` (M6)
- `Pi_star - Sigma_star*(4 - S_star) = -7.22456073765234902012354489563e-17`
  (residual still printed)

**Mathematica:** exit=0. All seven explicit PASS lines present
(`outlet-consistent reduction`, `0 < S_q(Pi_star) < 1`,
`Sigma_m^* numeric check`, `M_s^* numeric check`, `M_q^* numeric check`,
`mixed-lane correction numeric check`, `closure residual`,
`3 Sigma_m^* < Pi_* < 4 Sigma_m^*`).

**Output freshness:** confirmed. SymPy script mtime 1779925397 < output mtime
1779925552; Mathematica script mtime 1779924055 < output mtime 1779925669.
Both `.txt` outputs were regenerated after their corresponding script edits.

## Engine cross-check (numerical agreement, unchanged from audit)

| Quantity | SymPy | Mathematica | Agreement |
|---|---|---|---|
| `S_q(Pi_*)` | 0.658075937605428494... | 0.658075937605429486... | ~14 sf |
| `Sigma_m^*` | 0.451485277739089696... | 0.451485277739089808... | ~15 sf |
| `M_s^*` | 1.80594111095635878... | 1.80594111095635923... | ~14 sf |
| `M_q^*` | -0.451485277739089696... | -0.451485277739089808... | ~15 sf |
| mixed corr | -0.297111597463198745... | -0.297111597463199267... | ~14 sf |
| residual | -7.22e-17 | 0`` (precision-limited) | both ~0 |

Targets in both engines come from the notes file at the same digit string, and
both pass their respective tolerances by ~3+ orders of magnitude.

## Material-change assessment

`material_change`: false. No derived numerical values changed; the SymPy
script's numerical pipeline was untouched. The edits added assertions and
removed a tautological gate; they do not alter any quantity that downstream
units consume. No additional downstream concern beyond the routine
orchestrator `upstream_stale` flag.

## Side observations (non-blocking)

- The SymPy script declares `Sigma_var = sp.symbols("Sigma_var", positive=True, real=True)`
  at line 48 to solve for `Sigma_star`. This shadows the original
  one-parameter idea (since `Sigma_var` is now a fresh symbol) but is
  consistent with how the closed form is recovered numerically; not a finding.
- The `M_s_sym, M_q_sym` symbols at line 34 are declared `real=True` (not
  `positive=True`) — correct per the directive, because `M_q^* < 0`.

## Verdict justification

Both findings are resolved by a single coherent SymPy edit batch. The
tautological `if abs(residual) > 1e-12` gate is gone (F1); five substantive
assertions covering all six claim-manifest items (M1-M6) are now present and
each prints a non-trivial value in the output (F2). The Mathematica script is
unchanged (correctly, per directive). Both engines exit 0, agree to ≥14 sig
figs on every load-bearing quantity, and the outputs are fresh relative to
the scripts. Cluster A `STAGE 135` rename is verified across `.py`, `.wl`,
and notes H1.
