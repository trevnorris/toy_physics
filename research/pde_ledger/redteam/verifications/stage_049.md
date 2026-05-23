---
unit_id: 049
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T17:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 049

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py:65-68`
  ```python
  expect_zero(
      "k_n satisfies D/N Neumann boundary",
      sp.cos(k_n * L),
  )
  ```
- `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:50`
  ```wolfram
  expectZero["k_n satisfies D/N Neumann boundary", Cos[kN l]];
  ```

**Assessment:**
Both engines now check `cos(k_n · L)` rather than `k_n − (n+1/2)π/L`. Because `n` is declared an integer in both scripts (`n = sp.symbols("n", integer=True, nonnegative=True)` in sympy; `Element[{n}, Integers] && n >= 0` in `$Assumptions` in Mathematica), `cos(k_n·L) = cos((n+1/2)π)` simplifies to `0` non-trivially — the assertion now exercises the Neumann boundary condition rather than a definition. The fresh outputs confirm the residual is `0` (sympy line 7; mathematica lines 7-8 with `PASS`). No occurrence of `k_n - (n+1/2) pi / L` (or `kN - (n + 1/2) Pi/l`) remains in either script. The surrounding `chi_n`/`chiN` definitions and prints are untouched, matching the directive exactly.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`: helper definition `uniformDnOverlap[n_, l_] := ...` deleted (was at line 27).
- Line 53 now reads `overlapFormula = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];`
- Line 58 now reads `i0 = FullSimplify[overlapFormula /. n -> 0];`

**Assessment:**
The Mathematica overlap derivation now routes through `Integrate[chiN, {s, 0, l}]` with an explicit integer assumption, rather than echoing the SymPy closed form via a parallel helper. `i0` is now a substitution `n -> 0` on the integral-derived `overlapFormula`. The SymPy file's helper `uniform_dn_overlap` is correctly left intact (only the `.wl` mirror was to change). Output shows `I_n closed form = (2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` — algebraically equivalent to the SymPy `2*sqrt(2)*sqrt(L)/(pi*(2*n + 1))` — and the `uniform overlap integral` residual is `0` with `PASS`. `grep`-equivalent inspection confirms `uniformDnOverlap` no longer appears anywhere in the `.wl`. `overlapRatio`, `twinSupportRatio`, the zeta-twin assertions, and the carry-forward print block are all untouched, as required.

## Exec log assessment

Per the verifier prompt's path overrides, the canonical fresh transcripts are the saved `.txt` outputs (the `$RT exec-*` MANIFEST race meant `redteam/exec_logs/stage_049_*.log` are absent or stale; the saved outputs are the authoritative source).

**SymPy:** exit=0 (inferred — script ran to completion; no `AssertionError`; all seven `expect_zero` calls printed `= 0` and the final "Carry-forward formulas" block was emitted). Notable lines:
- `k_n satisfies D/N Neumann boundary = 0`
- `uniform overlap integral = 0`
- `lowest twin half-wave = 0`

**Mathematica:** exit=0 (script's terminal `Exit[0]` reached; "Stage 32 Mathematica audit passed." printed). Notable lines:
- `k_n satisfies D/N Neumann boundary = 0` / `PASS: k_n satisfies D/N Neumann boundary`
- `I_n closed form = (2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` / `PASS: uniform overlap integral`
- `PASS: lowest twin half-wave`

**Output freshness:** confirmed. Script mtimes vs output mtimes:
- `.py` 2026-05-22 16:55 → `.txt` 2026-05-22 16:56 (output newer)
- `.wl` 2026-05-22 17:24 → `.txt` 2026-05-22 17:25 (output newer)

Both outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false.

No derived results (closed-form `I_n`, ratio `I_n/I_0`, `zeta_n^(phys)`, twin-stiffness identity, `zeta_n^(twin)`, n=0 limit) changed. F1 only replaced a tautological check with a substantive one; the value of `k_n` itself is unchanged. F2 only swapped the Mathematica derivation path for the overlap (helper → integrator), and the resulting closed form is algebraically identical to the prior helper output. Downstream units that consume `I_n`, `I_n/I_0`, `zeta_n^(phys)`, or `zeta_n^(twin)` see no change.

## Side observations (non-blocking)

- The SymPy script's `overlap_from_integral` is computed from `sp.integrate(chi_n, (s, 0, L))` (line 70) and compared against the helper closed form (line 74), so the SymPy side already has independent integrator-vs-closed-form coverage; the F2 change brings the Mathematica side to a comparable level (integrator-derived `overlapFormula`, with `overlapFromIntegral` still computed separately at wl:52). The Mathematica `overlapFromIntegral - overlapFormula` assertion at wl:56 is now structurally `Integrate[...] - Integrate[...]` (both with similar assumptions), so it largely re-checks idempotence of `FullSimplify` rather than integrator-vs-closed-form; this is acceptable given F2's primary intent (break the SymPy mirror), and the SymPy script retains the integrator-vs-closed-form cross-check at line 74. Flagging only for awareness — not a verification issue.
- The banner still reads "STAGE 32" in both scripts despite the file name being `stage049`. Pre-existing; not in scope for this directive.

## Verdict justification

Both findings (F1, F2) are fully resolved: the tautological half-wave momentum check is replaced by a non-trivial Neumann boundary residual in both engines, and the Mathematica overlap derivation now flows through `Integrate` with an integer assumption rather than echoing the SymPy helper closed form. Fresh outputs show `= 0` (and `PASS` in Mathematica) for every assertion, output mtimes post-date script mtimes, the SymPy script's untouched regions remain untouched, and no derived downstream quantity changed. Verdict: verified.
