# Measurements — O2 comparator decision-list review dispositions (generated 2026-10-07 13:37)

Generator: `_scratch/s9b_build/gen/o2_comparator_directive_lookups.sh` (grep/sed lookups only; lines cut at 300 characters). Directive as reviewed: sha256 below.

```
$ sha256sum directives/O2_comparator_build_directive.md
64b79f55898736dd0cda65adfd002d3951c9ebe2d57e4d6b670f27241e62162f  directives/O2_comparator_build_directive.md
```

## U1 — unjoin rule; mass residual in both streams; Wolfram heads not named in the format note
```
$ sed -n '19,21p;31,36p' directives/O2_comparator_build_directive.md
  - `research/pde_ledger_v3/scripts/out/O2_live_balance_sympy_audit.out` (SymPy: one line per tag, `srepr` payloads);
  - `research/pde_ledger_v3/mathematica/out/O2_live_balance_mathematica_audit.out` (Wolfram: associations with
    `OPEN`, `OpenAction` and `Inactive` heads).
1. **Object-level join, declared and injective.** The two streams cut the same physics at different grain, so
   tags cannot be joined by name. The comparator carries an explicit join table in its source. Each row pairs a
   SymPy object (a tag, or a component or entry inside one) with a Wolfram object (a tag plus key path). Each row
   cites the construction line in each engine that produces its object. The table is injective: no object is on
   two rows. Every emitted object in either stream is either on a row or listed as unjoined with its reason.
   Nothing is skipped silently.
```

```
$ grep -n 'mass_residual =\|.MASS_RESIDUAL.: mass_residual' scripts/O2_live_balance_sympy_audit.py
267:    mass_residual = mass_equation.lhs - mass_equation.rhs
380:        'MASS_RESIDUAL': mass_residual,
```

```
$ grep -n '"Residual" ->' mathematica/O2_live_balance_mathematica_audit.wl
247:  "Equation" -> (massDivergence == massRHS), "Residual" -> (massDivergence - massRHS),
```

```
$ grep -o 'OpenFirstVariation\|OpenNativeField\|UnrestrictedNativeSection\|UnrestrictedSection\|OpenAction\|Inactive' mathematica/out/O2_live_balance_mathematica_audit.out | sort | uniq -c
    254 Inactive
   1277 OpenAction
     93 OpenFirstVariation
    326 OpenNativeField
   1309 UnrestrictedNativeSection
    488 UnrestrictedSection
```

## U2 — the named reusable residual; the SymPy payload printer
```
$ sed -n '1,8p;41p;628,630p;670,672p;816,820p' scripts/S11b_cross_engine_comparator.py
#!/usr/bin/env python3
"""Compare PY/WL S11b tag streams without interpreting their physics.

The comparator joins non-local tags by emitted object name, parses each
payload, applies only injective mechanical symbol/function-head transliteration, and prints
the operands, computed residual, and one of AGREE, DISAGREE, UNDECIDED, or
UNCOMPARED.  A disagreement is a reported result; only operational failures
make the process exit nonzero.
DEFAULT_RESIDUAL_LEAF_BUDGET_SECONDS = 5.0
def _canonical_basic_same(left: sp.Basic, right: sp.Basic) -> bool:
    """Compare non-residualable SymPy structures without object equality."""
    return type(left) is type(right) and sp.srepr(left) == sp.srepr(right)

    signal.signal(signal.SIGALRM, exceed_budget)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    if isinstance(py_value, TextAtom) or isinstance(wl_value, TextAtom):
        if isinstance(py_value, TextAtom) and isinstance(wl_value, TextAtom):
            if py_value.value == wl_value.value:
                return sp.S.Zero
            return Mismatch("TOKEN_DISAGREE", py_value, wl_value)
```

```
$ grep -n 'class LosslessReprPrinter' -A4 scripts/O2_live_balance_sympy_audit.py
60:class LosslessReprPrinter(ReprPrinter):
61-    """Preserve the undefined-function constructor metadata omitted by srepr.
62-
63-    In particular, a real native immersion function must not revive as an
64-    assumption-free function. This changes encoding, not any constructed object
```

```
$ grep -o "Function('o2_X_face_0', \*\*{" scripts/out/O2_live_balance_sympy_audit.out | wc -l
4440
```

## U3 — controls and the no-target rule
```
$ sed -n '59,79p' directives/O2_comparator_build_directive.md
   Then print the differences. A difference in representation (for example, native face geometry as a general
   immersion in one engine and an OPEN `𝒥_map` action in the other) is printed as a difference. It is never
   reconciled inside the comparator.

5. **No target.** Nothing in the comparator, its tests or its report states what any residual or comparison is
   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
   and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output belongs
   to the record (sub-step 7).

6. **Controls that can fail.** The tests use synthetic fixtures only and never load the measured streams.
   They show:
   - **One-sided corruption:** altering one engine's operand moves that row's residual away from what the
     unaltered pair prints.
   - **Repoint:** item 2.
   - **Form:** replacing an operand by one of a different functional form, not a rescaling, moves the residual.
   - **Injectivity:** a duplicate join or name row is rejected.
   - **Three-valued handling:** a native boolean is rejected as an operand while a sibling algebraic leaf is
     still subtracted.

   Run against the real streams, the comparator also prints per-row accounting: joined, unjoined with reason,
   and the count of extracted leaves. An under-measured row stays visible.
```

## Fold applied (generated 2026-10-07 13:39; generator `_scratch/s9b_build/gen/o2_comparator_directive_fold.sh`)

```
$ sha256sum _scratch/s9b_build/O2_comparator_build_directive_reviewed_v0.md research/pde_ledger_v3/directives/O2_comparator_build_directive.md
64b79f55898736dd0cda65adfd002d3951c9ebe2d57e4d6b670f27241e62162f  _scratch/s9b_build/O2_comparator_build_directive_reviewed_v0.md
97260120baa6b61269ac466919099246b5ec7016a08d8c33aa78240cd528ed6e  research/pde_ledger_v3/directives/O2_comparator_build_directive.md
```

```
$ diff _scratch/s9b_build/O2_comparator_build_directive_reviewed_v0.md research/pde_ledger_v3/directives/O2_comparator_build_directive.md | grep -c '^[<>]'
77
```

