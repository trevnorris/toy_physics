# Measurements behind `S11_sympy_build_directive.md`

Commands and their literal output. Rule 2: a claim about an artifact carries the command that produced
it. Regenerate with the commands as written; do not transcribe.

## The two inherited facts the directive cites about the S10 exporter

```
$ sed -n '1906,1910p;2087p' scripts/S10_brane_mode_spectrum_sympy_audit.py
def exact_reconstruction_match(live_value: object, reconstructed_value: object) -> bool:
    return (
        type(live_value) is type(reconstructed_value)
        and sp.srepr(live_value) == sp.srepr(reconstructed_value)
    )
            not record.overwrites_upstream or upstream["step"] != "S9"
```

## `§7` names the primary package — the directive points here rather than naming it

```
$ grep -n "primary result of the step" directives/S11_SHARED_PHYSICS.md
1036:⭐ **`MAIN` swept over `D` is the primary result of the step, ⛔ not a control.**
```

## `F1` is live and unsuperseded — it is what constrains the tag→key spelling

```
$ sed -n '16,20p' directives/S11_export_chain_decisions_v2.md
**F1 · Storage keys stay FLAT. `D5` is unchanged; ⛔ there is no producer prefix.**
⭐ A step writes the object's name. A later step re-deriving that object writes **the same key**.
⚠⚠ **⛔ AMENDED BY `F9`:** a producer prefix **is** written when the override would be wrong. ⛔ Read `F9`.

**F2 · Before writing a key that exists in the imported `LEDGER`, the writer compares the OBJECT.**
```

## `§Q6r`'s lookup is two-step — what "the next step can import" has to mean

```
$ grep -n "dimension_key" directives/S11_SHARED_PHYSICS.md | head -4
567:`LEDGER[name]['dimension_key']`, then `LEDGER[that key]['value']` lookup. Once per `(package, D)`, emit
```

## What S10's export population ACTUALLY is — two sources, not one

```
$ python3 -c "<census of key shapes and class carriers in S10_exports.py>"
total rows      : 617
tag-shaped keys : 0
KNOB rows       : ['rho_br', 'mu_R']
STRUCTURAL rows : ['D', 'length_dimension', 'time_dimension', 'mass_dimension']
```

```
$ grep -n "def add_symbol_records" -A 4 scripts/S10_brane_mode_spectrum_sympy_audit.py
2005:def add_symbol_records(merged: dict[str, dict[str, object]]) -> None:
2006-    occurrence_classes: dict[str, set[str]] = {}
2007-    occurrence_objects: dict[str, dict[str, sp.Symbol]] = {}
2008-    for record in merged.values():
2009-        value = record["value"]
```

⇒ the export is tag-derived rows PLUS a row per free symbol occurring in them.
⇒ the KNOB rows and `D` are bare-symbol rows from the second source; the three unit-basis markers are
  matrices carried from the predecessor's TAG source. ⛔ Corrected 2026-08-12: an earlier line here said
  "every KNOB and STRUCTURAL row" and was over-general — ⚠ the numbers above were exact, the inference
  was not, and a review leg caught it. ⭐ An inference is not a measurement.
