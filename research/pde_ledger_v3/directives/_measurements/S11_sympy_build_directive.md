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
