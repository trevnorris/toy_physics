# _measurements — S11_ledger_fold_repair_directive.md

Folds the two module-review legs (fresh Claude + Grok) against the committed baseline; each claim carries its
command (rule 2). Physics-free infra repair.

## baseline committed 05f86947 (reviewed module before repair)
```
$ git -C /var/projects/toy_physics log --oneline -1 05f86947
05f86947 Ledger fold+guard module (reviewed baseline before repair)
```

## both legs' evidence saved beside the review prompt
```
$ ls -la research/pde_ledger_v3/directives/_legs/ledger_fold_review*.* | awk '{print $5"  "$9}'
1049  research/pde_ledger_v3/directives/_legs/ledger_fold_review_agent.out
10346  research/pde_ledger_v3/directives/_legs/ledger_fold_review_grok.out
4527  research/pde_ledger_v3/directives/_legs/ledger_fold_review.md
```

## the Symbol-only edge extraction the repair generalizes (Fix 1) — ledger_fold.py:146-161
```
$ sed -n '158,161p' research/pde_ledger_v3/scripts/ledger_fold.py
            names.update(_free_symbol_names(item))
        return names
    free_symbols = getattr(value, "free_symbols", ())
    return {symbol.name for symbol in free_symbols if isinstance(symbol, sp.Symbol)}
```

## the Symbol-value F9c index the repair generalizes (Fix 2) — ledger_fold.py:169-174
```
$ sed -n '169,174p' research/pde_ledger_v3/scripts/ledger_fold.py
    index: dict[str, set[str]] = defaultdict(set)
    for key, row in fold.items():
        value = row.get("value")
        if isinstance(value, sp.Symbol):
            index[value.name].add(key)
    return {name: frozenset(keys) for name, keys in index.items()}
```

## the recursion the tests do not exercise (Fix 3) — ledger_fold.py:228 and :239
```
$ sed -n '228p;239p' research/pde_ledger_v3/scripts/ledger_fold.py
                pending.append(target)
                pending.append(target)
```

## the ledger's field-as-function representation that motivates Fix 1
```
$ grep -oE "Function\('[a-zA-Z0-9_]+'\)" research/pde_ledger_v3/scripts/S11c_b_exports.py | sort | uniq -c | sort -rn | head -5
   3284 Function('O_window')
    304 Function('u2')
    304 Function('u1')
    263 Function('u3')
    188 Function('u4')
```
