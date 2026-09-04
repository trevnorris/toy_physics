# _measurements — S11_ledger_fold_repair2_directive.md

Grounds the guard repair in the real-ledger verification + the F9c scan; each claim carries its command (rule 2).

## authority — design §D3 refined (edge resolution by identity), committed fd8c89d0
```
$ git -C /var/projects/toy_physics log --oneline -1 fd8c89d0
fd8c89d0 Design §D3 refined: edge resolution by full symbol identity (F9c scan)
```

## the false-positives the repair must revert (verified on the real base)
## O_window is a STRUCTURAL function: 0 rows, thousands of Function() refs
```
$ grep -cE "^    'O_window':     \{" research/pde_ledger_v3/scripts/S11c_b_exports.py; grep -oE "Function\('O_window'\)" research/pde_ledger_v3/scripts/S11c_b_exports.py | wc -l
0
3284
```

## f_hold_e_W_0 is a distinct PREMISE object (a USER of W_0), not a producer
```
$ awk "/^    'f_hold_e_W_0':     \{/{p=1} p{print} p&&/^    \},/{exit}" research/pde_ledger_v3/scripts/S11c_b_exports.py | grep -E "'class'|'step'"
        'class': 'PREMISE',
        'step': 'S11c-b',
```

## the omega F9c pair resolves by identity (real=True vs not) — the one critical-path pair
```
$ awk "/^    's11b_omega':     \{/{p=1} p{print} p&&/^    \},/{exit}" research/pde_ledger_v3/scripts/S11c_b_exports.py | grep value; awk "/^    'omega':     \{/{p=1} p{print} p&&/^    \},/{exit}" research/pde_ledger_v3/scripts/S11c_b_exports.py | grep -m1 value
        'value': _restore("Symbol('omega')"),
        'value_kind': 'COMPUTED_OBJECT',
        'value': _restore("Symbol('omega', real=True)"),
```

## the mandatory new test: check_consumer on the real critical roots must resolve (repair directive)
```
$ grep -nE 'MANDATORY real-ledger test|resolves without raising' research/pde_ledger_v3/directives/S11_ledger_fold_repair2_directive.md | head
53:5. ⭐ **MANDATORY real-ledger test** (`test_ledger_fold.py`, guarded by presence of the file): import
57:   `closure_shape_deriv` — and assert **each resolves without raising** (the `omega` edge resolves by identity;
```
