# S15 reading comparison inputs for items 1–3, 5 and 6

The lookup indexes the accepted reading run at `cc62dc9f` and the original S15 run at `3192123b`.
It verifies the full SymPy reading file, the Wolfram reading file through its base-pass marker, and
both original files against the accepted hashes. The Wolfram manifest records the base pass as a
separately hashed prefix. The runner appends the item-4 passes after that prefix. This lookup
makes no statement about the current bytes beyond it. [Wolfram manifest](S15_separation_departure_wolfram_run_manifest.json)
· [Wolfram output writer](../../mathematica/S15_separation_departure.wl) · [Original run record](S15_build_run_record.md)

The 336 names in each inventory include the surface, group, and Cartesian component where applicable.
The 309 identical payloads per engine link each reading expression to its same-engine original.
The 27 sums combine rest, transport or crossing terms on the inner, outer or both surfaces, in
each Cartesian direction. The lookup locates each original operand. Checking those sums and comparing
mathematical values across engines belong to the comparator. For item 3 the status counts below
are the engines’ own labels. For item 2 a difference in printing literal zero can arise from an
unevaluated integral and requires interpretation of the full expressions. [SymPy sum definitions](../../scripts/S15_separation_departure.py)
· [Wolfram sum definitions](../../mathematica/S15_separation_departure.wl)
· [Comparison remit](../S15_reading_decision_list.md#for-the-comparator-and-the-record-not-in-the-build-directive)

Run the command below from the repository root. To obtain the object-level index, add
`--index <new-local-file>`. Each selected record carries its engine-native tag (object name and field),
byte offset and length. The index also supplies the original operand tags and offsets for each sum.
It covers the principal records of the retained items. Additional native records, including domains,
dimension operands, engine-specific branch details and diagnostics, still need accounting in the formal comparison.

```text
$ python3 research/pde_ledger_v3/scripts/S15_reading_input_inventory.py
sympy original bytes and sha256: 24641539 4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c
sympy END_OF_TRANSCRIPT bytes and sha256: 192065031 10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6
manifest objects with selected tags: 336
item 1 same-engine lookup: {'identical payload': 309, 'new sum with all original operand tags located': 27}
item 3 status tag text: {'Exact family object retained unevaluated: integrand contains a general field, as permitted by item 3.': 192, 'Closed-form return from the exact moment, distance and finite-interval integration construction. The branch inventory is printed separately.': 144}
wolfram original bytes and sha256: 9797390 811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6
wolfram END_OF_BASE_PASS bytes and sha256: 333550018 75622128d04da2d934c2a52468505b4438fff9f17a902d7a6cd35993d769b6c4
manifest objects with selected tags: 336
item 1 same-engine lookup: {'identical payload': 309, 'new sum with all original operand tags located': 27}
item 3 status tag text: {'"Exact family integral with general fields retained, as permitted by item 3."': 192, '"Closed-form return from the angular, distance and finite-interval integration construction. The exponent-branch and retained-function inventories are printed separately."': 144}
common object names after label spelling normalization: 336
SymPy-only names: []
Wolfram-only names: []
E_geo literal-zero spelling differs across engines: 3
E_resp literal-zero spelling differs across engines: 72
sympy changed native payloads with same tag: []
sympy sum rows: 27 original operand references: 144
wolfram changed native payloads with same tag: []
wolfram sum rows: 27 original operand references: 144
```
