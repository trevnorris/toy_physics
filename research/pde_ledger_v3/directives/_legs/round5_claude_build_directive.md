# Independent physics review — S11 SymPy build directive (fresh-agent leg)

**Artifact:** `S11_sympy_no_ledger_build_directive.md` (66 lines). Read order honoured: `CLAUDE.md` → spec
at `cf4a21a4` → `S10_exports.py` → `S10_brane_mode_spectrum_sympy_audit.py` → directive last. Neither
replaced engine opened. No Mathematica. Every run wrapped in `timeout 600`.

Evidence: `/tmp/s11_sympy_directive_leg_A/`.

## BLOCKING

### B1 · Four SHARED tag suffixes fixed in a PY-only document; a name mismatch DELETES the comparison row

L30–31 fixes `LAGRANGIAN`, `EULER_LAGRANGE_SYSTEM`, `M_RESIDUAL`, `M_RATIO` — four **shared** Q1/Q2 objects
the spec orders but never names. L28–29 asserts the parity obligation, a **joint property of two engines a
single-engine directive cannot enforce**. Spec `§8:1084` affirmatively permits the other builder to choose
its own name. Same for L25–26's decomposition rule (the fix for the previously measured failure at
`§8:1076-1079`), stated only here.

**Mechanism.** The `§8` join is on the prefix-stripped suffix. On a mismatch the object yields **zero**
comparison rows and two orphan tags; `§8:1100` says orphan-gap reporting *"trains its reader to ignore it."*
A genuine physics difference is then never compared and both engines exit 0.

Injecting the relative-weight error `§3:103-106` warns about:

```
SCENARIO 1 - both directives happened to choose the same name
  COMPARISON ROWS PRODUCED: 1
      residual = B_comp*(a1*k1 + a2*k2)**2/4      is zero? False
SCENARIO 2 - the WL directive chose its own name
  COMPARISON ROWS PRODUCED: 0
  ORPHAN (PY only): ['S11_MAIN_D2_LAGRANGIAN']
  ORPHAN (WL only): ['S11_MAIN_D2_LAGRANGIAN_EXPANDED']
```

Control: identical payloads → residual `0`, so the row is not tautological.
Measured: `LAGRANGIAN`, `EULER_LAGRANGE_SYSTEM`, `M_RESIDUAL`, `M_RATIO` all `appears_in_spec=False`; no
S11 Wolfram file postdates the directive, so nothing carries the names.

**Must become true.** The four names and the one-payload-per-ordered-family rule fixed at a level **both
engines read — the shared spec**. Not fixable in this directive.

### ⛔⛔ B2 · **RETRACTED BY THE LEG ITSELF** — ⛔ do not fold it

⚠ **The leg re-checked its own blocking finding under rule 13 and withdrew it.** It had claimed the sweep
cannot run under `timeout 600`. ⛔ **Not established.** The cost is a property of the **route**, ⛔ not of
the obligation — on the *same* matrix:

```
D=5 root1  DomainMatrix.rank()                          0.02s  OK
D=5 root1  DomainMatrix rref_den() pivots               0.52s  OK
D=5 root1  N3_STACKED_RANK via DomainMatrix.rank()      0.01s  OK
D=5 root1  DomainMatrix nullspace [S10:1041]            0.66s  OK
D=5 root1  Matrix.rank(simplify=False)                120.00s  TIMEOUT
D=5 root1  Matrix.rank(iszerofunc,simplify=False)[S10:1033] 143.81s TIMEOUT
D=5 root1  Matrix.rank() default                      220.00s  TIMEOUT
```

⭐ Scale, which stands: **22 declared cells**, 252–1390 tags per cell, **5,584–30,620 tags for the sweep**.

**⇒ N6 (non-blocking), what survives:** the directive imposes **elapsed time** where `§7:1041` governs by
**observable progress**, and rule 11 forbids cost as a reason to narrow. ⭐ Since the directive says the
spec wins every conflict, the cap is **void** — ⛔ a textual inconsistency to remove, ⛔ not a blocker.

⚠⚠ **The aggravating fact is real and belongs to `B3`:** three of four exact routes exceed 100–220 s on one
object at D=5 — ⛔ **including `S10:1033`, the nearest built precedent the directive points the builder at**
— while a fourth costs `0.02 s`. ⇒ a builder who writes the precedent's route, hits a cap, and reaches for
*"leave the cell incomplete"* sheds cells and exits 0. ⭐ **The fix is `B3`'s**, plus removing the cap —
⛔ never raising it.

### ⚠ A NEW RISK THE LEG LOOKED FOR AND COULD NOT CLEAR — ⭐ carry it into the BUILD review

⛔ **Neither spec nor directive pins the rank route**, and `§Q8a` takes ρ to be *"the rank your own
computation returned"*, building **all** ρ×ρ minors from it ⇒ ρ drives the minor set, the rank-drop locus,
the strata and **every component-scoped object**. ⇒ ⛔ a route-dependent ρ would be a blocking defect.

```
D=3 root1/root2: routes_completed=4/4   ALL COMPLETED ROUTES AGREE? True
D=4 root1/root2: routes_completed=4/4   ALL COMPLETED ROUTES AGREE? True
D=5: two routes time out ⇒ agreement UNMEASURED
```

⭐ **D=5 is a coverage gap, ⛔ not a clear.**

### B3 · "REPORTING IS SUCCESS" makes emission CONDITIONAL ON THE PAYLOAD — corollary 4's exact prohibition

L61–64's *"leave the cell incomplete"* means **omit the tag**, on a condition about the payload's content.
Spec `§5` corollary 4 `:193-197`: *"EMISSION MUST NEVER BE CONDITIONAL ON A PAYLOAD'S VALUE … a value absent
is indistinguishable from never computed … No tag for a named object may be emitted 'if' anything."*
Repeated at `Q5:446` and `§5:294`.

**Mechanism.** The trigger is a **physics** condition, and the spec's answer to that condition is always
**emit the underdetermined object** — that is the measurement. Sharpest instance: `§Q11` supplies no
interface condition and forbids inventing one (`:893`); `C1_EQUATIONS` is ordered *"emit the list as built,
whatever length it has"* (`:919`); `C4_DIFFERENCE` exists to *"measure the deficiency of the content
supplied"* (`:924`). **No pinned undecided/failure token exists for C1–C4.** So on the one object whose
purpose is to report that the material does not determine it, the clause reads: no pinned form → do not
emit → cell incomplete. Q11 runs in every cell ⇒ **every cell in `SKIPPED_PAIRS`, exit 0.**

**Must become true.** Distinguish *"the object IS an empty/underdetermined/undecided object — emit it, that
is the answer"* from *"the engine failed operationally."* ⛔ Omission must never respond to payload content.

### B4 · `PD_TERM` is outside the scope sentence, and "Emit no additional shared tags" closes the door

The end-state sentence covers *"every **Q1–Q11** object"*. `PD_TERM` is a **`§7`** object. The directive
rescues other `§7`/`§9` objects it needs (`RUN_PAIRS`/`SKIPPED_PAIRS`, `PREMISE_INVENTORY`) and rescues
`PD_TERM` nowhere; L31 then forbids additional shared tags.

Spec `§7:996`: *"DO NOT TYPE `P_D` FOR ANY `D`. It must be read out of the `V6` object the engine computed,
so that corrupting the census moves the package's action"* — corollary 5's wiring requirement. Emitted once
per `(package, D)` for all eight packages; for the seven carrying no `P_D` the repetition is itself the
result (`§7:937-939`). Omitting it removes the only per-cell evidence that Q9's census reached the action.

```
PD_TERM   in_spec=1  in_directive=False
P_D       in_spec=1  in_directive=False
V6_BASIS  in_spec=1  in_directive=False
```

**Must become true.** The end-state sentence spans `§7` and `§9`, or `PD_TERM` is named.

### B5 · "Symbolic tests remain CAS objects rather than Python booleans" is UNVERIFIABLE as specified

L24 states the rule; L22–23 requires only *"re-parseable SymPy expressions"*, and `§8:1086-1087` permits
**either** `sympy.srepr` **or** `str()`. Under `str()` the rule cannot be checked.

```
=== rendered with str()  -- section 8 explicitly permits this ===
  compliant : STATUS_TOKEN: PROVED_TRUE, TEST_OBJECT: True, OPERANDS: (x**2 + 1, x**2 + 1)
  forbidden : STATUS_TOKEN: PROVED_TRUE, TEST_OBJECT: True, OPERANDS: (x**2 + 1, x**2 + 1)
  BYTE-IDENTICAL? True

=== rendered with sympy.srepr ===
  compliant : ... TEST_OBJECT: true, ...
  forbidden : ... TEST_OBJECT: True, ...
  BYTE-IDENTICAL? False
```

This is `§4:162`'s test verbatim — *"if you deleted the computation, would this tag change?"* — and
`§5:226-229` records it as the defect class the rebuild exists to remove.

**Must become true.** Every boolean-valued payload (`_IDENTICALLY_SATISFIED`, `_INCONSISTENT`,
`_REAL_ADMISSIBLE`'s `TEST_OBJECT`) rendered so a CAS boolean is distinguishable from a host boolean.
⚠ **The identical trap exists in Wolfram** (`InputForm[True]` is `True` either way) ⇒ the Wolfram directive
needs the same pin.

### B6 · `RUN_PAIRS` accumulation is not checkable by the engine

L57's condition cannot be evaluated without a declared per-cell required-tag manifest, which neither the
directive nor the spec asks for. The natural implementation is *"control flow reached the end of the cell."*
Both implementations below satisfy L57–59 and neither copies the declaration:

```
IMPL_A  inject_defect=True
  RUN_PAIRS: [('MAIN',2),('MAIN',3),('XFORM_EXTRA',3)]   SKIPPED_PAIRS: []
  tags emitted / required: 24 / 27
IMPL_B  inject_defect=True
  RUN_PAIRS: []   SKIPPED_PAIRS: [('MAIN',2),('MAIN',3),('XFORM_EXTRA',3)]
  tags emitted / required: 24 / 27
```

`IMPL_A`'s run record **affirmatively certifies an incomplete cell**.

## NON-BLOCKING

**N1 · Q6r's resolution rule adds a condition the spec does not have, deleting the build's only cross-step
comparison.** Spec `Q6r:566-568` keys resolution on **the lookup**; L42–44 keys it on the lookup **and** the
dimension value **and** every provenance field on both rows. Baseline agrees, so the difference is
invisible; form control — strip one metadata field, leave the dimension value intact:

```
  spec rule       RESOLVED=['mu_R', 'rho_br']
  directive rule  RESOLVED=['rho_br']
  COMPARISONS DELETED by the directive rule: ['mu_R']
```

Also: six of eight mapped names fail at step 1, so a literal `LEDGER[name]['dimension_key']` raises
`KeyError`; three distinct failure shapes collapse to one token. **Must become true:** resolution keyed on
the lookup; missing provenance reported *as* provenance.

**N2 · Three Q6r tags not marked engine-local** — `Q6R_RESOLVED_COEFFICIENTS`,
`Q6R_UNRESOLVED_COEFFICIENTS`, `Q6R_RESIDUAL_SCOPE`. If emitted shared: **3 × 22 = 66 orphan rows**.

**N3 · `M_ROUTE_USED` is corollary 5's named wiring obligation and the directive endorses the wrong
property.** Spec `Q2:344-346` requires its payload *"read from the same route-selection object that supplies
the matrix actually consumed downstream, ⛔ not from a second literal beside it"*; the corollary-5 exemption
list is closed at four entries and `M_ROUTE_USED` is not one. L24 *"Tokens remain symbolic tokens"* is about
**rendering, not sourcing** — `Str('M_B')` typed beside the code satisfies it.

**N4 · The directive orders a report against a Wolfram stream that does not exist.** The tempting response
is to open the old WL engine, which L8 rules out.

**N5 · The decomposition rule contradicts its own Q6r indexing** — L25–26 says only spec-indexed scopes
become separate tags; L47–50 introduces `Q6R_RESOLVED<q>_…`, a coefficient index the spec does not define.

## What ran and found nothing

- **Rule-5 leak scan**: 10 of 66 lines flagged, **all benign** (section numbers, `timeout 600`, "nonzero
  exit"). States no measured value, count, membership, sign or spectrum; names none of the eight Q6r
  coefficients and states no resolution outcome. **Clean on rule 5.**
- **Completeness of the four fixed names**: walking every "Emit" instruction, the shared objects the spec
  orders without naming are exactly those four. The directive's list is **complete**.
- **Q9 cost probe**: could have found the census infeasible; **it did not** — exact `V1`/`V2` nullspace plus
  the `V6` reflection eigenspace across D=2..5 totals **1.78 s**. The runtime problem is Q4/Q8a ranks.
- **`_CANONICAL_LOCUS` lex Gröbner probe**: ≤ 0.22 s. Not the bottleneck.
- **Dangling-`dimension_key` probe**: 0 of 50 rows point at a missing key.
