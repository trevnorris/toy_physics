# Build the harness

You are the builder. Build it; do not rewrite the specification.

## What to build from

`research/pde_ledger_v3/directives/S10_harness_rebuild_directive.md` (working tree, uncommitted) is the
specification. You wrote it against a decision list and it passed a measurement audit: every headline
count in it is true of HEAD. Build what it says, with the corrections below applied.

Two review legs then read it. One found nothing blocking; the other found five items, all verified
independently by the orchestrator. Those five are corrections to the specification, not new scope.

Do not open another specification round. Apply these while building.

## Correction 1 — B5's shape vocabulary cannot express any live PY dimension source

`family` and `vector` do not cover the shape the PY engine actually emits:

```
WL configured source  {{-3, 0, 1}, {-3, 0, 1}, {-3, 0, 1}}   -> seq[3]<seq[3]<scalar>>   = family
PY per-package source (((-D, 0, 1),),)                        -> seq[1]<seq[1]<seq[3]<scalar>>>
                                                                 matches neither
```

So B5/B6's per-(engine, package) requirement is unsatisfiable on the PY side, and acceptance items 11
and 12 cannot discriminate — item 11 passes without measuring anything because every PY table is already
unassessable, and item 12 has no valid PY vector to give.

The only PY tags carrying these symbols as a clean 3-vector are
`PY_S10_LOCAL_REGISTRY_RHO_BR_DERIVED_DIMENSION_SYMBOLIC` and its `MU_R` twin, which carry no package and
no dimension in their names — keying the table on those reinstates the one-table-across-the-sweep defect
B6 exists to kill.

Add the third shape, or name the exact PY source tags and their shapes for both symbols. Whichever you
choose, both acceptance items must become performable.

## Correction 2 — B8's probe list is wrong, and acceptance item 14 contradicts B7

Measured:

```
S10 table   '2'                 checked=1 homogeneous=1 unknown=[]
S10 table   'x'                 checked=1 homogeneous=0 unknown=['x']
S10 table   'Element[x, Reals]' checked=1 homogeneous=1 unknown=[]
S9  table   'x'                 checked=1 homogeneous=1 unknown=[]
```

`x` is a vacuous pass only under S9's table, which declares `x` a primitive; S10 declares `x1..x6`, so
there `x` is an UNKNOWN_SYMBOL, which already fails operationally. B8's point stands on `2` and
`Element[x, Reals]`.

Name the config each probe runs under. Then either drop `x` from acceptance item 14 or move it to
`unassessable` — as written, item 14 requires `x` to increase `no_comparison` while B7 defines a visibly
failed lookup as `unassessable`, so a builder greening item 14 destroys the partition B7 exists to create.

## Correction 3 — the acceptance fixture cannot fail for the defect A2 names

A2 forbids a zero-locus convention because it makes `x` and `2*x` equal. Acceptance item 3 says only
"a same-shape unequal vector", and `symbolic_equal([a],[b])` is already False, so a scale-erasing
normalisation passes.

Require the `[x]` vs `[2*x]` fixture by name. Also state whether `mode: multiset` is permitted on an
action row; positional comparison is what makes the reorder half of item 3 meaningful.

## Correction 4 — the `unpaired` bucket drops rows silently

`engine_output_checks.py:856-858, 879, 1468`. Measured on S9, the only config where controls are alive:

```
CONTROL_RESPONSE: compared=170 responsive=150 invariant=20 unparsed=1 unpaired=8
main tags dropped as 'unpaired', never compared, only counted: 8
'unpaired' appears in operational_failures: False
```

A main tag whose control counterpart set is incomplete is dropped from comparison entirely. B2 covers
only the all-or-nothing case, and acceptance item 4 deletes every control tag at a cell — deleting some
leaves the row silently uncompared with only a count moving.

Give the control layer its own coverage formula and a named per-row status for a partially paired declared
cell, and add a partial-deletion ablation. Acceptance item 4 also refers to "the target engine's declared-tag
coverage", which no section defines; define it.

## Correction 5 — the symbol alias table, which is the important one

`engine_output_checks.py:576-587` holds `_SYMBOL_ALIASES`, named nowhere in the specification. Measured:

```
raw WL text:  rhoBr 1102,  muR 2440,   rho_br 0,     mu_R 0
raw PY text:  rhoBr 0,     muR 0,      rho_br 1729,  mu_R 4222
registry (quantities.yaml / relations.yaml): declares mu_R and rho_br; declares neither
             muR nor rhoBr; declares neither omega2 nor omegaSquared
```

The engines share no spelling for these constants. Every `AGREE` the harness reports — 446 today — depends
on this table, which lives inside the comparator. The module's own docstring disclaims exactly that,
saying comparison targets come from another emission, a control emission, or the registry.

**The cause is a language constraint, not a defect.** Wolfram Language reads `_` as `Blank`, so the
registry's snake_case names are unrepresentable in the `.wl` engine. A transliteration is therefore
required and legitimate. What is wrong is where it lives and what has been mixed into it.

Two of the ten entries are not transliteration:

```
omega2       -> omega**2
omegaSquared -> omega**2      =>  symbolic_equal(omega2, omegaSquared) = True
```

That asserts two differently-named quantities are the same object, and neither name is declared anywhere
in the registry.

Required:

1. **Standardise the mapping as a declared rule**, not a hardcoded list: canonical registry name in
   snake_case, engine spelling in camelCase, applied mechanically and symmetrically. Declare it in config.
   Print the mapping actually applied alongside the verdict, so a reader can see what was reconciled.
2. **Any name the rule does not cover is a declared exception** with a stated reason, and is reported.
   An undeclared spelling mismatch is a comparison failure, not something to absorb.
3. **Remove the two algebraic entries from the naming layer.** An identity between two differently-named
   quantities belongs in the registry as a declared relation, or it does not hold.
4. **Do not iterate the resulting verdicts back to agreement.** Removing those entries may turn some
   `AGREE` rows into `DISAGREE`. That is a finding and must be reported as one. Report the before/after
   count and name the rows that changed.

## Constraints

- Do not commit. Do not touch `steps/`, `paper/`, `REBUILD_HANDOFF.md`, either engine, or any committed
  engine output.
- Do not change `symbolic_equal`'s equality semantics; `DEFECT_REGISTER.md#f7` is deferred by scope
  decision. Correction 5 changes the naming layer, not the equality kernel.
- A physics disagreement must not exit non-zero; only operational failure does.
- No script may run longer than 600 seconds. A timeout means reformulate, never raise the limit.
- `reduction/test_engine_output_checks.py` invokes Mathematica when `/tmp/s9_wl.txt` is absent. That file
  exists now. This machine has a two-seat licence; do not launch a kernel.
- Report what the rebuilt instrument reads. Do not hardcode, predict, or iterate toward any counter value.

## Deliverables, then stop

1. The rebuilt `engine_output_checks.py`, both configs, the test suite, and `harness_ablation.py`.
2. The literal acceptance output the specification requires.
3. A summary under 40 lines: one line per `A1`-`A3` and `B1`-`B13`, one line per correction above, then
   anything you measured that contradicts either document.

Stop at the deliverables. Do not commit, and do not begin a further revision of the specification.
