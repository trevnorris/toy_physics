# S10 harness repair, part A — make the ABLATION layer measure the sweep

**Primary file:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py`
**Config:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S10.yaml`
**Tests:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/test_engine_output_checks.py`
**One engine edit, item A5 only:** `mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`

⛔ **Do not commit.** ⛔ **Do not edit** `scripts/S10_brane_mode_spectrum_sympy_audit.py`, ⛔ anything
under `steps/` or `paper/`, ⛔ or `checks_S9.yaml`.
⚠ **Read `directives/S10_SHARED_PHYSICS.md` §7 and §8 before starting** — §8 defines the tag grammar this
checker parses, and §7 defines which packages run at which `D`.

⚠ **A second part (the parser and dimension-table repairs) follows separately. ⛔ Do not attempt it here**
— in particular ⛔ leave `derived_dimensions`, `_as_dimension`, and the Mathematica expression reader
alone except where an item below names them.

---

## ⭐ WHAT THIS IS

The checker was written for **S9**: one fixed spatial dimension, control packages named `X1…X8`.
**S10 sweeps `D = 2,3,4,5` and names its controls in words** (`XFORM_FULLGRAD`, `XFORM_DIVONLY`,
`XFORM_SIGNFLIP`, `XFORM_ANISO`, `XCOEF_SCALE`), which run at only **some** of those `D`.

⛔⛔ **The consequence is that the layer which tells a computed output from a typed one has been
comparing nothing, and saying so in a line that reads like health.**

⚠ **An earlier version of this directive proposed a fix, and a review leg demonstrated that the fix
would have produced a layer that is 79% dead while passing the acceptance test the directive asked for.**
⭐ Everything below is the corrected version. ⛔ Where this document contradicts anything you were told
before, this document wins.

---

## ⛔⛔ A1 — the layout must be keyed on `(package, D)`, not on package

```python
CONTROL_TAG_PATTERN = re.compile(r"^(?P<base>.+)_X(?P<control>[1-9][0-9]*)_(?P<suffix>.+)$")
```

⚠ **Measured: zero S10 tags match this.** `_package_layout` returns empty; `CONTROL_RESPONSE` reports
`compared=0 responsive=0 invariant=0` and `TAG_PARITY` reports `packages=0 gaps=0`.

**A1a — declare the packages in the config.** Add to `checks_S10.yaml`:

```yaml
main_package: MAIN
control_packages: [XFORM_FULLGRAD, XFORM_DIVONLY, XFORM_SIGNFLIP, XFORM_ANISO, XCOEF_SCALE]
```

⭐ Keep the existing `_X<digits>_` behaviour as the fallback when the config declares nothing, so
`checks_S9.yaml` keeps working unchanged. ⚠ **Match by LONGEST declared prefix** — package names contain
underscores, so a shortest-match rule splits them in the wrong place.

**A1b — ⛔⛔ and this is the part the first version got wrong. Key the layout on `(package, D)`.**
`check_control_response` requires a suffix to be present in **every** control package
(`len(control_tags) != len(control_packages)` ⇒ `unpaired`). But §7 runs `MAIN` at `D = 2,3,4,5` and the
controls at only some of those. ⚠ **Measured under the package-only fix:**

```
row statuses     : UNPAIRED 661 | RESPONSIVE 105 | INVARIANT 59 | UNPARSED 11
per-D [compared, unpaired] : D2 [0,207]   D3 [175,40]   D4 [0,207]   D5 [0,207]
```

⇒ ⛔ **the ablation would cover `D3` and nothing else**, while `compared = 164 ≠ 0` so no guard fires.
⭐ **Fix: a main tag at `D` is compared against the controls that ran at that same `D`.** A `D` with no
controls is **reported as uncovered**, ⛔ not silently folded into `unpaired` alongside genuine gaps.

**A1c — ⛔ exclude the engine-local family from the layout.** Tags are spelled
`WL_S10_LOCAL_MAIN_D2_…` — `_LOCAL_` sits **before** the package name. ⇒ package matching invents a
second base `WL_S10_LOCAL` with its own main and five controls, pulling the engine-local tags that §8
exists to keep **out** of comparison **into** the ablation layer. ⭐ Exclude them from `_package_layout`
itself. ⚠ 132 such tags in engine 1, 111 in engine 2.

---

## ⛔⛔ A2 — the guards. ⭐ Coverage, ⛔ not `compared > 0`

**A2a — a layer ASKED to compare something that compared nothing is an OPERATIONAL FAILURE.**
If a declared package contributes **zero** tags, fail and name it.

⛔⛔ **Scope it precisely.** The rule is **declared-but-matched-nothing**, ⛔ never "any count of zero".
⚠ `REGISTRY_RESIDUAL: configured=0` and `CROSS_ENGINE: configured=0` mean *nothing was asked for* — a
legitimate state that both shipped configs are in today. ⛔ Failing on those breaks `checks_S9.yaml`,
which you may not edit.

**A2b — ⭐⭐ guard on COVERAGE, and this is the guard that would have caught the original defect.**
`compared > 0` is satisfied by a layer that is inert over most of its input. ⭐ Report, and fail on,
the **uncovered fraction**: how many main tags were compared against a full control set, how many were
not, and **which `D` values got no ablation at all**.

**A2c — ⛔ `responsive` and `unpaired` have no failure condition anywhere.** `operational_failures`
ignores the entire `ControlReport`. ⇒ a run reporting `compared=677 responsive=0` — the ablation's
**total** failure — exits 0. ⭐ Add a floor: if a layer compared a substantial number of tags and **none**
responded, that is a failure, not a finding.

**A2d — ⛔ a base with NO control packages currently scores every main tag `INVARIANT`.** With no
controls, `len(control_tags) == len(control_packages)` is trivially true, so the row is recorded as
"the controls did not move it" **when no control was consulted** — and it inflates `compared`, defeating
A2a directly. ⭐ Such rows are **uncovered**, ⛔ not invariant.

---

## ⛔ A3 — `TAG_PARITY` does not measure what its name says

⚠ **It is INTRA-engine**: each control package against **that same engine's** main package. It never
compares one engine to the other. ⭐ Demonstrated by the review leg: an engine carrying a tag the other
lacks entirely still yields `GAPS: 0`. The layer's name, the header comment in `checks_S10.yaml`, and
§8 all describe **cross-engine** parity.

⭐ **Fix the description, ⛔ not (yet) the layer.** Rename the report line and its config comment to say
plainly that it is per-engine main-vs-control parity. ⭐ **Then add a genuine cross-engine tag-set
comparison** as a separate reported line: tags present in one engine and absent in the other, both
directions, counted and named.
⚠ Expect that count to be large — it is a known naming divergence, ⛔ not necessarily a defect — so
report it, ⛔ do not fail on it.

**A3b — ⛔ `PARITY_EXCLUDED` reports an exclusion that is not in force.** A pattern matching zero tags
prints `excluded=0 by_pattern={'_LOCAL_': 0}`, which reads as an exclusion working. ⚠ **Measured today:
both engines emit `_LOCAL_` tags and the line reports `excluded=0`**, because the tag set it draws from
comes from the empty layout. ⭐ Once A1 lands this changes on its own; ⭐ **add a guard anyway** — a
configured pattern that excludes nothing must say so as a warning.

---

## ⛔⛔ A4 — the `homogeneous` count is mostly not a measurement

⚠ **Measured by the review leg: 1699 of 2139 "homogeneous" tags had ZERO comparisons performed.**
`checked` increments for any tag that does not raise, and `homogeneous` increments when no complaint was
recorded — so a payload of `{}`, of `2`, or of a bare symbol scores **homogeneous**. By kind:
`list=1289, integer=403, relation=7`.

⛔⛔ **S9's `checked=1219 homogeneous=1219 non_homogeneous=0` is the terminal form of this**, and that
number is cited in a step record.

⭐ **Fix: count and report separately (i) tags where at least one dimensional comparison was actually
performed and (ii) tags that passed vacuously.** ⛔ Do not merge them. ⛔ Do not change the dimension
arithmetic — that is part B. ⚠ This item is **reporting only**, and it is the highest-value line in the
whole report because it is the one a reader currently over-trusts.

---

## ⛔ A5 — Mathematica diagnostics reach stdout (⭐ the only engine edit here)

`Solve::svars` prints **10 message lines plus 10 blank lines** into engine 1's output; the harness parses
them as tags.

⭐ **Fix: `Quiet` the message at the solve site — ⛔ only that message, ⛔ only at that site.** The engine
already emits `..._SOLVE_SVARS_MESSAGE` tags; ⛔ do not remove, weaken, or stop emitting those.
⛔ No blanket `Quiet`/`Off` around anything larger than the individual `Solve`. ⛔ Do not change what is
solved or how.

⚠ **Then re-run engine 1** and overwrite `mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`.
⭐ Launch it in the **background**, ⛔ never foreground under `timeout`. Budget **600 s**; ⛔ if it exceeds
that, stop and report. ⚠ An orphaned `WolframKernel` leaks memory without limit — confirm with
`ps -eo pid,rss,etime,comm | grep WolframKernel` that none survives.

⛔⛔ **A5b — the re-run rewrites all 2983 payloads and the obvious acceptance cannot see a changed
value.** ⭐ **Diff before against after**: the sorted tag-name set must be identical, and every payload
must be unchanged except for the removal of the stray lines. ⚠ `Quiet` must not change which branch
`Solve` returns. ⭐ Report the diff summary.

⚠ **A5c — the class, not the instance.** `TAG_PATTERN` accepts `Solve::svars: …` as a tag named `Solve`.
Only the **duplicate-key** check caught these ten; a *single* stray CAS message would have entered the
tag stream as data and been dimensioned and compared. ⭐ Reject a "tag" whose name is not in §8's
grammar, and report any line rejected.

---

## ⛔⛔ ACCEPTANCE — ⭐ COVERAGE ablation. ⛔ Corrupting one payload is not enough.

⚠ **The review leg showed that corrupt-one-payload passes a layer that is inert over 79% of its input,
because the original defect was exactly unreachability.** Under the package-only fix:

```
corrupt MAIN_D3_PREMISE_BACKGROUND_STATE : INVARIANT -> RESPONSIVE   <- acceptance PASSES
corrupt MAIN_D2_Q3_ROOT_COUNT            : UNPAIRED  -> UNPAIRED     <- layer is dead
corrupt MAIN_D5_Q3_DETERMINANT           : UNPAIRED  -> UNPAIRED     <- layer is dead
```

⭐⭐ **Required instead: a coverage ablation.** On a **copy** of an engine output, corrupt payloads
**one at a time across a large sample spanning every package and every `D`**, re-run, and report:
- the **fraction** of corrupted tags on which the layer moved, and
- **the names of the tags that did not move.**

⛔ A layer that moves on one tag and not on a whole `D` is not fixed. ⛔ Corrupt the copy, ⛔ never the
committed output.

⭐ **Four further requirements, each of which a corruption test alone does not establish:**

1. ⭐ **Ablate every control package, one at a time** — a layer that only ever consults the first control
   passes a single-corruption demo.
2. ⭐ **A SYMBOLIC-EQUIVALENCE probe.** The same object in two equivalent renderings — `x + y` against
   `y + x`, `1/2` against `Rational[1,2]` — must report **AGREE**. ⛔ A raw-string comparison passes every
   corruption test and fails this.
3. ⭐ **An EMPTY-INPUT probe.** Feed a layer an output with none of its tags; it must complain.
4. ⭐ **Corrupt into a WRONG VALUE, ⛔ not into something unparsable** — that only moves the UNPARSED
   count without exercising the comparison.

### ⭐ The commands

```
cd /var/projects/toy_physics/research/pde_ledger_v3 && \
timeout 900 python3 reduction/engine_output_checks.py \
  --config reduction/checks_S10.yaml \
  --output wl=mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out \
  --output py=scripts/out/S10_brane_mode_spectrum_sympy_audit.out
```

⚠ **The run will still exit 2 on 220 UNPARSED payloads — that is part B's work, ⛔ not a failure of
yours.** ⭐ Judge every acceptance by the **reported counters**, ⛔ not by the exit code.

### ⭐ The unit tests, and their baseline — ⛔ this is measured, do not re-derive it

`python3 -m pytest reduction/test_engine_output_checks.py -q` on the untouched tree gives
**`2 failed, 15 passed`** (105 s). ⛔ **Both failures are pre-existing and are NOT yours.** The cause:

```python
REAL_S9_OUTPUTS = {"wl": Path("/tmp/s9_wl.txt"), "py": Path("/tmp/s9_py.txt")}
```

⚠ Those `/tmp` fixtures are **1487 and 590 lines**; the engines' committed output is **1559 and 635**.
The pair `WL_S9_X7_ROOT1_SCALING_RESIDUAL` / `PY_S9_X7_ROOT_SCALING_RESIDUAL` **exists in the committed
output and is absent from the fixtures**, which is the reported `MISSING`.

⭐ **Fix: point the tests at the committed outputs**, `mathematica/out/S9_light_requires_shear_*.out` and
`scripts/out/S9_light_requires_shear_*.out`. ⚠ Both tests should then pass; ⭐ report what they do.
⛔ Do not change any other test to match your work — if one fails, ⭐ tell me it failed and why.
⚠ `test_absent_parity_exclude_preserves_report_bytes` freezes the report text byte-for-byte and your new
report lines will break it. ⭐ Report it, ⛔ do not quietly re-baseline it.

⭐ **Finally run the S9 config as a regression check** and report its summary line:

```
timeout 900 python3 reduction/engine_output_checks.py --config reduction/checks_S9.yaml \
  --output wl=mathematica/out/S9_light_requires_shear_mathematica_audit.out \
  --output py=scripts/out/S9_light_requires_shear_sympy_audit.out
```

---

## Report back — ⛔ under 35 lines, plus the coverage table

1. One line per `A1a`–`A1c`, `A2a`–`A2d`, `A3`, `A3b`, `A4`, `A5`–`A5c`: fixed / partially / not, with
   file and line numbers.
2. **The coverage-ablation table**: layer, sample size, moved fraction, and the tags that did not move.
   ⛔ This is the deliverable, not an appendix.
3. The four extra probes: pass / fail, one line each.
4. Exit code and counters of the acceptance command; the pytest result against the stated baseline; the
   S9 regression line.
5. ⛔ **Do not report what any physics value came out to be.** Counts are fine; ⛔ values are not yours
   to summarise.
6. ⭐ Anything here you believe is **wrong**, and anything you changed that I did not ask for. ⭐ Wanted.
