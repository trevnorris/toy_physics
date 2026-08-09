# S9/S10 export chain — fix round 3. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: `scripts/S9_light_requires_shear_sympy_audit.py` · `scripts/S10_brane_mode_spectrum_sympy_audit.py` ·
`scripts/extract_knob_inventory.py` · both generated modules · both regenerated `scripts/out/*.out` ·
a new directive under `directives/`. ⛔ Out of scope: every `.wl`, every step record, every `.tex`.

⭐ Two independent legs measured round 2. **All five repairs hold under mutation** — each rejects its
weakest change, and **no computed physics value moved.** ⛔ Every item below is a defect the legs found in
the material round 2 added. ⚠ S10 is slow: wrap runs in `timeout 1800`.

---

## F1 · ⛔⛔ One quantity, one name, one object — across the whole chain.

⚠ **The most serious finding, and round 2 created the appearance of it being solved.** S9 declares the
squared frequency `omega2`; S10 declares it `omegaSquared`. **Both are now ledger keys, and no record
carries both.** Measured: binding `LEDGER['omegaSquared']` and substituting moves **53** records and
**leaves `factored_determinant_d3` — the light cone — untouched.** Same split for `lambda_scale` vs
`lambdaScale` and for the dimension-unknown alphabet (`dim_mu_R_L…` vs `mu_R_dim_length…`).

⛔ `symbol_binding_residual` reads 0 **by construction**: it polices name→object, ⛔ never quantity→name.

**Decision: a quantity carries one name across the chain, and a consumer binding it reaches every record
that uses it.**

⚠ **Two things you must NOT decide silently:**
- The **assumptions differ** (`lambda_scale` is `positive`, `lambdaScale` is `real`). ⛔ That is a physics
  question, not a spelling one ⇒ **name it and stop** if the answer is not forced.
- ⭐ Where a third engine already spells the object, that spelling is evidence about the standard name.
  ⛔ Do not assume the newer engine is right.

## F2 · The corroboration must survive the merge.

The three overwritten records are the **only cross-step corroboration the ledger has** — S10 re-derives
them from its own general-`D` action and lands on S9's values. Measured: the merge stamps `step: 'S10'`
and **the agreement survives only in the transcript row**, not in the ledger.

**Decision: an entry that two steps independently computed and agreed on says so in the ledger.**

⭐ Record, ⛔ do not try to fix: `overwrites_upstream` is an author-set flag. A leg built a record that
declares the overwrite and copies S9's value verbatim — **exit 0, all residuals 0.** ⛔ Nothing can check
that a value was derived rather than copied; **state that limit.**

## F3 · The traversability guard must reach inside.

`traversable_export_value` returns any `sp.Basic` **unchanged without inspecting it**. Measured: a value
whose `.args[0]` is a raw python `list` is accepted, published and round-trips. ⛔ Not realized in the
artifact — a hole in the check, not a live defect.

**Decision: the guard inspects the whole object, not its top type.**

## F4 · An authored word is not a bindable quantity.

Round 2 made the distinction structural (`Str`) at **two** sites; **~30 sites still emit authored words as
`sp.Symbol`**, and **11 of them are now ledger records** indistinguishable from a physics quantity
(`roots_returned`, `solution_returned`, `M_B`, `exactly_determined`, …). ⇒ **24 records in 6 families hold
a word as their whole value**, so a downstream difference across `D` is **zero by construction** — a word
cannot vary with `D`.

**Decision: apply the round-2 distinction everywhere, not at two sites.** ⚠ Three of those sites are
**unconditional**, not branch-selected — a typed word where a readout was intended.

## F5 · A counter that cannot be anything but zero must go.

The `dim` field was deleted (50 → 0) but the counters reporting on it survive; the transcript now reads
`COMPUTED_DIMENSION_COUNT: 3 → 0` and `50 → 0`. ⇒ a reader sees dimension coverage collapse when in fact
the field was removed.

**Decision: delete the dead counters, and make an exported object's dimension DISCOVERABLE.** ⚠ The join
existed for **53** records and is gone — a consumer can no longer ask what units an object carries.
⛔ A sibling key the consumer has to guess the spelling of is not discoverable.

## F6 · Validate, then write — on both engines.

`S9:981` writes the module; the asserts are at `S9:1023-1026`. **S10 does it in the correct order.**
Measured on both legs: a failed identity scan leaves a rejected `S9_exports.py` on disk carrying both
`Symbol('mu_R', positive=True)` and `Symbol('mu_R')`.

**Decision: neither engine leaves a generated module behind when its own guards fail.**

## F7 · An exported value must not be mutable.

15 (S9) / 35 (S10) records hold `MutableDenseMatrix`. Measured: `LEDGER['length_dimension']['value'][0,0] =
Symbol('ZZ')` succeeds — **the ledger's record changes with no run.**

**Decision: what the ledger hands a consumer cannot be edited in place.**

---

## Constraints

- ⛔ The derivation, action, ansatz and every computed physics value stay untouched. ⭐ Report the complete
  `.out` diff for both engines and say whether anything beyond names, containers and export bookkeeping moved.
- ⭐ Every S9 record S10 does not overwrite stays identical between the two generated modules.
- ⛔ Do not add any in-run check on the export writer beyond repairing the ones named here.
- ⛔ Do not state anywhere what a count, tally or partition size comes out as. Emit it and let it be read.

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** For each item, construct the weakest change that
should be rejected and show whether it is.
⛔⛔ **A FORM ablation is mandatory** — change the structure of a load-bearing object, not a coefficient.

## Deliverables

The changed files, both literal stdouts, all diffs, and every ablation script with its stdout at named
absolute paths **outside the repository**.
