# W2 fix round 1 — seven decisions from two independent review legs

⚠ **Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

You built the registry-dimension witness. Two legs reviewed it and exhibited defeating cases with literal
output. ⛔ Do not re-litigate any decision below — each is backed by an executed counterexample.

**Files you may touch:** `reduction/registry_dimension_witness.py`, `.yaml`,
`registry_dimension_witness_able_to_fail.py`, `test_registry_dimension_witness.py`, `reduction/README.md`.
⛔ No engine, no step record, no committed `.out`. ⛔ Do not commit.

**Report:** ≤40 lines — one line per decision, plus the new committed-output run and the able-to-fail runs
with **literal stdout**.

---

## D1 · ⛔ The `D_brane` row is registry-vs-registry — stop counting it as agreement

`PY_S10_LOCAL_REGISTRY_D_BRANE_DIMENSION` is the engine reading `quantities.yaml` and emitting it back
(`S10_..._sympy_audit.py:1683,1688`). A leg proved it cannot fail: with the registry mutated to `[7,0,0]`
the tag emits `(7,0,0)` and the residual is zero. ⚠ Your own docstring says the script *"does not consume
an engine's own … declared-dimension output"* — the manifest contradicts it.

⭐ **Do not simply delete the row.** Classify sources: a source is **DERIVED** (the engine computed the
vector) or **ECHOED** (the engine re-emitted a registry declaration). ⭐ An ECHOED source gets its own
status — ⛔ never `AGREEMENT` — and is excluded from any agreement count. ⭐ Fix the docstring to match
whatever the code then does.

## D2 · ⛔⛔ Specialise with the BRANCH's dimension, not the registry's

`_specialise_emission` resolves `braneDimension`/`D` from `Q.brane.D_brane` unconditionally, discarding the
branch the tag came from. A leg pointed the manifest at S10's **D=5** cell and got:

```
WITNESS|artifact=S10-wl|quantity=Q.brane.rho_br|status=AGREEMENT|
  emitted_raw={[-braneDimension,0,1]@WL_S10_MAIN_D5_Q6_INERTIAL_COEFFICIENT_DIMENSIONS}|
  emitted={[-3,0,1]}|difference={[0,0,0]}          ← engine's actual D=5 value is [-5,0,1]
```

⛔ **S10 exists to run D = 2,3,4,5. The check is blind along the axis the step varies**, and the only thing
hiding it is the hand-written `dimensions: [3]` selector — load-bearing physics disguised as a filter.

⭐ Specialise with the branch dimension the tag was selected for (it is already the loop variable at
`:215-218`), **or** read the engine's own `*_SPECIALIZED` companion tags. ⛔ When the branch dimension
differs from the registry's declared value, that is an **explicit status**, ⛔ never a silent substitution.
⚠ Also: `specialisations` is keyed by bare symbol name over a whole artifact, so any free symbol spelled
`D` is replaced regardless of meaning. ⭐ Bind the substitution to the source, not the artifact.

## D3 · ⛔ Replace the pinned integers with a coverage invariant

`test_registry_dimension_witness.py:101-103` pins `len(rows) == 9`, `all(status == "AGREEMENT")`, and
`47`. ⛔⛔ **That is a sampled outcome asserted back** — the failure class this project has already paid
four rounds for. It also breaks for reasons unrelated to instrument health the moment a quantity or
artifact is added.

⭐ Replace with invariants: **every `(artifact, quantity)` pair named by a selected `dimension_sources`
entry produces a row** · **an artifact yielding zero rows is a guarded status** · **a `tag_template` that
expands to zero tags is a guarded status**, ⛔ not the silent `if not tags: continue` at `:279-280`.

## D4 · ⭐⭐ Compare the SPEEDS — they are emitted, as squares

⛔ `Q.brane.c_gamma` and `Q.brane.c_L` currently report `NOT_EMITTED`, and that reads as *"the engines said
nothing."* They did:

```
c_gamma declared [1,-1,0] ; emitted  WL_S9_MU_RHO_DIMENSION_DIFFERENCE          {2,-2,0}
                                     PY_S9_MAIN_DIM_STIFFNESS_MINUS_INERTIA     [2,-2,0]
                                     WL_S9_VELOCITY_SQUARED_DIMENSION
c_L     declared [1,-1,0] ; emitted  WL_S11_MAIN_D3_ROOT1_DIM_OVER_KSQ          {2,-2,0}
```

⇒ `[2,-2,0] = 2 × [1,-1,0]`. ⭐ **Add a DECLARED multiplier (power) on a source**, so a source may state
that its tag carries the dimension of *quantity^n*. ⛔ The multiplier is declared in the manifest and
**printed on every row beside the operands and the residual** — ⛔ never inferred from the numbers.

⚠ **This is where a wrong declaration could force agreement**, so the row must print: declared vector ·
multiplier · emitted vector · residual. ⭐ Add an able-to-fail case that perturbs **the multiplier** and
confirms the row flips.

⭐ Add `Q.brane.D_brane` at S9-wl, where `WL_S9_SPEED_SQUARED_SYMBOL_DIMENSION_MAP` carries `D -> {0,0,0}`
— ⚠ but first determine whether that entry is **derived or echoed** and classify it per D1. ⛔ Do not
assume.

## D5 · ⛔ Stop describing the convention as verified

`README.md:170` says *"The emitted convention is checked before an exponent residual is formed."* ⛔ It is
**declared** in the manifest and compared to another declaration; the witness never reads an axis label,
and `ordered_bases` is never read at all. A leg permuted an engine's axis order to `M,L,T` and the witness
reported `DISAGREEMENT`, ⛔ never `CONVENTION_MISMATCH`. ⚠ And a permutation-invariant declaration
(`[0,0,0]`) is blind to axis order entirely.

⭐ **Prefer the measurement**: S9's SymPy engine already emits axis-labelled components
(`dim_rho_br_L`/`_T`/`_M` in `DIM_SOLUTION`) — read them where they exist and report `UNDETERMINED` where
they do not. ⭐ If you keep the declared form anywhere, restate `README.md` and the report field so
`emitted_convention=` reads as an **input**, ⛔ not a verification.

## D6 · ⛔ Calibrate every KIND of source, not one row

All three cases target `S9-py` / `Q.brane.rho_br` / one tag (`:21-23`), so **1 of 9 rows** is exercised —
and ⛔ not the hand-added one. ⭐ Cover each source **kind** (inherited-config, multiplier, echoed), and
⭐ add a case that perturbs **the branch dimension** rather than an exponent.

## D7 · ⭐ State the scope on the report line

`NOT_EMITTED=47` is coherent slot accounting (14 quantities × 4 artifacts − 9) ⛔ but it reads as a
project-wide claim. ⭐ Either add the S11 artifacts, or print the artifact set the count is scoped to.

---

## Standing rules

- ⭐ **Print computed objects; ⛔ never state conclusions.** Both operands, the multiplier, and the
  residual — then guard.
- ⛔ No `assert` before a value it guards.
- ⛔ **Do not tune anything to make the committed run green.** ⭐ If the corrected witness reports a
  disagreement on committed outputs, **that is the deliverable** — report it and stop.
- ⭐ Re-run the full `reduction/` suite and report the literal result.
