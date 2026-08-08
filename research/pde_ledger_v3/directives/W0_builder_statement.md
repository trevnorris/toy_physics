# W0 — builder statement: emit the transverse-branch speed object

⚠ Do not read `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md`.

⛔ **LAUNCH HOLD.** Do not issue this statement until the work described by `W3_fix_round2.md` has landed
and the orchestrator has recorded that all four unchanged engines can regenerate their committed baseline
projections. This statement authorizes no engine work while that hold is active.

After release, this statement is shared by four independent builds: S9-WL, S9-Python, S10-WL, and
S10-Python. Each builder changes only its assigned engine and does not read another builder's change.

## Required objects

Let $r$ be a member of the selected set and $K=k\mathbin{\cdot}k$. For every such member, the final
emission contains one root-local speed bundle with these computed objects:

1. The selected-root value $r$, $K$, $v_T^2$, the two quotient-link operands $K v_T^2$ and $r$,
   and the quotient-link residual
   \[
   Q_T=K v_T^2-r.
   \]
2. A non-negative scalar $c_T$, the two square-link operands $c_T^2$ and $v_T^2$, and the square-link
   residual
   \[
   L_T=c_T^2-v_T^2.
   \]
3. Dimension witnesses for all of those objects. In particular, $c_T$ has LTM dimension `[1,-1,0]`
   and $v_T^2$ has dimension `[2,-2,0]`.
4. The premises establishing $K>0$, the required sign of $v_T^2$, and the non-negativity of $c_T$.

Emit both operands and each residual before emitting its zero guard. A label, a literal assertion, or a
guard emitted before its operands is not a computed certificate. If a required object, dimension, sign,
or positivity premise cannot be established, emit a pinned failure object at that placement; do not emit
a substitute scalar.

## Required transverse certificate and selection

The certificate covers every distinct squared-frequency root of the cell's solved spectrum. Each root
record contains:

1. the root and its algebraic multiplicity;
2. the dynamical matrix evaluated at that root;
3. the nullity of that evaluated matrix;
4. the transverse nullity \(\dim(\ker M(r)\cap k^\perp)\);
5. the computed Boolean predicate that the transverse nullity is positive; and
6. the root's computed membership decision, emitted together with the predicate value.

The selected set is exactly the truth set of the emitted per-root predicate: a root is a member if and
only if that root's emitted predicate is true. Emit the complete selected set and its computed cardinality.
For every member, emit the root-local speed bundle above. The builder has no required cardinality and must
not tune the predicate, membership, or set to a target count.

Existing computed fields may satisfy individual requirements, but the final emission must expose the
predicate-to-membership equality, the complete selected set, and the association from every selected root
to its speed bundle.

## Cells and placement grain

| Build | Required cells | Existing reusable object |
|---|---|---|
| S9-WL | one cell: MAIN with the generic three-component wavevector | squared-speed candidate |
| S9-Python | one cell: MAIN with the generic three-component wavevector | squared-speed candidate |
| S10-WL | four cells: MAIN at `D=2,3,4,5` | per-root squared-speed quotient |
| S10-Python | four cells: MAIN at `D=2,3,4,5` | per-root squared-speed quotient |

Placement is fixed as follows; tag spelling is not fixed.

- Emit one cell-level certificate and one cell-level selection summary in each table cell.
- Emit one root record per distinct algebraic root; multiplicity remains a field, not duplicated placement.
- Emit one scalar root-local speed bundle for every root whose membership is true. Do not emit a speed
  bundle for a false member and do not collapse bundles from different roots into one unassociated list.
- S9's direction-specialisation blocks remain evidence inside its one generic MAIN cell. They do not get
  separate selection summaries or separate speed bundles.

No S9 or S10 control package and no other S10 dimension is in scope. S11 and every longitudinal-speed
object are deferred.

## Required able-to-fail mutations

Preserve each mutation artifact and its literal output. A prose claim that a dependency moved is not
evidence.

1. **Upstream quotient sentinel.** In a mutation-only run, replace the $K$ input to the new quotient path
   by a fresh algebraically independent sentinel of the same dimension. Hold the solved spectrum, the
   transverse certificate, every coefficient emission, and every pre-existing emitted record fixed.
   Recompute $v_T^2$, both defining links, and $c_T$ from the mutated input. The required observable is
   that the production $c_T$ changes and depends on the sentinel while no pre-existing emitted record
   changes. Do not construct or emit an engine-side comparison target from named coefficients.
2. **Selection-only mutation.** Hold the solved spectrum, every root value, and every coefficient emission
   fixed. Change only the selection input so that a different root from that same emitted spectrum becomes
   a predicate-true member. The predicate/membership pair and selected set must remain equal by the rule
   above. The emitted selected-root value, $v_T^2$, and $c_T$ must each differ from their unmutated
   counterparts. Changing an action, matrix, root, or upstream coefficient does not satisfy this mutation.

## Satisfiable regression and manifest

After the launch hold is released and before editing, run the assigned unchanged engine to a temporary
output and establish its deterministic baseline against the committed `.out`. Never overwrite the
committed output. Compare the post-change production output with that baseline using the same deterministic
projection of physics tag/value records; exclude runtime metadata and inventory/count metadata.

- Every projected pre-existing record retains exactly the same name and value.
- The projected additions equal the manifest, with no unlisted addition, deletion, or modification.
- The manifest gives the exact emitted name **and exact emitted value** of every new production record.
- State the baseline count, final count, and numerical projection delta; the delta equals the manifest's
  cardinality.
- State separately the before/after value and numerical delta of every inventory/count record; each delta
  equals the total number of new emissions. Runtime metadata is identified but not compared.

The manifest is a review input, not an acceptance criterion. Agreement between a builder-authored output
and builder-authored manifest does not make the build pass; review must establish that every manifested
name/value is one of the required objects above and that the mutation artifacts establish its dependencies.

## Boundaries and handoff

Do not change a derivation, action, ansatz, solved spectrum, existing control package, harness
configuration, any file under `reduction/`, or any committed `.out`. Do not commit. In the production run,
do not change any existing emitted name, value, or shape except exactly the regression-named runtime and
inventory/count records: runtime may vary, and inventory/count values change only by their verified deltas.
The mutation-only artifacts may differ only as required in the mutation section.

Deliver the assigned engine change, the name-and-value manifest, the deterministic regression comparison
with counts, both mutation artifacts with literal output, and a report of at most 40 lines. This emission
work does not make any registry relation tested; binding and coverage are a separate reviewed change.
