# O2 comparator decision list: review dispositions (orchestrator)

**Artifact:** `directives/O2_comparator_build_directive.md` as reviewed. Its sha256 is in
`O2_comparator_directive_review_disposition_lookups.md`, and a frozen copy is
`_scratch/s9b_build/O2_comparator_build_directive_reviewed_v0.md`. It is orchestrator-written, so the legs are
Codex + Grok (G2): one pass, findings verified, folded once.

**Legs (identical prompt `_scratch/s9b_build/o2_comparator_directive_review_prompt.md`).** Both reported before
adjudication.
- **Codex:** NEEDS REVISION, 2 findings (`_scratch/s9b_build/o2_comparator_directive_review_codex.txt`).
- **Grok:** NEEDS REVISION, 3 findings (`_scratch/s9b_build/o2_comparator_directive_review_grok.txt`).
- **Agreed sound by both:** the oracle exclusion. Grok checked the transcripts for the oracle domains' terms and
  found the scoping §5 limits absent from both streams. Both also agree the list states no physics target, sign or
  frozen profile.

Each verification is a mechanical lookup. Commands and literal output are in
`O2_comparator_directive_review_disposition_lookups.md`.

| # | Finding (legs) | Disposition | What must be true after the fold |
|---|---|---|---|
| U1 | Item 1 lets an object be "unjoined with its reason" even when its counterpart exists under a different tag, key or head. For example, the mass residual is `PY_O2_MASS_RESIDUAL` in SymPy and the `"Residual"` key of `WL_O2_MASS_INPUT` in Wolfram. The format note names only `OPEN`/`OpenAction`/`Inactive`, but the Wolfram stream also carries `OpenFirstVariation`, `OpenNativeField`, `UnrestrictedSection` and `UnrestrictedNativeSection`. (Codex 1, Grok 1) | **ACCEPT.** Lookups: PY lines 267 and 380; WL line 247; Wolfram head counts. | An object is unjoined only when the other stream has no component, entry or key path for the same spec §9 object. A different tag name, key or head is a join row. Content that cannot be subtracted stays on its joined row as a coverage finding. No list of heads stands in for the engine sources. |
| U2 | The machinery named for reuse, `S11b_cross_engine_comparator.residual`, joins by name, returns zero for equal text atoms, compares other atoms by stock `srepr`, has a 5-second per-leaf alarm, and prints `AGREE`. With a name table, equal mapped strings would become a zero by name. The SymPy payloads are the engine's lossless constructor form, not stock `srepr`. (Grok 2) | **ACCEPT.** Lookups: S11b lines 1–8, 41, 628–630, 670–672, 816–820; PY `LosslessReprPrinter` at line 60; 4,440 lossless `o2_X_face_0` constructors in the transcript. | Reuse is limited to machinery consistent with items 1–8. Equal mapped names are a substitution, never by themselves a zero: agreement comes only from the compared content beneath them. No per-leaf time budget. SymPy payloads are read in the engine's lossless form. |
| U3 | The item 6 controls need not fail. Item 5 forbids test expectations, synthetic algebraic fixtures can bypass the serialized reader, extraction and join path, and several defect classes have no control: a stripped applied-function argument, a profile replaced by a constant, an unjoined existing counterpart, a lost OPEN dependency, orientation or sibling, and derivative or binder structure. "Under-measured" is undefined. (Codex 2, Grok 3) | **ACCEPT.** Lookup: directive lines 59–79. | Each control is a test that fails when its defect goes undetected. The assertion is value-free: the mutation changes the printed output. Fixtures are serialized in each engine's own format and pass through the same reader, extraction, join and comparison path as the production run. Each listed defect class has such a control. Coverage is an independently enumerated leaf count per object from each stream's parse, printed beside the compared count. |

**Fold.** U1–U3 were folded once into the directive, giving the folded version (sha `97260120…`, in the lookups'
"Fold applied" section). Where each lands:
- **U1:** the Inputs (no head list; the engine sources are the authority) and item 1 (the unjoin rule).
- **U2:** the Inputs (lossless SymPy form; S11b `residual` not reusable as is), item 2 (equal mapped names are a
  substitution) and item 3 (no per-leaf time budget).
- **U3:** item 5 (the no-target rule covers the measured streams, so value-free test assertions are allowed) and
  item 6 (serialized end-to-end fixtures, a control for each listed defect class, and independent leaf counts).

The routed r3 items are carried: item 1 injectivity and content-based pairing, item 4 for the representational
differences, and item 7 keeping ablation triples (K6, K7, K9) unjoined. Under G2 there is no second pass. The
builder may start.
