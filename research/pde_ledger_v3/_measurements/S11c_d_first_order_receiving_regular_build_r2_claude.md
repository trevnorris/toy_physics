NEEDS REVISION

This is a source and build assessment only. I didn't run anything, and the verdict is not a runtime acceptance. I read worker.py, build-notes.txt, method-plan.txt and guide.txt in full. From the evidence tree (about 2,000 files) I opened only the force-column files, `complete-transformed-operator.json`, the geometric-dual operand file and a few others. I did not audit the remaining files, including the 400 local cells and the individual C3 entry views. The hand arithmetic below is a preview only.

## Blocker: the live chart and C3 are not joined to the saved proofs

The worker takes `E5`, `D5`, `transformed` and `physicalMatrix` from `complete-transformed-operator.json`. It then builds `C`, `R3` and `E3` from them. Three of the proofs that should tie these together are not bound to the live objects:

- **`old_matrix('receiving','geometric-dual',right=sp.eye(5))`** (worker.py:141) leaves `left=None`. The saved left operand is an unsimplified D5·E5-type product, with the `l²/20+1/400` denominators visible in the file. Nothing requires it to be the live `D5*E5`.
- **`old_matrix('receiving','LEFT-raw-vs-receiving',right=physicalMatrix)`** (line 148) leaves the left operand unchecked. Line 147 checks only the right side against `full['physicalMatrix']`. Nothing ties `transformed`, and therefore `C`, to `D5 · physicalMatrix · E5` or to the native equation basis.
- **Line 157** (`full['chart']==E5 and full['dual']==D5`) is a tautology, because both sides were assigned from `full` at line 134. It binds nothing.
- **`FULL-offwave-source-correspondence`** (line 146) binds only the left side, to `physicalMatrix`. The source-side operand is not joined to the live force or pressure rows.

As written, a saved dual or transformation proof about different operands would pass. The checks on `R3`, `E3`, the eta table and C3 would then rest on an unproven basis. The guide requires "same original source/receiving equation basis" and proof input/return joins.

**Fix:** bind each proof's left operand to the live object:

- For the dual, compare `left.applyfunc(sp.cancel)` to `(D5*E5).applyfunc(sp.cancel)`.
- For raw-vs-receiving, compare the left operand to `transformed`.
- For the full off-wave correspondence, bind its source-side operand to the live force and pressure rows.

This binds the proof to the live operands. It does not replay the proof or require its two sides to be structurally equal. Emit the comparisons and delete the line-157 check.

## What I checked and found sound

- **Row selection:**
  - `R3` and `E3` are checked literally against the plan's indices (lines 154–155), and `R3·U` is checked against `iQ/(5K2)`.
  - The hand projection table is compared against the full 5×2 `D5*Feta` before any rows are selected.
  - I recomputed `R3(p)·9U = 0` by hand from the saved `localRight`, and U_B's `R3·U` is zero.
- **Step and cusp:**
  - `A'=(9/2+4T)(1-T²)/L` and the endpoint jump of 9 are verified.
  - The Q identity is applied before the cusp.
  - The delta/delta-prime expansion follows `f δ' = f(p)δ' − f'(p)δ`, with the T1 column checked separately.
  - The dropped-`R3'U` control is declared at entry [0,0]. By hand, `R3(p)B'_A = −i/(5K2)`, so it should respond.
- **Pressure and sigma:**
  - There is no double count: `physicalPathForce` is documented as excluding the pressure addition, and the pressure force is added separately.
  - The chemical tag is transported once, as `10·(A_mem/rho)` times the physical tag.
  - Both normal signs are checked per face against the saved flat factor.
  - By hand, the saved sigma row 0 vanishes at both T=±1, so it is localized as the worker requires.
- **Domains and transport:**
  - Base transport and the even-sheet map are consistent with the data. The saved `matrix` entries contain only `l**2` powers; the `l**3` terms sit in `physicalMatrix`.
  - Every negative-power base is traversed before cancel, and K2 is checked separately.
  - The raw and cancelled determinant denominators are checked on both rays.
  - The real-or-imaginary-only case (`b=0`) is handled correctly by `gcdex` and the Sturm step, and the endpoints at 0, p and ∞ are covered.
- **Current:**
  - The six old assignments are matched by AST and never executed.
  - The bra/ket momentum slots are consistent with the G0/Gref frames.
  - Both cross contractions and the conjugate relation are tested, and the wrong-bra control uses the same contraction.
- **Growth ledger:** R3, the eta multiplier, the pressure factors and E3 add positive powers with no decay credit. Nondecaying transverse components are emitted as not claimed.

## Non-blocking notes

- `localized()` is applied to all five raw sigma rows. That is stricter than the selected three-row class needs, and it fails if any row has a nonzero endpoint. I verified row 0 only, so a refusal from another row is possible.
- The controls will refuse rather than report if the cross-current contraction is silent or the declared entries come out zero. That behavior is intended.

I have no other substantive blockers.