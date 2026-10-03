# Method assessment: corrected weak ends versus the selected uniform subspace

## Packet coverage

**Read in full:**
- `method.md`
- `evidence-guide.md`
- `packet-index.json`
- `uniform/plan.md`
- `uniform/opaque-interface-receipts.json`
- `ends/method.md`
- `controls/translated-weak-end-conclusion.json`
- `ends/depth-domain.json`

**Read in part:**
- `source/uniform-worker.py`: the function index, `selected_lift`, `EndBinding`, `end_restriction`, `on_wave`, `uniform_applicability` and `exact_limits` (about lines 205–405 and 620–720).
- `source/ends-worker.py`: the function index plus a keyword grep.

**Not read:**
- The 200 cells in `ends/symbols.json`, the control records and `exact-evidence.jsonl`.
- `uniform/selected-result.json`.
- `engine.py`, `c1.py` and `c2.py`.
- `uniform-continuation.py` and the rest of `ends-worker.py`.
- Any opaque payload. I did not decode any, and I make no claim about their contents.

I did not read the review prompt file separately. Its text was in the request.

Findings about the E cells below come from the method text and the worker source, not from inspecting cell contents.

## Blocking defects

### B1. A nonzero A or R cannot be interpreted, and the method calls it "substantive"

- **Where:** `method.md` lines 89–113 and 99–106 (the "Two distinctions" section).
- **Why it matters:** E keeps grades eta^a sigma^b with a,b ∈ {0,1}. The excluded eta², sigma² and higher terms are preserved only as a note in `ends/method.md`. `EndBinding.retained` in `uniform-worker.py` (lines 305–310) also truncates to a,b < 2. `EndBinding.bind` substitutes the numeric `origin` into the saved closed pencil, so P_old is a finite law with all grades.
- **Failure mechanism:** Generically Ephys − P_old carries an O(eta²) excluded remainder. A and R can then be nonzero purely from truncation, with the retained grades agreeing. The method says nonzero A or R "is a substantive applicability finding". It also makes the grade-attributed comparison optional and forbids "regrading". Under those rules a truncation remainder is indistinguishable from a real mismatch. The method then says to stop. It also gives no way to report a retained-grade agreement.
- **Required revision:**
  - Say that a nonzero A or R is first a finding about the finite law.
  - Require the grade attribution: a coefficient-by-coefficient comparison of E_ab L_old with the eta^a sigma^b coefficients of P_old L.
  - Require it whenever raw unbound operands of P_old are restorable. Otherwise report `NONZERO_UNATTRIBUTED`.
  - Do not call a nonzero result an applicability failure unless a retained-grade entry differs.
  - Equality of the retained grades alone must never be reported as "A=R=0".

### B2. The row, leg and momentum-sign join does not allow source-derived weak-form factors

- **Where:** `method.md` lines 51–66, against the bilinear form on lines 17–24.
- **What E is:** E is the symbol of B(v,u) = 2π∫ v̂(−p)ᵀ E(p) û(p) dp. Its cells are built from consumer coefficient × source coefficient × wave × normal × response (`ends-worker.py` around lines 525–560). That is a test-row-labelled weak form.
- **What P_old is:** P_old is `CLOSED_PENCIL_LEGS[0]` or `[1]`, a signed-leg physical pencil bound with the leg convention (k, q) versus (−k, −q) (`EndBinding.native_leg`).
- **Failure mechanism:** The consumer rows may differ from the strong rows by source-derived factors. These could include integration-by-parts factors such as (ip)^n, consumer coefficient weights, transposition, the p→−p reversal on v, or conjugation. The method permits only a permutation or a constant field/row unit conversion. It explicitly forbids a "momentum-dependent row operation". It also never says which leg, which sign of p, or which time and phase convention is used for Ephys(p) versus P_old. A correct but p-dependent native weak-row factor would produce a false `NONZERO`. A wrong leg or sign choice could produce a false or accidental zero.
- **Required revision:**
  - Name the leg and the sign of p, and derive them from the `ends/phase-arguments.json` phase convention and the uniform leg map. Do not leave them to the worker.
  - Allow a left row map only if it is derived from the actual native weak-form definition. It may be p-dependent if the source dictates that.
  - Preserve that map and its inverse, and bar fitting it.
  - Make "no join, stop" the explicit outcome when no source-derived map reproduces the row labels.

### B3. A and R are not independent tests, and the inherited targets are grazing points

- **Where:** `method.md` lines 94–98 and 115–147 (the A/R definitions and the grazing paragraph).
- **A versus R:** A = (E − P)L and R = E L − L D. So A − R = −(P L − L D), which is the old invariant residual. The method restores that residual as zero on the wave. On the nongrazing surface A and R are then the same test, not two. The method should say so, and certify A − R against the restored residual as a consistency check. Treating them as independent evidence overstates the coverage.
- **Grazing targets:** The listed targets are LEFT cs = √(3/2), p² = 595/100 and RIGHT cs² = 150/101, p² = 601/100. They satisfy 9/cs² − 1/20 = p² exactly, so q = 0. Every inherited correspondence point is therefore at grazing. The old limit operands in `exact_limits` (lines 651–715) are evaluated only along two paths at the branch momentum k = ±√k23 (RADIATING and EVANESCENT). The nongrazing rational certificate on a (p, cs) surface says nothing about these points unless it is joined through continuity of the closed E and regularity of L and D. D is q-free by the worker check at line 358.
- **Required revision:**
  - State that the physical content of the selected branch values lives entirely in the grazing comparison.
  - State that the grazing claim is limited to those inherited paths and points, and is not a surface claim.
  - Prefer evaluating R at q = 0 directly from closed E, L and D when L and D are regular there. R does not involve the singular P. Use P-limits only for A.

## Nonblocking suggestions

1. **Certificates (lines 123–130).**
   - The existing `on_wave` uses `Poly(..., domain='EX')`. That is the broad domain the method warns against, and it must not be inherited unchanged.
   - Use a small exact domain over Q(i), with p and cs as rational generators, and q as the only reduced variable.
   - Check that every actual denominator, including p- and cs-dependent factors and the Gram determinant, is nonzero on the stated domain. q+β is nonzero on the closed first quadrant by `depth-domain.json`.
2. **Control baseline (lines 168–193).**
   - A control that perturbs E and expects a "responsive exact residual" needs a recorded unperturbed baseline for the same point. If the baseline A is nonzero, movement is ambiguous.
   - Controls 1 and 3 may be unavailable on the decoupled transverse lift, as the method already notes. The coverage-gap path is acceptable.
3. **Cs scope.**
   - E is a symbol over cs ∈ [1,2]. The old selected branch is independent of cs by the worker check at lines 368–372, so the cs and branch-point join is a statement at specific points.
   - Do not phrase a successful result as holding over the cs interval unless the certificate is an identity in cs.
4. **Conditionality.**
   - `uniform_applicability` (lines 641–643) records the old selected results as `CONDITIONAL_ON_UNSUPPLIED_NONUNIFORM_AND_DIRECT_MIXED_GRADE_COMPOSITION`.
   - Say explicitly what a pass does to that status: it can remove the condition only for the selected subspace at real ω = 3, and nowhere else.
5. **Opaque blobs.**
   - The method already correctly refuses to decode them. Add that the worker must verify the receipt hashes before restoring any blob, and that it must bind omega = 3 and cs from the actual input pickles before reuse. The method says this loosely at lines 60–61.

## Answers to the five requested questions

1. **Joins.** The sufficiency of the proposed joins fails on one point, B2. The sections on ends (0 and W/2), face and end distinction, memory β = (30+9i)/109, p→k, physical q versus algebraic Q, and actual maps are good. They do not specify the leg, the sign of p or the weak-row map. Fitted scales are correctly barred, but the barred class is too narrow.
2. **Eta/sigma assembly, then finite binding.** This is correctly specified, and Delta, A and R are correctly allowed to be nonzero. The gap is that nonzero results are not attributable to grades (B1), and A and R are the same test (B3).
3. **Rational certificates and grazing.** Preserving denominators and refusing raw singular substitution is correct. The domain and Poly caveats are in item 1 above. The grazing extension is restricted by B3.
4. **Scope of success.** The disclaimers at lines 149–166 are accurate and sufficient. They exclude a frequency neighborhood, a new current identity, a mode census, plane-wave action and loss. A pass shows only that span(L) is E-invariant with D_old at real ω = 3 on the stated domain and in the stated coordinate frame, subject to the row map found under B2.
5. **Controls and stop handling.**
   - The three controls are plausible. Their applicability is uncertain and the gap rule covers it.
   - Preserving completed returns is correct. The late SymPy Zero serialization failure is correctly separated from method acceptance.
   - Stopping on a substantive mismatch is correct once B1 defines what "substantive" means. The stop must still persist the E_ab L_old diagnostics and A, R and Delta.

NEEDS REVISION