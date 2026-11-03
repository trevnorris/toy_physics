# S15 reading item 4: an outside assessment of how it is computed (2026-11-03)

A read-only assessment the orchestrator obtained on 2026-11-03, after the user asked whether item 4 is computed correctly
and whether there is a faster way. It is advisory: no CAS run, kernel or computation stood behind it, every cost in it is
an estimate it labels as one, and nothing in it has been reviewed. It describes methods only and names no condition.

- Source: a fresh `codex exec` session (model gpt-6.1-sol, reasoning xhigh, sandbox read-only), session id `01a1a280-56e3-7da0-8d0b-f9924ffe67de`, 2026-11-03 17:40–17:52Z.
- Prompt sha256 `69f980d6ce9ae9086e311241cc9f7760a8c915071af84705a5fc2fb6434cb69b`. It handed the decision list, the r0–r3 dispositions, the repair-4 notes, the engines and their run harness, the item-4 timings of the current repair, and the round-4 review legs' derivation scripts; it asked five questions (is the construction computing the object item 4 names; why is it slow; is there a much cheaper exact construction; does that need a decision-list change; is the question itself wrong).
- The sandbox refused the answer file, so the answer below is the session's final message, copied from its log. Two edits: the path of the builder's SymPy working draft (outside the tree) is shortened to its file name, and the path of the refused answer file is dropped.

---

codex
The new constructions **can return exponent conditions**, so they have moved beyond the incapable attempts described in Y2. They do not yet establish delivery of the complete item-4 classification. The dominant concern is excessive symbolic representation and parameter branching; exactness itself does not require this approach.

I could not save `ANSWER.md` because the filesystem sandbox rejected the write. The assessment follows.

References use these source abbreviations:

- **WL**: [S15_separation_departure.wl](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl)
- **PY**: [S15_separation_departure.py](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py)
- **RUN**: [S15_separation_departure_run.py](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure_run.py)
- **DRAFT**: `sympy_six_hour_candidate.py` (the builder's SymPy working draft)
- **DL**: [S15_reading_decision_list.md](/var/projects/toy_physics/research/pde_ledger_v3/directives/S15_reading_decision_list.md)
- **R2/R3**: the corresponding build-review disposition files in `research/pde_ledger_v3/directives/_measurements/`.

**1. The target is right, but the present construction is not yet a secured complete answer.**

The target is an identity in separation magnitude **and direction**, at fixed instantaneous parameter values. It must preserve cancellations among amplitudes, rates and velocities. It is not an identity in independently variable amplitudes, and it concerns the geometric departure, not the departure with responses live (`DL:35–60,75–88`; `R3:90–101`).

Read from source: both engines eliminate direction identities on the unit sphere (`WL:1309–1323`; `PY:2085–2114`). Their preliminary sufficient constructions deliberately solve a stronger problem: annihilating physical monomials or matching independent formal coefficient coordinates. Exponents and radii are excluded from those solve variables (`WL:1324–1369`; `PY:39–45,104–157`). Those returns are sufficient subsets, not the complete requested classification. The sources describe that limitation correctly.

The additional constructions are materially different. Wolfram lifts Euler integrals, reconstructs the departure, generates exponent-sensitive necessary factors and recursively classifies parameter branches (`WL:660–682,1195–1264,1375–1388`). SymPy now has an additional exponent-sensitive construction too (`PY:797–817,1084–1206,1271–1349,2132–2138`). The historical Y2 diagnosis therefore should not simply be repeated against these paths.

Starting from closed `Q` is not inherently a change of object: item 3 itself is computed as the closed geometric derivative of that same `Q`, plus twice `Q` (`WL:1533–1538`; `PY:2267–2271`). Lifting and differentiating it is legitimate if all transformations preserve its branches and domain.

There is, however, a **conditional endpoint-substitution problem in SymPy**:

- Its factor solver can return a substitution for the upper radius (`PY:1094–1114`).
- Classification substitutes accumulated rules into `q` and inherited radial integrand/boundary expressions (`PY:1277–1287`).
- Radial integration still uses `self.lo,self.hi`; the cached radial moment also contains those original endpoints (`PY:831–847,850–870`).
- The finite-Laurent joint solver can consume those integrals and return before the later substitution of accumulated rules into `evaluated` (`PY:1235–1270,1307–1309`).
- Its polynomial walker assumes incoming equations are already restricted: it carries accumulated rules toward the answer but does not apply them to the incoming equations at entry (`PY:444–462`).

The inference is that this path solves an intermediate problem over the old symbolic interval while carrying a radius relation separately, rather than consistently restricting the whole integral. Subsequent solving can impose additional constraints or produce incompatible substitutions. I have not run a counterexample or established that a particular final object exercises this path; this is a source-level gap, not a measured false condition.

The draft retains the pattern (`DRAFT:838–853,1278–1294,1315`). Wolfram explicitly substitutes endpoints and establishes a branch-local moment cache (`WL:1190–1196`).

There is also a status-contract issue. Both engines can print `NOT OBTAINED EXACTLY` merely because their mathematical traversal returns unresolved branches (`WL:1390–1393`; `PY:2143–2148`; `DRAFT:2231–2234`). Examples include repeated input or failure to obtain a necessary polynomial (`WL:1189,1252–1256`). Those events are not themselves memory exhaustion or a six-hour stop, the only currently authorized reasons for such labels (`resume_s15_reading_repair_4g.md:7–17`).

The explicitly preserved historical failure labels are different and authorized: the launcher carries them into a reconstruction pass that skips the additional attempt (`RUN:33–64,212–217`; `WL:1370–1374`; `resume_s15_reading_repair_4g.md:11`).

**2. What grows is primarily branch count, expression size and artificial algebraic dimension. Series order can grow too.**

These mechanisms are visible in the source. Their individual shares of runtime were not measured in this consultation.

- **Primitive strata:** Wolfram recursively takes both zero and nonzero sides of every primitive boundary, pruning only when simplification proves a side false (`WL:1274–1288`). With \(p\) boundaries this admits up to \(2^p\) leaves before pruning. Further necessary-factor branches reclassify the restricted object (`WL:954–1010,1259–1264`). SymPy likewise partitions primitive boundaries and candidate substitutions recursively (`PY:722–744,1322–1338`).

- **Repeated expansion and rational normalization:** every Wolfram normalization expands the direction polynomial, reduces monomials, reconstructs it, then performs rational algebra in power coordinates (`WL:601–658`). Euler elimination repeats differentiation and normalization at every weight (`WL:753–763`). SymPy similarly expands powers, factors terms, combines denominators and cancels coefficients (`PY:663–715,876–890`). Common denominators and differentiation can greatly enlarge expressions without increasing the underlying number of physical scalar kernels.

- **Larger parameter problems than the physical problem:** Wolfram introduces coordinates for parameter powers, logarithms, amplitude monomials and then products of these coordinates (`WL:1081–1119`). Its final real reduction includes exponents, radius ratio/scale and product coordinates, followed by `LogicalExpand` (`WL:1120–1129,519–529`). This can be mathematically valid after exact pullback, but it solves a much larger ambient problem first. SymPy also branches on products and degenerate linear/quadratic pivots, with a Gröbner fallback (`PY:486–595`).

- **Increasing series order and matrix size:** Wolfram generates rows until there are as many pivots as columns (`WL:880–891`). SymPy repeatedly differentiates and computes rank toward the same full-column-rank target (`PY:1000–1019`). Neither loop has a separate rank-deficiency completion certificate. Thus the source provides no fixed order bound; a dependent retained basis would need detection. I have not established such dependence on the running object. SymPy subsequently constructs coefficient-block inverses and a dense weight determinant whose size grows with the selected indices across groups (`PY:1053–1081,608–626`).

The running stdout tail I read was on `reading/region/component_3/rate_labels_{drain1, drain2}` and contained these literal suffixes:

```text
/ item4 exponent-sensitive / PolynomialReduce normalization direction monomial
/ item4 exponent-sensitive / rational direction reconstruction residual
/ item4 exponent-sensitive / rational algebra in positive power coordinates
/ item4 exponent-sensitive / eliminate Euler integral 2 + item4Parameter55 + item4Parameter58
/ item4 exponent-sensitive / eliminate Euler integral 3 + item4Parameter55 + item4Parameter58
```

They correspond to `WL:607,615–616,656–657,758–760`. They show repeated normalization during Euler elimination. They do **not** establish that the present bottleneck is high-order series work.

The launcher preserves the earlier last step `eliminate series pivot 1` (`RUN:54–63`) and descriptions of coefficient-extraction and parameter-branch memory failures (`RUN:33–49`). These are recorded historical descriptions, not measurements I repeated.

**3. A substantially cheaper exact strategy is credible, but no complete-build speed claim has yet been established.**

The construction I recommend is:

1. **Compress the exact analytic object first.** Reduce the printed departure into a sparse direction basis once. Derive rotational/tensor structure from the actual object, collect identical scalar kernels and keep them factored before differentiating. Reuse exact kernel templates across Cartesian components within each engine.

   Reduce directions only; do not match velocity or amplitude coefficients as though those quantities were independently variable. For rates, construct the full instantaneous time derivative before imposing candidate relations. The present closed-rate constructors already retain that chain rule (`WL:440–455`; `PY:1990–2008`).

2. **Generate inexpensive exact necessary constraints.** Derive coefficient recurrences for the small scalar kernels instead of repeatedly differentiating large closed expressions.

   For example, derive the recurrence of
   \[
   G(z,\mu)=(1-2\mu z+z^2)^\gamma
   \]
   from
   \[
   (1-2\mu z+z^2)G'=\gamma(-2\mu+2z)G.
   \]
   Angular coefficients then use exact polynomial moments. Region coefficients use the exact analytic radial moment
   \[
   M(s)=\int_{\ell}^{h}r^s\,dr,
   \]
   retaining its logarithmic branch at the removable divisor. On a positive nonempty shell this moment is strictly positive.

   This avoids manufacturing many pole cases from quotient representations of regular integrals. Establish an exact bridge from these formulas to the printed closed departure in each engine.

3. **Certify the candidate cover, including resonances.** Low-order necessary coefficients can narrow the candidate parameter set rapidly. They cannot be treated as sufficient merely because several vanish.

   Mixed exponent families can cancel after index shifts. Where such possibilities remain, construct a small scalar differential system from the actual kernels. Lift scalar integrals once, use their homogeneous differential relations, and retain the integration constants removed by elimination. Endpoint powers satisfy first-order rational differential equations; logarithms give triangular extensions.

   Collect identical kernels and integer shifts before constructing cyclic rows. Compute necessary determinant/minor coefficient loci fraction-free, preserving every vanishing-pivot stratum. Solve those loci in the actual exponent coordinates and radius ratio. This supplies a finite necessary candidate cover without guessing a list of exponent resonances or using an enormous universal Taylor bound.

4. **Solve exponent strata before expanding amplitude alternatives.** On each exponent/radius stratum, solve the remaining small relations among amplitude, rate and velocity products. Substitute every product definition back into the named quantities. Keep final amplitude equations factored where possible; the contract does not require enumerating every overlapping zero-coordinate alternative.

   This must actually produce parameter relations. Printing an unsolved coefficient system, jet sum or aliased matrix condition would repeat X1 (`R2:102–115`).

5. **Verify every surviving branch on the full printed departure.** Then restrict the actual closed `Q` on that same branch. An expression that does not simplify to zero is not grounds to discard an entire parameter branch: refine it with another necessary equation or resolve it through the scalar system.

   Restrict integration endpoints and cached moments along with integrands. Preserve all divided-factor branches, pivot degeneracies and integration constants. Do not substitute a reference-separation template for the restricted `Q`.

The completeness argument is straightforward: every target parameter point lies in the exact necessary candidate cover; that cover is solved exhaustively; surviving branches are verified on the original identity; excluded branches have exact contradictions. This allows finite coefficients to serve as a filter **without asserting finite-truncation equivalence**.

There are also inexpensive exact improvements inside the existing machinery. SymPy’s weight matrix is a Vandermonde matrix in the displayed weight-minus-index values (`PY:612–625`). Construct its determinant as the product of pairwise differences, with an exact generic determinant identity, instead of computing a dense symbolic determinant. Likewise, the Wolfram spectral resultants are resultants of products of linear factors (`WL:893–896`); their factored products can be constructed from pairwise root differences. These operations compute the same necessary factors without supplying physical exponent conditions.

**Unmeasured cost estimate:** after kernel setup, seconds to a few minutes for small surface/single-family objects, and minutes to tens of minutes for harder mixed-rate objects, **if** the compact scalar ranks and candidate sets remain small. The work becomes short coefficient recurrences, low-rank elimination and exact candidate checks.

I have not measured those compact ranks or completed a candidate cover. A stubborn degenerate stratum could still dominate. Therefore I cannot responsibly promise a sub-hour runtime for every object. The sample review timings do not establish a complete-build runtime.

Conditions could be missed by dropping zero amplitudes or pivots, canceling exponent divisors, ignoring shifted exponent supports, matching velocities too strongly, losing integration constants, mishandling endpoints, or accepting numerical coverage. Exact substitution proves sufficiency; the exhaustive candidate-cover argument supplies necessity.

**4. This need not change the decision list.**

My reading of `DL:86–88` is that it forbids an expansion **in place of** an exact object. Exact coefficients used as necessary filters, followed by exhaustive exact parameter solving and verification of the full identity, do not replace the object.

An all-orders convergent representation with a proved tail identity also does not impose a small-radius regime: throughout the domain, \(r/d<1\). If the intended rule instead forbids any internal use of series even within an exact proof, that requires explicit clarification; it is stronger than the quoted wording. The scalar differential-system construction avoids depending on that interpretation.

A numerical step may locate candidates or suggest factors to investigate. I read that as allowed provided every printed condition is reconstructed and certified exactly, and completeness is proved independently of numerical coverage. Finding roots numerically and checking them exactly does not exclude missed isolated or multiple roots. No numerical locator is necessary for the construction above.

**5. Item 4 is meaningful; its wording does not force the current cost.**

It asks for a complete identity classification for a finite power-profile family. That is substantially stronger than finding representative exponent relations, but it does not require the present Cartesian expansion, exhaustive primitive traversal, large product-coordinate reduction or dense determinant construction.

Conversely, exactness does not guarantee a cheap classification. The proposed change should be described as a representation-based way to reduce cost until its complete coverage and runtime are demonstrated—not as an already measured seconds-per-object solution.

No CAS run, numerical computation, kernel, guarded job, repository edit or commit was performed. The unresolved empirical questions are whether particular final SymPy objects exercise the endpoint path and what compact ranks, candidate counts and complete runtimes the replacement achieves.

The attempted write of the answer file was rejected by the read-only sandbox and disabled approval settings; the file was not created.
hook: Stop
hook: Stop Completed
