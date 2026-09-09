1. **The missed issue: an energy-basis identity is not automatically a source-operand rewrite.**

   WL differentiates only atoms registered as field jets. Its `energyCoefficient<i>` symbols are declared coefficients, so their spatial derivatives are zero. It then computes `muE = el[energy["DENSITY"]]`. SymPy instead inserts the background-dependent local thickness and differentiates it inside the energy construction. See [WL differential algebra and EL](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:75), [WL coefficient allocation](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:222), and [SymPy local thickness](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1626).

   Map 4 sends some of those constant WL coefficients to **spatially varying expressions**. Consequently,
   \[
   \operatorname{EL}(T L)\ne T(\operatorname{EL}L)
   \]
   in general. Applying the correct density-level table to an already differentiated source misses product-rule terms.

   Here is a concrete constructor-level counterexample using **E04/E14**, with the independent E14 gamma coefficient set to zero solely to isolate the κ contribution. In one spatial dimension, write \(W=W_{\rm bg}\), \(R=W_0/W\), \(e=e_W\), and \(k=\kappa_{\theta W}\):
   \[
   L_{\rm WL}=a\,\theta'e'+b\,eW'\theta',
   \qquad T(a)=kR,\quad T(b)=-kR/W.
   \]
   The mapped density is exactly \(k\theta'(Re)'\), the SymPy local-thickness expression. But differentiating before versus after mapping gives
   \[
   \operatorname{EL}_\theta(TL)-T(\operatorname{EL}_\theta L)
   =\frac{kW_0}{W^2}W'e'
     -\frac{2kW_0}{W^3}e(W')^2.
   \]
   **The first term survives at \(\sigma_W^1\).** I verified this identity with a small symbolic calculation. The row definitions are in [map 4](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:105).

   This is **not** a finding that the measured cross-engine leaves genuinely disagree; I have not run their collapse. It establishes that the proposed post-construction bridge is insufficient even with grading repaired. Rank/invertibility over \(\mathbb Q(W_0,W_{\rm bg})\) checks an algebraic basis change, not compatibility with differentiation.

2. **The proposed grading fixes need a stronger formulation too.**

   Replacing \(R_W\) by its background value separately in each grade is unsound. E02 already demonstrates why:
   \[
   W_0(1+\eta w_1)c,\qquad c\mapsto C/(1+\eta w_1).
   \]
   Mapping before grading gives \(W_0C\), whose \(\eta^1\) coefficient vanishes. Replacing \(c\mapsto C\) in the existing \(\eta^1\) leaf leaves \(W_0Cw_1\).

   A sound graded bridge generally acts on the **whole retained coefficient tuple**. For a simple multiplicative map,
   \[
   (TF)_{1j}=T_{00}F_{1j}+T_{10}F_{0j}
   \]
   when \(T\) depends only on η. These contributions belong to the same final grade. This is legitimate coefficient convolution, not identifying η with σ or cancelling different final grades.

   Thus “per-grade images” is sound only if it means the correctly induced tuple transformation, including lower-grade contributions. The current prohibition on accessing another grade would need revision. Merely expanding a newly inserted \(W_{\rm bg}\) on isolated leaves cannot implement that transformation. See the [current grading and tripwire requirements](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:227).

   **Even that repair addresses grading only; it does not restore the missing coefficient derivatives above.**

3. **The cleaner construction is upstream of EL, with an explicit coefficient-field convention.**

   First decide what “the same constitutive law” means: which coefficients are constant, and which acquire background dependence under the thickness convention change. Two constant-coefficient parameterizations related by a spatially varying table need not describe the same restricted family of operators.

   Then apply the primitive thickness/basis transformation **before constitutive differentiation**, retaining the live background and its jets through all required derivatives. Let each engine independently construct its material amplitude, predicted amplitude, velocity, and source; expand profiles and extract grades afterward. Equivalently, a bridge could carry coefficient jets and their chain rules, but the current frozen outputs were constructed without those coefficient derivatives.

   There is an essential claim boundary: **a canonicalized replay validates the replayed constructions. It does not retroactively reconcile the original emitted streams** unless their relation to that replay is also established. “Absorb thickness upstream” is a cleaner experiment, but it changes more than the comparator.

   If further work is funded, I would start with a bounded diagnostic:

   - The E04/E14 sector, checking both density and EL transport.
   - Φ and its domain/control census, which do not require map 4 at all—the [typed scope explicitly excludes Φ from map 4](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:196).

   That would answer whether a larger instrument is well posed before another full build/review cycle. C’s coefficient rows plus contraction matching alone would miss the independent Φ question.

   A **joint numeric/PIT channel** could make the eventual comparison cheaper, particularly by evaluating circuits without expanding giant residuals. It cannot supply the missing convention bridge: samples must represent the same coefficient fields **and their derivatives**, points, and retained grades. Sampling frozen coefficient values alone reproduces the product-rule problem.

   A convention-invariant response comparison is another possible endpoint, but requires transforming inputs and output rows consistently and confronting the open omitted-block issue. It is unlikely to be the cheap route here.

4. **Reading B is sufficient for its stated claim; it is not an unconditional downstream correctness certificate.**

   The committed resolution establishes faithful implementation of the declared Φ, with PIT qualification. It explicitly leaves Φ’s physical correctness, velocity transformation, and extracted-block leakage open, and explicitly says **cross-engine agreement remains owed**. See [RESOLVED, caveats and downstream guidance](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md:36).

   S11c-d can perform calculations conditional on a nominated closed operator and its conventions. Binding profiles to that operator does not validate its coefficients. The intended downstream objects include gradient-driven mixing and leakage, so a retained gradient-term discrepancy cannot be dismissed simply because covariance holds. See [S11c-d’s scope and order counting](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:83).

   I therefore would **not** justify B by saying “c2 already has everything it needs.” Justify it as an explicit decision to carry unresolved corroboration while preserving the conditional result.

5. **The traps differ across the paths.**

   - **A:** The promised “definitive channel-(b) answer” is overstated. Besides the differential-stage problem, agreement remains limited to matched operands under the bridge; unmatched families and the three premise caveats survive. The [vetted question already imposes those limits](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:183).
   - **B:** “Representational-difference-unadjudicated” must not quietly become “known representational difference.” The unresolved alternatives include inconsistent constitutive coefficient-field conventions and actual implementation errors. Two internally covariant constructions can still produce different physical responses.
   - **C:** Verifying the normalization table and contraction identities validates the dictionary’s algebra. It does not establish source transport through differentiation or corroborate Φ.

   Conversely, B does not erase the independent per-engine constructions and their controls. Those remain useful evidence. The *additional cross-engine subtraction of their vanishing residuals* supplies no operand-agreement evidence.

I would record the disposition approximately as:

> Per-engine Reading B remains adopted, with its PIT qualification and three premise caveats. Matched cross-engine carrier, source, and Φ operand agreement remains UNDECIDED. The proposed collapse bridge has unresolved differential and grading-stage compatibility; current residuals establish neither representational reconciliation nor genuine disagreement. Cross-engine corroboration remains a carried obligation.

Choosing B also explicitly revises the [previously adopted A scope](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:219). Record the reason as **an inadequately specified comparison and a scope decision**, rather than computational intractability. The next worthwhile expenditure is the small upstream differential-compatibility diagnostic—not another grading-only v4.
