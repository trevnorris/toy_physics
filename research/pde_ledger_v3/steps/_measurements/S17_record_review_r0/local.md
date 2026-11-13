# S17 independent record review, round 0

**Verdict: revisions required for source fidelity.** I reviewed the S17 step record with SHA-256 `1e4802b348ab34fb9baf1fafdffecc1babbf9900c47d131f1909c67b821e8c1a`. I could inspect all seven local main TeX files, the Papers List, Stage 029, the S17 plan and the accepted S15/S16 records at `438d270c`. I read the decisive primary sections, including the theorem statements, matching calculations and limitations, and consulted the accompanying local source lookup. Findings concern the record only. This was a comparison of written sources and included no rederivation or Zenodo lookup.

1. **2.5PN, line 45.** The premise list ends with “and nonzero static orbital overlap. Under those conditions”. The theorem explicitly requires the overlap to “the static orbital/source overlap remains on the natural branch” with `m̂₀ = 1 + O(a²/r²)` and `m̂₀ ≠ 0` (2.5PN TeX 2469–2475, literal text below). An overlap tending to zero with the size could be nonzero at each finite size, yet would not satisfy this premise for a universal point-particle term. Give the natural-branch scaling. Preserve line 47's separate test of the absolute canonical strength `m̂₀²Γ₅ = 2G/(5c⁵)`.
2. **3PN, line 53.** “complete the static term against the imported geometry target” assigns the wrong provenance to the geometry coefficient. The paper imports the **full generic-frame static target**. The grouped response fit already predicts a static contribution. “Their exact difference is the unique geometry coefficient” (3PN TeX 8533–8545). Say that subtracting the already-fixed grouped contribution from the imported static target assigns the remaining part to geometry within this construction. The throat's physical geometry response is still open.
3. **4PN, lines 57–59.** “It uses the Schwarzschild one-body condition and reconstructs the local interaction against the comparable-mass reference” captures the reduced local result. The open-input paragraph ends with “A measured 2.5PN amplitude alone leaves that transport test open.” This accounts for the hereditary measurement but omits the paper's separate physical requirement for the **conservative** response. In the authors' words, “the relevant grouped bundle / still has to be computed as an actual branch output” (4PN TeX 3541–3553). The three quadrupole components must be obtained and their isotropy tested. Failure of that test rules out the assumed branch before the amplitude comparison (3579–3598). Add this to the moving-drain handoff for the 4PN paper. Line 53's 3PN response debt and line 45's assumed rotational invariance are consistent with it. The 4PN canonical coefficient bridge and the separate transport test as stated for the canonical branch are supported.
4. **S15 handoff, line 63.** “the classification of profiles for which each departure vanishes at every separation” expands the deferred item to all departures. Item 4 covers only the **closed-form geometric** departures on the power-law profile family. S15 says “The three **net** family departures retain general-field integrals. The restricted closed-piece results cannot classify their sum” (S15 record 25). The accepted record expressly excludes 192 general-field cases from item 4. The live field-response contribution to the total crossing is a separate physical question. Limit the description of the deferral to the closed geometric pieces, or simply state that item 4 stays deferred.

The other principal distinctions survive this comparison. The Newtonian potential has a source calibration and the receiving body has a supplied mass and a boundary-flux assumption. S15/S16 leave the relation to free-mouth acceleration open. The first paper matches the pressure–volume remainder to the orbit. The full 1PN paper applies the required response slope before predicting an internal energy partition. The EIH velocity targets and the 2PN and 3PN reference matches are identified as inputs. The two conditional results retain different observables, local dissipative reaction at 2.5PN and a conservative local-plus-hereditary result at 4PN. The seven listed title/DOI pairs match the pinned README, and none of the main TeX files contains its own DOI.

## Reproducible text lookups

Commands run from the repository root. Primary sources use the pinned revision. Candidate excerpts come from the working record identified by its SHA-256. The text after each command is its literal stdout.

```bash
sha256sum research/pde_ledger_v3/steps/S17_post_newtonian_ladder.md
```

```text
1e4802b348ab34fb9baf1fafdffecc1babbf9900c47d131f1909c67b821e8c1a  research/pde_ledger_v3/steps/S17_post_newtonian_ladder.md
```

```bash
nl -ba research/pde_ledger_v3/steps/S17_post_newtonian_ladder.md | sed -n '45p;47p;53p;57p;59p;63p'
```

```text
    45	**Result and assumptions.** The dissipative result concerns two nonspinning bodies in the strict small-body limit. Its channel argument assumes no additional direct scalar source, compact outgoing scalar and vector channels with the stated low-frequency thresholds, rotational invariance, a passive outgoing quadrupole, and nonzero static orbital overlap. Under those conditions the lower scalar and dipole terms are suppressed. The orbital symmetric trace-free quadrupole supplies the local Burke–Thorne radiation-reaction form in the Iyer–Will family. [Channel argument](../../4d_2_5pn/paper/4d_2_5pn.tex#L35) · [Explicit theorem assumptions](../../4d_2_5pn/paper/4d_2_5pn.tex#L2456)
    47	**Matching and open input.** The absolute coefficient in canonical source units is `γ_quad^eff = m̂₀² Γ₅`. The GR reaction strength requires `γ_quad^eff = 2G/(5c⁵)`. The actual moving drains must realize the branch and its amplitude. The reduced-model normalization comparison belongs to S18. [Burke–Thorne target](../../4d_2_5pn/paper/4d_2_5pn.tex#L791) · [Canonical coefficient](../../4d_2_5pn/paper/4d_2_5pn.tex#L2552) · [Assigned DOI](https://doi.org/10.5281/zenodo.19492270)
    53	**Matching and open input.** The Schwarzschild one-body limit fixes three coefficients. For comparable masses the authors invert the grouped response map on nine known GR residual coefficients and complete the static term against the imported geometry target. The physical response represented by that map and the throat geometry and compressibility that would supply the static contribution remain to be derived. [Nine-coefficient match](../../4d_3pn/paper/4d_3pn.tex#L8331) · [Imported static coefficient](../../4d_3pn/paper/4d_3pn.tex#L8533) · [Open physical origins](../../4d_3pn/paper/4d_3pn.tex#L8969) · [Assigned DOI](https://doi.org/10.5281/zenodo.19501724)
    57	**Result and assumptions.** The local, instantaneous calculation gives an ordinary two-body Lagrangian in a generic frame within the paper's chosen representation. It uses the Schwarzschild one-body condition and reconstructs the local interaction against the comparable-mass reference. For the history-dependent piece, the paper adopts the standard conservative GR tail functional. Its canonical coefficient satisfies `C_tail = (GM/2c³) γ_quad^eff`. [Local reconstruction](../../4d_4pn/paper/4d_4pn.tex#L3420) · [Choice of tail functional](../../4d_4pn/paper/4d_4pn.tex#L2169)
    59	**Matching and open input.** The canonical result gives `C_tail = G²M/(5c⁸)` if the 2.5PN branch supplies `γ_quad^eff = 2G/(5c⁵)`. This is the same quadrupole normalization. For a toy-medium measurement the paper separately writes `α_tail^toy = Θ_tail (GM/2c_s³) γ_toy`. Agreement with the canonical bridge requires `c_s = c`, the same canonically normalized source, and a monopole-scattering transport factor `Θ_tail = 1`. A measured 2.5PN amplitude alone leaves that transport test open. [Transport factor and branch](../../4d_4pn/paper/4d_4pn.tex#L2368) · [Canonical normalization](../../4d_4pn/paper/4d_4pn.tex#L2406) · [GR tail coefficient](../../4d_4pn/paper/4d_4pn.tex#L3486) · [Assigned DOI](https://doi.org/10.5281/zenodo.19561056)
    63	For the Newtonian orbit benchmark, compute the motion of each mouth from the evolving fields and compare its companion-induced acceleration to the reduced result. For the PN benchmarks, record which coefficients the simulation supplies before comparing with each paper's target. This distinction matters even for an exact algebraic match to the published Hamiltonian. S21 retains the conditional shell-crossing comparison. S15 item 4, the classification of profiles for which each departure vanishes at every separation, stays deferred. S18 owns the reduced radiation-reaction normalization check. S22 inventories the imported `κ_add`, `κ_PV`, `κ_ρ` and `n = 5` with their individual provenances. This citation step adds no register entry and closes no physical requirement. [S15 record](S15_two_drain_region.md#L47) · [S18 remit](../V3_STEP_PLAN.md#L760) · [S22 inputs](../V3_STEP_PLAN.md#L1303)
```

```bash
git show 438d270c:research/4d_2_5pn/paper/4d_2_5pn.tex | nl -ba | sed -n '2467,2476p;2487,2490p;2552,2561p'
```

```text
  2467	  \item the low-frequency defect/worldtube sector is rotationally invariant,
  2468	  \item the surviving STF quadrupole branch is passive and outgoing,
  2469	  \item the strict small-body hierarchy remains valid,
  2470	  \item and the static orbital/source overlap remains on the natural branch,
  2471	  \[
  2472	  \hat m_0 = 1 + O\!\left(\frac{a^2}{r^2}\right),
  2473	  \qquad
  2474	  \hat m_0\neq 0.
  2475	  \]
  2476	\end{enumerate}
  2487	where $\mathbf a_{2.5\rm PN}^{\rm quad}$ is the standard local
  2488	Burke--Thorne / Iyer--Will quadrupolar reaction term generated by the
  2489	orbital/worldtube STF quadrupole, while all other currently derived odd
  2490	channels are pushed above universal point-particle $2.5$PN.
  2552	and the retarded STF matching coefficient on the same branch is $\Gamma_5$, then
  2553	the canonically normalized odd coefficient is controlled only by the product
  2554	\[
  2555	\gamma_{\rm quad}^{\rm eff}=\hat m_0^2\Gamma_5.
  2556	\]
  2557	Matching to Burke--Thorne is therefore equivalent to the scalar condition
  2558	\[
  2559	\boxed{
  2560	\hat m_0^2\Gamma_5=\frac{2G}{5c^5}.
  2561	}
```

```bash
git show 438d270c:research/4d_3pn/paper/4d_3pn.tex | nl -ba | sed -n '8333,8338p;8533,8545p'
```

```text
  8333	
  8334	Now substitute the solved GR middle-block residual
  8335	\[
  8336	(L_6,\dots,L_{14})=(\Delta l_6,\dots,\Delta l_{14}).
  8337	\]
  8338	The inverse map then gives one unique grouped coefficient set:
  8533	The same conclusion can be read directly on the target coefficients. The imported generic-frame static target coefficient is
  8534	\[
  8535	c_U^{\rm target}=-\frac{227}{24}+\frac{21\pi^2}{32},
  8536	\]
  8537	while the richer grouped-$P_2$ middle-block fit already forces the static companion
  8538	\[
  8539	c_U^{(P_2\,\rm pred)}=\frac{293-308\nu-102\nu^2}{24}.
  8540	\]
  8541	Their exact difference is the unique geometry coefficient
  8542	\[
  8543	\boxed{
  8544	c_U^{(g)}=
  8545	\frac{408\nu^2+1232\nu-2080+63\pi^2}{96}.
```

```bash
git show 438d270c:research/4d_4pn/paper/4d_4pn.tex | nl -ba | sed -n '2381,2397p;2406,2415p;3541,3553p;3579,3585p;3594,3598p'
```

```text
  2381	On the branch with $c_s=c$ and
  2382	\begin{equation}
  2383	\gamma_{\rm toy}=
  2384	\gamma_{\rm quad}^{\rm eff},
  2385	\label{eq:4pn-tail-toy-canonical-branch}
  2386	\end{equation}
  2387	Eq.~\eqref{eq:4pn-tail-toy-transport} reduces to the exact bridge
  2388	\eqref{eq:4pn-tail-exact-bridge} when
  2389	\begin{equation}
  2390	\Theta_{\rm tail}=1.
  2391	\label{eq:4pn-tail-theta-one}
  2392	\end{equation}
  2393	This is the sharpest way to say what is and is not still open.  The
  2394	representation/source problem has already been solved by the $2.5$PN quadrupole
  2395	audit.  The local $4$PN existence problem has already been solved by the local
  2396	referee chain.  The only possible remaining hereditary mismatch is one scalar
  2397	monopole-scattering transport factor.
  2406	  \item the unique compatible conservative hereditary coefficient is the bridge
  2407	  coefficient \eqref{eq:4pn-tail-exact-bridge},
  2408	  \item therefore the $4$PN tail sector introduces no new independent
  2409	  quadrupole-normalization datum,
  2410	  \item and full GR matching of the hereditary coefficient is equivalent to the
  2411	  already-known $2.5$PN normalization target
  2412	  \begin{equation}
  2413	  \gamma_{\rm quad}^{\rm eff}
  2414	  =
  2415	  \frac{2G}{5c^5}.
  3541	Third, the deeper microscopic derivation of the grouped conservative quadrupole
  3542	data also remains open.  On the moving-throat side, the relevant grouped bundle
  3543	still has to be computed as an actual branch output,
  3544	\begin{equation}
  3545	D_A(\omega)=D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+O(\omega^6),
  3546	\qquad
  3547	N_A(\omega)=N_{A0}+N_{A2}\omega^2+N_{A4}\omega^4+O(\omega^6),
  3548	\label{eq:conclusions-grouped-bundle-open}
  3549	\end{equation}
  3550	with $A\in\{20,21,22\}$ and with the isotropic branch determined rather than
  3551	assumed.  These quantities are not leftover algebraic residues in the present
  3552	$4$PN theorem.  They are dynamical response data of the microscopic
  3553	mouth/support/mixed-sector system.
  3579	On the conservative side, the physical branch must first determine whether the
  3580	natural grouped real $P_2$ isotropy gate actually holds,
  3581	\begin{equation}
  3582	\boxed{a_2=0,
  3583	\qquad
  3584	b_2=0,}
  3585	\label{eq:conclusions-isotropy-gate-4pn}
  3594	If these first gates fail, then the minimal isotropic passive/outgoing
  3595	quadrupole branch is already dead before the final normalization test is even
  3596	attempted.  If they pass, then the next task is simply to evaluate the single
  3597	invariant normalization product in
  3598	\eqref{eq:conclusions-normalization-target-4pn}.
```

```bash
git show 438d270c:research/pde_ledger_v3/steps/S15_two_drain_region.md | nl -ba | sed -n '20p;25p;47p'
```

```text
    20	E_Q = 𝒟Q + 2Q = E_Q^geo + E_Q^resp .
    25	For the family `u_k=A_k(t)r^{a_k}`, `ρ_k=B_k(t)r^{b_k}`, the amplitudes remain time-dependent and the exponents are real constants. One closed example is the drain-2-only contribution to centre transport through the inner surface, Cartesian component 1. Its value is `−(4π/3) A₂B₂ r_in^(a₂+b₂+2) Ẋ₂¹`. The corresponding geometric departure is `−(8π/3) A₂B₂ r_in^(a₂+b₂+2) Ẋ₂¹`. Here the distance is the local inner radius. Neither value contains the inter-drain distance. This piece depends on drain 2's own motion and can be present without a companion-induced response. The three **net** family departures retain general-field integrals. The restricted closed-piece results cannot classify their sum. [Both native expressions, M2, and object 080, M3](_measurements/S15_record_measurements.md)
    47	**Item 4 remains deferred.** Wolfram reports 96 Reduce returns, 42 exponent-sensitive returns and six sufficient-only cases, with 192 general-field cases inapplicable. The six concern the two-drain-labelled region momentum and its rate, in each Cartesian component. Their printed sufficient conditions and corresponding restricted quantities survive as one-engine evidence. SymPy attempted no item-4 case in the accepted build. The stopped independent trial assigns no final classification, and the 144-object independent survey never started. The all-separation profile classification and the actual field responses remain open. [Per-object Wolfram returns and restrictions, M3–M4](_measurements/S15_record_measurements.md) · [Accepted reading build](../directives/_measurements/S15_reading_build_review_r5_disposition.md) · [Deferral](../directives/S15_reading_comparator_decision_list.md)
```

```bash
git show 438d270c:research/pde_ledger_v3/directives/S15_reading_decision_list.md | nl -ba | sed -n '75,83p;105,107p'
```

```text
    75	3. `E_Q^geo` on the profile family, with its dependence on `d` explicit and exact. Where the integrand holds no general
    76	   field besides the family's amplitudes, the region's quantities and the velocities, the integrals are evaluated in
    77	   closed form. Otherwise the object is printed with the family substituted and left unevaluated.
    78	4. For each `Q` whose item-3 object is in closed form:
    79	   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
    80	     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
    81	     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
    82	   - `Q` on the family under each condition.
    83	5. For each drain `k`, the part of `Q` that depends on that drain's velocity: `Q` minus `Q` at `Ẋ_k = 0`.
   105	  - for each `Q`, the printed conditions under which `E_Q^geo` vanishes. These are the conditions under which `Q`
   106	    falls off as the inverse square when no quantity other than `X_1` changes along the family. Beside them it reads
   107	    `E_Q^resp`, what any change adds. Neither is read as a condition for `E_Q = 0` while the changes are live;
```
