# S17 decision list, local independent source-fidelity pass 1

**Verdict. Findings on the pre-record S17 decision list.**

I read the candidate and adjacent source retrieval in the working tree, and consulted the seven main TeX manuscripts, README Papers List, Stage 029, the S15 and S16 records and the step plan at `438d270c0fa5576167ab237c61301f1435ef7800`. I inspected the claims and the outcomes they describe, including manuscript sections below the abstracts. This pass concerns decisions for a citation-only record. It supplies no independent calculation of the orbital physics and makes no disposition of the future S17 record.

I checked all seven DOI assignments against the pinned README and searched each paired main TeX source for its assigned DOI. Each assignment matches and each literal count is zero. The adjacent retrieval reproduced its literal stdout byte for byte. The first paper’s test-mass static limit and orbital fit, the bridge paper’s energy coefficient, the full 1PN paper’s protocol-dependent breathing coefficient, the 2PN residual matching, the 2.5/4PN normalization link, Stage 029’s bibliographic role and S15 item 4’s deferral have textual support. The following changes concern inputs or the physical interpretation a record could carry.

1. **Full 1PN Newtonian normalization and free-mouth target (candidate lines 7–9, 19).** “It integrates over a compact defect in a smooth external field” describes the response half. The same paper first constructs a quasi-static scalar for a localized source and makes a controlled constitutive identification with the Newtonian potential. The monopole *calibrates* that identification and `G` is its measured normalization. Equality of passive and inertial monopoles in the receiving body does not set the source strength of its companion. State the calibration and coherent-translation assumption beside the reduced force, and keep the physical drain-rate/mass map open. Line 9 says “S15 measures the signed momentum entering” for a freely draining mouth. S15 supplies a signed **net crossing rate** on prescribed shell surfaces with unspecified fields and compares expressions, leaving which flow is realized open. Carry its paired-versus-isolated free-mouth acceleration under the same ambient conditions, with direction and magnitude open. Otherwise a simulation author can use a free-mouth measurement containing self or background acceleration, or mistake an unspecified crossing for a measured flow.

2. **Full 1PN matching provenance (candidate line 19).** “Mass dressing ... and the wake construction” permits the full EIH agreement to be recorded without identifying two uses of targets in that paper. Newtonian matching fixes the scalar dressing `κ_ρ=1`. The velocity cross coefficients are imposed as EIH targets on the wake family and the thermodynamic wake closure selects a member. Include both directions of inference in this row, by cross-reference to the bridge where appropriate. Retain the exact assembled EIH equality and the separate `κ_PV` protocol statement.

3. **2.5PN branch hypotheses (candidate line 21).** “A compact, passive, outgoing quadrupole branch” describes the surviving channel. The uniqueness of that channel among the odd terms audited also assumes no new direct scalar source, outgoing compact lower scalar/vector outlets with the stated thresholds, low-frequency rotational invariance, the strict small-body hierarchy and the natural nonzero static orbital overlap. A passive outgoing quadrupole could coexist with an unsuppressed scalar or dipole if those hypotheses fail. Point the author to the paper’s explicit theorem hypotheses and distinguish structural survival and the remaining absolute coefficient. The normalization stays open in the moving-throat PDE.

4. **3PN comparable-mass inputs (candidate line 22).** “They inherit ... grouped quadrupole response” conflates the inherited *choice of channel* with its 3PN constitutive coefficients. Nine comparable-mass interaction coefficients are obtained by inverting the grouped map on the solved GR residual. The geometry term supplies the rest of an imported static target. The row names the one-body fit but omits these comparable-mass matches. Record the direction of the matching alongside the fixed-chart equality and the still-open moving-throat origins.

5. **4PN hereditary transport (candidate line 23).** The paper freezes the standard GR tail functional. Its canonical coefficient bridge and its explicit parameterization of toy transport occupy distinct equations. The latter gives `α_tail^toy = Θ_tail (GM/2c_s³) γ_toy` and becomes the former for `c_s=c`, `γ_toy=γ_quad^eff` and `Θ_tail=1`. The text calls possible mismatch a scalar monopole-scattering transport factor, while the final conditional theorem counts one shared quadrupole normalization. Cite both passages. Keep the single shared normalization in the canonical statement and attach the transport condition to its identification with the simulated toy tail. Recording only a shared physical input can imply that measuring the 2.5PN amplitude alone measures the toy tail.

6. **Type of result (candidate lines 3, 27).** “Which orbital observable to test at each order” asks for an orbital observable even when the cited result is an effective coefficient, a coordinate Hamiltonian or the existence of an ordinary Lagrangian representative in a fixed chart. Qualify the instruction to give each paper’s reported quantity and a direct orbital observable where the paper supplies one. An effective-dynamics comparison should identify its chart/representation. This preserves the cite-only remit and avoids asking the author to derive an unreported gauge-independent observable.

## Source retrievals

Numbered lines are literal command output. Candidate line numbers refer to the reviewed hash below.


```bash
sha256sum research/pde_ledger_v3/directives/S17_record_decision_list.md research/pde_ledger_v3/directives/_measurements/S17_record_decision_list.md
```

```text
f55a0eac7259ea020bff8eaaaaa3d38bb10959f20ec710cd69e626c6cb457b99  research/pde_ledger_v3/directives/S17_record_decision_list.md
4c8c4ba17168ccf6ee63e8cc56e73ae8dc114f1f7a64c63d89e6e06741a1d609  research/pde_ledger_v3/directives/_measurements/S17_record_decision_list.md
```


### Finding 1


```bash
git show 438d270c:research/4d_1pn_full/paper/4d_1pn_full.tex | nl -ba | sed -n '582,586p;673,680p;684,690p;812,820p'
```

```text
   582	It is important not to overstate the claim at this stage. Equation
   583	\eqref{eq:dictionary-longitudinal-identity} is exact. Equation
   584	\eqref{eq:dictionary-poisson-hook} is a regime statement. The identification of the
   585	resulting near-zone scalar with the Newtonian potential used by a reduced point particle
   586	is a further controlled constitutive step carried out later in the paper.
   673	Operationally, we therefore define the measured Newtonian potential \(\PhiN\) as the
   674	properly normalized near-zone scalar channel on the brane,
   675	\begin{equation}
   676	\PhiN(\x,t)\equiv \lambda_\Phi\,\phibr(\x,t),
   677	\label{eq:newtonian-potential-definition}
   678	\end{equation}
   679	with \(\lambda_\Phi\) fixed by monopole calibration. For an isolated compact source \(B\),
   680	the near-zone solution outside the source support is then
   684	-\frac{Gm_B}{|\x-\mathbf X_B(t)|}
   685	+ O\!\left(\frac{Gm_B a_B^2}{|\x-\mathbf X_B|^3}\right),
   686	\label{eq:isolated-source-newtonian-potential}
   687	\end{equation}
   688	where \(G\) is the measured brane normalization of the scalar channel and \(a_B\) is the
   689	source size. The remaining task in this section is to show that a second compact defect
   690	moving in this field reduces, at leading order, to a Newtonian point particle.
   812	At Newtonian order, the second term is an internal energy contribution. It can renormalize
   813	the defect rest energy, but it does not alter the leading external center-of-mass
   814	kinematics so long as the coherent-translation closure is maintained.
   815	
   816	This is the first place where the reduced description already displays an important
   817	feature that later becomes the Newtonian matching condition of the scalar sector:
   818	the same total monopole \(M\) appears in the kinetic term and in the leading coupling to
   819	the external potential below. In the reduced Newtonian description, inertial and
   820	gravitational mass therefore coincide automatically.
```


```bash
git show 438d270c:research/pde_ledger_v3/steps/S15_two_drain_region.md | nl -ba | sed -n '3p;5p;65p;67p'
```

```text
     3	**Result.** S15 supplies a test a two-drain flow must pass **if its signed net shell crossing is to scale as the inverse square of separation**. All three components still contain unspecified fields, even on the power-law profile family. The two engines agree exactly on the general crossing and its geometric departure. Their expressions for the additional field-response departure agree at 30 tested points per component. Which response an evolving pair of drains produces remains open. [Measurements M1–M3](_measurements/S15_record_measurements.md) · [Accepted comparator](../directives/_measurements/S15_reading_comparator_run_record.md)
     5	The Newtonian target is the companion-induced acceleration of a freely moving mouth. A finite shell can store momentum while material crosses its boundary. S15 neither identifies how much reaches the core nor supplies the core's relation between momentum and motion. Whether inverse-square shell crossing is necessary or sufficient for Newtonian acceleration remains unresolved. S21 receives the conditional shell test. Throat acceleration, its direction and its normalization remain open. [Moving-region account](S16_moving_region.md) · [Requirements register](../SUBSTRATE_REQUIREMENTS.md)
    65	The first population route files the user's flow-comparison decision against the existing `R-S21-01`, alongside its distinct light-sector conditions. The comparison stays OPEN. If S21 tests inverse-square signed shell crossing on the stated family, it must carry the geometric departure, live field response, pressure and viscosity inputs, both surfaces and their motion. If it tests Newtonian gravity, it needs companion-induced **free acceleration** relative to an isolated drain under the same ambient conditions. Separate the surroundings from self-field (`R-S16-01`), account for the momentum stored in the shell and exchanged at the sleeve/core and at each slab face (`R-S16-03`), and identify the core contents (S22) and inertial response (`R-S16-04`, S8). The compactness assumption, where used, stays with `R-S16-02`. The physical suppliers of the passive surroundings field, face-resolved momentum and compactness have no owner assigned in S16. [Register entries and S15 pass](../SUBSTRATE_REQUIREMENTS.md) · [S16 record](S16_moving_region.md)
    67	The slow limit and far behaviour remain conditional on `R-S14a-01`, jointly owed by S12, S8/S22 and the unnamed face/interior laws. S22 carries the core, the source map `{J, Θ_Q, N_∞,3, Q}` and the coefficient `I_F,12`. C3 remains the gap between drain rate and mass. The full sign retains quantum, confinement, Maxwell and profile residuals. The forms of `Π_o` and `f` have no assigned owner here. Extending steady P6 and resolving possible stress-carried momentum under P5 remain S8 material questions. A Newtonian orbit, its corrections, brane feedback, local thinning, and mouth support lie beyond this step. [Retained questions](../directives/S15_decision_list.md#retained-open) · [S22 debt rows](../V3_STEP_PLAN.md) · [Defect C3](../DEFECT_REGISTER.md)
```


```bash
git show 438d270c:research/pde_ledger_v3/steps/S16_moving_region.md | nl -ba | sed -n '394,402p'
```

```text
   394	The plan's response-side approximation, supplied mass/multipoles, compactness, smooth
   395	external field, calibrated potential, discarded boundary fluxes and reduced-force equation
   396	are **not printed as an S16 reduction** here. The list implements the later user decisions
   397	with live count, count-moment and boundary-flux objects. The list author's filing leaves
   398	`M`, `a`, multipoles of `M`, compactness, source normalization and any flux suppression
   399	as debts for S22's table.
   400	`a_WT` remains distinct from `a_mouth`; neither control radius nor a count second moment
   401	delivers their A11 bridge. S16 supplies no smooth surroundings field to expand, no
   402	Gaussian profile width and no source normalization. (Retrieve M01/M04/M07.)
```


### Finding 2


```bash
git show 438d270c:research/4d_1pn_full/paper/4d_1pn_full.tex | nl -ba | sed -n '203,204p;2074,2085p;2120,2123p;2146,2154p'
```

```text
   203	    Newtonian scalar coupling & \(\krho=1\) & Newtonian matching of the scalar-dressed particle action (historically \(q=1\)) & \Controlled \\
   204	    EOS exponent & \(n=5\) & weak-field optical/refraction coefficient & \Controlled \\
  2074	Imposing the EIH targets \eqref{eq:wake-eih-target} on
  2075	\eqref{eq:wake-C-parallel-long} gives
  2076	\begin{align}
  2077	-\frac72
  2078	&=
  2079	K_{\rm vec}\pi^2\bigl(-1+a_H^2-\alpha^2\bigr),
  2080	\label{eq:wake-match-para}
  2081	\\[3pt]
  2082	-\frac12
  2083	&=
  2084	K_{\rm vec}\pi^2\bigl(-1+a_H^2+\alpha^2\bigr).
  2085	\label{eq:wake-match-long}
  2120	Thus EIH matching alone does \emph{not} yet pick a unique point. It fixes two invariants
  2121	and leaves one apparent freedom, which may be parameterized by \(a_H^2\). Put
  2122	differently, the cross-tensor structure is already highly constrained by the wake basis,
  2123	but one more physical input is needed to collapse the family.
  2146	The minimal thermodynamic closure used here identifies the longitudinal wake fraction with
  2147	the non-internal part of the barotropic budget:
  2148	\begin{equation}
  2149	\alpha^2_{\rm thermo}
  2150	\equiv
  2151	1-\frac{U}{P}
  2152	=
  2153	1-\frac{1}{n-1}.
  2154	\label{eq:wake-thermo-alpha}
```


### Finding 3


```bash
git show 438d270c:research/4d_2_5pn/paper/4d_2_5pn.tex | nl -ba | sed -n '2456,2475p;2537,2542p;2552,2561p'
```

```text
  2456	The cleanest theorem statement is now conditional rather than unconditional.
  2457	Assume that the completed moving-throat PDE satisfies all of the following on
  2458	its physical low-frequency branch:
  2459	\begin{enumerate}
  2460	  \item there is no new direct scalar source beyond the continuity-mediated and
  2461	  reciprocal-boundary ontology isolated in
  2462	  Sections~\ref{sec:scalar-audit-demotion} and
  2463	  \ref{sec:benchmark-framework},
  2464	  \item the scalar and vector outlets are compact and outgoing at low
  2465	  frequency, so that the surviving scalar and dipole channels begin at the
  2466	  derivative/high-order thresholds established earlier,
  2467	  \item the low-frequency defect/worldtube sector is rotationally invariant,
  2468	  \item the surviving STF quadrupole branch is passive and outgoing,
  2469	  \item the strict small-body hierarchy remains valid,
  2470	  \item and the static orbital/source overlap remains on the natural branch,
  2471	  \[
  2472	  \hat m_0 = 1 + O\!\left(\frac{a^2}{r^2}\right),
  2473	  \qquad
  2474	  \hat m_0\neq 0.
  2475	  \]
  2537	and that the outgoing $\ell=2$ branch has the correct positive
  2538	$+i\omega^5$ parity/sign.  In other words, the surviving universal branch is
  2539	already fixed as a source-content, parity, and sign statement.  What is no
  2540	longer open is whether the model can support \emph{some} quadrupolar odd term.
  2541	What remains open is whether the actual moving-throat PDE realizes the branch
  2542	with the correct normalization.
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


### Finding 4


```bash
git show 438d270c:research/4d_3pn/paper/4d_3pn.tex | nl -ba | sed -n '8331,8338p;8533,8537p;8969,8971p'
```

```text
  8331	\subsection{Exact fit to the solved GR middle block}
  8332	\label{app:richer_grouped_p2_compiler_inverse_map:gr_fit}
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
  8969	In that precise sense, the conservative $3$PN two-body ledger is now fully assigned in the fixed ADM chart. What remains open is no longer algebraic coefficient matching.
  8970	
  8971	The open problems are now dynamical and constitutive. First, the richer grouped real $P_2$ constitutive map should be derived directly from the moving-throat problem rather than by inverse matching. Second, the unique geometry completion should be derived directly from the throat-side geometry/compressibility dynamics. Third, those two conservative lanes must be carried into the already narrowed $2.5$PN bridge, where the live universal route has already been reduced to the grouped real quadrupole bundle and its passive/outgoing normalization problem.
```


### Finding 5


```bash
git show 438d270c:research/4d_4pn/paper/4d_4pn.tex | nl -ba | sed -n '2169,2175p;2368,2397p;2406,2417p'
```

```text
  2169	side is therefore sharply narrower than the local side ever was.  The purpose of
  2170	this section is to make that narrowing explicit.  First, we freeze the standard
  2171	conservative GR tail functional and its universal local logarithmic shadow.
  2172	Second, we show that the hereditary coefficient is tied exactly to the same
  2173	canonically normalized STF quadrupole coefficient that already controls the
  2174	$2.5$PN Burke--Thorne channel.  Third, we isolate the only remaining model-side
  2175	ambiguity as a single scalar transport factor rather than a new tensor channel,
  2368	be the canonically normalized STF quadrupole coefficient inherited from the
  2369	$2.5$PN audit.  Then the hereditary side may be parameterized by one scalar
  2370	transport factor $\Theta_{\rm tail}$ through
  2371	\begin{equation}
  2372	\boxed{
  2373	\alpha_{\rm tail}^{\rm toy}
  2374	=
  2375	\Theta_{\rm tail}
  2376	\frac{GM}{2c_s^3}
  2377	\gamma_{\rm toy}.
  2378	}
  2379	\label{eq:4pn-tail-toy-transport}
  2380	\end{equation}
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
  2416	  \label{eq:4pn-tail-gr-recovery-condition}
  2417	  \end{equation}
```


### Finding 6


```bash
git show 438d270c:research/4d_2pn/paper/4d_2pn.tex | nl -ba | sed -n '110,115p'
```

```text
   110	\(\lambda_\rho=\tfrac12\).  The remaining comparable-mass mismatch against the
   111	generic-frame \(\ADM\) target collapses to a finite algebraic residual with a
   112	unique solution for the added cross-sector coefficients \(\{q_i,t_i,s_1\}\).
   113	After assembly, the full reduced two-body Lagrangian Legendre-transforms
   114	exactly to the standard generic-frame conservative \(\ADM\) Hamiltonian through
   115	\(\TwoPN\).  The paper therefore gives a full conservative \(\TwoPN\)
```


```bash
git show 438d270c:research/4d_4pn/paper/4d_4pn.tex | nl -ba | sed -n '3254,3258p'
```

```text
  3254	Together they verify the full Hamiltonian-first local theorem chain: existence
  3255	of a generic-frame Hamiltonian representative, exact ordinary sign flip,
  3256	structured liftability of the aligned-seed correction, and final reconstruction
  3257	of one exact generic-frame ordinary representative whose COM reduction
  3258	reproduces the complete fixed-chart local ordinary target.
```


```bash
git show 438d270c:research/pde_ledger_v3/V3_STEP_PLAN.md | nl -ba | sed -n '747,749p'
```

```text
   747	### S17 · The PN ladder — cite-only
   748	Seven separately-published papers, Zenodo-DOI'd (`ledger_stage029_pn_corpus_doi_cite.md:32-40`).
   749	⛔ **Do not re-derive.** Record what each establishes and its inputs.
```
