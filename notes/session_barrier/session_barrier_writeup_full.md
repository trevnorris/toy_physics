# Session barrier thread write-up — Sections 0 through 8

## 0. Scope, status labels, and what this write-up is actually about

This write-up covers the **barrier-crossing / stability / damping / materials** thread developed in this session. It does **not** try to retell the whole toy-model program, and it does not treat the session’s reduced barrier results as if they were already a theorem of the fully solved moving-throat PDE. The point here is narrower and more useful: to reconstruct what framework was carried in, what assumptions were relaxed, what the reduced computations actually established, and where the session still sat inside the project’s own claim-status hierarchy.

The status language should stay the same as the compact moving-throat program. At the bottom are **exact** parent-theory statements: the action, Euler–Lagrange equations, the corrected charge ontology, the projection identities, and the localized Maxwell equations. Above that are **exact within closure** statements, where the algebra is exact once a reduced family or branch is declared. Then come **reduced / controlled reductions**, where one explicitly assumes zero-mode dominance, weak gradients, quasi-staticity, or a low-frequency wall/support truncation. The session’s barrier calculations also used **effective closures**, meaning deliberate phenomenological choices added to stress-test a branch that the completed PDE has not yet uniquely fixed. Finally, some objects are still **open**, because they depend on the actual completed moving-throat branch rather than on already-frozen algebra.

That classification matters because these first three sections mix all of those layers. The **parent 4+1 theory** and the **projection / Maxwell identities** are part of the exact backbone. The moving-throat wall/support/outgoing stack is presently strong enough to organize reduced theorem work inside its carried closures, but the compact program still states plainly that the main remaining gap is **branch realization in the completed PDE**, not more algebraic packaging. The barrier session therefore lives mostly in the zone of **controlled reduction plus effective closure**, built on top of an exact action-level base.

Within that framing, this session actually did five things. It first relaxed three stationary constraints and recomputed a reduced same-charge effective potential. It then promoted that reduced potential into a simple time-evolution problem, added the magnetic/helicity sector, and extracted a turning point, WKB multiplier, and a preferred spin orientation. It next introduced a proton-proxy inertia scale and scanned for a crossing-versus-collapse window. It then damped the dressing leg to test whether an external lattice heatsink could defeat the earlier “speed trap.” Finally, it mapped the resulting phenomenological damping / gradient thresholds onto condensed-matter quantities such as
\(\lambda_{\rm ep}\omega_D\), an effective lattice stiffness, and the Korringa-limited operating temperature. Those results were recorded in the session artifacts rather than in the original source stack. [Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt), [dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt), [proton-proxy stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt), [damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt), and [material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt).

Just as important is what this session did **not** do. It did not solve the full nonlinear moving-throat PDE. It did not prove an unconditional same-charge fusion theorem. It did not derive spin, chirality, or electroweak structure from the defect ontology. It also did not prove that any specific reactor matrix—Palladium, Titanium, or anything else—already satisfies the derived threshold conditions. The compact program still treats the moving-throat branch, the full outgoing normalization, and the remaining endgame packet as the live theorem bottlenecks, and the earlier particle-target handoff likewise treated spin and magnetic-moment questions as requirement-audit territory rather than completed derivations.

Because of that, the right way to read the session is to separate **derived reduced-model statements** from **higher-level interpretation**. The reduced-model statements are things like: activating leakage produces explicit \(S_{\rm leak}\) and \(J^wE_w\) channels; a dynamic \(U/V\) coupling can lower the reduced \(V_{\rm eff}(r)\); aligned spins export more subscale helicity than anti-aligned ones in the chosen closure; a damped dressing leg opens a cold-survival window; and the resulting thresholds can be written as material-design equations. By contrast, claims such as “the reactor is proved to run at room temperature” or “the lattice naturally guarantees safe cold fusion” are one interpretive step higher than the session actually established. The barrier audit is especially useful here, because it had already narrowed the static same-charge mixed corridor to **short-range coefficient engineering**, not a new long-range attractive law, and had already shown that the first outgoing linear dynamic correction is **phase lag / pumping**, not conservative barrier lowering.

## 1. Source stack and baseline model

The starting point for everything in this session was the action-based **4+1-dimensional bulk theory** summarized in the corrected `4d` paper stack. In that baseline, the matter sector is a gauged 4D-spatial GNLS medium with a frozen stiff-polytropic equation of state, the electromagnetic sector is a localized 4+1 Maxwell theory with transverse profile \(Z(w)\), and the geometry sector is represented at the parent-paper level by effective throat variables \(a(t)\) and \(L(t)\). The same core paper is also where the project’s corrected charge ontology was frozen: electric sign is the topological puncture orientation \(\eta_Q=\pm1\), the microscopic branch charge is \(q_*=\eta_Q e_*\), and the observable brane charge strength is controlled by localization thickness rather than by circulation. The mixed core sector is suppressed only in the far-field zero-mode brane limit, not removed from the microscopic ontology.

A second foundational ingredient is the project’s insistence that a 3D brane observer sees an **operational projection**, not a hard boundary condition. With a projection weight \(W(w)\), the brane subsystem is generically **open**: projected continuity acquires the leakage source \(S_{\rm leak}\), and the projected longitudinal velocity obeys an exact identity that becomes a Poisson-like equation only in a controlled quasi-static regime. This was the gravity-side route by which inverse-square longitudinal behavior first appeared on the brane, and it is also why later barrier work can legitimately talk about energy or density leaving the resolved brane subsystem without violating bulk conservation.

The EM papers added the other half of the baseline. Starting from the localized 4+1 Maxwell action, they derive the exact bulk field equations and then show that ordinary 3+1 Maxwell on the brane is recovered only under explicit reduction assumptions such as zero-mode dominance and suppression of the mixed sector. In that controlled limit, \(\mu_0\) is replaced by \(\mu_0^{\rm eff}=\mu_0/Z_{\rm int}\), while canonical zero-mode normalization gives \(q_{\rm eff}=q_*/\sqrt{Z_{\rm int}}\) and \(e_{\rm eff}=e_*/\sqrt{Z_{\rm int}}\). For a Gaussian localization profile, the transverse Sturm–Liouville problem closes into a Hermite / KK tower with masses \(m_n^2=2n/\lambda^2\), and the brane-to-brane static field becomes Coulomb plus Yukawa corrections. Odd modes decouple from a centered brane source at leading order by parity, but they are not thereby deleted from the underlying ontology.

The plasma extension is the place where the “hidden channels” that mattered for this session become explicit. It reformulates the same parent ontology as a 4+1 bulk plasma observed by a 3D brane subsystem, and it emphasizes that relaxing the standard reduction assumptions opens four beyond-MHD channels: genuine mixed-sector electromagnetism through \(A_w\) and \(F_{\mu w}\), explicit brane–bulk exchange through \(j^w\), finite-localization transverse towers, and projected unresolved-stress channels. In that language, \(E_w=-\partial_tA_w-\partial_wA_0\) and \(C_a=F_{aw}\) are real observables, not artifacts, and the bulk EM energy ledger contains the explicit work term \(J^wE_w\) together with explicit mixed Poynting transport along \(w\). Those facts are exactly what later made it legitimate to monitor leakage, scalar-photon work, and energy export into unresolved or transverse structure when the barrier constraints were relaxed.

The moving-throat program then lifts the old finite-dimensional geometry closure into a distributed wall/support/outgoing language. In its compact form, the program now treats the parent 4+1 field theory, the geometry lift, the reduced wall/BdG/Maxwell/mixed hierarchy, and the grouped real \(P_2\) finish line as a stable conceptual backbone, while still marking the actual completed moving-throat branch as open. The non-negotiable notation firewall is repeated there: electric charge is carried by \(\eta_Q,q_*,q_{\rm eff}\), not by circulation; the historical gravity-side bare \(q=1\) is really \(\kappa_\rho=1\); grouped labels \(20/21/22\) are grouped real \(P_2\) lanes; and the mixed channels \(A_w,J^w,F_{\mu w},E_w,C_a\) remain microscopic degrees of freedom outside the strict far-field brane limit. This compact program is the cleanest source for understanding what the later barrier scripts were allowed to relax and what they were not.

Finally, before the session’s own computations began, the same-charge **barrier audit** had already sharply constrained what a mixed-sector corridor could mean. At the static one-port level, it found that the mixed bundle does **not** generate a qualitatively new long-range attractive family. It only renormalizes short-range structures, giving the bundle access to families such as \(x^{-6}\), \(e^{-2\kappa x}/x^4\), and \(e^{-4\kappa x}/x^2\). The audit also showed that the first linear outgoing dynamic correction is not a conservative barrier-softening term: it enters as phase lag / pumping / leakage. In other words, by the time this session started, the remaining same-charge corridor was already narrow. Any survival route would have to work through short-range coefficient engineering, energy export into transverse structure, or more genuinely dynamical effects—not by discovering a new long-range attraction hidden in the mixed sector.

## 2. Session I — Stationary relaxed-constraint barrier test

The first session block asked a very specific question: **what happens if the previous “standard-physics recovery” constraints are relaxed rather than imposed by hand?** The reason this was a sensible stress test is built into the source stack itself. The EM and plasma papers repeatedly state that suppressions such as \(A_w\approx0\), \(F_{\mu w}\approx0\), and especially \(J^w\approx0\) are part of a **controlled far-field brane reduction**, not a microscopic theorem about the core. So the session deliberately broke three assumptions that had previously kept the reduced same-charge problem close to standard brane physics: transverse current suppression, rigid-mouth transfer/dressing factorization, and positivity of the normalized axial source profile.

The first relaxation was to **activate transverse leakage** by allowing \(J^w\neq0\). Instead of forcing the reduced current to live only in the brane directions, the session let the transverse electric field \(E_w\) drive a nonzero \(v^w\), and therefore a nonzero bulk current component \(j^w=\rho v^w\). Conceptually, that move was not imported from nowhere; it was exactly the beyond-MHD scalar-photon channel the plasma paper had already identified, and it immediately turned \(S_{\rm leak}\) and \(J^wE_w\) into observables rather than vanishing identities. In the stationary reduced closure used by the script, the leakage source and scalar-photon work channel were computed explicitly and then monitored as the same-charge separation decreased. The session report shows that they stay finite and grow toward the short-range region, matching the interpretation that bulk-conserved energy can be exported out of the resolved 3D subsystem as the charges are pushed together. [Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt).

In that run, the leakage source and work channel were not left schematic. For the chosen Gaussian-style closure, the report gives
\[
S_{\rm leak}^{(s)}
=
\frac{\sqrt{2}\,E_0\mu_w\rho_0}{2\sqrt{\pi}\,\lambda^3},
\qquad
J^wE_w
=
\frac{2E_0^2\mu_w q\rho_0}{\lambda^2}.
\]
At the representative evaluation point \(r\approx\lambda\), these became approximately \(S_{\rm leak}\approx0.0366\) and \(J^wE_w\approx0.0211\), and the same report notes that both channels strengthen further as \(r\) is decreased toward the contact region. The significance was not the specific numbers by themselves, but the fact that the lowered barrier branch could now export energy into the bulk through an explicitly conserved work channel rather than pretending the brane subsystem was closed. [Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt).

The second relaxation was to **break the rigid-mouth assumption**. In the barrier-audit language this meant allowing the transfer-shape variable and the dressing variable to stop factorizing cleanly when the gradient became steep. The session introduced the logarithmic variables \(U\) and \(V\), then added a \(\chi_\lambda\)-scaled cross coupling so that the mouth geometry was allowed to flex as same-charge gradients sharpened. In the chosen reduced closure, this immediately destroyed the invariance of the dressing leg: \(V\) became nonzero, \(\epsilon_\eta\) stopped behaving as a frozen bystander, and the first-order dynamic barrier scalar \(\Xi_1\) became a directly trackable function on the relaxed branch. At the same representative point \(r\approx\lambda\), the report records \(\Xi_1\approx0.1431\), \(V\approx-0.0362\), and \(\epsilon_\eta\approx0.2893\), together with a positive \(U/V\) drain term that quantified transfer of energy into the dressing leg. In plain terms, the rigid-mouth closure had hidden an active channel: once relaxed, the dressing geometry began to absorb energy from the repulsive side of the interaction rather than simply watching it. [Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt).

The third relaxation was to allow **sign-changing multimode mouth/core sources**. Earlier positive-source Family-1 closures assumed a nonnegative normalized axial source profile \(\sigma(z)\). The session instead promoted \(\sigma\) to a small multimode superposition,
\[
\sigma(x)=1+a\cos(\pi x)+b\cos(2\pi x),
\]
and then tracked the corresponding mouth-bias functional \(\mathfrak g[\sigma]\) and mixed-to-shell loading ratio \(\mathcal R\). In the branch actually sampled by the stationary script, the report found an explicit sign change, with \(\sigma_{\min}<0\) and `sign_change = True`. That mattered because it meant the run was no longer living inside the positive-source Family-1 corridor at all; it was probing a compensated branch in which different source lobes could oppose one another while still altering the mouth bias and short-range susceptibility. Again, this was not presented as a theorem of the full mouth PDE. It was a controlled stress test of a closure that the prior audit had not allowed. [Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt).

Once those three relaxations were turned on together, the session recomputed the reduced same-charge effective potential \(V_{\rm eff}(r)\). The key result was that the potential **did** drop substantially below the naive Coulomb branch in the near-contact region, but in a way that remained consistent with the barrier audit’s earlier structural restrictions. The report’s strongest sampled softening occurred at \(r=0.18\), where
\[
\frac{V_{\rm eff}}{V_{\rm Coul}}\approx 0.31446203.
\]
So the barrier seen by the reduced same-charge coordinate was cut to about thirty-one percent of the pure Coulomb value in that window. But that lowering was not credited to a newly discovered long-range attraction. The report itself explicitly interprets the reduction as the combined effect of short-range families, transverse leakage, scalar-photon work, dressing-leg drain, and compensated mouth response—exactly the sort of short-range/open-system reshaping that the barrier audit still allowed. [Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt).

That final point is important enough to state plainly. Session I did **not** overturn the barrier audit. It did not show that the mixed sector secretly generates an attractive \(1/r\)-like law. Instead, it demonstrated that once the far-field suppression \(J^w\approx0\), the rigid-mouth factorization, and the positive-source assumption are relaxed together, the reduced same-charge branch can lower its near-contact effective barrier by exporting energy into the bulk and into internal dressing / mouth structure, while still respecting the audit’s basic message that the admissible correction families are short-ranged and tightly constrained by the same one-port bundle that later feeds the normalization chain.

The Session I artifacts are:
[SymPy stress-test script](sandbox:/mnt/data/stress_test_relaxed_constraints_sympy.py)
[Stationary relaxed-constraint report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt)
[Effective-potential plot](sandbox:/mnt/data/relaxed_constraints_veff.png)
[Leakage / Xi₁ / mouth-response diagnostics](sandbox:/mnt/data/relaxed_constraints_diagnostics.png)

## 3. Session II — Dynamic scattering, turning point, and helicity export

Session II took the stationary lowered-barrier branch from Session I and asked the next obvious question: does the reduction merely soften a plotted potential, or does it open a genuine **dynamical** corridor in which an approaching same-charge particle can get farther than the pure Coulomb model would allow? The answer was sought without changing the basic structural discipline. The run kept the same reduced short-range family structure already allowed by the barrier audit, rather than introducing a new long-range interaction by hand. That mattered because the barrier audit had already shown that the mixed bundle’s first linear dynamic correction is phase lag / pumping rather than a new conservative attraction, so any improvement in approach had to come from the lowered near-contact branch and the new open-system channels, not from inventing a different asymptotic law.

The dynamic script therefore promoted the reduced effective potential into a simple one-dimensional scattering problem for the separation coordinate,
\[
m_s\ddot r=-\partial_r V_{\rm eff}(r),
\]
and then compared that motion directly against the pure Coulomb reference using the same contact radius. On the lowered branch, the report found a barrier peak at
\[
r_{\rm peak}\approx 0.23944389,
\qquad
V_{\rm peak}\approx 3.42933112,
\]
with a finite reduced classical threshold speed
\[
v_{\rm crit,new}\approx 2.54139063.
\]
For the same contact radius, the Coulomb comparison needed about
\[
v_{\rm contact,Coul}\approx 3.27278339.
\]
That is the first clean dynamic result: the reduced branch creates a real window
\[
v_{\rm crit,new}<v_0<v_{\rm contact,Coul}
\]
in which the lowered same-charge barrier is classically traversable while the pure Coulomb model still cannot reach the same contact scale. [Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)

The subbarrier test then asked a softer question: even before classical crossing, does the lowered branch improve under-barrier penetration? On the default slice
\[
E_{\rm sub}=2.5,
\]
the new-model outer turning point moved slightly inward to
\[
r_{\rm turn,new}\approx 0.39096144,
\]
versus
\[
r_{\rm turn,Coul}\approx 0.40000141.
\]
More importantly, the WKB action fell from about
\[
I_{\rm Coul}\approx 0.30222297
\]
to
\[
I_{\rm new}\approx 0.19744614,
\]
which raised the reduced tunneling factor by
\[
\frac{T_{\rm new}}{T_{\rm Coul}}\approx 1.23312756.
\]
In the report’s own wording, that corresponds to a fusion-probability increase of about **23.31%** on that reduced slice. The same run also tracked the dynamic barrier scalar at the new turning point and found
\[
\Xi_1(r_{\rm turn})\approx 0.34437471.
\]
So the session did not stop at “the potential looks lower”; it carried that lowering into an explicit reduced scattering exponent. [Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)

The most concrete classical demonstration came just above threshold. At
\[
v_0\approx 2.59221845,
\]
the lowered-branch trajectory reached the chosen contact radius, while the pure Coulomb comparison still turned back at about
\[
r\approx 0.28091705.
\]
This is worth stating plainly because it answers the motivating dynamic question in the narrowest possible way: in the reduced closure actually used here, there exists a speed range where the softened same-charge branch reaches contact and the Coulomb branch does not. That is a stronger statement than the original stationary ratio plot, even though it still lives inside a reduced model rather than the full moving-throat PDE. [Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)

The same dynamic continuation also opened the magnetic / vortical sector. That move was not ad hoc. The parent theory already carries the exact minimal-coupling vorticity identity
\[
\Omega_{ij}= -\frac{q_s}{m_s}F_{ij},
\]
and the plasma extension already insists that the mixed fields \(A_w\), \(F_{\mu w}\), and \(J^w\) remain microscopic degrees of freedom outside the strict far-field brane reduction. So adding the Lorentz-force contribution and asking whether unresolved helicity is exported into the higher / transverse sector is consistent with the existing ontology, not a separate theory bolt-on.

Within that closure, two initial spin states were compared: aligned and anti-aligned. The outcome was not that one branch could cross and the other categorically could not. In the representative run, **both** branches reached contact. The difference lay in how efficiently each one exported repulsive magnetic or topological structure into the unresolved sector. The aligned case reached a peak helicity-export rate
\[
\max(\partial_t h_{\rm sub})\approx 281.7983,
\]
while the anti-aligned case reached only
\[
\max(\partial_t h_{\rm sub})\approx 56.9688.
\]
The peak export ratio was therefore about
\[
4.95,
\]
and the integrated helicity-export ratio was still about
\[
4.11.
\]
By the end of the run, the aligned branch carried
\[
h_{\rm sub,final}\approx 20.5807,
\]
compared to only
\[
5.0084
\]
for the anti-aligned branch. So the correct reading is not “anti-aligned fails absolutely.” It is that, on this reduced mixed-sector closure, **aligned spins are the branch that most effectively unload unresolved repulsive structure into the higher transverse sector**, while anti-aligned spins partially self-cancel that export channel. [Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)

The last task in Session II translated the reduced turning-point data into the gradient trigger
\[
\chi_\lambda\equiv \lambda\,\bigl|\partial_r\ln V_{\rm eff}(r)\bigr|.
\]
Using the operational definition adopted in the run,
\[
\lambda_{\rm th}(r_{\rm turn})
=
\left|\frac{V_{\rm eff}(r_{\rm turn})}{V'_{\rm eff}(r_{\rm turn})}\right|,
\]
the representative turning point gave
\[
\lambda_{\rm th}\approx 0.42826825.
\]
Scanning the turning-point branch produced a threshold-width range of roughly
\[
0.40673709\le \lambda_{\rm th}\le 1.06949146.
\]
This mattered because the session now had a direct way to say how steep the confinement had to be before the beyond-MHD / transverse-bypass closure even turned on. What it still did **not** have was a fully fixed physical unit map from reduced \(r\) and \(\lambda\) to a specific laboratory material. That later mapping step only became meaningful because Session II had already extracted this explicit trigger variable. [Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)

So Session II changed the status of the barrier thread in three ways. It demonstrated a reduced classical-crossing window that the pure Coulomb comparison did not share, it quantified a subbarrier enhancement on the same branch, and it identified aligned-spin helicity export as the preferred magnetic route into the pure-transfer subcorridor. But none of that overrode the broader claim-status firewall. The compact moving-throat program still treats branch realization and the final outgoing normalization as open, so the dynamic scattering and helicity results belong to the same “reduced / controlled reduction plus effective closure” category as the stationary barrier lowering that came before them.

The Session II artifacts are:
[Dynamic SymPy script](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_sympy.py)
[Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)
[Potential plot](sandbox:/mnt/data/dynamic_scattering_potential.png)
[Trajectory plot](sandbox:/mnt/data/dynamic_scattering_trajectories.png)
[Helicity-export plot](sandbox:/mnt/data/dynamic_helicity_export.png)
[Lambda-threshold curve](sandbox:/mnt/data/lambda_threshold_curve.png)

## 4. Session III — Structural stability and the Goldilocks window

Once Session II showed that the lowered branch could be crossed dynamically, the next question was no longer “can a particle get through?” but “can it get through **without its own dressing geometry failing first?**” That is a different question, because the earlier relaxed rigid-mouth closure had already shown that energy can be drained into the dressing leg \(V\). Session III therefore reframed the problem as a competition between two characteristic timescales: the barrier-region transit time and the dressing-leg collapse time. In spirit that move was natural, because the moving-throat lift already treats geometry through effective wall inertias and stiffnesses, and the distributed-wall reduction explicitly identifies \(\mu_\eta\) as the wall inertia density feeding the low-mode geometry equations. What remained closure-dependent was the exact choice of reduced “collapse” metric used for the stress test.

The script therefore introduced the characteristic estimate
\[
t_{\rm cross}(E)=\lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}},
\qquad
t_{\rm collapse}=\sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}},
\]
and studied the stability ratio
\[
\mathcal S(E)=\frac{t_{\rm cross}}{t_{\rm collapse}}.
\]
The simple but important algebraic observation built into the report is that if the heavy-throat proxy scales mass and wall inertia together,
\[
\mu_\eta=\alpha m_s,
\]
then
\[
\mathcal S(E)
=
\frac{\sqrt{2}\,\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}}}{2\sqrt{\alpha}}\,(E-V_{\rm peak})^{-1/2},
\]
and the corresponding lower edge becomes
\[
E_{\rm edge}=V_{\rm peak}+\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
\]
So the heavy-throat scaling does **not** move the edge by itself when \(m_s\) and \(\mu_\eta\) are scaled together; what matters is the confinement-width / steepness choice entering \(\chi_\lambda\). [Proton-proxy stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt)

The specific run then imposed the requested proton-proxy scaling
\[
m_s=\mu_\eta=1836.15267343,
\]
and, for its main branch, kept the previously derived trigger width
\[
\lambda_{\rm eff}=\lambda_{\rm th}(r_{\rm turn})\approx 0.42826825
\]
as the active barrier width. On that branch, the steepest logarithmic gradient was recorded at the contact-side edge and gave
\[
\chi_\lambda^{\rm peak}\approx 21.73204372,
\]
which in turn produced
\[
t_{\rm collapse}\approx 9.43066476.
\]
The analytic lower edge of the stable region then came out to
\[
E_{\rm safe,min}\approx 5.32265943,
\]
with corresponding speed
\[
v_{\rm safe,min}\approx 0.07469791.
\]
Measured in units of the proton-proxy classical threshold speed, that is about
\[
1.25948\,v_{\rm crit,p}.
\]
The scan found no upper collapse edge before the top of the sampled band, so the safe region on this closure was **one-sided** rather than a closed island:
\[
5.32265943\lesssim E_{\rm inc}\lesssim 80.93332737.
\]
[Proton-proxy stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt)

The aligned-spin dynamic cross-check then made the conclusion less abstract. Across the sampled over-barrier aligned trajectories, every run reached contact, and the actual barrier-region transit times lay between about
\[
0.204
\quad\text{and}\quad 4.054,
\]
which remained below the characteristic collapse time throughout the sampled band. So the analytic stability ratio was not flattering the branch; if anything, it was the more conservative criterion. That is why Session III could honestly claim a “Goldilocks zone”: not because every conceivable trajectory survives, but because on the chosen reduced branch there is a definable lower incident-energy edge above which crossing outruns collapse. [Proton-proxy stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt)

The section’s most important caveat came from the width sensitivity test. If one threw away the trigger-width interpretation and instead used the raw model width
\[
\lambda=1,
\]
then the same reduced branch gave
\[
\chi_\lambda^{\rm peak}\approx 50.74399964,
\qquad
t_{\rm collapse}\approx 6.17163516,
\qquad
E_{\rm safe,min}\approx 27.53273095.
\]
So the location of the safe region was far more sensitive to the geometric steepness choice than to the proton mass factor itself. In other words, Session III did not just locate a stability window; it also identified what the window was *most* sensitive to. The dominant lever was the confinement-width / gradient choice inherited from Session II, not simply the substitution of a heavier particle mass. [Proton-proxy stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt)

That sensitivity result is also what keeps the interpretation honest. Session III did **not** show that a proton-proxy defect is generically stable whenever it is heavy. It showed something narrower: once the same reduced branch is equipped with the earlier trigger width and the chosen aligned-spin closure, there is a one-sided over-barrier regime in which transit beats collapse, and the survival edge is controlled primarily by \(\chi_\lambda\) and \(\lambda_{\rm eff}\). The moving-throat hierarchy had already prepared that conclusion by making wall inertia and geometry response explicit; the session’s added value was to turn that structure into a concrete timescale inequality and then sweep it numerically.

The Session III artifacts are:
[Proton-proxy stability script](sandbox:/mnt/data/proton_proxy_stability_sweep_sympy.py)
[Proton-proxy stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt)
[Timescale plot](sandbox:/mnt/data/proton_proxy_stability_timescales.png)

## 5. Session IV — Damped dressing leg and lattice heatsink

Session III still left an obvious weakness: the safest region it found lived above the characteristic collapse edge, while a slower “cold” crossing remained unsafe if the dressing leg was treated as a perfect energy trap. Session IV addressed exactly that omission. Rather than pretending the \(V\)-leg could only store incoming drain energy, it added a Langevin-style shedding term
\[
\gamma_{\rm tot}\,\dot V,
\qquad
\gamma_{\rm tot}=\gamma_{\rm vac}+\gamma_{\rm lattice},
\]
and asked whether sufficiently rapid shedding could stabilize a slow aligned-spin event. This move was explicitly phenomenological: the source stack does allow effective damping completions in the geometry sector, but it does not yet derive this specific dressing-leg dissipation law from the completed PDE. So Session IV should be read as an open-system **effective closure** layered on top of the already reduced barrier branch, not as an exact theorem of the parent action.

Within that envelope closure, the report used the characteristic formulas
\[
\gamma_{\rm crit}=\sqrt{\frac{g_{UV}\chi_{\rm peak}}{\mu_\eta}},
\qquad
t_{\rm collapse}^{\rm damped}=
\begin{cases}
(\gamma_{\rm crit}-\gamma_{\rm tot})^{-1}, & \gamma_{\rm tot}<\gamma_{\rm crit},\\[4pt]
\infty, & \gamma_{\rm tot}\ge \gamma_{\rm crit},
\end{cases}
\]
and for a specific cold event defined
\[
\gamma_{\rm safe}(v_0)=\gamma_{\rm crit}-\frac{1}{t_{\rm cross}(v_0)}.
\]
So there are really two thresholds in this reduced model: one that is sufficient for a **particular** event to survive, and a larger one that makes the characteristic collapse time diverge in the envelope description. [Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt)

The cold run itself was fixed at
\[
v_0=2.6
\]
on the aligned-spin branch. The lowered barrier still reached contact, while the pure Coulomb comparison still turned back at about
\[
r\approx 0.27933174.
\]
The actual lowered-branch time to contact was about
\[
2.1129,
\]
while the characteristic crossing-time estimate, using the raw model width \(\lambda=1\) for the requested timescale law, was
\[
t_{\rm cross}\approx 1.82169718.
\]
The undamped collapse time on that same branch was only
\[
t_{\rm collapse,0}\approx 0.14402764,
\]
which gave an undamped stability ratio
\[
\mathcal S(0)\approx 12.64824697.
\]
So in the absence of shedding, the cold event was very much **not** safe: the dressing leg collapsed on the characteristic estimate long before the crossing could be treated as controlled. [Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt)

The central result of Session IV is that a finite shedding window does exist once the V-leg is allowed to dissipate. The report found
\[
\gamma_{\rm crit}\approx 6.94311167,
\]
for unconditional stability in the envelope closure, and
\[
\gamma_{\rm safe}(v_0=2.6)\approx 6.39417302,
\]
for the weaker requirement that this *particular* cold event satisfy
\[
\mathcal S<1.
\]
This is exactly the point where the earlier “cold crossing fails” story changed. The session did not claim that vacuum shedding alone would do the job. It showed instead that once total shedding is strong enough, the slow event can survive on the same lowered barrier branch that had previously looked like a trap. [Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt)

The run then split the dissipated dressing-leg energy between vacuum and lattice channels using a default 3:1 lattice-to-vacuum partition,
\[
\gamma_{\rm vac}:\gamma_{\rm lattice}=1:3.
\]
At the event-survival threshold, that produced
\[
\gamma_{\rm vac}\approx 1.59854325,
\qquad
\gamma_{\rm lattice}\approx 4.79562976.
\]
The total dissipated dressing-leg energy along the actual \(v_0=2.6\) crossing was then
\[
E_{\rm diss,total}\approx 0.01033460,
\]
of which
\[
E_{\rm vac}\approx 0.00258365,
\qquad
E_{\rm lattice}\approx 0.00775095.
\]
At the unconditional-stability threshold the lattice channel still carried essentially the same order of magnitude,
\[
E_{\rm lattice}\approx 0.00770830.
\]
So on this closure, most of the shed dressing-leg energy ends up as **lattice heat**, not as the vacuum-side component. That is the numerical fact that made the later condensed-matter mapping meaningful in the first place. [Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt)

What Session IV therefore changed was not the parent theory, but the practical interpretation of the earlier cold-crossing failure. In the undamped picture, the V-leg behaved like an energy trap and the slow aligned-spin event died on the dressing side before it could count as a controlled crossing. In the damped picture, the question becomes quantitative rather than absolute: how large must \(\gamma_{\rm lattice}\) be, relative to the weaker vacuum channel, for the cold event to survive, and how much energy is dumped into the lattice per event once that happens? The report answered both questions in explicit reduced units. That still falls short of a materials theorem, but it is already much more informative than the earlier undamped “speed trap” verdict. It converts a failure mode into a solvable design inequality.

The Session IV artifacts are:
[Damped-shedding SymPy script](sandbox:/mnt/data/damped_shedding_cold_sweep_sympy.py)
[Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt)
[Timescale plot](sandbox:/mnt/data/damped_shedding_timescales.png)
[Heat-partition plot](sandbox:/mnt/data/damped_shedding_heat_partition.png)

## 6. Session V — Condensed-matter/material mapping

By the time Session V began, the barrier thread had already reached a very specific engineering-style question. The reduced barrier branch no longer failed because “same charges can never meet” inside the chosen closure. It failed, or survived, according to explicit threshold inequalities: a confinement steepness condition from the \(\chi_\lambda\) trigger, a shedding condition from the damped dressing-leg audit, and a spin-survival condition from the aligned-helicity branch. Session V therefore did not introduce a new dynamical corridor of its own. Its job was to translate those already-derived reduced quantities into standard condensed-matter language so the project could eventually ask whether any real lattice even sits in the right regime. The resulting section is best read as a **parameter map**, not as a materials verdict. It turned phenomenological reduced variables into experimentally recognizable targets without yet proving that Palladium, Titanium, or any other specific host actually hits them. [Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

The first mapping was the lattice-shedding requirement. Session IV had already shown that the cold aligned-spin event at \(v_0=2.6\) becomes survivable only if the total dressing-leg shedding exceeds
\[
\gamma_{\rm safe}\approx 6.39417302,
\]
with the default 3:1 split giving a required lattice contribution
\[
\gamma_{\rm lattice}\approx 4.79562976.
\]
Session V then rewrote that demand as a standard electron-phonon turnover condition,
\[
\gamma_{\rm lattice}^{\rm phys}=\zeta_{\rm ep}\,\lambda_{\rm ep}\,\omega_D,
\qquad
\gamma_{\rm lattice}^{\rm phys}=\frac{\gamma_{\rm lattice}^{\rm red}}{t_*},
\]
which yields the exact threshold
\[
(\lambda_{\rm ep}\omega_D)_{\min}
=
\frac{4.79562976}{\zeta_{\rm ep} t_*}.
\]
Using the already-carried reduced crossing time \(t_{\rm cross}^{\rm red}\approx 1.82169718\), the same condition becomes
\[
(\lambda_{\rm ep}\omega_D)_{\min}
=
\frac{8.73618521}{\zeta_{\rm ep}\,t_{\rm cross}^{\rm phys}},
\qquad
(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys})_{\min}
=
\frac{8.73618521}{\zeta_{\rm ep}}.
\]
The report’s plain-language interpretation is cautious and useful: for the normalization \(\zeta_{\rm ep}=1\), the lattice drain has to turn over about **8.74 times during one physical crossing event**. That is the clean database-comparison target produced by the session. It is not yet a statement that any real metal already satisfies it. [Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

The second mapping revisited the geometric trigger \(\chi_\lambda\). Session II had already found that the representative lowered-barrier turning point sits at
\[
r_{\rm turn}\approx 0.39096144,
\qquad
\lambda_{\rm th}(r_{\rm turn})\approx 0.42826825.
\]
Session V asked what sort of interstitial trap could realize such a steepness. In a purely harmonic lattice model,
\[
V_{\rm lattice}(r)=\frac12 k_{\rm eff}r_{\rm phys}^2,
\qquad
\partial_r \ln V_{\rm lattice}=\frac{2}{r_{\rm phys}},
\]
so the \(\chi_\lambda\) condition by itself only fixes a geometry ratio,
\[
\chi_{\lambda,\rm lattice}=\frac{2\lambda_{\rm phys}}{r_{\rm phys}},
\qquad
\chi_\lambda\ge 1 \iff r_{\rm phys}\le 2\lambda_{\rm phys}.
\]
Evaluated on the earlier reduced turning point, the report found
\[
r_{\rm turn}^{\rm phys}=0.9128891530\,\lambda_{\rm phys},
\qquad
\chi_{\lambda,\rm lattice}(r_{\rm turn})\approx 2.19084649.
\]
That is an important clarification. The session did **not** prove a stiffness directly from \(\chi_\lambda\). It proved that the reduced turning-point geometry already lies inside the formal \(\chi_\lambda\gtrsim 1\) regime if one identifies the lattice trap with a harmonic interstitial void. A stiffness appears only after adding one more assumption, namely force matching between the harmonic trap and the reduced-model barrier force at the turning point. Under that extra step, the report obtained
\[
k_{\rm eff,req}
=
2.73855812\,
\frac{E_*[\mathrm{eV}]}{\lambda_{\rm phys}^2[\mathrm{\AA}^2]}
\quad \mathrm{eV/\AA^2},
\]
or, if \(\lambda_{\rm phys}=a_{\rm int}/2\),
\[
k_{\rm eff,req}
=
10.95423247\,
\frac{E_*[\mathrm{eV}]}{a_{\rm int}^2[\mathrm{\AA}^2]}
\quad \mathrm{eV/\AA^2}.
\]
So the clean outcome was an explicit stiffness formula, together with a warning that it does **not** come from \(\chi_\lambda\) alone. [Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

The third mapping concerned thermal spin survival. Session II had already identified aligned-spin helicity export as the preferred branch, so the material question was whether an ordinary metal scrambles that spin alignment before the barrier crossing finishes. Session V used the Korringa relation,
\[
T_1 T = \mathcal K_{\rm corr},
\]
and imposed the obvious survival condition
\[
T_1\ge t_{\rm cross}^{\rm phys}.
\]
That gave the ceiling
\[
T_{\max} = \frac{\mathcal K_{\rm corr}}{t_{\rm cross}^{\rm phys}}.
\]
Using the reduced crossing time carried from the cold branch,
\[
t_{\rm cross}^{\rm red}\approx 1.82169718,
\]
the report rewrote the same result as
\[
T_{\max}
=
0.548938655\,
\frac{\mathcal K_{\rm corr}}{t_*}.
\]
The important interpretive boundary here is that \(T_{\max}\) remains symbolic until the reduced time unit \(t_*\) and the material-specific Korringa constant are actually calibrated. So Session V produced a clean formula, not yet a laboratory temperature claim. [Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

Taken together, Session V changed the barrier thread in a subtler way than the earlier sections. It did not make the reduced corridor stronger. It made it **testable** in a condensed-matter sense. After this point, one can ask concrete questions such as: does a candidate host provide fast enough lattice turnover, enough effective confinement stiffness once force-matched, and a large enough Korringa constant relative to the actual physical crossing time? But the session very deliberately stopped short of saying that any named material already passes those tests. That restraint matters, because the compact moving-throat program still treats the real branch realization as open, and the physical unit map \((t_*,\lambda_{\rm phys},E_*,\zeta_{\rm ep})\) remained unfixed by the reduced barrier calculations themselves. See [Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt) and [moving_throat_pde_program_compact.md](sandbox:/mnt/data/moving_throat_pde_program_compact.md).

The Session V artifacts are:
[Material-mapping script](sandbox:/mnt/data/material_mapping_condensed_matter_sympy.py)
[Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

## 7. Synthesis of what this session actually established

The cleanest way to summarize the whole session is to separate what was **shown inside the reduced closure** from what remained **closure-dependent or interpretive**. Inside the reduced closure, the session established a coherent same-charge bypass story with five linked pieces. First, relaxing the far-field suppressions \(J^w\approx 0\), rigid-mouth factorization, and positive-source closure lowered the near-contact effective barrier significantly; the representative stationary run reached
\[
\frac{V_{\rm eff}}{V_{\rm Coul}}\approx 0.31446203
\]
at the smallest sampled separation. Second, the dynamic continuation promoted that lowered branch into an actual scattering problem, where the new branch admitted a finite classical threshold
\[
v_{\rm crit,new}\approx 2.54139063
\]
and a reduced WKB enhancement of about **23.31%** on the default subbarrier slice. Third, the magnetic/helicity audit found that aligned spins exported far more subscale helicity than anti-aligned ones on the chosen mixed-sector closure, with a peak export ratio of about **4.95**. Fourth, the proton-proxy sweep found a one-sided over-barrier stability window once transit was compared against dressing-leg collapse. Fifth, the damped V-leg closure showed that even a slower cold event can survive if sufficiently rapid total shedding is present, with
\[
\gamma_{\rm safe}(v_0=2.6)\approx 6.39417302
\]
and
\[
\gamma_{\rm crit}\approx 6.94311167.
\]
The material mapping then translated that whole reduced corridor into threshold equations for lattice turnover, confinement stiffness, and depolarization temperature. [Stationary report](sandbox:/mnt/data/stress_test_relaxed_constraints_report.txt) [Dynamic report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt) [Stability report](sandbox:/mnt/data/proton_proxy_stability_report.txt) [Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt) [Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

Equally important is what the session **did not** establish. It did not derive a new long-range attractive law between same charges. The barrier audit already ruled out that picture at the static one-port level, and the dynamic audit kept the same structural conclusion: the mixed sector reshapes short-range families and open-system channels, while the first outgoing linear correction is phase lag / pumping rather than conservative barrier lowering. So the live corridor after this session should still be described as a **short-range/open-system bypass**, not as a replacement of Coulomb by some hidden long-range attraction. See [barrier_audit_full.md](sandbox:/mnt/data/barrier_audit_full.md).

Several of the session’s most useful results were therefore conditional on declared closures. The reduced \(V_{\rm eff}(r)\) was built from a specific phenomenological \(U/V\) cross-coupling and a specific sign-changing mouth-source completion. The helicity comparison used a reduced mixed-sector closure rather than a fully solved spin-support theorem. The collapse and damped-collapse times were characteristic estimates rather than exact outputs of a solved moving-throat geometry PDE. The material mapping depended on unfixed physical-unit bridges such as \(t_*\), \(\lambda_{\rm phys}\), \(E_*\), and \(\zeta_{\rm ep}\). The compact program’s status firewall is therefore still the correct one: these results are strongest as **reduced / controlled reductions** and **effective closures**, not yet as a finished moving-throat theorem. See [moving_throat_pde_program_compact.md](sandbox:/mnt/data/moving_throat_pde_program_compact.md).

That distinction also separates direct consequences from higher-level “reactor blueprint” extrapolations. Directly supported consequences include the existence of an explicitly lowered reduced barrier, the presence of a finite dynamic crossing window relative to pure Coulomb, the aligned-spin preference for helicity export on the chosen closure, the existence of a cold-survival threshold once sufficient shedding is included, and the exact condensed-matter threshold equations extracted in Session V. By contrast, statements such as “room-temperature operation is already proved,” “ordinary lattice geometry automatically guarantees the bypass,” or “a specific metal like Pd or Ti already satisfies the thresholds” go beyond what the session actually computed. The write-up should keep those as hypotheses or next-step engineering questions, not as established results.

The most compact honest summary of the live corridor after this session is therefore this: **the same-charge bypass remained alive only as a narrow, short-range, mixed-sector/open-system branch.** It required transverse leakage or unresolved export, allowed the dressing geometry to absorb and then shed energy, preferred aligned-spin helicity export on the chosen closure, and demanded sufficient lattice-side dumping once the cold branch was tested. It did **not** survive as a purely static Maxwell trick, a new asymptotic attraction, or a theorem independent of branch realization. In that sense the session was highly productive: it narrowed a vague loophole into a concrete set of coupled gates.

## 8. Immediate next questions

The first next question is **physical-unit calibration**. The material mapping formulas are useful only once the reduced units are tied to laboratory scales. That means determining or constraining the time unit \(t_*\), the physical localization width \(\lambda_{\rm phys}\), the force-matching scale \(E_*\), and any conversion factor such as \(\zeta_{\rm ep}\) that enters the phenomenological damping map. Without that step, the threshold equations
\[
\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys}\ge \frac{8.73618521}{\zeta_{\rm ep}},
\qquad
k_{\rm eff,req}=2.73855812\,\frac{E_*}{\lambda_{\rm phys}^2},
\qquad
T_{\max}=\frac{\mathcal K_{\rm corr}}{t_{\rm cross}^{\rm phys}}
\]
are structurally informative but not yet falsifiable against a real material.

The second next question is **candidate-material screening**. Once the unit map is fixed tightly enough, the session’s formulas can be compared against actual metals, hydrides, deuterides, stressed alloys, or lattice-engineered hosts. The clean screen would not ask vaguely whether a material is “good.” It would ask whether its effective interstitial geometry, electron-phonon turnover, and Korringa-limited spin lifetime can simultaneously satisfy the three threshold conditions carried out of Session V. The current session intentionally stopped one step before that database comparison.

The third next question is **true moving-throat branch realization**. The compact program is explicit that the main remaining theorem gap is no longer algebraic compression but realization of the actual branch in the completed PDE. For the barrier thread, that means asking whether the real branch actually returns the kind of isotropic one-port / grouped-\(P_2\) / mixed-sector data that the reduced barrier scripts assumed, or whether the static placement / orbit-lock side kills the corridor earlier. The same-charge barrier audit is especially sharp on this point: after all of the static and dynamic narrowing, the first unresolved kill condition remains the **static placement / orbit-lock side**, not the wall-like dynamic window. See [barrier_audit_full.md](sandbox:/mnt/data/barrier_audit_full.md).

The fourth next question is **which effect is the next falsifier**. At the level of the current reduced chain, the most valuable falsifier would probably be one that couples back the full localized Maxwell/mixed block, the wall/support dynamics, and the aligned-helicity branch strongly enough to test whether the reduced damping-and-export picture survives once the hidden channels are no longer packaged phenomenologically. In practical terms, that means pushing the moving-throat reduction one step closer to the actual PDE-selected branch rather than inventing another symbolic closure layer. The compact program already points in that direction: the active bottleneck is branch realization, not yet another residual algebra package. See [moving_throat_pde_program_compact.md](sandbox:/mnt/data/moving_throat_pde_program_compact.md).

Finally, the session’s own internal logic suggests a good order for the next round. First, calibrate units and screen candidate materials. Second, use that calibration to decide whether the threshold equations are even plausibly reachable in a real host. Third, if the condensed-matter side remains plausible, return to the moving-throat PDE and test whether the actual branch delivers the required placement/orbit-lock data instead of the reduced surrogate closures used here. That sequence respects the project’s own methodology: derive the narrowest falsifier first, then go deeper only if the target remains alive.

The main artifacts relevant to those next steps are:
[Barrier audit](sandbox:/mnt/data/barrier_audit_full.md)
[Moving-throat compact program ledger](sandbox:/mnt/data/moving_throat_pde_program_compact.md)
[Material-mapping script](sandbox:/mnt/data/material_mapping_condensed_matter_sympy.py)
[Material-mapping report](sandbox:/mnt/data/material_mapping_condensed_matter_report.txt)

---

# Appendix — Barrier / Stability / Materials Thread

This appendix belongs to the same session write-up as Sections 0–8. Its purpose is practical rather than rhetorical: it collects the symbol dictionary, reduced-unit bookkeeping, artifact inventory, source-document roles, and the explicit non-claims / open theorem gaps that keep the thread readable after the fact.

The appendix is scoped to the **barrier-crossing / stability / damping / materials** chain developed in this session. It is not a complete glossary for the entire toy-model program.

---

## Appendix A. Symbol glossary

### A.1 Parent-theory and brane-reduction symbols

| Symbol | Meaning in this thread |
|---|---|
| \(x^M=(t,x,y,z,w)\) | Bulk spacetime coordinates of the parent \(4+1\) theory. |
| \(\mathbf X=(x,y,z,w)\) | Bulk spatial coordinate. |
| \(\mathbf x=(x,y,z)\) | Brane spatial coordinate. |
| \(w\) | Extra transverse / bulk coordinate. |
| \(\psi(\mathbf X,t)\) | Gauged GNLS matter/order parameter field. |
| \(\rho=|\psi|^2\) | Bulk matter density. |
| \(A_M\) | Localized \(4+1\) gauge potential. |
| \(F_{MN}\) | Gauge-field strength tensor. |
| \(Z(w)\) | EM localization profile appearing in the action. |
| \(Z_{\rm int}=\int Z(w)\,dw\) | Integrated localization weight controlling brane-effective EM coupling. |
| \(W(w)\) | Projection kernel defining what the brane observer measures. |
| \(\mu_0^{\rm eff}=\mu_0/Z_{\rm int}\) | Controlled brane-effective Maxwell coupling. |
| \(\eta_Q\) | Fixed topological puncture orientation carrying electric-charge sign. |
| \(e_*>0\) | Positive microscopic charge scale. |
| \(q_*=\eta_Q e_*\) | Signed microscopic branch charge. |
| \(q_{\rm eff}=q_*/\sqrt{Z_{\rm int}}\) | Observable brane charge after canonical zero-mode normalization. |
| \(j^A\) | Bulk number-current components, including \(j^w\). |
| \(S_{\rm leak}\) | Exact projected leakage source in brane continuity. |
| \(\varphi\) | Longitudinal brane velocity potential appearing in the exact longitudinal identity. |

### A.2 Mixed-sector / beyond-MHD symbols

| Symbol | Meaning in this thread |
|---|---|
| \(A_w\) | Extra gauge component that is suppressed only in the strict far-field brane limit. |
| \(F_{\mu w}\) | Mixed field-strength components. |
| \(E_w=-\partial_t A_w-\partial_w A_0\) | Transverse electric field. |
| \(C_a=F_{aw}\) | Mixed spatial field component. |
| \(J^w\) | Transverse electric current. |
| \(J^wE_w\) | Scalar-photon work channel that tracks energy export into the hidden sector. |
| \(v^w\) | Transverse velocity component driven by \(E_w\) when \(J^w\neq 0\). |
| \(h_{\rm sub}\) | Subscale / unresolved helicity, defined from projected minus resolved helicity. |
| \(\partial_t h_{\rm sub}\) | Helicity-export rate used as the spin-orientation diagnostic in Session II. |

### A.3 Moving-throat and bundle symbols

| Symbol | Meaning in this thread |
|---|---|
| \(R(\Omega,w,t)\) | Moving-throat shape field in the geometry lift. |
| \(\Sigma=r-R(\Omega,w,t)\) | Level-set representation of the throat surface. |
| \(a(t),L(t)\) | Old collective geometry variables, reinterpreted as low moments of the distributed throat. |
| \(q_{2m}\) | Grouped real \(P_2\) wall/support amplitudes in the moving-throat program. |
| \(K_*=K-C^2/\varpi^2\) | Effective wall stiffness after integrating out a stable support mode. |
| \(\Delta=\Omega_U^2\Omega_W^2-R^2\) | Internal one-port Maxwell/mixed determinant. |
| \(Q=G_U^2\Omega_W^2+2G_UG_WR+G_W^2\Omega_U^2\) | Static one-port numerator entering the conservative wall operator. |
| \(P=\Omega_U^2G_W+RG_U\) | Static outgoing-load factor in the one-port mixed bundle. |
| \(D_0=K_*-Q/\Delta\) | Static conservative wall operator in the one-port bundle. |
| \(N_0=P^2/\Delta^2\) | Static outgoing-transfer coefficient. |
| \(P_0=N_0/D_0\) | Static normalized outgoing prefactor that also enters the quadrupole-normalization chain. |
| \(m_{\hat 0}\) | Source-map factor in the outgoing-normalization product \(m_{\hat 0}^2P_0\). |

### A.4 Session-I stationary barrier symbols

| Symbol | Meaning in this thread |
|---|---|
| \(r\) | Reduced same-charge separation coordinate. |
| \(r_{\rm reg}\) | Short-distance regularization length used in the reduced barrier ansatz. |
| \(\lambda\) | Reduced localization / width parameter in the stress-test closures. |
| \(U\) | Logarithmic transfer-shape coordinate introduced when rigid-mouth factorization is relaxed. |
| \(V\) | Logarithmic dressing-leg coordinate introduced together with \(U\). |
| \(\epsilon_\eta\) | Dressing variable whose invariance is broken once \(U/V\) coupling is allowed. |
| \(\Xi_1\) | First-order dynamic barrier scalar; in Session I it matched the reduced \(U(r)\) branch. |
| \(\widetilde\Xi\) | Finite nonlinear rigid-mouth proxy carried in the stationary script. |
| \(a_U,a_V\) | Quadratic free-energy coefficients for the reduced \(U,V\) closure. |
| \(g_{UV}\) | Cross-coupling strength between \(U\) and \(V\). |
| \(s_0\) | Amplitude scale in the reduced transfer-source profile. |
| \(\sigma(z)\) | Axial mouth/core source profile, relaxed to a sign-changing multimode branch. |
| \(\mathfrak g[\sigma]\) | Mouth-bias functional of the axial source profile. |
| \(\mathcal R[\sigma]\) | Mixed-to-shell loading ratio used in the Family-1 source stress test. |
| \(\beta_Q,\beta_U,\beta_W\) | Primitive short-range source amplitudes entering the reduced barrier kernel. |
| \(\kappa\) | Yukawa inverse-range parameter appearing in \(e^{-2\kappa x}\)-type families. |
| \(V_{\rm eff}(r)\) | Reduced same-charge effective potential after the three constraint relaxations are turned on. |
| \(V_{\rm Coul}(r)=1/r\) | Pure Coulomb comparison branch. |

### A.5 Session-II dynamic and helicity symbols

| Symbol | Meaning in this thread |
|---|---|
| \(m_s\) | Reduced particle mass in the dynamic scattering problem. |
| \(\hbar_{\rm eff}\) | Reduced Planck-like scale entering the WKB formulas. |
| \(r_{\rm turn}\) | Classical turning point on the reduced potential branch. |
| \(V_{\rm peak}\) | Height of the reduced barrier peak. |
| \(r_{\rm peak}\) | Location of the reduced barrier peak. |
| \(E_{\rm sub}\) | Representative subbarrier energy used for the WKB comparison. |
| \(I_{\rm WKB}\) | WKB action integral on the reduced branch. |
| \(T_{\rm new},T_{\rm Coul}\) | Reduced tunneling weights for the lowered branch and pure Coulomb comparison. |
| \(\Omega_{ij}=-(q_s/m_s)F_{ij}\) | Canonical-vorticity identity used to define the magnetic/vortical sector. |
| aligned / anti-aligned | The two reduced spin closures compared in the helicity-export diagnostic. |
| \(\chi_\lambda=\lambda |\partial_r\ln V_{\rm eff}|\) | Gradient trigger for beyond-MHD / transverse-bypass behavior. |
| \(\lambda_{\rm th}(r_{\rm turn})\) | Threshold width extracted from the turning-point gradient. |

### A.6 Session-III stability symbols

| Symbol | Meaning in this thread |
|---|---|
| \(t_{\rm cross}\) | Characteristic crossing time through the lowered barrier region. |
| \(t_{\rm collapse}\) | Characteristic dressing-leg collapse time. |
| \(\mathcal S(E)=t_{\rm cross}/t_{\rm collapse}\) | Stability ratio used to define the Goldilocks window. |
| \(\mu_\eta\) | Wall-inertia scale for the dressing / geometry leg. |
| \(\chi_\lambda^{\rm peak}\) | Steepest logarithmic gradient of the reduced barrier branch. |
| \(E_{\rm safe,min}\) | Lower edge of the stable over-barrier window in the proton-proxy sweep. |
| \(v_{\rm safe,min}\) | Equivalent lower-edge speed. |
| \(\alpha\) | Ratio \(\mu_\eta/m_s\) in the heavy-throat scaling discussion. |

### A.7 Session-IV damping / heat-sink symbols

| Symbol | Meaning in this thread |
|---|---|
| \(\gamma_{\rm tot}\) | Total V-leg shedding rate. |
| \(\gamma_{\rm vac}\) | Vacuum-side shedding / radiation component. |
| \(\gamma_{\rm lattice}\) | Lattice-side shedding / heat-sink component. |
| \(\gamma_{\rm safe}(v_0)\) | Minimum total shedding that stabilizes a specific cold event. |
| \(\gamma_{\rm crit}\) | Envelope-level unconditional-stability threshold where the characteristic collapse time diverges. |
| \(E_{\rm diss}\) | Total dissipated dressing-leg energy during a crossing. |
| \(E_{\rm vac}\) | Vacuum-side share of dissipated energy. |
| \(E_{\rm lattice}\) | Lattice-heat share of dissipated energy. |
| \(v_0\) | Initial inward speed for the cold-event sweep; default value \(v_0=2.6\). |

### A.8 Session-V condensed-matter mapping symbols

| Symbol | Meaning in this thread |
|---|---|
| \(\lambda_{\rm ep}\) | Dimensionless electron-phonon coupling constant. |
| \(\omega_D\) | Debye frequency of the host material. |
| \(\zeta_{\rm ep}\) | Phenomenological proportionality constant in the damping map \(\gamma_{\rm lattice}^{\rm phys}=\zeta_{\rm ep}\lambda_{\rm ep}\omega_D\). |
| \(t_*\) | Physical time-unit conversion between reduced and laboratory time. |
| \(\lambda_{\rm phys}\) | Physical localization width entering the lattice mapping. |
| \(a_{\rm int}\) | Interstitial spacing used when identifying \(\lambda_{\rm phys}=a_{\rm int}/2\). |
| \(k_{\rm eff}\) | Effective harmonic-trap stiffness in the lattice mapping. |
| \(E_*\) | Physical energy scale used in the force-matching stiffness map. |
| \(T_1\) | Spin-lattice relaxation time. |
| \(\mathcal K_{\rm corr}\) | Korringa constant. |
| \(T_{\max}\) | Maximum temperature below which spin alignment survives the crossing time. |

---

## Appendix B. Reduced-unit and parameter ledger

### B.1 Carry-forward structural constants from the source stack

These constants and structural outputs were already frozen before the barrier session and were treated as carry-forward background rather than as newly refitted session parameters.

| Quantity | Carry-forward value / meaning | Source role |
|---|---|---|
| \(n\) | \(5\) | Frozen stiff-polytropic EOS exponent. |
| \(\kappa_\rho\) | \(1\) | Newtonian mass-dressing coefficient in the gravity-side ledger. |
| \(\kappa_{\rm add}\) | \(1/2\) | Added-mass topology coefficient for the \(w\)-uniform throat. |
| \(\kappa_{\rm PV}\) | \(3/2\) | Adiabatic one-DOF pressure-volume response coefficient. |
| \(\beta_{\rm 1PN}\) | \(3\) | Conservative 1PN precession ledger. |
| \(\alpha^2\) | \(3/4\) | Wake-mixing weight in the vector cross sector. |
| \(a_H\) | \(0\) | Helical admixture set to zero on the carried parity-even branch. |
| \(K_{\rm vec}\) | \(2/\pi^2\) | Vector-sector normalization fixed by the 1PN chain. |
| \(L/a\) | \(\approx 1.85\) | Preferred throat aspect ratio from the EM cavity/throat geometry branch. |
| \(\mu_0^{\rm eff}\) | \(\mu_0/Z_{\rm int}\) | Controlled brane-effective Maxwell coupling. |
| \(q_{\rm eff}\) | \(q_*/\sqrt{Z_{\rm int}}\) | Controlled brane-effective charge. |
| One-port admissibility | \(\Delta>0,\;D_0>0\) | Same bundle conditions used by the barrier audit and normalization chain. |
| Static mixed-kernel families | \(x^{-6},\;e^{-2\kappa x}/x^4,\;e^{-4\kappa x}/x^2\) | Barrier-audit restriction on the static same-charge corridor. |
| First linear dynamic mixed correction | phase lag / pumping, not conservative barrier lowering | Barrier-audit restriction on the linear outgoing channel. |

### B.2 Session-I stationary relaxed-constraint run

#### B.2.1 Parameter set

| Parameter | Value |
|---|---:|
| \(\lambda\) (`lam`) | 1.0 |
| \(r_{\rm reg}\) | 0.25 |
| \(E_{0,\rm amp}\) | 0.18 |
| \(\rho_0\) | 1.0 |
| \(\mu_w\) | 0.8 |
| \(q\) | 1.0 |
| \(a_U\) | 2.5 |
| \(a_V\) | 3.0 |
| \(g_{UV}\) | 0.95 |
| \(s_0\) | 0.9 |
| \(\epsilon_{\rm ref}\) | 0.3 |
| \(K_*\) | 4.0 |
| \(\Omega_U^2\) | 9.0 |
| \(\Omega_W^2\) | 16.0 |
| \(G_U\) | 1.0 |
| \(G_W\) | 1.25 |
| \(R_{\rm mix}\) | 1.35 |
| \(\beta_Q\) | 0.03 |
| \(\beta_{U0}\) | 0.15 |
| \(\beta_{W0}\) | 0.20 |
| \(\kappa\) | 1.0 |
| \(a_0\) | 2.2 |
| \(b_0\) | -0.6 |
| \(r_\sigma\) | 0.8 |
| \(\xi_R\) | 0.9 |
| \(\eta_{\rm leak}\) | 0.03 |
| \(\eta_{UV}\) | 0.22 |
| \(kk_{\rm amp}\) | 0.0 |
| \(r_{\rm F1}\) | 1.77799353547498 |

#### B.2.2 Headline outputs

| Observable | Value |
|---|---:|
| \(\Delta\) | 142.17750000 |
| \(D_0\) | 3.76481862 |
| \(r_{\rm eval}\) | 1.00217028 |
| \(\Xi_1(r_{\rm eval})\) | 0.14313458 |
| \(\widetilde\Xi(r_{\rm eval})\) | 0.14352690 |
| \(U(r_{\rm eval})\) | 0.14313458 |
| \(V(r_{\rm eval})\) | -0.03619791 |
| \(\epsilon_\eta(r_{\rm eval})\) | 0.28933482 |
| \(S_{\rm leak}(r_{\rm eval})\) | 0.03663918 |
| \(J^wE_w(r_{\rm eval})\) | 0.02108684 |
| \(\mathfrak g[\sigma](r_{\rm eval})\) | 0.82823667 |
| \(\mathcal R[\sigma](r_{\rm eval})\) | 0.21677037 |
| \(\sigma_{\min}\) | -0.08979545 |
| sign-changing branch? | True |
| strongest softening point \(r_{\rm soft}\) | 0.18000000 |
| \(V_{\rm eff}(r_{\rm soft})\) | 1.74701126 |
| \(V_{\rm Coul}(r_{\rm soft})\) | 5.55555556 |
| \(V_{\rm eff}/V_{\rm Coul}\) | 0.31446203 |
| UV energy drop at \(r_{\rm soft}\) | 0.21064278 |
| bulk work at \(r_{\rm soft}\) | 1.51632107 |

### B.3 Session-II dynamic scattering / helicity run

#### B.3.1 Additional dynamic parameters

| Parameter | Value |
|---|---:|
| \(m_s\) | 1.0 |
| \(\hbar_{\rm eff}\) | 1.0 |
| \(r_0\) | 5.0 |
| \(r_{\rm contact}\) | 0.18 |
| \(E_{\rm sub}\) | 2.5 |
| `cross_factor` | 1.02 |
| \(\xi_{\rm spin}\) | 0.4 |
| \(C_{0,\rm spin}\) | 0.05 |
| \(r_{\rm spin}\) | 0.8 |
| \(r_{\rm core}\) | 0.15 |
| \(\eta_h\) | 0.25 |
| \(\Omega_{0,\rm sub}\) | 0.7 |
| \(dt\) | 0.0001 |
| \(t_{\max}\) | 10.0 |

#### B.3.2 Headline outputs

| Observable | Value |
|---|---:|
| \(r_{\rm peak}\) | 0.23944389 |
| \(V_{\rm peak}\) | 3.42933112 |
| \(V_{\rm eff}(r_0=5)\) | 0.19999794 |
| Coulomb \(V(r_0)\) | 0.20000000 |
| \(v_{\rm crit,new}\) | 2.54139063 |
| Coulomb contact speed to \(r_{\rm contact}\) | 3.27278339 |
| \(v_{0,\rm sub}\) | 2.14476202 |
| \(r_{\rm turn,new}\) | 0.39096144 |
| \(r_{\rm turn,Coul}\) | 0.40000141 |
| inner turning point | 0.19039548 |
| \(I_{\rm new}\) | 0.19744614 |
| \(I_{\rm Coul}\) | 0.30222297 |
| \(T_{\rm new}\) | \(6.73752615\times 10^{-1}\) |
| \(T_{\rm Coul}\) | \(5.46377065\times 10^{-1}\) |
| \(T_{\rm new}/T_{\rm Coul}\) | 1.23312756 |
| fusion-probability increase | 23.3128% |
| \(\Xi_1(r_{\rm turn})\) | 0.34437471 |
| \(\lambda_{\rm th}(r_{\rm turn})\) | 0.42826825 |
| above-threshold demo \(v_{0,\rm cross}\) | 2.59221845 |
| Coulomb turning point at same \(v_0\) | 0.28091705 |
| aligned \(\max(\partial_t h_{\rm sub})\) | 281.79830789 |
| anti-aligned \(\max(\partial_t h_{\rm sub})\) | 56.96878122 |
| peak helicity-export ratio | 4.94653917 |
| aligned \(h_{\rm sub,final}\) | 20.58070146 |
| anti-aligned \(h_{\rm sub,final}\) | 5.00843357 |
| integrated helicity-export ratio | 4.10920923 |
| scanned \(\lambda_{\rm th}\) range | [0.40673709, 1.06949146] |
| scanned \(\Xi_1\) range | [0.25095422, 0.53817934] |
| scanned WKB-multiplier range | [1.18016972, 1.31627906] |

### B.4 Session-III proton-proxy stability run

#### B.4.1 Stability-sweep additions

| Parameter | Value |
|---|---:|
| \(E_{\rm sub,reference}\) | 2.5 |
| proton mass ratio | 1836.15267343 |
| \(v\)-multiplier min | 1.001 |
| \(v\)-multiplier max | 5.0 |
| number of energy points | 350 |
| number of dynamic samples | 24 |
| \(dt_{\rm dynamic}\) | 0.002 |
| \(t_{\max,\rm dynamic}\) | 500.0 |
| `lambda_choice` | trigger |

#### B.4.2 Headline outputs

| Observable | Value |
|---|---:|
| trigger width \(\lambda_{\rm th}(r_{\rm turn})\) | 0.42826825 |
| chosen \(\lambda_{\rm eff}\) | 0.42826825 |
| proton-proxy \(m_s\) | 1836.15267343 |
| proton-proxy \(\mu_\eta\) | 1836.15267343 |
| proton-proxy threshold speed \(v_{\rm crit,p}\) | 0.05930851 |
| steepest-gradient location | 0.18000000 |
| \(\chi_\lambda^{\rm peak}\) | 21.73204372 |
| \(t_{\rm collapse}\) | 9.43066476 |
| analytic \(E_{\rm safe,min}\) | 5.32265943 |
| analytic \(v_{\rm safe,min}\) | 0.07469791 |
| lower edge in threshold units | \(1.25948037\,v_{\rm crit,p}\) |
| numeric lower edge on scan | 5.36393605 |
| numeric upper edge on scan | 80.93332737 |
| aligned min transit | 0.20400000 |
| aligned max transit | 4.05400000 |
| raw-width \(\chi_\lambda^{\rm peak}\) | 50.74399964 |
| raw-width \(t_{\rm collapse}\) | 6.17163516 |
| raw-width \(E_{\rm safe,min}\) | 27.53273095 |
| reported main safe window | [5.32265943, 80.93332737] |
| equivalent speed window | [0.07469791, 0.29654256] |

### B.5 Session-IV damped V-leg run

#### B.5.1 Damping-specific additions

| Parameter | Value |
|---|---:|
| \(v_{0,\rm cold}\) | 2.6 |
| \(\mu_\eta\) | 1.0 |
| vacuum fraction | 0.25 |
| `gamma_scan_factor_max` | 1.4 |
| number of \(\gamma\) points | 321 |
| `lambda_choice` | model |

#### B.5.2 Headline outputs

| Observable | Value |
|---|---:|
| \(r_{\rm peak}\) | 0.23944389 |
| \(V_{\rm peak}\) | 3.42933112 |
| \(V_{\rm eff}(r_0)\) | 0.19999794 |
| \(v_{\rm crit}\) | 2.54139063 |
| chosen cold speed \(v_0\) | 2.60000000 |
| cold-event energy \(E_{\rm cold}\) | 3.57999794 |
| lowered-branch status | contact |
| pure Coulomb status | turn |
| pure Coulomb turning radius at \(v_0=2.6\) | 0.27933174 |
| actual lowered-branch time to contact | 2.11290000 |
| width used in \(t_{\rm cross}\) | 1.00000000 |
| trigger-width cross-check | 0.42826825 |
| \(\chi_\lambda^{\rm peak}\) | 50.74399964 |
| undamped \(t_{\rm collapse,0}\) | 0.14402764 |
| characteristic \(t_{\rm cross}(v_0=2.6)\) | 1.82169718 |
| undamped stability ratio \(\mathcal S(0)\) | 12.64824697 |
| \(\gamma_{\rm crit}\) | 6.94311167 |
| \(\gamma_{\rm safe}(v_0=2.6)\) | 6.39417302 |
| \(\gamma_{\rm vac}/\gamma_{\rm tot}\) | 0.25 |
| \(\gamma_{\rm lattice}/\gamma_{\rm tot}\) | 0.75 |
| \(\gamma_{\rm vac}\) at \(\gamma_{\rm safe}\) | 1.59854325 |
| \(\gamma_{\rm lattice}\) at \(\gamma_{\rm safe}\) | 4.79562976 |
| \(E_{\rm diss}\) at \(\gamma_{\rm safe}\) | 0.01033460 |
| \(E_{\rm vac}\) at \(\gamma_{\rm safe}\) | 0.00258365 |
| \(E_{\rm lattice}\) at \(\gamma_{\rm safe}\) | 0.00775095 |
| \(E_{\rm store,final}\) at \(\gamma_{\rm safe}\) | 0.00960157 |
| \(\gamma_{\rm vac}\) at \(\gamma_{\rm crit}\) | 1.73577792 |
| \(\gamma_{\rm lattice}\) at \(\gamma_{\rm crit}\) | 5.20733375 |
| \(E_{\rm lattice}\) at \(\gamma_{\rm crit}\) | 0.00770830 |

### B.6 Session-V material-mapping inputs and outputs

#### B.6.1 Inputs carried from earlier reduced runs

| Quantity | Value |
|---|---:|
| required \(\gamma_{\rm lattice}\) | 4.79562976 |
| reduced \(t_{\rm cross}\) | 1.82169718 |
| reduced \(\lambda_{\rm th}\) | 0.42826825 |
| reduced \(r_{\rm turn}\) | 0.39096144 |
| reduced \(V_{\rm eff}(r_{\rm turn})\) | 2.50000000 |

#### B.6.2 Headline mapping formulas and benchmark numbers

| Mapping output | Result |
|---|---|
| \((\lambda_{\rm ep}\omega_D)_{\min}\) | \(4.79562976/(t_*\zeta_{\rm ep})\) |
| equivalent threshold | \(8.73618521011608/(t_{\rm cross}^{\rm phys}\zeta_{\rm ep})\) |
| turnover product | \(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys}\ge 8.73618521011608/\zeta_{\rm ep}\) |
| harmonic-trap geometry ratio | \(\chi_{\lambda,\rm lattice}=2\lambda_{\rm phys}/r_{\rm phys}\) |
| threshold geometry condition | \(r_{\rm phys}\le 2\lambda_{\rm phys}\) |
| \(r_{\rm turn}^{\rm phys}/\lambda_{\rm phys}\) | 0.912889153001653 |
| \(\chi_{\lambda,\rm lattice}(r_{\rm turn})\) | 2.19084649 |
| force-matched stiffness | \(k_{\rm eff,req}=2.7385581171381\,E_*[\mathrm{eV}]/\lambda_{\rm phys}^2[\mathrm{\AA}^2]\) |
| interstitial-spacing form | \(k_{\rm eff,req}=10.9542324685524\,E_*[\mathrm{eV}]/a_{\rm int}^2[\mathrm{\AA}^2]\) |
| \(T_{\max}\) | \(\mathcal K_{\rm corr}/t_{\rm cross}^{\rm phys}\) |
| reduced-unit form of \(T_{\max}\) | \(0.548938655106224\,\mathcal K_{\rm corr}/t_*\) |

### B.7 Reduced-unit conventions that stayed unresolved

The session intentionally left several physical conversion factors unresolved. They are recorded here because many later interpretation mistakes come from forgetting them.

| Quantity | Status in this session |
|---|---|
| \(t_*\) | Unfixed reduced-to-physical time conversion. |
| \(\lambda_{\rm phys}\) | Unfixed physical interpretation of the reduced width \(\lambda\). |
| \(E_*\) | Unfixed physical energy scale used in the force-matching stiffness map. |
| \(\zeta_{\rm ep}\) | Unfixed proportionality between phenomenological lattice damping and \(\lambda_{\rm ep}\omega_D\). |
| \(\mathcal K_{\rm corr}\) | Material-specific Korringa constant, not inserted numerically in the session. |
| candidate-host elastic / electronic data | Not yet screened against the derived thresholds. |

---

## Appendix C. Artifact index

### C.1 Session write-up files

| File | Role |
|---|---|
| `session_barrier_thread_toc.md` | Session table of contents used to organize the thread. |
| `session_barrier_sections_0_2.md` | Draft write-up of Sections 0–2. |
| `session_barrier_sections_3_5.md` | Draft write-up of Sections 3–5. |
| `session_barrier_sections_6_8.md` | Draft write-up of Sections 6–8. |
| `session_barrier_writeup_0_8_clean.md` | Clean merged main write-up without appendix. |

### C.2 Core source-stack files used in the write-up

| File | Role in the barrier thread |
|---|---|
| `4d_summary.md` | Parent \(4+1\) action, projection vs reduction, leakage identity, charge ontology, brane Maxwell hook. |
| `4d_em_fields_summary.md` | Localized Maxwell sector, \(q_{\rm eff}\), \(\mu_0^{\rm eff}\), KK/Yukawa corrections, current-conservation structure. |
| `4d_plasma_summary.md` | Beyond-MHD mixed channels, \(E_w\), \(C_a\), \(J^w\), projected leakage, open-system interpretation. |
| `moving_throat_pde_program_compact.md` | Claim-status firewall, grouped real \(P_2\) program status, open branch-realization bottleneck. |
| `barrier_audit_full.md` | Same-charge static/dynamic mixed-kernel restrictions and corridor narrowing. |

### C.3 Session-I computational artifacts

| File | Role |
|---|---|
| `stress_test_relaxed_constraints_sympy.py` | SymPy stress test for relaxed stationary constraints. |
| `stress_test_relaxed_constraints_report.txt` | Stationary run report with formulas, parameters, and outputs. |
| `relaxed_constraints_veff.png` | Plot of the reduced \(V_{\rm eff}(r)\) against Coulomb. |
| `relaxed_constraints_diagnostics.png` | Leakage / \(\Xi_1\) / mouth-response diagnostics for Session I. |

### C.4 Session-II computational artifacts

| File | Role |
|---|---|
| `dynamic_scattering_helicity_lambda_sympy.py` | Dynamic scattering and helicity continuation. |
| `dynamic_scattering_helicity_lambda_report.txt` | Turning-point, WKB, helicity-export, and \(\lambda_{\rm th}\) report. |
| `dynamic_scattering_potential.png` | Potential plot for the dynamic continuation. |
| `dynamic_scattering_trajectories.png` | Trajectory comparison for lowered branch vs Coulomb. |
| `dynamic_helicity_export.png` | Helicity-export comparison for aligned vs anti-aligned closures. |
| `lambda_threshold_curve.png` | \(\lambda_{\rm th}\) / turning-point trigger curve. |

### C.5 Session-III computational artifacts

| File | Role |
|---|---|
| `proton_proxy_stability_sweep_sympy.py` | Proton-proxy Goldilocks-window sweep. |
| `proton_proxy_stability_report.txt` | Stability-ratio report and safe-window summary. |
| `proton_proxy_stability_timescales.png` | Crossing-vs-collapse timescale plot. |

### C.6 Session-IV computational artifacts

| File | Role |
|---|---|
| `damped_shedding_cold_sweep_sympy.py` | Damped V-leg / cold-sweep script. |
| `damped_shedding_report.txt` | Damped-shedding report with \(\gamma_{\rm safe}\), \(\gamma_{\rm crit}\), and heat partition. |
| `damped_shedding_timescales.png` | Timescale comparison under damping. |
| `damped_shedding_heat_partition.png` | Vacuum-vs-lattice heat partition plot. |

### C.7 Session-V computational artifacts

| File | Role |
|---|---|
| `material_mapping_condensed_matter_sympy.py` | Condensed-matter parameter map script. |
| `material_mapping_condensed_matter_report.txt` | Threshold equations for \(\lambda_{\rm ep}\omega_D\), \(k_{\rm eff}\), and \(T_{\max}\). |

### C.8 Workspace note on expected-but-missing plotting artifacts

The current workspace snapshot did **not** contain the three material-mapping plot files that were provisionally named in the early table of contents (`material_mapping_lambdaep_vs_Debye.png`, `material_mapping_keff_vs_lambda.png`, `material_mapping_Tmax_vs_Korringa.png`). The session thread therefore has script-plus-report coverage for Session V, but not companion plot files under those exact names in the present directory snapshot.

---

## Appendix D. Source-document map

### D.1 Parent-theory and reduction backbone

| Source file | What it contributed to this session |
|---|---|
| `4d_summary.md` | Provided the exact parent GNLS + localized Maxwell + geometry action, the projection/reduction distinction, the exact projected leakage source, the exact longitudinal identity, and the corrected charge ontology used throughout the session. |
| `4d_em_fields_summary.md` | Supplied the localized \(4+1\) Maxwell equations, the controlled brane reduction to \(\mu_0^{\rm eff}=\mu_0/Z_{\rm int}\), the canonical \(q_{\rm eff}\) formula, and the KK/Yukawa packaging that underlies the Coulomb-plus-short-range view of the EM sector. |
| `4d_plasma_summary.md` | Made the beyond-MHD mixed sector explicit: \(A_w\), \(F_{\mu w}\), \(E_w\), \(C_a\), \(J^w\), projected leakage, and the interpretation of brane non-ideality as conservative export into hidden or higher-mode channels. |

### D.2 Moving-throat and barrier-audit backbone

| Source file | What it contributed to this session |
|---|---|
| `moving_throat_pde_program_compact.md` | Supplied the session’s status firewall: what counts as exact, exact within closure, reduced/controlled, effective closure, or open. Also stated explicitly that the main remaining theorem gap is **branch realization** on the completed moving-throat PDE. |
| `barrier_audit_full.md` | Narrowed the same-charge corridor before the session began: the static one-port mixed bundle only renormalizes short-range families, and the first linear outgoing correction is phase lag / pumping rather than conservative long-range barrier lowering. |

### D.3 Session-generated reports

| Source file | What it contributed to this session |
|---|---|
| `stress_test_relaxed_constraints_report.txt` | Recorded the stationary relaxed-constraint closure, the explicit \(S_{\rm leak}\) and \(J^wE_w\) channels, the \(U/V\) drain, the sign-changing source branch, and the \(V_{\rm eff}/V_{\rm Coul}\approx 0.314\) near-contact softening. |
| `dynamic_scattering_helicity_lambda_report.txt` | Recorded the dynamic continuation: barrier peak, threshold speed, turning points, WKB enhancement, aligned-vs-anti-aligned helicity export, and the \(\lambda_{\rm th}\) trigger curve. |
| `proton_proxy_stability_report.txt` | Recorded the Goldilocks timescale sweep under proton-proxy mass/inertia scaling, including \(t_{\rm cross}\), \(t_{\rm collapse}\), and the one-sided safe window. |
| `damped_shedding_report.txt` | Recorded the damped V-leg closure, \(\gamma_{\rm safe}\), \(\gamma_{\rm crit}\), and the split of dissipated energy into vacuum and lattice channels. |
| `material_mapping_condensed_matter_report.txt` | Recorded the condensed-matter thresholds for lattice turnover, harmonic-trap stiffness, and Korringa-limited temperature ceiling. |

### D.4 How the source map is supposed to be read

The session used the core source stack to constrain **what kind of branch was even allowed** and then used the session reports to explore one reduced barrier/stability/materials branch numerically. The source stack did not “prove” the final session outputs directly. It defined the allowed ontology, the admissible closures, and the already-existing restrictions that the session outputs were expected to respect.

---

## Appendix E. Non-claims and open theorem gaps

This appendix is deliberately blunt. It records what the session **did not** show, because many later misunderstandings come from reading a reduced closure as if it were a theorem of the full moving-throat PDE.

### E.1 The session did not solve the full moving-throat PDE

The compact moving-throat program still treats the actual branch realization problem as open. The session’s barrier/stability/materials chain therefore lives inside a hierarchy of reduced closures and phenomenological completions. It is strongest as a **reduced / controlled reduction** plus **effective closure** story, not as a completed theorem of the fully solved bulk/interface PDE.

### E.2 The session did not discover a new long-range same-charge attraction

The barrier audit already constrained the static same-charge mixed bundle to short-range families. Nothing in Sessions I–V overturned that. The stationary reduction lowered the near-contact barrier, but it did so through short-range families, leakage/export channels, dressing-leg drain, and source compensation. The dynamic audit likewise did not turn the first linear outgoing correction into a new conservative long-range force. The live corridor remained a **short-range/open-system bypass**, not a replacement of Coulomb by a hidden attractive \(1/r\)-law.

### E.3 The stationary relaxed-constraint branch is still closure-dependent

The Session-I result depends on a specific phenomenological choice for the \(U/V\) free energy, a specific scaling of the leakage channel, and a specific sign-changing source completion. The fact that the reduced \(V_{\rm eff}(r)\) dropped to about \(0.314\) of Coulomb in the plotted near-contact window is therefore a result **inside that declared closure**, not a branch-independent theorem.

### E.4 The helicity-export preference is not yet a spin theorem

The aligned-vs-anti-aligned comparison in Session II established a strong preference for aligned spins on the chosen reduced mixed-sector closure, because aligned spins exported far more subscale helicity. That is a useful reduced diagnostic, but it is not yet a full theorem of the microscopic spin-support structure of the parent theory.

### E.5 The collapse and damping times are characteristic closures, not solved geometry times

The \(t_{\rm cross}\), \(t_{\rm collapse}\), and \(t_{\rm collapse}^{\rm damped}\) formulas were characteristic estimates chosen to stress-test the branch. In particular, the finite \(\gamma_{\rm crit}\) of Session IV comes from the adopted **envelope** closure. A raw damped inverted-oscillator equation would not, by itself, generate the same finite unconditional-stability threshold. So \(\gamma_{\rm crit}\) should be read as a property of the adopted damped-envelope reduction, not as a closure-free theorem.

### E.6 The proton-proxy window did not prove generic heavy-particle stability

Session III showed that a one-sided stable window opens when the heavy-throat scaling is combined with the trigger-width choice inherited from Session II. It did **not** show that proton-like defects are generically stable simply because they are heavy. The safe lower edge remained highly sensitive to the confinement-width / \(\chi_\lambda\) choice.

### E.7 The condensed-matter map did not yet prove any specific candidate material

Session V stopped at threshold equations. It did not insert real Palladium, Titanium, hydride, deuteride, alloy, or stressed-lattice data and did not prove that any particular host already satisfies all three thresholds simultaneously. The material screen is still a next step, not a completed result.

### E.8 The session did not prove a room-temperature reactor theorem

A few higher-level interpretations suggested that very fast crossing could make spin scrambling irrelevant in ordinary laboratory conditions. The session itself did **not** prove that. It produced the symbolic ceiling
\[
T_{\max}=\frac{\mathcal K_{\rm corr}}{t_{\rm cross}^{\rm phys}},
\]
but it left both \(t_*\) and the material-specific Korringa constant unfixed. So any room-temperature or hot-operation claim remains interpretive until the physical unit map and material constants are inserted.

### E.9 The stiffness map contains an extra force-matching assumption

The \(\chi_\lambda\) criterion alone only fixes the geometry ratio \(r_{\rm phys}/\lambda_{\rm phys}\) for a harmonic trap. The explicit \(k_{\rm eff}\) formula required one more assumption: matching the absolute harmonic-trap force to the reduced-model barrier force at the turning point. That formula is therefore useful, but it is not a pure consequence of \(\chi_\lambda\) by itself.

### E.10 The most important remaining theorem gaps

The remaining gaps are now fairly narrow.

1. **Physical-unit calibration:** determine \(t_*\), \(\lambda_{\rm phys}\), \(E_*\), and \(\zeta_{\rm ep}\).
2. **Candidate-material screening:** test real materials against the derived lattice-turnover, stiffness, and Korringa thresholds.
3. **True branch realization:** determine whether the completed moving-throat PDE actually returns the isotropic / grouped-\(P_2\) / mixed-sector branch that the reduced barrier scripts assumed.
4. **Placement/orbit-lock kill test:** determine whether the actual branch survives the static placement / orbit-lock constraints highlighted by the barrier audit.
5. **Higher-fidelity damping/export closure:** replace the phenomenological V-leg envelope law with a more microscopic wall/support/mixed-sector export calculation.

### E.11 The cleanest falsifier after this appendix

The narrowest next falsifier is probably not “build a full reactor model.” It is:

1. calibrate reduced units tightly enough to turn the Session-V threshold equations into real laboratory inequalities,
2. test whether any plausible host material satisfies them,
3. and, only if that screen survives, return to the moving-throat PDE to see whether the actual branch realization preserves the same corridor once the hidden channels are no longer packaged phenomenologically.

That order keeps the project faithful to its established methodology: derive the tightest kill test first, then go deeper only if the target remains alive.
