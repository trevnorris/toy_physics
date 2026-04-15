# Session barrier thread write-up — Sections 3, 4, and 5

## 3. Session II — Dynamic scattering, turning point, and helicity export

Session II took the stationary lowered-barrier branch from Session I and asked the next obvious question: does the reduction merely soften a plotted potential, or does it open a genuine **dynamical** corridor in which an approaching same-charge particle can get farther than the pure Coulomb model would allow? The answer was sought without changing the basic structural discipline. The run kept the same reduced short-range family structure already allowed by the barrier audit, rather than introducing a new long-range interaction by hand. That mattered because the barrier audit had already shown that the mixed bundle’s first linear dynamic correction is phase lag / pumping rather than a new conservative attraction, so any improvement in approach had to come from the lowered near-contact branch and the new open-system channels, not from inventing a different asymptotic law. fileciteturn17file15

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
and the plasma extension already insists that the mixed fields \(A_w\), \(F_{\mu w}\), and \(J^w\) remain microscopic degrees of freedom outside the strict far-field brane reduction. So adding the Lorentz-force contribution and asking whether unresolved helicity is exported into the higher / transverse sector is consistent with the existing ontology, not a separate theory bolt-on. fileciteturn14file1turn14file2

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

So Session II changed the status of the barrier thread in three ways. It demonstrated a reduced classical-crossing window that the pure Coulomb comparison did not share, it quantified a subbarrier enhancement on the same branch, and it identified aligned-spin helicity export as the preferred magnetic route into the pure-transfer subcorridor. But none of that overrode the broader claim-status firewall. The compact moving-throat program still treats branch realization and the final outgoing normalization as open, so the dynamic scattering and helicity results belong to the same “reduced / controlled reduction plus effective closure” category as the stationary barrier lowering that came before them. fileciteturn14file16

The Session II artifacts are:
[Dynamic SymPy script](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_sympy.py)  
[Dynamic scattering / helicity report](sandbox:/mnt/data/dynamic_scattering_helicity_lambda_report.txt)  
[Potential plot](sandbox:/mnt/data/dynamic_scattering_potential.png)  
[Trajectory plot](sandbox:/mnt/data/dynamic_scattering_trajectories.png)  
[Helicity-export plot](sandbox:/mnt/data/dynamic_helicity_export.png)  
[Lambda-threshold curve](sandbox:/mnt/data/lambda_threshold_curve.png)

## 4. Session III — Structural stability and the Goldilocks window

Once Session II showed that the lowered branch could be crossed dynamically, the next question was no longer “can a particle get through?” but “can it get through **without its own dressing geometry failing first?**” That is a different question, because the earlier relaxed rigid-mouth closure had already shown that energy can be drained into the dressing leg \(V\). Session III therefore reframed the problem as a competition between two characteristic timescales: the barrier-region transit time and the dressing-leg collapse time. In spirit that move was natural, because the moving-throat lift already treats geometry through effective wall inertias and stiffnesses, and the distributed-wall reduction explicitly identifies \(\mu_\eta\) as the wall inertia density feeding the low-mode geometry equations. What remained closure-dependent was the exact choice of reduced “collapse” metric used for the stress test. fileciteturn17file3turn14file16

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

That sensitivity result is also what keeps the interpretation honest. Session III did **not** show that a proton-proxy defect is generically stable whenever it is heavy. It showed something narrower: once the same reduced branch is equipped with the earlier trigger width and the chosen aligned-spin closure, there is a one-sided over-barrier regime in which transit beats collapse, and the survival edge is controlled primarily by \(\chi_\lambda\) and \(\lambda_{\rm eff}\). The moving-throat hierarchy had already prepared that conclusion by making wall inertia and geometry response explicit; the session’s added value was to turn that structure into a concrete timescale inequality and then sweep it numerically. fileciteturn17file3turn14file16

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
and asked whether sufficiently rapid shedding could stabilize a slow aligned-spin event. This move was explicitly phenomenological: the source stack does allow effective damping completions in the geometry sector, but it does not yet derive this specific dressing-leg dissipation law from the completed PDE. So Session IV should be read as an open-system **effective closure** layered on top of the already reduced barrier branch, not as an exact theorem of the parent action. fileciteturn17file3turn14file16

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

What Session IV therefore changed was not the parent theory, but the practical interpretation of the earlier cold-crossing failure. In the undamped picture, the V-leg behaved like an energy trap and the slow aligned-spin event died on the dressing side before it could count as a controlled crossing. In the damped picture, the question becomes quantitative rather than absolute: how large must \(\gamma_{\rm lattice}\) be, relative to the weaker vacuum channel, for the cold event to survive, and how much energy is dumped into the lattice per event once that happens? The report answered both questions in explicit reduced units. That still falls short of a materials theorem, but it is already much more informative than the earlier undamped “speed trap” verdict. It converts a failure mode into a solvable design inequality. fileciteturn17file3turn14file16

The Session IV artifacts are:
[Damped-shedding SymPy script](sandbox:/mnt/data/damped_shedding_cold_sweep_sympy.py)  
[Damped-shedding report](sandbox:/mnt/data/damped_shedding_report.txt)  
[Timescale plot](sandbox:/mnt/data/damped_shedding_timescales.png)  
[Heat-partition plot](sandbox:/mnt/data/damped_shedding_heat_partition.png)
