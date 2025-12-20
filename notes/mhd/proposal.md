**Research Program Proposal**  
**Topology-Controlled Energy Release and Confinement Management**  
**using the Superfluid Defect Toy Model**

*Fluid Universe Project (internal)*  
Date: 19 Dec 2025 (America/Denver)  
Prepared for: future focused sessions and iterative compilation

# **Table of Contents**

1. 1\. Executive summary  
2. 2\. Motivation and scope boundaries  
3. 3\. What the toy model uniquely contributes  
4. 4\. Fusion: current blockers to reliable power  
5. 5\. Research thrusts and work packages  
6. 6\. Simulation program and diagnostics  
7. 7\. Validation plan against fusion-relevant phenomena  
8. 8\. Tritium and advanced-fuel pathways (what we can and cannot change)  
9. 9\. Milestones, resources, and collaboration needs  
10. 10\. Risk register and guardrails  
11. 11\. References (fusion constraints and reconnection literature)

# **1\. Executive summary**

This document proposes a focused research program to leverage the Superfluid Defect Toy Model as a practical ‘conceptual and computational lens’ on fusion-relevant physics—especially the control of topology-changing events (reconnection, tearing, bursty relaxations) that limit reliability in magnetically confined plasmas.

The program does not assume the toy model is a literal plasma model, nor does it attempt to redesign a reactor. Instead, it targets portable insights:  
• a minimal, equation-level framework for when topology is frozen-in vs when reconnection is permitted;  
• prediction of thresholds and phase diagrams for burst onset (gentle relaxation vs runaway cascade);  
• energy-partition rules (what fraction becomes waves/flows/heat) during topology changes;  
• a control vocabulary that connects external ‘pulses’ to non-ideal windows that enable controlled reconnection.

These outputs can inform (i) disruption prediction/mitigation, (ii) ELM/edge-burst management, (iii) heat-load shaping, and (iv) reducing operational margins needed for advanced fuels (D–D, D–He3, p–B11) by improving stability and control.

# **2\. Motivation and scope boundaries**

Fusion reactors struggle not only with achieving fusion conditions, but with doing so continuously, safely, and economically. Two recurring themes in modern fusion engineering are:  
(a) plasmas store large magnetic and thermal energy that can be released abruptly by instabilities; and  
(b) the wall/divertor must survive both steady and transient heat loads.

Your toy model already developed a ‘Rosetta stone’ mapping between hydrodynamic vorticity evolution and the MHD induction equation (including the role of a small diffusive/non-ideal term as the gate that permits topology change). This mapping is valuable even if it is only a kinematic correspondence, because reconnection/tearing theory is largely organized around that induction-law structure.

Scope boundaries (guardrails):  
• We do not claim new nuclear reactions or altered fusion cross-sections.  
• We do not propose ‘vacuum energy extraction’. Energy release is treated strictly as conversion of stored field/flow energy.  
• We do not provide step-by-step reactor construction guidance. The focus is on physics insight, diagnostics, and control concepts.  
• When discussing fuels (e.g., p–B11), we treat them as constraints that amplify the need for stability/control, not as targets of the toy model.

# **3\. What the toy model uniquely contributes**

## **3.1 A kinematic dictionary that makes topology-control questions portable**

The toy model provides an explicit mapping (developed in the project’s EM/MHD bridge notes) in which vector potential A is proportional to the flow velocity v, and magnetic field B corresponds to vorticity ω \= ∇×v. With these definitions, inviscid Euler dynamics implies the ideal ‘frozen-in’ condition E \+ v×B \= 0 and therefore the ideal induction law. Adding a small viscous/diffusive term produces the resistive/diffusive induction slot that permits reconnection. This is the central mechanism that turns “switching” into a precise statement: a pulse toggles access to a non-ideal window.

## **3.2 A principled interpretation of ‘switching’ (fusion/fission) as gating the non-ideal channel**

Experimental systems that show externally controlled fusion vs fission of topological defects can be interpreted through a common lens: (i) an energy landscape with multiple metastable topological states, and (ii) a temporary non-ideal pathway that allows transitions between them. The toy model’s main contribution is to separate these two roles cleanly and to provide a quantitative way to test which one dominates in a given protocol.

## **3.3 A simulation philosophy optimized for topology-change events**

The project’s existing CUDA-style ‘one-off’ Python scripts (e.g., GP reconnection prototypes and energy-partition diagnostics) are already tailored to the hardest part of these studies: distinguishing physical emission (sound/phonons, jets, cascades) from numerical artifacts. This program proposes standardizing those scripts into an experiment runner that supports fast parameter sweeps and produces comparable outputs.

## **3.4 A cross-domain bridge: from defect topology in soft matter to reconnection in plasmas**

Because the toy model is organized around topology change and energy release rather than material microphysics, it can act as a bridge between soft-matter topological solitons (where switching is demonstrable and repeatable) and fusion plasmas (where switching is dangerous). This provides a safe pathway for hypothesis generation: demonstrate control strategies in an analog platform, then translate to fusion diagnostics.

# **4\. Fusion: current blockers to reliable power**

## **4.1 Disruptions and runaway-energy events**

Disruptions are abrupt losses of confinement that can dump stored plasma thermal and magnetic energy onto plasma-facing components, leading to severe melting and large electromagnetic forces; runaway electrons are an additional hazard in some scenarios. These risks force conservative operating margins and drive complex mitigation systems (e.g., shattered pellet injection).

## **4.2 Edge localized modes (ELMs) and transient heat loads**

High-performance confinement regimes often produce transient edge events (ELMs) that eject energy and particles toward the divertor. Even if a device can operate at high fusion gain, repeated transient loads threaten divertor lifetime and availability.

## **4.3 Heat exhaust and plasma-facing materials**

Heat exhaust remains a central engineering bottleneck. Public ITER-facing discussions commonly cite divertor components designed to tolerate \~10 MW/m² steady and \~20 MW/m² for short transients, emphasizing how narrow the design window is. DEMO-class plants must manage these loads with high uptime and remote maintainability.

## **4.4 Tritium fuel cycle constraints**

D–T is the near-term ‘easiest’ fusion reaction, but tritium is scarce and decays, requiring a breeding and recovery fuel cycle for power plants. DEMO studies often use a tritium breeding ratio (TBR) criterion above unity (commonly quoted \>1.1 to cover uncertainties) to achieve self-sufficiency. This is a coupled system constraint: plasma performance, blanket design, and availability all affect fuel viability.

## **4.5 Advanced fuels (D–D, D–He3, p–B11): attractive but less forgiving**

More abundant or less neutron-heavy fuel cycles are often discussed, but they generally require higher temperatures and/or better confinement, and can be constrained by radiation losses (notably bremsstrahlung in p–B11 concepts). Therefore, any pathway toward reduced tritium dependence is, in practice, also a pathway toward tighter stability and control requirements.

# **5\. Research thrusts and work packages**

The program is organized as work packages (WPs). Each WP is written so it can become a dedicated future session (or mini-report) with concrete outputs.

| Work package | Title | Purpose / main deliverable |
| :---- | :---- | :---- |
| WP1 | Formalize the topology-switching framework | Produce a clear mathematical and computational statement of the ‘frozen-in’ vs ‘reconnecting’ regimes under the toy model’s MHD dictionary. |
| WP2 | Pulse-gated reconnection phase diagrams | Simulate controlled reconnection/fusion/fission outcomes as a function of pulse amplitude, duration, and localization; build a regime map. |
| WP3 | Invariant and bookkeeping diagnostics (helicity/Hopf-like measures) | Define and measure robust invariants/proxies that track topology through reconnection and quantify where they ‘go’ (linking→twist→dissipation). |
| WP4 | Sheet/tearing/plasmoid onset analogs | Establish thresholds and scaling for ‘single-event’ vs ‘cascade’ reconnection, using current-sheet analog initial conditions. |
| WP5 | Fusion-relevant translation: disruptions, ELMs, and mitigation shaping | Translate WP1–WP4 results into hypotheses for prediction/mitigation: early-warning markers, ‘pacing’ strategies, and load-shaping control ideas. |
| WP6 | Tritium reduction and advanced-fuel margins | Quantify how improved stability/control widens viability margins for non-DT fuels; identify which control levers matter most. |
| WP7 | Analog-platform bridge (soft matter defect switching) | Use experimental defect switching (e.g., electrically driven knot fusion/fission) as a controllable testbed for the same control logic. |
| WP8 | Software unification and reproducible pipelines | Refactor one-off CUDA scripts into a single experiment runner with standardized outputs, postprocessing, and sweep automation. |

## **WP1. Formalize the topology-switching framework**

### **Rationale:**

Reconnection research is organized around a simple but powerful dichotomy: ideal evolution conserves field-line topology; non-ideal effects (diffusion, resistivity, mutual friction, viscosity) open the door to topology change. In the toy model, the MHD dictionary makes this dichotomy explicit at the level of the induction/vorticity equation. WP1 turns that insight into a reusable set of definitions, dimensionless control numbers, and diagnostic criteria.

### **Key questions:**

* What is the most stable, minimal set of equations (toy-level) that still reproduces frozen-in vs reconnecting behavior?  
* Which ‘non-ideal’ term is the best proxy for an externally applied switch (time-dependent coefficient, spatial localization, both)?  
* Which dimensionless numbers best organize behavior (diffusion length vs core/sheet thickness; advective vs diffusive time scales)?

### **Methods and simulations:**

* Extract the core correspondence from mhd.md and the supporting Mathematica scripts; rewrite as a compact ‘model card’.  
* Define a standardized pulse model: η(t) or γ(t) windows; optionally η(x,t) localized to a reconnection region.  
* Define a regime classifier that labels outcomes: no reconnection, single reconnection (gentle), multi-event cascade (violent), breakup.

### **Deliverables:**

* 2–3 page ‘Topology Switching Primer’ with equations and definitions in plain language.  
* Implementation-ready pseudocode for the pulse term and for regime classification.  
* A shortlist of diagnostic plots that must be produced for every run (energy vs time; event markers; topology proxy).

### **Suggested next-session prompt:**

***Let’s finalize WP1: define the minimal equations, the pulse model, and the regime classifier, then write a ‘model card’ we can reuse.***

## **WP2. Pulse-gated reconnection phase diagrams**

### **Rationale:**

Many ‘switching’ experiments (including defect-knot fusion/fission work in soft matter) can be reinterpreted as operating in different regions of a control phase space: pulse strength/duration vs geometry. Fusion devices also use pacing/perturbations to avoid worse events. WP2 builds a concrete regime map from simulations, turning intuition into an actionable control diagram.

### **Key questions:**

* For a fixed initial topology, what pulse amplitude and duration reliably produces a single controlled reconnection?  
* Where is the boundary between ‘nothing happens’ and ‘too much happens’ (cascade)?  
* How sensitive are outcomes to geometry (separation, relative orientation, twist)?

### **Methods and simulations:**

* Use GP-based prototypes (quantum\_vortex\_reconnection.py) and the diagnostic harness pattern in length\_loss.py.  
* Implement pulsed dissipation (γ(t) or η(t)) and optionally localize it in space (Gaussian window around expected sheet).  
* Run parameter sweeps (grid search or Latin hypercube) and auto-classify outcomes to build phase diagrams.

### **Deliverables:**

* Phase diagram plots: outcome class vs (pulse amplitude, pulse duration) for 2–3 canonical initial conditions.  
* Scaling hypothesis: optimal pulse when diffusion length sqrt(η τ) is comparable to core/sheet thickness.  
* A curated set of 5–10 ‘representative runs’ with saved fields for later analysis.

### **Suggested next-session prompt:**

***For WP2, pick the first canonical initial condition (e.g., two linked rings). Define pulse parameters and outputs, then design the sweep grid.***

## **WP3. Invariant and bookkeeping diagnostics (helicity/Hopf-like measures)**

### **Rationale:**

Claims of ‘topological protection’ are only meaningful if we can identify what is approximately conserved and under what conditions. In MHD, magnetic helicity is a standard topological measure; in soft-matter defects, a Hopf index plays a similar role. WP3 creates practical proxies that can be computed robustly in simulations and used as constraints on what reconnection can do.

### **Key questions:**

* Which invariant/proxy is numerically stable in our discretization (global helicity, relative helicity, spectral helicity density)?  
* During reconnection, does the ‘topological currency’ transfer across scales (linking → twist) before dissipating?  
* Can invariants predict whether a pulse will lead to fusion vs fission outcomes?

### **Methods and simulations:**

* Compute phase-derived velocity v ∝ ∇arg(ψ) with masking near cores; decompose into incompressible component.  
* Compute helicity-like integrals (A·B) or spectral helicity proxies; track time evolution and spectra.  
* Correlate invariant changes with event times and with energy-partition bursts.

### **Deliverables:**

* A robust diagnostic module that computes one or more invariant proxies per run.  
* Empirical rules: what remains conserved in gentle reconnection vs what breaks in cascades.  
* A ‘bookkeeping narrative’ to use in future papers: where the topology ‘goes’ when it changes.

### **Suggested next-session prompt:**

***For WP3, let’s decide on a helicity proxy we can compute reliably from ψ, then add it to the length\_loss.py diagnostic pipeline.***

## **WP4. Sheet/tearing/plasmoid onset analogs**

### **Rationale:**

Fusion-relevant reconnection often involves current sheets that become unstable, leading to multi-island/plasmoid cascades. The toy model can serve as a simplified environment to study the onset boundary between single reconnection and explosive multi-event behavior. This is valuable for control: avoid the region where cascades become likely.

### **Key questions:**

* What initial conditions best represent a ‘sheet’ in the toy model (shear layer, anti-parallel tube array, driven boundary)?  
* Is there a sharp onset boundary in a dimensionless number (e.g., effective Lundquist-like ratio)?  
* How does pulse localization affect whether a cascade occurs?

### **Methods and simulations:**

* Construct a thin shear/vorticity layer initial condition; add small diffusion; monitor growth of secondary structures.  
* Measure event statistics: number of reconnections, peak energy release rate, fragmentation metrics.  
* Test localized vs global diffusion pulses; compare outcomes.

### **Deliverables:**

* An onset curve separating single-event from cascade regimes.  
* Qualitative scaling comparisons with known reconnection/plasmoid heuristics (without over-claiming identity).  
* Recommendations for ‘operate near threshold’ control strategies.

### **Suggested next-session prompt:**

***For WP4, define the simplest sheet initial condition we can implement on your grid and specify what we’ll measure to detect a cascade.***

## **WP5. Fusion-relevant translation: disruptions, ELMs, and mitigation shaping**

### **Rationale:**

Tokamaks already use mitigation strategies (pellet injection, pacing, magnetic perturbations), but prediction/control remains challenging. WP5 translates toy-model phase diagrams and diagnostics into fusion-relevant hypotheses: early-warning indicators, and control waveforms that trade rare catastrophic events for frequent benign ones.

### **Key questions:**

* Which toy-model diagnostics map to quantities measurable in tokamaks (magnetic fluctuation spectra, mode amplitudes, energy decay rates)?  
* Can we propose a generic ‘pacing’ criterion: intervene when a sheet-thickness proxy crosses a threshold?  
* How can control aim to ‘radiate’ or distribute energy rather than deposit it locally?

### **Methods and simulations:**

* Create a ‘translation table’ from toy-model diagnostics to MHD/tokamak observables.  
* Use WP2–WP4 results to propose control waveforms (pulse shapes) that prefer gentle relaxation.  
* Compare qualitative predictions with published mitigation approaches (e.g., shattered pellet injection aims to convert energy into radiation).

### **Deliverables:**

* 1–2 page hypothesis list: ‘If you see X, intervene with Y’ (high level, not engineering instructions).  
* A prioritized set of measurable markers to evaluate against fusion data (public shots where possible).  
* A draft outline for a cross-disciplinary paper or note.

### **Suggested next-session prompt:**

***For WP5, let’s build a translation table: our diagnostics → tokamak observables, and then draft 3–5 testable mitigation hypotheses.***

## **WP6. Tritium reduction and advanced-fuel margins**

### **Rationale:**

The toy model cannot change nuclear cross sections, but it can help in the place advanced fuels fail most often: tight margins. If better control reduces the frequency/severity of disruptive events, it increases availability and lowers fuel throughput needs. Similarly, if control strategies improve energy partition (less into damaging wall loads), it makes more demanding fuels less implausible.

### **Key questions:**

* Which improvements (reduced disruptions, reduced transient loads, better confinement stability) most directly reduce tritium inventory/throughput needs?  
* Can we define a ‘control margin’ metric: how far from an instability boundary a device can operate while maintaining uptime?  
* For p–B11 and D–He3 concepts, which loss channels are most sensitive to bursty events (e.g., radiative losses, electron heating)?

### **Methods and simulations:**

* Define an abstract plant-level metric: availability × safe operating window; relate to event rate predicted by WP2–WP4.  
* Use simplified loss models (not nuclear design) to show how reduced burst frequency improves net power prospects.  
* Write clear ‘what we can/cannot change’ statements to avoid overstated claims.

### **Deliverables:**

* A framing note: ‘Toy-model contributions to tritium challenge’ (focus: control and availability).  
* A list of advanced-fuel ‘margin multipliers’ most affected by reconnection control.  
* Concrete next targets: which WP results to extend if we later choose a specific fuel pathway.

### **Suggested next-session prompt:**

***For WP6, let’s write a clean tritium/advanced-fuel section: what constraints are nuclear/engineering vs what control theory can genuinely improve.***

## **WP7. Analog-platform bridge (soft matter defect switching)**

### **Rationale:**

Soft-matter systems demonstrate topological defect switching (fusion/fission) under controlled fields in a way fusion plasmas cannot easily permit. WP7 uses those systems as a safe experimental sandbox to test control concepts (pulse timing, localization, outcome classification) before translating the logic to fusion.

### **Key questions:**

* What control parameters in defect-switching experiments correspond to the toy model’s non-ideal gate (diffusion window) and energy tilts?  
* Which observables in those experiments correspond to our energy-partition diagnostics (waves, dissipation bursts)?  
* Can we propose a minimal ‘switching model’ that predicts fusion vs fission outcomes from pulse protocol?

### **Methods and simulations:**

* Extract the experiment’s control protocol and outcomes; build a reduced model mirroring WP2 phase diagrams.  
* Use toy-model simulations to reproduce qualitatively similar switching maps (without matching microphysics).  
* Identify invariant/bookkeeping parallels (Hopf-like charges vs helicity proxies).

### **Deliverables:**

* A bridging memo: how to interpret defect switching as reconnection control.  
* A proposed set of ‘universal plots’ that both communities can compare.  
* A list of candidate collaborators / experimental groups to follow.

### **Suggested next-session prompt:**

***For WP7, let’s build the universal switching model: define state variables, control pulse, and predicted fusion/fission boundary.***

## **WP8. Software unification and reproducible pipelines**

### **Rationale:**

One-off CUDA scripts are an asset for speed, but without standard outputs it is hard to compare runs and build cumulative knowledge. WP8 standardizes the workflow so that every new idea becomes a comparable experiment, enabling systematic sweeps and a growing dataset.

### **Key questions:**

* What is the minimal common interface for all experiments (grid, time step, pulse schedule, IC generator)?  
* What outputs must every run save (energies, spectra, topology proxies, checkpoints)?  
* How do we ensure numerical trust (dealiasing, convergence checks, control runs)?

### **Methods and simulations:**

* Refactor length\_loss.py into an ‘experiment runner’ with subcommands for initial conditions and pulse models.  
* Add a postprocessing script that generates canonical figures and a run report.  
* Add a configuration schema (YAML/JSON) so sweeps can be queued reproducibly.

### **Deliverables:**

* A clean repo folder: /experiments with configs, /runs with .npz outputs, /reports with auto-generated summaries.  
* A run-report template (markdown/PDF) that captures settings and results in one page.  
* Codex-ready instructions for the refactor and for adding new experiments safely.

### **Suggested next-session prompt:**

***For WP8, let’s define the standard run schema and refactor plan: one config → one run → one report.***

# **6\. Simulation program and diagnostics**

## **6.1 Recommended model ladder**

To avoid overfitting to any one interpretation, we use a ladder of models, from cheapest to richest:  
• Reduced filament/Biot–Savart models (fast sweeps, qualitative outcome maps).  
• Gross–Pitaevskii / compressible superfluid PDE (captures reconnection, sound emission, core physics).  
• Viscous vorticity transport / induction-equation analogs (directly tests the MHD dictionary).  
• Optional: reduced MHD or Hall-MHD toy variants (only if needed for a specific translation).

## **6.2 Standard diagnostics (every run)**

* Total energy vs time, and partition into incompressible vs compressible components (or closest available decomposition).  
* Event markers: reconnection time(s), peak release rate max(-dE/dt), number of reconnection events.  
* Topology proxy/invariant (helicity-like) vs time; optionally spectral helicity density.  
* Spatial ‘where did the energy go’ maps: compressible burst fronts, jet formation, core density depletion volume.  
* Convergence/control checks: pulse-off control, resolution check for representative cases, energy conservation in ideal limit.

## **6.3 Outcome classification**

Each simulation should end with an automatic classification label that supports phase diagrams:  
• No reconnection / topology preserved.  
• Single controlled reconnection (gentle relaxation).  
• Multi-event reconnection cascade (plasmoid-like fragmentation).  
• Fission-like splitting (one object becomes multiple).  
• Fusion-like merging (multiple objects become one).  
• Loss of coherence / turbulent decay.  
Classification can be based on core-line tracking, event counts, and energy-release signatures.

# **7\. Validation plan against fusion-relevant phenomena**

Validation is framed as ‘does the control logic transfer?’, not ‘does the toy model reproduce a tokamak’. Key validation steps:  
1\) Internal validation: ideal runs preserve invariants and conserve energy to numerical tolerance.  
2\) Regime validation: non-ideal pulses produce topology change with diffusion-length scaling.  
3\) Qualitative fusion alignment: toy-model ‘violent vs gentle’ boundaries resemble the need for disruption avoidance and ELM pacing.  
4\) Literature alignment: compare onset thresholds and event statistics with reconnection/instability heuristics in MHD literature.

# **8\. Tritium and advanced-fuel pathways (what we can and cannot change)**

## **8.1 What is fixed by nuclear physics**

The relative ease of D–T ignition compared to D–D, D–He3, or p–B11 is primarily a nuclear cross-section and reaction-rate fact. The toy model cannot change that. Therefore, we treat alternative fuels as scenarios that increase the importance of stability, confinement, and loss-channel management.

## **8.2 Where control and topology management can reduce tritium dependence**

Even within D–T, improved reliability can reduce tritium throughput and inventory requirements by increasing availability and reducing unplanned dumps. For advanced fuels, improved control can widen the viable operating window by reducing disruptive bursts and shaping where energy is deposited.

## **8.3 Specific research questions for fuels**

* Can reconnection control reduce transient heat loads enough to relax divertor constraints (a prerequisite for any advanced-fuel plant)?  
* Can control strategies bias energy partition away from wall-damaging channels (e.g., toward distributed waves/radiation analogs)?  
* Which control knobs matter most for p–B11 concepts where radiative losses are dominant—e.g., avoiding electron overheating during bursts?  
* How does event-rate reduction translate into plant-level fuel-cycle feasibility metrics (availability × duty cycle)?

# **9\. Milestones, resources, and collaboration needs**

## **9.1 Suggested milestone plan**

| Milestone | Outcome |
| :---- | :---- |
| M1 (0–6 weeks) | WP1 complete \+ experiment runner skeleton; first pulse-gated reconnection phase map for one IC. |
| M2 (6–12 weeks) | WP2 expanded to 2–3 ICs; stable outcome classifier; first invariant proxy (WP3) integrated. |
| M3 (3–6 months) | WP4 sheet onset curve; WP5 translation table and 3–5 testable mitigation hypotheses. |
| M4 (6–12 months) | WP6 tritium/advanced-fuel margin note; WP7 analog-platform bridging memo; draft paper outline. |

## **9.2 Resources and skills**

* Compute: GPU-capable runs for parameter sweeps; storage for saved checkpoints.  
* Skills: numerical PDE/spectral methods; basic MHD/reconnection literacy; experimental analog reading.  
* External collaboration: plasma physicist feedback for WP5 translation; possibly soft-matter experimentalists for WP7.

# **10\. Risk register and guardrails**

| Risk | Mitigation |
| :---- | :---- |
| Model mismatch | Toy-model dynamics may omit key plasma terms (two-fluid, kinetic effects). Mitigation: focus on control logic transfer; validate qualitatively. |
| Numerical artifacts | Topology change is sensitive to resolution and dissipation. Mitigation: control runs, convergence checks, standardized diagnostics. |
| Overclaiming | Risk of implying new energy sources or reactor designs. Mitigation: explicit scope guardrails; quantify only what simulations support. |
| Translation ambiguity | Mapping diagnostics to tokamak observables may be non-unique. Mitigation: build multiple candidate mappings; consult domain experts. |
| Complexity creep | Too many one-off scripts and untracked parameters. Mitigation: WP8 standardization; freeze run schemas. |

# **11\. References (fusion constraints and reconnection context)**

12. ITER Organization, “Disruption mitigation” (webpage). [link](https://www.iter.org/machine/supporting-systems/disruption-mitigation)  
13. ITER Organization, “Final design review is a major step forward” (disruption mitigation system, shattered pellet injection). [link](https://www.iter.org/node/20687/final-design-review-major-step-forward)  
14. M. Lehnen et al., “Disruptions in ITER and strategies for their control and mitigation,” Journal of Nuclear Materials (2015). [link](https://www.sciencedirect.com/science/article/abs/pii/S0022311514007594)  
15. R. A. Pitts et al., “Physics basis for the first ITER tungsten divertor,” Nuclear Materials and Energy (2019). [link](https://www.sciencedirect.com/science/article/pii/S2352179119300237)  
16. CEA/IRFM, “Innovative design for plasma-facing divertor components” (heat-flux figures, 2025). [link](https://irfm.cea.fr/en/2025/10/an-innovative-design-for-plasma-facing-divertor-components/)  
17. S. Zhen/Zheng et al., “Study of impacts on tritium breeding ratio of a fusion DEMO …” (criterion TBR \> 1.1, preprint/PDF). [link](https://scientific-publications.ukaea.uk/wp-content/uploads/Preprints/pre-CCFE-PR1732.pdf)  
18. I. E. Ochs et al., “Improving the feasibility of economical proton–boron-11 fusion …,” Physical Review E (2022). [link](https://link.aps.org/doi/10.1103/PhysRevE.106.055215)  
19. V. Sizyuk et al., “Comprehensive analysis of disruption mitigation methods …” Scientific Reports (2025). [link](https://www.nature.com/articles/s41598-025-14407-z)
