# Can a bowed conversion surface keep a draining mouth inward at 0PN?

**Continuing the flat feeding field into a bowed opening produces outward crossing near the rim, even for one isolated drain.** The total intake stays inward. Carrying the current and density with the bend preserves inward crossing but fails the carrier phase check. Neither construction measures the force on a free mouth. [Guarded runs and independent reviews](../../../research/pde_ledger_v3/_measurements/mouth_responding_rim.json)

**Shape test.** Keep both radii at 1 and each inward number rate at 0.25. For a stationary graph \(w=h(\mathbf x)\) attached at its edge, steady source-free transport through the layer gives the number current per projected volume

\[
q_h=j_w(0)-\nabla_\parallel\!\cdot\int_0^h\mathbf j_\parallel\,dw .
\]

The first trial carries the measured in-plane current unchanged across that depth and lets the shape solve \(-K\nabla^2h=q_h\). The uniform conversion coefficient \(K\) could be calibrated from an isolated mouth’s intake and mean rim slope. We tested three trial values. A separate flat order-energy calculation gives \(\sigma=1/3\) for the trial material. Neither coefficient has been measured for the archived mouth. The reconstructed aperture current was normalized to its recorded intake. [Shape solver](../tools/mouth_advected_rim.py) · [Earlier tension and slope model](mouth_front_pressure_dipole.md)

On the fine feeding grid at \(K=0.4\), the companion-directed rim-slope dipole changed by **\(1.469\times10^{-5}\)** at \(8R\) and **\(3.660\times10^{-6}\)** at \(16R\). At the edge the current crossing the surface is \(j_w+s j_r\), where \(s=-\partial_r h>0\). The archived downward current tapers to zero while the brane flows radially inward, \(j_r<0\), around the entire rim. This sign test only requires continuity at the edge. In the isolated trial, crossing reversed over the outer 12–13% of the radius and the returned amount was **1.05%** of the net inward rate. It fell to **0.125%** at \(K=1.6\). The band widths depend on the trial law and interpolated feed. [Interior check](../tools/mouth_rim_reverse_extent.py) · [Independent edge-current check](../tools/mouth_rim_link_current.py)

**Bending the flow with the mouth.** Carrying density and current with the bend in the fixed laboratory geometry adds \(\mathbf j_\parallel\cdot\nabla h\) to its downstream component and preserves inward crossing pointwise. Its predicted dipole changes at \(K=0.4\) are \(1.250\times10^{-7}\) and \(3.122\times10^{-8}\). A smooth neutral carrier requires zero circulation around a contractible loop with nonzero density. On one loop the calculated circulations at \(8R\) and \(16R\) were **\(+2.232\times10^{-6}\)** and **\(+5.568\times10^{-7}\)**, versus \(-2.65\times10^{-18}\) for the isolated mouth. A second loop and a second velocity reconstruction confirmed the paired signal. This particular deformation violates the feeding phase law. [Material map](../tools/mouth_following_front.py) · [Circulation test](../tools/mouth_lifted_front_circulation.py)

**Independent bulk pressure.** Two outlets discharge into one initially quiescent, constant-density, irrotational 4D half-space. Supplying both reconstructed local current profiles gives receiving-side pressure moments of \(4.136\times10^{-6}\) and \(5.051\times10^{-7}\), respectively **82.43%** and **40.25%** of the corresponding feeding moments. Under full transfer of those flat-face pressures, \(\sigma s_1\) would be \(2.104\times10^{-7}\) and \(1.790\times10^{-7}\). Both pressures respond. This matches the number flux on the flat opening. The actual phase and pressure on the bowed surface remain open. [Redistributed discharge](../tools/mouth_open_bulk_archived_current.py) · [Radial control](../tools/mouth_open_bulk_pressure_moment.py) · [Prior phase-matching failure](one_mouth_outer_holding_force.md)

**What the next physical solve must do.** Let the carrier phase adjust on the same curved surface on which pressure acts. For surface phase \(\theta_s\), number density \(n_s\) and projected crossing \(q_h\), the required derivative into the bulk is

\[
\partial_w\phi=\frac{q_h/n_s+\nabla h\cdot\nabla\theta_s}{1+|\nabla h|^2}.
\]

This relation allows tangency at the lip and positive crossing inside. A local phase expansion also satisfies stationary number conservation at the surface for a trial phase. This allows the feeding trace to change. For the measured inward radial current, tangency requires a positive downstream component at the lip, \(j_w=-s j_r\). The bulk and feed equations must select compatible phase and pressure there. Test the isolated mouth first, then use the same law at both separations. Working base `49dca23e`, formal review pending. [Local phase calculation](../tools/mouth_rim_phase_trace.py)
