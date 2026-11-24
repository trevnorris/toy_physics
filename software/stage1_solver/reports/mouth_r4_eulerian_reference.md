# Radius 4 — the single Euclidean reference map fails at the bonded ring

**Stage 1 stopped before descent.** Checkpoint 10 and the post-shape material state in `e2e1d048` were tested at the held entrance radius 4, one quantum, fixed phase drop 0.4510084282857874 and \(U(n)=(n^5-5n+4)/20\). The core was registered stress-free in the shallow cavity. The collar retained its flat-brane reference. They share their current placement at the bond.

Independent review concluded “Stage 1 fails the single-reference-map strain gate on both saved states.”

**The same ring has two natural lengths.** A continuous spatial reference map has one tangential derivative along a bonded ring. Even if its normal derivative jumps, the magnitude of its azimuthal derivative fixes one natural circumference. At material depth fraction 1/48 and polar angle 0.65449847 radians, the core assigns natural radius **3.5261816063** and the collar assigns **4.5** in both snapshots.

| Saved state | Order χ at bond | Core hoop stretch | Collar hoop stretch |
|---|---:|---:|---:|
| Checkpoint 10 | 0.18527171 | 0.99904009 | 0.78284373 |
| After the shape and crystal steps | 0.20317028 | 0.99862911 | 0.78252169 |

For a single common hoop stretch, the smallest possible maximum relative error is \(|h_c-h_s|/(h_c+h_s)=0.1213302217\). This is a local one-sided strain bound under fixed material correspondence. It applies at all three specified gauge strengths because the ring and its immediate neighborhood have positive order. [Trace calculation](../tools/mouth_inverse_bond_trace.py)

**Energy and orientation controls.** Writing \(G=\nabla X\), \(j=\det G\), the spatial energy density on each oriented material reference patch is \(\chi j W(G^{-1})\). We integrated it at the mapped native quadrature points using the full four-dimensional volume, including the azimuthal stretch. Summed elementwise differences from reference-volume energy were \(+1.11505\times10^{-16}\) and \(+1.03260\times10^{-16}\), respectively. Maximum absolute element errors were \(2.30393\times10^{-18}\) and \(3.17807\times10^{-18}\). Total energies rounded to 92185.53614224629 and 92185.47335775776. All eight held-field current/balance readings reproduced their saved values exactly. Reflecting one inverse-gradient row at an active quadrature point on each state gave \(j=-0.98970857\) and \(-1.03300085\). The physical energy evaluator rejected both. These are two-reference integration controls. [Element results and signed controls](../../../research/pde_ledger_v3/_measurements/mouth_r4_eulerian_reference.json)

**Six-force control.** Forces below use the checkpoint-10 shape coordinates and sign. The analytic column uses both registered references and holds normalized spatial order knots fixed. The last column is the independently calculated signed crystal work with the production knot motion at the original step. The first column includes all other field and surface loads.

| Coordinate | Saved total | Analytic crystal, fixed knots | Crystal, production motion |
|---|---:|---:|---:|
| Log depth | +0.150565 | −0.441280 | −0.022267 |
| Scaled log exit | −0.641236 | +0.144850 | +0.144850 |
| Mean width | +0.705499 | +0.539551 | +0.539551 |
| Depth asymmetry | −0.666834 | −0.284254 | −0.284254 |
| Waist profile | +0.465603 | +0.404060 | +0.404060 |
| Polar–equatorial shape | +0.009007 | −0.040334 | −0.040334 |

Across all six directions the largest analytic versus half-step difference **at fixed knots** was \(8.014\times10^{-9}\) at checkpoint 10 and \(8.789\times10^{-9}\) after the shape step. Transport of the order knots and quadrature contributes +0.41901314 and +0.32986399 to the respective half-step depth forces. The equivalent single-map total-force residual and a free-wall kinematic correction are undefined because a continuously bonded single-map baseline already changes the ordered strain.

**Spatial resolution and next test.** At the post-step tail, the saved \((\chi,J,W)=(8.29976\times10^{-21},2.11022\times10^{-14},31.09361)\) requires \(J_{\rm inv}\approx4.74\times10^{13}\) and \(\|\nabla u\|_2\ge3.6785\times10^{13}-1\). Local 128× trilinear interpolation on the physical grid gave \(J=0.01102371\) there, using saved labels at four melted vertices. A regular checkpoint-10 control gave 0.24348121 against 0.24346730. The trace bound holds for ε = 0, 10⁻¹⁰ and 10⁻⁸ with a gauge supported in melted material. Finite-gauge whole-mouth energies and forces were not evaluated. Checkpoint 11, free release and the blind radius-3 prediction remain pending. The next formulation to test carries one spatial field of material labels and evaluates strain with the original core and collar reference metrics on their respective sides of the moving bond. In this formulation the inverse tensor on side \(m\) is \(G_m=(\partial_\eta X_m)\nabla_x\eta\), with the corresponding physical metric factors. Both sides share the material chart \(\eta\) at the bond. The compressed-tail resolution test still applies. [Local spatial interpolation](../tools/mouth_eulerian_r4_spatial_cell.py)
