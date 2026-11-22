# How a sideways push changes mean particle throughput

At radius 1, baseline intake 0.25 and angular frequency 0.01, the mesh-34 rigid mouth acquired an additional **0.0315034140 units of feed-to-bulk particle current per unit squared peak sideways speed**. Both distant reservoir phases were fixed on the translating chart. The result is the second-order coefficient about this draining state. It allows the separate reservoir energies to be converted to a common chemical reference while preserving their sum. [Guarded readings and source hash](../../../research/pde_ledger_v3/_measurements/mouth_push_mean_throughput.json)

Write the time-dependent carrier as \(f+A\mu(c\cos\omega t+s\sin\omega t)+A^2k+\cdots\), where \(A\) is peak speed, \(\mu=x/r\), and \(k\) is the angularly uniform, time-independent part of the second-order correction. The omitted angular and temporal harmonics of second order integrate to zero in the mean throughput at this order. The base field \(f\) and both velocity-driven dipoles \(c,s\) come from the reference rigid run. At the two fixed complex reservoir traces \(k=0\). Regularity on the symmetry axis and the moving-wall condition complete the mean problem. [Harmonic fields](../tools/mouth_rigid_push_energy.py)

Let \(n=|f|^2\) and \(v_h=2\operatorname{Re}(f^*h)\). Expanding the carrier interaction \( (|\psi|^8-1)\psi/4 \) gives the forcing

\[
Q=\tfrac16\sum_{h=c,s}\left[n^3(v_hh+|h|^2f)+\tfrac32 n^2v_h^2 f\right].
\]

The factor \(1/6\) is the product of the angular average of \(\mu^2\) and the time average of one squared sinusoid. With axisymmetric real test function \(u\), the static real-linear derivative \(L_0\) of the carrier equation gives

\[
\langle u,L_0 k\rangle =
\frac{i}{12}\int\!\left[(\partial_r u)c-u(\partial_r c+2c/r)\right]dV
-\int uQ\,dV .
\]

The first term is the projection of the moving-chart transport and reflecting-wall condition. The free-node residual in this solve was \(2.18\times10^{-15}\). Testing the residual on the fixed phase at each bath gives outward current coefficients of −0.0315034140250737 at the feed and +0.03150341402502623 at the bulk. Direct integration onto four horizontal sections of the conduit gives 0.0315174204, 0.0315038032, 0.0315033360 and 0.0315027527 at depths 0.7, 1.4, 2.4 and 3.0. Both the mean-field correction and the averaged dipole current contribute. With the reservoir size and frequency held fixed, the receiving current coefficient was 0.0315024732 on mesh 26 and 0.0315042243 on mesh 48. [Mean field and flux implementation](../tools/mouth_push_mean_throughput.py)

The canonical energy current used in the push audit subtracts a common chemical potential. Restoring a reference \(\mu_{\rm chem}\) adds \(\mu_{\rm chem}\) times these particle currents to the two outward energy readings. The combined correction vanishes to the accuracy of the conserved particle balance.

As a separate check of the energy observable, direct space and time integration on the dipole-plus-mean field gave feed power \(+3.22824608\times10^{-14}\) and bulk-side power \(-5.34570162\times10^{-14}\) at applied push \(10^{-5}\). Reversing the push gave the same values to the quoted precision. Half-strength readings were \(+8.07061520\times10^{-15}\) and \(-1.33642540\times10^{-14}\). These are powers in the shifted energy convention evaluated on a perturbative field. Independent fully nonlinear trajectories and the higher-order correction to the actuator impedance remain open. Direct port normal derivatives also carry discretization error. The conserved phase reactions and four interior sections above supply the stronger check on particle throughput.
