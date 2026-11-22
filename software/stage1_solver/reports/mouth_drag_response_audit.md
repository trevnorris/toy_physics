# Does the draining mouth experience drag at low speed?

**Step 1, 22 November 2026.** The saved rigid-mouth runs separate the lag caused by driving the incident stream from the resistance inferred by pushing the mouth through resting material. The resistance depends on where force is measured. The out-of-phase part of the reconstructed stream force grows approximately in proportion to frequency across the sampled band. The sign and size of steady physical drag remain open. All entries below use model units. [Complete complex responses, source hashes and guarded command](../../../research/pde_ledger_v3/_measurements/mouth_drag_response_audit.json)

The subsequent multi-box study completed the independent energy and momentum comparison on a coherent carrier. [Step 1 frequency table and finite-box result](mouth_rigid_drag_table.md)

### What the two drives measure

The same radius-1 cubic mouth translates rigidly in finite feeding and receiving reservoirs. In the stream drive a phase gradient oscillates on the distant brane while the far bulk rests. In the push drive both baths rest and a known force acts on the center. With the solver’s $e^{+i\omega t}$ convention, $a=i\omega V$. Write $\chi_a=a/F$ for the external push response and $H$ for the force a holder would exert on the stationary mouth in the stream. Eliminating the common linear response gives

$$
\gamma_{\rm app}=\Re(F/V)=-\omega\Im(1/\chi_a),\qquad
S=-H/U=\frac{a/U}{a/F}.
$$

The second equality for $S$ reconstructs the total stream force per unit incident speed. Identifying a Froude–Krylov contribution separately would require the force from the undisturbed incident pressure integrated over a specified surface. [Drive and force conventions](../tools/mouth_carrier_free_center.py)

| Mesh | Intake | $\omega$ | Push $\gamma_{\rm app}$ | Free stream $\Im(a/U)$ | Held load $\Re S$ | Held load $\Im S/\omega$ |
|---:|---:|---:|---:|---:|---:|---:|
| 26 | 0.25 | 0.01 | +0.0005017514 | −0.000207991 | 0.2502100 | −0.3573734 |
| 26 | 0.25 | 0.02 | +0.0005017430 | −0.000542575 | 0.2508019 | −0.3582144 |
| 26 | 0.25 | 0.04 | +0.0005017092 | −0.001158659 | 0.2531869 | −0.3616035 |
| 34 | 0.25 | 0.01 | +0.0005421093 | −0.000201243 | 0.2502081 | −0.3573725 |
| 34 | 0.25 | 0.02 | +0.0005421025 | −0.000539087 | 0.2508000 | −0.3582135 |
| 18 | 0 | 0.01 | 0 | −0.000310550 | 0 | −0.3813427 |
| 18 | 0 | 0.02 | 0 | −0.000622535 | 0 | −0.3822382 |

The radial and depth extents of these baths are 9 and 10. On mesh 34 at 0.01, $\Im(a/U)/\Re(a/U)=-0.00986085$. The stream load itself has $\Im S/\Re S\simeq-0.0142830$, while the push contributes the opposite phase ratio $+0.00442153$. Since $a/U=S\chi_a$, the phases add. The apparent one-model-time coefficient from the free response mixes these contributions. The signed time coefficient of the reconstructed stream force, $\Im S/(\omega\Re S)$, is −1.42830 model times on this preparation. The zero-intake mouth retains a stream lag with zero measured push resistance. A constant $\gamma$ would appear as a $1/\omega$ term in inverse push susceptibility. Expanding the phase of $a/U$ around an inertial response requires $|\gamma|/(M\omega)\ll1$. [Frequency plot](figures/mouth_drag_response_audit.png)

### Where the force measurement disagrees

The table gives the external force required per unit steady towing speed, in units of $10^{-6}$. Its sign is the reverse of the archive’s fluid force on the moving wall. The two interior readings integrate the momentum balance using different spatial windows around the mouth. The outer cut is a third momentum estimate. Both remote media are at rest in the laboratory frame. [Steady force integrals](../tools/mouth_carrier_lateral_inertia.py)

| Mesh | Bath radius / depth | Intake | Direct wall | Interior 1 | Interior 2 | Outer cut |
|---:|:---:|---:|---:|---:|---:|---:|
| 26 | 9 / 10 | 0.25 | +501.68 | −11.72 | −10.95 | +201.00 |
| 34 | 9 / 10 | 0.25 | +542.08 | −6.07 | −8.36 | +118.61 |
| 48 | 9 / 10 | 0.25 | +521.24 | −6.91 | −7.05 | +58.13 |
| 28 | 14 / 15 | 0.25 | −1515.95 | −1.20 | −4.28 | +221.68 |
| 36 | 20 / 21 | 0.25 | −2472.19 | +6.05 | −2.84 | +148.95 |
| 22, 26 | 9 / 10 | 0 | 0 | 0 | 0 | 0 |

At mesh 34 the steady wall resistance nearly equals the pushed wall’s harmonic resistance. They share a wall-force discretization. The momentum estimates disagree even on that same mesh. The enlarged rigid baths change mesh as well as outer boundary, so their differences cannot isolate a box-size effect or establish an error bound.

The rigid archive has no harmonic drives in an enlarged bath, no mesh-34 harmonic result at 0.04, and no rigid results at 0.005 or 0.002.

The occupied-mouth control uses two free bends, a trapped wave of unit action and two assumed tensions. Across 48 entrance-push responses (four mesh/bath preparations, two tensions and six frequencies from 0.00015625 to 0.005), $|\Re(F/V)|\le 9.58\times10^{-15}$. A Hermitian finite-bath operator forces a real collocated susceptibility away from poles at real drive frequency. The observed resistance is at roundoff. The separate occupied flux checks measure momentum. The archive lacks an independent cycle-averaged energy balance for the draining rigid push. [Occupied controls](mouth_occupied_slow_push.md)

**Next measurement.** Use the mesh-34, $\omega=0.01$ rigid push to compare cycle-mean actuator power $\Re(FV^*)/2$ with independent energy currents in the brane and bulk. Record the matched undriven flow. Average opposite drives before subtracting that baseline to retain the quadratic work. Halving the amplitude should quarter the leading excess. Retain the change in mean throughflow and energy supply at second order in the drive, or derive a balance that eliminates it. On each of two enclosing surfaces include the rate of change of stored energy, net outward power, and the work supplied by reservoirs or fixtures, counting each port once. The cycle-mean storage rate vanishes once periodicity is established. For lab-frame energy density $e$ and energy current $\mathbf J_E$, the outward flux on a moving cut is $(\mathbf J_E-e\,\mathbf v_{\rm cut})\cdot\mathbf n$. Repeat with mesh and outer boundary controlled separately. This can tell whether the wall resistance corresponds to net work, where the energy goes, and whether the numerical error can resolve it. A bound on steady drag still requires steady towing and the open-boundary, low-frequency limit.

Two independent local checks covered the signed reduction of the rigid drives and the steady and occupied controls, respectively. Each located the same missing energy measurement. All response data came from saved runs. [Archive reducer](../tools/mouth_drag_response_audit.py)
