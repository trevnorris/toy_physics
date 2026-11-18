# 0PN gravity from two inward drains

**Result.** In an ideal three-dimensional brane exterior, one drain sends material toward itself. Where that flow enters a second drain, it delivers momentum toward the first. The change in carried momentum and pressure gives an attractive load proportional to \(1/d^2\). The physical model still owes the step from that load to the acceleration of a free mouth. The calculation below uses continuity and momentum balance. Newton's law enters only at the end as the comparison. [Held-drain derivation](held_drain_pair.md) · [Free-fall response test](../software/stage1_solver/reports/free_fall_response_ratio.md)

### 1. From one open mouth to the incident flow

Assume both mouths stay open and drain inward. At fixed ambient conditions a physical mouth of radius \(R\) must select its inward material mass rate \(D(R)\). The existing exterior solves prescribe the intake. In the ideal benchmark take a steady, spherical exterior with mass density \(\rho\) and no material exchange between concentric spheres. The same source mass rate \(D_s\) crosses every sphere. At distance \(d\) its inward speed is therefore

\[
U_s(d)=\frac{D_s}{4\pi\rho d^2}.
\]

The numerical carrier test fixed base density at the distant reservoir and solved for the interior field. At radii 2, 3, 4, 6 and 8 the inward number current on closed grid surfaces rounded to the prescribed \(0.25\). This checks conservation in that solver. Relating the physical drain and any distributed replenishment to its far-field material flux remains a separate requirement. [One-drain and two-mouth measurements](drain_driven_initial_acceleration.md#one-drain-current-and-two-mouth-force-with-a-distant-lateral-reservoir) · [Physical flux conditions](../research/pde_ledger_v3/steps/S14a_2_slow_limit.md)

### 2. From incident flow to a load toward the companion

Compare the same target alone and with its companion. Enclose the target in a sphere and integrate the steady momentum current \(\rho\mathbf v\mathbf v+p\mathbf 1\), including the Bernoulli pressure. In the constant-density potential exterior the extra inward momentum per second toward the companion is

\[
\Delta I_x=D_t U_s=\frac{D_tD_s}{4\pi\rho d^2}.
\]

On that enclosing sphere, carriage gives \(4/3\) of the answer and pressure gives \(-1/3\). Doubling the separation from 10 to 20 reduced the net by exactly four in the analytic benchmark. The finite feeding-slab solver separately integrated carriage and stress over the three-dimensional mouth and found a companion-directed load. Attributing the entire exterior load to the mouth requires knowing what is transmitted onward. The S15 shell record establishes a conditional momentum crossing; its item-4 independent classification remains deferred. [Sphere calculation and checks](held_drain_pair.md#alone-and-with-the-companion) · [Finite-mouth integral](drain_driven_initial_acceleration.md#paired-minus-alone-load-across-the-3d-mouth-at-0pn) · [S15 scope](../research/pde_ledger_v3/steps/S15_two_drain_region.md)

### 3. From a held mouth to a free mouth

Momentum arriving on the brane side can accelerate bound material, pass into the receiving bulk, or be transmitted by stress. The full budget includes the mouth, the slab faces, and the change in stored momentum. The fixed-center runs measure the incoming contribution. The five matched histories for the radius-one drain at inward number rate \(0.25\) contain a commanded rim motion and feeding-side readings. They supply no independent receiving-side force measured together with local motion. An earlier release example assigned a solid mass to a permeable sphere; the physical throat's inertia must emerge from its own state. [Held-mouth balance](held_drain_pair.md#what-the-holder-would-measure) · [Signed-run inventory](_measurements/receiving_force_motion_check.txt) · [Physical mass comparison](drain_driven_initial_acceleration.md#physical-free-fall-comparison-at-0pn)

A trapped standing wave could carry momentum as the mouth moves. The last paper test held the radius, intake and steady pull \(F\) fixed. Two candidate laws let the same rim motion carry different amounts of field momentum. Both conserved momentum. With other co-moving field inertia \(M\), they gave slow accelerations \(F/(M+m)\) and \(F/(M+2m)\). Here \(m\) comes from integrating momentum in the responding field. SymPy and Mathematica agreed. The example shows why knowing the steady pull alone cannot select a physical slow acceleration. Sudden release can also excite internal motion. [Derivation and assumptions](../software/stage1_solver/reports/free_fall_response_ratio.md) · [Verification record](_measurements/free_fall_response_ratio.json)

### 4. What a Newtonian result still requires

Let \(C_t\) be the total pull on a stationary target per unit weak incident speed, measured by the opposite holding reaction. Let \(M_t\) be its independently determined slow translational inertia in resting brane. If the incident speed supplies the whole leading drive, other loads are accounted for, and the scalar slow response applies, the predicted companion-directed acceleration is \(a_t=(C_t/M_t)U_s(d)\). Equal fall toward the same source requires the same positive \(C_t/M_t\) across physical target radii. Reproducing the source-mass dependence also requires far-field source strength proportional to independently determined source inertia. Neither comparison has been made for freely responding physical mouths. [Mass and radius requirements](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#how-drain-rate-and-retained-momentum-must-scale-for-free-fall)

Energy supplies a separate check. In the ideal pair benchmark, holding one mouth still and prescribing approach speed \(s\) gives exterior interaction-energy growth \(Fs\) and mechanical work \(Fs\), where \(F\) is that benchmark's attractive load magnitude. The inner regions together owe incremental power \(2Fs\) relative to the isolated controls. The physical supply remains open. [Stationary-mouth test and review scope](_measurements/stationary_mouth_energy_review/round2/status.md)

**Where to resume.** A mouth-level law or an independent measurement must determine the pull transferred to the translating mouth and the inertia moving with it on the same inward branch. Compare the resulting acceleration across radii. An orbit calculation would then have an independent force law to test. [Release protocol](drain_driven_initial_acceleration.md#physical-free-fall-comparison-at-0pn)

**Research record.** Earlier local commits cover the unequal-mouth budget (`806ae7b6`), held-mouth accounting (`4b49fe61`), release requirements (`9c526a84`), one-mouth response (`fd667ea9`), and the stationary-mouth energy benchmark (`cda9a8eb`). The final response-ratio test and this synthesis are working results pending independent review. [Detailed calculation history](drain_driven_initial_acceleration.md)
