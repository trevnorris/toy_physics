# Does returning sideways momentum to the rim remove towing drag?

**Working result.** An elastic contact can return all the sideways momentum it acquired. If the rim takes that entire reaction, the contact also loses its steady pull on a held mouth. Sharing the reaction with independently supported material separates the two forces. Simultaneous impulse and work accounting constrains how far they can separate. [Numerical, symbolic and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_release_reaction.json)

### Fixed-intake comparison

We kept the [previous contact preparation](mouth_registering_rim.md). The mean time to the decision to leave is 128 feed residence times. The mean time to angular locking is 256. A departing parcel ordinarily takes its acquired sideways velocity with it. We tested a reversible alternative motivated by unloading an elastic contact. At the decision it stops exerting angular torque, but its center spring can finish the oscillation it started on entry. Switching anchors preserves stretch and parcel velocity. A fraction \(h\) of the unloading reaction goes to surrounding material moving at \(W\), and the rest goes to the rim moving at \(V\). The spring endpoint therefore moves at \((1-h)V+hW\). This ideal linkage and the release timing are hypotheses. Its extra residence increases material inventory at the same steady throughput. The angular preparation stays fixed by assuming no crowding or feedback from unloading.

The two coefficients below are forces per participating mass rate and per unit sideways speed. We held the rim in an incoming stream at speed \(U\), then towed it in resting material. The second support was stationary in both trials. Opposite drives and half-strength drives checked the signed linear response.

| Release at the rim | Held pull \(c\) | Towing drag \(g\) |
|---|---:|---:|
| Immediate departure | 0.999939 | 0.999939 |
| Locked contacts unload fully onto the rim | 0.666575 | 0.666575 |
| Locked contacts share the reaction equally | 0.833257 | 0.749916 |
| Every contact unloads fully onto the rim | 0 | 0 |

One third of departures were angularly locked. When every contact unloaded onto the rim, delaying final release by another \(1/12\) of an oscillation raised both zero coefficients to 0.133975. With the rim and support both stationary, different assignments of the unloading reaction gave the *same outgoing parcel momentum* and different rim forces. The total spring force includes the reaction on each support. Identifying the receiving support with part of the freely moving mouth would require including both reactions in that body's balance.

### What the energy balance requires at 0PN

For a stationary second support the actuator force is \(\dot M(gV-cU)\). If the only incoming sideways power is parcel kinetic energy, the power available for outgoing motion and dissipation is

\[
\frac{P_{\rm available}}{\dot M}=\frac{U^2}{2}+gV^2-cUV
=\frac{(U-cV)^2}{2}+\left(g-\frac{c^2}{2}\right)V^2 .
\]

Steady passive operation for all weak \(U,V\) thus requires \(c^2\le 2g\) under these assumptions. For the ideal benchmark \(c=1,\ g=0\), at \(V=2U\) even zero outgoing energy requires an additional \(3\dot M U^2/2\). A draining system can also exchange pressure work and inward-flow energy. The bound specifies what those channels would have to supply in this preparation. The [earlier ideal exterior](mouth_tow_and_companion.md) permits finite pull and zero towing drag. Its energy and stresses extend beyond this contact.

**Next test.** At fixed intake independently reverse the incident stream and rim speed. The part of the actuator power odd under both reversals is \(-\dot M cUV\). Outgoing energy and dissipation can also have a signed part. Close the balance with incoming and outgoing energy, storage, pressure work and each powered support's reaction and velocity at matched cuts. The shear-free bulk admits normal traction on a tilted surface. The five [archived histories](../../../docs/_measurements/receiving_force_motion_check.txt) supply feeding-side readings and prescribed rim motion. Release timing, receiving traction and whole-mouth acceleration remain open. Base 49dca23e. Independent local review is recorded in the measurement file. Formal review is pending.
