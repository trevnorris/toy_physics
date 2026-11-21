# Can a draining mouth retain an orientation while its material turns over?

**Working result.** A delay before a contact can exert torque is sufficient to transmit a prepared signed winding through the feeding layer in the seven-shell model. The same contacts export sideways momentum under a steady tow. This gives two independently testable consequences of one attachment law. [Angular field and preparations](mouth_collar_charge_transmission.md) · [Commands, outputs and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_registering_rim.json)

## What a contact remembers

Assume a contact attaches elastically on arrival, initially swivels, and locks to the local preferred orientation after an exponentially distributed time of mean \(T\). The prepared orientation must persist as individual contacts leave. A pair of opposite spring arms of lever length \(\ell\) and total stiffness \(k\) stores translational energy \(kx^2/2\). Once registered it also stores \(k\ell^2 v(1-\cos\Delta\theta)\), up to a constant, with \(v\) the prepared phase visibility. At arrival rate \(j\) per link and independent exponential residence time of mean \(\tau\), the mean angular stiffness before applying visibility is

\[
G=j k\ell^2\tau\,\frac{\tau}{\tau+T}.
\]

All comparison mouths have identical intake, geometry, center springs and residence times. We chose illustrative rim and feed residence times in the ratio 128 to 1 and set \(j k\ell^2\tau_f/K=2\), where \(K\) is the field stiffness. At \(T=256\tau_f\), one third of rim residents and one in 257 feed residents carry torque. The mirror mouths have exterior ground-state winding \(+1,-1\) in all eight archived finite-mesh comparisons. Immediate registration gives \(0,0\). So does a freely swiveling contact throughout its life. These are equilibrium predictions for the mean angular energy with the prepared reference. The physical origin and energy supply of registration, charge readout, radius dependence and lifetime remain open.

Discrete contacts can matter. For Poisson arrivals and fixed preferred phase and visibility per registered contact, we checked one frozen snapshot of the finer mesh statistically. At a mean of 16 feed residents per link, the sufficient bound is inconclusive. At 16,384, with weaker individual contacts chosen to preserve the same mean stiffness, the numerical upper bound on failing the sufficient transmission conditions is approximately \(3.37\times10^{-10}\). The bound describes one snapshot.

## Where the sideways momentum goes

For one representative patch with no separate center mass, incoming parcels of mass \(m\) enter unstretched at lateral velocity \(U\). Let \(V\) be the center velocity, \(u\) the mean resident velocity and \(\dot M=jm\) the participating mass rate. Momentum conservation gives the actuator force \(F_{\rm act}=d(j\tau m u)/dt+\dot M(u-U)\). For steady motion this yields

\[
F_{\rm act}=\Gamma(V-U),\qquad
\frac{\Gamma}{\dot M}=\frac{(k/m)\tau^2}{1+(k/m)\tau^2}.
\]

With \(\sqrt{k/m}\,\tau_f=1\), the rim coefficient is \(0.999938968569\). Opposite full and half pushes in resting material gave opposite velocities proportional to the applied force. A constant push approaches a terminal speed as departing material carries away momentum. Towing work balances outgoing kinetic energy and spring energy released at detachment. This is the ensemble response of one participating stream. Whole-mouth inertia and the reactions of surrounding brane and bulk remain open.

## Discriminating test

At fixed intake, label material entering the feed and rim, and measure torque per surviving contact versus age under opposite weak twists. Check whether the reference persists across replacements and what powers locking. Under opposite steady lateral motions measure both exit momentum and total actuator reaction. These readings distinguish export of the acquired impulse from transfer back to the rim or surrounding brane before release. The latter distinction directly constrains 0PN free fall. Base commit 49dca23e. Scoped local review is recorded with the measurements. Formal review is pending.
