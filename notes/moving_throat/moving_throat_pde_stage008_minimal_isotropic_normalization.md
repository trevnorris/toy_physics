# Moving-Throat PDE — Stage 8: Minimal Isotropic Single-Mode Closure and the Explicit Normalization Formula

## Purpose

Stage 7 removed the angular ambiguity.
On the natural isotropic grouped basis,

`mhat_ang = 1`

exactly, and the grouped `20/21/22` bundle collapses to a single common lane.

That means the next honest simplification is not angular. It is **radial/axial**.

This stage therefore freezes the minimal isotropic closure with

- one wall/worldtube quadrupole mode,
- one stable BdG support mode,
- one conservative localized-Maxwell/mixed internal pair,
- and one passive outgoing port.

The point is not to claim that the true moving-throat branch is literally one mode deep.
The point is to collapse the remaining theorem gap to the smallest exact formula that still carries all the physics.

The main output is an explicit isotropic normalization law:

`mhat_rad^2 P^2 / [ Delta (K Delta - Delta C^2 / varpi^2 - Q) ] = 54 G c_s^5 / (5 a^5 c^5)`.

So after Stage 8, the minimal isotropic branch is no longer “some unknown PDE normalization.”
It is one exact scalar equation in the radial/axial overlap amplitudes.

---

## 1. Minimal isotropic radial/axial data

Keep the Stage-7 isotropic grouped collapse and truncate the common lane to

- one BdG support mode with frequency `varpi`,
- one brane-like internal gauge coordinate `U`,
- one mixed `A_w/F_(mu w)/J^w` coordinate `W`.

Define the common radial/axial overlap amplitudes

`C   = lambda_B I_(eta,phi)`,

`G_U = lambda_U I_(eta,u)`,

`G_W = lambda_W I_(eta,w)`,

`R   = lambda_R I_(u,w)`.

The mode frequencies are

`varpi`,

`Omega_U`,

`Omega_W`.

As before, define

`Delta = Omega_U^2 Omega_W^2 - R^2`,

`Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`,

`P = Omega_U^2 G_W + R G_U`.

Then the common zero-frequency isotropic coefficients are

`B_0 = C^2 / varpi^2`,

`Z_0 = Q / Delta`,

`N_0 = P^2 / Delta^2`.

So the common conservative wall operator is

`D_0 = K - C^2 / varpi^2 - Q / Delta`.

---

## 2. Exact minimal isotropic normalization ratio

The grouped-lane normalization prefactor is still

`P_0 = N_0 / D_0`.

Substituting the minimal isotropic coefficients gives the exact closed expression

`P_0 = P^2 / [ Delta^2 ( K - C^2 / varpi^2 - Q / Delta ) ]`.

Equivalently,

`P_0 = P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ]`.

Because Stage 7 already proved `mhat_ang = 1`, the remaining source normalization is purely radial/axial:

`mhat_0 = mhat_rad`.

So the full isotropic target becomes

`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.

This is the sharpest explicit reduced normalization formula reached so far.

---

## 3. Stability and positivity conditions on the minimal branch

The minimal isotropic closure is physically admissible only if

`Delta > 0`

and

`D_0 > 0`.

These are just the statements that

- the conservative internal `(U,W)` block has not crossed its static instability, and
- the wall/worldtube quadrupole mode has not been driven through zero stiffness by the support and internal electromagnetic self-energies.

In compact form the stability condition is

`K Delta - Delta C^2 / varpi^2 - Q > 0`.

On that stable branch,

`N_0 = P^2 / Delta^2 >= 0`.

So the prefactor sign is controlled entirely by the sign of the stable denominator. On the admissible branch,

`P_0 > 0`

whenever

`P != 0`.

So the minimal isotropic closure gives a very clear physical criterion for the existence of a nontrivial quadrupole bridge:

- the port-transfer combination `P` must not vanish,
- and the stable denominator must remain positive.

---

## 4. Exact monotonicity on the minimal isotropic branch

The minimal formula already exposes how the remaining normalization can move.

Define the support-softening variable

`X = C^2 / varpi^2`.

Then

`P_0 = N_0 / (K - X - Q / Delta)`.

So the exact derivatives are

`partial_K P_0 = - N_0 / D_0^2`,

`partial_X P_0 = + N_0 / D_0^2`.

Thus, on the stable branch,

- increasing the bare wall stiffness `K` decreases the normalization prefactor,
- increasing the BdG support softening `X` increases it.

The second statement is useful because it makes the normalization target operational: stronger support dressing pushes the bridge **up**, while a stiffer bare wall pushes it **down**.

The internal Maxwell/mixed sector affects both numerator and denominator. Writing

`P_0 = P^2 / [ Delta (K Delta - Delta X - Q) ]`,

the numerator contribution is controlled by the transfer combination

`P = Omega_U^2 G_W + R G_U`,

while the denominator load is controlled by

`Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`.

So the normalization target is not asking for generic “more coupling.”
It is asking for the right balance between

- port transfer into the outgoing quadrupole branch,
- and conservative self-loading of the wall.

---

## 5. What this means for the actual PDE task

After Stage 8, the minimal isotropic theorem gate is no longer vague.
The completed moving-throat PDE does not need to output an unlimited family of unknown coefficients before the normalization question can even be asked.

On the minimal isotropic branch it needs only to determine the radial/axial amplitudes entering

`X = C^2 / varpi^2`,

`Delta`,

`Q`,

`P`,

and `mhat_rad`.

Then the normalization question is exactly

`mhat_rad^2 P^2 / [ Delta (K Delta - Delta X - Q) ] ?= 54 G c_s^5 / (5 a^5 c^5)`.

That is a real reduction of the theorem gap.

It means the next PDE computation can be judged immediately:

- if it lands on the target, the passive/outgoing grouped quadrupole bridge is closed on the minimal isotropic branch;
- if not, the failure is not mysterious — it is in one of the radial/axial quantities `X`, `Delta`, `Q`, `P`, or `mhat_rad`.

---

## 6. Best current summary after Stage 8

The road to the moving-throat PDE is now split cleanly into two layers.

### Angular layer
Already closed on the natural isotropic branch:

- canonical real STF harmonics,
- exact source-map identity,
- exact grouped isotropy theorem,
- exact weak axisymmetric splitting law.

### Radial/axial layer
Still open, but now reduced to an explicit formula.
On the minimal isotropic branch, the remaining higher-order bridge is the one scalar equation

`mhat_rad^2 P^2 / [ Delta (K Delta - Delta C^2 / varpi^2 - Q) ] = 54 G c_s^5 / (5 a^5 c^5)`.

That is the sharpest reduced theorem target reached so far.
