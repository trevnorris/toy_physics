(* ================================================================= *)
(* HAMILTONIAN DERIVATION OF STABLE 1PN EQUATIONS *)
(* ================================================================= *)
(* Goal: Derive the acceleration from a Hamiltonian H(r, p) to guarantee *)
(* energy conservation and stability (no anti-damping). *)

Print["--- DERIVING HAMILTONIAN EQUATIONS OF MOTION ---"];

ClearAll["Global`*"];

(* 1. Define Variables *)
(* r and p are vectors *)
rVec = {x[t], y[t], z[t]};
pVec = {px[t], py[t], pz[t]};
rMag = Sqrt[rVec . rVec];
pMagSq = pVec . pVec;

(* 2. Define the Hamiltonian H(r, p) *)
(* Standard Kinetic: p^2 / (2m) *)
(* Potential: -GmM/r - chi/r^2 *)
(* Interaction: We apply the superfluid correction (1 + alpha p^2) to the Potential *)
(* Note: This couples position and momentum, generating the 'lag' effect safely *)

HKinetic = pMagSq / (2 * m);
HPotential = -(G * M * m / rMag) - (chi / rMag^2);
HInteraction = (1 + alpha * pMagSq / (m^2 * c^2));

H = HKinetic + (HPotential * HInteraction);

(* 3. Hamilton's Equations *)
(* dr/dt = dH/dp *)
(* dp/dt = -dH/dr *)

Print["Calculating Hamilton's Equations..."];

(* Velocity Equation (v = dr/dt) *)
vx = D[H, px[t]];
vy = D[H, py[t]];
vz = D[H, pz[t]];

(* Force Equation (F = dp/dt) *)
Fpx = -D[H, x[t]];
Fpy = -D[H, y[t]];
Fpz = -D[H, z[t]];

(* 4. Solve for Acceleration *)
(* We have v as a function of p and r. We need a = dv/dt. *)
(* a = d/dt (dH/dp) *)
(* This involves dp/dt, which we know is -dH/dr. *)

(* Compute dvx/dt using the Chain Rule *)
(* dvx/dt = sum( d(vx)/d(coord) * d(coord)/dt ) *)
(* Variables are x, y, z, px, py, pz *)
(* Note: d(coord)/dt are vx, vy, vz for positions, and Fpx, Fpy, Fpz for momenta *)

axExact =
  D[vx, x[t]] * vx + D[vx, y[t]] * vy + D[vx, z[t]] * vz +
  D[vx, px[t]] * Fpx + D[vx, py[t]] * Fpy + D[vx, pz[t]] * Fpz;

(* 5. Series Expansion *)
(* We assume p ~ m*v for the small relativistic corrections to make the output readable. *)
(* This simplifies the momentum terms back to velocity terms. *)

Print["\n--- EXPANDING ACCELERATION TO 1PN ORDER ---"];

(* Replace momentum p with m*v for the expansion (valid to 1st order in 1/c^2 for the correction) *)
(* We use a symbolic vector vVec for readability in the output *)
readableAx = axExact /. {
   px[t] -> m * vxSym, py[t] -> m * vySym, pz[t] -> m * vzSym,
   x[t] -> x, y[t] -> y, z[t] -> z
};

(* Expand *)
axSeries = Series[readableAx, {c, Infinity, 2}];
axSimplified = Simplify[Normal[axSeries]];

Print["Acceleration x-component (a_x):"];
Print[axSimplified];

Print["\n--- CHECKING FOR DRAG TERMS ---"];
(* We check if there is a term proportional to velocity v *)
CoeffV = Coefficient[axSimplified, vxSym];
Print["Coefficient of Velocity (Drag Check):"];
Print[Simplify[CoeffV]];

Print["\n--- DIAGNOSIS ---"];
Print["If the Hamiltonian is correct, this acceleration derives from a conserved energy."];
Print["Even if velocity terms appear, they must be conservative (symplectic)."];
