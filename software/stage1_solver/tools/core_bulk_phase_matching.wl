(* Test the boundary terms by integration by parts on the two sides of a
   moving cut. Then use a frame following a tilted, moving order level to
   solve for the phase slope on that level. *)
ClearAll["Global`*"];
kin = hb^2/(2 m);
(* Coefficients of a test value f and its normal derivative fp in the
   distribution i hb psi_t + kin psi_zz. The side fields are regular. *)
weakBoundary = I hb speed (left - right) f +
  kin (left fp - leftSlope f - right fp + rightSlope f);
Print["smooth cut, coefficients of f and fp = ",
      InputForm[Expand[{Coefficient[weakBoundary, f], Coefficient[weakBoundary, fp]}]]];

(* chi(x,z,t)=F(z-b*x-c*t), carrier phase gradient (px,pz).
   The common surface moves with normal speed c/norm, pointing to +z.
   dF is its signed slope in the unnormalized coordinate. *)
norm = Sqrt[1 + b^2];
u = hb {px,pz}/m;
normal = {-b,1}/norm;
balance = n (Dot[u, {-b,1}] dF - c dF) + divJ - n rate;
pzRule = First[Solve[balance == 0, pz]];
normalPhase = Simplify[Dot[{px,pz}, normal] /. pzRule];
normalCurrent = Simplify[n (Dot[u, normal] - c/norm) /. pzRule];
Print["normal phase slope on tilted level = ", InputForm[normalPhase]];
Print["relative number current on tilted level = ", InputForm[normalCurrent]];
Print["signed order gradient along the normal = ", InputForm[Simplify[dF {-b,1}.normal]]];

(* A real, homogeneous ordinary Robin exit applied to a complex carrier
   in the ungauged sector, with the boundary at rest. *)
field = a + I b0;
Print["ordinary Robin, normal phase-current numerator = ",
      InputForm[ComplexExpand[Im[Conjugate[field] (-alpha field)]]]];
