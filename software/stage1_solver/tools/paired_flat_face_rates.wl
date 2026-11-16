(* Direct quadrature of the paired and isolated material currents on a flat
   face annulus. Sample the velocity from the gradient of two scalar
   potentials. Positions, velocities and number currents use slab-test units.
   Both sources are inward. Each face normal points out of the slab. *)
ClearAll["Global`*"];
depth = 2;
numberDensity = 3;
constituentMass = 1;
supplyGradient = 3/25;
distanceToSource = 10;
originRate = 24 Pi;
companionRate = 24 Pi;
xyz = {x, y, z};
radial = Sqrt[xyz.xyz];
offset = xyz - {distanceToSource, 0, 0};
potential[rate_] := (originRate/radial + rate/Sqrt[offset.offset])/
                     (4 Pi depth numberDensity);
(* The flat-slab benchmark's regular radial velocity is r/100 at t=0. *)
velocity[rate_] := xyz/100 + Grad[potential[rate], xyz];
verticalCurrent = -supplyGradient w;
axis = {1, 0, 0, 0};

Do[
  normal = {0, 0, 0, sign};
  graphRule = w -> sign depth/2;
  outwardMass = constituentMass Append[{0, 0, 0}, verticalCurrent].normal /. graphRule;
  flows = Table[
    component = Append[velocity[rate], verticalCurrent/numberDensity].axis;
    integrand = 2 Pi r^2 Sin[theta] outwardMass component /.
                   {x -> r Cos[theta], y -> r Sin[theta], z -> 0};
    NIntegrate[integrand, {r, 1, 3}, {theta, 0, Pi},
               WorkingPrecision -> 40, AccuracyGoal -> 28],
    {rate, {0, companionRate}}];
  Print[If[sign == 1, "upper", "lower"],
        " sample flow/pi [isolated, paired, change] = ",
        ToString[InputForm[N[Append[flows, Last[flows]-First[flows]]/Pi, 18]]]];
  Print["outward normal stress along axis per area = ", axis.(p normal)],
  {sign, {1, -1}}];
