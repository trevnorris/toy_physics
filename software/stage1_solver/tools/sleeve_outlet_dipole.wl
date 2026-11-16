(* Cartesian check on a sphere inside the flat outlet. Local radial and
   vertical Taylor coefficients are independent at the evaluation radius.
   After taking the first variation the surviving angular polynomials have
   degree at most two, so two Gauss-Legendre points integrate them exactly. *)
ClearAll["Global`*"];
radius = Sqrt[x^2 + y^2 + z^2];
cosine = x/radius;
offset = radius - r0;
logDensity = h0 + hr offset + hw w + hrw offset w +
             eps cosine (b0 + br offset + bw w + brw offset w);
velocityPotential = f0 + fr offset + fw w +
                    eps cosine (g0 + gr offset + gw w);
density = Exp[logDensity];
massFlux = m density D[velocityPotential, w];
carried = massFlux D[velocityPotential, x];
quantum = -hb^2/(4 m) density D[logDensity, x, w];

angularMean[term_] := FullSimplify[
  Total[Table[term /. {x -> t r0/Sqrt[3], y -> r0 Sqrt[2/3], z -> 0, w -> 0},
              {t, {-1, 1}}]]/2, r0 > 0];
Do[
  expr = item[[2]];
  Print[item[[1]], " [isolated mean, dipole derivative mean] = ",
        ToString[InputForm[{angularMean[expr /. eps -> 0],
                            angularMean[D[expr, eps] /. eps -> 0]}]]],
  {item, {{"mass", massFlux}, {"carried", carried}, {"density-gradient stress", quantum}}}];
