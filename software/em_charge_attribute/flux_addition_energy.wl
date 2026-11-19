(* Express energy in units with the prefactor e^2/(8 Pi eps) removed.
   Use the pair Coulomb energy, with the spherical average of the Green kernel. *)
ClearAll["Global`*"];
emit[label_, expr_] := Print[label, " = ", InputForm[expr]];
weight[r_] := 3 r^2/(b^3 - a^3);
(* For u<r the angular average of inverse distance is 1/r. *)
angle = Integrate[1/(2 Sqrt[r^2 + u^2 - 2 r u x]), {x, -1, 1},
                  Assumptions -> 0 < u < r];
kernel = FullSimplify[angle, 0 < u < r];
oneUnit = Factor[2 Integrate[weight[r] Integrate[weight[u] kernel, {u, a, r}], {r, a, b}]];
emit["angular Green kernel for nested shells", kernel];
emit["one-unit shell coefficient", oneUnit];
emit["thin-shell limit", Limit[oneUnit, b -> a]];
If[!TrueQ[kernel == 1/r], Exit[1]];

width = 1/20;
radii = {1/2, 1, 2};
increments[q_, sgn_, k_] := Expand[k ((q + sgn)^2 - q^2)];
Do[
   kval = oneUnit /. {a -> radius, b -> radius + width};
   emit[{"radius", radius, "step sign", dir},
        {increments[0, dir, kval], increments[dir, dir, kval]}],
   {radius, radii}, {dir, {-1, 1}}];
(* Determine whether a common linear binding cost can serve at every radius.
   Require the aligned unit to beat each neighboring sector. *)
trial[n_, radius_] := (oneUnit /. {a -> radius, b -> radius + width}) n^2 - bias n;
singleWindow[radius_] := Reduce[trial[1, radius] < trial[0, radius] &&
                                trial[1, radius] < trial[2, radius], bias, Reals];
emit["individual linear-reward windows", singleWindow /@ radii];
emit["common linear-reward window", Reduce[And @@ (singleWindow /@ radii), bias, Reals]];

(* Use two reflected depth samples rather than parity-mode coordinates.
   Any real reflection-invariant Green matrix on one such pair has this form. *)
green = {{g0, g1}, {g1, g0}};
v[n_] := {bg1, bg2} + n {u1, u2};
state[1, n_] := v[n];
state[-1, n_] := Reverse[v[n]];
energy[s_, n_, lam_] := state[s, n].green.state[s, n]/2 - lam s n;
oddContrast[lam_] := Expand[(energy[1, n, lam] - energy[1, -n, lam] -
                            energy[-1, n, lam] + energy[-1, -n, lam])/4];
emit["first-unit cost difference across sides with background", Simplify[
  energy[1, 1, 0] - energy[1, 0, 0] - energy[-1, 1, 0] + energy[-1, 0, 0]]];
emit["mirror/sign double difference, written field", oddContrast[0]];
emit["mirror/sign double difference, inserted reward", oddContrast[lam]];
If[!TrueQ[oddContrast[0] == 0 && oddContrast[lam] == -lam n], Exit[1]];
