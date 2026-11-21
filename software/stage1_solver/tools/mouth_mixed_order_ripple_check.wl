(* Homogeneous transverse velocities proportional to Exp[I (k x - omega t)].
   y = omega rho (1-q)/b, and the reference speed is Sqrt[mu/rho]. *)
ClearAll["Global`*"];
$HistoryLength = 0;
ordered = rho q;
disordered = rho (1 - q);
stress = mu q;
dragLaw = b -> rho q (1 - q)/tau;
emit[label_, value_] := Print[label, " = ", InputForm[value]];

(* Retain the two velocities separately, then eliminate the unforced one. *)
operator = {{b - I omega ordered + I stress k2/omega, -b},
            {-b, b - I omega disordered}};
dispersion = k2 /. First[Solve[Det[operator] == 0, k2]];
ratio = vd /. First[Solve[operator[[2]] . {1, vd} == 0, vd]];
scaled[expr_] := Factor[expr /. b -> omega disordered/y];
z = scaled[(mu/rho) dispersion/omega^2];
{zr, zi} = Factor /@ ComplexExpand[{Re[z], Im[z]}];
limits = {Limit[z, y -> 0], Limit[z, y -> Infinity]};
emit["disordered/ordered velocity", scaled[ratio]];
emit["squared slowness times mu/rho", z];
emit["its real and imaginary parts", {zr, zi}];
emit["low and high limits", limits];
emit["loss circle operands", Factor /@ {zi^2, (zr - Last[limits]) (First[limits] - zr)}];
emit["loss / excess real slowness squared", Factor[zi/(zr - Last[limits])]];

(* In the uniform release, relative momentum exchange is the only force. *)
accelerations = {-b (vo - vd)/ordered, b (vo - vd)/disordered};
emit["center momentum rate", Factor[{ordered, disordered} . accelerations]];
emit["relative acceleration / slip", Factor[(Subtract @@ accelerations)/(vo - vd) /. dragLaw]];
emit["y from the release time", Factor[omega disordered/b /. dragLaw]];

(* For unit local ordered velocity, -<stress velocity> = mu q Re[k]/(2 omega).
   Its decrease per length uses 2 Re[k] Im[k] = Im[k^2]. *)
heat = scaled[(b/2) ComplexExpand[Abs[1 - ratio]^2]];
fluxDecrease = Factor[stress omega (rho/mu) zi/2];
emit["drag heat and flux decrease per squared velocity", {heat, fluxDecrease}];

(* Evaluate the right-going principal square root; retain the attenuation sign. *)
Do[
  waveNumber = N[Sqrt[z /. {q -> 1/4, y -> probe}], 15];
  emit["{y, phase speed/reference speed, attenuation*reference speed/omega, attenuation/Re[k]}",
       {probe, 1/Re[waveNumber], Im[waveNumber], Im[waveNumber]/Re[waveNumber]}],
  {probe, {1/10, 1, 10}}];
