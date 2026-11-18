(* Fixed-width exterior and a conditional static thickness diagnostic. *)
ClearAll["Global`*"];
$HistoryLength = 0;
lap[f_] := D[r^2 D[f, r], r]/r^2;
(* At weak intake, build a complex radial wave with a phase gradient of
   order eps and a real amplitude correction of order eps^2. The sign of
   the phase gives an inward flow in the outward radial coordinate. *)
wave = (1 + eps^2 amp[r]) Exp[I eps q/r];
potential = ((1 + eps^2 amp[r])^8 - 1)/4;
row = Expand[SeriesCoefficient[
    Exp[-I eps q/r] (-lap[wave]/2 + potential wave), {eps, 0, 2}]];
Print["order-two amplitude equation = ", InputForm[row]];
If[Simplify[row + lap[amp[r]]/2 - 2 amp[r] - q^2/(2 r^4)] =!= 0, Exit[1]];
(* Solve successive algebraic coefficients of the density deficit. *)
tail = Sum[c[i]/r^i, {i, 1, 6}]/2;
equation = Expand[row /. {amp[r] -> tail, amp'[r] -> D[tail, r],
                         amp''[r] -> D[tail, {r, 2}]} /. r -> 1/z];
solution = First[Solve[Table[Coefficient[equation, z, i] == 0, {i, 1, 6}],
                       Table[c[i], {i, 1, 6}]]];
Print["number-density coefficients 1 through 6 = ", InputForm[Table[c[i], {i, 6}] /. solution]];
If[(Table[c[i], {i, 6}] /. solution) =!= {0,0,0,-q^2/2,0,-3 q^2/2}, Exit[1]];

(* Static variational slice from the recorded uniform quadratic energy.
   Width and density are independent here; a live draining state need
   not relax to this slice. The applied load has compact support. *)
energy = modulus theta^2/2 + coupling theta width + kWidth width^2/2 + stiff slope^2/2;
relaxed = First[Solve[D[energy, theta] == 0, theta]];
restoring = Simplify[D[D[energy, width] /. relaxed, width]];
Print["relaxed width stiffness = ", InputForm[restoring]];
radialTail = Exp[-mu r]/r;
widthResidual = FullSimplify[(restoring radialTail - stiff lap[radialTail])/radialTail];
Print["exterior exponential-over-radius residual divided by field = ", InputForm[widthResidual]];
If[Simplify[restoring - kWidth + coupling^2/modulus] =!= 0 ||
   Simplify[widthResidual - restoring + stiff mu^2] =!= 0, Exit[1]];
Exit[0];
