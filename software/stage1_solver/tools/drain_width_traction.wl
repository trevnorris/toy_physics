(* Virtual work and the portion of the normal momentum supplied by the neutral
   carrier on an intact reflecting face. All loads below are in the lab frame. *)
ClearAll["Global`*"];
$HistoryLength = 0;
pos = {mid + width/2, mid - width/2};
force = {pUp, -pLo};
readings = {force.D[pos,mid], force.D[pos,width]};
Print["center and width readings = ", InputForm[readings]];
If[Expand[readings - {pUp-pLo,(pUp+pLo)/2}] =!= {0,0}, Exit[1]];

dens = re[w]^2 + im[w]^2;
tww = D[re[w],w]^2 + D[im[w],w]^2 - D[dens,{w,2}]/4 + dens^5/5;
reflect = {re'[w] -> 0, im'[w] -> 0};
fromField = Simplify[(tww /. reflect) - dens^5/5];
Print["reflecting-face correction to pressure = ",InputForm[fromField]];
If[Simplify[fromField + (re[w] re''[w] + im[w] im''[w])/2] =!= 0, Exit[1]];

(* Thickness restoring coefficients depend on which other fields can respond.
   Compare the sealed uniform column with the fixed-chemical-potential slice. *)
hess = {{modulus,cross},{cross,kWidth}};
sealed = { -1/w0, 1}.hess.{-1/w0,1};
relaxed = 1/Inverse[hess][[2,2]];
Print["sealed minus relaxed = ",InputForm[Factor[sealed-relaxed]]];
If[Simplify[sealed-relaxed-modulus (1/w0-cross/modulus)^2] =!= 0, Exit[1]];

(* Exterior radial loading f=-a/r^4. Integrate flux balance in the case
   with zero local restoring stiffness; retain the unspecified monopole. *)
b = mono/r + particular/r^2;
row = Expand[-grad D[r^2 D[b,r],r]/r^2 + a/r^4];
coeff = First[Solve[row == 0,particular]];
trial = lead/r^4;
leading = Limit[r^4 (stiff trial-grad D[r^2 D[trial,r],r]/r^2),r->Infinity];
Print["zero stiffness particular = ",InputForm[particular /. coeff],
      " ; finite stiffness leading coefficient = ",InputForm[lead /. First[Solve[leading == -a,lead]]]];
If[Simplify[(particular /. coeff)-a/(2 grad)] =!= 0 || Simplify[leading - stiff lead] =!= 0, Exit[1]];
Exit[0];
