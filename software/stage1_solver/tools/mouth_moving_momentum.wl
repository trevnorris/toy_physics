(* Derive the interfacial balance from a finite interval with an internal
   moving boundary. The spatial derivative of each primitive is momentum
   density. Bulk conservation sets its time derivative to minus the flux
   up to a spatial constant, which cancels on each interval. *)
inventory = left[s[t],t]-left[lo,t]+stored[t]+right[hi,t]-right[s[t],t];
rate = D[inventory,t] /. {
  Derivative[0,1][left][pos_,t] :> -pl[pos],
  Derivative[0,1][right][pos_,t] :> -pr[pos]};
sheetRate = stored'[t] /. First[Solve[
  rate == f + pl[lo]-pr[hi]-edge,stored'[t]]];
names = {pl[s[t]]->pin, pr[s[t]]->pout,
  Derivative[1,0][left][s[t],t]->gin, Derivative[1,0][right][s[t],t]->gout,
  s'[t]->wn};
sheetRate = Expand[sheetRate /. names];
Print["stored rate from finite moving interval = ",InputForm[sheetRate]];
If[Expand[sheetRate - (f + pin-gin*wn - pout+gout*wn - edge)] =!= 0,Exit[1]];

(* Compare the total carrier budget in the lab and in a frame following
   the mouth. Include particle transport across the edge as well. *)
numberSupply = jin-jout-jedge;
boostedBudget = D[labP[t]-m*v[t]*num[t],t] -
    (f+(pin-m*v[t]*jin)-(pout-m*v[t]*jout)-(edge-m*v[t]*jedge));
change = Simplify[boostedBudget /. {labP'[t]->f+pin-pout-edge,
   num'[t]->numberSupply}];
Print["accelerating-frame remainder = ",InputForm[change]];
If[Simplify[change+m*num[t]*v'[t]] =!= 0,Exit[1]];

(* Use the mixed carrier stress in amplitude variables. Perturb only the
   receiving normal amplitude derivative. Both phase derivatives, the
   surface density and the normal speed of the cut are shared. *)
mixed[normalAmp_] := 2 beta*(u'[x]*normalAmp-u[x]*D[normalAmp,x]);
inc = u[x]^2*(vn-wn);
carriage = m*inc*vx;
sigmaIn = mixed[d[x]];
sigmaOut = mixed[d[x]+eps*h[x]*u[x]/2];
jump = Simplify[(carriage+sigmaIn)-(carriage+sigmaOut)];
Print["two-trace flux difference = ",InputForm[jump]];
If[Simplify[jump-beta*eps*u[x]^2*h'[x]] =!= 0,Exit[1]];

(* The archived binary stores the real face, imaginary face and inward
   current in C order. Use the adjoint derivative of masked density for
   the integral of the trial stress, retaining the full 3D aperture. *)
file = "_scratch/s11c/mouth_rim_push_sympy01/resting_face.f64";
If[FileHash[file,"SHA256","HexString"] =!=
   "e76812f953ed1af4fd216a5d7fb9d70da6801536dfc786a48fcea165c8ca8053",Exit[1]];
{re,im,current} = ArrayReshape[BinaryReadList[file,"Real64",ByteOrdering->-1],{3,32,32,32}];
dx = 1/8;
coords = Range[-31,31,2]/16;
radius2 = Table[i^2+j^2+k^2,{i,coords},{j,coords},{k,coords}];
mask = Map[Boole[#<1]&,radius2,{3}];
trial = N[Table[i*Max[1-4*(i^2+j^2+k^2),0]^3,{i,coords},{j,coords},{k,coords}]];
der[field_] := (RotateRight[field,2]-8 RotateRight[field]+8 RotateLeft[field]-RotateLeft[field,2])/(12 dx);
integral[field_] := dx^3 Total[Flatten[field]];
density = re^2+im^2;
rateNumber = integral[current*mask];
slopeCoefficient = -integral[trial*der[mask*density]]/4;
oneSide[velocity_] := integral[mask*current*((re*der[im]-im*der[re])/density-velocity)];
frameShift = oneSide[0.02]-oneSide[0];
Print["saved rate, two-face coefficient, single-face frame shift = ",
   ExportString[{rateNumber,slopeCoefficient,frameShift},"RawJSON","Compact"->True]];
Print["integrated trial differences at amplitudes 0, .01, -.01, .005 = ",
   ExportString[N[{0,.01,-.01,.005}*slopeCoefficient],"RawJSON","Compact"->True]];
If[!TrueQ[Max[Abs[Flatten[(1-mask)*der[trial]]]]==0 &&
     Abs[rateNumber-1/4]<10^-13 && Abs[slopeCoefficient]>10^-13 &&
     Abs[frameShift+.02*rateNumber]<10^-14],Exit[1]];
