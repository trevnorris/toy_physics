(* Independent route: the stream-function radial powers give the creeping
   branch. A unit normal dipole then reads the first convective force by
   the reciprocal identity, without solving the corrected pressure field. *)
ClearAll["Global`*"];
$Assumptions = a>0 && eta>0 && dens>0 && r>=a && q>0 && gam>0;
check[name_, left_, right_] := Module[{diff = FullSimplify[left-right]},
  Print[name, " | ", ToString[left,InputForm], " | ", ToString[right,InputForm],
    " | residual ", ToString[diff,InputForm]];
  If[diff =!= 0, Exit[1]]];
op[y_] := D[y,{r,2}] - 2 y/r^2;
powers = Sort[n /. Solve[Factor[op[op[r^n]] /. r->1] == 0, n]];
Print["stream function radial powers = ", ToString[powers,InputForm]];
If[powers =!= {-1,1,2,4}, Exit[1]];
stream = (w r^2/2 + c r + h/r) Sin[th]^2;
vr = FullSimplify[D[stream,th]/(r^2 Sin[th])];
vt = FullSimplify[-D[stream,r]/(r Sin[th])];
shear = FullSimplify[eta (r D[vt/r,r] + D[vr,th]/r)];
sol = First@Solve[{(vr/.{th->0,r->a})==b,(shear/.{th->Pi/2,r->a})==0},{c,h}];
(* Far stress integral of the c r mode gives the force on the opening.
   Its Stokeslet pressure is 2 eta c cos(theta)/r^2. *)
normal = 2 eta D[vr,r]-2 eta c Cos[th]/r^2;
force0 = FullSimplify[2 Pi r^2 Integrate[
    (normal Cos[th]-shear Sin[th]) Sin[th],{th,0,Pi}] /. sol];
z = 4 Pi a eta;
check["surface stress integral",force0,z (w-b)];
vr = FullSimplify[vr /. sol]; vt = FullSimplify[vt /. sol];
(* At fixed relative normal rim velocity the correction has homogeneous
   normal-velocity and tangential-traction data. The auxiliary solution
   has unit rim velocity and zero incident flow. The exterior's inner
   normal points toward the opening, giving the minus in the force. *)
testR = vr /. {w->0,b->1}; testT = vt /. {w->0,b->1};
base = -q/r^2;
advR = base D[vr,r] + vr D[base,r];
advT = base D[vt,r] + vt base/r;
shell = FullSimplify[2 Pi r^2 dens Integrate[
    (testR advR+testT advT) Sin[th],{th,0,Pi}]];
first = FullSimplify[-Integrate[shell,{r,a,Infinity},Assumptions->a>0]];
rate = 4 Pi dens q;
check["first-order reciprocal volume integral",first,rate (w/6-b/2)];
(* Determine each coefficient by balancing powers of q at the boundary. *)
orderEq = Expand[(force0+first-gam*b) /. b->b0+q b1];
bc = First@Solve[{Coefficient[orderEq,q,0]==0,Coefficient[orderEq,q,1]==0},{b0,b1}];
eff = FullSimplify[gam (b0+q b1)/w /. bc];
load[incident_,mouth_] := eff (incident-mouth);
check["flow slope and negative drift slope", D[load[u,v],u],-D[load[u,v],v]];
check["strong coupling through first order", Limit[eff, gam->Infinity], z+rate/6];
check["relative dipole changes intake by",
  2 Pi a^2 Integrate[b Cos[th] Sin[th],{th,0,Pi}],0];
