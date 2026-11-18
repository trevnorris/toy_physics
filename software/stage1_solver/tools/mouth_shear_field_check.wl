(* Independent frequency-domain calculation from the vorticity equation.
   On the inner side of the brane exterior, the normal for its traction is -e_r.
   Only tangential contact and outgoing curl-elastic shear enter this probe. *)
path = Environment["MOUTH_SHEAR_FIELDS"];
If[!StringQ[path], path = "_scratch/s11c/mouth_shear_field01/fields.json"];
packet = Import[path, "RawJSON"];
If[!AssociationQ[packet], Exit[1]];
rat[x_] := Rationalize[x, 10^-14];
{a, rho, c} = rat /@ Lookup[packet, {"radius", "density", "wave_speed"}];
mu = rho c^2;
k = s/c;
omega = amp Exp[-k (r-a)] (1+k r)/r^2;
residual = FullSimplify[D[r^2 D[omega,r],r] - (2+k^2 r^2) omega];
(* Recover both displacements from rho s^2 u = -mu curl(omega e_phi sin(theta)). *)
p = -2 omega/(r k^2);
d = -(D[omega,r]+omega/r)/k^2;
force = 8 Pi mu r^2 omega/3;
dragArea[g_] := 8 Pi a^2 g/3;
(* Unit acceleration switched on at t=0. *)
contact = (force - dragArea[g] (1/s^2 - s d)) /. r -> a;
rule = First[Solve[contact==0, amp]];
If[residual =!= 0 || Together[contact /. rule] =!= 0, Exit[1]];

(* Strip the travel delay before inversion and evaluate at the retarded time. *)
impulseResponse[drag_, gap_] := impulseResponse[drag,gap] = Module[{spectra},
  spectra = (Cancel[# Exp[k (r-a)]] & /@ ({p,d,force} /. rule)) /.
             {g->drag, r->a+gap};
  InverseLaplaceTransform[#,s,time] & /@ spectra
];
reference[drag_, gap_, t_, acc_] := N[acc impulseResponse[drag,gap] /. time->t-gap/c, 40];
relativeError[x_,ref_] := If[ref==0, If[x==0, 0, Infinity], Abs[(x-ref)/ref]];

runs = packet["runs"];
amps = {0,1/50,-1/50,1/100,-1/100};
expected = Tuples[{{1/2,1},amps}];
If[{a,rho,c} != {1,1,1} ||
   (rat /@ Lookup[#,{"slip_drag","acceleration"}]& /@ runs) != expected, Exit[1]];
errors = {0,0,0};
rimError = 0;
refinement = 0;
count = 0;
Do[
  {drag, acc} = rat /@ Lookup[run,{"slip_drag","acceleration"}];
  If[(rat /@ Lookup[run["samples"],"t"]) != {1/4000,1/2000,1/1000}, Exit[1]];
  Do[
    t = rat[row["t"]];
    atRim = reference[drag,0,t,acc];
    rimError = Max[rimError,relativeError[row["rim_d"],atRim[[2]]],
                   relativeError[row["rim_force"],atRim[[3]]]];
    If[(rat /@ Lookup[row["shells"],"gap"]) != {1/40000,1/10000}, Exit[1]];
    Do[
      gap = rat[shell["gap"]];
      ref = reference[drag,gap,t,acc];
      measured = Lookup[shell,{"p","d","shear_force"}];
      errors = MapThread[Max,{errors,MapThread[relativeError,{measured,ref}]}];
      refinement = Max[refinement,relativeError[shell["force_coarse"],shell["force_fine"]]];
      count++,
      {shell,row["shells"]}],
    {row,run["samples"]}],
  {run,runs}];
Print["Vorticity wave and tangential contact residuals: ", {residual,Together[contact /. rule]}];
Print["Off-rim shell samples compared: ",count];
Print["Maximum relative errors (p, d, force): ",InputForm[errors]];
Print["Maximum relative error on rim displacement and force: ",InputForm[rimError]];
Print["Maximum relative change on halving radial difference spacing: ",InputForm[refinement]];
If[!(count==60 && Max[errors]<10^-8 && rimError<10^-8 && refinement<10^-6),Exit[1]];
