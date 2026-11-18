(* Prescribe the signed equatorial traction q in the Laplace domain.
   Derive the outward solution from the radial wave operator. *)
trial = Exp[-s (r-a)/c] (b/r + j);
radial = Factor[Exp[s (r-a)/c] (D[trial,{r,2}] - (2/r^2+s^2/c^2) trial)];
constants = First[Solve[{radial==0,(trial /. r->a)==a q/(rho c^2)}, {b,j}]];
h = trial /. constants;
y = Factor[(-c^2 D[h,r]/(a s q)) /. r->a];
surfaceResidual = Factor[(rho c^2 h/a /. r->a)-q];
Print["Exterior wave and imposed traction residuals: ",{Factor[radial /. constants],surfaceResidual}];
Print["Velocity per net equatorial traction: ",InputForm[y]];
If[Factor[radial /. constants] =!= 0 || surfaceResidual =!= 0, Exit[1]];

path = Environment["BRANE_DIRECT_TRACTION_FIELDS"];
If[!StringQ[path], path = "_scratch/s11c/direct_brane_traction01/velocities.json"];
packet = Import[path,"RawJSON"];
If[!AssociationQ[packet], Exit[1]];
rat[x_] := Rationalize[x,10^-14];
{radius,density,speed,duration} = rat /@ Lookup[packet,{"radius","density","speed","duration"}];
(* q=q_peak (t/T)^3 over the observation window. Causality permits its
   cubic continuation beyond T when inverting the transform on [0,T]. *)
unitRamp = InverseLaplaceTransform[
    (y /. {a->radius,rho->density,c->speed}) 6/(duration^3 s^4),s,time];
levels = {"384","768","1536"};
errors = {0,0,0};
If[(rat /@ Lookup[packet["runs"],"peak"]) != {0,1/50000,-1/50000,1/100000,-1/100000},Exit[1]];
Do[
  peak = rat[run["peak"]];
  Do[
    rows = run["grids"][levels[[index]]];
    If[(rat /@ Lookup[rows,"t"]) != duration {1/4,1/2,1},Exit[1]];
    Do[
      tRead = rat[row["t"]];
      reference = N[peak unitRamp /. time->tRead,40];
      If[Abs[row["traction"]-peak(tRead/duration)^3]>10^-17,Exit[1]];
      err = If[reference==0,If[row["velocity"]==0,0,Infinity],
               Abs[(row["velocity"]-reference)/reference]];
      errors[[index]] = Max[errors[[index]],err],
      {row,rows}],
    {index,Length[levels]}],
  {run,packet["runs"]}];
Print["Full positive ramp continuum rim velocities at T/4, T/2, T: ",
      InputForm[N[(unitRamp/50000 /. time->#)& /@ (duration {1/4,1/2,1}),16]]];
Print["Maximum relative velocity errors at 384, 768, 1536 radial cells: ",InputForm[errors]];
If[!TrueQ[errors[[3]]<10^-4 && 3.5<errors[[1]]/errors[[2]]<4.5 &&
          3.5<errors[[2]]/errors[[3]]<4.5],Exit[1]];
