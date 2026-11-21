(* Independent laboratory potential gradient and time derivative at a translating
   spherical cut. Angular integration is exact for axisymmetric polynomials. *)
ClearAll["Global`*"];
$HistoryLength = 0;

meanOnSphere[polynomial_, d_] := Module[{terms, moment},
  terms = CoefficientList[Expand[polynomial], mu];
  moment[n_] := If[OddQ[n], 0, Pochhammer[1/2, n/2]/Pochhammer[d/2, n/2]];
  Together[Sum[terms[[i]] moment[i - 1], {i, Length[terms]}]]
];

surface[d_, q_, rho_, radius_, farPressure_, stream_, center_, dip_] :=
 Module[{area, potential, squareRadius, radial, along, transverse, partialTime,
   vx, vy, dt, normalVelocity, kinetic, pressure, sweep, integrate},
  area = 2 Pi^(d/2) radius^(d-1)/Gamma[d/2];
  squareRadius = (x-center t)^2 + y^2;
  radial = Sqrt[squareRadius];
  potential = stream x + q radius^(d-1) radial^(2-d)/((2-d) rho area)
             + dip (x-center t) / radial^d;
  along = D[potential, x];
  transverse = D[potential, y];
  partialTime = D[potential, t];
  {vx,vy,dt} = FullSimplify[
    {along,transverse,partialTime} /. {t -> 0,x -> radius mu,y -> radius Sqrt[1-mu^2]},
    -1 < mu < 1];
  normalVelocity = Expand[vx mu + vy Sqrt[1-mu^2]];
  sweep = normalVelocity - center mu;
  kinetic = Expand[rho (vx^2 + vy^2)/2];
  pressure = Expand[farPressure + rho(stream^2/2 - dt) - kinetic];
  integrate[expr_] := area meanOnSphere[expr,d] // Simplify;
  (* Store the two omitted pieces as ablation diagnostics. *)
  integrate /@ {rho sweep, rho vx sweep + pressure mu, kinetic sweep,
                pressure normalVelocity, rho dt normalVelocity, pressure center mu}
];

(* Use a lab-fixed far-pressure datum. Dipole is tied to relative incident
   speed with fixed susceptibility. Both reservoirs can have distinct radii,
   densities and dipoles. *)
readPair[rb_, rk_, sb_, sk_, u_, v_, w_] :=
 surface[3,-1/4,2,rb,3/5,u,v,sb (u-v)] +
 surface[4,1/4,3,rk,6/5,w,v,sk (w-v)];
oddBoth[pp_, pm_, mp_, mm_] := (pp - pm - mp + mm)/4;
u0=1/50;
v0=1/25;
q0=1/4;
Do[
  {rb,rk,sb,sk}=geometry;
  Do[
    u=scale u0; v=scale v0;
    matrix=Table[readPair[rb,rk,sb,sk,s u,tSign v,0],
                 {s,{1,-1}}, {tSign,{1,-1}}];
    cross=oddBoth[matrix[[1,1]],matrix[[1,2]],matrix[[2,1]],matrix[[2,2]]];
    raw=matrix[[1,1]];
    act=Table[tSign v matrix[[i,j,2]] /. tSign -> {1,-1}[[j]],{i,2},{j,2}];
    crossAct=Apply[oddBoth,Flatten[act]];
    deltaEnergy=raw[[3]]+raw[[4]]-q0 (6/5/3-3/5/2);
    needed=deltaEnergy-act[[1,1]];
    Print["cuts=", {rb,rk}, " susceptibilities=", {sb,sk}, " scale=", scale,
          " signed powers per mdot U V {advected,pressure,net,actuator,residual}=",
          InputForm[Together[{cross[[3]],cross[[4]],Total[cross[[3;;4]]],
                        crossAct,Total[cross[[3;;4]]]-crossAct}/(q0 u v)]],
          " positive quadrant powers per mdot U^2 {net_above_rest,actuator,residual}=",
          InputForm[Together[{deltaEnergy,act[[1,1]],needed}/(q0 u^2)]],
          " spurious mixed flux per mdot U V {omit_dt,subtract_pressure_sweep}=",
          InputForm[Together[{cross[[5]],-cross[[6]]}/(q0 u v)]]];
    If[!TrueQ[Together[cross[[3]] + cross[[4]]]==0 &&
             crossAct==-q0 u v && needed==3 q0 u^2/2 &&
             raw[[2]]==-q0 u], Exit[1]],
    {scale,{1,1/2}}],
  {geometry,{{1,1,0,0},{1,1,1/4,-3/7},{2,3,1/4,-3/7}}}];

(* Independently moving bulk reservoir tests the input/output speed signs and
   a common Galilean translation, while the pressure datum stays fixed. *)
Do[
  {u,v,w}=speeds;
  raw=readPair[2,3,1/4,-3/7,u,v,w];
  deltaEnergy=Total[raw[[3;;4]]]-q0 (6/5/3-3/5/2);
  required=Together[deltaEnergy-raw[[2]] v];
  Print["independent lab speeds {U,V,W}=",speeds," holder=",raw[[2]],
        " signed residual above rest=",required];
  If[!TrueQ[raw[[2]]==q0(w-u) && deltaEnergy==q0(w^2-u^2)/2 &&
            required==q0 ((w-v)^2-(u-v)^2)/2],Exit[1]],
  {speeds,{{3/100,1/100,-1/100},{3/50,1/25,1/50},{1/50,1/50,1/50}}}];
Exit[0];
