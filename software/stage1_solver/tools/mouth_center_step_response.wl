(* Cross-check the earlier SymPy displacement response and signed push readout. *)
other = y /. First@Solve[dyx*x + dyy*y == 0, y];
travel = x /. First@Solve[dxx*x + dxy*other == (f-hold)/s, x];
perStep = Factor[D[travel,f]];
heldDerivative = D[displacement[t,support[u],u],u] /. u->0 /. support[0]->0;
heldDerivative = heldDerivative /. {support'[0]->-flowCoefficient,
    Derivative[0,1,0][displacement][t,0,0]->response[t],
    Derivative[0,0,1][displacement][t,0,0]->flowSlope};
flow = flowSlope /. First@Solve[heldDerivative==0,flowSlope];
Print["displacement per force step = ",InputForm[perStep]];
Print["fixed-time displacement slope in flow = ",InputForm[flow]];

(* Integrate the calibrated normal traction with its actual axial projection.
   Use the original archived force only to select the four drive magnitudes. *)
net = Integrate[2 Pi radius^2 Sin[theta]*Cos[theta]*
    3 f Cos[theta]/(4 Pi radius^2), {theta,0,Pi}];
archive = Import["docs/_measurements/held_force_incident_flow.json","RawJSON"];
row = SelectFirst[archive["runs"], #["name"]=="L96C6_d16"&];
If[!TrueQ[row["guard_exit"]==0 && row["worker_exit"]==0 && row["limits_verified"]],Exit[1]];
scale = ImportString[Last[row["stdout"]],"RawJSON"]["force"];
If[!TrueQ[NumericQ[scale] && scale>0],Exit[1]];
pushes = (net /. f->scale*#)& /@ {-1,-1/2,1/2,1};
Print["integrated signed push inputs ",ExportString[pushes,"RawJSON","Compact"->True]];

(* Series of a general smooth acceleration, rather than a prescribed
   polynomial in the drive. Weights are for the four signed force samples. *)
estimate = (8*(g[h/2]-g[-h/2])-(g[h]-g[-h]))/(6 h);
expansion = Expand[Normal@Series[estimate,{h,0,4}]];
Print["four-push derivative estimate through fourth order = ",InputForm[expansion]];
If[!TrueQ[Simplify[s*(dxx-dxy*dyx/dyy)*perStep]===1 &&
    Simplify[travel /. f->hold]===0 && flow===flowCoefficient*response[t] && net===f &&
    Simplify[expansion-g'[0]+h^4*Derivative[5][g][0]/480]===0], Exit[1]];
