(* Stationary scalar momentum through a smooth aperture and the distant receiver. *)
ClearAll["Global`*"];
$HistoryLength = 0;
components = {u[r,w],v[r,w]};
dens = components.components;
potential = (dens^4-1)/4;
lap[f_] := D[f,{r,2}]+2 D[f,r]/r+D[f,{w,2}];
rw = D[components,r].D[components,w]-D[dens,r,w]/4;
ww = D[components,w].D[components,w]+(dens^5-1)/5-D[dens,{w,2}]/4;
div = D[r^2 rw,r]/r^2+D[ww,w];
equations = lap /@ components - 2 potential components;
(* Keep derivatives of the equations as well as the equations themselves:
   the density Hessian in the stress brings in third derivatives. *)
identity = Simplify[Expand[div -(D[components,w].equations-
                                      components.D[equations,w])/2]];
Print["off-shell normal momentum identity residual = ",InputForm[identity]];
(* At an open interior point the two sides have the same value, normal
   derivative, lateral Laplacian and EOS. Use the field equation for f_ww. *)
traceStress[am_,slope_,lat_] := With[{nn=am.am},
   (nn^5-1)/5+(slope.slope-am.((nn^4-1) am/2-lat))/2];
common = {a,b}; commonSlope={c,d}; lateral={e,f};
jump = Simplify[traceStress[common,commonSlope,lateral]-
                 (ww /. {Derivative[0,2][u][r,w]->((a^2+b^2)^4-1)a/2-e,
                         Derivative[0,2][v][r,w]->((a^2+b^2)^4-1)b/2-f,
                         Derivative[0,1][u][r,w]->c,Derivative[0,1][v][r,w]->d,
                         u[r,w]->a,v[r,w]->b})];
Print["smooth-aperture normal traction jump = ",InputForm[jump]];
(* Far-field control for a radially spreading ideal 4D receiver of unit density. *)
areaDensity=4 Pi rad^3 Sin[theta]^2;
speed=q/Integrate[areaDensity,{theta,0,Pi/2}];
export=Simplify[Integrate[areaDensity(speed^2-speed^2/2)Cos[theta],{theta,0,Pi/2}]];
Print["ideal hemispherical normal export = ",InputForm[export]];
If[!TrueQ[identity==0 && jump==0 && Limit[export,rad->Infinity]==0],Exit[1]];
Exit[0];
