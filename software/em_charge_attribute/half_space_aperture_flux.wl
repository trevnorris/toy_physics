(* Compute the normal momentum entering the ideal receiver directly on its
   aperture. Use the continuum radial potential on the inside of the rim. *)
ClearAll["Global`*"];
$HistoryLength = 0;
g = (1/4)/Integrate[4 Pi y^2 (1-y^2)^3, {y,0,1}];
(* Integrate the logarithmic kernel by parts. Subtract the residue at y=r
   before evaluating the remaining Cauchy principal value. *)
poly = Cancel[((1-y^2)^4-(1-r^2)^4)/(r^2-y^2)];
phi = -g/(4 Pi) (Integrate[poly,{y,0,1}] +
                       (1-r^2)^4 Log[(1+r)/(1-r)]/(2 r));
v = D[phi,r];
integrate[expr_] := NIntegrate[Evaluate[4 Pi r^2 expr],{r,0,1},
                               WorkingPrecision->40,AccuracyGoal->22,PrecisionGoal->18];
speed = g (1-r^2)^3;
carriage = Integrate[4 Pi r^2 speed^2,{r,0,1}];
pressure = -integrate[(speed^2+v^2)/2];
normal = carriage+pressure;
Print["volume rate = ", InputForm[Integrate[4 Pi r^2 speed,{r,0,1}]],
      "; aperture carriage = ",InputForm[carriage]];
Print["aperture excess pressure = ",InputForm[N[pressure,16]],
      "; aperture normal flux = ",InputForm[N[normal,16]]];
phaseContrast = (phi /. r->1/2) - (phi /. r->1/Sqrt[2]);
Print["aperture phase contrast = ",InputForm[N[phaseContrast,16]]];
If[!TrueQ[carriage>0 && pressure<0 && normal>0 &&
          Abs[phaseContrast-(-0.02830834579470561)]<10^-14],Exit[1]];
Exit[0];
