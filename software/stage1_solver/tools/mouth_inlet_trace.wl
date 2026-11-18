(* Mixed receiving stress on the full flat 3D aperture in the +w convention.
   A density fixed at the inlet leaves its normal logarithmic slope free. *)
ClearAll["Global`*"];
check[name_, lhs_, rhs_] := Module[{diff = Factor[lhs-rhs]},
  Print[name, " | ", ToString[lhs,InputForm], " | ", ToString[rhs,InputForm],
    " | residual ", ToString[diff,InputForm]];
  If[diff =!= 0, Exit[1]]];

(* mu = x/r. Compare the Cartesian derivative integrated on each shell to
   the radial derivative of the boundary moment, without fixing its profile. *)
slope = Sum[f[l][r] LegendreP[l,mu],{l,0,4}];
dxSlope = mu D[slope,r]+(1-mu^2) D[slope,mu]/r;
shellFlux = Factor[-2 Pi k n r^2 Integrate[dxSlope,{mu,-1,1}]];
boundaryFlux = Factor[-2 Pi k n r^2 Integrate[mu slope,{mu,-1,1}]];
check["shell integrand vs derivative of rim moment",shellFlux,D[boundaryFlux,r]];
check["integrated receiver stress for regular profiles",boundaryFlux,
      -(4 Pi/3) k n r^2 f[1][r]];
inlet = bath Exp[I phase];
check["sideways velocity of scalar inlet",D[inlet,x],0];
required = First[Solve[(boundaryFlux/.{r->1,k->1/4,n->1})==load, f[1][1]]];
check["unit aperture slope required to transmit incident load",
      f[1][1]/.required,-3 load/Pi];
