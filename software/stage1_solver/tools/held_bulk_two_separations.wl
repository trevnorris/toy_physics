(* Independent 4D Neumann solution for two identical unit spheres. Both local
   axes point toward the other center, so symmetry gives one set of coefficients.
   On-axis values U_l(1)=l+1 and the Taylor series of (d-x)^(-l-2) give
   the translation into regular harmonics r^n U_n. *)
ClearAll["Global`*"];
$HistoryLength = 0;
Print["Wolfram version = ", $Version];
solve[d_, degree_] := Module[
  {translation, scatter, change, matrix, coeff, slope},
  translation = Table[(ell+1)/(n+1) Binomial[ell+n+1,n]/d^(ell+n+2),
    {n,1,degree},{ell,0,degree}];
  (* d_r[C_n r^(-n-2) + A_n r^n] at unit radius must vanish for n>0. *)
  scatter = N[DiagonalMatrix[Range[degree]/(Range[degree]+2)].translation,50];
  matrix = IdentityMatrix[degree]-scatter[[All,2;;]];
  coeff = Prepend[LinearSolve[matrix, -scatter[[All,1]]/2],-1/2];
  (* Differentiate the same boundary problem at the midpoint without differencing
     the inferred B values at neighboring distances. *)
  change = Table[-(ell+n+2)/d scatter[[n,ell+1]],{n,1,degree},{ell,0,degree}];
  slope = LinearSolve[matrix,change.coeff];
  2 {coeff[[2]],slope[[1]]}
];
fmt[z_] := ToString[N[z,14],InputForm];
Do[
  reference = solve[4,degree];
  Print["L=",degree," midpoint {B,Bprime}=",fmt[reference]];
  slopes = Table[
    left = First@solve[4-h,degree]; right = First@solve[4+h,degree];
    secant = (right-left)/(2h);
    Print["h=",h," Bminus=",fmt[left]," Bplus=",fmt[right],
      " secant=",fmt[secant]," error vs midpoint derivative=",fmt[secant-Last[reference]]];
    secant,{h,{1/10,1/20,1/40}}];
  extrapolated = (4 Last[slopes]-slopes[[-2]])/3;
  Print["extrapolated derivative=",fmt[extrapolated],
    " difference from direct derivative=",fmt[extrapolated-Last[reference]]];
  If[Abs[extrapolated-Last[reference]] > 10^-9,Exit[1]];
  If[degree==24,
    Print["midpoint change with harmonic order=",fmt[Max[Abs[reference-coarse]]]];
    If[Max[Abs[reference-coarse]]>10^-10,Exit[1]]];
  coarse=reference,
  {degree,{14,24}}];
