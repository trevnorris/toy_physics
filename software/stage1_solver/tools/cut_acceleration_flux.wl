(* Check the same decaying ambiguity using the complex carrier linearization
   and the Cartesian momentum tensor, retaining the complete spherical cut. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, expr_, target_] := Module[{v=FullSimplify[expr]},
  Print[label, " = ", InputForm[v]];
  If[!TrueQ[FullSimplify[v==target]],Exit[1]]
];
hb = 1/Sqrt[6];
eta[w_] := (-3 Exp[-2 w]+8 Exp[-3 w])/2 + I (-2 Exp[-2 w]+3 Exp[-3 w])/hb;
etaReal[w_] := (-3 Exp[-2 w]+8 Exp[-3 w])/2;
residual = -hb^2/2 (eta''[w]+2 I eta'[w]/hb-eta[w]) + 2 (35/24) etaReal[w];
emit["linearized complex receiver equation", residual, 0];
(* At the cut psi=e^(iw/hb)(1+eps Y eta[w]), base density and speed are 1. *)
psi = Exp[I w/hb] (1+eps y eta[w]);
derX = D[psi,y] yx;
n = ComplexExpand[Conjugate[psi] psi] // FullSimplify;
currentN = hb ComplexExpand[Im[Conjugate[psi] D[psi,w]]];
currentX = hb ComplexExpand[Im[Conjugate[psi] derX]];
first[expr_] := FullSimplify[D[expr,eps] /. {eps->0, w->0}];
mixed = hb^2 (ComplexExpand[Re[Conjugate[derX] D[psi,w]]] - yx D[n,w,y]/4);
emit["same normal current",first[currentN],0];
emit["advective current change",first[currentN currentX/n],yx];
emit["total Cartesian mixed current change",first[mixed],7 yx/4];
(* Use radial/angular integration of d_x Y, with Y(r,theta)=j(r) Cos[theta].
   Regular j(0)=0; j(a)=1. *)
angular = Integrate[2 Pi r^2 (mu^2 j'[r]+(1-mu^2)j[r]/r),{mu,-1,1}];
emit["angular integral is a radial derivative",angular-D[4 Pi r^2 j[r]/3,r],0];
edge = (4 Pi r^2 j[r]/3 /. r->a) /. j[a]->1;
emit["total on the complete section", first[mixed] /. yx->edge,7 Pi a^2/3];
Exit[0];
