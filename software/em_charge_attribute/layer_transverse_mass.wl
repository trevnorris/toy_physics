(* Expand the input stiffness T/U looking inward from the traction-free
   edge, then check the phase against a finite-frequency field solve. *)
ClearAll["Global`*"]; $HistoryLength = 0;
emit[name_, v_] := Print[name, " = ", ToString[v, InputForm]];
check[name_, a_, b_] := Module[{err=Factor[a-b]},
  emit[name,{a,b,err}]; If[err=!=0, Exit[1]]];
(* x rises from the free edge toward the entrance, q runs oppositely.
   Traction is mu dU/dq. The wave equation gives dT/dx=h omega^2 rho U. *)
logDerivative = D[t[x]/u[x],x] /. {t'[x]->h omega^2 dens[x] u[x],
                                     u'[x]->-h t[x]/stiff[x]};
riccati = Expand[logDerivative /. t[x]->stiffIn[x] u[x]];
emit["input stiffness derivative",riccati];
ansatz = omega^2 mass[x]+omega^4 quad[x];
coef = Expand[riccati /. stiffIn[x]->ansatz];
check["input stiffness quadratic source",Coefficient[coef,omega,2], h dens[x]];
check["input stiffness quartic source",Coefficient[coef,omega,4], h mass[x]^2/stiff[x]];

(* In h=rho0=c=1 solve the regular displacement series at fixed nonzero k:
   (x y')' = -k^2 weight(x) y, y(0)=1. The second solution is logarithmic.
   Match the entrance traction -y'(1) and displacement y(1) to the incident
   and reflected waves. All Taylor coefficients are exact rationals. *)
phaseDepth[weight_, wave_, degree_] := Module[{a,weights,i,j,order},
  weights=CoefficientList[weight,x];
  a=ConstantArray[0,degree+1]; a[[1]]=1;
  Do[a[[i+1]]=-(wave^2/i^2) Sum[weights[[j+1]] a[[i-j]],
                              {j,0,Min[i-1,Length[weights]-1]}],{i,1,degree}];
  order=Range[0,degree];
  N[ArcTan[-(order.a)/(wave Total[a])]/wave,36]];
Do[
  massFn = Integrate[h rho0 (weight /. x->v),{v,0,x}];
  quadEnd = Integrate[Cancel[h massFn^2/(rho0 c^2 x)],{x,0,1}];
  edgeK = ansatz /. {mass[x]->(massFn /. x->1),quad[x]->quadEnd,omega->c k};
  depth = Normal[Series[ArcTan[edgeK/(rho0 c^2 k)]/k,{k,0,2}]];
  emit["{weight,column,phase depth to k^2}",{weight,massFn/.x->1,Factor[depth]}];
  unitDepth = Expand[depth /. {h->1,c->1,rho0->1}];
  terms = Table[
    first=phaseDepth[weight,probe,24]; second=phaseDepth[weight,probe,36];
    fitted=(second-(unitDepth/.k->0))/probe^2;
    emit["{kh,phase depth,scaled correction,series-degree difference}",
         N[{probe,second,fitted,Abs[first-second]},14]];
    If[Abs[first-second]>10^-20,Exit[1]];
    Abs[fitted-Coefficient[unitDepth,k,2]],{probe,{1/10,1/20,1/40}}];
  If[!TrueQ[terms[[3]] < terms[[2]]/3 < terms[[1]]/9 && Last[terms]<10^-5],Exit[1]],
  {weight,{x,x^2(3-2x)}}];
Exit[0];
