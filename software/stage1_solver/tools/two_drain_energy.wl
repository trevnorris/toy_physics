(* Gauss radial accumulation for the volume cross energy, followed by full
   finite-sphere quadrature of the moving-cut fluxes. The companion field is
   spatially varying in these quadratures. No numerical throat is constructed. *)
ClearAll["Global`*"];
$HistoryLength = 0;
Print["Wolfram version = ",$Version];
emit[label_,expr_] := Print[label," = ",ToString[Factor[expr],InputForm]];
check[label_,lhs_,rhs_] := Module[{residual=FullSimplify[lhs-rhs]},
 Print[label," | ",ToString[lhs,InputForm]," | ",ToString[rhs,InputForm],
       " | residual ",ToString[residual,InputForm]];
 If[residual=!=0,Exit[1]]];

Do[
 omega=2 Pi^(n/2)/Gamma[n/2];
 (* Center spheres about source 1. Gauss gives zero companion radial flux
    before enclosing source 2, and q2/rho afterwards. Removing either small
    disjoint ball leaves the cross energy unchanged by harmonic mean value. *)
 energy=Integrate[q1 q2/(rho omega radius^(n-1)),{radius,d,Infinity},
                  Assumptions->rho>0 && d>0];
 edot=D[energy,d](vel2-vel1);
 c21=-q2/(rho omega d^(n-1)); c12=q1/(rho omega d^(n-1));
 (* Tensor identity before integration: self radial speed g; exact companion
    components cr, cx. Unsteady pressure follows from each moving potential. *)
 crossTensorNormal=rho (g mu cr+g cx-g cr mu);
 force1=-FullSimplify[crossTensorNormal/.{g->q1/(rho area),cx->c21}] area;
 force2=-FullSimplify[crossTensorNormal/.{g->q2/(rho area),cx->c12}] area;
 mechanical=FullSimplify[force1 vel1+force2 vel2];
 Print["dimension = ",n];
 emit["common-domain interaction energy",energy];
 emit["interaction energy rate",edot];
 emit["mechanical work rate",mechanical];
 (* Pressure/current energy from each moving source, with Reynolds term. *)
 j1=q1 vel2 c21; j2=q2 vel1 c12;
 check["exterior balance",edot,j1+j2];
 emit["inner energy required at mouth 1",j1+force1 vel1];
 emit["inner energy required at mouth 2",j2+force2 vel2];
 emit["whole interaction energy requirement",j1+j2+mechanical],
 {n,{3,4}}];

(* Numeric inputs are controls only. Both signs of mass flow, unequal rates,
   different cut radii and two velocity allocations are tested. All quadrature
   operands are built from the full two potentials, not the force formula. *)
cutRead[n_,density_,qi_,qj_,radius_,separation_,vi_,vj_,direction_] := Module[
 {omega,measure,mu,dist,g,cx,cr,phiSelf,phiComp,derComp,
  momentum,energyCross,pressureTime,labEnergy,movingEnergy,green,integrands},
 omega=2 Pi^(n/2)/Gamma[n/2];
 mu=Cos[theta];
 measure=2 Pi^((n-1)/2)/Gamma[(n-1)/2] radius^(n-1) Sin[theta]^(n-2);
 dist=Sqrt[radius^2+separation^2-2 direction radius separation mu];
 g=qi/(density omega radius^(n-1));
 cx=qj (radius mu-direction separation)/(density omega dist^n);
 cr=qj (radius-direction separation mu)/(density omega dist^n);
 phiSelf=qi radius^(2-n)/((2-n) density omega);
 phiComp=qj dist^(2-n)/((2-n) density omega);
 derComp=cr;
 momentum=density (g mu cr+cx g-g cr mu);
 energyCross=density g cr;
 (* Cross part of (rho e + rho v^2/2 + p)v is
    -rho[dot(phi_i)*v_j+dot(phi_j)*v_i]. *)
 labEnergy=density (vi g mu cr+vj cx g);
 movingEnergy=labEnergy-energyCross vi mu;
 (* Symmetric Green boundary density with outward normal of exterior=-n_i. *)
 green=-density (phiSelf derComp+phiComp g)/2;
 integrands={green,-momentum,movingEnergy,
             -energyCross vi mu,labEnergy};
 Table[NIntegrate[Evaluate[measure integrand],{theta,0,Pi},
                  WorkingPrecision->40,AccuracyGoal->27,PrecisionGoal->25],
       {integrand,integrands}]
];
Do[
 Do[
  With[{density=If[n==3,2,3],qone=sign/4,qtwo=sign*3/8,
        separation=8,oneCut=cuts[[1]],twoCut=cuts[[2]],
        va=speeds[[1]],vb=speeds[[2]]},
   one=cutRead[n,density,qone,qtwo,oneCut,separation,va,vb,1];
   two=cutRead[n,density,qtwo,qone,twoCut,separation,vb,va,-1];
   omega=2 Pi^(n/2)/Gamma[n/2];
   referenceEnergy=qone qtwo/(density omega (n-2) separation^(n-2));
   referenceRate=-qone qtwo (vb-va)/(density omega separation^(n-1));
   observed={one[[1]]+two[[1]],one[[2]] va+two[[2]] vb,
             one[[3]]+two[[3]],one[[3]]+one[[2]] va,two[[3]]+two[[2]] vb};
   Print["finite-cut n,sign,cuts,velocities = ",InputForm[{n,sign,cuts,speeds}]];
   Print["energy,mechanical,field influx,inner1,inner2 = ",InputForm[N[observed,18]]];
   residual=Max[Abs[observed-{referenceEnergy,referenceRate,referenceRate,referenceRate,referenceRate}]];
   Print["largest absolute analytic/quadrature residual = ",InputForm[residual]];
   Print["omit unsteady pressure or Reynolds transport, energy defects = ",
         InputForm[N[{referenceRate-one[[4]]-two[[4]],referenceRate-one[[5]]-two[[5]]},18]]];
   If[residual>10^-23,Exit[1]]
  ],{cuts,{{1,2},{2,3}}},{speeds,{{1/100,-1/200},{0,-1/100}}}],
 {n,{3,4}},{sign,{-1,1}}];
