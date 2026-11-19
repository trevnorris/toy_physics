(* Shoot from the shear-free edge in s=-log(chi). The second state is
   chi du/dchi. No special-function solution enters the calculation. *)
ClearAll["Global`*"];
$HistoryLength = 0;
spectraOnly = MemberQ[$CommandLine,"--spectra-only"];
edge[k_?NumericQ,h_,m_,loss_,cut_:40] := Module[{pair, h2},
  h2 = With[{arg=SetPrecision[k,36]}, N[h^2 (arg^2+I loss arg),32]];
  pair = NDSolveValue[{u'[s]==-p[s], p'[s]==h2 Exp[-(m+1)s] u[s],
                       u[cut]==1, p[cut]==0}, {u[0],p[0]}, {s,0,cut},
                      WorkingPrecision->30, AccuracyGoal->21, PrecisionGoal->20];
  {First[pair], -Last[pair]/(h k)}];
rimNode[k_?NumericQ,y_,h_,m_] := Module[{a,b},
  {a,b}=edge[k,h,m,0]; a Cos[k (y-h/2)]-b Sin[k (y-h/2)]];
mode[y_,h_,m_,{lo_,hi_}] := k /. FindRoot[rimNode[k,y,h,m]==0,{k,lo,hi},
  Method->"Brent",WorkingPrecision->28,AccuracyGoal->16,PrecisionGoal->16];
branches[y_,h_,m_] := mode[y,h,m,#]& /@ {{1,2},{7/2,21/4}};
reflection[k_,y_,h_,m_,loss_,cut_:40] := Module[{a,b},
  {a,b}=edge[k,h,m,loss,cut];
  Exp[2 I k (y-h/2)] (I a-b)/(I a+b)];
spacingLength[{k0_,k1_}] := Sqrt[2] Pi/Sqrt[k1^2-k0^2];
emit[row_] := WriteString[First[$Output],ExportString[row,"RawJSON","Compact"->True],"\n"];

Do[
  probes=branches[1,1/10,m];
  Do[
    {name,y,h}=geometry;
    roots=If[name=="base",probes,branches[y,h,m]];
    spectrum=<|"modes"->N[roots],"spacing_length"->N[spacingLength[roots]],
               "spacing_shift"->N[spacingLength[roots]-spacingLength[probes]]|>;
    If[spectraOnly,
      emit[Join[<|"m"->m,"geometry"->name,"loss"->0|>,spectrum]];
      Continue[]];
    Do[
      refl=reflection[#,y,h,m,loss]& /@ probes;
      reference=reflection[#,1,1/10,m,loss]& /@ probes;
      phase=Arg[refl/reference];
      row=<|"m"->m,"geometry"->name,"loss"->N[loss],
            "power"->N[Abs[refl]^2],"relative_phase"->N[phase],
            "phase_shift"->N[phase/(2 probes)]|>;
      If[loss==0,row=Join[row,spectrum]];
      emit[row];
      If[Max[Abs[refl]]>1+10^-12 || Max[Abs[phase]]>=Pi/2,Exit[1]];
      If[loss==0 && Max[Abs[Abs[refl]-1]]>10^-12,Exit[1]];
      If[Max[Abs[refl-(reflection[#,y,h,m,loss,48]& /@ probes)]]>10^-13,Exit[1]],
      {loss,{0,1/10}}],
    {geometry,{{"base",1,1/10},{"shift",21/20,1/10},{"widen",1,1/5}}}],
  {m,{0,1}}];
Exit[0];
