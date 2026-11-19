(* Sharp single-time mass bounds from a signed velocity profile and net impulse. *)
ClearAll["Global`*"]; $HistoryLength=0;
emit[name_,value_]:=Print[name," = ",ToString[value,InputForm]];
(* p is total stored odd momentum and vLo,vHi are the observed odd velocity
   extrema. Extremal realizations use two nonnegative masses at those speeds. *)
weights=First[Solve[{mLo+mHi==mass,mLo vLo+mHi vHi==p},{mLo,mHi}]];
emit["masses at extrema",weights];
domain=0<vLo<vHi && p>0;
allowed=FullSimplify[(mLo>=0 && mHi>=0)/.weights,Assumptions->domain];
emit["positive masses allowed",allowed];
valid=Resolve[ForAll[{vLo,vHi,p,mass},Implies[domain,
  Equivalent[allowed,p/vHi<=mass<=p/vLo]]],Reals];
emit["same interval",valid]; If[valid=!=True,Exit[1]];
design={{vLo,vHi},{vLo/2,vHi/2}};
emit["full and half velocity rows at matched time: rank",MatrixRank[design]];
(* Independent signed readout from the raw carrier packets. *)
Do[
  pair=Table[Last[Import[FileNameJoin[{"_scratch","s11c",
          "mouth_prescribed_accel_"<>tag<>suffix<>"02","stdout"}],"RawJSON"]["32"]],
          {tag,{"positive","negative"}}];
  values=Table[(pair[[1]][key]-pair[[2]][key])/2,{key,{"carried","stress","feed_load"}}];
  emit["{suffix,t,intake,odd carried,stress,total}",Join[{suffix,pair[[1]]["t"],
       pair[[1]]["inward_rate"]},values]];
  If[Abs[Total[Most[values]]-Last[values]]>10^-20,Exit[1]],
  {suffix,{"","_half"}}];
exact=Integrate[Cos[Pi x/(2 len)]^2,{x,0,len/4},Assumptions->len>0];
emit["specified benchmark wave mass to cut",{exact,N[exact/.len->37/20,15]}];
Exit[0];
