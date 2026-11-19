(* Independent graph calculation from div(grad h / W) in radial and
   cos(polar angle) coordinates. Curvature is the average of three roots. *)
ClearAll["Global`*"];
gauss[n_]:=gauss[n]=Module[{a=ConstantArray[0.,{n,n}],v,q},
 Do[a[[i,i+1]]=a[[i+1,i]]=N[i/Sqrt[4i^2-1]],{i,n-1}];
 {v,q}=Eigensystem[a]; SortBy[Transpose[{v,2q[[All,1]]^2}],First]];
gate[s_]:=10s^3-15s^4+6s^5;
kernel[g_] := Module[{h,hr,hm,w,hmean},
 h=(r^2/2-b r m g)/5;
 hr=D[h,r]; hm=D[h,m];
 w=Sqrt[1+hr^2+(1-m^2)hm^2/r^2];
 hmean=(D[r^2 hr/w,r]+D[(1-m^2)hm/w,m])/(3r^2);
 Compile[{{r,_Real},{m,_Real},{b,_Real}},Evaluate[2Pi r^2 w {1,hmean^2}]]];
funcs={kernel[gate[r-1]],kernel[1]};
part[b_,region_,n_] := Module[{q=gauss[n],fn=funcs[[region]]},
 Total[Flatten[Table[fn[region+(1+x[[1]])/2,y[[1]],N[b]] x[[2]]y[[2]]/2,
   {x,q},{y,q}],1]]];
wall[b_,n_] := Table[part[b,i,n],{i,2}];
(* Use the independent archived fine Wolfram frequency difference at the
   crossing. The wider energy differences use the archived fine spline mesh. *)
manifest=Import["docs/_measurements/wave_at_level_patch.json","RawJSON"];
lines=manifest["runs"]["sympy"]["stdout_claim_lines"];
freq=Association[(N[ToExpression[First[#]]]->ToExpression[Last[#]])& /@
 Flatten[StringCases[lines,RegularExpression[
 "^cells=36 angular=14 tilt=([0-9.]+) omega=([0-9.]+).*"]:>{"$1","$2"}],1]];
If[Length[freq]!=7,Exit[1]];
wlLine=Last[manifest["runs"]["wolfram"]["stdout_claim_lines"]];
waveForce=ToExpression[Last[StringCases[wlLine,RegularExpression["0\\.[0-9]+"]]]];
base=wall[0,64]; {a0,b0}=Total[base];
Print["baseline join, conversion = ",ToString[base,InputForm]];
Do[Module[{parts=wall[b,64],fine=wall[b,48],da,db,dw},
 {da,db}=Total[parts-base]; dw=freq[N[b]]-freq[0.];
 Print["tilt, area gain, bending gain, frequency drop = ",
       ToString[{b,da,db,-dw},InputForm]];
 Print["scale balanced offset and slope = ",
       ToString[{da+3a0 dw/freq[0.],db+b0 dw/freq[0.]},InputForm]];
 If[Max[Abs[Flatten[parts-fine]]]>5.*^-8,Exit[1]]],{b,{2,2.5,3,4}}];
Do[Module[{der=(Total[wall[2.5+d,64]]-Total[wall[2.5-d,64]])/(2d)},
 Print["step, area slope, bend slope, force-balance action coefficients = ",
 ToString[{d,der[[1]],der[[2]],der[[1]]/waveForce,der[[2]]/waveForce},InputForm]]],
 {d,{0.01,0.001}}];
