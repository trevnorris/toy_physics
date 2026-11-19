(* Incremental scalar-wave energies in the archived Phase A fixed-action
   cylinder. All constants below are the exact rational versions of the
   printed Phase A inputs. The carrier action unit is tested as an explicit
   additional assumption for the scalar wave. The native I is continuous. *)
ClearAll["Global`*"];
emit[label_, values_] := WriteString["stdout", label, " = ",
    ToString[N[values, 17], InputForm], "\n"];
pressure = 3/1000; tension = 3/1000; gap = 1/5;
cut = 5/2; cap = 5/2; da = 3/100;
coords = {r, l};
wall = pressure (4 Pi/3) r^3 l + tension (4 Pi r^2 l + 8 Pi r^3/3);
wave = Sqrt[gap^2 + (Pi/r)^2 + (Pi/(2 l))^2];
energy = wall + act wave;
gradient = D[energy, {coords}]; hessian = D[gradient, {coords}];
waveSlope = D[wave, {coords}];
wp = 50;
fullPoint[i_] := fullPoint[i] = coords /. FindRoot[
    Evaluate[Thread[(gradient /. act -> i) == {0, 0}]],
    {r,N[19/10 i^(1/4),wp]}, {l,N[18/10 i^(1/4),wp]}, WorkingPrecision->wp];
reference = fullPoint[1];
fixedRPoint[i_] := fixedRPoint[i] = {First[reference],l /. FindRoot[
    Evaluate[(Last[gradient] /. {act->i, r->First[reference]}) == 0],
    {l,N[18/10 Sqrt[i],wp]},WorkingPrecision->wp]};
point[type_,i_] := Switch[type,"fixed",reference,"free length",fixedRPoint[i],"free",fullPoint[i]];
rules[type_,i_] := Join[Thread[coords -> point[type,i]],{act->i}];
level[type_,i_] := energy /. rules[type,i];
om[type_,i_] := wave /. rules[type,i];
labels = {"fixed", "free length", "free"};
residuals = {};
record[v_] := AppendTo[residuals, Max[Abs[Flatten[{v}]]]];
emit["Phase A {radius,length,omega,E} at reference", Join[reference, {om["free",1],level["free",1]}]];
record[reference - {1892410005412/10^12,181206770991/10^11}];
Do[
    eq = rules[k,1];
    compliance = Switch[k,
        "fixed", 0,
        "free length", (Last[waveSlope]^2/Last[Last[hessian]]) /. eq,
        "free", With[{v=waveSlope /. eq, h=hessian /. eq}, v.LinearSolve[h,v]]];
    eps = 1/100000;
    firstDifference = (level[k,1+eps] - level[k,1-eps])/(2 eps);
    secondDifference = (level[k,1+eps] + level[k,1-eps] - 2 level[k,1])/eps^2;
    emit[k <> " {E',E'' implicit,E'' centered,deltaI^2 E''}", {firstDifference,-compliance,secondDifference,-da^2 compliance}];
    record[{firstDifference - om[k,1],secondDifference+compliance}];
    energies = Table[level[k,1+n da],{n,0,2}];
    steps = Differences[energies];
    emit[k <> " {first,second,second-first,first/deltaI,second/deltaI}", Join[steps,Differences[steps],steps/da]];
    emit[k <> " {action,r,L,omega,curv_min,cap_margin,cut_margin,frozen_jump_relaxation_energy}",
        Table[Module[{rr=rules[k,1+n da],mincurv,loss},mincurv=Min[Eigenvalues[hessian/.rr]];
            loss = If[n==0,0,(energy/.Join[Thread[coords->point[k,1+(n-1)da]],{act->1+n da}])-energies[[n+1]]];
            {1+n da,r,l,wave,mincurv,cap-l,cut-wave,loss}/.rr], {n,0,2}]];
    Do[
        rr = rules[k,1+n da];
        record[Switch[k,"fixed",0,"free length",Last[gradient] /. rr, "free",gradient /. rr]];
        If[!TrueQ[(0<wave<cut && l<cap && r>0 && l>0)/.rr],Exit[2]],
        {n,0,2}];
    If[!TrueQ[k=="fixed" || Last[steps]<First[steps]],Exit[3]],
    {k,labels}];
threshold = {act,r,l} /. FindRoot[
    Evaluate[Join[Thread[gradient=={0,0}],{wave==cut}]],
    {act,1/4},{r,7/5},{l,7/5},WorkingPrecision->wp];
lengthThreshold = {act,r,l}/.{r->First[reference]}/.FindRoot[
    Evaluate[({Last[gradient]==0,wave==cut}/.r->First[reference])],
    {act,1/10},{l,9/10},WorkingPrecision->wp];
emit["unconstrained and fixed-radius cutoff {action,r,length}",{threshold,lengthThreshold}];
emit["unsupported low-action trials {action,r,L,omega,cut_margin}",
    Table[{i,Sequence@@fullPoint[i],om["free",i],cut-om["free",i]},{i,{da,2da}}]];
record[(Join[gradient,{wave-cut}] /.Thread[{act,r,l}->threshold])];
record[({Last[gradient],wave-cut} /. Thread[{act,r,l}->lengthThreshold])];
emit["maximum derivative, stationarity, cutoff and archived-dimension residual", Max[residuals]];
Exit[If[Max[residuals]<10^-9,0,1]];
