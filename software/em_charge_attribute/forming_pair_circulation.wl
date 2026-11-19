(* Independent reduction: exact oscillator modes in x, Fourier derivatives in
   w, classical RK4 in time. Hermite quadrature integrates the x moments. *)
ClearAll["Global`*"];
h = 1./16; dt = 1./128; n = 288; steps = 512;
w = N[Range[-144,143]/16];
freq = N[2 Pi/18 Join[Range[0,143],Range[-144,-1]]];
d2[v_] := Re[InverseFourier[-freq^2 Fourier[v, FourierParameters->{0,-1}],
                            FourierParameters->{0,-1}]];
d1[v_] := Re[InverseFourier[I freq Fourier[v, FourierParameters->{0,-1}],
                            FourierParameters->{0,-1}]];
depth[t_, split_] := If[split, 2 (t/4)^2 (3-2 t/4), 0.];
well[t_, split_] := (w^2-depth[t,split]^2)^2/4;
rhs[state_, t_, split_] := With[{pot=well[t,split]},
  {state[[3]],state[[4]], d2[state[[1]]]-(3+pot) state[[1]],
                           d2[state[[2]]]-(1+pot) state[[2]]}];
nodes = N[Flatten[Table[sign Sqrt[(3+branch Sqrt[6])/2],
                        {branch,{-1,1}},{sign,{-1,1}}]]];
weights = N[Flatten[Table[Sqrt[Pi](1-branch Sqrt[2/3])/4,
                          {branch,{-1,1}},{2}]]];

observe[state_, t_, split_, seed_, side_] := Module[
  {f,g,vf,vg,fw,gw,weight,pot,tot,xx,wc,px,pw,jj,band,bandEnergy,
   i,x,u,vel,ux,uz,en,mx,mw},
  {f,g,vf,vg}=state; fw=d1[f]; gw=d1[g]; pot=well[t,split];
  weight=N[(Sign[side w]+1)/2] h;
  band=N[Map[Which[#<1,1,#==1,1/2,True,0]&,Abs[Abs[w]-2]]];
  tot=xx=wc=px=pw=jj=0.; bandEnergy=0.;
  Do[
    x=nodes[[i]];
    u=x f+seed g; vel=x vf+seed vg;
    ux=(1-x^2)f-x seed g; uz=x fw+seed gw;
    en=(vel^2+ux^2+uz^2+(x^2+pot)u^2)/2;
    mx=-vel ux; mw=-vel uz;
    tot+=weights[[i]] weight.en;
    xx+=weights[[i]] x weight.en; wc+=weights[[i]] weight.(w en);
    px+=weights[[i]] weight.mx; pw+=weights[[i]] weight.mw;
    jj+=weights[[i]] weight.(x mw-w mx);
    bandEnergy+=weights[[i]] weight.(band en),
    {i,Length[nodes]}];
  <|"L"->jj-(xx pw-wc px)/tot,"E"->tot,"x"->xx/tot,"w"->wc/tot,
    "well_band"->bandEnergy/tot|>];

Do[
  state = {Exp[-w^2/2],0. w,0. w,(w^2-.5) Exp[-w^2/2]};
  Do[
    time=step dt;
    If[MemberQ[{0,256,512},step],
      Do[Print[ExportString[<|"forming"->split,"t"->time,"seed"->seed,
        "up"->observe[state,time,split,seed,1],
        "down"->observe[state,time,split,seed,-1]|>,"RawJSON","Compact"->True]],
        {seed,{.01,-.01,.005,0.}}]];
    If[step<steps,
      k1=rhs[state,time,split]; k2=rhs[state+dt k1/2,time+dt/2,split];
      k3=rhs[state+dt k2/2,time+dt/2,split]; k4=rhs[state+dt k3,time+dt,split];
      state+=dt (k1+2 k2+2 k3+k4)/6],
    {step,0,steps}],
  {split,{False,True}}];
