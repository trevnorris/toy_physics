(* Static mouth exterior and linear waves in a curl-elastic axial field.
   flux is the measured integral of D through an outward-oriented sphere.
   The radial coefficient multiplying n/r^2 is thus flux/(4 Pi).
   No condition on the allowed numerical values of flux enters here. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[value,InputForm]];
checks = {};
compare[label_, left_, right_] := Module[{residual},
  residual = Simplify[left-right];
  emit["check " <> label <> " {left,right,residual}",{left,right,residual}];
  AppendTo[checks,AllTrue[Flatten[{residual}],#===0&]]
];
xyz = {x,y,z};
curl[field_] := Curl[field,xyz];
div[field_] := Div[field,xyz];
outer[a_,bvec_] := Outer[Times,a,bvec];
tensorDiv[tensor_] := Table[Sum[D[tensor[[i,j]],xyz[[j]]],{j,3}],{i,3}];

(* Derive the Fourier stiffness by differentiating the mode energy.
   The longitudinal initial velocity is an independent canonical datum. *)
wv = {wx,wy,wz};
vel = {vx,vy,vz};
kv = {0,0,k};
curlAmp = Cross[kv,wv]; (* the missing I cancels in the Hermitian norm *)
modeLag = b vel.vel/2 - curlAmp.curlAmp/(2 eps);
hess[expr_,coords_] := Table[D[expr,coords[[i]],coords[[j]]],{i,3},{j,3}];
mass = hess[modeLag,vel];
stiff = -hess[modeLag,wv];
waveOp = stiff-omega^2 mass;
emit["one-mode {L,canonical momenta,stiffness,inertia}",{modeLag,D[modeLag,#]&/@vel,stiff,mass}];
emit["Fourier {operator,determinant,squared transverse frequencies}",
  {waveOp,Factor[Det[waveOp]],omega2 /. Solve[waveOp[[1,1]]==0 /. omega^2->omega2,omega2]}];
longEq = mass[[3,3]] wlong''[t]+stiff[[3,3]] wlong[t]==0;
longSol = DSolveValue[{longEq,wlong[0]==u0,wlong'[0]==v0},wlong[t],t];
emit["longitudinal {equation,initial-value solution,H_L,B_L,energy density}",
  {longEq,longSol,D[longSol,t],b D[longSol,t],b D[longSol,t]^2/2}];
emit["initial div B on mode with longitudinal velocity v0",I k b v0];
emit["k=0 stiffness",stiff /. k->0];

(* D=curl W, H=Wdot, E=D/eps and B=b H. The equation of motion and
   the curl identity give Bdot=-curl E and Ddot=curl H. Use independent
   first spatial derivatives to test the local conservation laws. *)
dfield = {d1[x,y,z],d2[x,y,z],d3[x,y,z]};
hfield = {h1[x,y,z],h2[x,y,z],h3[x,y,z]};
efield = dfield/eps;
bfield = b hfield;
energy = (dfield.efield+hfield.bfield)/2;
flow = Cross[efield,hfield];
mom = Cross[dfield,bfield];
stress = outer[efield,dfield]+outer[hfield,bfield]-IdentityMatrix[3] energy;
power = efield.curl[hfield]-hfield.curl[efield];
momdot = Cross[curl[hfield],bfield]-Cross[dfield,curl[efield]];
sym = {d1[x,y,z]->dx,d2[x,y,z]->dy,d3[x,y,z]->dz,
       h1[x,y,z]->hx,h2[x,y,z]->hy,h3[x,y,z]->hz};
emit["local {energy,energy flux,field momentum,stress zz}",{energy,flow,mom,stress[[3,3]]}/.sym];
emit["time derivative div B from curl E",div[-curl[efield]]];
compare["local energy continuity",Expand[power],Expand[-div[flow]]];
compare["exterior momentum balance including divergence sources",
  Expand[tensorDiv[stress]-momdot],Expand[efield div[dfield]+hfield div[bfield]]];
speed = 1/Sqrt[b eps];
waveW = {0,amp Cos[k z-speed k t],0};
waveD = curl[waveW];
waveH = D[waveW,t];
waveDensity = Simplify[(waveD.waveD/eps+b waveH.waveH)/2];
waveFlux = Simplify[Cross[waveD/eps,waveH]];
emit["rightward real plane wave {W,D,H,energy density,energy flux}",
  {waveW,waveD,waveH,waveDensity,waveFlux}];
emit["rightward wave {kinetic energy density,curl-elastic energy density}",
  Simplify[{b waveH.waveH/2,waveD.waveD/(2 eps)}]];
compare["wave transports energy at transverse phase speed",waveFlux[[3]],speed waveDensity];

(* On a single globally smooth chart the common equator is traversed
   with opposite boundary orientations. This applies even when the ball
   inside the sphere includes an excised core. *)
edge = {r Cos[ph],r Sin[ph],0};
edgeW = {w1@@edge,w2@@edge,w3@@edge};
northLine = edgeW.D[edge,ph];
southLinePulledBack = edgeW.(-D[edge,ph]);
emit["global Stokes cap integrands on the equator, common azimuth parameter",
  {northLine,southLinePulledBack}];
globalFlux = Integrate[Simplify[northLine+southLinePulledBack],{ph,0,2 Pi}];
emit["global Stokes flux through intact sphere",globalFlux];

(* Independently integrate the radial curl equation starting at each
   regular pole. W components below use the orthonormal spherical frame. *)
radialD = flux/(4 Pi r^2);
curlSource = r Sin[th] radialD;
capAzimuthal[end_] := Simplify[
  Integrate[curlSource /. th->ang,{ang,end,th}]/Sin[th]];
wn = capAzimuthal[0];
ws = capAzimuthal[Pi];
axisymCurl[wphi_] := Simplify[{D[Sin[th] wphi,th]/(r Sin[th]),
                                 -D[r wphi,r]/r,0}];
capCurlN = axisymCurl[wn];
capCurlS = axisymCurl[ws];
emit["{target radial D, d_theta(sin theta W_phi), north W_phi, south W_phi}",
  {radialD,curlSource,wn,ws}];
emit["{north value at north pole, south value at south pole}",
  {Limit[wn,th->0,Direction->"FromAbove"],Limit[ws,th->Pi,Direction->"FromBelow"]}];
emit["north and south curls in {r,theta,phi} frame",{capCurlN,capCurlS}];
fluxN = Integrate[capCurlN[[1]] r^2 Sin[th],{ph,0,2 Pi},{th,0,Pi/2}];
fluxS = Integrate[capCurlS[[1]] r^2 Sin[th],{ph,0,2 Pi},{th,Pi/2,Pi}];
overlapW = Simplify[wn-ws];
overlapCirc = Integrate[Simplify[r Sin[th] overlapW],{ph,0,2 Pi}];
chi = Integrate[Simplify[r Sin[th] overlapW],ph];
emit["{northern and southern outward flux,sum,overlap W_phi,positive-azimuth overlap circulation}",
  {fluxN,fluxS,fluxN+fluxS,overlapW,overlapCirc}];
emit["{local transition scalar chi,chi(phi+2 Pi)-chi(phi),gradient phi component}",
  {chi,Expand[(chi /. ph->ph+2 Pi)-chi],D[chi,ph]/(r Sin[th])}];
compare["north curl",capCurlN,{radialD,0,0}];
compare["south curl",capCurlS,{radialD,0,0}];
compare["two-cap Stokes",fluxN+fluxS,overlapCirc];

(* At the equator, lambda W_N and lambda W_S lie on the same oriented
   rotation axis. Equality as SO(3) rotations requires their signed angle
   difference to be an integral multiple of 2 Pi at each point. *)
angles = Simplify[lambda {wn,ws} /. th->Pi/2];
deltaAngle = Simplify[angles[[1]]-angles[[2]]];
emit["equatorial ordinary rotation {signed angles,angle difference,radial derivative,number of turns}",
  {angles,deltaAngle,D[deltaAngle,r],deltaAngle/(2 Pi)}];
rotationDelta = {{Cos[deltaAngle],-Sin[deltaAngle],0},
                 {Sin[deltaAngle],Cos[deltaAngle],0},{0,0,1}};
emit["relative rotation trace",Tr[rotationDelta]];

(* Separate control for identifying W with an ordinary displacement.
   Pure symmetric shear and infinitesimal rigid rotation have different
   symmetric strain and curl, even in the same Cartesian frame. *)
strain[field_] := Module[{grad=Table[D[field[[i]],xyz[[j]]],{i,3},{j,3}]},
  (grad+Transpose[grad])/2];
Do[emit["ordinary displacement control {displacement,symmetric strain,curl,curl energy}",
  {disp,strain[disp],curl[disp],curl[disp].curl[disp]/(2 eps)}],
  {disp,{{shear y,shear x,0},{-spin y,spin x,0}}}];

shellEnergy = Integrate[radialD^2/(2 eps) r^2 Sin[th],{r,inner,outerR},
  {ph,0,2 Pi},{th,0,Pi},Assumptions->0<inner<outerR && eps>0];
emit["static {density,shell energy,energy outside inner radius}",
  {radialD^2/(2 eps),Factor[shellEnergy],Limit[shellEnergy,outerR->Infinity]}];

(* Use the outward normal of a sphere around the fixed flux. No
   constitutive law for the rim itself has been included. *)
normal = {Sin[th] Cos[ph],Sin[th] Sin[ph],Cos[th]};
background = {ex,ey,ez}; (* prescribed uniform electric field *)
electricStress[dvec_] := (outer[dvec,dvec]-IdentityMatrix[3] (dvec.dvec)/2)/eps;
dRad = radialD normal;
isolatedTr = Simplify[electricStress[dRad].normal];
backgroundTr = Simplify[electricStress[eps background].normal];
crossTr = FullSimplify[
  (electricStress[dRad+eps background]-electricStress[dRad]-electricStress[eps background]).normal];
sphereIntegral[v_] := Map[Integrate[# r^2 Sin[th],{ph,0,2 Pi},{th,0,Pi}]&,v];
isoForce = sphereIntegral[isolatedTr];
bgForce = sphereIntegral[backgroundTr];
crossForce = sphereIntegral[crossTr];
emit["isolated {radial stress,integrated traction}",{Simplify[normal.isolatedTr],isoForce}];
emit["uniform-background {cross traction,integrated background traction,integrated cross traction}",
  {crossTr,bgForce,crossForce}];
unitAB = {nx,ny,nz}; (* unit vector from the remote source to this sphere *)
pairForce = crossForce /. Thread[background->otherFlux unitAB/(4 Pi eps separation^2)];
emit["exterior stress force in locally uniform radial background from other flux",pairForce];
compare["uniform exterior stress cross term",crossForce,flux background];

(* Reflection of the material frame. D, E, the oriented-outward unit
   normal and the area vector transform as polar vectors; W,H are axial.
   Ordinary reflection of the entire solution preserves the signed flux. *)
reflect = DiagonalMatrix[{-1,1,1}];
dSample = {dxx,dyy,dzz};
nSample = {nxx,nyy,nzz};
wSample = {wxx,wyy,wzz};
kSample = {kxx,kyy,kzz};
emit["reflection {axial W,polar D,reflected flux integrand,original flux integrand}",
  {Det[reflect] reflect.wSample,reflect.dSample,
   (reflect.dSample).(reflect.nSample),dSample.nSample}];
compare["curl of reflected axial field",Cross[reflect.kSample,Det[reflect] reflect.wSample],
  reflect.Cross[kSample,wSample]];
compare["reflected flux integrand",(reflect.dSample).(reflect.nSample),dSample.nSample];
emit["all consistency checks",checks];
Exit[If[And@@checks,0,1]];
