(* Recover the reaction from a translated Cartesian potential and from the
   carrier dipole's potential variation. Vr absorbs the confinement width. *)
radius = Sqrt[(x-shift)^2 + y^2 + z^2];
potential = v[radius-rim];
wallVariation = -D[potential,shift] /. shift->0;
carrierVariation = -D[potential,x] /. shift->0;
atPoint = {x->r Cos[theta],y->r Sin[theta],z->0};
dipole = FullSimplify[wallVariation /. atPoint, Assumptions->r>0 && Element[theta,Reals]];
drive = Integrate[2 Pi Sin[theta]*dn Cos[theta]*dipole, {theta,0,Pi}];
Print["wall force per radial volume measure and unit dipolar drive = ", InputForm[drive]];
If[!TrueQ[Simplify[wallVariation+carrierVariation]===0 &&
    Simplify[drive-4 Pi dn v'[r-rim]/3]===0], Exit[1]];

(* Complete the square in the harmonic energy, without an X velocity term.
   Force/displacement in Fourier space is the relaxed quadratic coefficient. *)
harmonic = (k*u^2+(om^2-w^2)*b^2)/2-g*u*b;
centered = Expand[harmonic /. b->mode+g*u/(om^2-w^2)];
kernel = Factor[D[centered,{u,2}]];
slowMass = Simplify[-(D[kernel,{w,2}] /. w->0)/2];
fast = Limit[kernel,w->Infinity];
Print["Fourier force/displacement = ", InputForm[kernel]];
Print["slow induced coefficient and fast limit = ", InputForm[{slowMass,fast}]];
If[!TrueQ[Simplify[D[centered,mode,u]]===0 && slowMass===g^2/om^4 && fast===k],Exit[1]];

(* In the closed rectangular path, the two amplitude legs cancel since the
   restoring derivative is phase independent. The two phase legs give the
   difference between their currents. Differentiate the local current first. *)
normalization = j+2 e*m1+e^2*m2;
intake = j*base*(1+e*h)^2/normalization;
dj = Expand[D[intake,e] /. e->0];
weighted = Expand[-h*dj/base];
curl = Factor[Coefficient[weighted,h,1]*m1 + Coefficient[weighted,h,2]*m2];
Print["oriented work per infinitesimal trace rectangle = ",InputForm[curl]];
If[!TrueQ[Simplify[curl+2*(m2-m1^2/j)]===0],Exit[1]];

(* Raw, checksummed checkpoint face shared with the previous paired field
   readout. Compute the work around a finite trace rectangle from the laws. *)
file = "_scratch/s11c/mouth_rim_push_sympy01/resting_face.f64";
If[FileHash[file,"SHA256","HexString"] =!=
    "e76812f953ed1af4fd216a5d7fb9d70da6801536dfc786a48fcea165c8ca8053",Exit[1]];
{re,im,flux} = ArrayReshape[BinaryReadList[file,"Real64",ByteOrdering->-1],{3,32,32,32}];
dv = (1/8)^3;
axis = N[Range[-31,31,2]/16];
weight = Table[Max[1-i*i-l*l-p*p,0.]^3,{i,axis},{l,axis},{p,axis}];
direction = MapThread[#1*#2&,{axis,weight}];
integral[field_] := dv Total[Flatten[field]];
rate = integral[flux];
amplitude = Sqrt[re^2+im^2];
speedShape = weight/amplitude^2;
If[!TrueQ[Max[Abs[Flatten[flux-(1/4)*weight/integral[weight]]]]<10^-13 &&
    Abs[rate-1/4]<10^-13],Exit[1]];
pairedCurrent[am_,test_] := Module[{raw=speedShape*am^2},
  rate*integral[raw*test]/integral[raw]];
loop[test_,delta_] := (pairedCurrent[amplitude*(1-delta*test),test]-
                       pairedCurrent[amplitude*(1+delta*test),test])/(2 delta);
variance = -2 integral[flux*(direction-integral[flux*direction]/rate)^2];
loops = loop[direction,#]& /@ {0.0001,0.00005};
uniform = loop[ConstantArray[1.,{32,32,32}],0.0001];
Print["rate, variance, rectangle values at full and half width, uniform control ",
  ExportString[{rate,variance,Sequence@@loops,uniform},"RawJSON","Compact"->True]];
If[!TrueQ[variance<0 && Max[Abs[loops/variance-1]]<2*10^-8 && Abs[uniform]<10^-10],Exit[1]];
