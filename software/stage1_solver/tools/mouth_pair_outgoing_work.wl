(* Isotropic three-dimensional continuum of the counterrotating branch.
   C = stiff, b = inertia, c^2=C/b.  Q is measured on ONE primitive cut.
   Coefficients may instead include both locked faces at this same flux.
   A prescribed compact change of natural curl strain P gives
      D = curl a - P, H = a_t, rho = -div P, J = P_t.
   Variation of the kinetic energy b |H|^2/2 and strain energy C |D|^2/2
   gives b H_t = -curl E, E = C D.  Consequently the outward torque power
   is E x H and the supplied power density is -E.J.  No postulated force
   on a moving center is used. For a general conserved current the same
   driven equations follow if D_t=curl H-J and b H_t=-curl E are supplied.
   All phasors in this calculation are PEAK phasors, exp(-i omega t).
   Calculate a single polarization along z. Rotational isotropy gives the
   identical scalar coefficient on each component of a vector dipole. *)
ClearAll["Global`*"];
assumptions = {r > 0, k > 0, omega > 0, c > 0, stiff > 0, radius > 0,
   s > 0, 0 < theta < Pi, Element[{p, q, amplitude, inertia}, Reals]};
reduce[expression_] := FullSimplify[expression, Assumptions -> And @@ assumptions];
emit[label_, value_] := WriteString["stdout", label, " = ", ToString[value, InputForm], "\n"];
checks = {};
check[label_, expressions_] := Module[{residuals},
  residuals = DeleteDuplicates[reduce /@ Flatten[{expressions}]];
  emit[label, residuals];
  AppendTo[checks, residuals === {0}]];

(* a = curl (p zhat G) in the dipole limit of P, where outgoing G solves
   (Laplacian+k^2)G=-delta. This construction works on the transverse
   dynamical branch and sets initial free longitudinal velocity to zero. *)
green = Exp[I k r]/(4 Pi r);
angle = -p Sin[theta] D[green, r];
dr = reduce[D[Sin[theta] angle, theta]/(r Sin[theta])];
dt = reduce[-D[r angle, r]/r];
er = stiff dr;
et = stiff dt;
hp = -I omega angle;
phasorRemoved = Exp[-I k r];
emit["outgoing {a_phi,D_r,D_theta,H_phi} without exp(i k r)",
   Factor /@ (phasorRemoved {angle, dr, dt, hp})];
check["exterior displacement divergence",
   D[r^2 dr, r]/r^2 + D[Sin[theta] dt, theta]/(r Sin[theta])];
check["exterior angular momentum balance on omega^2 b = C k^2",
   (D[r et, r] - D[er, theta])/r - I omega (stiff k^2/omega^2) hp];
(* The measured radial current contains all near-zone terms. *)
complexFlux = reduce[Pi r^2 Integrate[
    reduce[et Conjugate[hp]] Sin[theta], {theta, 0, Pi}]];
emit["sphere complex power at radius r", Collect[Expand[complexFlux], I]];
averagePower = reduce[Re[complexFlux] /. k -> omega/c];
cycleWork = reduce[2 Pi averagePower/omega];
emit["{mean outgoing power, work per period} at peak dipole p", {averagePower, cycleWork}];
check["mean current is independent of sphere radius", D[averagePower, r]];
check["electric dipole normalization and reactive current",
   complexFlux - stiff p^2 omega (k^3 + I/r^3)/(12 Pi)];

(* Independent real-time near-zone check. Values of f and its derivatives
   are evaluated at the retarded time u=t-r/c. The remainder after the
   outgoing positive acceleration power is a full derivative. *)
radialCurlTime = stiff Sin[theta] (f[u]/r^3+f'[u]/(c r^2)+f''[u]/(c^2 r))/(4 Pi);
angularSpeedTime = Sin[theta] (f'[u]/r^2+f''[u]/(c r))/(4 Pi);
powerTime = 2 Pi r^2 Integrate[radialCurlTime angularSpeedTime Sin[theta], {theta, 0, Pi}];
nearPrimitive = stiff (f[u]^2/(2 r^3)+f[u] f'[u]/(c r^2)+f'[u]^2/(c^2 r))/(6 Pi);
check["instantaneous sphere power minus positive retarded acceleration power and near derivative",
   powerTime - stiff f''[u]^2/(6 Pi c^3) - D[nearPrimitive, u]];

(* Compute cycle work on a finite source from its own retarded stress at
   the drive, independently of the flux integration. P = p zhat f, with
   f uniform on a ball, unit integral. Inside and outside, as distributions,
   D=(grad grad+k^2 I) G*P.  Only Im G contributes to -<E.J> for real p.
   In the isotropic double average the zz Hessian becomes one third of the
   Laplacian. Im G=sin(kr)/(4 Pi r) is regular even at the origin. *)
imagGreen = Sin[k r]/(4 Pi r);
laplace = D[r^2 D[imagGreen, r], r]/r^2;
check["regular absorptive Green trace equation", laplace+k^2 imagGreen];
emit["regular imaginary scalar Green at the origin", Limit[imagGreen, r -> 0]];
density = 3/(4 Pi radius^3);
overlap = 2 Integrate[Pi (radius^2-z^2), {z, r/2, radius}];
doubleMean = density^2 Integrate[4 Pi r^2 overlap imagGreen,
   {r, 0, 2 radius}, Assumptions -> k > 0 && radius > 0];
drivePower = reduce[(omega/2) stiff p^2 (2 k^2/3) doubleMean];
formFactor = reduce[4 Pi density Integrate[r^2 Sin[k r]/(k r), {r, 0, radius},
   Assumptions -> k > 0 && radius > 0]];
emit["finite ball dipole form factor at k", formFactor];
emit["finite ball drive power from integrated source work", Factor[drivePower]];
check["finite ball source work versus radiated current including size correction",
   drivePower - (averagePower /. c -> omega/k) formFactor^2];
check["point limit of ball work normalization",
   Limit[drivePower, radius -> 0] - (averagePower /. c -> omega/k)];

(* Point cores of fixed opposite signed flux, at z=+d/2 and -d/2.
   Both displacements along z, each of peak amplitude amplitude.
   Linearized convective current only. The Fourier transform at k n
   includes each propagation phase. s=k d, mu=n_z. Per unit solid angle
   the single-core power has the transverse weight 1-mu^2.
   Normalizing common or opposed by the radiation of ONE displaced core
   keeps the individual motion and the exterior coefficients the same. *)
relativeCurrent = Exp[-I s mu/2]+Exp[I s mu/2];
commonCurrent = Exp[-I s mu/2]-Exp[I s mu/2];
pattern[current_] := FullSimplify[ComplexExpand[current Conjugate[current]]];
weight = 1-mu^2;
singleNorm = Integrate[weight, {mu, -1, 1}];
pairPowers = (FullSimplify[Integrate[weight pattern[#], {mu, -1, 1},
   Assumptions -> s > 0]/singleNorm, Assumptions -> s > 0] &) /@
   {relativeCurrent,commonCurrent};
emit["parallel relative and common power / single-core power", pairPowers];
smallPair = Normal[Series[#, {s, 0, 2}]] & /@ pairPowers;
emit["parallel relative and common power / single-core power through (k d)^2", smallPair];
check["parallel relative and common interference limits", smallPair-{4-s^2/5,s^2/5}];
emit["single-core power at peak displacement amplitude", averagePower /. p -> q amplitude];
emit["all identity checks", And @@ checks];
Exit[If[And @@ checks, 0, 1]];
