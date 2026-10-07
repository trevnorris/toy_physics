(* Blind O2 engine. Sole physical input: O2_SHARED_PHYSICS.md, sections cited below.
   No runtime reads, imports, exports, selected constitutive laws or expected results.
   Construction conventions: Cartesian ambient components (x1,x2,x3,w), force
   positive INTO material, exchange positive OUTWARD. All reduced densities d^3x.
   OPEN first variations below denote derivatives of unrestricted functionals of
   complete field sections/history, NOT functions with a closed finite jet list. *)
Begin["O2`"];
ClearAll["O2`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};
emitted = {};
emit[name_String, payload_] := Module[{tag, stream},
  tag = "WL_O2_" <> name;
  If[MemberQ[emitted, tag], Quit[91]];
  AppendTo[emitted, tag];
  stream = First[$Output];
  WriteString[stream, tag <> ": " <>
    ToString[payload, InputForm, PageWidth -> Infinity] <> "\n"];
  Flush[stream];
];
origin[section_, contract_] := <|"Spec" -> section, "Contract" -> contract|>;
operand[name_, section_, contract_] := OPEN[name, origin[section, contract]];

(* SUPPLIED INPUT / ANSATZ construction, §§1--3,5--8. No historical response
   is installed into these live inputs. Cartesian r>0; no symbol has a number. *)
x = {x1, x2, x3};
radius = Sqrt[x.x];
profiles = <|"V_r" -> VR[radius], "rho_br" -> RhoBr[radius],
  "mu_perp" -> MuPerp[radius], "xi_w" -> XiW[radius], "h" -> H[radius],
  "delta" -> DeltaOpt[radius], "j_n" -> Jn[radius], "f" -> Fbulk[radius]|>;
vPlane = profiles["V_r"] x/radius;
fieldIdentity = profiles["xi_w"] == ell profiles["h"];
opticalInputs = {CGamma[radius]^2 == profiles["mu_perp"]/profiles["rho_br"],
  CGamma[radius] == c0 (1 + profiles["delta"])};
bulkInputs = {PBulk[radius] == eosK RhoBulk[radius]^eosN,
  CSound[radius]^2 == eosN eosK RhoBulk[radius]^(eosN - 1)/particleMass,
  profiles["f"] == RhoBulk[radius]/rho0 - 1};
massDensity = profiles["rho_br"]; (* SITE K1 *)
carriedVelocity = vPlane; (* SITE K3 *)
stressReference = operand[RRefStrainLive, "3.2,4", "4"];
momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)
sourceInventory = operand[S12LocalConversionReturnControllers, "3.3,5", "1,2,6,8,10"];
boundaryInventory = operand[S12MouthCollarReturnIRBulkBoundary, "3.3,5", "1,6,7,8,10"];
branch = operand[BA13, "3.2", "2"];
inertia = operand[IBrLive, "3.2,4", "3"];
consMomentum = operand[PBrCons, "3.2,4", "3"];
fullStress = operand[TBrLive, "3.2,4", "3"];
consStress = operand[TBrCons, "3.2,4", "3"];
normalResponse = operand[NBrLive, "3.2,4", "3"];
rotation = operand[ARotLive, "3.2,4,6", "3,8"];
mapOperand = operand[JMap, "3.3,5", "6"];
core = operand[HCore, "3.3", "5,7"];
stiffnessResponse = operand[MPerp, "3.2", "6"];
densityResponse = operand[RBr, "3.2", "6"];
embeddingRelation = operand[EhLive, "3.2", "5"];
holdOperand = operand[THold[s], "3.3,4,5", "7"];
bulkAmplitude = operand[TBulkNormalLive[s], "3.3,5", "7"];
exchangeOperand = operand[PiN, "3.3,5", "1,7,10"];
partnerOperand = operand[S12AdditionalMomentumAndReactionSystem, "3.3,5", "1,2,6,8,10"];
energyOperand = operand[EBrLive, "6", "8"];
energyCurrentOperand = operand[JELive, "6", "8"];
relaxationPower = operand[PRefRelaxLive, "6", "8"]; (* SITE K11 *)
conversionPower = operand[PConvertExchangeLive, "6", "8"];
boundaryPower = operand[PBoundaryLive, "6", "8"];
supplier = operand[SENet, "6", "8"];
supplyBudget = operand[PESupply, "6", "8"];
energyReference = operand[CRef, "8.6", "8"];
bodyForceEntries = {}; (* SITE K12 *)

(* The supplied brane graph and a CONDITIONAL native height-chart ansatz.
   The native chart domain is printed on every affected output and retained
   inside its reduction/velocity actions. No native-face identification with
   xi_w, h, centre or thickness is supplied by this representation. *)
faceRepresentationDomain = <|
  "Condition" -> "Every contributing native face admits a real differentiable single-valued height over far-field x on the represented region",
  "Coverage" -> "Restricted to that height-chart domain; faces without such a chart are outside this representation",
  "NativeGeometryAndIdentifications" -> mapOperand|>;
faceQualified[object_Association] := Append[object,
  "NativeFaceRepresentationDomain" -> faceRepresentationDomain];
graphEmbedding = Append[x,profiles["xi_w"]];
slope = D[profiles["xi_w"], #] & /@ x;
tangent = D[graphEmbedding, #] & /@ x;
(* Oriented codimension-one normal density from tangent minors. For the
   ordered ambient basis (x1,x2,x3,w), the final-column cofactor is positive.
   Normalization is Euclidean on real regular charts; no normal is supplied. *)
normalDensity[rows_List] := Table[
  (-1)^(Length[First[rows]] + j) Det[
    rows[[All,Delete[Range[Length[First[rows]]],j]]]],
  {j,Length[First[rows]]}];
unitNormal[rows_List] := With[{density = normalDensity[rows]},
  density/Sqrt[density.density]];
metric = tangent.Transpose[tangent]; (* SITE K4 *)
metricInverse = Simplify[Inverse[metric], Element[x,Reals] && x.x > 0];
metricDet = Factor[Det[metric]];
graphNormal = unitNormal[tangent];
vBulk = vPlane.slope; (* SITE K2 *)
vMaterial = Join[vPlane, {vBulk}];
faceHeight = qFace[s][x1,x2,x3,t];
faceEmbedding = Append[x,faceHeight];
faceTangents = D[faceEmbedding, #] & /@ x;
faceNormalDensity = normalDensity[faceTangents];
faceArea = Sqrt[faceNormalDensity.faceNormalDensity];
faceNormal = orientation[s] unitNormal[faceTangents]; (* SITE K7 *)
bulkTraction = bulkAmplitude faceNormal;
geometry = <|"g_ij" -> metric, "g_inverse" -> metricInverse,
  "det_g" -> metricDet, "tangent_vectors" -> tangent,
  "graph_normal" -> graphNormal, "native_face_normal" -> faceNormal,
  "native_face_tangents" -> faceTangents,
  "native_face_area_factor" -> faceArea|>;
emit["BASIS_MEASURES_GEOMETRY", faceQualified[<|"Basis" -> {ex1,ex2,ex3,ew},
  "Coordinates" -> Append[x,w], "Domain" -> (radius > 0),
  "DensityMeasure" -> CoordinateVolume[x], "NativeFaceMeasure" -> faceArea CoordinateVolume[x],
  "NativeFaceOrientation" -> (orientation[s]^2 == 1),
  "NativeFaceMap" -> mapOperand, "NativeChart" -> faceHeight,
  "FieldIdentity" -> fieldIdentity, "Geometry" -> geometry,
  "MaterialVelocity" -> vMaterial, "Origin" -> origin["1,3.1,5", "5,6,7"]|>]];

(* General section calculus. The explicit entries are evaluated coordinates on
   an UNRESTRICTED section, not a finite list of constitutive arguments. The
   OtherDependence operand contains arbitrary spatial jets, nonlocal dependence,
   material histories, formation/return, normal and rotational content. Its
   derivative is deliberately unevaluated. FirstVariation is the OPEN Gateaux
   action on the computed tangent to the ENTIRE section. No momentum density is
   identified with rho V and no transport current is identified with p V. *)
materialOperands = {inertia,consMomentum,fullStress,consStress,normalResponse,
  rotation,branch,mapOperand,densityResponse,stiffnessResponse,
  sourceInventory,boundaryInventory};
section = <|"Profiles" -> profiles, "Velocity" -> vMaterial,
  "Metric" -> metric, "OtherDependence" -> UnrestrictedSection[
    materialOperands, x, t, AllSpatialJets, EntireMaterialHistory,
    UnresolvedStressInertiaNormalIdentifications]|>;
sectionTangent[sec_, coordinate_] := <|
  "Profiles" -> Map[D[#, coordinate] &, sec["Profiles"]],
  "Velocity" -> D[sec["Velocity"], coordinate],
  "Metric" -> D[sec["Metric"], coordinate],
  "OtherDependence" -> Inactive[D][sec["OtherDependence"], coordinate],
  "FurtherDependence" -> Map[Inactive[D][#,coordinate] &,
     KeyDrop[sec,{"Profiles","Velocity","Metric","OtherDependence"}]]|>;
(* These are differential actions, not separate normal forces. *)
response[role_, operands_, state_] := OpenAction[role, operands, state];
sectionDerivative[object_, sec_, coordinate_] :=
  OpenFirstVariation[object, sectionTangent[sec, coordinate]];
momentumSection = Join[section,<|"MaterialReference" -> stressReference|>];
momentumDifferentiatedSection = momentumSection /. momentumDerivativeRules;
momentumDensity = Table[response[MomentumDensity[a],
  {inertia,consMomentum,normalResponse,rotation,branch,mapOperand,
   UnresolvedStressInertiaNormalIdentifications},momentumSection], {a,4}];
momentumCurrent = Table[response[MomentumCurrent[a,i],
  {inertia,consMomentum,normalResponse,rotation,branch,mapOperand,
   UnresolvedStressInertiaNormalIdentifications},momentumSection], {a,4},{i,3}];
storage = sectionDerivative[#, momentumDifferentiatedSection, t] & /@ momentumDensity;
transport = Table[Sum[sectionDerivative[momentumCurrent[[a,i]],
  momentumDifferentiatedSection,x[[i]]],{i,3}],{a,4}]; (* SITE K8 *)
materialDerivativeProfiles = Map[Function[value,
  Sum[vPlane[[i]] D[value,x[[i]]],{i,3}]], profiles];
emit["MATERIAL_MOMENTUM", <|"DensityAction" -> momentumDensity,
  "CurrentAction" -> momentumCurrent, "Storage" -> storage,
  "SpatialTransport" -> transport,
  "DifferentiatedSection" -> momentumDifferentiatedSection,
  "MaterialProfileDerivatives" -> materialDerivativeProfiles,
  "NormalAccounting" -> UnresolvedIdentification[normalResponse,inertia,fullStress,mapOperand],
  "Calculus" -> "OpenFirstVariation acts on the whole section tangent; no finite jet/history cutoff",
  "Origin" -> origin["1,3.2,4", "2,3,4,6"]|>];

stressInputs = {fullStress,consStress,normalResponse,rotation,stressReference,
  branch,UnresolvedStressInertiaNormalIdentifications}; (* SITE K5 *)
internalForce = Table[response[InternalForce[a],stressInputs,section],{a,4}];
emit["INTERNAL_MATERIAL_FORCE", <|"Components" -> internalForce,
  "GraphNormalProjection" -> graphNormal.internalForce,
  "Origin" -> origin["3.2,4", "3,4"]|>];

(* Full load is ONE joint action constrained by its bulk part. A support split
   is not selected. The native face integral and weighting are inside JMap;
   using a chart here imposes no sheet/slab material reduction. *)
nativeHold = Table[response[NativeCompleteHold[a],
  {holdOperand,BulkPartConstraint[bulkTraction],core,
   NoDeclaredExternalSupport,UnresolvedSupportPartition},section],{a,4}];
nativeFaceSet = response[NativeBoundingFaceSet,
  {mapOperand,boundaryInventory},section];
reduceNative[object_] := response[NativeToCoordinateDensity,
  {mapOperand,nativeFaceSet,boundaryInventory,core},
  faceQualified[<|"NativeObject" -> object, "AreaFactor" -> faceArea,
    "NativeNormal" -> faceNormal, "Measure" -> CoordinateVolume[x],
    "UnrestrictedDependence" -> section|>]];
(* Sum the complete per-face reduction over the OPEN O6 face set. Function
   binds s throughout the reduced object, including geometry and application
   data. Apply constructs that binding AFTER the per-face object is evaluated.
   Inactive Map/Total retain the unknown set without choosing its cardinality
   or face identifications. The same aggregation is used for face work below. *)
totalNativeFaces[object_] := Inactive[Total][
  Inactive[Map][Function @@ {{s},object},nativeFaceSet]];
mechanicalLoad = totalNativeFaces /@ (reduceNative /@ nativeHold);
emit["MECHANICAL_LOAD", faceQualified[<|"NativeBulkTraction" -> bulkTraction,
  "NativeCompleteLoad" -> nativeHold, "NativeFaceSet" -> nativeFaceSet,
  "CoordinateComponents" -> mechanicalLoad,
  "GraphNormalProjection" -> graphNormal.mechanicalLoad,
  "Origin" -> origin["3.3,4,5", "5,6,7"]|>]];

carriedBulk = response[OutwardCarriedBulkMomentum,
  {exchangeOperand,mapOperand,normalResponse,branch,sourceInventory,boundaryInventory,
   NativeRelativeMassCurrent,Premise3LocalMaterialVelocity,
   UnspecifiedFaceToMaterialVelocity},section]; (* SITE K13 *)
carriedMomentum = Append[profiles["j_n"] carriedVelocity,carriedBulk];
sourcePartners = Table[response[OutwardAdditionalMomentumPartner[a],
  {partnerOperand,exchangeOperand,mapOperand,branch,sourceInventory,
   OPENReactionSystem,NoDuplicateMappedConvectiveCurrent},section],{a,4}];
emit["EXCHANGE_MOMENTUM", <|"Orientation" -> "Outward loss",
  "Carried" -> carriedMomentum, "AdditionalPartners" -> sourcePartners,
  "MaterialIdentification" -> SameExchangedMaterial[profiles["j_n"],exchangeOperand,mapOperand],
  "Origin" -> origin["2,3.3,5", "1,2,6,7,10"]|>];

emit["DRIVE_PROVENANCE", faceQualified[<|"SeparateBodyForceEntries" -> bodyForceEntries,
  "Drive" -> DynamicalOrderConversionDrain[branch,sourceInventory],
  "LocalSourceInventory" -> sourceInventory, "BoundaryInventory" -> boundaryInventory,
  "EntryDependencies" -> {internalForce,mechanicalLoad,carriedMomentum,sourcePartners},
  "GMInterface" -> OPEN[S16ResponseMatching], "Origin" -> origin["2,3.3,4,5", "1,2,6,7,10"]|>]];

(* Generic oriented accounting: storage and outward flux enter positively;
   applied/internal forces enter negatively. Inputs above, never a supplied
   assembled balance, determine the component expression. *)
account[entries_] := Total[(#["Orientation"] #["Object"]) & /@ entries];
entry[role_, sense_, object_, trace_] := <|"Role"->role,"Orientation"->sense,
  "Object"->object,"Origin"->trace|>;
momentumEntries = Join[{
  entry[Storage,1,storage,origin["4","3"]],
  entry[SpatialTransport,1,transport,origin["4","3"]],
  entry[InternalMaterialForce,-1,internalForce,origin["4","3,4"]],
  entry[MechanicalFaceSupport,-1,mechanicalLoad,origin["4,5","7"]],
  entry[CarriedExchange,1,carriedMomentum,origin["5","1,7,10"]],
  entry[AdditionalExchange,1,sourcePartners,origin["5","1,2,7,10"]]},
  (entry[SeparateBodyForce,-1,#,origin["2","1"]]& /@ bodyForceEntries)];
holdBalance = account[momentumEntries];
emit["B_HOLD_LIVE", faceQualified[<|"Entries" -> momentumEntries,
  "InPlane" -> Take[holdBalance,3], "BulkCoordinate" -> Last[holdBalance],
  "GraphNormalProjection" -> graphNormal.holdBalance,
  "Relation" -> Thread[holdBalance == ConstantArray[0,4]],
  "Status" -> ConditionalNamedBalance, "O4Identity" -> UnresolvedIdentification[embeddingRelation,holdBalance],
  "Origin" -> origin["4,5,9", "1,3,5,7,10"]|>]];

massFlux = massDensity vPlane;
massDivergence = Total[MapThread[D,{massFlux,x}]];
massRHS = -profiles["j_n"];
emit["MASS_INPUT", <|"DensityOperand" -> massDensity, "VelocityOperand" -> vPlane,
  "Flux" -> massFlux, "Divergence" -> massDivergence, "RHS" -> massRHS,
  "Equation" -> (massDivergence == massRHS), "Residual" -> (massDivergence - massRHS),
  "Measure" -> CoordinateVolume[x],
  "Qualification" -> "Recorded relative O(epsilon) qualification on j_n when transferred to induced measure",
  "NativeIdentification" -> mapOperand, "Origin" -> origin["3.1,5,7", "6,9"]|>];

(* Native mechanical work uses the ACTUAL application velocity as an OPEN map
   of the SAME graph velocity. It is reduced as one scalar work action, never
   a product of independently averaged force and velocity. *)
pairedVelocity = vMaterial; (* SITE K9 *)
faceVelocity = Table[response[NativeApplicationVelocity[a],
  {normalResponse,mapOperand,UnspecifiedFaceToMaterialVelocity},
  faceQualified[<|"GraphVelocity" -> pairedVelocity,"NativeChart" -> faceHeight,
    "UnrestrictedDependence" -> section|>]],{a,4}];
nativeFaceWork = nativeHold.faceVelocity;
faceWork = totalNativeFaces[reduceNative[nativeFaceWork]];
stressWork = response[MaterialStressNormalRotationalWork,
  {stressInputs,inertia,normalResponse,rotation},
  <|"ForceAction" -> internalForce,"GraphVelocity" -> vMaterial,
    "GeneralizedRates" -> operand[UnspecifiedRotationalNormalRates,"6","3,8"],
    "UnrestrictedDependence" -> section|>];
energyInputs = {energyOperand,energyCurrentOperand,relaxationPower,
  conversionPower,boundaryPower,supplier,supplyBudget,energyReference,
  stressInputs,inertia,normalResponse,rotation,branch,sourceInventory,boundaryInventory};
energySection = Join[section,<|"EnergyDependence" -> energyInputs|>];
energyDensity = response[MaterialEnergyDensity,energyInputs,energySection];
energyCurrent = Table[response[MaterialEnergyCurrent[i],energyInputs,energySection],{i,3}];
energyStorage = sectionDerivative[energyDensity,energySection,t];
energyTransport = Sum[sectionDerivative[energyCurrent[[i]],energySection,x[[i]]],{i,3}]; (* SITE K10 *)
(* Joint energy occurrence resolves no overlaps. P_ref/relax is explicitly
   located WITHIN the material stress/internal-energy accounting; conversion,
   boundary and supplier occurrences are also jointly identified, not summed
   as independent powers on top of their existing work. *)
energyPower = response[JointNetPowerOccurrence,
  {energyInputs,UnresolvedEnergyOccurrenceIdentifications},
  <|"MaterialWork" -> stressWork,"MechanicalBoundaryWork" -> faceWork,
    "ReferencePowerOccurrence" -> IncludedOccurrence[relaxationPower,
       {MaterialStressWork,MaterialInternalEnergy}],
    "ConversionOccurrence" -> IncludedOccurrence[conversionPower,{CarriedEnergy,S12Partners}],
    "BoundaryOccurrence" -> IncludedOccurrence[boundaryPower,{faceWork,OtherBoundaryTransfer}],
    "SupplyOccurrence" -> IncludedOccurrence[supplyBudget,{supplier}],
    "UnrestrictedDependence" -> energySection|>];
energyEntries = {entry[EnergyStorage,1,energyStorage,origin["6","8"]],
  entry[EnergySpatialTransport,1,energyTransport,origin["6","8"]],
  entry[JointPower,-1,energyPower,origin["6","8"]]};
energyBalance = account[energyEntries];
emit["FORCE_POWER_PAIRINGS", faceQualified[<|"GraphVelocity" -> vMaterial,
  "NativeApplicationVelocity" -> faceVelocity, "NativeFaceWork" -> nativeFaceWork,
  "CoordinateFaceWork" -> faceWork, "MaterialWork" -> stressWork,
  "SharedMap" -> mapOperand,"Measure" -> CoordinateVolume[x],
  "Origin" -> origin["1,4,5,6", "3,6,7,8"]|>]];
emit["B_E_STEADY", faceQualified[<|"Storage" -> energyStorage,"Transport" -> energyTransport,
  "PowerOccurrence" -> energyPower,"Entries" -> energyEntries,
  "Object" -> energyBalance,"Relation" -> (energyBalance == 0),
  "Supplier" -> supplier,"Budget" -> supplyBudget,"EnergyReference" -> energyReference,
  "NonPassiveObligation" -> ConditionalObligation[NonPassiveClosure,
    {NamedReservoir,StatedPowerBudget,OPENLiveSuccessorOwner}],
  "Origin" -> origin["2,6,8.6", "4,8"]|>]];

emit["COUPLED_INPUTS_MODEL_POINT", <|"CoupledInputs" ->
  {branch,stiffnessResponse,embeddingRelation,core,mapOperand,densityResponse},
  "Premises" -> AdoptedConditionalSubstrateInputs[Range[4],"2026-10-06"],
  "Profiles" -> profiles,"OpticalIdentifications" -> opticalInputs,"BulkInputs" -> bulkInputs,
  "Anchoring" -> LABHELD[CGamma,NoMaterialReferenceLaw,NoPhysicalHolder],
  "Counting" -> <|"epsilon" -> GM/(c0^2 radius),
    "OpticalMonomialBox" -> Tuples[{Range[0,1],Range[0,2],Range[0,1]}],
    "Grades" -> {OpticalGrade[DeltaOpt,1],OpticalGrade[SlopeSquared,1],
      OpticalGrade[VelocityOverC0,1/2],SeparateBulkResponseDomain[FirstOrderInF]},
    "Truncation" -> None|>,
  "MissingGradesScales" -> OPEN[{RhoBr,MuPerp,IBrLive,TBrLive,NBrLive,
    RRefStrainLive,PRefRelaxLive,PiN,THold,HCore,EhLive,AllUnsuppliedDerivativeScales}],
  "Restrictions" -> {FarField,IsolatedSphericalMassAtRest,SteadyEulerianProfiles,
    LabTime,LinearOpticalWaves,LeadingEikonal,LiveRadialProfiles,
    NoConstitutiveIsotropyParityStressSymmetryOrCoupleRestriction},
  "TransferLimits" -> {NoStrongFieldMouthInteriorMovingRotatingMassTransfer,
    NoTimeVaryingDrainTransfer,NoAngularProfileOrSwirlTransfer,
    NoDirectionDependentOpticalOrPolarizationExtension,
    NoSharpSheetMaterialReduction,NoSteadyNoDissipationAssumption,
    NoReferenceFreezeFromLABHELD},
  "HistoricalDomains" -> <|"8.1" -> HomogeneousLinearDisplacementAndPostulatedOntology,
    "8.2" -> UniformQuadraticMaterialAndStaticSingleKink,
    "8.3" -> PostulatedLocalizedParentConstantCoefficientStaticExteriorHeldMouth,
    "8.4" -> UniformTransverseAnchorAndFiniteSlabKinematicsFrozenFirstShapeProjection,
    "8.5" -> RestAcousticsFrozenPerturbationTractionAndHeldBackground,
    "8.6" -> HistoricalRealFractionOrderWorkAndN5EnergyReference|>,
  "O4EquationIdentityCount" -> Unsettled,
  "RelaxationOwner" -> Unassigned,
  "LaterInterfaces" -> {S12,S14a,S14ConditionalOnS14a,S16InferredResponseMatching,
    S21,S1p5,S8,Q1,Q2,S22},"Origin" -> origin["1,2,3,7,8,10", "1,2,3,4,5,6,7,8,9,10"]|>];
emit["LOCAL_NAMES", {}];
End[];
