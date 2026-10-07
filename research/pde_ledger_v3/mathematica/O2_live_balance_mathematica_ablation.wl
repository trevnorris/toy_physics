(* Run this harness through s11c_guarded_run.py, 4 GiB, default TasksMax.
   Exactly one kernel: Get executes the live file and each single-mutated copy
   serially inside this kernel. The engine itself reads no file. All artifacts
   written by this harness are below the repository's ignored _scratch.
   Usage: math -script /absolute/path/O2_live_balance_mathematica_ablation.wl
   No result values, signs or expected physical outcomes are tested. *)
Begin["O2Harness`"];
ClearAll["O2Harness`*"];
$HistoryLength = 0;
engine = FileNameJoin[{DirectoryName[$InputFileName],"O2_live_balance_mathematica_audit.wl"}];
repo = Nest[DirectoryName,DirectoryName[$InputFileName],3];
scratch = CreateDirectory[FileNameJoin[{repo,"_scratch","o2_wl_ablation_" <>
  StringReplace[CreateUUID[],"-"->""]}]];
source = Import[engine,"Text"];
knives = {
 {"K1","massDensity = profiles[\"rho_br\"]; (* SITE K1 *)",
   "massDensity = RhoBrConstant; (* SITE K1 *)"},
 {"K2","vBulk = vPlane.slope; (* SITE K2 *)",
   "vBulk = 0; (* SITE K2 *)"},
 {"K3","carriedVelocity = vPlane; (* SITE K3 *)",
   "carriedVelocity = Array[IndependentCarriedVelocity,3]; (* SITE K3 *)"},
 {"K4","metric = tangent.Transpose[tangent]; (* SITE K4 *)",
   "metric = IdentityMatrix[3]; (* SITE K4 *)"},
 {"K5","stressInputs = {fullStress,consStress,normalResponse,rotation,stressReference,\n  branch,UnresolvedStressInertiaNormalIdentifications}; (* SITE K5 *)",
   "stressInputs = {fullStress,consStress,normalResponse,rotation,\n  branch,UnresolvedStressInertiaNormalIdentifications}; (* SITE K5 *)"},
 {"K6a","momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)",
   "momentumDerivativeRules = {VR -> Function[{z},VRConstant]}; (* SITE K6a K6b K6c *)"},
 {"K6b","momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)",
   "momentumDerivativeRules = {RhoBr -> Function[{z},RhoBrConstant]}; (* SITE K6a K6b K6c *)"},
 {"K6c","momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)",
   "momentumDerivativeRules = {XiW -> Function[{z},XiWConstant]}; (* SITE K6a K6b K6c *)"},
 {"K7","faceNormal = orientation[s] Join[-faceSlope, {1}]/faceArea; (* SITE K7 *)",
   "faceNormal = orientation[s] {0,0,0,1}; (* SITE K7 *)"},
 {"K8","transport = Table[Sum[sectionDerivative[momentumCurrent[[a,i]],\n  momentumDifferentiatedSection,x[[i]]],{i,3}],{a,4}]; (* SITE K8 *)",
   "transport = ConstantArray[0,4]; (* SITE K8 *)"},
 {"K9","pairedVelocity = vMaterial; (* SITE K9 *)",
   "pairedVelocity = Array[IndependentPowerVelocity,4]; (* SITE K9 *)"},
 {"K10","energyTransport = Sum[sectionDerivative[energyCurrent[[i]],energySection,x[[i]]],{i,3}]; (* SITE K10 *)",
   "energyTransport = 0; (* SITE K10 *)"},
 {"K11","relaxationPower = operand[PRefRelaxLive, \"6\", \"8\"]; (* SITE K11 *)",
   "relaxationPower = 0; (* SITE K11 *)"},
 {"K12","bodyForceEntries = {}; (* SITE K12 *)",
   "bodyForceEntries = {Array[IndependentBodyForce,4]}; (* SITE K12 *)"},
 {"K13","carriedBulk = response[OutwardCarriedBulkMomentum,\n  {exchangeOperand,mapOperand,normalResponse,branch,sourceInventory,boundaryInventory,\n   NativeRelativeMassCurrent,Premise3LocalMaterialVelocity,\n   UnspecifiedFaceToMaterialVelocity},section]; (* SITE K13 *)",
   "carriedBulk = profiles[\"j_n\"] vBulk; (* SITE K13 *)"}
};
writeTag[name_, value_] := (WriteString[First[$Output],name <> ": " <>
  ToString[value,InputForm,PageWidth->Infinity] <> "\n"]; Flush[First[$Output]];);
(* Strings and semantic records have a formal structural difference; algebraic
   leaves have corrupted-minus-baseline. Identical records, including strings,
   produce the CAS object 0, never suppression of their tag. *)
difference[a_,b_] /; SameQ[a,b] := 0;
difference[a_Association,b_Association] /; Keys[a] === Keys[b] :=
  AssociationThread[Keys[a],MapThread[difference,{Values[a],Values[b]}]];
difference[a_List,b_List] /; Length[a] === Length[b] := MapThread[difference,{a,b}];
difference[a_String,b_String] := StructuralDifference[b,a];
difference[a_List,b_List] := StructuralDifference[b,a];
difference[a_Association,b_Association] := StructuralDifference[b,a];
difference[a_Equal,b_Equal] := MapThread[difference,{List@@a,List@@b}];
difference[a_,b_] := Expand[b-a];
parseCapture[path_] := Module[{lines,parts},
  lines = Select[StringSplit[Import[path,"Text"],"\n"],StringLength[#]>0&];
  parts = StringSplit[#,": ",2]& /@ lines;
  If[!AllTrue[parts,Length[#]===2 && StringStartsQ[First[#],"WL_O2_"]&],Quit[92]];
  If[Length[DeleteDuplicates[First /@ parts]] != Length[parts],Quit[93]];
  Association[(First[#] -> Block[{$Context="O2Payload`",$ContextPath={"System`"}},
     ToExpression[Last[#],InputForm]])& /@ parts]
];
runFile[file_,capture_] := Module[{stream},
  stream = OpenWrite[capture,PageWidth->Infinity];
  Block[{$Output={stream}},Get[file]];
  Close[stream];
  parseCapture[capture]
];
writeTag["WL_LOCAL_O2_ABLATION_MANIFEST",<|"Engine"->engine,"Scratch"->scratch,
 "Knives"->knives,"DifferenceConvention"->"corrupted minus baseline; structural difference for semantic records"|>];
baseline = runFile[engine,FileNameJoin[{scratch,"baseline.stdout"}]];
writeTag["WL_LOCAL_O2_ABLATION_BASELINE_TAGS",Keys[baseline]];
Do[
  Module[{knife=First[spec],old=spec[[2]],new=spec[[3]],count,file,corrupted},
    count = StringCount[source,old];
    writeTag["WL_LOCAL_O2_ABLATION_"<>knife<>"_SITE",<|"Occurrences"->count,"Source"->old,"Replacement"->new|>];
    If[count != 1,
      writeTag["WL_LOCAL_O2_ABLATION_"<>knife<>"_CONSTRUCTION",Missing["ConstructionSite",old]];
      Quit[94]];
    file = FileNameJoin[{scratch,knife<>".wl"}];
    Export[file,StringReplace[source,old->new],"Text"];
    corrupted = runFile[file,FileNameJoin[{scratch,knife<>".stdout"}]];
    writeTag["WL_LOCAL_O2_ABLATION_"<>knife<>"_TAGS",Keys[corrupted]];
    If[Keys[corrupted] =!= Keys[baseline],Quit[95]];
    Do[writeTag["WL_LOCAL_O2_ABLATION_"<>knife<>"_"<>tag,
      <|"Baseline"->baseline[tag],"Corrupted"->corrupted[tag],
        "Difference"->difference[baseline[tag],corrupted[tag]]|>],{tag,Keys[baseline]}];
  ],{spec,knives}];
writeTag["WL_LOCAL_O2_ABLATION_LOCAL_NAMES",Join[
 {"WL_LOCAL_O2_ABLATION_MANIFEST","WL_LOCAL_O2_ABLATION_BASELINE_TAGS"},
 Flatten[Table[Join[{"WL_LOCAL_O2_ABLATION_"<>First[spec]<>"_SITE",
   "WL_LOCAL_O2_ABLATION_"<>First[spec]<>"_TAGS"},
   ("WL_LOCAL_O2_ABLATION_"<>First[spec]<>"_"<>#& /@ Keys[baseline])],{spec,knives}]],
 {"WL_LOCAL_O2_ABLATION_LOCAL_NAMES"}]];
End[];
