(* One guarded kernel. The canonical engine is loaded unchanged for baseline.
   Each copy has exactly one manifest replacement. No subprocess kernels.
   Run from this repository root; all harness artifacts stay in _scratch. *)
Begin["S9bHarness`"];
$HistoryLength = 0;
(* Optional scope changes execution selection only, never a knife or payload. *)
scope = Replace[Environment["S9B_ABLATION_SCOPE"], $Failed -> "ALL"];
If[!MemberQ[{"ALL", "K11"}, scope], Quit[94]];
root = Directory[];
engine = FileNameJoin[{root, "research", "pde_ledger_v3", "mathematica", "S9b_light_bending_mathematica_audit.wl"}];
work = FileNameJoin[{root, "_scratch", "s9b_wl_repair7", "ablation-final4"}];
If[!DirectoryQ[work], CreateDirectory[work, CreateIntermediateDirectories -> True]];
source = Import[engine, "Text"];
S9bHarness`progress = OpenWrite[FileNameJoin[{work, "progress.log"}]];
note[stage_, name_] := (WriteString[S9bHarness`progress,
  ToString[{stage, name, MemoryInUse[], MaxMemoryUsed[]}, InputForm] <> "\n"]; Flush[S9bHarness`progress]);
knives = <|
 "K1" -> {"kineticSymbol = (omega - vVector.kCov)^2;", "kineticSymbol = omega^2;"},
 "K2" -> {"metric2 = DiagonalMatrix[{aMetric, r^2}];", "metric2 = DiagonalMatrix[{1, r^2}];"},
 "K3" -> {"dispersionSpeed = cSquared;", "dispersionSpeed = c0^2;"},
 "K4a" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] /. HoldPattern[Derivative[j_][deltaProfile][arg_]] :> 0;"},
 "K4b" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] /. HoldPattern[Derivative[j_][velocityProfile][arg_]] :> 0;"},
 "K4c" -> {"xiSlope = Simplify[D[xiAnsatz, r], $Assumptions];", "xiSlope = Simplify[D[xiConstant, r], $Assumptions];"},
 "K5" -> {"returnDirection = -1;", "returnDirection = 1;"},
 "K6" -> {"radarSource = roundTrip;", "radarSource = oneWayER;"},
 "K7" -> {"quantifierMode = \"EVERY_B\";", "quantifierMode = \"FIXED_B\";"},
 "K8" -> {"massDensity = rhoBr[r];", "massDensity = rhoBrConstant;"},
 "K9a" -> {"localSoundSquared = bulkSoundSquared /. rho -> rho0 (1 + fLocal);", "localSoundSquared = bulkSoundSquared /. rho -> rho0;"},
 "K9b" -> {"localBulkDensity = rho0 (1 + fLocal);", "localBulkDensity = rho0;"},
 "K10" -> {"branchGateLocal = propagatingLocal;", "branchGateLocal = Not[propagatingLocal];"},
 "K11" -> {"vVector = {vLocal, 0};", "vVector = {vLocal, vLocal/r};"}|>;
print[name_, value_] := (WriteString[First[$Output], "WL_LOCAL_S9B_ABLATION_" <> name <> ": " <>
  ToString[value, InputForm, PageWidth -> Infinity] <> "\n"]; Flush[First[$Output]]);
(* BEGIN_DOMAIN_DIFFERENCE_METHOD
   Domain attachment is a property of the engine's payload schema, not of
   a predicate's head. In particular Closedness, ProfileEquation, Constraint,
   ConstraintOnOtherInputs and NoRealVelocityDomain remain computed VALUES.
   ConditionalExpression conditions and the named attachment fields below
   are collected at their complete structural paths, including nested cases.
   A predicate emitted alone is a value, not an attached domain. Coverage is
   never passed to this projection (K11's sample domain is not an attachment).
   Baseline and corrupted payloads are always printed in full separately. *)
attachedDomainKeys = {"Domain", "Domains", "BranchDomain", "AlongRayDomain",
  "ConstantDomain", "GeneralConstantDomain", "FunctionalDomain", "Premise",
  "Component", "ZeroSetJoining"};
splitPayload[payload_] := Module[{domains = <||>, paths = {}, walk, value},
  walk[z_Association, path_] := Association[KeyValueMap[
    Function[{key, item}, If[MemberQ[attachedDomainKeys,key],
      AssociateTo[domains, ToString[Append[path,key],InputForm] -> item];
      AppendTo[paths,{path,item}]; Nothing,
      key -> walk[item,Append[path,key]]]], z]];
  walk[z_List, path_] := MapIndexed[walk[#1,Join[path,#2]] &,z];
  walk[z_ConditionalExpression, path_] := (
    AssociateTo[domains, ToString[Append[path,"ConditionalExpression"],InputForm] -> z[[2]]];
    AppendTo[paths,{path,z[[2]]}]; walk[z[[1]],path]);
  (* Traverse other heads as well: nested ConditionalExpression attachments
     are separated without changing a relation into a domain. *)
  walk[z_, path_] := If[AtomQ[z] || FreeQ[z,_ConditionalExpression],z,
    MapIndexed[walk[#1,Join[path,#2]] &,z]];
  value = walk[payload,{}]; <|"Value" -> value,"Domain" -> domains,"DomainPaths" -> paths|>];
logicalObjectQ[z_] := MemberQ[{True,False},z] ||
  MemberQ[{Equal,Unequal,Less,LessEqual,Greater,GreaterEqual,Inequality,
    And,Or,Not,Xor,Equivalent,Implies,Element,Exists,ForAll},Head[z]] ||
  MemberQ[{Inactive[ForAll],Inactive[Exists]},Head[z]];
(* Mathematical comparison is three-valued. Normalize inside containers and
   quantifiers, then ask the CAS for equality/equivalence. Opaque functionals
   are replaced by shared algebraic indeterminates only during a proof:
   a proved identity in that abstraction is valid for the original objects.
   A nonzero structural/algebraic residual alone NEVER proves a change. *)
normalCache = <||>;
opaqueQ[z_] := Module[{head = Head[z]}, !AtomQ[z] &&
  (Head[head] =!= Symbol || Context[Evaluate[head]] =!= "System`")];
(* Exact polynomial identity proof. Every non-polynomial operand, including
   complete reciprocals and noninteger powers, is a shared indeterminate
   during the proof and is restored afterward. No relation between those
   operands is assumed. No unknown denominator is canceled. An unresolved
   equality stays UNDECIDED with complete original operands and residual.
   Every term participates. No size cutoff, sampling or timed fallback. *)
polynomialNormal[z_] := Module[{objects,slots,algebra},
  objects = DeleteDuplicates[Cases[z,q_ /; !AtomQ[q] &&
    (!MemberQ[{Plus,Times,Power},Head[q]] ||
     (Head[q] === Power && (!IntegerQ[Last[q]] || Last[q] < 0))),{0,Infinity}]];
  slots = Table[Unique["S9bHarness`proofAtom"],{Length[objects]}];
  algebra = z /. Thread[objects -> slots];
  Expand[algebra] /. Thread[slots -> objects]];
relationNormal[z_] := Module[{head=Head[z],d,factor},
  If[Length[z] != 2,Return[z]];
  d=polynomialNormal[z[[1]]-z[[2]]];
  factor=If[Head[d]===Times && NumberQ[First[d]],First[d],1];
  (* Positive numeric scaling preserves inequalities; any nonzero numeric
     scaling preserves equality. Never divide by an unknown profile. *)
  If[MemberQ[{Equal,Unequal},head] || TrueQ[factor>0],d=d/factor];
  head[d,0]];
booleanNormal[z_] := Module[{atoms,slots,algebra,result},
  atoms=DeleteDuplicates[Cases[z,q_ /; !MemberQ[{True,False},q] &&
    !MemberQ[{And,Or,Not,Xor,Equivalent,Implies},Head[q]],{0,Infinity}]];
  slots=Table[Unique["S9bHarness`proofBoolean"],{Length[atoms]}];
  algebra=z /. Thread[atoms->slots];
  result=BooleanMinimize[algebra];
  result /. Thread[slots->atoms]];
proofSimplify[z_] := Which[
  MemberQ[{Equal,Unequal,Less,LessEqual,Greater,GreaterEqual},Head[z]],relationNormal[z],
  MemberQ[{And,Or,Not,Xor,Equivalent,Implies},Head[z]],booleanNormal[z],
  logicalObjectQ[z],z,
  True,polynomialNormal[z]];
canonical[z_] := Module[{key = HoldComplete[z],result,head = Head[z]},
  If[KeyExistsQ[normalCache,key],Return[normalCache[key]]];
  result = Which[
    !FreeQ[z, HoldPattern[Inactive[Integrate][0, ___]]],
      canonical[z /. HoldPattern[Inactive[Integrate][0, ___]] -> 0],
    AtomQ[z] || StringQ[z] || head === Function || head === HoldForm || head === HoldComplete, z,
    AssociationQ[z], AssociationMap[canonical[z[#]] &,Keys[z]],
    ListQ[z], canonical /@ z,
    MemberQ[{ForAll,Exists,Inactive[ForAll],Inactive[Exists]},head],
      head @@ Prepend[canonical /@ Rest[List @@ z],First[z]],
    head === ImplicitRegion, With[{condition = canonical[z[[1]]],variables = z[[2]]}, ImplicitRegion[condition,variables]],
    MemberQ[{And,Or,Not,Xor,Equivalent,Implies},head],
      proofSimplify[head @@ (canonical /@ (List @@ z))],
    head === Rule, First[z] -> canonical[Last[z]],
    True,proofSimplify[z]];
  AssociateTo[normalCache,key -> result];result];
zeroTree[z_Association] := If[And @@ (SameQ[#,0] & /@ Values[z]),0,z];
zeroTree[z_List] := If[And @@ (SameQ[#,0] & /@ z),0,z];
(* Amendment 10: identity has priority over every evaluation. Exact
   existential proofs are sufficient symbolic change proofs. Otherwise try
   concrete exact points with a single fixed function for each radial profile;
   failed or agreeing evaluations establish nothing. Bound variables are
   alpha-renamed before point substitution. Sets are tested at an element,
   but statements containing binders are evaluated on free variables only. *)
profileHeads = {Global`deltaProfile,Global`velocityProfile,Global`xiProfile,
  Global`fProfile,Global`rhoBr,Global`muPerp,Global`cGamma};
changeCache = <||>;
realAlgebraicQ[z_] := Which[
  NumberQ[z], TrueQ[Im[z] == 0],
  Head[z] === Symbol, !NumericQ[z] && !MemberQ[profileHeads,z] && !MemberQ[{True,False,Reals,Complexes,Infinity,Indeterminate,Undefined},z],
  MemberQ[{Plus,Times},Head[z]], And @@ (realAlgebraicQ /@ List @@ z),
  Head[z] === Power, realAlgebraicQ[z[[1]]] && (IntegerQ[z[[2]]] || Head[z[[2]]] === Rational),
  True,False];
realFormulaQ[z_] := Which[
  MemberQ[{True,False},z],True,
  MemberQ[{And,Or,Not,Xor,Equivalent,Implies},Head[z]],And @@ (realFormulaQ /@ List @@ z),
  MemberQ[{Equal,Unequal,Less,LessEqual,Greater,GreaterEqual},Head[z]],And @@ (realAlgebraicQ /@ List @@ z),
  Head[z] === Element, z[[2]] === Reals && And @@ (realAlgebraicQ /@ Flatten[{z[[1]] /. Alternatives -> List}]),
  MemberQ[{ForAll,Exists,Inactive[ForAll],Inactive[Exists]},Head[z]],
    And @@ (realFormulaQ /@ Rest[List @@ z]),
  True,False];
realChart[z_] := And @@ DeleteDuplicates[Cases[z,
  Power[base_,exponent_] /; IntegerQ[exponent] || Head[exponent] === Rational :>
    Which[IntegerQ[exponent] && exponent < 0,base != 0,
      IntegerQ[exponent],True,exponent < 0,base > 0,True,base >= 0],{0,Infinity}]];
scalarVariables[z_] := DeleteDuplicates[Cases[z,q_Symbol /;
  !NumericQ[q] && !MemberQ[{True,False,Reals,Complexes,Infinity,Indeterminate,Undefined},q],{0,Infinity}]];
existentialChange[a_,b_,logical_,context_,equality_] := Module[
  {key=HoldComplete[a,b,logical,context],condition,chart,variables,obligation,result,supported,complete,aa,bb,elements},
  If[KeyExistsQ[changeCache,key],Return[changeCache[key]]];
  {aa,bb}=scoped /@ {a,b};
  If[Head[aa]===ImplicitRegion && Head[bb]===ImplicitRegion && Length[aa[[2]]]===Length[bb[[2]]],
    elements=Table[Unique["S9bHarness`element"],{Length[aa[[2]]]}];
    aa=aa[[1]] /. Thread[aa[[2]]->elements];
    bb=bb[[1]] /. Thread[bb[[2]]->elements]];
  condition = comparisonAssumptions[a,b,context] && scoped[context] && Which[
    logicalObjectQ[aa] && logicalObjectQ[bb],Xor[aa,bb],
    SameQ[equality,False],True,True,aa != bb];
  supported = realFormulaQ[condition] &&
    (logical || (Head[a] === ImplicitRegion && Head[b] === ImplicitRegion && Length[a[[2]]]===Length[b[[2]]]) ||
      (StringQ[a] && StringQ[b]) || (realAlgebraicQ[a] && realAlgebraicQ[b]));
  (* A power inside a binder needs its definedness condition inside that
     binder. Leave that full obligation unresolved instead of moving it out. *)
  If[!FreeQ[condition,q_ /; MemberQ[{ForAll,Exists,Inactive[ForAll],Inactive[Exists]},Head[q]] &&
      !FreeQ[q,Power[_,exponent_] /; !IntegerQ[exponent] || exponent < 0]],supported=False];
  (* Only an everywhere-defined polynomial formula can turn failure to
     find a witness into zero. A chart excluding any pole or root boundary
     is sufficient for a witness, never for a global equality proof. *)
  complete = supported && FreeQ[{a,b,condition},Power[_,exponent_] /;
    !IntegerQ[exponent] || exponent < 0];
  chart = If[supported,realChart[{aa,bb,condition}],True];
  variables = witnessVariables[condition && chart];
  obligation = With[{v=variables,c=condition && chart},Inactive[Exists][v,c]];
  result = If[supported,
    With[{v=variables,c=Activate[condition && chart,ForAll|Exists]},
      If[v === {},Resolve[c,Reals],Resolve[Exists[v,c],Reals]]],Missing["UNDECIDED"]];
  result = <|"ExistenceCondition" -> obligation,"ExistenceResult" -> result,
    "Method" -> "exact real quantifier elimination on the common defined domain",
    "SupportedFormula" -> supported,"CompleteRealDecision" -> complete|>;
  AssociateTo[changeCache,key -> result];result];
(* Capture-avoiding witness substitution. Fresh bound names are NEVER part
   of a witness point. Integrals and quantifiers retain their whole domains. *)
witnessCache = <||>; witnessEvaluationCache = <||>;
scoped[z_] := Module[{h=Head[z],v,fresh,rules,args},
  Which[
    AtomQ[z],z,
    MemberQ[{HoldForm,HoldComplete},h],z,
    h===Function && Length[z]>=2 && (ListQ[First[z]] || Head[First[z]]===Symbol),
      v=Flatten[{First[z]}];fresh=Table[Unique["S9bHarness`bound"],{Length[v]}];
      rules=Thread[v->fresh];args=scoped /@ (Rest[List@@z] /. rules);
      h@@Prepend[args,If[ListQ[First[z]],fresh,First[fresh]]],
    MemberQ[{ForAll,Exists,Inactive[ForAll],Inactive[Exists]},h],
      v=Flatten[{First[z]}]; fresh=Table[Unique["S9bHarness`bound"],{Length[v]}];
      rules=Thread[v->fresh]; args=scoped /@ (Rest[List@@z] /. rules);
      h@@Prepend[args,If[ListQ[First[z]],fresh,First[fresh]]],
    h===Inactive[Integrate] && Length[z]===2 && ListQ[z[[2]]] && Length[z[[2]]]===3,
      v=z[[2,1]]; fresh=Unique["S9bHarness`bound"];
      Inactive[Integrate][scoped[z[[1]] /. v->fresh],
        {fresh,scoped[z[[2,2]]],scoped[z[[2,3]]]}],
    h===Inactive[Limit] && Head[z[[2]]]===Rule,
      v=z[[2,1]]; fresh=Unique["S9bHarness`bound"];
      Inactive[Limit][scoped[z[[1]] /. v->fresh],fresh->scoped[z[[2,2]]]],
    h===ImplicitRegion,
      v=z[[2]];fresh=Table[Unique["S9bHarness`bound"],{Length[v]}];
      With[{c=scoped[z[[1]] /. Thread[v->fresh]],q=fresh},ImplicitRegion[c,q]],
    AssociationQ[z],AssociationMap[scoped[z[#]]&,Keys[z]],
    True,Map[scoped,z]]];
witnessVariables[z_] := Select[scalarVariables[z],
  !(Context[#]==="S9bHarness`" && StringStartsQ[SymbolName[#],"bound"]) &&
  Context[#]=!="System`" && !MemberQ[profileHeads,#] &];
comparisonAssumptions[a_,b_,context_] := Module[{variables,clauses},
  variables=witnessVariables[scoped[{a,b,context}]];
  clauses=If[Head[$Assumptions]===And,List@@$Assumptions,{$Assumptions}];
  And@@Select[clauses,Complement[scalarVariables[#],variables]==={} &]];
(* Four finite test configurations; no fitted jets or independent derivative
   values. These entire rational functions are smooth on the real line.
   They are witness-only and never enter an engine run or a baseline payload. *)
witnessProfiles[index_] := Module[{d,v,x,f,rho},
  {d,v,x,f,rho}=Switch[index,
    1,{Function[{q},0],Function[{q},0],Function[{q},0],Function[{q},0],Function[{q},1]},
    2,{Function[{q},1/(10(1+q^2))],Function[{q},1/(5(1+q^2))],
       Function[{q},q/(10(1+q^2))],Function[{q},1/(3(1+q^2))],Function[{q},1+1/(1+q^2)]},
    3,{Function[{q},-1/(5(1+q^2))],Function[{q},-1/(3(1+q^2))],
       Function[{q},1/(7(1+q^2))],Function[{q},-2],Function[{q},2+q^2]},
    4,{Function[{q},1/4],Function[{q},1/3],Function[{q},q/5],
       Function[{q},-1],Function[{q},1+q^2]}];
  {Global`deltaProfile->d,Global`velocityProfile->v,Global`xiProfile->x,
    Global`fProfile->f,Global`rhoBr->rho,
    Global`muPerp->Function[{q},Evaluate[4 rho[q](1+d[q])^2]],
    Global`cGamma->Function[{q},Evaluate[2(1+d[q])]]}];
witnessPoint[variables_,index_] := Map[Function[var,var->Switch[SymbolName[var],
  "c0",2,"bFar",1,"b"|"bFixed"|"impact"|"impactRadar"|"rTurnFlyby"|"rTurnRadar",2,
  "r"|"rr"|"rAnchor",3,"ZE",4,"ZR",5,"rho0"|"m"|"ell",1,
  "n",2,"s",If[index==3,-1,2],_,{0,1,-1,2}[[index]]]],variables];
(* Exact evaluation only. Rational integrals are attempted in full, with
   their endpoints; unsupported analytic integrals stay inert. There is no
   quadrature, endpoint cutoff, tolerance or numerical error estimate. *)
exactWitnessEvaluation[z_] := Module[{key=HoldComplete[z],result,integrate},
  If[KeyExistsQ[witnessEvaluationCache,key],Return[witnessEvaluationCache[key]]];
  integrate[q_,{v_,lo_,hi_}] := Module[{t=Together[q]},
    If[PolynomialQ[Numerator[t],v] && PolynomialQ[Denominator[t],v] &&
       FreeQ[t, _Real] && FreeQ[{lo,hi},Indeterminate|Undefined|ComplexInfinity],
      Integrate[t,{v,lo,hi},GenerateConditions->True,Assumptions->True],
      Inactive[Integrate][q,{v,lo,hi}]]];
  result=Check[
  result=z /. HoldPattern[Inactive[Integrate][q_,range:{_,_,_}]] :> integrate[q,range];
  result=result /. HoldPattern[Inactive[Integrate][0,___]]->0;
  If[realFormulaQ[result],result=Resolve[Activate[result,ForAll|Exists],Reals]];
  If[FreeQ[result,Inactive|Integrate|ForAll|Exists|ImplicitRegion] &&
     (NumericQ[result] || logicalObjectQ[result]),result=FullSimplify[result,Assumptions->True]];
  result,$Failed,{Power::infy,Infinity::indet,Integrate::idiv}];
  If[result===$Failed,Return[$Failed]];
  AssociateTo[witnessEvaluationCache,key->result];result];
exactValueQ[z_] := MemberQ[{True,False},z] || StringQ[z] ||
  (NumericQ[z] && Precision[z]===Infinity && FreeQ[z,Indeterminate|Undefined|DirectedInfinity]);
(* Undefined arithmetic rejects that candidate only. Keep the candidate's
   work local so a failed Check cannot leave a previously assigned witness.
   Only mathematical singularity/divergence messages are caught; programming
   errors remain visible. No result carrying such a message enters the cache. *)
tryWitness[aa_,bb_,cc_,variables_,elements_,index_] := Module[
  {profiles=witnessProfiles[index],point=witnessPoint[variables,index],
   sa,sb,sc,condition,findCondition,instance,prepared,av,bv,cv,witness},
  {sa,sb,sc}={aa,bb,cc}/.profiles;
  condition=sc && If[logicalObjectQ[sa] && logicalObjectQ[sb],Xor[sa,sb],sa!=sb];
  findCondition=condition && realChart[{sa,sb,sc}];
  If[variables=!={} && !MemberQ[{True,False},findCondition] &&
     realFormulaQ[findCondition] &&
     Complement[scalarVariables[findCondition],variables]==={} &&
     FreeQ[findCondition,_DirectedInfinity|Indeterminate|Undefined|_Real] &&
     FreeQ[findCondition,
      q_ /; MemberQ[{ForAll,Exists,Inactive[ForAll],Inactive[Exists]},Head[q]]],
    instance=FindInstance[findCondition,variables,Reals,1];
    If[ListQ[instance] && Length[instance]>0,point=Join[First[instance],
      Select[point,!MemberQ[First/@First[instance],First[#]]&]]]];
  prepared={sa,sb,sc}/.point;
  cv=exactWitnessEvaluation[prepared[[3]]];
  If[TrueQ[cv],
    av=exactWitnessEvaluation[prepared[[1]]];bv=exactWitnessEvaluation[prepared[[2]]];
    If[exactValueQ[av] && exactValueQ[bv] &&
        TrueQ[If[MemberQ[{True,False},av] && MemberQ[{True,False},bv],Xor[av,bv],
          If[StringQ[av] && StringQ[bv],av=!=bv,FullSimplify[av!=bv]]]],
      witness=<|"Point"->Select[point,!MemberQ[elements,First[#]]&],
        "SetElement"->(elements/.point),"Profiles"->profiles,
        "CorruptedEvaluation"->av,"BaselineEvaluation"->bv,"CommonDomainEvaluation"->cv,
        "Arithmetic"->"exact","ErrorBound"->0,
        "Admissibility"->"fixed rational radial functions, smooth on the real line; all derivatives obtained from those functions"|>;
      Return[<|"Witness"->witness,"Attempt"-><|"Configuration"->index,"Result"->"witness"|>|>]]];
  <|"Witness"->Missing["NO_WITNESS_FOUND"],"Attempt"-><|"Configuration"->index,
    "Point"->point,"CommonDomainEvaluation"->cv,"Result"->"no certified difference"|>|>];
witnessSearch[a_,b_,logical_,context_] := Module[
  {key=HoldComplete[a,b,logical,context],aa,bb,cc,variables,
   attempts={},record,result,elements={},candidate},
  If[KeyExistsQ[witnessCache,key],Return[witnessCache[key]]];
  {aa,bb,cc}=scoped /@ {a,b,comparisonAssumptions[a,b,context] && context};
  If[Head[aa]===ImplicitRegion && Head[bb]===ImplicitRegion && Length[aa[[2]]]===Length[bb[[2]]],
    elements=Table[Unique["S9bHarness`element"],{Length[aa[[2]]]}];
    aa=aa[[1]] /. Thread[aa[[2]]->elements];
    bb=bb[[1]] /. Thread[bb[[2]]->elements]];
  variables=witnessVariables[{aa,bb,cc}];
  result=Missing["NO_WITNESS_FOUND"];
  Do[
    candidate=Quiet[Check[tryWitness[aa,bb,cc,variables,elements,index],
      $Failed,{Power::infy,Infinity::indet,Integrate::idiv}],{Power::infy,Infinity::indet,Integrate::idiv}];
    If[candidate===$Failed,
      AppendTo[attempts,<|"Configuration"->index,
        "Result"->"undefined candidate arithmetic or divergent integral"|>],
      AppendTo[attempts,candidate["Attempt"]];
      If[AssociationQ[candidate["Witness"]],result=candidate["Witness"];Break[]]],{index,4}];
  record=<|"Witness"->result,"Attempts"->attempts,
    "Search"->"four fixed smooth profile configurations and exact scalar points; algebraic FindInstance where supported"|>;
  AssociateTo[witnessCache,key->record];record];
valueContext[paths_,path_] := And @@ Cases[paths,{anchor_,condition_} /;
  Length[anchor] <= Length[path] && Take[path,Length[anchor]] === anchor :> condition];
leafDifference[a_,b_,asPredicate_:False,context_:True] := Module[{na,nb,logical,proof,residual,change,search},
  na = canonical[a]; nb = canonical[b];
  If[SameQ[na,nb],Return[0]];
  logical = (logicalObjectQ[na] && logicalObjectQ[nb]) ||
    (asPredicate && (logicalObjectQ[na] || Head[na] === Symbol) &&
      (logicalObjectQ[nb] || Head[nb] === Symbol));
  proof = proofSimplify[If[logical,Equivalent[na,nb],na == nb]];
  If[TrueQ[proof],Return[0]];
  change = existentialChange[na,nb,logical,context,proof];
  If[TrueQ[change["CompleteRealDecision"]] && SameQ[change["ExistenceResult"],False],Return[0]];
  search = witnessSearch[a,b,logical,context];
  residual = If[logical,Boole[a]-Boole[b],
    If[NumericQ[a] && NumericQ[b],a-b,Inactive[Subtract][a,b]]];
  <|"Status" -> If[TrueQ[change["ExistenceResult"]] || AssociationQ[search["Witness"]],"CHANGED","UNDECIDED"],
    "Equality" -> proof,"ChangeProof" -> change,"WitnessSearch" -> search,
    "Residual" -> residual,"Corrupted" -> a,"Baseline" -> b|>];
partDifference[a_,b_,path_:{},paths_:{}] := Which[
  SameQ[a,b],0,
  AssociationQ[a] && AssociationQ[b],zeroTree[AssociationMap[
    partDifference[Lookup[a,#,Missing["Unemitted"]],Lookup[b,#,Missing["Unemitted"]],Append[path,#],paths] &,
    Union[Keys[a],Keys[b]]]],
  ListQ[a] && ListQ[b] && Length[a] == Length[b],zeroTree[Table[
    partDifference[a[[i]],b[[i]],Append[path,i],paths],{i,Length[a]}]],
  Head[a] === Rule && Head[b] === Rule && SameQ[First[a],First[b]],
    With[{d = partDifference[Last[a],Last[b],Append[path,2],paths]},If[SameQ[d,0],0,First[a] -> d]],
  True,leafDifference[a,b,False,valueContext[paths,path]]];
attachedDifference[a_,b_] := If[SameQ[a,b],0,
  If[AssociationQ[a] || ListQ[a] || AssociationQ[b] || ListQ[b],partDifference[a,b],
    leafDifference[a,b,True]]];
domainDifference[a_Association,b_Association] := If[SameQ[a,b],0,
  zeroTree[AssociationMap[attachedDifference[Lookup[a,#,Missing["Unemitted"]],
    Lookup[b,#,Missing["Unemitted"]]] &,Union[Keys[a],Keys[b]]]]];
difference[a_,b_] := Module[{partsA = splitPayload[a],partsB = splitPayload[b]},
  <|"ValueDifference" -> partDifference[partsA["Value"],partsB["Value"],{},
      Join[partsA["DomainPaths"],partsB["DomainPaths"]]],
    "DomainDifference" -> domainDifference[partsA["Domain"],partsB["Domain"]]|>];
(* END_DOMAIN_DIFFERENCE_METHOD *)
run[path_, name_] := Module[{capture, result},
  ClearSystemCache[]; S9bHarness`normalCache = <||>; S9bHarness`changeCache = <||>; S9bHarness`witnessCache = <||>; S9bHarness`witnessEvaluationCache = <||>; note["engine begin", name];
  capture = OpenWrite[FileNameJoin[{work, name <> ".stdout"}]];
  Block[{$Context = "Global`", $ContextPath = {"Global`", "System`"}, $Output = {capture}}, Get[path]];
  Close[capture]; result = Global`emitted; note["engine end", name]; result];
mutate[name_] := Module[{pair = knives[name], target},
  print[name <> "_SITE", <|"Construction" -> First[pair], "Mutation" -> Last[pair],
    "Occurrences" -> StringCount[source, First[pair]]|>];
  If[StringCount[source, First[pair]] != 1, Quit[92]];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, StringReplace[source, Rule @@ pair], "Text"]; target];
baseline = run[engine, "BASELINE"];
print["MANIFEST", knives];
Do[copy = mutate[knife]; corrupted = run[copy, knife];
  Do[note["difference begin", knife <> "_" <> tag]; print[knife <> "_" <> tag, <|"Baseline" -> Lookup[baseline, tag, Missing["Unemitted"]],
    "Corrupted" -> Lookup[corrupted, tag, Missing["Unemitted"]],
    "Difference" -> difference[Lookup[corrupted, tag, Missing["Unemitted"]],
      Lookup[baseline, tag, Missing["Unemitted"]]]|>];
    note["difference end", knife <> "_" <> tag], {tag, Union[Keys[baseline], Keys[corrupted]]}],
  {knife, If[scope == "K11", {}, DeleteCases[Keys[knives], "K11"]]}];
(* Compact evaluation is exact symbolic prefix evaluation. No points, seeds,
   finite precision or extra grade truncation. Both copies have the same
   extraction boundary. Later tags are explicitly unreached. *)
copy = mutate["K11"];
(* Fully qualified harness state survives ClearAll["Global`*"]. Prepare
   and validate BOTH extracts before executing either. No fallback to a full
   source is possible when the marker is missing or extraction is malformed. *)
S9bHarness`compactBoundary = "(* COMPACT_K11_BOUNDARY *)";
S9bHarness`prepareCompact[path_, name_] := Module[{text, target, position, prefix},
  text = Import[path, "Text"];
  If[!StringQ[S9bHarness`compactBoundary] ||
    StringCount[text, S9bHarness`compactBoundary] != 1, Quit[93]];
  position = First[First[StringPosition[text, S9bHarness`compactBoundary]]];
  prefix = StringTake[text, position - 1];
  If[!StringQ[prefix] || StringLength[prefix] == 0 ||
    StringLength[prefix] >= StringLength[text] ||
    StringContainsQ[prefix, S9bHarness`compactBoundary], Quit[93]];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, prefix, "Text"];
  <|"Path" -> target, "Text" -> prefix|>];
S9bHarness`compactInputs = {
  S9bHarness`prepareCompact[engine, "K11_COMPACT_BASELINE"],
  S9bHarness`prepareCompact[copy, "K11_COMPACT_CORRUPTED"]};
If[S9bHarness`compactInputs[[2]]["Text"] =!=
  StringReplace[S9bHarness`compactInputs[[1]]["Text"], Rule @@ knives["K11"]], Quit[93]];
compactBaseline = run[S9bHarness`compactInputs[[1]]["Path"], "K11_COMPACT_BASELINE"];
compactCorrupted = run[S9bHarness`compactInputs[[2]]["Path"], "K11_COMPACT_CORRUPTED"];
coverage = <|"Method" -> "exact symbolic mechanical prefix extraction",
  "Boundary" -> S9bHarness`compactBoundary, "RetainedGrades" -> Global`grades,
  "LocalConstructions" -> "untruncated dispersion/Fermat and retained-ring nonreciprocal one-form",
  "AdditionalTruncation" -> None, "SampleDomain" -> "symbolic constructor domains",
  "Arithmetic" -> "exact", "Samples" -> {}, "Seeds" -> {}|>;
Do[note["difference begin", "K11_" <> tag]; print["K11_" <> tag, If[KeyExistsQ[compactBaseline, tag] && KeyExistsQ[compactCorrupted, tag],
  <|"Baseline" -> baseline[tag], "CompactBaseline" -> compactBaseline[tag],
    "CompactCorrupted" -> compactCorrupted[tag],
    "Difference" -> difference[compactCorrupted[tag], compactBaseline[tag]],
    "MethodResidual" -> difference[compactBaseline[tag], baseline[tag]], "Coverage" -> coverage|>,
  <|"Baseline" -> baseline[tag], "CompactBaseline" -> Missing["NOT_EVALUATED"],
    "CompactCorrupted" -> Missing["NOT_EVALUATED"], "Difference" -> Missing["NOT_EVALUATED"],
    "MethodResidual" -> Missing["NOT_EVALUATED"], "Coverage" -> coverage|>]]; note["difference end", "K11_" <> tag], {tag, Keys[baseline]}];
(* Harness dead-path ablation: both inputs are the actual baseline. *)
If[scope == "ALL", Do[note["difference begin", "HARNESS_DEAD_PATH_" <> tag]; print["HARNESS_DEAD_PATH_" <> tag, <|"Baseline" -> baseline[tag],
 "Corrupted" -> baseline[tag], "Difference" -> difference[baseline[tag], baseline[tag]]|>]; note["difference end", "HARNESS_DEAD_PATH_" <> tag], {tag, Keys[baseline]}]];
note["complete", scope]; Close[S9bHarness`progress];
End[];
