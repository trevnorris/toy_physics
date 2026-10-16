(* Executable syntax/structural demonstrations from the actual harness.
   No expected physics value, canonicalization or verdict. *)
$HistoryLength=0;
$Messages={OutputStream["stderr",2]};
src=Import["research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl","Text"];
syntax=Block[{$Context="S9bSyntax`",$ContextPath={"System`"}},ToExpression[src,InputForm,HoldComplete]];
Print["HARNESS_PARSE_HEAD: ",Head[syntax]];
method=StringSplit[StringSplit[src,"(* BEGIN_LITERAL_METHOD *)"][[2]],"(* END_LITERAL_METHOD *)"][[1]];
Begin["S9bHarness`"];
Block[{$Context="S9bHarness`",$ContextPath={"System`"}},ToExpression[Global`method]];
End[];
show[name_,a_HoldComplete,b_HoldComplete] := Print[name,": ",
 ToString[<|"Operands"->{a,b},"Comparison"->S9bHarness`comparison[a,b]|>,InputForm,PageWidth->Infinity]];
show["SAME_SYNTAX",HoldComplete[x+1],HoldComplete[x+1]];
show["POLYNOMIAL_FORMS",HoldComplete[(x+1)^2],HoldComplete[x^2+2 x+1]];
show["INERT_ZERO_INTEGRAL",HoldComplete[Inactive[Integrate][0,{r,1,Infinity}]],HoldComplete[0]];
show["DIFFERENT_APPROXIMATE_NUMBERS",HoldComplete[1.],HoldComplete[1.0000000000000002]];
show["DIFFERENT_PRECISIONS",HoldComplete[1.0`20],HoldComplete[1.0`30]];
show["DIFFERENT_DOMAINS",HoldComplete[ConditionalExpression[x,x>0]],HoldComplete[ConditionalExpression[x,x>1]]];
show["ALPHA_RENAMED_BINDERS",HoldComplete[Inactive[ForAll][r,r>0]],HoldComplete[Inactive[ForAll][x,x>0]]];
show["NO_EVALUATION",S9bHarness`freeze["1+1"],S9bHarness`freeze["2"]];
show["ORDERED_CASES",HoldComplete[Piecewise[{{1,x>0},{2,x>1}},3]],
 HoldComplete[Piecewise[{{2,x>1},{1,x>0}},3]]];
(* Arithmetic fixtures exercise the actual witness bookkeeping definitions.
   These are synthetic objects, not additional engine physics inputs. *)
witnessMethod=StringSplit[StringSplit[src,"(* BEGIN_WITNESS_METHOD *)"][[2]],"(* END_WITNESS_METHOD *)"][[1]];
Block[{$Context="S9bHarness`",$ContextPath={"System`"}},ToExpression[Global`witnessMethod]];
$Assumptions=True;
fixture1=<|"State"-><|"RadialMetric"->1|>|>;
fixture2=<|"State"-><|"RadialMetric"->2|>|>;
Print["COMMON_FIXED_SOLUTION_SETS: ",ToString[S9bHarness`fixedReading[{Global`fLocal},fixture1,fixture1],InputForm,PageWidth->Infinity]];
Print["DISJOINT_FIXED_SOLUTION_SETS: ",ToString[S9bHarness`fixedReading[{Global`aMetric},fixture1,fixture2],InputForm,PageWidth->Infinity]];
Print["ORDERED_CASE_SELECTION: ",ToString[S9bHarness`caseReadings[
 Piecewise[{{11,Global`n==0},{12,Global`n==1}},13],fixture1],InputForm,PageWidth->Infinity]];
