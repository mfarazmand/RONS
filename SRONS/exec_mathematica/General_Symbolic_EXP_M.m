(* ::Package:: *)

(* ::Input::Initialization:: *)
(**define some values**)
param1 = 2 + dim; (*number of parameters in one mode*)
param2 = 2*param1;

(**define parameters for two different modes**)
qgen[t_]=Table[Symbol["q"<>ToString[i]][t],{i,1,param2}]; 
(*setup spatial coordinates and centers in each dimmension*)
xmat=Table[Symbol["x"<>ToString[i]],{i,1,dim}];
xc[t_]=Table[Symbol["xc"<>ToString[i]][t],{i,1,dim}];
xtilde=Table[Symbol["xtilde"<>ToString[i]],{i,1,dim}];
xtildevals=Table[Symbol["xtilde"<>ToString[i]],{i,1,dim}];
xtildevalsF=Table[Symbol["xtildeF"<>ToString[i]],{i,1,dim}];

(**Assumptions on paramter values, this can make integration significantly faster**)

(*Params and derivs are real*)
Do[$Assumptions=$Assumptions&&qgen[t][[i]]\[Element]Reals&&D[qgen[t][[i]],t]\[Element]Reals,{i,1,param2}] 
(*weights are strictly positive*)
Do[$Assumptions=$Assumptions&&qgen[t][[i]]>0,{i,2,2+param1,param1}]
(*spatial coordinates are real*)
Do[$Assumptions=$Assumptions&&xmat[[i]]\[Element]Reals&&xc[t][[i]]\[Element]Reals&&xtilde[[i]]\[Element]Reals,{i,1,dim}]
(*Assumptions for COV we will make later*)
$Assumptions=$Assumptions&&c\[Element]Reals&&c>=1&&a\[Element]Reals&&a>0;

assumps = $Assumptions;

(**define approximate solution, use qi[t] for ith parameter**)
(*First define |x - x_i(t)|^2 for both modes*)
centersl = 0;
centersj = 0;
Do[centersl = centersl + (xmat[[i]] - qgen[t][[i+2]])^2;
centersj = centersj + (xmat[[i]] - qgen[t][[i+2+param1]])^2,{i,1,dim}]
xcentl = xmat - Table[qgen[t][[2+i ]],{i,1,dim}];xcentj = xmat - Table[qgen[t][[2+i + param1]],{i,1,dim}];
If[posamp ==0,
(
uhatl = qgen[t][[1]]*Exp[-qgen[t][[2]]^-2*(centersl)];uhatj = qgen[t][[1+param1]]*Exp[-qgen[t][[2+param1]]^-2*(centersj)];
),
(
uhatl = qgen[t][[1]]^2*Exp[-qgen[t][[2]]^-2*(centersl)];uhatj = qgen[t][[1+param1]]^2*Exp[-qgen[t][[2+param1]]^-2*(centersj)];
)
]
uhat = uhatl + uhatj;

(**Initializing things for Metric tensor**)
numints = (param1+1)*(param1)/2; (*number of ints for metric tensor*)
Mblock=Table[ij,{i,param1},{j,param1}]; (*one block of metric tensor*)
dpdqgen = Table[0,{i,1,param2}]; (*store all derivatives*)Do[dpdqgen[[i]] = D[uhat,qgen[t][[i]]],{i,1,param2}]

(**COV so thats x_i - c_i(t) \[Rule] \tilde{x}_i for each i. Looks messy, but makes integration and final outputs nicer. Doesn't need to be touched**)
xreplace = Table[xmat[[i]] - xc[t][[i]]->xtilde[[i]],{i,1,dim}];
xreplace = Join[xreplace,Table[xmat[[i]] ->xtilde[[i]]+xc[t][[i]],{i,1,dim}] ];

(*Makes COV easier*)
CompleteSquare[f_,x_]:=Module[{a,b,c},{c,b,a}=CoefficientList[f,x];
a* (x+Simplify[b/2/a])^2+Simplify[(c-b^2/4/a)]]
CompleteSquareCoeffs[f_,x_]:=Module[{a,b,c},{c,b,a}=CoefficientList[f,x];
{a,b,c}]

DistributeDefinitions[qgen,dpdqgen,xtilde,xmat,xc,dim,param1,xreplace,uhatl,uhatj,uhat,assumps];
SetSharedVariable[Mblock];
ParallelEvaluate[$Assumptions = assumps];

(* Calculate metric tensor *)
Print["Calculating "<>ToString[numints]<>" integrals for M"]
ParallelDo[integrand = Collect[Expand[Expand[dpdqgen[[i]]*dpdqgen[[j+param1]]]/.{Exp[p_]->c*Exp[-a*Sum[(xmat[[k]] -  xc[t][[k]])^2 ,{k,1,dim}]]}/.xreplace],Join[{Exp[-a*Sum[xtilde[[k]]^2,{k,1,dim}]]},xtilde],Simplify];
(*integrand = Together[integrand/Exp[-a*(Simplify[Norm[xtilde]]^2)]];*)
integrand = integrand/.Exp[p_]->1;
integrand = Simplify[(Pi/a)^(dim/2)*Expectation[integrand,xtilde\[Distributed]MultinormalDistribution[Table[0,{l,1,dim}],DiagonalMatrix[Table[1/a/2,{l,1,dim}]]]]];Mblock[[i,j]] = Collect[integrand,a,Simplify];
Print[ToString[i]<>","<>ToString[j]],{i,1,param1},{j,1,i}]

(**Build upper triangular part by flipping all parameter values**)
genvarreplace = Table[qgen[t][[i]]->Symbol["newqgen"<>ToString[i]][t],{i,1,param1}];
genvarreplace = Join[genvarreplace,Table[qgen[t][[i]]->qgen[t][[i-param1]],{i,param1+1,param2}]];
genvarreplace2= Table[Symbol["newqgen"<>ToString[i]][t]->qgen[t][[i+param1]],{i,1,param1}];

DistributeDefinitions[genvarreplace,genvarreplace2];

ParallelDo[ Part[Mblock,j,i ]=Part[Mblock,i,j ]/.genvarreplace/.genvarreplace2,{i,1 ,param1},{j,1,i-1},DistributedContexts->None]

(*pop up in COV*)
prod = -centersl*qgen[t][[2]]^-2 -centersj*qgen[t][[2+param1]]^-2;
exppart = 0;
Do[{acs,bcs,ccs} = CompleteSquareCoeffs[prod,xmat[[i]]];
xtildevals[[i]] = Simplify[bcs/2/acs];
exppart =exppart +  Simplify[acs]* (xmat[[i]]+Simplify[bcs/2/acs])^2;
prod =  Simplify[(ccs-bcs^2/4/acs)],{i,1,dim}]

(*DumpSave["GeneralInts_Symbolic_EXP_M_"<>pdename<>".mx","Global`"]*)
