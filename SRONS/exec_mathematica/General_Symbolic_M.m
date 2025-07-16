(* ::Package:: *)

(* ::Input::Initialization:: *)
(**define some values**)
param1 = 2 + dim; (*number of parameters in one mode*)
param2 = 2*param1;

(**define parameters for two different modes**)
qgen[t_]=Table[Symbol["q"<>ToString[i]][t],{i,1,param2}]; 
(*setup spatial coordinates in vector*)
xmat=Table[Symbol["x"<>ToString[i]],{i,1,dim}];

(**Assumptions on paramter values, this can make integration significantly faster**)
(*Params and derivs are real*)
Do[$Assumptions=$Assumptions&&qgen[t][[i]]\[Element]Reals&&D[qgen[t][[i]],t]\[Element]Reals,{i,1,param2}] 
(*weights are strictly positive*)
Do[$Assumptions=$Assumptions&&qgen[t][[i]]>0,{i,2,2+param1,param1}]
(*spatial coordinates are real*)
Do[$Assumptions=$Assumptions&&xmat[[i]]\[Element]Reals,{i,1,dim}]

assumps = $Assumptions

(**define approximate solution, use qi[t] for ith parameter**)
(*First define |x - x_i(t)|^2 for both modes*)
centersl = 0;
centersj = 0;
Do[centersl = centersl + (xmat[[i]] - qgen[t][[i+2]])^2;
centersj = centersj + (xmat[[i]] - qgen[t][[i+2+param1]])^2,{i,1,dim}]
If[posamp ==0,
(
uhatl =qgen[t][[1]]*ToExpression[phi<>"["<>ToString[qgen[t][[2]]^-2*(centersl),InputForm]<>"]"];uhatj = qgen[t][[1+param1]]*ToExpression[phi<>"["<>ToString[qgen[t][[2+param1]]^-2*(centersj),InputForm]<>"]"];
),
(
uhatl = qgen[t][[1]]^2*ToExpression[phi<>"["<>ToString[-qgen[t][[2]]^-2*(centersl),InputForm]<>"]"];uhatj = qgen[t][[1+param1]]^2*ToExpression[phi<>"["<>ToString[-qgen[t][[2+param1]]^-2*(centersj),InputForm]<>"]"];
)
]

uhat = uhatl + uhatj;

(**Initializing stuff for metric tensor**)
numints = (param1+1)*(param1)/2; (*number of ints for metric tensor*)
Mblock=Table[ij,{i,param1},{j,param1}];dpdqgen = Table[0,{i,1,param2}]; (*store all derivatives*)
Do[dpdqgen[[i]] = D[uhat,qgen[t][[i]]],{i,1,param2}]

DistributeDefinitions[qgen,dpdqgen,xmat,dim,param1,uhatl,uhatj,uhat,assumps];
SetSharedVariable[Mblock]
ParallelEvaluate[$Assumptions = assumps]


(* Calculate metric tensor *)
Print["Calculating "<>ToString[numints]<>" integrals for M"]
ParallelDo[integrand =Collect[Expand[dpdqgen[[i]]*dpdqgen[[j]]],{uhatl,uhatj},Simplify];
Do[integrand = Expand[integrand,xmat[[k]]],{k,1,dim}];
Do[integrand = Integrate[integrand,{xmat[[k]],lb,ub}],{k,1,dim}];
Mblock[[i,j]] = Simplify[integrand];
Print[ToString[i]<>","<>ToString[j]],{i,1,param1},{j,1,i}]
intMTable = Simplify[intMTable];

(**Build upper triangular part by flipping all parameter values**)
genvarreplace = Table[qgen[t][[i]]->Symbol["newqgen"<>ToString[i]][t],{i,1,param1}];
genvarreplace = Join[genvarreplace,Table[qgen[t][[i]]->qgen[t][[i-param1]],{i,param1+1,param2}]];
genvarreplace2= Table[Symbol["newqgen"<>ToString[i]][t]->qgen[t][[i+param1]],{i,1,param1}];

DistributeDefinitions[genvarreplace,genvarreplace2];

ParallelDo[ Part[Mblock,j,i ]=Part[Mblock,i,j ]/.genvarreplace/.genvarreplace2,{i,1 ,param1},{j,1,i-1},DistributedContexts->None]

(*SetDirectory[NotebookDirectory[]]
DumpSave["GeneralInts_Symbolic_M_"<>pdename<>".mx","Global`"]*)





