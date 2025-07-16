(* ::Package:: *)

(* ::Input::Initialization:: *)
(**Build RHS vector f**)

(*Parameters which appear on RHS of PDE*)
params = {alpha1, alpha2,alpha3,nu};

(*Add assumptions to parameters for PDE, e.g. if parameters are real valued or strictly positive. This makes integration easier*)
$Assumptions=$Assumptions&&nu\[Element]Reals&& nu>=0 &&sigma\[Element]Reals;
$Assumptions=$Assumptions&&alpha1\[Element]Reals&&alpha2\[Element]Reals&&alpha3\[Element]Reals;
ordernonlin = 1; (*Highest order nonlinearity of the RHS. For example:
 1th order (i.e. linear) -- u or D[u,{x,1}] or Laplacian[u]
2st order -- u^2 or u*D[u,{x,1}] or u*D[u,{x,2}]
3nd order -- u^3 or u*(D[u,{x,1}])^2 or u*D[u,{x,2}]*D[u,{y,1}]
etc. *)

(** Build extra terms for nonlinearities, don't need to change this**)
qgenk[t_]=Table[Symbol["qk"<>ToString[i]<>ToString[j]][t],{i,1,ordernonlin-1},{j,1,param1}];

centersk = Table[0,{i,1,ordernonlin-1}];
Do[centersk[[i]] = centersk[[i]] + (xmat[[j]] - qgenk[t][[i,j+2]])^2,{i,1,ordernonlin-1},{j,1,dim}]

uhatk = Table[0,{i,1,ordernonlin-1}];
Do[If[posamp ==0,
(
uhatk[[i]] = qgenk[t][[i,1]]*Exp[-qgenk[t][[i,2]]^-2*(centersk[[i]])];
),
(
uhatk[[i]] = qgenk[t][[i,1]]^2*Exp[-qgenk[t][[i,2]]^-2*(centersk[[i]])];
)
],{i,1,ordernonlin-1}]
(****)

(** Enter RHS of PDE. Some examples are provided below which appear in Anderson and Farazmand 2023 preprint: 
Fisher information and shape-morphing modes for solving the Fokker-Planck equation in higher dimensions **)

(*1D Bistable Potential*)
v = alpha1*xmat[[1]]^4 + alpha2*xmat[[1]]^3 + alpha3*xmat[[1]]^2;
vprime = D[v,xmat[[1]]];
F0=D[vprime*uhatj,xmat[[1]]]+nu*D[uhatj,{xmat[[1]],2}];

(*Stochastic Duffing Oscillator*)
(*F0=-( xmat[[2]] *D[uhatj,{xmat[[1]],1}] + alpha2*uhatj + (alpha1*xmat[[1]] + alpha2*xmat[[2]]+alpha3*xmat[[1]]^3)* D[uhatj,{xmat[[2]],1}]) + nu*D[uhatj,{xmat[[2]],2}]*)

(*****************************)

(* In general, this is how we define terms to apply S-RONS:

F0 terms -- uhatj
F1 terms -- uhatl and uhatj
F2 terms -- uhatl and uhatj and uhatk[[1]]
F3 terms -- uhatl and uhatj and uhatk[[1]] and uhatk[[2]]
F4 terms -- uhatl and uhatj and uhatk[[1]] and uhatk[[2]] and uhatk[[3]]
...;

For example:
;

If RHS of PDE is nu*Laplacian[u], then
F0 = nu*Simplify[Laplacian[uhatj,xmat]]; 

If RHS of PDE is u^2, then
F1 = uhatl*uhatj; 
or, if RHS of PDE is u*D[u,x], then
F1 = uhatl*D[uhatj,xmat[[1]]]; or we could also define
F1 = uhatj*D[uhatl,xmat[[1]]];

If RHS of PDE is u*D[u,x]*D[u,y], then
F2 = uhatl*D[uhatj,xmat[[1]]]*D[uhatk[[1]],xmat[[2]]]; (*u^3*)

If RHS of PDE is u^4, then
F3 = uhatl*uhatj*uhatk[[1]]*uhatk[[2]];
etc.*)

(*NOTE: must enter Fi = 0 if there are higher order nonlinear terms with no lower order terms. For example: if RHS is F(u) = u^2, then F0 = 0, F1 = uhatl*uhatj*)


(* ::Input::Initialization:: *)
(**Don't have to change anything in the code below this**)

(**Initializing stuff**)
(*store ints for RHS in table*)
intFTable = Table[ToString[i]<>ToString[j],{i,1,ordernonlin},{j,1,param1}]; 
(*derivs*)
dpdqgenk = Table[0,{i,1,ordernonlin-1},{j,1,param1}]; (*storage for derivatives*)
Do[dpdqgenk[[i,j]] = D[uhatk[[i]],qgenk[t][[i,j]]],{i,1,ordernonlin - 1},{j,1,param1}]


(** Update assumptions for new parameters **)
(*Params and derivs are real*)
Do[$Assumptions=$Assumptions&&qgenk[t][[i,j]]\[Element]Reals&&D[qgenk[t][[i,j]],t]\[Element]Reals,{i,1,ordernonlin - 1},{j,1,param1}] 
(*weights are strictly positive*)
Do[$Assumptions=$Assumptions&&qgenk[t][[i,2]]>0,{i,1,ordernonlin - 1}]
assumps = $Assumptions;

SetSharedVariable[intFTable,F0,F1,F2,F3,F4,F5,F6,F7]; (*Need to add more F values if higher than 6th order nonlinearity (gross)*)
ParallelEvaluate[$Assumptions = assumps];

(** Calculate RHS vector **)

ivals = Table[i,{i,1,param1}];
ivals =  Delete[ivals,2]; (* we will calculate the derivatives with respect to length scales seperately because it is faster*)

(*Linear terms*)
Print["Calculating "<>ToString[param1]<>" integrals for linear terms of F"]
ParallelDo[integrand = Collect[Expand[Expand[dpdqgen[[i]]*F0]/.{Exp[p_]->c*Exp[-a*Sum[(xmat[[k]] -  xc[t][[k]])^2 ,{k,1,dim}]]}/.xreplace],Join[{Exp[-a*Sum[xtilde[[k]]^2,{k,1,dim}]]},xtilde],Simplify];
(*integrand = Together[integrand/Exp[-a*(Simplify[Norm[xtilde]]^2)]];*)
integrand = integrand/.Exp[p_]->1;
integrand = Simplify[(Pi/a)^(dim/2)*Expectation[integrand,xtilde\[Distributed]MultinormalDistribution[Table[0,{l,1,dim}],DiagonalMatrix[Table[1/a/2,{l,1,dim}]]]]];
intFTable[[1,i]] = integrand;
Print[i],{i,ivals}]

(*integrate derivatives w.r.t. length scales*)
lints = Table[ToString[j],{j,1,dim}]; 
SetSharedVariable[lints]; 
ParallelDo[integrand = Collect[Expand[Expand[uhatl*(xmat[[i]]-qgen[t][[2+i]])^2*F0]/.{Exp[p_]->c*Exp[-a*Sum[(xmat[[k]] -  xc[t][[k]])^2 ,{k,1,dim}]]}/.xreplace],Join[{Exp[-a*Sum[xtilde[[k]]^2,{k,1,dim}]]},xtilde],Simplify];
(*integrand = Together[integrand/Exp[-a*(Simplify[Norm[xtilde]]^2)]];*)
integrand = integrand/.Exp[p_]->1;
integrand = Simplify[(Pi/a)^(dim/2)*Expectation[integrand,xtilde\[Distributed]MultinormalDistribution[Table[0,{l,1,dim}],DiagonalMatrix[Table[1/a/2,{l,1,dim}]]]]];
lints[[i]] = integrand,{i,1,dim}]

intFTable[[1,2]] = Simplify[Sum[lints[[i]],{i,1,dim}]]*2/qgen[t][[2]]^3;
Print["2"]

(*ParallelDo[integrand = Collect[Expand[Expand[dpdqgen[[2+i]]*(xmat[[i]]-qgen[t][[2+i]])/qgen[t][[2]]*F0]/.{Exp[p_]->c*Exp[-a*Sum[(xmat[[k]] -  xc[t][[k]])^2 ,{k,1,dim}]]}/.xreplace],Join[{Exp[-a*Sum[xtilde[[k]]^2,{k,1,dim}]]},xtilde],Simplify];
(*integrand = Together[integrand/Exp[-a*(Simplify[Norm[xtilde]]^2)]];*)
integrand = integrand/.Exp[p_]->1;
integrand = Simplify[(Pi/a)^(dim/2)*Expectation[integrand,xtilde\[Distributed]MultinormalDistribution[Table[0,{l,1,dim}],DiagonalMatrix[Table[1/a/2,{l,1,dim}]]]]];
intFTable[[1,2]] = intFTable[[1,2]]+integrand;
Print[i],{i,1,dim}]*)


(*Nonlinear terms*)
Do[Which[iord==1,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"st order nonlinear terms of F"],
iord==2,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"nd order nonlinear terms of F"],
iord==3,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"rd order nonlinear terms of F"],
iord>3,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"th order nonlinear terms of F"]
];ParallelDo[integrand = Collect[Expand[Expand[dpdqgenk[[iord,i]]*ToExpression["F"<>ToString[iord]]]/.{Exp[p_]->c*Exp[-a*Sum[(xmat[[k]] -  xc[t][[k]])^2 ,{k,1,dim}]]}/.xreplace],Join[{Exp[-a*Sum[xtilde[[k]]^2,{k,1,dim}]]},xtilde],Simplify];
(*integrand = Together[integrand/Exp[-a*(Simplify[Norm[xtilde]]^2)]];*)
integrand = integrand/.Exp[p_]->1;
integrand = Simplify[(Pi/a)^(dim/2)*Expectation[integrand,xtilde\[Distributed]MultinormalDistribution[Table[0,{l,1,dim}],DiagonalMatrix[Table[1/a/2,{l,1,dim}]]]]];
intFTable[[iord+1,i]] = integrand;
Print[i],{i,ivals}];

ParallelDo[integrand = Collect[Expand[Expand[uhatk[[iord]]*(xmat[[i]]-qgenk[t][[iord,2+i]])^2*ToExpression["F"<>ToString[iord]]]/.{Exp[p_]->c*Exp[-a*Sum[(xmat[[k]] -  xc[t][[k]])^2 ,{k,1,dim}]]}/.xreplace],Join[{Exp[-a*Sum[xtilde[[k]]^2,{k,1,dim}]]},xtilde],Simplify];
integrand = integrand/.Exp[p_]->1;
integrand = Simplify[(Pi/a)^(dim/2)*Expectation[integrand,xtilde\[Distributed]MultinormalDistribution[Table[0,{l,1,dim}],DiagonalMatrix[Table[1/a/2,{l,1,dim}]]]]];
lints[[i]] = integrand,{i,1,dim}];
intFTable[[iord+1,2]] = Simplify[Sum[lints[[i]],{i,1,dim}]]*2/qgenk[t][[2]]^3;
Print["2"],{iord,1,ordernonlin - 1}]


intFTable = Collect[intFTable,a,Simplify];

(*SetDirectory[NotebookDirectory[]];
DumpSave["GeneralInts_Symbolic_EXP_f_"<>pdename<>".mx","Global`"];*)
