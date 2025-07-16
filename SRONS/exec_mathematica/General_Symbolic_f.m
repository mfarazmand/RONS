(* ::Package:: *)

(* ::Input::Initialization:: *)
(**RHS of PDE**)

(*Parameters which appear on RHS of PDE*)
params = {nu, sigma};

(*Add assumptions to parameter values for PDE, makes integration easier*)
$Assumptions=$Assumptions&&nu\[Element]Reals&& nu>=0 &&sigma\[Element]Reals;

ordernonlin = 1; (*Highest order nonlinearity of the RHS. For example:
 1th order (i.e. linear) -- u or D[u,{x,1}] or Laplacian[u]
2st order -- u^2 or u*D[u,{x,1}] or u*D[u,{x,2}]
3nd order -- u^3 or u*(D[u,{x,1}])^2 or u*D[u,{x,2}]*D[u,{y,1}]
etc. *)

(** Build extra terms for nonlinearities, don't need to change this**)
qgenk[t_]=Table[Symbol["qk"<>ToString[i]<>ToString[j]][t],{i,1,ordernonlin - 1},{j,1,param1}];

centersk = Table[0,{i,1,ordernonlin - 1}];
Do[centersk[[i]] = centersk[[i]] + (xmat[[j]] - qgenk[t][[i,j+2]])^2,{i,1,ordernonlin - 1},{j,1,dim}]

uhatk = Table[0,{i,1,ordernonlin - 1}];
Do[If[posamp ==0,
(
uhatk[[i]] =qgenk[t][[i,1]]*ToExpression[phi<>"["<>ToString[qgenk[t][[i,2]]^-2*(centersk[[i]]),InputForm]<>"]"];
),
(
uhatk[[i]] =qgenk[t][[i,1]]^2*ToExpression[phi<>"["<>ToString[qgenk[t][[i,2]]^-2*(centersk[[i]]),InputForm]<>"]"];
)
],{i,1,ordernonlin - 1}]
(** Build extra terms for nonlinearities, don't need to change this**)

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
(**Don't have to change below this**)

(**Initializing stuff**)
(*store ints for RHS*)
intFTable = Table[ToString[i]<>ToString[j],{i,1,ordernonlin},{j,1,param1}]; 
(*derivs*)
dpdqgenk = Table[0,{i,1,ordernonlin - 1},{j,1,param1}]; (*store all derivatives*)
Do[dpdqgenk[[i,j]] = D[uhatk[[i]],qgenk[t][[i,j]]],{i,1,ordernonlin - 1},{j,1,param1}]


(** Update assumptions for new parameters **)
(*Params and derivs are real*)
Do[$Assumptions=$Assumptions&&qgenk[t][[i,j]]\[Element]Reals&&D[qgenk[t][[i,j]],t]\[Element]Reals,{i,1,ordernonlin - 1},{j,1,param1}] 
(*weights are strictly positive*)
Do[$Assumptions=$Assumptions&&qgenk[t][[i,2]]>0,{i,1,ordernonlin - 1}]
assumps = $Assumptions;

SetSharedVariable[intFTable,F0,F1,F2,F3,F4,F5,F6,F7]; (*Need to add more F values if higher than 6th order nonlinearity (ew)*)
ParallelEvaluate[$Assumptions = assumps]

(** Calculate RHS vector **)
(*Linear terms*)
Print["Calculating "<>ToString[param1]<>" integrals for linear terms of F"]
ParallelDo[
integrand =Collect[Expand[dpdqgen[[i]]*F0],{uhatl,uhatj},Simplify];
Do[integrand = Expand[integrand,xmat[[k]]],{k,1,dim}];
Do[integrand = Integrate[integrand,{xmat[[k]],lb,ub}],{k,1,dim}];
intFTable[[1,i]] = integrand;
Print[i],{i,1,param1}]

intFTable = Simplify[intFTable];

(*Nonlinear terms*)
Do[Which[iord==1,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"st order nonlinear terms of F"],
iord==2,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"nd order nonlinear terms of F"],
iord==3,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"rd order nonlinear terms of F"],
iord>3,
Print["Calculating "<>ToString[param1]<>" integrals for " <>ToString[iord] <>"th order nonlinear terms of F"]
];
ParallelDo[
integrand =Collect[Expand[dpdqgenk[[iord,i]]*ToExpression["F"<>ToString[iord]]],{uhatl,uhatj},Simplify];
Do[integrand = Expand[integrand,xmat[[k]]],{k,1,dim}];
Do[integrand = Integrate[integrand,{xmat[[k]],lb,ub}],{k,1,dim}];
intFTable[[iord+1,i]] = integrand;
Print[i],{i,1,param1}],{iord,1,ordernonlin - 1}]

(*SetDirectory[NotebookDirectory[]]
DumpSave["GeneralInts_Symbolic_f_"<>pdename<>".mx","Global`"]*)



