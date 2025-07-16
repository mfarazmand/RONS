(* ::Package:: *)

(* ::Input::Initialization:: *)
(*SetDirectory[NotebookDirectory[]];*)
Get["ToMatlab_SRONS.m"];
Print["Writing files to Matlab"]
dir = NotebookDirectory[]<>pdename;
If[DirectoryQ[dir],SetDirectory[dir],CreateDirectory[dir];
SetDirectory[dir]];
DumpSave["GeneralInts_Symbolic_"<>pdename<>".mx","Global`"]
(*DeleteDirectory[pdename,DeleteContents\[Rule]True]
CreateDirectory[pdename]
SetDirectory[NotebookDirectory[]<>pdename];*)

(** Write M **)
varreplace = Table[qgen[t][[k]]->Symbol["q"<>ToString[k]],{k,1,param2}];
varreplace = Join[varreplace,Table[xc[t][[k]]->Symbol["xc"<>ToString[k]],{k,1,dim}]];

centerstring = "";Do[centerstring = StringJoin[centerstring,"xc"<>ToString[i]<>" = "<>ToMatlab[-xtildevals[[i]]/.varreplace]<>"\n"],{i,1,dim}]

qstring = "";
Do[qstring = qstring<>"q"<>ToString[k]<>"=q("<>ToString[k]<>");",{k,1,param2}]

filename = OpenWrite["Mblock_"<>pdename<>".m"];
WriteMatlabMEXP[Mblock/.varreplace,qstring,prod/.varreplace,centerstring,acs/.varreplace,pdename,filename]; 
Close[filename];

(** Write f **)
intFTableprint = intFTable;

Do[If[i==1,
(
prodk =-centersl*qgen[t][[2]]^-2 -centersj*qgen[t][[2+param1]]^-2;
exppartk = 0;
Do[{acsk,bcsk,ccsk} = CompleteSquareCoeffs[prodk,xmat[[k]]];
xtildevalsF[[k]] = Simplify[bcsk/2/acsk];
exppartk =exppartk +  Simplify[acsk]* (xmat[[k]]+Simplify[bcsk/2/acsk])^2;
prodk =  Simplify[(ccsk-bcsk^2/4/acsk)],{k,1,dim}];

qstring = "";Do[qstring = qstring<>"ql"<>ToString[k]<>"=ql("<>ToString[k]<>");",{k,1,param1}];
qstring = qstring<>"\n \nqj = reshape(q,[param1,r])'; \n \n";
Do[qstring = qstring<>"qj"<>ToString[k]<>"=qj(:,"<>ToString[k]<>");",{k,1,param1}];

qinputs = "ql,q";

varreplaceF = Table[qgen[t][[k]]->"ql"<>ToString[k],{k,1,param1}];
varreplaceF = Join[varreplaceF,Table[qgen[t][[k]]->Symbol["qj"<>ToString[k-param1]],{k,param1+1,2*param1}]];
varreplaceF = Join[varreplaceF,Table[xc[t][[k]]->Symbol["xc"<>ToString[k]],{k,1,dim}]];

),
(

prodk =-centersl*qgen[t][[2]]^-2 -centersj*qgen[t][[2+param1]]^-2 + Sum[centersk[[k]]*-qgenk[t][[k,2]]^-2,{k,1,i-1}];
exppartk = 0;
Do[{acsk,bcsk,ccsk} = CompleteSquareCoeffs[prodk,xmat[[k]]];
xtildevalsF[[k]] = Simplify[bcsk/2/acsk];
exppartk =exppartk +  Simplify[acsk]* (xmat[[k]]+Simplify[bcsk/2/acsk])^2;
prodk =  Simplify[(ccsk-bcsk^2/4/acsk)],{k,1,dim}];

qstring = "";Do[qstring = qstring<>"qk"<>ToString[i-1]<>ToString[k]<>"=qk("<>ToString[i-1]<>","<>ToString[k]<>");";
If [k ==param1,qstring = qstring<>"\n"],{k,1,param1}];
qstring = qstring<>"\n\nq=reshape(q,[param1,r])';\n\n";
Do[qstring = qstring<>"ql"<>ToString[k]<>" = reshape(q(:,"<>ToString[k]<>"),[1,r]);\n",{k,1,param1}];
Do[qstring = qstring<>"qj"<>ToString[k]<>" = ql"<>ToString[k]<>"';",{k,1,param1}];
qstring = qstring<>"\n\n";

If[i>=3,
Do[qstring = qstring<>"qk"<>ToString[j]<>ToString[k]<>"=reshape(q(:,"<>ToString[k]<>"),["<>StringRepeat["1,",j+1]<>"r]);";
If [k ==param1,qstring = qstring<>"\n"],{j,1,i-2},{k,1,param1}]];
(*If[i\[GreaterEqual]3,
Do[qstring = qstring<>"qk"<>ToString[j]<>ToString[k],{j,1,i-2},{k,1,param1}]];*)


(*qstring = "";Do[qstring = qstring<>"qk"<>ToString[j]<>ToString[k]<>"=qk("<>ToString[j]<>","<>ToString[k]<>");",{j,1,i-1},{k,1,param1}];
qstring = qstring<>"\n\n";
qstring = qstring<>"q=reshape(q,[param1,r])'; \nq=repmat(q,[r,1]);";
qstring = qstring<>"\n\n";
Do[qstring = qstring<>"ql"<>ToString[k]<>" = reshape(q(:,"<>ToString[k]<>"),[r,r]);\n",{k,1,param1}];
qstring = qstring<>"\n";
Do[qstring = qstring<>"qj"<>ToString[k]<>" = ql"<>ToString[k]<>"';",{k,1,param1}];
*)
qinputs ="qk,q";

varreplaceF = Table[qgen[t][[k]]->"ql"<>ToString[k],{k,1,param1}];
varreplaceF = Join[varreplaceF,Table[qgen[t][[k]]->Symbol["qj"<>ToString[k-param1]],{k,param1+1,2*param1}]];
varreplaceF = Join[varreplaceF,Table[xc[t][[k]]->Symbol["xc"<>ToString[k]],{k,1,dim}]];
qrep = Table[qgenk[t][[k,j]]->Symbol["qk"<>ToString[k]<>ToString[j]],{k,1,i-1},{j,1,param1}];
qrep = Flatten[qrep];
varreplaceF = Join[varreplaceF,qrep];

)
];

centerstring = "";Do[centerstring = StringJoin[centerstring,"xc"<>ToString[k]<>" = "<>ToMatlab[-xtildevalsF[[k]]/.varreplaceF]<>"\n"],{k,1,dim}];

Do[intFTableprint[[i,k]]="sum("<>StringReplace[ToMatlab[intFTable[[i,k]]/.varreplaceF],";"~~___->""]<>",'all');\n",{k,1,param1}];

filename = OpenWrite["F"<>ToString[i-1]<>"block_"<>pdename<>".m"];WriteMatlabFEXP[Transpose[{intFTableprint[[i,;;]]}]/.varreplaceF,ToString[i-1],qinputs,qstring,StringRiffle[params,","],prodk/.varreplaceF,centerstring,acsk/.varreplaceF,pdename,filename]; Close[filename],{i,1,ordernonlin }]
(*filename = OpenWrite["F"<>ToString[i-1]<>"block_"<>pdename<>".m"];WriteMatlabFEXP[Transpose[{intFTable[[i,;;]]}]/.varreplaceF,ToString[i-1],qinputs,qstring,StringRiffle[params,","],prodk/.varreplaceF,centerstring,acsk/.varreplaceF,pdename,filename]; Close[filename],{i,1,ordernonlin - 1+1}]*)





frhsstring = "function dqdt = frhs_"<>pdename<>"(t,q,"<>"r,param1,delta,"<>StringRiffle[params,","]<>")";

frhsstring = frhsstring<>"\n\nparams = length(q);\nM = zeros(params);\nf = zeros(params, 1);\n\n";
(*ordernonlin - 1=4;*)
frhsstring = frhsstring<>"%build M using blocks
for i = 1:r
    for j = 1:i
        Mij = Mblock_"<>pdename<>"(...
            [q(param1*(i-1)+1:param1*i),q(param1*(j-1)+1:param1*j)]);
        M(param1*(i-1)+1:param1*i,param1*(j-1)+1:param1*j) = Mij;
        M(param1*(j-1)+1:param1*(j),param1*(i-1)+1:param1*(i)) = Mij';
    end
end

qk = reshape(q,[param1,r])';

";
iord = 0;
frhsstring = frhsstring<>StringRepeat["\t",iord]<>"for k = 1:r\n\tf((1:param1) + param1*(k-1)) = ...";
Do[frhsstring = frhsstring<>"\n"<>
StringRepeat["\t",2]<>"F"<>ToString[iord]<>"block_"<>pdename<>"(qk(k,:),q,r,param1,"<>StringRiffle[params,","]<>")+...";
,{iord,0,ordernonlin - 1}];
frhsstring = StringDrop[frhsstring,-4]<>";";
frhsstring = frhsstring<>"\nend";

frhsstring = frhsstring<>"\n\ndqdt = (M + delta * eye(params)) \\ f;
 
end";

filename = OpenWrite["frhs_"<>pdename<>".m"];
WriteMatlabFRHS[frhsstring,filename]; 
Close[filename];


filename = OpenWrite["frhs_"<>pdename<>".m"];
WriteMatlabFRHS[frhsstring,filename]; 
Close[filename];


