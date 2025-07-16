(* ::Package:: *)

(* ::Input::Initialization:: *)
Get["ToMatlab_SRONS.m"];
Print["Writing files to Matlab"]
dir = NotebookDirectory[]<>pdename;
If[DirectoryQ[dir],SetDirectory[dir],CreateDirectory[dir];
SetDirectory[dir]];
DumpSave["GeneralInts_Symbolic_"<>pdename<>".mx","Global`"]

(** Write M **)
varreplace = Table[qgen[t][[k]]->Symbol["q"<>ToString[k]],{k,1,param2}];
qstring = "";
Do[qstring = qstring<>"q"<>ToString[k]<>"=q("<>ToString[k]<>");",{k,1,param2}]

filename = OpenWrite["Mblock_"<>pdename<>".m"];
WriteMatlabM[Mblock/.varreplace,qstring,pdename,filename]; 
Close[filename];

(** Write f **)

Do[If[i==1,
(
qstring = "";Do[qstring = qstring<>"ql"<>ToString[k]<>"=ql("<>ToString[k]<>");",{k,1,param1}];
qstring = qstring<>"\n \nqj = reshape(q,[param1,r])'; \n \n";
Do[qstring = qstring<>"qj"<>ToString[k]<>"=qj(:,"<>ToString[k]<>");",{k,1,param1}];

qinputs = "ql,q";

varreplaceF = Table[qgen[t][[k]]->"ql"<>ToString[k],{k,1,param1}];
varreplaceF = Join[varreplaceF,Table[qgen[t][[k]]->Symbol["qj"<>ToString[k-param1]],{k,param1+1,2*param1}]];
),
(
qstring = "";Do[qstring = qstring<>"qk"<>ToString[j]<>ToString[k]<>"=qk("<>ToString[j]<>","<>ToString[k]<>");",{j,1,i-1},{k,1,param1}];
qstring = qstring<>"\n\n";
qstring = qstring<>"q=reshape(q,[param1,r])'; \nq=repmat(q,[r,1]);";
qstring = qstring<>"\n\n";
Do[qstring = qstring<>"ql"<>ToString[k]<>" = reshape(q(:,"<>ToString[k]<>"),[r,r]);\n",{k,1,param1}];
qstring = qstring<>"\n";
Do[qstring = qstring<>"qj"<>ToString[k]<>" = ql"<>ToString[k]<>"';",{k,1,param1}];

qinputs ="qk,q";

varreplaceF = Table[qgen[t][[k]]->"ql"<>ToString[k],{k,1,param1}];
varreplaceF = Join[varreplaceF,Table[qgen[t][[k]]->Symbol["qj"<>ToString[k-param1]],{k,param1+1,2*param1}]];
qrep = Table[qgenk[t][[k,j]]->Symbol["qk"<>ToString[k]<>ToString[j]],{k,1,i-1},{j,1,param1}];
qrep = Flatten[qrep];
varreplaceF = Join[varreplaceF,qrep];

)
];

filename = OpenWrite["F"<>ToString[i-1]<>"block_"<>pdename<>".m"];WriteMatlabF[Transpose[{intFTable[[i,;;]]}]/.varreplaceF,ToString[i-1],qinputs,qstring,StringRiffle[params,","],pdename,filename]; Close[filename],{i,1,ordernonlin}]





frhsstring = "function dqdt = frhs_"<>pdename<>"(t,q,"<>"r,param1,delta,"<>StringRiffle[params,","]<>")";

frhsstring = frhsstring<>"\n\nparams = length(q);\nM = zeros(params);\nf = zeros(params, 1);\n\n";
(*ordernonlin=4;*)
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


