(* ::Package:: *)

(*SRONS*)

BeginPackage["SRONSUtils`SRONSintegrate`"]

SRONS::usage 

SRONSEXP::usage 

(*SetMargin::usage = "SetMargin[margin]"
RestoreMargin::usage = "RestoreMargin[]"*)


Begin["`Private`"]
(*SetDirectory[NotebookDirectory[]];*)

SRONS[dim_, phi_,posamp_] :=
    (Get["General_Symbolic_M.m"]; 
    Get["General_Symbolic_f.m"];
    Get["M_f_loadints.m"];
    Get["Write_Files.m"];)

SRONSEXP[dim_, phi_,posamp_] :=
    (Get["General_Symbolic_EXP_M.m"]; 
    Get["General_Symbolic_EXP_f_general.m"];
    Get["M_f_loadints_EXP.m"];
    Get["Write_Files_EXP.m"];)



End[];
EndPackage[];
