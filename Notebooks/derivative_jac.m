#!/usr/local/bin/WolframScript -script
upfunc = "d_up";
downfunc = "d_down";
cc = 65;
nmax=7;

derivative[st_, i_] := 
  Module[{splitpos, splitfield, split, q0, stringup, stringdown, der, 
    res, macroname, totres, index},
   macroname = 
    If[i <= st[[3]], 
     StringJoin["d", FromCharacterCode[i + cc], "_", st[[1]]], 
     StringJoin["d", FromCharacterCode[i + cc], "_", st[[1]], ",", 
      FromCharacterCode[i + cc]]];
   index = If[i <= st[[3]], st[[3]], st[[3]] + 1];
   split = StringSplit[st[[2]], "*+*"];
   der = Table["", {q0, 1, Length[split]}];
   For[q0 = 1, q0 <= Length[split], q0++,
    stringup = StringJoin[split[[q0]], "_u:", ToString[i]];
    stringdown = StringJoin[split[[q0]], "_d:", ToString[i]];
    der[[q0]] = {"*+*", stringup, "*+*", stringdown}
    
    ];
   res = StringJoin@der;
   totres = {macroname, res, index};
   totres
   ];

simplify[st_] := 
  Module[{q0, q1, split, up, down, der, res, n, toreplace1, 
    toreplace2, pos1, pos2, totres},
   split = StringSplit[st[[2]], "*+*"];
   der = Table["", {q0, 1, Length[split]}];
   For[q0 = 1, q0 <= Length[split], q0++,
    For[q1 = 1, q1 <= Length[split], q1++,
     up = StringCases[split[[q0]], StringJoin["u:", ToString[q1]]];
     down = StringCases[split[[q0]], StringJoin["d:", ToString[q1]]];
     n = If[Length[up] >= Length[down], Length[down], Length[up]];
     toreplace1 = StringJoin["d:", ToString[q1]];
     pos1 = StringPosition[split[[q0]], toreplace1, n];
     toreplace2 = StringJoin["u:", ToString[q1]];
     pos2 = StringPosition[split[[q0]], toreplace2, n];
     split[[q0]] = 
      If[n > 0, StringReplacePart[split[[q0]], "n:-", pos1], 
       split[[q0]]];
     split[[q0]] = 
      If[n > 0, StringReplacePart[split[[q0]], "n:+", pos2], 
       split[[q0]]];
     ];
    der[[q0]] = {"*+*", split[[q0]]};
    ];
   res = StringJoin@der;
   totres = {st[[1]], res};
   totres
   ];


prova={"phi(A", "A", 0};

recursivesum[st_,index_]:=Module[{q0,temp},
				 temp=st;
				 For[q0 = 1, q0 <= Length[index], q0++,
				     temp = sumEPE[temp, index[[q0]]];
				     ];
				 temp
				 ];




recursive[st_, index_] := Module[{q0, temp},
				 temp = st;
				 For[q0 = 1, q0 <= Length[index], q0++,
				     temp = derivative[temp, index[[q0]]];
				     ];
				 toEPE[temp]
				 ];

generateindex[imax_] := Module[{a, q0, q1, qi, part, nindex, i, l},
			       l = {};
			       For[part = 1, part <= imax, part++,
				   For[nindex = 1, nindex <= part, nindex++,
				       a = IntegerPartitions[part, {nindex}];
				       For[i = 1, i <= Length[a], i++,
					   index = Table[0, {i, 1, part}];
					   count = 1;
					   For[q0 = 1, q0 <= Length[a[[i]]], q0++,
					       For[q1 = 1, q1 <= a[[i, q0]], q1++,
						   index[[count]] = q0;
						   count++;
						   ];
					       ];
					   AppendTo[l, index]
					   ];
				       ];
				   ];
			       l
			       ];




doubleindex[index_]:=Module[{count,listcount,di},
			    listcount={};
			    For[q0 = 1, q0 <= Max[index], q0++,
				count=Count[index,q0];
				AppendTo[listcount,count];
				];
			    di={};
			    For[q0 = 1, q0 <= Max[index], q0++,
				If[listcount[[q0]]==2,AppendTo[di,q0]];
				];
			    di
			    ];



generatemacro[imax_,ifieldname_] := Module[{macro,macroC,count,name,dupl,sumind, index, q0,q1,q2, prova,subsets,temp},
					     macro = {};
					     index = generateindex[imax];
					     name=StringJoin[ToString[ifieldname],"(tflow,A"];
					     prova = {name, "A", 0};
					     For[q0 = 1, q0 <= Length[index], q0++,
						 temp=recursive[prova, index[[q0]]];
						 AppendTo[macro,temp];
						 subsets=DeleteDuplicates[Subsets[doubleindex[index[[q0]]]]];

						 For[q1 = 2, q1 <= Length[subsets], q1++,
						    AppendTo[macro,recursivesum[temp,subsets[[q1]]]];
						     ];
						 ];
					     macroC={};
					     For[q0 = 1, q0 <= Length[macro], q0++,
						 AppendTo[macroC,EPEtoC[macro[[q0]]]];
						 ];
					   macroC
					     ];

 


outermatch[string_]:=Module[{open,closed,q0,match,ordop,ordcl,fmatch,fmatchpos,matchop,matchcl},
			    open=StringPosition[string,"["];
			    closed=StringPosition[string,"]"];
			    match={};
			    ordop=Table[0,{i,1,Length[open]}];
			    ordcl=Table[0,{i,1,Length[open]}];

			    For[q0=1,q0<=Length[open],q0++,
				ordop[[q0]]=open[[q0]][[1]];
				ordcl[[q0]]=closed[[q0]][[1]];
				];
			    
			    count=0;
			    sist=0;
			     For[q0=1,q0<=Length[ordop],q0++,
				fmatch=Select[ordop, # > ordcl[[q0]]  &, 1];
				 fmatchpos=If[Length[fmatch]==1,Position[ordop,fmatch[[1]]][[1,1]],Length[ordop]+1];
				 If[fmatchpos-q0==1,AppendTo[match,{ordop[[fmatchpos-q0+sist]],ordcl[[q0]]}],fmatch=0];
				 If[fmatchpos-q0==1,sist=q0,fmatch=0];
				];
			    match
			    ]



dim=3;





toEPE[st_] :=
  Module[{simp, split, subsplit, q0, q1, der, res, arg, dir, coeff,
	  sign, duplicates, dersimp, field},
    simp = simplify[st];
    split = StringSplit[simp[[2]], "*+*"];
    der = Table["", {j, 1, Length[split]}];
    subsplit = StringSplit[st[[1]], "_"];
    field = StringCases[subsplit[[Length[subsplit]]], __ ~~ "("];
    field = StringTrim[StringReplace[field, "(" -> " "]];
    For[q0 = 1, q0 <= Length[split], q0++,
	subsplit = StringSplit[split[[q0]], "_"];
	arg = subsplit[[1]];
	coeff = 1;
	For[q1 = 2, q1 <= Length[subsplit], q1++,
	    dir = StringSplit[subsplit[[q1]], ":"];
	    coeff *= If[dir[[1]] == "d", -1, If[dir[[2]] == "-", -1, 1]];
     arg = If[dir[[2]]!="S",If[dir[[1]] == "n", arg,
	       If[dir[[1]] == "d",StringJoin[downfunc, "[", arg, "*d+",
			     FromCharacterCode[ToExpression[dir[[2]]] + cc], "]"],
		  StringJoin[upfunc, "[", arg, "*d+",
			     FromCharacterCode[ToExpression[dir[[2]]] + cc], "]"]]],If[dir[[1]] == "n", arg,
	       If[dir[[1]] == "d",StringJoin[downfunc, "[", arg, "*d+S]"],
		  StringJoin[upfunc, "[", arg, "*d+S]"]]]];
	  

	    ];
	sign = If[coeff < 0, "-", "+"];
	der[[q0]] = StringJoin[sign, "/[", arg, "]"];
	];

    dersimp = Table["", {j, 1, Length[der]}];

    For[q0 = 1, q0 <= Length[der], q0++,
	duplicates = Count[der, der[[q0]]];
	subsplit = StringSplit[der[[q0]], "/"];
    dersimp[[q0]] =StringJoin["*$*",subsplit[[1]], ToString[duplicates], "$$", field,"$$",
			      subsplit[[2]]];

	];

    dersimp = DeleteDuplicates[dersimp];
    res = StringJoin@dersimp;

    StringJoin["#define ", st[[1]], " ",
	       ToString[N[1/Length[split]]], " ", res]
    ];





sumEPE[st_,index_]:=Module[{splitblank,split,subsplit,q0,q1,sub,coeff,newcoeff,arg,summed,match,res,totres,resstring},
			   splitblank=StringSplit[st];
			   splitblank[[2]]=StringDelete[splitblank[[2]],StringJoin[",",FromCharacterCode[ToExpression[index]+cc]]];
			   splitblank[[2]]=StringJoin["s",FromCharacterCode[ToExpression[index]+cc],"_",splitblank[[2]]];
			   match=StringJoin["*",FromCharacterCode[ToExpression[index]+cc],"*"];
			   split=StringSplit[splitblank[[4]],"*$*"];
			   summed={};
			   For[q0=1,q0<=Length[split],q0++,
			       subsplit=StringSplit[split[[q0]],"$$"];
			       coeff=If[StringMatchQ[subsplit[[3]],match],subsplit[[1]],ToString[+3*ToExpression[subsplit[[1]]]]];
			       newcoeff=If[ToExpression[coeff]>1,StringJoin["+",ToString[ToExpression[coeff]]],coeff];
			       If[StringMatchQ[subsplit[[3]],match],
				  For[q1=0,q1<3,q1++,
				      AppendTo[summed,StringJoin["*$*",newcoeff,"$$",subsplit[[2]],"$$",StringReplace[subsplit[[3]],FromCharacterCode[ToExpression[index]+cc]->ToString[q1]]]];


				      ],AppendTo[summed,StringJoin["*$*",newcoeff,"$$",subsplit[[2]],"$$",subsplit[[3]]]]
				  ];		       
																    
			       ];
			   summed=StringDelete[summed,"+0"];
			   res=StringJoin@summed;
			   totres={splitblank[[1]],StringJoin[" ",splitblank[[2]]],StringJoin[" ",splitblank[[3]],StringJoin[" ",res]]};
			   resstring=StringJoin@totres;
			   resstring
]

  EPEtoC[st_]:=Module[{splitblank,split,coeff,subsplit,cvers,q0},
		      splitblank=StringSplit[st];
		      splitblank[[1]]=StringJoin[splitblank[[1]]," "];
		      splitblank[[2]]=StringJoin[splitblank[[2]],") "];
		      splitblank[[3]]=StringJoin["(",splitblank[[3]],"*("];
		      split=StringSplit[splitblank[[4]],"*$*"];
		      cvers={};
		      For[q0=1,q0<=Length[split],q0++,
			  subsplit=StringSplit[split[[q0]],"$$"];
			  coeff=Which[subsplit[[1]]=="+1","+",subsplit[[1]]=="-1","-",True,StringJoin[subsplit[[1]],"*"]];


			  AppendTo[cvers,StringJoin[coeff,subsplit[[2]],"[tflow]",subsplit[[3]]]];

			  ];





			  splitblank[[4]]=StringJoin@cvers;
			  splitblank[[4]]=StringJoin[cvers,"))"];


		      res=StringJoin@splitblank;
		      res

]

Print["#ifndef __MACRO_JAC_D"];
Print["#define __MACRO_JAC_D"];
For[qi=3,qi<=Length[$ScriptCommandLine],qi++,

    Print[StringJoin["#define ",$ScriptCommandLine[[qi]],"(tflow,A) ",$ScriptCommandLine[[qi]],"[tflow][A]"]];

    macrolist=If[ToExpression[$ScriptCommandLine[[2]]]<=nmax,generatemacro[ToExpression[$ScriptCommandLine[[2]]],$ScriptCommandLine[[qi]]],{StringJoin["N too large it is set to be less than ",ToString[nmax]]}];

    For[q0=1,q0<=Length[macrolist],q0++,
	Print[macrolist[[q0]]];
	];
    macrolist={};
    ];
Print["#endif"];
