#!/usr/local/bin/WolframScript -script

GenerateOp[nd_]:=Module[{list={},q0},
			For[q0=0,q0<=nd,q0++,
			    
			    AppendTo[list,Join[Table[dd[lor],{i,1,q0}],{Phi}]];
			    ];
			list

			];


Generate[nd_,nf_]:=Module[{list={},finallist={}},
			  list=Tuples[GenerateOp[nd],nf];
			  For[q0=1,q0<=Length[list],q0++,
			      count=0;
			      For[q1=1,q1<=Length[list[[q0]]],q1++,
				  count+=Count[list[[q0]][[q1]],dd[lor]];
				  ];
			      If[count==nd,AppendTo[finallist,list[[q0]]],count=0];
			      ];

			  For[q0=1,q0<=Length[finallist],q0++,
			      finallist[[q0]]=Sort[finallist[[q0]]];
			      ];
			  finallist=DeleteDuplicates[finallist];
			  
			  finallist
			  ];



ExtLorentz[el_,mu_]:=Module[{res={},idx,newidx,nidx,all,count=0},
			    For[q0=1,q0<=Length[el],q0++,
				count+=Count[el[[q0]],dd[lor]];
				];

			    idx=Table[lor,{i,1,count}];
			   
If[EvenQ[Length[idx]],Print["Error: wrong number of lorentz indices"];];
Do[
   all=Union[Permutations[Flatten[Table[{a[j],a[j]},{j,1,(Length[idx]-1)/2}]]]];
   all=Union[all/.{
       {A1___,a[e1_],A2___,a[e1_],A3___,a[e2_],A4___,a[e2_],A5___}:> {A1,a[e2],A2,a[e2],A3,a[e1],A4,a[e1],A5}/;e2<e1, 
       {A1___,a[e1_],A2___,a[e2_],A3___,a[e1_],A4___,a[e2_],A5___}:> {A1,a[e2],A2,a[e1],A3,a[e2],A4,a[e1],A5}/;e2<e1,{A1___,a[e1_],A2___,a[e2_],A3___,a[e2_],A4___,a[e1_],A5___}:> {A1,a[e2],A2,a[e1],A3,a[e1],A4,a[e2],A5}/;e2<e1}]; 
   all=Union[all/.{ 
       {A1___,a[e1_],A2___,a[e1_],A3___,a[e2_],A4___,a[e2_],A5___}:> {A1,a[e2],A2,a[e2],A3,a[e1],A4,a[e1],A5}/;e2<e1, 
       {A1___,a[e1_],A2___,a[e2_],A3___,a[e1_],A4___,a[e2_],A5___}:> {A1,a[e2],A2,a[e1],A3,a[e2],A4,a[e1],A5}/;e2<e1,{A1___,a[e1_],A2___,a[e2_],A3___,a[e2_],A4___,a[e1_],A5___}:> {A1,a[e2],A2,a[e1],A3,a[e1],A4,a[e2],A5}/;e2<e1}]; 
   all=Union[all/.{ 
       {A1___,a[e1_],A2___,a[e1_],A3___,a[e2_],A4___,a[e2_],A5___}:> {A1,a[e2],A2,a[e2],A3,a[e1],A4,a[e1],A5}/;e2<e1, 
       {A1___,a[e1_],A2___,a[e2_],A3___,a[e1_],A4___,a[e2_],A5___}:> {A1,a[e2],A2,a[e1],A3,a[e2],A4,a[e1],A5}/;e2<e1,{A1___,a[e1_],A2___,a[e2_],A3___,a[e2_],A4___,a[e1_],A5___}:> {A1,a[e2],A2,a[e1],A3,a[e1],A4,a[e2],A5}/;e2<e1}]; 
   all=Union[all//. {{A1___, a[e1_], A2___, a[e1_], A3___, a[e2_], A4___, a[e2_], A5___} :> {A1, a[e2], A2, a[e2], A3, a[e1], A4, a[e1], A5} /;  e2 < e1,
	     {A1___, a[e1_], A2___, a[e2_], A3___, a[e1_], A4___, a[e2_], A5___} :> {A1, a[e2], A2, a[e1], A3, a[e2], A4, a[e1], A5} /; e2 < e1, {A1___, a[e1_], A2___, a[e2_], A3___, a[e2_], A4___, a[e1_], A5___} :> {A1, a[e2], A2, a[e1], A3, a[e1], A4, a[e2], A5} /; e2 < e1}];
			    Do[ 
			       newidx=Insert[all[[k]],mu,i]; 
			       nidx=1; 
			       new=el; 
			       
			       Do[ 
				  Do[
				     If[Not[FreeQ[el[[j,q]],dd]],new[[j,q]]=el[[j,q]]//. dd[a1_]-> dd[newidx[[nidx]]];nidx=nidx+1];
					, {q, 1, Length[el[[j]]]}];
				     ,{j,1,Length[el]}]; 
				  res = Union[res, {new}];
				  ,{k,1,Length[all]}];
			       ,{i,1,Length[idx]}];
			    Do[
			       Do[
				  res[[i,k]]=Sort[res[[i,k]]];

				  If[Length[res[[i,k]]]>1,res[[i,k]]=Permute[res[[i,k]],Cycles[{{1,Length[res[[i,k]]]}}]]];

				  ,{k,1,Length[res[[i]]]}];
			       res[[i]]=Sort[res[[i]]];
			       ,{i,1,Length[res]}];


			    res=DeleteDuplicates[res];

			    res];



VectorialProbes[ndim_]:=Module[{q0,q1,q2,list={},probes}, 
			       For[q0=2,q0<=2 ndim,q0+=2,
				   For[q1=1,q1+q0/2<=ndim,q1+=2,
				       probes=Generate[q1,q0];
				       For[q2=1,q2<=Length[probes],q2++,
					   
					   list=Join[list,ExtLorentz[probes[[q2]],rho]];
					   
					   ];
				       ];
				   ];

			       list


			       ];

FunctionalD[list_,index_]:=Module[{listder={},q0,obj,listtemp={},totlist={},count,countlist={}},
				 
			   For[q0=1,q0<=Length[list],q0++,
			       obj=list[[q0]];
			       obj[[Length[obj]]]=Phiboundary;
			       AppendTo[obj,dd[index]];
			       listtemp=list;
			       obj=Sort[obj];
			       If[Length[obj]>1,obj=Permute[obj,Cycles[{{1,Length[obj]}}]]];
			       listtemp[[q0]]=obj;
			       listtemp=Sort[listtemp];
			       AppendTo[totlist,listtemp];
			       ];

				  For[q0=1,q0<=Length[totlist],q0++,
				      count=Count[totlist,totlist[[q0]]];
				      AppendTo[countlist,count];
				      ];

				  For[q0=1,q0<=Length[totlist],q0++,
				      PrependTo[totlist[[q0]],{countlist[[q0]]}];

				      ];

				  totlist=DeleteDuplicates[totlist];
				  totlist

			   ];

cc = 65;
toCprobe[list_,n_]:=Module[{q0,q1,q2,stringtot,string,idx,stringheader,idxrip={},countlist={},extidx={},opensq,closesq,opidx},
			For[q0=1,q0<=Length[list],q0++,
			   			    

			    string=" ";
			    string=ToString[list[[q0,Length[list[[q0]]]]]];
			    
			    idx=list[[q0]]//. {e1___,dd[a1_],e2___}-> {e1,a1,e2};
			    idx=Drop[idx,{Length[idx]}];
			    countlist={};
			   
			    For[q1=1,q1<=Length[idx],q1++,
				AppendTo[countlist,Count[idx,idx[[q1]]]];
				];
			    idxrip={};
			    
			    For[q1=1,q1<=Length[idx],q1++,
				AppendTo[idxrip,{idx[[q1]],countlist[[q1]]}];
				];
			    idxrip=DeleteDuplicates[idxrip];
			    idxrip=SortBy[idxrip,Last];
                            idxrip=Reverse[idxrip];

			    For[q1=1,q1<=Length[idxrip],q1++,
				For[q2=1,q2<=idxrip[[q1,2]],q2++,
				    string=StringJoin["d",FromCharacterCode[q1 + cc],"_",string];
				    ];
				];

			    For[q1=1,q1<=Length[idxrip],q1++,
				If[idxrip[[q1,2]]==2,string=StringJoin["s",FromCharacterCode[q1 + cc],"_",string]];
				];
			    
			    string=StringJoin[string,"(x"];


			    
			    For[q1=1,q1<=Length[idxrip],q1++,
				If[idxrip[[q1,2]]==1,string=StringJoin[string,",",ToString[idxrip[[q1,1]]]]];
				If[idxrip[[q1,2]]==1,AppendTo[extidx,idxrip[[q1,1]]]];
				];
			    string=StringJoin[string,")"];

			    stringtot=If[q0==1,string,StringJoin[string,"*",stringtot]];
			    ];

			

			countlist={};
			For[q0=1,q0<=Length[extidx],q0++,
			    AppendTo[countlist,Count[extidx,extidx[[q0]]]];

			    ];

			idxrip={};

			For[q0=1,q0<=Length[extidx],q0++,
			    AppendTo[idxrip,{extidx[[q0]],countlist[[q0]]}];

			    ];
			   idxrip=DeleteDuplicates[idxrip];


			   stringtot=StringJoin["sum+=",stringtot,";\n"];


			For[q0=1,q0<=Length[idxrip],q0++,
			    
			    If[idxrip[[q0,2]]==2,
			       stringtot=StringJoin["\n","for(",ToString[idxrip[[q0,1]]],"=0;",ToString[idxrip[[q0,1]]],"<d;",ToString[idxrip[[q0,1]]],"++){\n",stringtot];


			       ];
			    ];



			For[q0=1,q0<=Length[idxrip],q0++,
			    If[idxrip[[q0,2]]==2,
			       stringtot=StringJoin["int ",ToString[idxrip[[q0,1]]],";\n",stringtot];
			       ];
			    ];


			For[q0=1,q0<=Length[idxrip],q0++,
			    If[idxrip[[q0,2]]==2,
			       stringtot=StringJoin[stringtot,"\n}"];

			       ];
			    ];


			   For[q0=1,q0<=Length[idxrip],q0++,
			       If[idxrip[[q0,2]]==1,
				  opidx=ToString[idxrip[[q0,1]]];
				  ];
			       ];


			   factor=Factorial[Length[list]-1];
			   
			   stringtot=StringJoin["\ndouble sum=0.;\n",stringtot];
			   If[factor!=1, stringtot=StringJoin[stringtot,"\nreturn sum/",ToString[factor],".;\n"],stringtot=StringJoin[stringtot,"\nreturn sum;\n"]];
			   stringtot=StringJoin["double probe_",ToString[n],"(int x, int ",opidx,",double *Phi) {\n",stringtot];
			   stringheader=StringJoin["double probe_",ToString[n],"(int x, int ",opidx,",double *Phi);\n"];
			   stringtot=StringJoin[stringtot,"\n}\n\n"];
			   stringtot=StringReplace[stringtot,"["->""];
			   stringtot=StringReplace[stringtot,"]"->""];
			   {stringheader,stringtot}
			   ];

probinteststring="
#ifndef __TWIPROB
#define __TWIPROB

#include <math.h>
#include \"macro.h\"
#include <stdio.h>
#include \"observables.h\"
#include \"geometry.h\"
#include \"lattice_config.h\"
\n"


  cPrgen[listpr_,filename_,filein_]:=Module[{q0,filec,fileh,probenames={},filenamec,filenameh,cprobesc={},cprobesh={},temp},
				   
				   filenamec=StringJoin[filename,".c"];
				   filenameh=StringJoin[filename,".h"];
				   
				  				   
				   For[q0=1,q0<=Length[listpr],q0++,
				       temp=toCprobe[listpr[[q0]],q0];
				       AppendTo[cprobesc,temp[[2]]];
				       AppendTo[cprobesh,temp[[1]]];
				       AppendTo[probenames,StringTrim[StringReplace[StringCases[temp[[1]]," "~~__~~"("],"("->" "]]];
				       ];
				   
				   
				  filec=OpenWrite[filenamec];
				    WriteString[filec,"#include \"",filenameh,"\"\n"];
				    
				    For[q0=1,q0<=Length[cprobesc],q0++,
				      WriteString[filec,cprobesc[[q0]]];
				      ];

				  fileh=OpenWrite[filenameh];
				    WriteString[fileh,probinteststring];
				    WriteString[fileh,"\n#define NUM_PROB ",ToString[Length[listpr]],"\n"];
				    WriteString[fileh,"\ntypedef double (*prb_ptr_t)(int,int,double *);\n"];
				    
				    
				    For[q0=1,q0<=Length[cprobesh],q0++,
					WriteString[fileh,cprobesh[[q0]]];
					];
				    
				    WriteString[filein,StringJoin["\nstatic prb_ptr_t PROB_OP[",ToString[Length[listpr]],"]={"]];
				    WriteString[filein,probenames[[1,1]]];
				    For[q0=2,q0<=Length[probenames],q0++,
					WriteString[filein,StringJoin[",",probenames[[q0,1]]]];
					];
				    WriteString[filein,"};\n\n"];
				    WriteString[fileh,"\n#endif"];
				    
				    Close[filec];
				    Close[fileh];
				    
				    ];





toCprobered[list_,n_]:=Module[{q0,q1,q2,stringtot,string,idx,idxrip={},countlist={},extidx={},opensq,closesq,opidx},
			   For[q0=1,q0<=Length[list],q0++,
			       string=" ";
			       string=ToString[list[[q0,Length[list[[q0]]]]]];
			       
			       idx=list[[q0]]//. {e1___,dd[a1_],e2___}-> {e1,a1,e2};
			       idx=Drop[idx,{Length[idx]}];
			       countlist={};
			       
			       For[q1=1,q1<=Length[idx],q1++,
				   AppendTo[countlist,Count[idx,idx[[q1]]]];
				   ];
			       idxrip={};
			       
			       For[q1=1,q1<=Length[idx],q1++,
				   AppendTo[idxrip,{idx[[q1]],countlist[[q1]]}];
				   ];
			       idxrip=DeleteDuplicates[idxrip];
			       idxrip=SortBy[idxrip,Last];
			       idxrip=Reverse[idxrip];
			       
			       For[q1=1,q1<=Length[idxrip],q1++,
				   For[q2=1,q2<=idxrip[[q1,2]],q2++,
				       string=StringJoin["d",FromCharacterCode[q1 + cc],"_",string];
				       ];
				   ];
			       
			       For[q1=1,q1<=Length[idxrip],q1++,
				   If[idxrip[[q1,2]]==2,string=StringJoin["s",FromCharacterCode[q1 + cc],"_",string]];
				   ];
			       
			       string=If[ToString[list[[q0,Length[list[[q0]]]]]]=="Phiboundary",StringJoin[string,"(x"],StringJoin[string,"(y"]];
			       
			       
			       
			       For[q1=1,q1<=Length[idxrip],q1++,
				   If[idxrip[[q1,2]]==1,string=StringJoin[string,",",ToString[idxrip[[q1,1]]]]];
				   If[idxrip[[q1,2]]==1,AppendTo[extidx,idxrip[[q1,1]]]];
				   ];
			       string=StringJoin[string,")"];
			       
			       stringtot=If[q0==1,string,StringJoin[string,"*",stringtot]];
			       ];
			   
			   
			   
			   countlist={};
			   For[q0=1,q0<=Length[extidx],q0++,
			       AppendTo[countlist,Count[extidx,extidx[[q0]]]];
			       
			       ];
			   
			   idxrip={};
			   
			   For[q0=1,q0<=Length[extidx],q0++,
			       AppendTo[idxrip,{extidx[[q0]],countlist[[q0]]}];
			       
			       ];
			      idxrip=DeleteDuplicates[idxrip];
			      idxrip=SortBy[idxrip,Last];
			      idxrip=Reverse[idxrip];			   
			   
			      stringtot=StringJoin["sum",ToString[n],"+=",stringtot,";\n"];
			   
			      For[q0=1,q0<=Length[idxrip],q0++,
				  If[idxrip[[q0,2]]==2,

				     stringtot=StringReplace[stringtot,ToString[idxrip[[q0,1]]]->StringJoin[ToString[idxrip[[q0,1]]],ToString[n]]];

				     ];
				  ];
			   For[q0=1,q0<=Length[idxrip],q0++,
			       
			       If[idxrip[[q0,2]]==2,
				  stringtot=StringJoin["\n","for(",StringJoin[ToString[idxrip[[q0,1]]],ToString[n]],"=0;",StringJoin[ToString[idxrip[[q0,1]]],ToString[n]],"<d;",StringJoin[ToString[idxrip[[q0,1]]],ToString[n]],"++){\n",stringtot];
				  
				  
				  ];
			       ];
			   
			   
			   
			      For[q0=1,q0<=Length[idxrip],q0++,
				  If[idxrip[[q0,2]]==2,
				     stringtot=StringJoin["int ",StringJoin[ToString[idxrip[[q0,1]]],ToString[n]],";\n",stringtot];
				     ];
				  ];
			      
			      
			      For[q0=1,q0<=Length[idxrip],q0++,
				  If[idxrip[[q0,2]]==2,
				     stringtot=StringJoin[stringtot,"\n}"];
				     
				     ];
				  ];
			      

			      stringtot=StringJoin["\ndouble sum",ToString[n],"=0.;\n",stringtot];			   

			      factor=Factorial[Length[list]-1];			      
			      If[factor!=1, stringtot=StringJoin[stringtot,"\nsum",ToString[n],"/=",ToString[factor],".;\n"]];

			      stringtot=StringReplace[stringtot,"["->""];
			      stringtot=StringReplace[stringtot,"]"->""];
			      {stringtot}
			      ];



jacobianfuncheader="double TWIJac(int x,int y,int tflow);\n"


jacobianfunc="\ndouble TWIJac(int x,int y,int tflow){ \n   
  int jp;\n
  int xJ[d];\n
  int yJ[d];\n
  coordinates(x,xJ);\n
  coordinates(y,yJ);\n
  int xyJ[d];\n
  for(jp=0;jp<d;jp++) xyJ[jp]=abs(xJ[jp]-yJ[jp]);\n
  int iJ=get_index(xyJ[0],xyJ[1],xyJ[2]);\n
 return jac_v[iJ+Vol*tflow];\n\n}\n\n";

jacobianfuncname="TWIJac(x,y,tflow)";

jacobian="jac_v[tflow]"



toCvariation[list_,n_]:=Module[{q0,q1,idxfree={},temp,temp2,countlist={},idxrip={},idx={},fact,stringlist={},element,tempstring,stringtot,stringheader=""},
					  For[q0=1,q0<=Length[list],q0++,

					      element=list[[q0]];
					      
					      fact=element[[1,1]];
					      element=Drop[element,1];
					      temp=element//. {e1___,dd[a1_],e2___}-> {e1,a1,e2};
					      idx={};
					      For[q1=1,q1<=Length[temp],q1++,
						  If[Length[temp[[q1]]]>1,
						     temp2=Drop[temp[[q1]],{Length[temp[[q1]]]}];
						     AppendTo[idx,temp2];
						     ];
						  ];
					      idx=Flatten[idx];
					      tempstring=toCprobered[element,q0];
					      If[fact!=1,tempstring=StringJoin[tempstring,"\nsum",ToString[q0],"*=",ToString[fact],";\n"];];
					      AppendTo[stringlist,tempstring];
					      
					      ];
			    

			    For[q0=1,q0<=Length[idx],q0++,
				AppendTo[countlist,Count[idx,idx[[q0]]]];

				];

			    For[q0=1,q0<=Length[idx],q0++,
				AppendTo[idxrip,{idx[[q0]],countlist[[q0]]}];

				];

			       idxrip=DeleteDuplicates[idxrip];
			       
			       idxrip=SortBy[idxrip,Last];
			       idxrip=Reverse[idxrip];
			    For[q0=1,q0<=Length[idxrip],q0++,
				If[idxrip[[q0,2]]==1,
				   AppendTo[idxfree,idxrip[[q0,1]]];
				   ];
				];
			       
			       
			       idxfree=Sort[idxfree];
			       stringtot=StringJoin@stringlist;
			       

			    stringtot=StringJoin[stringtot,"\nreturn "];
			    stringtot=StringJoin[stringtot,"(sum1"];
			    For[q0=2,q0<=Length[list],q0++,
				stringtot=StringJoin[stringtot,"+sum",ToString[q0]];

				];
			       stringtot=StringJoin[stringtot,")*",jacobian,";\n}\n\n"];
			       
			       stringtot=StringJoin["){\n",stringtot];
			       stringheader=StringJoin[");\n",stringheader];
			       stringtot=StringJoin[",int tflow,double *Phi,double *Phiboundary",stringtot];
			       stringheader=StringJoin[",int tflow,double *Phi,double *Phiboundary",stringheader];
			       
			       For[q0=1,q0<=Length[idxfree],q0++,
				   stringtot=StringJoin[",int ",ToString[idxfree[[q0]]],stringtot];
				   stringheader=StringJoin[",int ",ToString[idxfree[[q0]]],stringheader];
				   ];

			    
			       stringtot=StringJoin["double delta_pr_",ToString[n],"(int x,int y",stringtot];
			       stringheader=StringJoin["double delta_pr_",ToString[n],"(int x,int y",stringheader];

			       {stringheader,stringtot}
			    ];


varinteststring="
#ifndef __TWIVAR
#define __TWIVAR

#include <math.h>
#include \"macro.h\"
#include <stdio.h>
#include \"observables.h\"
#include \"geometry.h\"
#include \"lattice_config.h\"
\n"

  cVargen[listpr_,index_,filename_,filein_]:=Module[{q0,filec,fileh,filenamec,filenameh,cprobesc={},tempfunc,cprobesh={},temp,probenames={}},
							    
							    filenamec=StringJoin[filename,".c"];
							    filenameh=StringJoin[filename,".h"];
							    
							    
							    For[q0=1,q0<=Length[listpr],q0++,
								temp=toCvariation[FunctionalD[listpr[[q0]],index],q0];
								
								AppendTo[cprobesc,temp[[2]]];
								AppendTo[cprobesh,temp[[1]]];
								AppendTo[probenames,StringTrim[StringReplace[StringCases[temp[[1]]," "~~__~~"("],"("->" "]]];
								];
							    
							    filec=OpenWrite[filenamec];
							    
							    WriteString[filec,StringJoin["#include \"",filenameh,"\"\n"]];
							    
						    (*WriteString[filec,jacobianfunc];*)
							    For[q0=1,q0<=Length[cprobesc],q0++,                                                                                                                                               
								WriteString[filec,cprobesc[[q0]]];
								];
							    
							    fileh=OpenWrite[filenameh];
							    WriteString[fileh,varinteststring];
							    (*WriteString[fileh,"#include \"",filenameprobes,".h\""];*)
							    WriteString[fileh,"\ntypedef double (*dprb_ptr_t)(int,int,int,int,int,double *, double *);\n\n"];
							    
						    (* WriteString[fileh,jacobianfuncheader];*)

							    
							    For[q0=1,q0<=Length[cprobesh],q0++,
								WriteString[fileh,cprobesh[[q0]]];
								];
							    
							    Close[filec];
							    WriteString[filein,StringJoin["\nstatic dprb_ptr_t Delta_PROB_OP[",ToString[Length[listpr]],"]={"]];
							    WriteString[filein,probenames[[1,1]]];
							    For[q0=2,q0<=Length[probenames],q0++,
								WriteString[filein,StringJoin[",",probenames[[q0,1]]]];
								];
							    WriteString[filein,"};\n\n"];
							    WriteString[fileh,"\n#endif"];
							    Close[fileh];
							    
							    ];


probes=VectorialProbes[ToExpression[$ScriptCommandLine[[2]]]];
filein=OpenWrite[StringJoin[$ScriptCommandLine[[3]],".h"]];
WriteString[filein,"#ifndef __TWIINCLUDE \n#define __TWIINCLUDE\n"];
WriteString[filein,StringJoin["\n#include \"",$ScriptCommandLine[[4]],".h\"\n"]];
WriteString[filein,StringJoin["#include \"",$ScriptCommandLine[[5]],".h\"\n"]];
WriteString[filein,"\n#endif \n"];
Close[filein];
filevec=OpenWrite["probevec"];
cPrgen[probes,$ScriptCommandLine[[4]],filevec];
cVargen[probes,sigma,$ScriptCommandLine[[5]],filevec];
Close[filevec];
