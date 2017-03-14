(* ::Package:: *)

(* ::Section::Closed:: *)
(*Preamble*)


BeginPackage["SR`"];

$SRLibVersion = "SR-v0.4.1";


(* ::Subsubsection::Closed:: *)
(*ChangeLog & TODO*)
(**)


(* ::Text:: *)
(*Changes (0.4.1): *)
(*- Introduced variable data trim for SRProfileB.*)
(*	- Specify TrimLength Option = # of points to drop from upper and lower data limits (default = 2)*)
(*- Dropped obsolete generic parser spec.*)
(**)
(**)
(*Changes (0.4.0): *)
(*- Developed SRSimulateD: a generalization of the very successful SRSimulateC:*)
(*	- User-defined summarizing functions with matching log-files and export formats*)
(*	- Better restart tracking by using running index of completed iterations*)
(*	- Changed noise input from generating function form to simple scalar. In line with typical use*)
(**)


(* ::Text:: *)
(*Todo:*)
(*- Error reporting can only be compared to the true value of the estimand. But EM is*)
(*converges only to the ML estimate. Usually MLE != true value. Therefore error reporting*)
(*is useless for relative error calculation*)
(*	- Temp Solution: For large training set, the zero-noise case should average to the MLE.*)
(*	 => Use avg zero-noise estimate as surrogate for MLE.*)
(*		- Not good enough: The point is to see if noise improves zero-noise results*)
(*	=> should not baseline on zero-noise results. *)
(*	=> Estimation error comparison still unresolved... *)
(*	=> Temp bypass: compare to true value instead of ML estimate*)
(**)


(* ::Subsubsection::Closed:: *)
(*Function Usage Description*)


SRSimulateC::usage = 
			"SRSimulateC[testfxn_,dataGen_,noiseGen_,train_,\[Sigma]1_,\[Sigma]2_,d\[Sigma]_, opts___] is a higher-level 
testing routine that simulates a given iterative statistical algorithm with varied amounts of noise 
injected into the test function. The test function is a 2-argument pure-function. The 1st argument 
is the data generator and the 2nd is the noise generator. This procedure allows for exotic data/noise
 dependencies including non-additive noise models. The magic bit is in how you define testfxn to use 
the noise parameter in its second argument.";

SRSimulateD::usage = 
			"SRSimulateD generalizes SRSimulateC. We lose the noise generator and just tell the test 
function what noise stddev to use. This version also has generic data loggers and better simulation recovery ";


(* ::Subsubsection::Closed:: *)
(*Plotter Usage Description*)


SRProfile::usage = 
			"SRProfile[data_List,index_List,n_Integer:20] plots profile of output. 
data is numerical record of output of an SR procedure. 
Index is noise dispersion indexing list. 
n is the # number of training results to omit from upper and lower extremes.";

SRProfileB::usage = 
			"SRProfileB[data_List,index_List,bn_Integer:20] plots profile of output. 
It generates estimates for the data average by bootstrapping the results at each noise level. 
Useful for converting integer-numbered results into rational-numbered results
data is the numerical record of output of an SR procedure. 
Index is noise dispersion indexing list. 
bn = {bn1, bn2}: bn1 is the # number of bootstrap estimates to acquire per noise level. 
bn2 is the size of each ensemble used to calculate each bootstrap estimates.";


(* ::Subsubsection::Closed:: *)
(*Options Usage Description*)


Train::usage = 
			"Option that sets the # of iterations to run";

TimeLimit::usage = 
			"Option that sets the max time (in seconds) allowed per iteration";

RecoveryFile::usage = 
			"String Option that sets up disk file for temporary results backup. Best specify relative to Research dir. ";

InterpolationOrd::usage = 
			"Integer Option for SR plotter. Determines curve fit type.";

MedianSummary::usage = 
			"Boolean Option for SR plotter. Choose to use median instead of trimmed mean to analyze data.";

Origin::usage = 
			"(x,y)-Coordinates for Axes origin.";

RecoveryFile::usage = 
			"Specifies prefix for intermediate result storage files...";

SaveSuffix::usage = 
			"Specifies suffic for final result storage files. Set to empty string if you don't want final results saved to disk.";

PrimaryColor::usage="";
BandColor::usage="";
TrendColor::usage="";
TrimLength::usage="";

ParseFunctions::usage = 
			"Option that specifies parsing functions for each iteration";

RecordLabels::usage = 
			"Option that specifies labels for output files";

RecordFormats::usage = 
			"Options that specifies export formats for each record";



(* ::Section:: *)
(*Implementation*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Setup*)


(* ::Text:: *)
(*Messages and Options definitions...*)


SRSimulateC::convstop = "Convergence Stopped for \[Sigma] = `1`";
SRSimulateD::convstop = "Convergence Stopped for \[Sigma] = `1`";

Options[SRSimulateC]:={TimeLimit->360,Train->32,RecoveryFile->"$SRC-Recovery",SaveSuffix->"-Final"};
Options[SRSimulateD]:={TimeLimit->360,Train->32,RecoveryFile->"$SRC-Recovery",SaveSuffix->"-Final",
ParseFunctions->{Length,Last},RecordLabels->{"Speeds","Estimates"},RecordFormats->{"Table","MAT"}
};

Options[SRProfile]:={InterpolationOrd->1,MedianSummary->False, Origin->Automatic,PrimaryColor->Blue,BandColor->Green,TrendColor->Red};
Options[SRProfileB]:={InterpolationOrd->1,MedianSummary->False, Origin->Automatic,PrimaryColor->Blue,BandColor->Green,TrendColor->Red,TrimLength->2};


(* ::Subsection:: *)
(*SR Procedures*)


(* ::Subsubsection::Closed:: *)
(*SRSimulateC :*)


(* ::Text:: *)
(*This procedure allows for exotic data/noise dependencies including non-additive noise models. It is going to give some freaking sweet results!*)


(*  *)
SRSimulateC[testfxn_,dataGen_,noiseGen_,\[Sigma]1_,\[Sigma]2_,d\[Sigma]_,opts___?OptionQ]:=
Module[
{\[Sigma],testResult,counter,tofill,rcounter,update,train,timelimit,storage,final,speedsR,estimatesR},

{train,timelimit,storage,final}={Train,TimeLimit,RecoveryFile,SaveSuffix}/.{opts}/.Options[SRSimulateC];
DistributeDefinitions[dataGen,noiseGen,testfxn,timelimit,train,\[Sigma]1,\[Sigma]2,d\[Sigma]];

\[Sigma]=\[Sigma]1;counter=1;rcounter=1;
estimatesR=ConstantArray[0,{Length@Range[\[Sigma]1,\[Sigma]2,d\[Sigma]],train}];(* 2- or 3-Dimensions *)
speedsR=ConstantArray[0,{Length@Range[\[Sigma]1,\[Sigma]2,d\[Sigma]],train}]; (* Always 2-Dimensions *)

(*update[temp_]:=(speedsR[[rcounter, counter]]=Length@temp;estimatesR[[rcounter, counter(*++*)]]=Last@temp;);*)
tofill[]:=Flatten@Position[speedsR[[rcounter]], 0];
update[temp_]:=With[{mcounter = If[Length@tofill[]!=0, First@tofill[],Mod[counter,train,1]]},	
	speedsR[[rcounter, mcounter]]=Length@temp;
	estimatesR[[rcounter, mcounter]]=Last@temp;
];
(* Solves 2 problems:
- Makes the results row vector a circular buffer to fix Set::partw errors after convstop  
- Fixes 0-speed logging errors
*)

SetSharedVariable[\[Sigma],rcounter,counter];
SetSharedFunction[update,tofill];

If[FileExistsQ[storage<>"Speeds"],
speedsR = Import[storage<>"Speeds","Table"];
estimatesR = Import[storage<>"Estimates","XLS"];
If[Length@estimatesR==1, estimatesR = First@estimatesR];
rcounter=First@Select[Range@Length@speedsR,Not@FreeQ[speedsR[[#]],0]&];
\[Sigma]=Part[Range[\[Sigma]1,\[Sigma]2,d\[Sigma]],rcounter]
Print["Resuming from \[Sigma] = ", \[Sigma]];
];

Monitor[
While[\[Sigma]<=\[Sigma]2,
counter=1;
While[(*(counter<=train)&&*)(Length@tofill[]>0),(* Run sim until min train size reached *)
ParallelDo[
Off[ParallelSum::"subpar"];SRSimulateC::convstop = "Convergence Stopped for \[Sigma] = `1`";

testResult=If[\[Sigma]!=0,
TimeConstrained[testfxn[dataGen[],noiseGen[\[Sigma]]],timelimit],
TimeConstrained[testfxn[dataGen[],(0*noiseGen[1])],timelimit]
];

If[testResult=!=$Aborted,
(update[testResult];counter++;),
Message[SRSimulateC::convstop,\[Sigma]];
],
{train}
];
];
\[Sigma]+=d\[Sigma];
rcounter++;
(* Order of export matters in case of abort *)
Export[storage<>"Estimates",estimatesR,"XLS"];
Export[storage<>"Speeds",speedsR,"Table"]
],
Grid[{
{"Completed Tests",counter-1,ProgressIndicator[counter-1,{0,train-1}]},
{"@ Noise Level",N[\[Sigma]],ProgressIndicator[\[Sigma],{\[Sigma]1,\[Sigma]2}]}
}]
];
UnsetShared[\[Sigma],rcounter,counter,update];

If[And[final!="", Or[FileExistsQ[storage<>"Estimates"<>final], FileExistsQ[storage<>"Speeds"<>final]] ],
	final = "-"<>DateString[{"Day","MonthNameShort","YearShort","-","Hour24","Minute"}]<>"h"<>final
];

If[final!="",
Export[storage<>"Estimates"<>final,estimatesR,"XLS"];
Export[storage<>"Speeds"<>final,speedsR,"Table"]
];
DeleteFile[storage<>"Speeds"];
DeleteFile[storage<>"Estimates"];

{speedsR,estimatesR}
]


(* ::Subsubsection:: *)
(*SRSimulateD :*)


(*  *)
SRSimulateD[testfxn_,dataGen_,\[Sigma]1_,\[Sigma]2_,d\[Sigma]_,opts___?OptionQ]:=
Module[
{\[Sigma],testResult,counter,tofill,rcounter,update,train,
timelimit,storage,final,speedsR,estimatesR,parsers,
reclabels,index,records,recforms},

{train,timelimit,storage,final}={Train,TimeLimit,RecoveryFile,SaveSuffix}/.{opts}/.Options[SRSimulateD];
{parsers,reclabels,recforms} = {ParseFunctions,RecordLabels,RecordFormats}/.{opts}/.Options[SRSimulateD];
DistributeDefinitions[dataGen,testfxn,timelimit,train,\[Sigma]1,\[Sigma]2,d\[Sigma]];

\[Sigma]=\[Sigma]1;counter=1;rcounter=1;
records = ConstantArray[0,{Length@parsers,Length@Range[\[Sigma]1,\[Sigma]2,d\[Sigma]],train}];
index = ConstantArray[0,{Length@Range[\[Sigma]1,\[Sigma]2,d\[Sigma]],train}];

tofill[]:=Flatten@Position[index[[rcounter]], 0];
update[temp_]:=With[{mcounter = If[Length@tofill[]!=0, First@tofill[],Mod[counter,train,1]]},	
	index[[rcounter, mcounter]]=1;
	Table[records[[k,rcounter, mcounter]]=(parsers[[k]])@temp,
		{k,Length@parsers}
	];
];

SetSharedVariable[\[Sigma],rcounter,counter];
SetSharedFunction[update,tofill];

If[FileExistsQ[storage<>"Index"],
	index = Import[storage<>"Index","Table"];
	Do[records[[k]]=Import[storage<>reclabels[[k]], recforms[[k]] ],{k,Length@parsers}];
	rcounter=First@Select[Range@Length@index,Not@FreeQ[index[[#]],0]&];
	\[Sigma]=Part[Range[\[Sigma]1,\[Sigma]2,d\[Sigma]],rcounter];
	Print["Resuming from \[Sigma] = ", \[Sigma]];
];

Monitor[
While[\[Sigma]<=\[Sigma]2,
counter=1;
While[(Length@tofill[]>0),(* Run sim until min train size reached *)
ParallelDo[
Off[ParallelSum::"subpar"];SRSimulateD::convstop = "Convergence Stopped for \[Sigma] = `1`";

testResult=TimeConstrained[testfxn[dataGen[],\[Sigma]],timelimit];

If[testResult=!=$Aborted,
(update[testResult];counter++;),
Message[SRSimulateD::convstop,\[Sigma]];
],
{train}
];
];

(* Order of export matters in case of abort *)
Do[
Export[storage<>reclabels[[k]],records[[k,1;;rcounter]], recforms[[k]]],{k,Length@parsers}];
Export[storage<>"Index",index,"Table"];

\[Sigma]+=d\[Sigma];rcounter++;
],
Grid[{
{"Completed Tests",counter-1,ProgressIndicator[counter-1,{0,train-1}]},
{"@ Noise Level",N[\[Sigma]],ProgressIndicator[\[Sigma],{\[Sigma]1,\[Sigma]2}]}
}]
];
UnsetShared[\[Sigma],rcounter,counter,update];

If[And[final!="", FileExistsQ[storage<>"Index"<>final]],
	final = "-"<>DateString[{"Day","MonthNameShort","YearShort","-","Hour24","Minute"}]<>"h"<>final
];

If[final!="",
Table[Export[storage<>reclabels[[k]]<>final,records[[k]], recforms[[k]]],{k,Length@parsers}];
Export[storage<>"Index"<>final,index,"Table"]
];
DeleteFile[storage<>"Index"];
Do[DeleteFile[storage<>reclabels[[k]]],{k,Length@parsers}];

records
]


(* ::Subsection::Closed:: *)
(*Output Plotter & Parsers*)


(* ::Subsubsection::Closed:: *)
(*Plotter: With Bootstrapped Confidence Band Estimation*)


SRProfileB[data_List, index_List, bn_List:{200,53}, opts___?OptionQ]:=
Module[
	{fulldata,bdata,trend,interpOrd,median,orig,pcol,bcol,tcol,droppts},
	{interpOrd,median,orig,pcol,bcol,tcol,droppts}={InterpolationOrd,MedianSummary,Origin,PrimaryColor,BandColor,TrendColor,TrimLength}/.{opts}/.Options[SRProfileB];
	fulldata = (Sort/@data)[[All,(droppts+1);;-(droppts+1)]]; (*fulldata = (Sort/@data)[[All,3;;-3]];*) 
	(* Drop off #droppts [U,L]-endpoints to deal with extreme outliers *)

	bdata=Table[
		Table[N@Mean@RandomChoice[fulldata[[k]],bn[[2]]],{bn[[1]]}],
		{k,Length@fulldata}
		];
	bdata = Transpose[Sort/@bdata];
	(* ~ 95% confidence band *)
	bdata = bdata[[Ceiling[0.025*bn[[1]]];;Floor[0.975*bn[[1]]],;;]];  

	trend=If[median,
			Median[Transpose[fulldata]],
			Mean[Transpose[fulldata]]
		];

	Show[
		Table[
			ListPlot[Transpose@{index,bdata[[k]]},
				PlotStyle->{bcol,PointSize[Small]}],
			{k,Length@bdata}
		],(* Dots *)
		ListPlot[Transpose@{index[[1;;Length@trend]],trend},
			PlotStyle->{pcol,Thickness[0.004]},Joined->True,InterpolationOrder->interpOrd],
		(* Trend Line *)
		ListPlot[Table[{index[[k]],trend[[1]]},{k,Length@trend}],
			PlotStyle->{tcol,Dashed},Joined->True], 
		(* Noise-free Threshold *)
		AxesOrigin->orig,
		PlotRange->Automatic
	]
]


(* ::Subsubsection::Closed:: *)
(*Plotter: Non-Bootstrap Aggregation*)


(* ::Text:: *)
(*Trimmed minmax SR profile plotter:*)


SRProfile[data_List,index_List,n_Integer:1,opts___?OptionQ]:=
Module[{fullset=Min[Length/@data],fulldata,trend,median,interpOrd,orig,pcol,bcol,tcol},
	{interpOrd,median,orig,pcol,bcol,tcol}={InterpolationOrd,MedianSummary,Origin,PrimaryColor,BandColor,TrendColor}/.{opts}/.Options[SRProfile];
	fulldata=Transpose[Sort/@data[[All,1;;fullset]]];
	trend=If[median,
			Median[fulldata],
			TrimmedMean[fulldata,Min[0.49,n/Length@fulldata[[All,1]]]]
			];
	Show[
		Table[
			ListPlot[Transpose@{index[[1;;Length@fulldata[[k]]]],fulldata[[k]]},
				PlotStyle->{bcol,PointSize[Small]}],
			{k,1+n,Length@fulldata-n}
		],
		ListPlot[Transpose@{index[[1;;Length@trend]],trend},
			PlotStyle->{pcol,Thickness[0.004]},Joined->True,InterpolationOrder->interpOrd],
		ListPlot[Table[{index[[k]],trend[[1]]},{k,Length@trend}],
			PlotStyle->{tcol,Dashed},Joined->True], 
		AxesOrigin->orig,
		PlotRange->Automatic
	]
] 


(* ::Subsection::Closed:: *)
(*Parallelization Prep*)


(* ::Text:: *)
(*- If I played my cards right, no prep necessary... These routines are supposed "wrap" around test-functions. *)
(*So SR routines are technically all overhead i.e. Master Kernel processing. *)
(*Parallelization is just for running multiple instances of the test. *)
(*- Parsing functions might require distribution if they are executed outside the locking-block.*)


(* ::Section::Closed:: *)
(*End*)


End[];
(*Protect[SRSimulateD,SRSimulateC];*)
EndPackage[];
