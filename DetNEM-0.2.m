(* ::Package:: *)

(* ::Section:: *)
(*Preamble*)


BeginPackage["EM`"];

$EMLibVersion = "DetEM--(EM-5.7 fork)";

$MStepIterLim = 1000;
$EMDebug = False; (* Turns on diagnostic messages for iterative routines *)


(* ::Subsubsection:: *)
(*Details*)


(* ::Text:: *)
(*Deterministic Perturbation Variation on NEM Algorithm*)
(*GMM only for now...*)
(**)


(* ::Subsection:: *)
(*Usage Description*)


(* ::Subsubsection:: *)
(*GMM-(N)EM in 1-D and n-D*)


GmmEM::usage = 
			"GmmEM[data_List,tol_:6,maxIter_:25,options] runs the EM algorithm for a 
mixed Gaussian population model. This version works with 2 populations and a Bernoulli mixing process.
			data_List: the sample dataset under investigation
			tol: the algorithm converged if successive estimations differ by < \!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"tol\"}]]\)
			maxIter: maximum # of iterations allowed
			options: specify Init and Estimand options";


GmmNEM::usage = 
			"GmmNEM[data_List,\[Sigma]n_:0,tol_:6,maxIter_:25,options] runs the EM algorithm for a mixed Gaussian population model. New noise samples added at each iteration. This version works with 2 populations and a Bernoulli mixer.
			data_List: the sample dataset under investigation
			tol: the algorithm converged if successive estimations differ by < \!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"tol\"}]]\)
			maxIter: maximum # of iterations allowed
			options: specify Init and Estimand options";


GmmChaosNEM::usage = "GmmChaosNEM: GMM-NEM with logistic chaotic noise samples";


GmmMultNEM::usage = "GmmMultNEM: GMM-NEM with multiplicative noise samples";


GmmMultNEMneg::usage = "GmmMultNEMneg...";


GmmMultNEMInd::usage = "GmmMultNEMInd: GMM-NEM with *i.i.d* multiplicative noise samples";


GmmNDEMFLaw::usage = "GmmNDEMLaw with full covariance matrices. Uses ...";


(* ::Subsubsection::Closed:: *)
(*Options Declarations*)


CoolingSchedule::usage = 
			"Option to set the noise std.dev's rate of decay in an Annealed-Noise EM procedure. Default is 2 i.e. the noise std.dev \[Sigma]n becomes \[Sigma]n/k^2 on the kth iteration.";
Init::usage = 
			"Option to set the initial parameters \[Theta] for mixture model with 2 populations. For MMs \[Theta]={\[Alpha],\[Theta]_0,\[Theta]_1]} where \[Theta]_k={m_k,d_k} and \[Alpha] is the mixing probability";
Estimand::usage = 
			"Option to set the parameter to estimate in a mixture model.";
Distribution::usage = 
			"Option for generalized mixture model to set the location-scale population distribution. Must be a distribution object...";
NoiseDistribution::usage = 
			"Option for Noisy EM procedures to set additiive noise distribution. Must be a distribution object parametrized by scale parameter...";
Init\[Alpha]::usage = "Initial cluster weights";
Initw::usage = "Initial cluster weights for routines supporting K>2 subpopns";
Init\[Mu]::usage = "Initial cluster centroids";
Init\[Sigma]::usage = "Initial cluster dispersions";

ParameterLock::usage = "3-bit Mask (for {\[Alpha],\[Mu],\[Sigma]}) specifying which parameters to leave alone";
ClusterCount::usage = "# of clusters";
ToleranceLevel::usage = "";
MaxIteration::usage = "";

LStoppingRule::usage = "Use \[CapitalDelta]L (True) for stopping criterion. Default is False i.e. use \[CapitalDelta]\[Theta] stopping rule";
ForceDiagonal::usage = "Assume (force) diagonal covariance matrices for subpopns";


(* ::Section:: *)
(*Implementation*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Gaussian Mixture EM Model*)


Begin["`GmmEM`"];


(* ::Subsubsection::Closed:: *)
(*Auxiliary Functions*)


(* ::Text:: *)
(*Helper functions for gaussian mixture EM*)
(*Reformulated EM routines. Different from versions in pre-parallel-exec*)


(* Parameter Specification:  
\[Theta]={\[Alpha],\[Theta]_ 0,\[Theta]_ 1}
where \[Theta]_k={m_k,d_k}   
K is the number of different popns  *)

pzi[zi:(0|1), \[Alpha]_/;0<=\[Alpha]<=1]:=1-PDF[BernoulliDistribution[\[Alpha]],zi];
fyi[yi_,\[Theta]k_List/;VectorQ[\[Theta]k]]:=PDF[NormalDistribution[\[Theta]k[[1]], \[Theta]k[[2]]],yi];
(* f_y^total (y|\[Theta])=(Overscript[Subscript[\[Sum], k=0], K]\[Alpha]_k*f_k (y|\[Theta]_k)) *)
Fy[yi_,\[Theta]:{_,{_,_},{_,_}}]:=PDF[MixtureDistribution[{\[Theta][[1]],1-\[Theta][[1]]},Table[NormalDistribution[\[Theta]k[[1]], \[Theta]k[[2]]],{\[Theta]k,\[Theta][[2;;]]}]],yi];
FyOLD[yi_,\[Theta]:{_,{_,_},{_,_}}]:=ParallelSum[pzi[k,\[Theta][[1]]]*fyi[yi,\[Theta][[k+2]]],{k,0,1}];

(* In fyzi[]: k+2 index locates the correct popn parameters for the k^th popn *)
fyzi[yi_,zi_,\[Theta]:{_,{_,_},{_,_}}]:=Exp@Sum[DiscreteDelta[zi-k]*(Log@pzi[k,\[Theta][[1]]]+Log@fyi[yi,\[Theta][[k+2]]]),{k,0,1}]; 
fyziNEW[yi_,zi_,\[Theta]:{_,{_,_},{_,_}}]:= (pzi[zi,\[Theta][[1]]] * fyi[yi,\[Theta][[k+2]]]); 

pzjyi[zj:(0|1),yi_,\[Theta]:{_,{_,_},{_,_}}]:=fyzi[yi,zj,\[Theta]]/Fy[yi,\[Theta]];

Lyizi[\[Theta]:{_,{_,_},{_,_}},yi_,zi_]:=ParallelSum[
(Log@pzi[k,\[Theta][[1]]]+Log@fyi[yi,\[Theta][[k+2]]])*DiscreteDelta[zi-k],
{k,0,1}
];
(* E[Lyz[\[Theta];y,Z] | Z,\[Theta]t ] = E[\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\(Lyizi[\[Theta]; 
\*SubscriptBox[\(y\), \(i\)], 
\*SubscriptBox[\(z\), \(i\)]]\)\)|Z,\[Theta]t]
					   =\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\(E[Lyizi[\[Theta]; 
\*SubscriptBox[\(y\), \(i\)], 
\*SubscriptBox[\(z\), \(i\)]] | 
\*SubscriptBox[\(z\), \(i\)], \[Theta]t]\)\)
					  =\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\(
\*SubscriptBox[\(\[Sum]\), \(z = j\)]Lyizi[\[Theta]; 
\*SubscriptBox[\(y\), \(i\)], j]*p[z = j | 
\*SubscriptBox[\(y\), \(i\)], \[Theta]t]\)\);
This summation is the basis for Qver2. Qv1 is a simplification based on delta eliminations.
*)
Qv1[\[Theta]:{_,{_,_},{_,_}}, \[Theta]t:{_,{_,_},{_,_}}, y_List] :=Sum[
(Log@pzi[j,\[Theta][[1]]]+Log@fyi[y[[i]],\[Theta][[j+2]]])*pzjyi[j,y[[i]],\[Theta]t],
{i,Length@y},{j,0,1}
];
Qv2[\[Theta]:{_,{_,_},{_,_}}, \[Theta]t:{_,{_,_},{_,_}}, y_List] :=Sum[
Lyizi[\[Theta],y[[i]],j]*pzjyi[j,y[[i]],\[Theta]t],
{i,Length@y},{j,0,1}
];
(* Main routine: more robust since no weird data structure assumed for \[Theta] *)
Q[\[Theta]_List, \[Theta]t_List,y_List] :=Sum[
Lyizi[{\[Theta][[1]],{\[Theta][[2]],\[Theta][[3]]},{\[Theta][[4]],\[Theta][[5]]}},y[[i]],j]*pzjyi[j,y[[i]],{\[Theta]t[[1]],{\[Theta]t[[2]],\[Theta]t[[3]]},{\[Theta]t[[4]],\[Theta]t[[5]]}}],
{i,Length@y},{j,0,1}
];

(* Clean-Conditioning version of Q() *)
QN[\[Theta]_List, \[Theta]t_List,y_List,y0_List] := Sum[
Lyizi[{\[Theta][[1]],{\[Theta][[2]],\[Theta][[3]]},{\[Theta][[4]],\[Theta][[5]]}},y[[i]],j]*pzjyi[j,y0[[i]],{\[Theta]t[[1]],{\[Theta]t[[2]],\[Theta]t[[3]]},{\[Theta]t[[4]],\[Theta]t[[5]]}}],
{i,Length@y},{j,0,1}
];


(* ::Subsubsection::Closed:: *)
(*GMM-EM: Classical Noiseless Version*)


(* ::Text:: *)
(*Static data implementation*)


(* Current: FixedPoint EM implementation *)
PlaceParam[param_List,new\[Theta]_,loc_]:=ReplacePart[param,MapThread[Rule,{Flatten@{loc},Flatten@{new\[Theta]}}]];
(* Ensures that each objective function changes only parameters it is spec'd to change *)

Options[GmmEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all"};
GmmEM[data_List,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta]},
{init\[Theta],choose\[Theta]}={Init,Estimand}/.{opts}/.Options[GmmEM];
Off[General::"unfl"];
obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{Q[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#,data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{Q[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#,data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{Q[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#,data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{Q[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{Q[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];

(*Quiet@*)FixedPointList[(*obj*)
Check[ obj[#], (Abort[]), {NArgMax::nnum,NMaximize::cvmit,Power::infy,\[Infinity]::indet} ]&,
init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)]
]


(* ::Subsubsection::Closed:: *)
(*Noise Functions: Gaussian Mixture*)


detNEM[data_List,means_List,\[Sigma]_]:=Module[
{valints,nind,n=0*Range@Length@data},
valints = Table[
	IntervalUnion[Interval[{Min[0,Max[ 2(means-data[[j]]) ]],0}],
		Interval[{0,Max[0,Min[ 2(means-data[[j]]) ]]}]
	],
	{j,Length@data}
];
(* Can only positively identify  [0,0]-intervals *)
nind=DeleteCases[Range@Length@data, j_/;valints[[j]]==0];
Do[ 
	n[[j]]= 0.5*\[Sigma]*If[Min[valints[[j]]]==0,
				Max[valints[[j]]],
				Min[valints[[j]]]
			],
	{j,nind}
];

data+n
];


chaosNEM[data_List,means_List,\[Sigma]_]:=Module[
{valints,nind,n=0*Range@Length@data,nbucket,len},
valints = Table[
	IntervalUnion[Interval[{Min[0,Max[ 2(means-data[[j]]) ]],0}],
		Interval[{0,Max[0,Min[ 2(means-data[[j]]) ]]}]
	],
	{j,Length@data}
];

nind=DeleteCases[Range@Length@data, j_/;valints[[j]]==0];
len=Length@data;
nbucket = \[Sigma]*NestList[4#(1-#)&,0.123456789,3*len];
nbucket = nbucket[[RandomSample[Range@Length@nbucket, len]]];

Do[ 
	n[[j]]= 0.5*nbucket[[j]]*If[Min[valints[[j]]]==0,
				Max[valints[[j]]],
				Min[valints[[j]]]
			],
	{j,nind}
];

data+n
];


(* ::Subsubsection::Closed:: *)
(*Noise Functions: Gaussian Mixture + Multiplicative Noise*)


testQ[y_,n_, m_List]:=Apply[And,(y^2 (n^2-1)-2*y*#*(n-1))<=0&/@m]

Sampler[y_,\[Mu]_List,dist_]:=Module[
{ncand=RandomVariate@dist,t=0,tries=100},
While[ And[Not[testQ[y,ncand,\[Mu]]], t++<tries], ncand=RandomVariate@dist ];
Return@If[t<tries,ncand,1]
];(*Default return is 1 since multiplicative*)


multNEM[data_List,means_List,dist_]:=Module[
{valints,n=ConstantArray[1,Length@data]},
	n=Table[
		Sampler[data[[j]], means, dist],
		{j, Length@data}
	];
	data*n
];



negSampler[y_,\[Mu]_List,dist_]:=Module[
{ncand=RandomVariate@dist,t=0,tries=100},
While[ 
	And[testQ[y,ncand,\[Mu]], t++<tries], 
	ncand=RandomVariate@dist 
];
Return@If[t<tries,ncand,1]
];

multNEMneg[data_List,means_List,dist_]:=Module[
{valints,n=ConstantArray[1,Length@data]},
	n=Table[
		negSampler[data[[j]], means, dist],
		{j, Length@data}
	];
	data*n
];


multNEMInd[data_List,means_List,dist_]:=Module[
{valints,n=ConstantArray[1,Length[data]]},
	n=RandomVariate[dist,Length[data]];
	data*n
];


multNEMold[data_List,means_List,dist_]:=Module[
{valints,n=ConstantArray[1,Length@data]},
	n=Table[ 
		Which[
(data[[j]]<0),RandomVariate[TruncatedDistribution[{(2-data[[j]])/data[[j]],1},dist]],
(0<data[[j]]<1),RandomVariate[TruncatedDistribution[{1,(2-data[[j]])/data[[j]]},dist]],
(data[[j]]>1),RandomVariate[TruncatedDistribution[{(2-data[[j]])/data[[j]],1},dist]],
True, 1
]
		,{j,Length@data}
	];

data*n
];


(* ::Subsubsection::Closed:: *)
(*Noisy GMM-EM*)


Options[GmmNEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

GmmNEM[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNEM];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#1, data, data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#1, data, data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#1,detNEM[data,m,#2], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


Options[GmmChaosNEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

GmmChaosNEM[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmChaosNEM];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#1, data, data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#1, data, data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#1,chaosNEM[data,m,#2], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection::Closed:: *)
(*Multiplicative Noise GMM-EM*)


Options[GmmMultNEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[1,#]&)
}; (* Changed: Subscript[\[Mu], n]=1*)

GmmMultNEM[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmMultNEM];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#1, data, data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#1, data, data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#1,multNEM[data,m,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


Options[GmmMultNEMneg]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[1,#]&)
}; (* Changed: Subscript[\[Mu], n]=1*)

GmmMultNEMneg[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmMultNEMneg];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#1, data, data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#1, data, data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#1,multNEMneg[data,m,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection::Closed:: *)
(*IID Multiplicative Noise GMM-EM*)


Options[GmmMultNEMInd]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[1,#]&)
}; (* Change: Subscript[\[Mu], n]=1 *)

GmmMultNEMInd[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmMultNEMInd];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#1, data, data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#1, data, data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#1,multNEMInd[data,m,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#1,data, data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection::Closed:: *)
(*Parallelization Prep*)


(* ::Text:: *)
(*The routines are annoyingly slow in sequential evaluation mode*)


DistributeDefinitions[pzi,fyi,fyzi,Fy,pzjyi,Lyizi,Qv1,Qv2,Q,QN,PlaceParam,
						GaussNEMConditionQ,OutsideMeanClusterQ,
						advNemCond
]
DistributeDefinitions[GmmEM,GmmNEM,GmmNECM,GmmNEMInd]
End[];


(* ::Subsection::Closed:: *)
(*Gaussian N-D Mixture EM with Non-diagonal Covariance*)


Begin["`GmmNDEMFC`"];
(* CAVEAT: 1-D NormalDistribution specs std.dev, N-D MultinormalDistribution specs \[CapitalSigma] matrix i.e. \[Sigma]^2 *)


(* ::Subsubsection::Closed:: *)
(*Auxillary Functions*)


pzi[zi_Integer, w_List/;VectorQ[w]]:=PDF[EmpiricalDistribution[w->Range@Length[w]], zi];

fyi[yi_?VectorQ,\[Theta]k_List]:=PDF[MultinormalDistribution[\[Theta]k[[1]], \[Theta]k[[2]]],yi];

Fy[yi_?VectorQ, \[Theta]:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta]))]:=PDF[
		MixtureDistribution[
			\[Theta][[1]],
			Table[MultinormalDistribution[\[Theta][[2,k]], \[Theta][[3,k]]],{k,Length@\[Theta][[1]]}]
	],yi
]; (* This is the observed likelihood at \[Theta] for a single sample yi*)

(* f(y,z|\[Theta]): *)
fyzi[yi_?VectorQ, zi_Integer, \[Theta]:{{__},{__List},{__List}}] := Product[
(pzi[k,\[Theta][[1]]]*fyi[yi,{\[Theta][[2,k]],\[Theta][[3,k]]}])^DiscreteDelta[zi-k],
{k,Length[\[Theta][[1]]]}
]/;(Equal@@(Length/@\[Theta]))&&(1<=zi<=Length[\[Theta][[1]]]) ;

(* p(z|y,\[Theta]) *)
pzjyi[zj_Integer, yi_,\[Theta]:{{__},{__List},{__List}}]:=
fyzi[yi,zj,\[Theta]]/Fy[yi,\[Theta]]/;(Equal@@(Length/@\[Theta]))&&(1<=zj<=Length[\[Theta][[1]]]);

Ly[samples_List, \[Theta]:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta]))]:=Sum[N@Fy[yt, \[Theta]],{yt,samples}];


(* ::Subsubsection::Closed:: *)
(*Deterministic Noise Function: Gaussian Mixture*)


detNEM[data_List,means_List,\[Sigma]_]:=Module[
{valints,nind,n=0*Range@Length@data},
valints = Table[
	IntervalUnion[Interval[{Min[0,Max[ 2(means-data[[j]]) ]],0}],
		Interval[{0,Max[0,Min[ 2(means-data[[j]]) ]]}]
	],
	{j,Length@data}
];
(* Can only positively identify  [0,0]-intervals *)
nind=DeleteCases[Range@Length@data, j_/;valints[[j]]==0];
Do[ 
	n[[j]]= 0.5*\[Sigma]*If[Min[valints[[j]]]==0,
				Max[valints[[j]]],
				Min[valints[[j]]]
			],
	{j,nind}
];

data+n
];(* This version is a bit faster. But it uses a lot more memory in the interval struct *)


(* *)
ndDetNEM[data:{__List},means:{__List},\[Sigma]_/;\[Sigma]<=0]:=data;
ndDetNEM[data:{__List},means:{__List},\[Sigma]_/;\[Sigma]>0]:=Module[{dim = Length@data[[1]],d=ConstantArray[0,Dimensions[data]]},
Table[d[[;;,col]]=detNEM[data[[;;,col]],means[[;;,col]],\[Sigma]],{col,dim}];
d
]


(* ::Subsubsection:: *)
(*N-D Learning Law EM  with Full Covariance Matrix*)


Options[GmmNDEMFLaw]:={
Init\[Alpha]->ConstantArray[1/3,3], Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->Table[IdentityMatrix[2],{3}],
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
ClusterCount->3,ParameterLock->0,
ForceDiagonal->False,LStoppingRule->False
};

GmmNDEMFLaw[data_List,\[Sigma]n_?NonNegative,tol_:6,maxIter_:\[Infinity],opts___?OptionQ]:=
Module[{obj,init\[Alpha],init\[Mu],init\[Sigma],K,dim=Length[data[[1]]],cs,ndist,k=1,plock,recalc\[Sigma]=\[Sigma]n,srule,white\[Sigma]},
{init\[Alpha],init\[Mu],init\[Sigma],K,cs,ndist,white\[Sigma]}={
	Init\[Alpha],Init\[Mu],Init\[Sigma],ClusterCount,CoolingSchedule,NoiseDistribution,ForceDiagonal
}/.{opts}/.Options[GmmNDEMFLaw];

{plock}={ParameterLock}/.{opts}/.Options[GmmNDEMFLaw];
plock = (#==1)&/@IntegerDigits[plock,2,3]; (* Mask-style parameter lock {a,m,s} *)
If[And@@plock, Return[{init\[Alpha],init\[Mu],init\[Sigma]}]];  (* Nothing to update *)

If[Length[init\[Alpha]]<K, init\[Alpha] = ConstantArray[1,K]/K];

obj[cur\[Theta]_,\[Sigma]_]:=Module[{tmp=cur\[Theta],pzjsums,smptmp,zjs=Range[Length@First@cur\[Theta]]},

		If[$EMDebug,
			Print[
				"Likelihood @iter#",(k),"= ", Ly[data, cur\[Theta]],"\n",
				Grid[Transpose@{cur\[Theta][[1]],cur\[Theta][[2]],(*MatrixForm/@*)cur\[Theta][[3]]},Frame->All]
			]
		];
			pzjsums=Total@Table[(pzjyi[#,sample,cur\[Theta]])&/@zjs,{sample,data}];
			If[plock[[1]],
				tmp[[1]]=init\[Alpha],
				tmp[[1]]=(pzjsums/Length@data);
				tmp[[1]]=tmp[[1]]/Total@tmp[[1]];
			];

			tmp[[2]]=If[plock[[2]],
						init\[Mu],
						Total@Table[
							(sample*pzjyi[#,sample,cur\[Theta]])&/@zjs,
							{sample,data}]/pzjsums
						];
			tmp[[3]]=If[plock[[3]],
						init\[Sigma],
						smptmp=data;
						smptmp = ndDetNEM[data,init\[Mu],\[Sigma]];
						(Total@Table[
							(pzjyi[#,sample,cur\[Theta]]*
							Outer[Times,
								(sample-cur\[Theta][[2,#]]),
								(sample-cur\[Theta][[2,#]])
								]
							)&/@zjs,
						{sample,smptmp}])/pzjsums
					];
			If[white\[Sigma],tmp[[3]]=(DiagonalMatrix/@Diagonal/@tmp[[3]])];
			tmp
	];

If[$EMDebug,
	Print["Starting Law FP iters.\nCluster Count = ", K,
		"\nInitial Parameters:\n",
		Grid[Transpose@{init\[Alpha],init\[Mu],(*MatrixForm/@*)init\[Sigma]},Frame->All]
	]
];

srule = If[Not[(LStoppingRule/.{opts}/.Options[GmmNDEMFLaw])],
			((Max[Abs[Flatten@#1-Flatten@#2]]<10^-tol)&),
			((Abs[Ly[data,#1]-Ly[data,#2]]<10^-tol)&)
];

Quiet@FixedPointList[
Check[ obj[#,\[Sigma]n/(k++)^cs], (Abort[]), 
	{General::unfl, Power::infy,\[Infinity]::indet}(*MultinormalDistribution ::posdefprm,DiagonalMatrix::mindet,*)
]&,
{init\[Alpha],init\[Mu],init\[Sigma]},
maxIter,
SameTest->srule
]
]


(* ::Subsubsection::Closed:: *)
(*End*)


End[];


(* ::Section::Closed:: *)
(*End*)


End[];
EndPackage[];
