(* ::Package:: *)

(* ::Section:: *)
(*Preamble*)


BeginPackage["EM`"];

$EMLibVersion = "EM-v0.5.7.1";

$MStepIterLim = 1000;
$EMDebug = False; (* Turns on diagnostic messages for iterative routines *)


(* ::Subsubsection:: *)
(*Changes in version 5.6 & 5.7*)


(* ::Text:: *)
(*CAVEAT: 1-D NormalDistribution specs std.dev, N-D MultinormalDistribution specs \[CapitalSigma] matrix i.e. (\[Sigma]^2)!!!*)
(*Also: GmmNDEM behaves erratically. Law version works ok though.*)
(**)
(*0. - Old GMM-NEM condition left out factor of 2. => Old NEM sets half as small than they should be...*)
(*1. - GmmNDEMFLaw works with full covariance matrices \[CapitalSigma] [Done]*)
(*2. - Law EM routines: \[CapitalSigma]-updates use \[Mu](t) instead of \[Mu](t+1) [Done]*)
(**)


(* ::Subsection:: *)
(*Usage Description*)


(* ::Subsubsection::Closed:: *)
(*EM functions declarations for Gamma distributions*)


GammaEM::usage = 
			"GammaEM[data_List,\[Alpha]_,T_,m_:0,tol_:6, maxIter_:25] runs the censored-EM algorithm for a 
gamma distributed data set using a Fixed-Point implementation. The estimand is the \[Theta] of the gamma pdf. 
The last 2 parameters determine how to stop the algorithm.
			data_List: the sample dataset under investigation
			\[Alpha]: the known alpha parameter of the gamma distribution
			T: data censorship point
			m: # of uncorrupted points at the beginning of the dataset.
			tol: The algorithm converged if successive estimations differ by < \!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"tol\"}]]\)
			maxIter: maximum # of iterations allowed";


GammaNEM::usage = 
			"GammaNEM[data, \[Alpha], T, m, \[Sigma]n, tol, maxIter] runs the censored-EM algorithm for a 
gamma distributed data set using a Fixed-Point implementation. The estimand is the \[Theta] of the gamma pdf. New noise samples added at each iteration
The last 2 parameters determine how to stop the algorithm.
			data_List: the sample dataset under investigation
			\[Alpha]: the known alpha parameter of the gamma distribution
			T: data censorship point
			m: # of uncorrupted points at the beginning of the dataset.
			tol: The algorithm converged if successive estimations differ by < \!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"tol\"}]]\)
			maxIter: maximum # of iterations allowed";

GammaMultiplicativeNoiseNEM::usage = "";


(* ::Subsubsection:: *)
(*EM functions declarations for Mixtures with Gaussian subpopulations*)


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

GmmNEMInd::usage = "<<fill this later>>";

GmmNEMneg::usage = "<<fill this later>>";

GmmNECM::usage = "Noisy Expectation Conditional Maximization <<fill out details later>>";


(* ::Subsubsection::Closed:: *)
(*EM functions declarations for Mixtures with Cauchy subpopulations*)


CmmEM::usage = 
			"CmmEM[data_List,tol_:6,maxIter_:25,options] runs the EM algorithm for a 
mixed Cauchy population model. This version works with 2 populations and a Bernoulli mixing process.
			data_List: the sample dataset under investigation
			tol: the algorithm converged if successive estimations differ by < \!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"tol\"}]]\)
			maxIter: maximum # of iterations allowed
			options: specify Init and Estimand options";

CmmNEM::usage = "CmmNEM[data_List,\[Sigma]n_:0,tol_:6,maxIter_:25,options] runs the EM algorithm for a 
mixed Cauchy population model. New noise samples added at each iteration. This version works with 2 populations and a Bernoulli mixer.
			data_List: the sample dataset under investigation
			tol: the algorithm converged if successive estimations differ by < \!\(\*SuperscriptBox[\"10\", RowBox[{\"-\", \"tol\"}]]\)
			maxIter: maximum # of iterations allowed
			options: specify Init and Estimand options";


(* ::Subsubsection::Closed:: *)
(*EM functions declarations for Multi-K Mixtures of Multi-Dimensional subpopulations*)


GmmNDEM::usage = "GmmNDEM runs the EM algorithm for a multidimensional Gaussian mixture of an arbitrary number of populations. Uses numerical maximizer instead of analytic GMM-EM update equations.";
GmmNDECM::usage = "GmmNDECM runs the ECM algorithm for a multidimensional Gaussian mixture of an arbitrary number of populations.";
GmmNDEMLaw::usage = "Full (N)EM routine for multidimensional, multi-cluster-count Gaussian Mixture models. ParameterLock bitvector mask ({\[Alpha],\[Mu],\[Sigma]}) specifies which parameters to estimate. Users must specify number of clusters and appropriately-sized initializations for all 3 parameter types. Assumes diagonal covariance (indep. subpopns). Not necessarily white...";
GmmNDEMFLaw::usage = "GmmNDEMLaw with full covariance matrices. Uses ...";


(* ::Text:: *)
(*Functions for Clustering Applications:*)


EMClassifier::usage = "Assigns cluster labels to samples in a data list based on supplied GMM parameter";
MatchLabelsToRef::usage = "Returns most charitable cluster-label translation for comparison to a reference label list.
				Necessary because clustering algorithms are 'invariant' under a permutations of the cluster-labels";
NEMClusterError::usage = "NEM cluster-error measurement procedure. Return Hamming mismatch of classfication and EM params used in comparisons..";
NEMLClusterError::usage = "Same as above. But uses seemingly fast learning law NEM version.";
NEMClusterErrorGTref::usage = "Same as NEMClusterError except compare to ground truth instead of final EM classfication. returns Hamming mismatch p-value KS-GOF test";

FindSoftClusters::usage = "Uses EMClassifier cluster labels to partition data set";
EMClustering::usage = "Partitions data set into clusters after finding EM estimate for GMM parameters";



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


(* ::Subsection::Closed:: *)
(*Gamma EM Model*)


Begin["`GammaEM`"];


(* ::Subsubsection::Closed:: *)
(*Auxillary Functions*)


(* ::Text:: *)
(*Regular and Clean Conditioning helper functions*)


Elnz[\[Alpha]_,\[Theta]_,T_]:=MeijerG[{{},{1,1}},{{0,0,\[Alpha]},{}},T/\[Theta]]/Gamma[\[Alpha],T/\[Theta]]+Log[T]
Lz[\[Alpha]_,\[Theta]_,z_List]:=ParallelSum[(\[Alpha]-1) Log[z[[i]]]-z[[i]]/\[Theta]-\[Alpha] Log[\[Theta]]-Log[Gamma[\[Alpha]]],{i,1,Length[z]}];
\[Mu]j[xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,xj,xj>=T,(\[Theta]t Gamma[\[Alpha]t+1,T/\[Theta]t])/Gamma[\[Alpha]t,T/\[Theta]t]];
\[Mu]jLn[xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,Log[xj],xj>=T,Elnz[\[Alpha]t,\[Theta]t,T]];

Q[\[Theta]_,\[Theta]t_,x_List,\[Alpha]i_,m_,T_]:= Lz[\[Alpha]i,\[Theta],x[[1;;m]]]+ParallelSum[-(\[Mu]j[x[[i+m]],\[Alpha]i,\[Theta]t,T]/\[Theta])+(\[Alpha]i-1) \[Mu]jLn[x[[i+m]],\[Alpha]i,\[Theta]t,T]-\[Alpha]i Log[\[Theta]]-Log[Gamma[\[Alpha]i]],{i,1,Length[x[[m+1;;All]]]}];

(* Clean Conditioning Procedures *)
\[Mu]jN[xnj_,xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,xnj,xj>=T,(\[Theta]t Gamma[\[Alpha]t+1,T/\[Theta]t])/Gamma[\[Alpha]t,T/\[Theta]t]];
\[Mu]jLnN[xnj_,xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,Log[xnj],xj>=T,Elnz[\[Alpha]t,\[Theta]t,T]];

QN[\[Theta]_,\[Theta]t_,xn_List,x0_List,\[Alpha]i_,m_,T_]:= Lz[\[Alpha]i,\[Theta],x0[[1;;m]]] + ParallelSum[-(\[Mu]jN[xn[[i+m]],x0[[i+m]],\[Alpha]i,\[Theta]t,T]/\[Theta])+(\[Alpha]i-1) \[Mu]jLnN[xn[[i+m]],x0[[i+m]],\[Alpha]i,\[Theta]t,T]-\[Alpha]i Log[\[Theta]]-Log[Gamma[\[Alpha]i]],{i,1,Length[x0[[m+1;;All]]]}];


(* ::Subsubsection::Closed:: *)
(*EM Routine: Classical Noiseless Version*)


(* ::Text:: *)
(*Fixed-Point implementation*)


(* Current: FixedPoint EM implementation *)
Options[GammaEM]:={Init:>Abs@RandomVariate@NormalDistribution[0,0.5]};

GammaEM[data_List,\[Alpha]_,T_,m_:0,tol_:6, maxIter_:25,opts___?OptionQ]:=Module[{\[Theta],init},
{init}={Init}/.{opts}/.Options[GammaEM];
FixedPointList[ 
	Check[
		ArgMax[{Q[\[Theta],#,data,\[Alpha],m,T],\[Theta]>0},\[Theta]], 
		(PrintTemporary["M-Step Error..."];Abort[]), 
		{NMaximize::cvmit, ArgMax::objv}
	]&,
	init,
	maxIter,
	SameTest->(Abs[#1-#2]<10^-tol&)
]
]



(* ::Subsubsection::Closed:: *)
(*Noise Auxillary Functions*)


(* Can't clip at Zero!!! *)
Options[NoiseModel]:={NoiseDistribution->(NormalDistribution[0,#]&)};

NoiseModel[data_List,\[Sigma]n_:0.1,opts___?OptionQ]:= Module[{ndist},
{ndist} = {NoiseDistribution}/.{opts}/.Options[NoiseModel];
Max[10000*$MachineEpsilon,#]&/@(data+RandomVariate[ndist[\[Sigma]n],Length@data])
]


(* ::Subsubsection::Closed:: *)
(*Noisy Gamma EM*)


Options[GammaNEM]:={CoolingSchedule->2, NoiseDistribution->(NormalDistribution[0,#]&), Init:>Abs@RandomVariate@NormalDistribution[0,0.5]};

GammaNEM[data_List,\[Alpha]_,T_,m_:0,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=Module[
{\[Theta],k=1,cs,ndist,init},
{cs,ndist,init}={CoolingSchedule,NoiseDistribution,Init}/.{opts}/.Options[GammaNEM];

If[\[Sigma]n!=0,
FixedPointList[
	Check[ 
		NArgMax[{QN[\[Theta],#,NoiseModel[data,\[Sigma]n/(k++)^cs,NoiseDistribution->ndist],data,\[Alpha],m,T],\[Theta]>0},\[Theta]], 
		(PrintTemporary["M-Step Error..."];Abort[]), 
		{NArgMax::cvmit,NMaximize::cvmit(*, ArgMax::objv*)} 
	]&,
	init,
	maxIter,
	SameTest->(Abs[#1-#2]<10^-tol&)
],
GammaEM[data,\[Alpha],T,m,tol, maxIter]
]
]



(* ::Subsubsection::Closed:: *)
(*Noisy Gamma EM using Multiplicative Noise*)


Options[MultNoiseModel]:={NoiseDistribution->(GammaDistribution[#^-2,#^2]&)};

MultNoiseModel[data_List,\[Sigma]n_:0.1,opts___?OptionQ]:= Module[{ndist},
{ndist} = {NoiseDistribution}/.{opts}/.Options[MultNoiseModel];
Inner[Times,data,RandomVariate[ndist[\[Sigma]n],Length@data],List]
]



Options[GammaMultiplicativeNoiseNEM]:={CoolingSchedule->2, NoiseDistribution->(GammaDistribution[#^-2,#^2]&),Init:>1};

GammaMultiplicativeNoiseNEM[data_List,\[Alpha]_,T_,m_:0,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=Module[
{\[Theta],k=1,cs,ndist,init},
{cs,ndist,init}={CoolingSchedule,NoiseDistribution,Init}/.{opts}/.Options[GammaMultiplicativeNoiseNEM];

If[\[Sigma]n!=0,
FixedPointList[
	Check[ 
		ArgMax[{QN[\[Theta],#,MultNoiseModel[data,\[Sigma]n/(k++)^cs,NoiseDistribution->ndist],data,\[Alpha],m,T],0<\[Theta]<\[Infinity]},\[Theta]], 
		(PrintTemporary["M-Step Error..."];Abort[]), 
		{NArgMax::cvmit,NMaximize::cvmit, ArgMax::objv} 
	]&,
	init, 
	maxIter,
	SameTest->(Abs[#1-#2]<10^-tol&)
],
GammaEM[data,\[Alpha],T,m,tol, maxIter,Init->init]
]
]


(* ::Subsubsection::Closed:: *)
(*Parallelization Prep*)


(* ::Text:: *)
(*The routines are annoyingly slow in sequential evaluation mode*)


DistributeDefinitions[Elnz,Lz,\[Mu]j,\[Mu]jLn,Q];
DistributeDefinitions[\[Mu]jN,\[Mu]jLnN,QN];
DistributeDefinitions[GammaEM];
End[]; 


(* ::Subsection::Closed:: *)
(*Gamma Full-EM Model*)


Begin["`GammaFullEM`"];


(* ::Subsubsection:: *)
(*Auxillary Functions*)


(* ::Text:: *)
(*Regular and Clean Conditioning helper functions*)


Elnz[\[Alpha]_,\[Theta]_,T_]:=MeijerG[{{},{1,1}},{{0,0,\[Alpha]},{}},T/\[Theta]]/Gamma[\[Alpha],T/\[Theta]]+Log[T]
Lz[\[Alpha]_,\[Theta]_,z_List]:=ParallelSum[(\[Alpha]-1) Log[z[[i]]]-z[[i]]/\[Theta]-\[Alpha] Log[\[Theta]]-Log[Gamma[\[Alpha]]],{i,1,Length[z]}];
\[Mu]j[xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,xj,xj>=T,(\[Theta]t Gamma[\[Alpha]t+1,T/\[Theta]t])/Gamma[\[Alpha]t,T/\[Theta]t]];
\[Mu]jLn[xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,Log[xj],xj>=T,Elnz[\[Alpha]t,\[Theta]t,T]];

Q[\[Theta]_,\[Theta]t_,x_List,\[Alpha]i_,m_,T_]:= Lz[\[Alpha]i,\[Theta],x[[1;;m]]]+ParallelSum[-(\[Mu]j[x[[i+m]],\[Alpha]i,\[Theta]t,T]/\[Theta])+(\[Alpha]i-1) \[Mu]jLn[x[[i+m]],\[Alpha]i,\[Theta]t,T]-\[Alpha]i Log[\[Theta]]-Log[Gamma[\[Alpha]i]],{i,1,Length[x[[m+1;;All]]]}];

\[Mu]jN[xnj_,xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,xnj,xj>=T,(\[Theta]t Gamma[\[Alpha]t+1,T/\[Theta]t])/Gamma[\[Alpha]t,T/\[Theta]t]];
\[Mu]jLnN[xnj_,xj_,\[Alpha]t_,\[Theta]t_,T_]:=Which[xj<T,Log[xnj],xj>=T,Elnz[\[Alpha]t,\[Theta]t,T]];

QN[\[Theta]_,\[Theta]t_,xn_List,x0_List,\[Alpha]i_,m_,T_]:= Lz[\[Alpha]i,\[Theta],x0[[1;;m]]] + ParallelSum[-(\[Mu]jN[xn[[i+m]],x0[[i+m]],\[Alpha]i,\[Theta]t,T]/\[Theta])+(\[Alpha]i-1) \[Mu]jLnN[xn[[i+m]],x0[[i+m]],\[Alpha]i,\[Theta]t,T]-\[Alpha]i Log[\[Theta]]-Log[Gamma[\[Alpha]i]],{i,1,Length[x0[[m+1;;All]]]}];


(* ::Subsubsection:: *)
(*EM Routine: Classical Noiseless Version*)


(* ::Text:: *)
(*Fixed-Point implementation*)


(* Current: FixedPoint EM implementation *)
Options[GammaFullEM]:={Init:>Abs@RandomVariate@NormalDistribution[0,0.5]};

GammaFullEM[data_List,\[Alpha]_,T_,m_:0,tol_:6, maxIter_:25,opts___?OptionQ]:=Module[{\[Theta],init},
{init}={Init}/.{opts}/.Options[GammaFullEM];
FixedPointList[ 
	Check[
		ArgMax[{Q[\[Theta],#,data,\[Alpha],m,T],\[Theta]>0},\[Theta]], 
		(PrintTemporary["M-Step Error..."];Abort[]), 
		{NMaximize::cvmit, ArgMax::objv}
	]&,
	init,
	maxIter,
	SameTest->(Abs[#1-#2]<10^-tol&)
]
]



(* ::Subsubsection::Closed:: *)
(*Noise Auxillary Functions*)


(* Can't clip at Zero!!! *)
Options[NoiseModel]:={NoiseDistribution->(NormalDistribution[0,#]&)};

NoiseModel[data_List,\[Sigma]n_:0.1,opts___?OptionQ]:= Module[{ndist},
{ndist} = {NoiseDistribution}/.{opts}/.Options[NoiseModel];
Max[10000*$MachineEpsilon,#]&/@(data+RandomVariate[ndist[\[Sigma]n],Length@data])
]


(* ::Subsubsection:: *)
(*Noisy Gamma EM*)


Options[GammaFullNEM]:={CoolingSchedule->2, NoiseDistribution->(NormalDistribution[0,#]&), Init:>Abs@RandomVariate@NormalDistribution[0,0.5]};

GammaFullNEM[data_List,\[Alpha]_,T_,m_:0,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=Module[
{\[Theta],k=1,cs,ndist,init},
{cs,ndist,init}={CoolingSchedule,NoiseDistribution,Init}/.{opts}/.Options[GammaFullNEM];

If[\[Sigma]n!=0,
FixedPointList[
	Check[ 
		NArgMax[{QN[\[Theta],#,NoiseModel[data,\[Sigma]n/(k++)^cs,NoiseDistribution->ndist],data,\[Alpha],m,T],\[Theta]>0},\[Theta]], 
		(PrintTemporary["M-Step Error..."];Abort[]), 
		{NArgMax::cvmit,NMaximize::cvmit(*, ArgMax::objv*)} 
	]&,
	init,
	maxIter,
	SameTest->(Abs[#1-#2]<10^-tol&)
],
GammaFullEM[data,\[Alpha],T,m,tol, maxIter]
]
]



(* ::Subsubsection::Closed:: *)
(*Parallelization Prep*)


(* ::Text:: *)
(*The routines are annoyingly slow in sequential evaluation mode*)


DistributeDefinitions[Elnz,Lz,\[Mu]j,\[Mu]jLn,Q];
DistributeDefinitions[\[Mu]jN,\[Mu]jLnN,QN];
DistributeDefinitions[GammaFullEM,GammaFullNEM];
End[]; 


(* ::Subsection:: *)
(*Gaussian Mixture EM Model*)


Begin["`GmmEM`"];


(* ::Subsubsection::Closed:: *)
(*Auxillary Functions*)


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


(* ::Subsubsection:: *)
(*Noise Conditioning Functions: Gaussian Mixture*)


GaussNEMConditionQ[data_, noise_, means_List] := And@@Thread@(noise^2 <= 2*noise*(means-data));
OutsideMeanClusterQ[data_,means_List]:=Not/@IntervalMemberQ[Interval[{Min@means, Max@means}],data];

(*OutsideMeanClusterQv2[data_,m_List]:=Xnor@@Thread[(data<m)] <----This version is not Listable unlike the main one*)


indNoiseModel[data_List,dist_]:=Module[{n=RandomVariate[dist,Length@data]},data+n];


advNemCond[data_List,means_List,dist_]:=Module[
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
	n[[j]]=RandomVariate[TruncatedDistribution[(#[valints[[j]]])&/@{Min,Max},dist]],
	{j,nind}
];

data+n
];(* This version is a bit faster. But it uses a lot more memory in the interval struct *)


advNemCondv2[data_List,means_List,dist_]:=Module[
{modsamp=OutsideMeanClusterQ[data,means], n=0*Range@Length@data, nind, ps, ns},
nind = Select[Range@Length@data,modsamp[[#]]&];
ps=Select[nind,(data[[#]]<Min[means])&];
ns=Select[nind,(data[[#]]>Max[means])&];

Do[n[[j]]=RandomVariate[TruncatedDistribution[{Max[ 2(means-data[[j]]) ],0},dist]],
{j,ns}];
Do[n[[j]]=RandomVariate[TruncatedDistribution[{0,Min[ 2(means-data[[j]]) ]},dist]],
{j,ps}];

data+n
];


negNemCond[data_List,means_List,dist_]:=Module[
{valints,nind,n=0*Range@Length@data},
valints = Table[
	IntervalUnion[Interval[{Min[0,Max[ 2(data[[j]]-means) ]],0}],
		Interval[{0,Max[0,Min[ 2(data[[j]]-means) ]]}]
	],
	{j,Length@data}
];
(* Can only positively identify  [0,0]-intervals *)
nind=DeleteCases[Range@Length@data, j_/;valints[[j]]==0];
Do[ 
	n[[j]]=RandomVariate[TruncatedDistribution[(#[valints[[j]]])&/@{Min,Max},dist]],
	{j,nind}
];

data+n
];(* This version is a bit faster. But it uses a lot more memory in the interval struct *)


(* ::Subsubsection::Closed:: *)
(*Noisy GMM-EM*)


Options[GmmNEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

GmmNEM[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNEM];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#, indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#, indNoiseModel[data,ndist[#2]], data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#,advNemCond[data,m,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	(*Quiet@*)FixedPointList[(*obj*)
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection:: *)
(*Noise Harm (Neg) GMM-EM*)


Options[GmmNEMneg]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

GmmNEMneg[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNEMneg];
m = init\[Theta][[{2,4}]];
(*advNemCond[data_List,means_List,dist_]*)

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#, indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#, indNoiseModel[data,ndist[#2]], data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#,negNemCond[data,m,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];
(*Off[General::"unfl"];*)
If[\[Sigma]n!=0,	
	(*Quiet@*)FixedPointList[(*obj*)
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection::Closed:: *)
(*Noisy GMM-EM using Independent Noise (no NEM condition)*)


Options[GmmNEMInd]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

GmmNEMInd[data_List,\[Sigma]n_:0,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNEMInd];

obj=Switch[choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#,indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"means",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#,indNoiseModel[data,ndist[#2]], data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"sigmas",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];

Off[General::"unfl"];
If[\[Sigma]n!=0,	
	FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	GmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection::Closed:: *)
(*Noisy Expectation-Conditional-Maximization using modified NEM condition*)


Options[GmmNECM]:={Init:>{RandomReal[],0,1,3,1},CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

GmmNECM[data_List,\[Sigma]n_:0,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,init\[Theta],choose\[Theta],cs,ndist,k=1},
{init\[Theta],cs,ndist}={Init,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNECM];

obj[cur\[Theta]_,\[Sigma]_/;(\[Sigma]>0)]:=Module[{tmp\[Theta]=cur\[Theta],\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},
			(*Expectation Conditional-Maximization*)
			(* Maximize Q conditioned on current estimates of \[Sigma]1 and \[Sigma]2 *)
			tmp\[Theta][[1]]=NArgMax[{Q[{\[Theta]\[Alpha]}~Join~Flatten@tmp\[Theta][[2;;]],cur\[Theta],data],0<\[Theta]\[Alpha]<1},\[Theta]\[Alpha]];
			(*tmp\[Theta][[{2,4}]]=NArgMax[{QN[{tmp\[Theta][[1]],\[Theta]m1,tmp\[Theta][[3]],\[Theta]m2,tmp\[Theta][[5]]},tmp\[Theta],indNoiseModel[data,ndist@\[Sigma]],data]},{\[Theta]m1,\[Theta]m2}];*)
			tmp\[Theta][[{2,4}]]=NArgMax[{Q[{tmp\[Theta][[1]],\[Theta]m1,tmp\[Theta][[3]],\[Theta]m2,tmp\[Theta][[5]]},tmp\[Theta],data]},{\[Theta]m1,\[Theta]m2}];
			(* Maximize Q conditioned on updated estimates of \[Alpha], m1, and m2 *)
			tmp\[Theta][[{3,5}]]=NArgMax[{QN[{tmp\[Theta][[1]],tmp\[Theta][[2]],\[Theta]s1,tmp\[Theta][[4]],\[Theta]s2},tmp\[Theta],advNemCond[data,tmp\[Theta][[{2,4}]],ndist@\[Sigma]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}];
			tmp\[Theta]
		];

obj[cur\[Theta]_,\[Sigma]_/;(\[Sigma]<=0)]:=Module[{tmp\[Theta]=cur\[Theta],\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},
			(*Expectation Conditional-Maximization*)
			(* Maximize Q conditioned on current estimates of \[Sigma]1 and \[Sigma]2 *)
			(*tmp\[Theta]=(*N*)ArgMax[{Q[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},tmp\[Theta],d],0<\[Theta]\[Alpha]<1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}];*)
			tmp\[Theta][[1]]=NArgMax[{Q[{\[Theta]\[Alpha]}~Join~Flatten@tmp\[Theta][[2;;]],cur\[Theta],data],0<\[Theta]\[Alpha]<1},\[Theta]\[Alpha]];
			tmp\[Theta][[{2,4}]]=NArgMax[{Q[{tmp\[Theta][[1]],\[Theta]m1,tmp\[Theta][[3]],\[Theta]m2,tmp\[Theta][[5]]},tmp\[Theta],data]},{\[Theta]m1,\[Theta]m2}];
			(* Maximize Q conditioned on updated estimates of \[Alpha], m1, and m2 *)
			tmp\[Theta][[{3,5}]]=NArgMax[{Q[{tmp\[Theta][[1]],tmp\[Theta][[2]],\[Theta]s1,tmp\[Theta][[4]],\[Theta]s2},tmp\[Theta],data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}];
			tmp\[Theta]
		];

FixedPointList[
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,NArgMax::cvmit,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	]
(*,GmmEM[data,tol, maxIter,opts,Estimand->"popns"]]*)
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
(*Cauchy Mixture EM Model*)


Begin["`CmmEM`"];


(* ::Subsubsection::Closed:: *)
(*Auxillary Functions*)


(* ::Text:: *)
(*Helper functions for Cauchy mixture EM*)


(* Parameter Specification:  
\[Theta]={\[Alpha],\[Theta]_ 0,\[Theta]_ 1}
where \[Theta]_k={m_k,d_k}   
K is the number of different popns  *)

pzi[zi:(0|1), \[Alpha]_/;0<=\[Alpha]<=1]:=1-PDF[BernoulliDistribution[\[Alpha]],zi];
fyi[yi_,\[Theta]k_List/;VectorQ[\[Theta]k]]:=PDF[CauchyDistribution[\[Theta]k[[1]], \[Theta]k[[2]]],yi];
(* In fyzi[]: k+2 index locates the correct popn parameters for the k^th popn *)
fyzi[yi_,zi_,\[Theta]:{_,{_,_},{_,_}}]:=Exp@Sum[DiscreteDelta[zi-k]*(Log@pzi[k,\[Theta][[1]]]+Log@fyi[yi,\[Theta][[k+2]]]),{k,0,1}]; 
(* f_y^total (y|\[Theta])=(




Overscript[Subscript[\[Sum], k=0], K]\[Alpha]_k*f_k (y|\[Theta]_k)) *)
Fy[yi_,\[Theta]:{_,{_,_},{_,_}}]:=Sum[pzi[k,\[Theta][[1]]]*fyi[yi,\[Theta][[k+2]]],{k,0,1}];
pzjyi[zj:(0|1),yi_,\[Theta]:{_,{_,_},{_,_}}]:=fyzi[yi,zj,\[Theta]]/Fy[yi,\[Theta]]
Lyizi[\[Theta]:{_,{_,_},{_,_}},yi_,zi_]:=Sum[
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
QN[\[Theta]_List, \[Theta]t_List,y_List,y0_List] :=Sum[
Lyizi[{\[Theta][[1]],{\[Theta][[2]],\[Theta][[3]]},{\[Theta][[4]],\[Theta][[5]]}},y[[i]],j]*pzjyi[j,y0[[i]],{\[Theta]t[[1]],{\[Theta]t[[2]],\[Theta]t[[3]]},{\[Theta]t[[4]],\[Theta]t[[5]]}}],
{i,Length@y},{j,0,1}
];


(* ::Subsubsection::Closed:: *)
(*CMM-EM: Classical Noiseless Version*)


(* ::Text:: *)
(*Fixed-Point implementation*)


(* Current: FixedPoint EM implementation *)
PlaceParam[param_List,new\[Theta]_,loc_]:=ReplacePart[param,MapThread[Rule,{Flatten@{loc},Flatten@{new\[Theta]}}]];
(* Ensures that each objective function changes only parameters it is spec'd to change *)

Options[CmmEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all"};

CmmEM[data_List,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta]},
{init\[Theta],choose\[Theta]}={Init,Estimand}/.{opts}/.Options[CmmEM];
Off[General::"unfl"];
obj=Switch[(*OptionValue[Estimand]*)choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{Q[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#,data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"locs",PlaceParam[init\[Theta],NArgMax[{Q[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#,data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"disps",PlaceParam[init\[Theta],NArgMax[{Q[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#,data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{Q[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{Q[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];

(*Quiet@*)FixedPointList[(*obj*)
Check[ obj[#], (Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)]
]



(* ::Subsubsection:: *)
(*Noise Conditioning Functions: Cauchy Mixture*)


CauchyNEMConditionQ[data_, noise_, means_List] := And@@Thread@(noise^2 <= 2*noise*(means-data));
OutsideMeanClusterQ[data_,means_List]:=Not/@IntervalMemberQ[Interval[{Min@means, Max@means}],data];


indNoiseModel[data_List,dist_]:=Module[{n=RandomVariate[dist,Length@data]},data+n];


advNemCond[data_List,means_List,dist_]:=Module[
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
	n[[j]]=RandomVariate[TruncatedDistribution[(#[valints[[j]]])&/@{Min,Max},dist]],
	{j,nind}
];

data+n
];


(* ::Subsubsection::Closed:: *)
(*Noisy CMM-EM*)


Options[CmmNEM]:={Init:>{RandomReal[],0,1,3,1},Estimand->"all",CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)};

CmmNEM[data_List,\[Sigma]n_?NonNegative,tol_:6, maxIter_:25,opts___?OptionQ]:=
Module[{obj,\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2,init\[Theta],choose\[Theta],cs,ndist,k=1,m},
{init\[Theta],choose\[Theta],cs,ndist}={Init,Estimand,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[CmmNEM];
m = init\[Theta][[{2,4}]];
Off[General::"unfl"];
obj=Switch[(*OptionValue[Estimand]*)choose\[Theta],
"alpha",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha]}~Join~Flatten@init\[Theta][[2;;]],#, indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1},\[Theta]\[Alpha]],1]&,
"locs",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,init\[Theta][[3]],\[Theta]m2,init\[Theta][[5]]},#,indNoiseModel[data,ndist[#2]], data]},{\[Theta]m1,\[Theta]m2}],{2,4}]&,
"disps",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],init\[Theta][[2]],\[Theta]s1,init\[Theta][[4]],\[Theta]s2},#,advNemCond[data,m,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]s1,\[Theta]s2}],{3,5}]&,
"popns",PlaceParam[init\[Theta],NArgMax[{QN[{init\[Theta][[1]],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{2,3,4,5}]&,
"all",PlaceParam[init\[Theta],NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2},#,indNoiseModel[data,ndist[#2]], data],0<=\[Theta]\[Alpha]<=1&&\[Theta]s1>0&&\[Theta]s2>0},{\[Theta]\[Alpha],\[Theta]m1,\[Theta]s1,\[Theta]m2,\[Theta]s2}],{1,2,3,4,5}]&
];

Off[General::"unfl"];
If[\[Sigma]n!=0,	
	(*Quiet@*)FixedPointList[(*obj*)
		Check[ obj[#,\[Sigma]n/(k++)^cs], (PrintTemporary["M-Step Error..."];Abort[]), {NArgMax::nnum,Power::infy,\[Infinity]::indet} ]&,
		init\[Theta],maxIter,SameTest->(Max@Abs[#1-#2]<10^-tol&)
	],
	CmmEM[data,tol, maxIter,opts]
]
]


(* ::Subsubsection::Closed:: *)
(*Parallelization Prep*)


(* ::Text:: *)
(*The routines are annoyingly slow in sequential evaluation mode*)


DistributeDefinitions[pzi,fyi,fyzi,Fy,pzjyi,Lyizi,Qv1,Qv2,Q,QN,PlaceParam,
					CauchyNEMConditionQ,OutsideMeanClusterQ,
					advNemCond
]
DistributeDefinitions[CmmEM,CmmNEM]
End[];


(* ::Subsection::Closed:: *)
(*Gaussian N-D Mixture EM Model*)


Begin["`GmmNDEM`"];
(* Assume independent components for now. K is an option now... *)
(* CAVEAT: 1-D NormalDistribution specs std.dev, N-D MultinormalDistribution specs \[CapitalSigma] matrix i.e. \[Sigma]^2 *)


(* ::Subsubsection::Closed:: *)
(*Function Declarations*)


Options[GmmNDEM]:={
Init\[Alpha]:>ConstantArray[1/3,3], Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],
ClusterCount->3,ParameterLock->0,
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
LStoppingRule->False
};


Options[GmmNDECM]:={
Init\[Alpha]->ConstantArray[1/3,3], Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],
ClusterCount->3,ParameterLock->0,
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
LStoppingRule->False
};


Options[GmmNDEMLaw]:={
Init\[Alpha]->ConstantArray[1/3,3], Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],
CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
ClusterCount->3,ParameterLock->0,
LStoppingRule->False
};


Options[NEMClusterError]:={Init\[Alpha]->ConstantArray[1/3,3],
Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],
ClusterCount->3,CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
ParameterLock->4,MaxIteration->\[Infinity]
};


Options[NEMLClusterError]:={Init\[Alpha]->ConstantArray[1/3,3],
Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],
ClusterCount->3,CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
ParameterLock->4,MaxIteration->\[Infinity]
};


Options[EMClustering]:={ToleranceLevel->2, MaxIteration->25,Init\[Alpha]:>ConstantArray[1/3,3], 
Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],Estimand->"popn",
ClusterCount->3,CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&)
};


Options[NEMClusterErrorGTref]:={(*Init\[Alpha]:>ConstantArray[1/3,3],*)
Init\[Mu]->2{{-1,1},{-1,-1},{1,1}},Init\[Sigma]->ConstantArray[1,{3,2}],
ClusterCount->3,CoolingSchedule->2,NoiseDistribution->(NormalDistribution[0,#]&),
ParameterLock->0
};


(* ::Subsubsection::Closed:: *)
(*Auxillary Functions*)


(* Switching to one-based index instead of zero-based index *)
pzi[zi_Integer, w_List/;VectorQ[w]]:=PDF[EmpiricalDistribution[w->Range@Length[w]], zi];

fyi[yi_?VectorQ,\[Theta]k_List]:=PDF[MultinormalDistribution[\[Theta]k[[1]], DiagonalMatrix@\[Theta]k[[2]]],yi];

Fy[yi_?VectorQ, \[Theta]:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta]))]:=PDF[
		MixtureDistribution[\[Theta][[1]],Table[MultinormalDistribution[\[Theta][[2,k]], DiagonalMatrix@\[Theta][[3,k]]],{k,Length@\[Theta][[1]]}]
	],yi
]; (* This is the observed likelihood at \[Theta] for a single sample yi*)

(* f(y,z|\[Theta]): *)
fyzi[yi_?VectorQ, zi_Integer, \[Theta]:{{__},{__List},{__List}}] := Exp@Sum[
DiscreteDelta[zi-k]*(Log@pzi[k,\[Theta][[1]]]+Log@fyi[yi,{\[Theta][[2,k]],\[Theta][[3,k]]}]),
{k,Length[\[Theta][[1]]]}
]/;(Equal@@(Length/@\[Theta]))&&(1<=zi<=Length[\[Theta][[1]]]) ;

Ly[samples_List, \[Theta]:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta]))]:=Sum[N@Fy[yt, \[Theta]],{yt,samples}]


(* Implement Log-Sum-Exp Trick in case of excessive underflow *)

(* Log[f(y,z|\[Theta])]: Complete likelihood function. Nec. for Q-fxns. Not very useful for diagnostics *)
Lyizi[\[Theta]:{{__},{__List},{__List}}, yi_?VectorQ, zi_Integer]:=
Sum[
	(Log@pzi[k,\[Theta][[1]]]+Log@fyi[yi,{\[Theta][[2,k]],\[Theta][[3,k]]}])*DiscreteDelta[zi-k],
	{k,Length@\[Theta][[1]]}
]/;(Equal@@(Length/@\[Theta]))&&(1<=zi<=Length[\[Theta][[1]]]);

(* p(z|y,\[Theta]) *)
pzjyi[zj_Integer, yi_,\[Theta]:{{__},{__List},{__List}}]:=
fyzi[yi,zj,\[Theta]]/Fy[yi,\[Theta]]/;(Equal@@(Length/@\[Theta]))&&(1<=zj<=Length[\[Theta][[1]]]);
(*Exp[lfyzi[yi,zj,\[Theta]] - Log@Fy[yi,\[Theta]]]/;(Equal@@(Length/@\[Theta]))&&(1<=zj<=Length[\[Theta][[1]]]);*)


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
*)
Q[\[Theta]:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta])), 
\[Theta]t:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta]t)),
y_List] := Sum[
Lyizi[\[Theta],y[[i]],j]*pzjyi[j,y[[i]],\[Theta]t],
{i,Length@y},{j,Range@Length@\[Theta][[1]]}
];

(* Clean-Conditioning version of Q() *)
QN[\[Theta]:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta])), 
\[Theta]t:{{__},{__List},{__List}}/;(Equal@@(Length/@\[Theta]t)),
y_List,y0_List] := Sum[
Lyizi[\[Theta],y[[i]],j]*pzjyi[j,y0[[i]],\[Theta]t],
{i,Length@y},{j,Range@Length@\[Theta][[1]]}
];



Q[\[Theta]_?VectorQ, \[Theta]t_?VectorQ,y_List,K_Integer] := Q[shape[\[Theta],K,Length@y[[1]]],shape[\[Theta]t,K,Length@y[[1]]],y];

QN[\[Theta]_?VectorQ, \[Theta]t_?VectorQ,y_List,y0_List,K_Integer] := QN[shape[\[Theta],K,Length@y[[1]]],shape[\[Theta]t,K,Length@y[[1]]],y,y0];


shape[t_?VectorQ,K_Integer,Dim_Integer]:={t[[;;K]]/Total[t[[;;K]]]}~Join~Partition[Partition[t[[K+1;;]],Dim],K];

fillResults[results_List,inits__List, pl:{Repeated[(False|True),{3}]}, K_Integer,Dim_Integer]:=Module[
{fin=inits, res = Partition[results,K]},
	If[Not[pl[[1]]],fin[[1]] = res[[1]];
		res = res[[2;;]]
	];
	If[Not[pl[[2]]],fin[[2]] = First@Partition[Partition[Flatten[res[[;;Dim]]],Dim],K];
		res = res[[-Dim;;]]
	];
	If[Not[pl[[3]]], fin[[3]] = First@Partition[Partition[Flatten@res[[-Dim;;]],Dim],K] ];
	fin[[1]] = fin[[1]]/Total[fin[[1]]];
	fin
]


(* ::Subsubsection::Closed:: *)
(*Noise Conditioning Functions: Gaussian Mixture*)


advNemCond[data_List,means_List,dist_]:=Module[
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
	n[[j]]=RandomVariate[TruncatedDistribution[(#[valints[[j]]])&/@{Min,Max},dist]],
	{j,nind}
];

data+n
];(* This version is a bit faster. But it uses a lot more memory in the interval struct *)


(* *)
ndNemCond[data:{__List},means:{__List},dist_,\[Sigma]_/;\[Sigma]<=0]:=data;
ndNemCond[data:{__List},means:{__List},dist_,\[Sigma]_/;\[Sigma]>0]:=Module[{dim = Length@data[[1]],d=ConstantArray[0,Dimensions[data]]},
Table[d[[;;,col]]=advNemCond[data[[;;,col]],means[[;;,col]],dist@\[Sigma]],{col,dim}];
d
]


(* ::Subsubsection::Closed:: *)
(*N-D EM Routine: Merged Noiseless & Noisy*)


GmmNDEM[data_List,\[Sigma]n_?NonNegative,tol_:6,maxIter_:25,opts___?OptionQ]:=
Module[{obj,init\[Alpha],init\[Mu],init\[Sigma],K,dim=Length[data[[1]]],cs,ndist,k=1,plock,
\[Theta]\[Alpha],\[Theta]m,\[Theta]s,\[Alpha],m,s,upd\[Theta]s={},
recalc\[Sigma]=\[Sigma]n,srule
},

{init\[Alpha],init\[Mu],init\[Sigma],K,cs,ndist}={Init\[Alpha],Init\[Mu],Init\[Sigma],ClusterCount,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNDEM];
{plock}={ParameterLock}/.{opts}/.Options[GmmNDEM];
plock = (#==1)&/@IntegerDigits[plock,2,3]; (* Mask-style parameter lock {a,m,s} *)
If[And@@plock, Return[{init\[Alpha],init\[Mu],init\[Sigma]}]];  (* Nothing to update *)
If[Length[init\[Alpha]]<K, init\[Alpha] = ConstantArray[1,K]/K];

\[Theta]\[Alpha]=If[plock[[1]],init\[Alpha], AppendTo[upd\[Theta]s,Array[Subscript[\[Alpha], #]&,K]]; Array[Subscript[\[Alpha], #]&,K] ];
\[Theta]m=If[plock[[2]],init\[Mu], AppendTo[upd\[Theta]s,Array[(Subscript[m, ##])&,{K,dim}]]; Array[(Subscript[m, ##])&,{K,dim}]];
\[Theta]s=If[plock[[3]],init\[Sigma], AppendTo[upd\[Theta]s,Array[(Subscript[s, ##])&,{K,dim}]]; Array[(Subscript[s, ##])&,{K,dim}]];
upd\[Theta]s=Flatten@upd\[Theta]s;
If[$EMDebug, Print@upd\[Theta]s ];

obj[cur\[Theta]_,\[Sigma]_]:=Module[{tmp=cur\[Theta],smptmp},
		smptmp=data;
If[$EMDebug, Print[k,Grid[Transpose@cur\[Theta],Frame->All], "\n", Ly[data, cur\[Theta]] ] ];
		If[And[(k>=3),Not[plock[[2]]]], (* No \[Sigma]-noise on 1st iter. when \[Mu]-upd. is active. Since 1st iter is random init *)
			recalc\[Sigma]=\[Sigma]*((k-1)/(k-2))^cs;
			smptmp = ndNemCond[data,cur\[Theta][[2]],ndist,recalc\[Sigma]];
			If[$EMDebug,Print@{k,recalc\[Sigma]}]
		];
		If[plock[[2]], smptmp = ndNemCond[data,init\[Mu],ndist,\[Sigma]]];
		
		tmp=NArgMax[{QN[{\[Theta]\[Alpha],\[Theta]m,\[Theta]s},cur\[Theta],smptmp,data],
				And@@Thread[0.001<Flatten@{\[Theta]\[Alpha],\[Theta]s}] }, Flatten@upd\[Theta]s(*,MaxIterations->$MStepIterLim*)
		];
		(*shape[Flatten@tmp,K,dim]*)
		fillResults[Flatten@tmp,{init\[Alpha],init\[Mu],init\[Sigma]},plock,K,dim]
	];

srule = If[Not[(LStoppingRule/.{opts}/.Options[GmmNDECM])],
			((Max[Abs[Flatten@#1-Flatten@#2]]<10^-tol)&),
			((Abs[Ly[data,#1]-Ly[data,#2]]<10^-tol)&)
];

FixedPointList[
Check[ obj[#,\[Sigma]n/(k++)^cs], (Abort[]), {NArgMax::nnum,NMaximize::cvmit,NMaximize::nrnum,Power::infy,\[Infinity]::indet} ]&,
{init\[Alpha],init\[Mu],init\[Sigma]},
maxIter,
SameTest->srule
]
]


(* ::Subsubsection::Closed:: *)
(*N-D ECM Routine*)


GmmNDECM[data_List,\[Sigma]n_?NonNegative,tol_:6,maxIter_:25,opts___?OptionQ]:=
Module[{obj,init\[Alpha],init\[Mu],init\[Sigma],K,dim=Length[data[[1]]],cs,ndist,k=1,plock, recalc\[Sigma]=\[Sigma]n, srule},
{init\[Alpha],init\[Mu],init\[Sigma],K,cs,ndist}={Init\[Alpha],Init\[Mu],Init\[Sigma],ClusterCount,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNDECM];
{plock}={ParameterLock}/.{opts}/.Options[GmmNDECM];

plock = (#==1)&/@IntegerDigits[plock,2,3]; (* Mask-style parameter lock {a,m,s} *)
If[And@@plock, Return[{init\[Alpha],init\[Mu],init\[Sigma]}]];  (* Nothing to update *)
If[Length[init\[Alpha]]<K, init\[Alpha] = ConstantArray[1,K]/K];

obj[cur\[Theta]_,\[Sigma]_]:=Module[{tmp=cur\[Theta],\[Theta]\[Alpha],\[Theta]m,\[Theta]s,\[Alpha],m,s,smptmp},

If[$EMDebug, Print[k,Grid[Transpose@cur\[Theta],Frame->All], "\n", Ly[data, cur\[Theta]] ] ];
			\[Theta]\[Alpha]=Array[Subscript[\[Alpha], #]&,K];\[Theta]m=Array[(Subscript[m, ##])&,{K,dim}];\[Theta]s=Array[(Subscript[s, ##])&,{K,dim}];
			tmp[[1]]=If[plock[[1]],
						init\[Alpha],
						ArgMax[{Q[{\[Theta]\[Alpha],cur\[Theta][[2]],cur\[Theta][[3]]},cur\[Theta],data],And@@Thread[0.001<\[Theta]\[Alpha]]}, \[Theta]\[Alpha]]
					];
			tmp[[1]]=tmp[[1]]/Total[tmp[[1]]];

			tmp[[2]]=If[plock[[2]],
						init\[Mu],
						NArgMax[{Q[Flatten@{tmp[[1]],\[Theta]m,cur\[Theta][[3]]},Flatten@cur\[Theta],data,K]}, Flatten@\[Theta]m, MaxIterations->$MStepIterLim]
					];
			tmp[[3]]=If[plock[[3]],
						init\[Sigma],
						smptmp=data;
						If[And[(k>=3),Not[plock[[2]]]], (* No \[Sigma]-noise on 1st iter. when \[Mu]-upd. is active. Since 1st iter is random init *)
							recalc\[Sigma]=\[Sigma]*((k-1)/(k-2))^cs;
							smptmp = ndNemCond[data,Partition[tmp[[2]],dim],ndist,recalc\[Sigma]];
							If[$EMDebug,Print@{k,recalc\[Sigma]}]
						];
						If[plock[[2]], smptmp = ndNemCond[data,init\[Mu],ndist,\[Sigma]]];

						NArgMax[{QN[shape[Flatten@{tmp[[1]],tmp[[2]],\[Theta]s},K,dim], cur\[Theta], smptmp, data],
								And@@Thread[0.001<Flatten@\[Theta]s]}, Flatten@\[Theta]s,MaxIterations->$MStepIterLim]
					];
			shape[Flatten@tmp,K,dim]
	];

srule = If[Not[(LStoppingRule/.{opts}/.Options[GmmNDECM])],
			((Max[Abs[Flatten@#1-Flatten@#2]]<10^-tol)&),
			((Abs[Ly[data,#1]-Ly[data,#2]]<10^-tol)&)
];

FixedPointList[
Check[ obj[#,\[Sigma]n/(k++)^cs], (Abort[]), {NArgMax::nnum,NMaximize::cvmit,NMaximize::nrnum,Power::infy,\[Infinity]::indet} ]&,
{init\[Alpha],init\[Mu],init\[Sigma]},
maxIter,
SameTest->srule
]
]


(* ::Subsubsection::Closed:: *)
(*N-D EM by "Learning Law"*)


GmmNDEMLaw[data_List,\[Sigma]n_?NonNegative,tol_:6,maxIter_:25,opts___?OptionQ]:=
Module[{obj,init\[Alpha],init\[Mu],init\[Sigma],K,dim=Length[data[[1]]],cs,ndist,k=1,plock, recalc\[Sigma]=\[Sigma]n,srule},
{init\[Alpha],init\[Mu],init\[Sigma],K,cs,ndist}={Init\[Alpha],Init\[Mu],Init\[Sigma],ClusterCount,CoolingSchedule,NoiseDistribution}/.{opts}/.Options[GmmNDEMLaw];
{plock}={ParameterLock}/.{opts}/.Options[GmmNDEMLaw];

plock = (#==1)&/@IntegerDigits[plock,2,3]; (* Mask-style parameter lock {a,m,s} *)
If[And@@plock, Return[{init\[Alpha],init\[Mu],init\[Sigma]}]];  (* Nothing to update *)
If[Length[init\[Alpha]]<K, init\[Alpha] = ConstantArray[1,K]/K];

obj[cur\[Theta]_,\[Sigma]_]:=Module[{tmp=cur\[Theta],\[Theta]\[Alpha],\[Theta]m,\[Theta]s,\[Alpha],m,s,pzjsums,smptmp},

If[$EMDebug,Print[k,Grid[Transpose@cur\[Theta],Frame->All] , "\nLikelihood=", Ly[data, cur\[Theta]] ] ];
			pzjsums=Total@Table[(pzjyi[#,sample,cur\[Theta]])&/@Range[Length@First@cur\[Theta]],{sample,data}];
			If[plock[[1]],
				tmp[[1]]=init\[Alpha],
				tmp[[1]]=(pzjsums/Length@data);
				tmp[[1]]=tmp[[1]]/Total@tmp[[1]];
			];

			tmp[[2]]=If[plock[[2]],
						init\[Mu],
						Total@Table[
							(sample*pzjyi[#,sample,cur\[Theta]])&/@Range[Length@First@cur\[Theta]],
						{sample,data}]/pzjsums
					];
			tmp[[3]]=If[plock[[3]],
						init\[Sigma],
						smptmp=data;
						If[And[(k>=3),Not[plock[[2]]]], (* No \[Sigma]-noise on 1st iter. when \[Mu]-upd. is active. Since 1st iter is random init *)
							recalc\[Sigma]=\[Sigma]*((k-1)/(k-2))^cs;
							smptmp = ndNemCond[data,Partition[tmp[[2]],dim],ndist,recalc\[Sigma]];
							If[$EMDebug,Print@{k,recalc\[Sigma]}]
						];
						If[plock[[2]], smptmp = ndNemCond[data,init\[Mu],ndist,\[Sigma]]];
						Flatten@(Diagonal/@(Total@Table[
							(Outer[Times,
									(sample-cur\[Theta][[2,#]]),
									(sample-cur\[Theta][[2,#]])]*pzjyi[#,sample,cur\[Theta]]
							)&/@Range[Length@First@cur\[Theta]],
						{sample,smptmp}]/pzjsums))
					];
			shape[Flatten@tmp,K,dim]
	];

If[$EMDebug,
	Print["Starting Law FP iters.\nCluster Count = ", K, "\nInitial Parameters: " Grid[Transpose@{init\[Alpha],init\[Mu],init\[Sigma]},Frame->All]]
];

srule = If[Not[(LStoppingRule/.{opts}/.Options[GmmNDEMLaw])],
			((Max[Abs[Flatten@#1-Flatten@#2]]<10^-tol)&),
			((Abs[Ly[data,#1]-Ly[data,#2]]<10^-tol)&)
];

Quiet@FixedPointList[
Check[ obj[#,\[Sigma]n/(k++)^cs], (Abort[]), {DiagonalMatrix::mindet, MultinormalDistribution::posdefprm,General::unfl,Power::infy,\[Infinity]::indet} ]&,
{init\[Alpha],init\[Mu],init\[Sigma]},
maxIter,
SameTest->srule
]
]


(* ::Subsubsection::Closed:: *)
(*Membership Classifier for EM-Clustering*)


EMClassifier[data__List,\[Theta]:{{__},{__List},{__List}}]:=Table[
	Last@Ordering[(pzjyi[#,sample,\[Theta]])&/@Range[Length@\[Theta][[1]]]],
	{sample,data}
]/;(Equal@@(Length/@\[Theta]));


(*Comparison in terms of ref_List (ground truth)!!!
Generate every possible translation into reference "space"
Then pick minimum-dissimilarity label translation
*)
MatchLabelsToRef[new_List,ref_List]:=Block[{fullset=Table[new/.MapThread[Rule,{Union[new],p}],
{p,Permutations@Union[ref]}]},
fullset[[First@Ordering[HammingDistance[ref,#]&/@fullset],All]]
]



NEMClusterError[data_List,\[Sigma]n_,compAt_,opts:OptionsPattern[]]:=Module[
{\[Theta]comp,\[Theta]em,refclasses,compclasses},

\[Theta]em=GmmNDECM[data, \[Sigma]n, 2, OptionValue[MaxIteration],
Init\[Alpha]->OptionValue[Init\[Alpha]], Init\[Mu]->OptionValue[Init\[Mu]],Init\[Sigma]->OptionValue[Init\[Sigma]],
ParameterLock->OptionValue[ParameterLock],
ClusterCount->OptionValue[ClusterCount]
];
\[Theta]comp = \[Theta]em[[Ceiling[compAt*Length@\[Theta]em]]];

refclasses=EMClassifier[data,Last@\[Theta]em];
compclasses=EMClassifier[data,\[Theta]comp];

If[$EMDebug,
Print[{Ceiling[compAt*Length@\[Theta]em], Length@\[Theta]em}];
Print[TableForm@\[Theta]em];
];

{
N[HammingDistance[refclasses,MatchLabelsToRef[compclasses,refclasses]]/Length@compclasses],
Flatten[{\[Theta]comp,Last@\[Theta]em}]
}
]



NEMLClusterError[data_List,\[Sigma]n_,compAt_,opts:OptionsPattern[]]:=Module[
{\[Theta]comp,\[Theta]em,refclasses,compclasses},

\[Theta]em=GmmNDEMLaw[data, \[Sigma]n, 2, OptionValue[MaxIteration],
Init\[Alpha]->OptionValue[Init\[Alpha]], Init\[Mu]->OptionValue[Init\[Mu]],Init\[Sigma]->OptionValue[Init\[Sigma]],
ParameterLock->OptionValue[ParameterLock],
ClusterCount->OptionValue[ClusterCount]
];
If[$EMDebug,Print["End EM-Law"]];
\[Theta]comp = \[Theta]em[[Ceiling[compAt*Length@\[Theta]em]]];

refclasses=EMClassifier[data,Last@\[Theta]em];
compclasses=EMClassifier[data,\[Theta]comp];
If[$EMDebug,Print[refclasses]];

If[$EMDebug,
Print[{Ceiling[compAt*Length@\[Theta]em], Length@\[Theta]em}];
Print[TableForm@\[Theta]em];
];

{
N[HammingDistance[refclasses,MatchLabelsToRef[compclasses,refclasses]]/Length@compclasses],
Flatten[{\[Theta]comp,Last@\[Theta]em,Length@\[Theta]em}]
}
]


EMClustering[data__List,opts___?OptionQ]:=Module[{clusters,labls,\[Theta]fin,dim=Length@data[[1]],tol,maxiter},
{tol,maxiter}={ToleranceLevel,MaxIteration}/.{opts}/.Options[EMClustering];

\[Theta]fin = Last@GmmNDEM[data,0,tol,maxiter,opts];
labls = EMClassifier[data,\[Theta]fin];
clusters = (Pick[data,labls,#])&/@Range[Length@First[\[Theta]fin]]
];


(* ::Subsubsection::Closed:: *)
(*Less Useful EM-Clustering fxns*)


FindSoftClusters[data__List,\[Theta]:{{__},{__List},{__List}}]:=Map[
	(Pick[data,EMClassifier[data,\[Theta]],#])&,
	Range[Length@First[\[Theta]]]
]/;(Equal@@(Length/@\[Theta]));

(*giveClass:=Compile[
{{sample,_Real,1},{\[Theta]c, _Real, 1},{numClasses, _Integer}},
Last@Ordering@Chop@((pzjyi[#,sample,shape[\[Theta]c,numClasses,Length@sample]])&/@Range[numClasses])
];

EMClassifierCompiled[data__List,\[Theta]:{{__},{__List},{__List}}]:=Block[
{num=Length@First@\[Theta],thc = Flatten@\[Theta]},
Table[
	giveClass[sample,thc,num],
	{sample,data}
]
]/;(Equal@@(Length/@\[Theta]));*)


NEMClusterErrorGTref[sdata_List,\[Sigma]n_,numIters_Integer,w_List,opts___?OptionQ]:=Module[
{data = sdata[[All,2]], ref = sdata[[All,1]],\[Theta]em,classes},
\[Theta]em=Last@GmmNDECM[
data,\[Sigma]n,2,numIters,
Init\[Alpha]->w,ParameterLock->4,Init\[Mu]->OptionValue[Init\[Mu]](*,
ClusterCount->OptionValue[ClusterCount],Init\[Mu]->OptionValue[Init\[Mu]],Init\[Sigma]->OptionValue[Init\[Sigma]],
CoolingSchedule->OptionValue[CoolingSchedule],NoiseDistribution->OptionValue[NoiseDistribution]*)
];
classes=EMClassifier[data,\[Theta]em];
{
N[HammingDistance[ref,MatchLabelsToRef[classes,ref]]/Length@classes],
KolmogorovSmirnovTest[data,MixtureDistribution[\[Theta]em[[1]],
Table[MultinormalDistribution[\[Theta]em[[2,k]],DiagonalMatrix@\[Theta]em[[3,k]]],
{k,Length@\[Theta]em[[1]]}]] ]
}
]


(* ::Subsubsection::Closed:: *)
(*End*)


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
(*Noise Conditioning Functions: Gaussian Mixture*)


advNemCond[data_List,means_List,dist_]:=Module[
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
	n[[j]]=RandomVariate[TruncatedDistribution[(#[valints[[j]]])&/@{Min,Max},dist]],
	{j,nind}
];

data+n
];(* This version is a bit faster. But it uses a lot more memory in the interval struct *)


(* *)
ndNemCond[data:{__List},means:{__List},dist_,\[Sigma]_/;\[Sigma]<=0]:=data;
ndNemCond[data:{__List},means:{__List},dist_,\[Sigma]_/;\[Sigma]>0]:=Module[{dim = Length@data[[1]],d=ConstantArray[0,Dimensions[data]]},
Table[d[[;;,col]]=advNemCond[data[[;;,col]],means[[;;,col]],dist@\[Sigma]],{col,dim}];
d
]


(* ::Subsubsection::Closed:: *)
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
						If[And[(k>=3),Not[plock[[2]]]], (* No \[Sigma]-noise on 1st iter. when \[Mu]-upd. is active. Since 1st iter is random init *)
							recalc\[Sigma]=\[Sigma]*((k-1)/(k-2))^cs;(*recalc\[Sigma]=\[Sigma];*)
							smptmp = ndNemCond[data,Partition[tmp[[2]],dim],ndist,recalc\[Sigma]];
							If[$EMDebug,Print["{EM Iter, Noise Iter} = ", {k,recalc\[Sigma]}] ]
						];
						If[plock[[2]], smptmp = ndNemCond[data,init\[Mu],ndist,\[Sigma]]];
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
