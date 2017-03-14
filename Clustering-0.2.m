(* ::Package:: *)

(* ::Section::Closed:: *)
(*Preamble*)


BeginPackage["Clustering`"];

$ClusteringLibVersion = "Clustering-v0.2";


(* Free-Mode Clustering Tools/Options *)
AssignCluster::usage = "AssignCluster[data_List,\[Mu]i_List] assigns data samples to clusters based on proximity to supplied centroids \[Mu]i.";

(* Available Options *)
InitCentroids::usage = "Specifies initial cluster centroids. Use empty list to indicate a random data sample initialization.";

NoiseDistribution::usage = "Specify Noise distribution";
Quenching::usage = "";


(* K-means Clustering *)
KMeansClusters::usage = "K-Means Clustering Procedure. Takes a data list and number of clusters. Returns a k-means partition of the data list. Extra options are InitCentroids and Medoid.";

KMeansEvolution::usage = "K-Means Clustering Procedure. Takes a data list and number of clusters. Returns a k-means partition of the data and the ordered evolution list of k-means cluster centroids; the last entry gives the k-means optimal cluster centroids.  Extra options are InitCentroids and Medoid.";

KMeansNoisyEvolution::usage = "";

KMeansNEMNoisyEvolution::usage = "Uses GMM-NEM condition for noise injection";
(* Available Options *)
Medoid::usage = "Specfies measure of central tendency for summarizing test clusters. Default is the sample mean function(Mean).";


(* Competitive Learning Adaptive Clustering *)

UCLearning::usage = "Unsupervised Competitive Learning procedure. Takes a data generator and number of clusters. Returns the evolution of tuned cluster-centroids. Extra options are InitCentroids, LearningScale and NumLearningEpochs";

DCLearning::usage = "Differential Competitive Learning procedure. Takes a data generator and number of clusters. Returns the evolution of tuned cluster-centroids. Extra options are InitCentroids, LearningScale and NumLearningEpochs";

SCLearning::usage = "Supervised Competitive Learning procedure. Takes a data generator and number of clusters. Returns the evolution of tuned cluster-centroids. Extra options are InitCentroids, LearningScale and NumLearningEpochs";

SCLearningFreeRun::usage = "Same as SCL except don't stop until converged to within tolerance stationarity.";
(* Available Options *)
LearningScale::usage = "Initial value for slowly-decaying learning coefficients.";
NumLearningEpochs::usage = "# of Learning Iterations.";


(* ::Section:: *)
(*Implementation*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Auxilliary Functions*)


AssignCluster[data_List,\[Mu]i_List(*,opts___?OptionQ*)]:=Module[{clusters = ConstantArray[{},Length@\[Mu]i]},
Do[
AppendTo[clusters[[First[Ordering[Norm/@((data[[j]]-#)&/@\[Mu]i)]]]],data[[j]]]
,{j,Length@data}];
clusters
];


(* ::Subsection:: *)
(*K-Means Clustering*)


Begin["`Kmeans`"];


(* ::Subsubsection:: *)
(*Regular K-Means*)


Options[KMeansClusters]:={InitCentroids->{}, Medoid->Mean};

KMeansClusters[data_List, K_Integer,opts___?OptionQ]:=Module[{\[Mu]i, icent, med},(* K <<n *)
{icent, med}={InitCentroids, Medoid}/.{opts}/.Options[KMeansClusters];
If[Length@icent==K, \[Mu]i=icent, \[Mu]i=RandomSample[data,K]];

AssignCluster[data, FixedPoint[(med/@AssignCluster[data, #])&, \[Mu]i]]
];


Options[KMeansEvolution]:={InitCentroids->{}, Medoid->Mean};

KMeansEvolution[data_List, K_Integer,opts___?OptionQ]:=Module[{\[Mu]i, evolve, icent, med},(* K <<n *)
{icent, med}={InitCentroids, Medoid}/.{opts}/.Options[KMeansEvolution];
If[Length@icent==K, \[Mu]i=icent, \[Mu]i=RandomSample[data,K] ];

evolve=FixedPointList[(med/@AssignCluster[data, #])&, \[Mu]i];
{evolve[[;;-2]],Last@evolve,AssignCluster[data, Last@evolve]}
];


(* ::Subsubsection::Closed:: *)
(*K-Means + Blind Noise*)


Options[KMeansNoisyEvolution]:={InitCentroids->{}, Medoid->Mean, 
NoiseDistribution->(NormalDistribution[0,#]&),Quenching->2};

KMeansNoisyEvolution[data_List,K_Integer,\[Sigma]n_:0,opts:OptionsPattern[]]:=Module[
{\[Mu]i, evolve, icent, med, ndist,t=1,cs},(* K <<n *)
{icent,med,ndist,cs}={InitCentroids,Medoid,NoiseDistribution,Quenching}/.{opts}/.Options[KMeansNoisyEvolution];
If[Length@icent==K, \[Mu]i=icent, \[Mu]i=RandomSample[data,K] ];

evolve=FixedPointList[
(
med/@AssignCluster[
data+If[\[Sigma]n!=0,RandomVariate[ndist[\[Sigma]n/(t++)^cs], Dimensions@data],0], 
#
])&,
\[Mu]i,
SameTest->((Max[Abs[Flatten@#1-Flatten@#2]]<10^-2)&)];
{evolve[[;;-2]],Last@evolve,AssignCluster[data, Last@evolve]}
];



(* ::Subsubsection:: *)
(*K-Means + NEM Condition*)


Options[KMeansNEMNoisyEvolution]:={InitCentroids->{}, Medoid->Mean, 
NoiseDistribution->(NormalDistribution[0,#]&),Quenching->2};

KMeansNEMNoisyEvolution[data_List,K_Integer,\[Sigma]n_:0,opts:OptionsPattern[]]:=Module[
{\[Mu]i, evolve, icent, med, ndist,t=1,cs},(* K <<n *)
{icent,med,ndist,cs}={InitCentroids,Medoid,NoiseDistribution,Quenching}/.{opts}/.Options[KMeansNEMNoisyEvolution];
If[Length@icent==K, \[Mu]i=icent, \[Mu]i=RandomSample[data,K] ];
\[Mu]i = med/@AssignCluster[data, \[Mu]i]; (* Do 1st cluster assignment & get 1st good set of centroids *)

evolve=FixedPointList[
	(med/@AssignCluster[ndNemCond[data,#,ndist,(\[Sigma]n/(t++)^cs)],#])&,
	\[Mu]i,
	SameTest->((Max[Abs[Flatten@#1-Flatten@#2]]<10^-2)&)
];

{evolve,Last@evolve,AssignCluster[data, Last@evolve]}
];



(* *)
ndNemCond[data:{__List},means:{__List},dist_,\[Sigma]_/;\[Sigma]<=0]:=data;
ndNemCond[data:{__List},means:{__List},dist_,\[Sigma]_/;\[Sigma]>0]:=Module[{dim = Length@data[[1]],d=ConstantArray[0,Dimensions[data]]},
Table[d[[;;,col]]=advNemCond[data[[;;,col]],means[[;;,col]],dist@\[Sigma]],{col,dim}];
d
]


advNemCond[data_List,means_List,dist_]:=Module[
{valints,nind,n=0*Range@Length@data},
valints = Table[
	IntervalUnion[Interval[{Min[0,Max[means-data[[j]]]],0}],
		Interval[{0,Max[0,Min[means-data[[j]]]]}]
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
];(* Using N-Dim GMM data NEMifier *)


End[];


(* ::Subsection::Closed:: *)
(*Competitive Learning*)


Begin["`CL`"];


Options[UCLearning]:={InitCentroids->{}, LearningScale->0.3, NumLearningEpochs->1000};

UCLearning[datagen_,p_Integer,opts___?OptionQ]:=Module[{xt,c,m,stop,winner},
{m,c,stop}={InitCentroids, LearningScale, NumLearningEpochs}/.{opts}/.Options[UCLearning];
(*m = Partition[Partition[m,2],p];*)
(*Print[m];*)
If[Length@m!=p, m =Table[datagen[],{p}] ];

First@Last@Reap[
Sow@m;
Do[
xt=datagen[];
winner = First@Ordering[Norm[(xt-#)]&/@m];
m[[winner]] =m[[winner]] +c*(1-t/stop)*(xt-m[[winner]]);
Sow@m,
{t,1,stop}];
]
]



Options[SCLearning]:={InitCentroids->{}, LearningScale->0.3, NumLearningEpochs->1000};

SCLearning[sdatagen_,p_Integer,opts___?OptionQ]:=Module[
{xt,rxt,c,m,stop,rm,winner,super},
{m,c,stop}={InitCentroids, LearningScale, NumLearningEpochs}/.{opts}/.Options[SCLearning];
If[Length@m!=p, m =sdatagen[InitMode->True] ];
rm = Range@p;

First@Last@Reap[
Sow@m;
Do[
{rxt,xt}=sdatagen[];
winner = First@Ordering[Norm[(xt-#)]&/@m];
super=Boole[rxt==rm[[winner]]]-Boole[rxt!=rm[[winner]]];
m[[winner]] =m[[winner]] +c*(1-t/stop)*super*(xt-m[[winner]]);
Sow@m,
{t,1,stop}];
]
]


SCLearningFreeRun[sdatagen_,p_Integer,tol_:2,opts___?OptionQ]:=Module[
{xt,rxt,c,m,stop,rm,winner,super,delta,t},
{m,c,stop}={InitCentroids, LearningScale, NumLearningEpochs}/.{opts}/.Options[SCLearning];
If[Length@m!=p, m =sdatagen[InitMode->True] ];
rm = Range@p;

First@Last@Reap[
Sow@m;
t=1;delta=10;
While[Not@And[t>stop, Abs[delta]<10^-tol],
{rxt,xt}=sdatagen[];
winner = First@Ordering[Norm[(xt-#)]&/@m];
super=Boole[rxt==rm[[winner]]]-Boole[rxt!=rm[[winner]]];
delta=(c/t)*super*(xt-m[[winner]]);
m[[winner]] =m[[winner]] + delta;
Sow@m; t++
];
]
]


Options[DCLearning]:={InitCentroids->{}, LearningScale->0.3, NumLearningEpochs->1000};

DCLearning[datagen_,p_Integer,opts___?OptionQ]:=Module[{xt,yt,c,m,stop,winner,\[CapitalDelta]y,W=3*IdentityMatrix[p]+ConstantArray[-1,{p,p}]},
{m,c,stop}={InitCentroids, LearningScale, NumLearningEpochs}/.{opts}/.Options[DCLearning];
If[Length@m!=p, m =Table[datagen[],{p}] ];

First@Last@Reap[
Sow@m;
Do[
xt = datagen[];
yt = xt.Transpose[m];
winner = First@Ordering[Norm[(xt-#)]&/@m];
\[CapitalDelta]y = (yt.W[[;;,winner]])+(xt.m[[winner]]);
m[[winner]] = m[[winner]] + c*Sign[\[CapitalDelta]y]*(1-t/stop)*(xt-m[[winner]]);
Sow@m,
{t,1,stop}];
]
]


End[];


(* ::Section::Closed:: *)
(*End*)


End[];
EndPackage[];
