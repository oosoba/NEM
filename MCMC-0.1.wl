(* ::Package:: *)

BeginPackage["MCMC`"];

$MCMCLibVersion = "MCMC-v0.1";


MetropolisHastingsSampling::usage = "";
nDimMetropolisHastingsSampling::usage = "";

SimulatedAnnealing::usage = "";


SampleDimensions::usage="";
CoolingSchedule::usage="";
TrendColor::usage="";
InitTemp::usage="";
ProposalDistribution::usage="";
Iterations::usage="";


Begin["`Private`"];


(* ::Subsection:: *)
(*MCMC for Optimization*)


(* 
q: 1-dim conditional sampling distribution; 
Temp0: Initial temperature;
Niter: Number of iterations;
opts: Miscellaneous options list
*)
Options[SimulatedAnnealing]:={
SampleDimensions->2,
CoolingSchedule-> ((#^-0.5)&),
InitTemp->10,
ProposalDistribution-> (CauchyDistribution[#,1]&),
Iterations->1000
};


(* cost: objective function for optimization; *)
SimulatedAnnealing[cost_,opts___?OptionQ]:=Module[
{\[Alpha]=1,n=1,
x0,xs,xprev,xcand,T,
ndims,cs,Temp0,q,Niter},

{ndims, cs, Temp0,q, Niter}={SampleDimensions,CoolingSchedule, InitTemp,ProposalDistribution,Iterations}/.{opts}/.Options[SimulatedAnnealing];

x0=RandomReal[CauchyDistribution[0,5],ndims];
xprev=x0;
xs = Partition[
Flatten@Last@Reap@While[
n<=Niter,
xcand =Table[RandomReal@q[xprev[[k]]],{k,ndims}];
(* q(.|xprev) *)
T=Temp0*cs[n];
\[Alpha] =Exp[(-cost[xcand]+cost[xprev])/T] ;
If[(RandomReal[]<=\[Alpha]), 
{Sow[{xcand}], xprev=xcand},
 n--
];
n++;],
ndims
];
Return@xs
]


(* ::Subsection:: *)
(*MCMC for Sampling*)


(* 
ProposalDistribution: candidate proposal distribution object;
Iterations: # of samples from p to generate
*)
Options[MetropolisHastingsSampling]:={
ProposalDistribution->((CauchyDistribution[#,1])&),
Iterations->1000
};


(*p: pdf of target distribution;*)
MetropolisHastingsSampling[p_, opts___?OptionQ]:=Module[
{t=0,\[Alpha]=1, n=1, 
xold=RandomVariate@CauchyDistribution[],xcand,
Niter,q
},
{q,Niter}={ProposalDistribution, Iterations}/.{opts}/.Options[MetropolisHastingsSampling];

While[
p[xold]<=0,
xold = RandomVariate@CauchyDistribution[]
];
Return@Flatten@Last@Reap@While[n<=Niter,
xcand = RandomReal[q[xold]];(* q(.|xold) *)
\[Alpha] = Min[(p[xcand]/p[xold]*PDF[q[xcand],xold]/PDF[q[xold],xcand]),1];
If[(RandomReal[]<=\[Alpha]), 
{Sow[xcand], xold=xcand},
 n--
];
n++;
]
]


(* 
ProposalDistribution: candidate proposal distribution object;
Iterations: # of samples from p to generate
*)
Options[nDimMetropolisHastingsSampling]:={
SampleDimensions->2,
ProposalDistribution->((CauchyDistribution[#,1])&),
Iterations->1000
};


(*p: n-D pdf of target distribution;*)
nDimMetropolisHastingsSampling[p_, opts___?OptionQ]:=Module[
{t=0,\[Alpha]=1, n=0, 
xold,xcand,
Niter,q,ndims
},
{ndims, q,Niter}={SampleDimensions,ProposalDistribution, Iterations}/.{opts}/.Options[nDimMetropolisHastingsSampling];


xold=RandomVariate[CauchyDistribution[],ndims];
While[
p[xold]<=0,
xold = RandomVariate[CauchyDistribution[], ndims]
];
Return@Partition[Flatten@Last@Reap@While[n<Niter,
xcand = Table[RandomReal@q[xold[[k]]],{k,ndims}];(* q(.|xold) *)
\[Alpha] = Min[(p[xcand]/p[xold]*Product[PDF[q[xcand[[k]]],xold[[k]]],{k,ndims}]/Product[PDF[q[xold[[k]]],xcand[[k]]],{k,ndims}]),1];
If[(RandomReal[]<=\[Alpha]), 
{Sow[xcand], xold=xcand},
 n--
];
n++;
],ndims]
];


(* ::Subsection:: *)
(*Miscellaneous Hidden Functions*)


CompareHistogram3D[dist_,{xmin_,xmax_,dx_},{ymin_,ymax_,dy_},df_:PDF,opts___?OptionQ]:=
Module[{hist,dplot},
hist=Histogram3D[RandomVariate[dist,10^4],{{xmin,xmax,dx},{ymin,ymax,dy}},ToString[df],PlotRange->{{xmin,xmax},{ymin,ymax},All},ChartBaseStyle->Opacity[0.5]];
dplot=Plot3D[df[dist,{x,y}],{x,xmin,xmax},{y,ymin,ymax},Mesh->{Range[xmin,xmax,dx],Range[ymin,ymax,dy]},PlotRange->All,MeshStyle->Gray,PlotStyle->Hue[.15,.7,.8]];
Show[hist,dplot]
]


(* ::Subsection:: *)
(*Optimization Test Functions*)


(* Rastrigin: Subscript[x, i]^opt=0 *)
rastrigin[x_]:= 10*Length@x+Total[(#^2-10Cos[2\[Pi]*#])&/@x]


(* Schwefel: Subscript[x, i]^opt=420.9687 *)
schwefel[x_]:= -418.9829*Length@x - Total[(#*Sin[Sqrt[Abs[#]]])&/@x]


(* Holder Table: x^opt=(+/-8.05502, +/-9.66459) *)
holder[x_]:= -Abs[Sin[x[[1]]]*Cos[x[[2]]]*Exp[Abs[1-Sqrt[x[[1]]^2+x[[2]]^2]/\[Pi]]]]

(* Levy #13 : x^opt=(1,1) *)
levy13[x_]:= (Sin[3\[Pi]*x[[1]]])^2+(x[[1]]-1)^2 (1+(Sin[2\[Pi]*x[[2]]])^2) +(x[[2]]-1)^2 (1+(Sin[2\[Pi]*x[[2]]])^2)

(* damavandi: Subscript[x, i]^opt=2 *)
damavandi[x_]:= (1-Abs[(Sin[\[Pi](x[[1]]-2)]*Sin[\[Pi](x[[2]]-2)])/(\[Pi]^2 (x[[1]]-2)(x[[2]]-2))])*(2+(x[[1]]-7)^2+2(x[[2]]-7)^2)


(* griewank: Subscript[x, i]^opt=0 *)
griewank[x_]:= Total[x^2]/4000-Times@@(Cos[x]/Sqrt[Range@Length@x])+1


(*
Plot3D[rastrigin[{x,y}], {x,-5,5},{y,-5,5}]

Plot3D[schwefel[{x,y}], {x,-500,500},{y,-500,500}]

Plot3D[holder[{x,y}], 
{x,-10,10},{y,-10,10},
PlotRange\[Rule]Full
]

Plot3D[levy13[{x,y}], 
{x,-10,10},{y,-10,10}
]

Plot3D[damavandi[{x,y}], 
{x,0,14},{y,0,14},
PerformanceGoal\[Rule]"Quality"
]

GraphicsRow[{Plot3D[griewank[{x,y}], 
{x,-600,600},{y,-600,600}
],Plot3D[griewank[{x,y}], 
{x,-5,5},{y,-5,5},PlotRange\[Rule]Full,PerformanceGoal\[Rule]"Quality"
]}]
*)


(* ::Subsection::Closed:: *)
(*Close Package*)


End[];
EndPackage[];
