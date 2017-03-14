(* ::Package:: *)

BeginPackage["Bootstrap`"];

$BSLibVersion = "BS-v0.0.1";


BootstrapCI::usage = "Create Confidence Intervals for each aggregated data point";

BootstrapCIBands::usage = "Plot profile with Confidence Bands for all aggregated data points ";


PrimaryColor::usage="";
BandColor::usage="";
TrendColor::usage="";


Begin["`Private`"];


BootstrapCI[samples_List, bn_List:{200,53}]:=Module[{trimData,bdata},
trimData = (Sort/@samples)[[All,3;;-3]]; 
(* Sorts each row.
1 row per param value.
*)
(* Drop off 4 endpoints to deal with extreme outliers *)

bdata=Table[
Table[N@Mean@RandomChoice[trimData[[k]],bn[[2]] ],{bn[[1]]}],
{k,Length@trimData}
];

bdata = Transpose[Sort/@bdata]; (*Converts to column vector per param *)
bdata = bdata[[Ceiling[0.025*bn[[1]]];;Floor[0.975*bn[[1]]],;;]] 
(* ~ 95% confidence band *)
]


Options[BootstrapCIBands]:={PrimaryColor->Blue,BandColor->Green,TrendColor->Red};


BootstrapCIBands[bdata_,index_,opts___?OptionQ]:=Module[{trend,pcol,bcol,tcol},
{pcol,bcol,tcol}={PrimaryColor,BandColor,TrendColor}/.{opts}/.Options[BootstrapCIBands];

trend=Mean[bdata];
Show[
Table[
ListPlot[Transpose@{index,bdata[[k]]},
PlotStyle->{bcol,PointSize[Small]}],
{k,Length@bdata}
],(* Dots *)
ListPlot[Transpose@{index[[1;;Length@trend]],trend},
PlotStyle->{pcol,Thickness[0.004]},Joined->True,
InterpolationOrder->1],(* Trend Line *)
ListPlot[Table[{index[[k]],trend[[1]]},{k,Length@trend}],
PlotStyle->{tcol,Dashed},Joined->True],(* Noise-free Threshold *)
(*AxesOrigin->orig,*)
PlotRange->Automatic
]
]

(*Protect[BootstrapCI,BootstrapCIBands];*)


End[];
EndPackage[];
