(* ::Package:: *)

SetDirectory["/home/paolo/Desktop/waveform/case_h/NFW_2/amps"];


Date[]


(* ::Section:: *)
(*constants*)


c=2.99792459*10^8 ;(*m/s*)
G=6.67408*10^(-11) ;(*m^3 kg^-1 s^-2*)
smtokg=1.98847*10^30 ;
(*m = 5*10^1*smtokg; (*mass source*)
m = 1*10^5*smtokg;
xx[f_]:= (Pi*G*m*f/c^3)^(2/3); (*this formula requires directly GW frequency *)
*)
(*Cosmological Parameters *)
convpc2m=3.08567758130573*10^(16);
(*convkpc2m=convpc2m*10^3;*)
convMpc2m=convpc2m*10^6;

(*H0=68.*10^3/convMpc2m;*)
H0=74.*10^3/convMpc2m;
om=0.3061;
w0=-1.;
w1=0.;

Da[zin_,zout_]:=1/(1+zout)*NIntegrate[1/Sqrt[om*(1+z)^3+(1-om)*(1+z)^(-3(1+w0+w1))*Exp[-3*w1*(z/(1+z))]],{z,zin,zout}];

zL=0.25;
zS=0.5;

DL=c/H0*Da[0.,zL]; (* m *)
DS=c/H0*Da[0.,zS]; (* m *)
DLS=c/H0*Da[zL,zS]; (* m *)


(* ::Section:: *)
(*NFW params*)


(*
(*M=10^11*)
\[Rho]s=6.56*10^-21;
rs=3.72*10^20;M*)

\[Rho]s = 7.954495462000312*10^-22;
rs = 3.4312734704119716*10^20; (*M=10^9*)
(*rs = 2.882022860939552*10^20; (*M=5*10^8*) *)

h[x_]:=-((2 Sqrt[-1+2/(1+x)] ArcTanh[Sqrt[-1+2/(1+x)]])/(-1+x))+Log[x/2];
M[\[Theta]_]:=Re[4*\[Pi]*\[Rho]s*rs^3*h[DL*\[Theta]/rs]]; (* theta in rad ; mass in kg *)(* 2\[Pi]*Integrate[x*\[CapitalSigma][x],x] *)
Mx[x_]:=Re[4*\[Pi]*\[Rho]s*rs^3*h[x]];


(* ::Section:: *)
(*Amplification factor*)


\[Lambda] = 1.;
y = 0.02804;
Print["\[Lambda] = "<>ToString[\[Lambda]]]
Print["y = "<>ToString[y]]
ystr = StringReplace[ToString[y], "."->""];
Print["amps_case_h_H74_y"<>ystr<>"_M10-9_zL05_FP"]


le = y-x+(16*\[Pi]*G/c^2*\[Rho]s*rs*DLS*DL/DS*h[x]/x); (*lens equation*)
xm1 = NSolve[le==0 && 0<x<1, x];
xm2 = NSolve[le==0 && 1<x<10^2, x];
If[Length[xm1]==1, xm=Values[xm1][[1,1]],xm=Values[xm2][[1,1]]];
(*xm*)


(*M= 10^2*smtokg/(1+zL);*)
M= 10^9*smtokg; (*mass lens rest (lens) frame*)

(*NFW*)
g[x_]:=1/8 (-(\[Pi]-2 I Log[2])^2+4 (2 Log[x]-Log[4 I-4 Sqrt[-1+x^2]]) (-I ArcTan[Sqrt[-1+x^2]]+Log[I-Sqrt[-1+x^2]])+4 I ArcTan[Sqrt[-1+x^2]] Log[I+Sqrt[-1+x^2]]);
gM[x_]:=1/2*Log[x/2]^2-(2*ArcTanh[Sqrt[(1-x)/(1+x)]]^2);

\[Psi]hat1[x_]:=16*\[Pi]*G/c^2*\[Rho]s*rs*DLS*DL/DS*g[x];
(*\[Psi]hat1[x_]:=g[x];*)
\[Psi]hatx[x_]:=Re[\[Psi]hat1[x]]; (* Imaginary term is much smaller *)
\[Phi]m = -(1/2*(xm-y)^2-\[Psi]hatx[xm]);

amp2[f_]:=Parallelize[(
(*w=DS*DL/(c*DLS)*\[Theta]E^2*(1+zL)*2\[Pi]*f;*)
w=(1+zL)/c*DS*rs^2/(DL*DLS)*2\[Pi]*f;
-I*w*Exp[I*w*\[Lambda]^2*y^2/2]*NIntegrate[x*BesselJ[0.,w*\[Lambda]*y*x]*Exp[I*w*\[Lambda]*(x^2/2-\[Psi]hatx[x]+\[Phi]m)],{x,0,Infinity},Method->{"LevinRule"},WorkingPrecision->10]
)];


(*
ff = Range[10^-5, 2*10^-4, 2*10^-6];
Length[ff]
(*ww=(1+zL)/c*DS*rs^2/(DL*DLS)*2\[Pi]*ff;*)
af = amp2[ff];
(*Re[-I*Log[af/Abs[af]]];
LogLinearPlot[Re[-I*Log[af/Abs[af]]]]*)
ListLogLinearPlot[Transpose@{ff, Re[-I*Log[af/Abs[af]]]}, Joined\[Rule]True, PlotRange\[Rule]All]
*)


n = 513; (*complete waveform in f domain*)
df = 10^-6;
prefs=Range[0, df*(n-1), df]; 

Length[prefs]
amps =amp2[prefs]; (*Exp[-I*Arg[amp2[10.]]]*)
Length[amps]


\[Lambda]str = StringReplace[ToString[\[Lambda]], "."->""];
(*Export["amps_case_h_H74_y"<>ystr<>"_M10-9_zL05_L"<>\[Lambda]str<>"_new", amps, "Binary", "DataFormat"->"Complex64"];*)
Export["amps_case_h_H74_y"<>ystr<>"_M10-9_zS05_zL025_FP", amps, "Binary", "DataFormat"->"Complex64"];


(* ::Code::Initialization::Bold:: *)
Date[]


(* ::Subsection:: *)
(*prove*)


(*
n = 513; (*complete waveform in f domain*)
df = 10^-6;
prefs=Range[0, df*(n-1), df]; 
*)
