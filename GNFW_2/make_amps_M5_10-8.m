(* ::Package:: *)

SetDirectory["/home/paolo/Desktop/waveform/case_h/GNFW_2/amps"];


Date[]


(* ::Section:: *)
(*constants*)


c=2.99792459*10^8 ;(*m/s*)
G=6.67408*10^(-11) ;(*m^3 kg^-1 s^-2*)
smtokg=1.98847*10^30 ;

(*Cosmological Parameters *)
convpc2m=3.08567758130573*10^(16);
convkpc2m=convpc2m*10^3;
convMpc2m=convpc2m*10^6;

H0=74.*10^3/convMpc2m;
om=0.3061;
w0=-1.;
w1=0.;

Da[zin_,zout_]:=1/(1+zout)*NIntegrate[1/Sqrt[om*(1+z)^3+(1-om)*(1+z)^(-3(1+w0+w1))*Exp[-3*w1*(z/(1+z))]],{z,zin,zout}];

zL=0.5;
zS=1.;

DL=c/H0*Da[0.,zL]; (* m *)
DS=c/H0*Da[0.,zS]; (* m *)
DLS=c/H0*Da[zL,zS]; (* m *)


(* ::Section:: *)
(*lens parameter*)


ML = 5*10^8*smtokg; (*mass lens rest (lens) frame*)

\[Rho]s = 4.938784307737533*10^-22;
rs = 1.271299163497*10^20; (*M=10^9*)
(*rs = 1.2712991634979609*10^20; (*M=5*10^8 - CHANGE ALSO ML*) *)

\[Theta]E = Sqrt[4*G*ML/c^2*DLS/(DS*DL)];
\[Alpha][R_, rhos_, rs_]:=(2 rhos rs (\[Pi] R-2 (Sqrt[-1+R^2] ArcTan[Sqrt[-1+R^2]]+Log[2/R])))/R;
M2D[rs_]:=\[Pi]*Re[\[Theta]E*DL/rs*rs^2*\[Alpha][\[Theta]E*DL/rs, \[Rho]s, rs]];


(* ::Section:: *)
(*Amplification factor*)


\[Lambda] = 1.;
y = 5.;
Print["\[Lambda] = "<>ToString[\[Lambda]]]
Print["y = "<>ToString[y]]
ystr = StringReplace[ToString[y], "."->""];
Print["amps_case_h_H74_y"<>ystr<>"_M10-9_zL05"]


(*GNFW 2 *)
g[R_]:=1/8 (-\[Pi]^2+4 \[Pi] R+4 (-2+Log[2]) Log[2]-4 Log[I-Sqrt[-1+R^2]]^2+Log[R] (8-4 Log[4]+8 Log[I-Sqrt[-1+R^2]])-4 I ArcTan[Sqrt[-1+R^2]] (2 Log[R]-Log[I-Sqrt[-1+R^2]]-Log[I+Sqrt[-1+R^2]])+4 I Sqrt[-1+R^2] (Log[I-Sqrt[-1+R^2]]-Log[I+Sqrt[-1+R^2]]));
\[Psi]hat1[x_]:=16*\[Pi]*G/c^2*\[Rho]s*rs*DLS*DL/DS*g[x];
\[Psi]hatx[x_]:=Re[\[Psi]hat1[x]]; (* Imaginary term is much smaller *)

amp2[f_]:=Parallelize[(
(*w=DS*DL/(c*DLS)*\[Theta]E^2*(1+zL)*2\[Pi]*f;*)
w=(1+zL)/c*DS*rs^2/(DL*DLS)*2\[Pi]*f;
-I*w*Exp[I*w*\[Lambda]^2*y^2/2]*NIntegrate[x*BesselJ[0.,w*\[Lambda]*y*x]*Exp[I*w*\[Lambda]*(x^2/2-\[Psi]hatx[x])],{x,10^-9,Infinity},Method->{"LevinRule"},WorkingPrecision->10]
)];


(* ::Section:: *)
(*calcolo AF*)


dt = 1;
n = 3157032; (*complete waveform*)
(*prefs=Range[(n-1)/2]/(dt*n);*) (*n dispari*)
prefs=Range[n/2]/(dt*n); (*n pari*)
df = 1/(dt*n)//N;
indfin = IntegerPart[4*10^-4/df];
prefs = Take[prefs, {1, indfin}];
(*
len = 220201;
ind=1; 
fs=Table[prefs[[i]],{i,(ind-1)*len+1,ind*len}];
*)
Length[prefs]
amps =amp2[prefs]; (*Exp[-I*Arg[amp2[10.]]]*)
Length[amps]


\[Lambda]str = StringReplace[ToString[\[Lambda]], "."->""];
(*Export["amps_case_h_H74_y"<>ystr<>"_M10-9_zL05_L"<>\[Lambda]str<>"_new", amps, "Binary", "DataFormat"->"Complex64"];*)
Export["amps_case_h_H74_y"<>ystr<>"_M5_10-8_zL05", amps, "Binary", "DataFormat"->"Complex64"];


(* ::Code::Initialization::Bold:: *)
Date[]


(* ::Subsection:: *)
(*prove*)


(*dt = 1;
n = 3157032; (*complete waveform*)
prefs=Range[n/2]/(dt*n);
df = 1/(dt*n)//N
indfin = IntegerPart[4*10^-4/df];
prefs = Take[prefs, {1, indfin}];
Print["Number of points : "<>ToString[Length[prefs]]]
Length[prefs]/4//N
Print["Final frequency : "<>ToString[prefs[[-1]]//N]]*)
