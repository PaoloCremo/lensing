(* ::Package:: *)

SetDirectory["/home/paolo/Desktop/waveform/case_h/SIS/amps"];


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

zL=0.5;
zS=1.;

DL=c/H0*Da[0.,zL]; (* m *)
DS=c/H0*Da[0.,zS]; (* m *)
DLS=c/H0*Da[zL,zS]; (* m *)


(* ::Section:: *)
(*Amplification factor*)


\[Lambda] = 1.;
y = 5.;
Print["\[Lambda] = "<>ToString[\[Lambda]]]
Print["y = "<>ToString[y]]
ystr = StringReplace[ToString[y], "."->""];
Print["amps_case_h_H74_y"<>ystr<>"_M5_10-8_zL05"]


M= 5*10^8*smtokg; (*mass lens rest (lens) frame*)

(*SIS*)
\[Sigma]=((M*c^2*G)/(4*Pi^2)*(DS/(DL*DLS)))^(1/4);
\[Theta]E = 4*Pi*(\[Sigma]/c)^2*DLS/DS;
\[Psi]hatx[x_]:=x; 
(*
(*point mass*)
\[Theta]E = Sqrt[4*G*M/c^2*DLS/(DS*DL)];
\[Psi]hatx[x_] := Log[x];
(*\[Psi]hatx[x_] := \[Lambda]*Log[x]-(1/2*(\[Lambda]-1)*(x)^2);*)

amp1[f_]:=(
w=DS*DL/(c*DLS)*\[Theta]E^2*(1+zL)*2\[Pi]*f;
Exp[(Pi*w/4)+(I*w/2*Log[w/2])]*Gamma[1-(I*w/2)]*Hypergeometric1F1[I*w/2,1,(I*w*y^2)/2]
);
*)

amp2[f_]:=Parallelize[(
w=DS*DL/(c*DLS)*\[Theta]E^2*(1+zL)*2\[Pi]*f;
-I*w*Exp[I*w*\[Lambda]^2*y^2/2]*NIntegrate[x*BesselJ[0.,w*\[Lambda]*y*x]*Exp[I*w*\[Lambda]*(x^2/2-\[Psi]hatx[x])],{x,0,Infinity},Method->{"LevinRule"},WorkingPrecision->10]
)];

fhx[x_]:=Abs[x];
amp3[f_]:=Parallelize[(
w=DS*DL/(c*DLS)*\[Theta]E^2*(1+zL)*2\[Pi]*f;
-I*w*Exp[I*w*\[Lambda]^2*y^2/2]*NIntegrate[x*BesselJ[0.,w*\[Lambda]*y*x]*Exp[I*w*\[Lambda]*(x^2/2-fhx[x])],{x,0,Infinity},Method->{"LevinRule"},WorkingPrecision->10]
)];


(*xx = Range[0.001, 0.1, 0.001];
a = amp2[xx];
b = amp3[xx];*)


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


\[Lambda]str = StringReplace[ToString[NumberForm[\[Lambda],{2,1}]], "."->""];
(*Export["amps_case_h_H74_y"<>ystr<>"_M10-9_zL05_L"<>\[Lambda]str<>"_new", amps, "Binary", "DataFormat"->"Complex64"];*)
Export["amps_case_h_H74_y"<>ystr<>"_M5_10-8_zL05", amps, "Binary", "DataFormat"->"Complex64"];


(* ::Code::Initialization::Bold:: *)
Date[]


(* ::Subsection:: *)
(*prove*)


(*dt = 1;
n = 2157032; (*complete waveform*)
prefs=Range[n/2]/(dt*n);
df = 1/(dt*n)//N
indfin = IntegerPart[10^-3/df];
prefs = Take[prefs, {1, indfin}];
Length[prefs]
Length[prefs]/4//N
prefs[[-1]]//N*)
