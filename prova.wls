(*```Mathematica*)

#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*You can do a lot of things with ".wls" files*)

Plot[Sin[x],{x,0,10}]

(* ::Text:: *)
(*Extra information can be found reading the official Mathematica doc: http://reference.wolfram.com/lang
uage/tutorial/WolframLanguageScripts.html*)

(* ::Subchapter:: *)
(*The big advantage is that graphics and outputs are not saved:*)

Print[ReadString["~/Temp/demo.wls"]]

(* ::Subsection:: *)
(*So you can use git as usual!*)
