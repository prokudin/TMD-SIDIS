(******************************************************************)
(* :Package:                                                      *)
(* :Title: NNPDF                                                  *)
(* :Version: 2.2                                                  *)
(* :Author: Emanuele R. Nocera (emanuele.nocera@edu.unige.it)     *)
(* :Summary: Tools and Tables                                     *)
(* :Terms of use: GPL                                             *)
(* :Mathematica version 8.0                                       *)
(* :History: Version 1.0 written by Emanuele R. Nocera May 2012   *)
(*           Version 1.1 revised Sep 2012 (added alphas function) *)     
(*           Version 2.0 modified March 2013 (updated the package *)
(*           for handling polarized NNPDF grids)                  *)
(*           Version 2.1 modified March 2014 (updated the package *)
(*           for handling NNPDF grids with photon and new NNPDF   *)            
(*           polarized grids)                                     *)
(*           Version 2.2 modified Sep 2014 (updated th package    *)
(*           for handling NNPDF3.0 parton sets)                   *)
(******************************************************************)

Print["Welcome in the Mathematica Package for NNPDF PDFs. The following Functions are available:"];
Print["- Initializing and interpolating PDFs"];
Print["  InitializePDFGrid;"];
Print["- Using NNPDF PDFs"];
Print["  xPDFcv; xPDFEnsemble; xPDFRep; xPDF; xPDFCL; xMin; xMax; Q2min; Q2Max; NumberPDF; HasPhoton;"];
Print["- Parameters for the QCD analysis"];
Print["  alphaSMZ; mCharm; mBottom; mTop; MZ; alphas; Infoalphas; Lam4; Lam5;"];
Print["For the usage of each of these Functions, please read the README file or type ?Function."];

BeginPackage["NNPDF`"];

InitializePDFGrid::usage="The function InitializePDFGrid[path, namegrid] reads the PDF grid into memory specified by namegrid from the location specified by path. It also performs the needed interpolation among the grid points.";

xPDFcv::usage="The function xPDFcv[x, QSq, f] returns x times the central value of the PDF with flavor f at a given momentum fraction x and scale QSq in Gev^2. Note that f must be an integer, and x and QSq must be numeric quantities. For the unpolarized case, the LHAPDF convention is used for the flavor f , that is, f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 corresponds to tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t; f=7 returns the photon PDF (if available). Only for NNPDFpol1.0, the following convention is used for the flavor f, that is, f = 0, 1, 2, 3, 4 corresponds to polarized g, u+ubar, d+dbar, s+sbar.";

xPDFEnsemble::usage="The function xPDFEnsemble[x, QSq, f] returns x times the vector of PDF replicas of flavor f at a given momentum fraction x and scale QSq in GeV^2. Note that f must be an integer, and x and QSq must be numeric quantities. For the unpolarized case, the LHAPDF convention is used for the flavor f, that is, f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 corresponds to tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t; f=7 returns the photon PDF (if available). Only for NNPDFpol1.0, the following convention is used for the flavor f, that is, f = 0, 1, 2, 3, 4 corresponds to polarized g, u+ubar, d+dbar, s+sbar.";

xPDFRep::usage="The function xPDFRep[x,QSq,f,irep] returns x times the irep PDF replica of flavor f at a given momentum fraction x and scale QSq in GeV^2. Note that f must be an integer, and x and QSq must be numeric quantities. For the unpolarized case, the LHAPDF convention is used for the flavor f, that is, f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 corresponds to tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t; f=7 returns the photon PDF (if available). Only for NNPDFpol1.0, the following convention is used for the flavor f, that is, f = 0, 1, 2, 3, 4 corresponds to polarized g, u+ubar, d+dbar, s+sbar.";

xPDF::usage="The function xPDF[x,QSq,f] returns x times the value of the PDF of flavor f at a given momentum fraction x and scale QSq in GeV^2 with its standard deviation. Note that f must be an integer, and x and QSq must be numeric quantities. For the unpolarized case, the LHAPDF convention is used for the flavor f, that is, f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 corresponds to tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t; f=7 returns the photon PDF (if available). Only for NNPDFpol1.0, the following convention is used for the flavor f, that is, f = 0, 1, 2, 3, 4 corresponds to polarized g, u+ubar, d+dbar, s+sbar.";

xPDFCL::usage="The function xPDFCL[ensemble,x,QSq,f,CL] returns x times the value of the PDF of flavor f at a given momentum fraction x and scale Q2 in GeV^2 with its error.  Note that f must be an integer, and x and QSq must be numeric quantities. For the unpolarized case, the LHAPDF convention is used for the flavor f, that is, f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 corresponds to tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t; f=7 returns the photon PDF (if available). Only for NNPDFpol1.0, the following convention is used for the flavor f, that is, f = 0, 1, 2, 3, 4 corresponds to polarized g, u+ubar, d+dbar, s+sbar. The variable ensemble is the PDF ensemble, as a function of x, QSq and f. It can be given by the function xPDFEnsemble (LHAPDF basis) or by any PDF ensemble (in a different basis) defined by the user. The variable CL denotes the confidence interval and must be a real number between 0 and 1.";

ForcepolPDFtoLHAPDF::usage="Polarized PDFs are given in the flavour basis, with quark antiquark flavors not disentangled. This function allows to separate the quark and antiquark contributions using dummy assumptions.";

UnForcepolPDFtoLHAPDF::usage="Polarized PDFs are correctly resetted to use only quark and antiquark combinations.";

HasPhoton::usage="Returns information about the availability of the photon PDF in the initialized grid.";

alphaSMZ::usage="Returns the List of alphaS values at the Z mass used in the QCD analysis for each replica.";
mCharm::usage="Charm quark mass in GeV.";
mBottom::usage="Bottom quark masss in GeV.";
mTop::usage="Top quark mass in GeV.";
MZ::usage="Z boson mass in GeV.";

alphas::usage="The function alphas[Q2,ipo,imodev] returns the QCD strong coupling constant alpha_s. The inputs are:
Q2: the energy scale (in GeV^2) at which alpha_s is computed;
ipo: the perturbative order at which alpha_s is computed (0=LO; 1=NLO; 2=NNLO);
imodev: the evolution mode with which alpha_s is computed (0: alpha_s is computed as a function of alpha_s at the Z mass reference scale as given in the .LHgrid file; 1: exact solution of the QCD beta function equation using Runge-Kutta algorithm)";

xMin::usage="The function xMin[] returns the minimum value of the x-grid.";
xMax::usage="The function xMax[] returns the maximum value of the x-grid.";

Q2Min::usage="The function Q2Min[] returns the minimum value of the Q2-grid.";
Q2Max::usage="The function Q2Max[] the maximum value of the Q2-grid.";

NumberPDF::usage="The function NumberPDF[] returns the number of PDF membes in the set.";

Lam4::usage="The function Lam4[] returns the lambdaQCD4 used in the QCD analysis.";
Lam5::usage="The function Lamm5[]returns the lambdaQCD5 used in the QCD analysis.";
Infoalphas::usage="The function Infoalphas[] returns info on alpha_s evolution.";

Begin["`Private`"];

InitializePDFGrid[path_String,namegrid_String]:=Module[{pdffile,pdfdat},
nPDF=13;
pdffile=path<>namegrid<>".LHgrid";
pdfdat=OpenRead[pdffile];

(* Initialize label of polarized set *) 
polflag=0;
photonflag=0;
lbl=Characters[namegrid];
If[lbl[[6]]=="p"&&lbl[[7]]=="o"&&lbl[[8]]=="l"&&lbl[[9]]=="1"&&lbl[[10]]=="0",polflag=1];

(* Read earliest version of LHAPDF*)
line=Read[pdfdat,String];

(* Read description of PDF set*)
While[line!=" 'Alphas:'",{Print[line],line=Read[pdfdat,String]}];

(* Read information on alphaS *)
alphainfo=ReadList[pdfdat,String,1];
idum=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!="{,}",space=ReadList[pdfdat,Character,1]];
MZ=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];
mCharm=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];
mBottom=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];
mTop=Read[pdfdat,Number];					
					 
(* Read information on interpolation grid  *)
While[line!=" 'MinMax:'",line=Read[pdfdat,String]];
NRep=Read[pdfdat,Number];
space=ReadList[pdfdat,Character,2];
idum=Read[pdfdat,Number];
xmin=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];
xmax=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];
Q2min=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];
Q2max=Read[pdfdat,Number];
			
(* Read QCD parameters *)
While[line!=" 'QCDparams:'",line=Read[pdfdat,String]];
NRep=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];					 
idum=Read[pdfdat,Number]; 
lambdaQCD4=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];	
lambdaQCD5=Read[pdfdat,Number];

(* Read alphaS parameters *)
While[line!=" 'Parameterlist:'",line=Read[pdfdat,String]];
space=Read[pdfdat,Word];
NRep=Read[pdfdat,Number];
space=ReadList[pdfdat,Character,2];
idum=Read[pdfdat,Number];
alphaSMZ=ReadList[pdfdat,Number,NRep+1];
			
(* Read the interpolation grid *)
While[line!=" 'Evolution:'",line=Read[pdfdat,String]];
space=Read[pdfdat,Word];
Q20=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!=",",space=Read[pdfdat,Character]];	
idum=Read[pdfdat,Number];
space=Read[pdfdat,Character];
While[space!="'",space=Read[pdfdat,Character]];	
interpolation=Read[pdfdat,Word];
If[interpolation=="NNPDF20intqed'",{photonflag=1;nPDF=14}];
line=Read[pdfdat,String];
				
npx=Read[pdfdat,Number];
xgrid=ReadList[pdfdat,Number,npx]; 
Which[Length[xgrid]!=npx,{Print["Mismatch in the number of points of the x grid."],Print[Length[xgrid],npx],Abort[]}];
Which[xgrid[[1]]!=xmin,{Print["Mismatch in the minimum extremum of the x grid."],Abort[]}];
Which[xgrid[[npx]]!=xmax,{Print["Mismatch in the maximum extremum of the x grid."],Abort[]}];

npQ2=Read[pdfdat,Number]; 
Q20=Read[pdfdat,Number];
Q2grid=ReadList[pdfdat,Number,npQ2];
Which[npQ2!=ptQ2,{Print["Mismatch in the number of points in the QSq grid."],Abort[]}];
Which[Length[Q2grid]!=npQ2,{Print["Mismatch in the number of points of the QSq grid."],Abort[]}];
(* Which[Q2grid[[1]]!=Q2min,{Print["Mismatch in the minimum extremum of the QSq grid."],Print[{Q2grid[[1]],Q2min}],Abort[]}]; *)
Which[Q2grid[[npQ2]]!=Q2max,{Print["Mismatch in the maximum extremum of the QSq grid."],Print[{npQ2,Q2max}],Abort[]}];

For[idx=0,idx<npQ2-1,idx++,If[Q2grid[[idx+1]]==Q2grid[[idx+2]],{Q2grid[[idx+1]]=Q2grid[[idx+1]]-0.000001;Q2grid[[idx+2]]=Q2grid[[idx+2]]+0.000001}]];

(* Output to screen some useful information *)
Print["Grid division: "];
Print["x axis: ", npx," points; ", ScientificForm[xmin,1], " <= x <= ", xmax];
Print["Q^2 axis: ", npQ2," points; ", Q2min, " <= Q^2 <= ", ScientificForm[Q2max,1]];
				 
(* Read PDFs *)	
Print["Reading PDF set from file. Please wait."];				 
NRep=Read[pdfdat,Number]; 
allPDF=ReadList[pdfdat,Number,nPDF*npQ2*npx*(NRep+1)];
space=Read[pdfdat,Character];
While[space!="'",space=Read[pdfdat,Character]];			 
end=Read[pdfdat,Word];
If[end=="End:'",{Close[pdfdat],Print["PDF set successfully read from file."]},{Print["an error occurred while processing LHPDF grid"],Abort[]}]; 
	 
(* Interpolating PDFs *)
Print["Interpolating PDF set on the grid. Please wait."];					 
xPDFfv=Partition[Partition[Partition[allPDF,nPDF],npQ2],npx];
Do[xPDFev[irep,ipdf]=ListInterpolation[Table[Table[xPDFfv[[irep]][[ix]][[iQ2]][[ipdf+7]],{iQ2,1,npQ2}],{ix,1,npx}],{xgrid,Q2grid},InterpolationOrder->2],{irep,1,NRep+1},{ipdf,-6,nPDF-7}];
Print["PDF grid interpolated. Ready to use the NNPDF set."];

];


xPDFcv[x_,Q2_,f_]:=If[polflag==0,Which[photonflag==0&&(f==0||f==1||f==2||f==3||f==4||f==5||f==6||f==-1||f==-2||f==-3||f==-4||f==-5||f==-6),xPDFev[1,f][x,Q2],photonflag==1&&(f==0||f==1||f==2||f==3||f==4||f==5||f==6||f==7||f==-1||f==-2||f==-3||f==-4||f==-5||f==-6),xPDFev[1,f][x,Q2],photonflag==0&&f!=0&&f!=1&&f!=2&&f!=3&&f!=4&&f!=5&&f!=6&&f!=-1&&f!=-2&&f!=-3&&f!=-4&&f!=-5&&f!=-6,Print["Incorrect flavor. Please enter the correct flavor"],photonflag==1&&f!=0&&f!=1&&f!=2&&f!=3&&f!=4&&f!=5&&f!=6&&f!=7&&f!=-1&&f!=-2&&f!=-3&&f!=-4&&f!=-5&&f!=-6,Print["Incorrect flavor. Please enter the correct flavor"]],Which[f==1,xPDFev[1,2][x,Q2]+xPDFev[1,-2][x,Q2],f==2,xPDFev[1,1][x,Q2]+xPDFev[1,-1][x,Q2],f==3,xPDFev[1,3][x,Q2]+xPDFev[1,-3][x,Q2],f==0,xPDFev[1,0][x,Q2],f!=0&&f!=1&&f!=2&&f!=3,Print["Incorrect polarized flavor. Please enter the correct flavor: ","0: gluon; ","1: u+ubar; ","2: d+dbar; ","3: s+sbar"]]];

xPDFRep[x_,Q2_,f_,irep_]:=If[polflag==0,Which[photonflag==0&&(f==0||f==1||f==2||f==3||f==4||f==5||f==6||f==-1||f==-2||f==-3||f==-4||f==-5||f==-6),xPDFev[irep+1,f][x,Q2],photonflag==1&&(f==0||f==1||f==2||f==3||f==4||f==5||f==6||f==7||f==-1||f==-2||f==-3||f==-4||f==-5||f==-6),xPDFev[irep+1,f][x,Q2],photonflag==0&&f!=0&&f!=1&&f!=2&&f!=3&&f!=4&&f!=5&&f!=6&&f!=-1&&f!=-2&&f!=-3&&f!=-4&&f!=-5&&f!=-6,Print["Incorrect flavor. Please enter the correct flavor"],photonflag==1&&f!=0&&f!=1&&f!=2&&f!=3&&f!=4&&f!=5&&f!=6&&f!=7&&f!=-1&&f!=-2&&f!=-3&&f!=-4&&f!=-5&&f!=-6,Print["Incorrect flavor. Please enter the correct flavor"]],Which[f==1,xPDFev[irep+1,2][x,Q2]+xPDFev[irep+1,-2][x,Q2],f==2,xPDFev[irep+1,1][x,Q2]+xPDFev[irep+1,-1][x,Q2],f==3,xPDFev[irep+1,3][x,Q2]+xPDFev[irep+1,-3][x,Q2],f==0,xPDFev[irep+1,0][x,Q2],f!=0&&f!=1&&f!=2&&f!=3,Print["Incorrect polarized flavor. Please enter the correct flavor: ","0: gluon; ","1: u+ubar; ","2: d+dbar; ","3: s+sbar"]]];

xPDFEnsemble[x_,Q2_,f_]:=If[polflag==0,Which[photonflag==0&&(f==0||f==1||f==2||f==3||f==4||f==5||f==6||f==-1||f==-2||f==-3||f==-4||f==-5||f==-6),Table[xPDFev[irep+1,f][x,Q2],{irep,1,NRep}],photonflag==1&&(f==0||f==1||f==2||f==3||f==4||f==5||f==6||f==7||f==-1||f==-2||f==-3||f==-4||f==-5||f==-6),Table[xPDFev[irep+1,f][x,Q2],{irep,1,NRep}],photonflag==0&&f!=0&&f!=1&&f!=2&&f!=3&&f!=4&&f!=5&&f!=6&&f!=-1&&f!=-2&&f!=-3&&f!=-4&&f!=-5&&f!=-6,Print["Incorrect flavor. Please enter the correct flavor"],photonflag==1&&f!=0&&f!=1&&f!=2&&f!=3&&f!=4&&f!=5&&f!=6&&f!=7&&f!=-1&&f!=-2&&f!=-3&&f!=-4&&f!=-5&&f!=-6,Print["Incorrect flavor. Please enter the correct flavor"]],Which[f==1,Table[xPDFev[irep+1,2][x,Q2]+xPDFev[irep+1,-2][x,Q2],{irep,1,NRep}],f==2,Table[xPDFev[irep+1,1][x,Q2]+xPDFev[irep+1,-1][x,Q2],{irep,1,NRep}],f==3,Table[xPDFev[irep+1,3][x,Q2]+xPDFev[irep+1,-3][x,Q2],{irep,1,NRep}],f==0,Table[xPDFev[irep+1,0][x,Q2],{irep,1,NRep}],f!=0&&f!=1&&f!=2&&f!=3,Print["Incorrect polarized flavor. Please enter the correct flavor: ","0: gluon; ","1: u+ubar; ","2: d+dbar; ","3: s+sbar"]]];

(* Table[xPDFev[irep+1,f][x,Q2],{irep,1,NRep}]; *)

xPDF[x_,Q2_,f_]:={Mean[xPDFEnsemble[x,Q2,f]],StandardDeviation[xPDFEnsemble[x,Q2,f]]};

xMin[]:=Print[xmin];
xMax[]:=Print[xmax];
Q2Min[]:=Print[Q2min];
Q2Max[]:=Print[Q2max];

NumberPDF[]:=Print[NRep];
Lam4[]:=Print[lambdaQCD4];
Lam5[]:=Print[lambdaQCD5];

Infoalphas[]:=Print[alphainfo];

ForcepolPDFtoLHAPDF[]:={polflag=0;Print["Warning: you are using dummy quark-antiquark PDFs which are dummy separated. Only q+qbar combinations have physical meaning."]};

UnForcepolPDFtoLHAPDF[]:={polflag=0;Print["Polarized PDFs properly reactivated."];};
HasPhoton[]:=If[photonflag==0,Print["The photon PDF is not available in this NNPDF grid."],Print["The photon PDF is available in this grid."],Print["An error occurred while processing the photon PDF."]];

xPDFCL[ensemble_,x_,Q2_,f_,CL_]:={sortedEnsemble=Sort[ensemble[x,Q2,f]];min=Round[((1-CL)/2)*Length[ensemble[x,Q2,f]]];max=Length[ensemble[x,Q2,f]]-Round[((1-CL)/2)*Length[ensemble[x,Q2,f]]];meanvalue=(sortedEnsemble[[max]]+sortedEnsemble[[min]])/2;error=(sortedEnsemble[[max]]-sortedEnsemble[[min]])/2;meanvalue,error};

(* Constants *)
PIGRECO=N[Pi,20];
ZETA2=N[Zeta[2],10];
ZETA3=N[Zeta[3],10];
ZETA4=N[Zeta[4],10];
c2=14/3;

(* Beta coefficient functions *)
BETA0=Table[(33-2*i)/3,{i,3,6,1}];
BETA1=Table[102-38/3*i,{i,3,6,1}];
BETA2=Table[2857/2-5033/18*i+325/54*i^2,{i,3,6,1}];
B1=Table[BETA1[[i]]/BETA0[[i]],{i,1,4,1}];
B2=Table[BETA2[[i]]/BETA0[[i]],{i,1,4,1}];

(* Implementation of AlphaMZ function *)
AlphaMZ[nfi_,mz2_,asz_,q2_,ipt_]:=
{
  q2ref=mz2;
  asref=asz/4/PIGRECO;
  asi=asref;
  t=Log[q2/q2ref];
  den=1+BETA0[[nfi-2]]*asi*t;
  alo=asi/den;

  (* LO *)
  as=alo;

  (* NLO *)
  If[ipt>=1,as=alo*(1-B1[[nfi-2]]*alo*Log[den])];

  (*NNLO *)
  If[ipt>=2,as=alo*(1+(alo*(alo-asi)*(B2[[nfi-2]]-B1[[nfi-2]]^2)+as*B1[[nfi-2]]*Log[as/asi]))];

  4*PIGRECO*as
};

(* Implementation of AlphaExactBeta *) 

(* fbeta *)
fbeta[a_,nf_,ipt_]:=
{
  fb=0;
  If[ipt==1,fb=-a^2*(BETA0[[nf-2]]+a*BETA1[[nf-2]]),If[ipt==2,fb=-a^2*(BETA0[[nf-2]]+a*(BETA1[[nf-2]]+a*BETA2[[nf-2]])),{Print["Invalid perturbative order"],Abort[];}]];
  fb
}

(* AlphaExactBeta *)
AlphaExactBeta[nf_,r20_,as0_,r2_,ipt_]:=
{
  nstep=50;
  sxth=1/6;
  nnf=nf;
  aspluto=as0;
  aspluto=aspluto/4/PIGRECO;
  as=aspluto;
  lrrat=Log[r2/r20];
  dlr=lrrat/nstep;

  If[ipt==0,as=aspluto/(1+BETA0[[nf-2]]*aspluto*lrrat),Indeterminate];
  
  If[ipt==1,For[k1=1,k1<=nstep,k1++,
      {
	xk0=dlr*fbeta[as,nnf,ipt];
	as1=Flatten[as+0.5*xk0];
	xk1=dlr*fbeta[as1,nnf,ipt];
	as2=Flatten[as+0.5*xk1];
	xk2=dlr*fbeta[as2,nnf,ipt];
	as3=Flatten[as+xk2];
	xk3=dlr*fbeta[as3,nnf,ipt];
	as=as+sxth*(xk0+2*xk1+2*xk2+xk3);
      }],Indeterminate];

  If[ipt==2,For[k1=1,k1<=nstep,k1++,
      {
	xk0=dlr*fbeta[as,nnf,ipt];
	as1=Flatten[as+0.5*xk0];
	xk1=dlr*fbeta[as1,nnf,ipt];
	as2=Flatten[as+0.5*xk1];
	xk2=dlr*fbeta[as2,nnf,ipt];
	as3=Flatten[as+xk2];
	xk3=dlr*fbeta[as3,nnf,ipt];
	as=as+sxth*(xk0+2*xk1+2*xk2+xk3);
      }]];
   
  as*4*PIGRECO

};

(* Main alphas function as VFN subroutine *)
alphas[qq2_,ipt_,imodev_]:=
  {
    q2th={mCharm^2,mBottom^2,mTop^2};
    q2=qq2;
    q2ref=MZ^2;
    If[q2>q2th[[3]],nff=6,If[q2>q2th[[2]],nff=5,If[q2>q2th[[1]],nff=4,nff=3]]];
    If[q2ref>q2th[[3]],nfi=6,If[q2ref>q2th[[2]],nfi=5,If[q2ref>q2th[[1]],nfi=4,nfi=3]]];
    alphasref=alphaSMZ[[1]];
    qq2ref=q2ref;

    While[nff!=nfi,
	     {
	       If[nff>nfi,{dnf=1;snf=1},{dnf=-1;snf=0}];
	       If[imodev==0,asi=AlphaMZ[nfi,qq2ref,alphasref,q2th[[nfi-2+snf-1]],ipt]];
	       If[imodev==1,asi=AlphaExactBeta[nfi,qq2ref,alphasref,q2th[[nfi-2+snf-1]],ipt]]
	       If[ipt>=2,{If[nff>nfi,asi=asi+(c2/(4*PIGRECO)^2)*asi^3];If[nff<nfi,asi=asi-(c2/(4*PIGRECO)^2)*asi^3]}];
	       alphasref=asi;
	       qq2ref=q2th[[nfi-2+snf-1]];
	     };nfi=nfi+dnf];

    If[nff==nfi,{If[imodev==0,alphastrong=Flatten[AlphaMZ[nfi,qq2ref,alphasref,q2,ipt]]];If[imodev==1,alphastrong=Flatten[AlphaExactBeta[nfi,qq2ref,alphasref,q2,ipt]]]},{Print["An error occurred while processing alphas"],Abort[]}];
    Flatten[alphastrong]

  };





End[];


EndPackage[];
