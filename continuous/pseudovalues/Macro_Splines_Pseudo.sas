/* #############################################################################
** ORIGINAL AUTHOR: Jeong-Gun Park
** EDITED: Michael Palazzolo
** DATE: July 24, 2023
** ----------------------------------------------------------------
** PURPOSE: GLM Pseudo-values Splines Macro
** ----------------------------------------------------------------
** This Macro also provides 95% confidence bands;
** ############################################################################## */


*** ===============================;
*** Description of MACRO variables ;
*** ===============================;
*** INDAT = SAS dataset to be used for analysis;
*** RESPVAR = Response(pseudo-value) variable;
*** ATTIME = Time at which the pseudo-values were derived at;
*** TIMEUNIT = Time unit for time (e.g. days,years,ect.);
*** PREDVAR = Covariate to be splined;
*** RESPINTEREST = Output of interest: Failure Probability (1) or Survival Probability (2);
*** PTKNOTS = Locations of knots used for restricted cubic spline function for the covariate ('rangefractions' method);
*** BYGRP = 1 or 2 groups to be analyzed by: assigned value must be (1) for one group or (1, 2) for two groups; 
*** BYGRPNM1 = The first group (BYGRP = 1) name to be used for plots: If one goup, then it can be any name;
*** BYGRPNM2 = The second group (BYGRP = 2) name to be used for plots: If one goup, then it can be any name;
*** XLABEL = Label of X-axis;
*** YLABELVAR = Label of Y-axis;
*** ENDPTNM = Name of endpoint;
*** FILENMSAVED = Name of graphic file that is going to be saved if uncommented below;
*** ================================;


*** ###################$$$$$$$$$$$$$$############################;
*** MACRO starts here;
*** #############################################################;
%macro ANALY_PSEUDO_RCS(indat=, respvar=, attime=, timeunit=, predvar=, respinterest=, ptknots=, bygrp=, bygrpnm1=, bygrpnm2=, xlabel=, ylabel=, endptnm=, filenmsaved=);
*** ------------------------------------;
*** Determine the number of groups: 1 or 2 groups;
proc freq data=&indat noprint;
  tables &bygrp / out=frqoutD outcum;
run;
proc means data=frqoutD noprint;
  var &bygrp;
  output out=frqoutED n=n;
run;
data _NULL_; set frqoutED;
  call symput("ngroups", n);
run;
%put &ngroups.;

*** -----------------------------;

*** ==============================================;
proc sort data=&indat.; by &bygrp.; run;
ods exclude all;
proc glimmix data=&indat.;
	by &bygrp.;
	effect splxvar = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.));
	model &respvar. = splxvar;
	output out=bloutD predicted=estimate lcl=lower ucl=upper;
run;
ods exclude close;
/*proc print data=bloutD(obs=1000); run;*/

*** ==================================;
title2 "Estimates of Event Probability of Failure or Survival";
data bloutED; set bloutD;
  if &respinterest. eq 1 then do;
    call symput("ylabeltitle", "Probability of Failure");
    if &bygrp. eq 1 then do;
      *** Prob of failure;
      estimA1 = estimate;
      lowerestimA1 = upper;
   	  upperestimA1 = lower;
    end;
    if &bygrp. eq 2 then do;
      *** Prob of failure;
      estimA2 = estimate;
      lowerestimA2 = upper;
      upperestimA2 = lower;
    end;
  end;

  if &respinterest. eq 2 then do;
    call symput("ylabeltitle", "Probability of Survival");
    if &bygrp. eq 1 then do;
      *** Prob of survival;
      estimA1 = 1 - estimate;
      lowerestimA1 = 1 - upper;
   upperestimA1 = 1 - lower;
    end;
    if &bygrp. eq 2 then do;
      *** Prob of survival;
      estimA2 = 1 - estimate;
      lowerestimA2 = 1 - upper;
      upperestimA2 = 1 - lower;
    end;
  end;
  label estimA1="&bygrpnm1" estimA1="&bygrpnm1" estimA2="&bygrpnm2" estimA2="&bygrpnm2";
run;
/*proc print data=bloutED(obs=10); run;*/

%if &ngroups.=1 %then %do;
*** ------------------------------------;
*** Cubic spline smooting;
proc univariate data=bloutED noprint;
  class &predvar.;
  var estimA1 lowerestimA1 upperestimA1;
  output out=bloutFD n=nv mean=Failure_m lower_fail_m upper_fail_m;
run;

ods exclude all;
proc glmselect data=bloutFD;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.)); 
  model Failure_m = spl / selection=none ;  
  output out=predmD predicted;
run;
proc sort data=predmD; by &predvar.; run;

proc glmselect data=bloutFD;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.));
  model lower_fail_m = spl / selection=none ;
  output out=predlD predicted;
run;
proc sort data=predlD; by &predvar.; run;

proc glmselect data=bloutFD;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.)); 
  model upper_fail_m = spl / selection=none ;
  output out=preduD predicted;
run;
proc sort data=preduD; by &predvar.; run;
ods exclude close;
*** -----------------------------------;
data statsZD; merge predmD predlD preduD; by &predvar.;  
run;
/*proc print data=statsZD(obs=10); run;*/
%end;

%if &ngroups.=2 %then %do;
proc univariate data=bloutED noprint;
  where &bygrp. = 1;
  class &predvar.;
  var estimA1 lowerestimA1 upperestimA1;
  output out=bloutFD1 n=nv1 mean=Failure_m1 lower_fail_m1 upper_fail_m1;
run;

proc univariate data=bloutED noprint;
  where &bygrp. = 2;
  class &predvar.;
  var estimA2 lowerestimA2 upperestimA2;
  output out=bloutFD2 n=nv2 mean=Failure_m2 lower_fail_m2 upper_fail_m2;
run;

*** Group = 1 smoothing;
ods exclude all;
proc glmselect data=bloutFD1;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.)); 
  model Failure_m1 = spl / selection=none ;  
  output out=predmD1 predicted;
run;
proc sort data=predmD1; by &predvar.; run;

proc glmselect data=bloutFD1;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.));
  model lower_fail_m1 = spl / selection=none ;
  output out=predlD1 predicted;
run;
proc sort data=predlD1; by &predvar.; run;

proc glmselect data=bloutFD1;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.)); 
  model upper_fail_m1 = spl / selection=none ;
  output out=preduD1 predicted;
run;
proc sort data=preduD1; by &predvar.; run;
ods exclude close;

*** Group = 2 smoothing;
ods exclude all;
proc glmselect data=bloutFD2;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.)); 
  model Failure_m2 = spl / selection=none ;  
  output out=predmD2 predicted;
run;
proc sort data=predmD2; by &predvar.; run;

proc glmselect data=bloutFD2;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.));
  model lower_fail_m2 = spl / selection=none ;
  output out=predlD2 predicted;
run;
proc sort data=predlD2; by &predvar.; run;

proc glmselect data=bloutFD2;
  effect spl = spline(&predvar. / naturalcubic knotmethod=rangefractions(&ptknots.)); 
  model upper_fail_m2 = spl / selection=none ;
  output out=preduD2 predicted;
run;
proc sort data=preduD2; by &predvar.; run;
ods exclude close;

data statsZD; merge predmD1 predlD1 preduD1 predmD2 predlD2 preduD2; 
by &predvar.; 
label p_Failure_m1="&bygrpnm1" p_Failure_m1="&bygrpnm1" p_Failure_m2="&bygrpnm2" p_Failure_m2="&bygrpnm2"; 
run;
%end;

*** #####################################################;
*** PROC TEMPLATE procedure;
*** #####################################################;
proc template;
  define statgraph curveplot;
  beginGraph / border=FALSE; * designwidth=750px designheight=1000px;
    entrytitle "&ylabeltitle. of &endptnm. and 95% CI at &attime &timeunit";
 layout lattice / border=FALSE rows=1 columns=1;
      layout overlay / xaxisopts=(label=("&xlabel") 
                            labelattrs=(family="arial" size=10 weight=bold) tickvalueattrs=(family="arial" size=9) griddisplay=on)
                       yaxisopts=(label="&ylabel"  
                            labelattrs=(family="arial" size=10 weight=bold) tickvalueattrs=(family="arial" size=9) griddisplay=on);

   %if &ngroups eq 1 %then %do;
      bandplot x=&predvar limitlower=p_lower_fail_m limitupper=p_upper_fail_m / fillattrs=(color=blue transparency=0.6);
           seriesplot x=&predvar y=p_Failure_m / lineattrs=(pattern=1 thickness=3 color=red);
   %end; 
   %else %if &ngroups eq 2 %then %do;
     bandplot x=&predvar limitlower=p_lower_fail_m1 limitupper=p_upper_fail_m1 / fillattrs=(color=red transparency=0.6);
           seriesplot x=&predvar y=p_Failure_m1 / name="grp1" lineattrs=(pattern=1 thickness=3 color=red);
     bandplot x=&predvar limitlower=p_lower_fail_m2 limitupper=p_upper_fail_m2 / fillattrs=(color=blue transparency=0.6);
           seriesplot x=&predvar y=p_Failure_m2 / name="grp2" lineattrs=(pattern=1 thickness=3 color=blue);
         discretelegend "grp1" "grp2" / location=inside across=1 halign=left valign=top valueattrs=(family="Arial" size=9);
   %end;

      endlayout;
 endlayout;
  endGraph;
  end;
run;
*** ------------------------------;

*** =================================================;
*** ODS Graphics;
*** =================================================;
options orientation=landscape;

ods graphics on / reset height=650px width=850px;
ods escapechar='^';

/*ods rtf file="&filenmsaved" bodytitle; */
proc sgrender data=statsZD template=curveplot;
  title;
  footnote;
run;
/*ods rtf close;*/

ods graphics off;
*** End of Graphics;

%mend ANALY_PSEUDO_RCS;
*** ##############################################################;
*** MACRO ends here;
*** ##############################################################;


