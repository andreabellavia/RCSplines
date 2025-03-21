/* #############################################################################
** ORIGINAL AUTHOR: Jeong-Gun Park
** EDITED: Michael Palazzolo
** DATE: July 24, 2023
** ----------------------------------------------------------------
** PURPOSE: This macro estimates probability of study events (1 - survival probability) and cumulative hazard of study event at a specified time when Cox regression model 
**             is fitted (using PROC PHREG) including a continuous covariate smoothed using restricted cubic spline. 
** This macro creates plots of study events against the continuous covariate by the specified group using SAS PROC TEMPLATE procedure.
** Please see below information about the baseline covariate dataset.
** ----------------------------------------------------------------
** This Macro also provides 95% confidence bands;
** ############################################################################## */

*** ===============================;
*** Baseline Covariate Dataset     ;
*** ===============================;
*** NOTES: The macro calls for a baseline covariate dataset to be specified. This should be a SAS dataset.
*** OPTIONS:
***      (1) Use PROC IML to specify the range of the predictor variable you would like to make estimations.
***      (2) Use the input dataset and keep only the predictor variable (there will be a predicted probability for each observation in the dataset).
***
*** PROC IML EXAMPLE CODE:
***                 proc iml;
***                 x = do(28,38,0.1); *will create a dataset that has values for 'x' of 28 to 38 in 0.1 increments;
***                 create inbase;
***                 append;
***                 quit;
*** 'x' should be the predictor variable | 'inbase' is the output dataset that you will use for the macro parameter 'BASELINECOV'
*** ================================;

*** ===============================;
*** Description of MACRO variables ;
*** ===============================;
*** INDAT = SAS dataset to be used for analysis;
*** BASELINECOV = SAS dataset of baseline covariates;
*** RESPVAR = Response(time) variable in survival analysis;
*** TIMEUNIT = Time unit for response(time) variable;
*** CENSVAR = Censoring indicator variable;
*** PREDVAR = Covariate to be splined (predictor variable);
*** RESPINTEREST = Output of interest: Event Probability (1) or Cummulative Hazard (2);
*** PTKNOTS = Locations of knots used for restricted cubic spline function for the covariate ('rangefractions' method);
*** ATTIME = Time at which event estimation will be made; 
*** TIMEATLM = time at landmarked  (This is for landmark analysis.  If not landmarked, use zero (0));
*** BYGRP = 1 or 2 groups to be analyzed by: assigned value must be (1) for one group or (1, 2) for two groups; 
*** BYGRPNM1 = The first group (BYGRP = 1) name to be used for plots: If one goup, then it can be any pseudo name;
*** BYGRPNM2 = The second group (BYGRP = 2) name to be used for plots: If one goup, then it can be any pseudo name;
*** XLABEL = Label of X-axis;
*** YLABELVAR = Label of Y-axis;
*** ENDPTNM = Name of endpoint;
*** FILENMSAVED = Name of graphic file that is going to be saved if uncommented below;
*** ================================;


*** ###################$$$$$$$$$$$$$$############################;
*** MACRO starts here;
*** #############################################################;
%macro ANALY_PHREG_RCS(indat=, baselinecov=, respvar=, timeunit=, censvar=, predvar=, respinterest=, ptknots=, attime=, timeatlm=, bygrp=, bygrpnm1=, bygrpnm2=, xlabel=, ylabel=, endptnm=, filenmsaved=);

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
proc phreg data=&indat.;
  by &bygrp.;
  effect splxvar = spline(&predvar / naturalcubic knotmethod=rangefractions(&ptknots.));
  model &respvar*&censvar(0) = splxvar;
  baseline out=bloutD covariates=&baselinecov. survival=survival lower=lower upper=upper cumhaz=cumhaz lowercumhaz=lowercumhaz uppercumhaz=uppercumhaz;
  hazardratio &predvar.;
  ods output ParameterEstimates=estD;
  ods output HazardRatios=hrD;
run;
/*proc print data=bloutD(obs=20); run;*/
ods exclude close;

*** ==================================;
title2 "Estimates of Event Probability or Cumulative Hazard";
data bloutED; set bloutD;
  if &respinterest. eq 1 then do;
    call symput("ylabeltitle", "Probability");
    if &bygrp. eq 1 then do;
      *** Prob of failure;
      estimA1 = 1 - survival;
      lowerestimA1 = 1 - upper;
   upperestimA1 = 1 - lower;
    end;
    if &bygrp. eq 2 then do;
      *** Prob of failure;
      estimA2 = 1 - survival;
      lowerestimA2 = 1 - upper;
      upperestimA2 = 1 - lower;
    end;
  end;

  if &respinterest. eq 2 then do;
    call symput("ylabeltitle", "Cummulative Hazard");
      if &bygrp. eq 1 then do;
      *** Cummulative hazard;
      estimA1 = cumhaz;
      lowerestimA1 = lowercumhaz;
      upperestimA1 = uppercumhaz;
    end;

    if &bygrp. eq 2 then do;
      *** Cummulative hazard;
      estimA2 = cumhaz;
      lowerestimA2 = lowercumhaz;
      upperestimA2 = uppercumhaz;
    end;
  end;

 label estimA1="&bygrpnm1" estimA1="&bygrpnm1" estimA2="&bygrpnm2" estimA2="&bygrpnm2";
run;
proc sort data=bloutED; by &bygrp &predvar &respvar; run;
*proc print data=bloutED(obs=20); run;
*** ======================================;
title2 "Choose the time of interest";
data timeD(keep=&bygrp &respvar); set bloutED;
  if &respvar le (&attime. - &timeatlm.);
run;
proc sort data=timeD; by &bygrp &respvar; run;
*** --------;
data _NULL_; set timeD; by &bygrp;
  if last.&bygrp then do; 
    if &bygrp eq 1 then call symput("attime1", &respvar);
    if &bygrp eq 2 then call symput("attime2", &respvar);
  end;
run;

*** ======================================;
title2 "Final data for plots";
%if &ngroups eq 1 %then %do;
data bloutZD; set bloutED;
    if &bygrp eq 1 &  round(&respvar, 0.000001) eq round(&attime1, 0.000001);
run;
%end; 
%else %if &ngroups eq 2 %then %do;
data bloutZD; set bloutED;
    if &bygrp eq 1 &  round(&respvar, 0.000001) eq round(&attime1, 0.000001)
       | &bygrp eq 2 &  round(&respvar, 0.000001) eq round(&attime2, 0.000001);
run;
%end;

%if &ngroups.=1 %then %do;
*** ------------------------------------;
*** Cubic spline smooting;
proc univariate data=bloutZD noprint;
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
proc univariate data=bloutZD noprint;
  where &bygrp. = 1;
  class &predvar.;
  var estimA1 lowerestimA1 upperestimA1;
  output out=bloutFD1 n=nv1 mean=Failure_m1 lower_fail_m1 upper_fail_m1;
run;

proc univariate data=bloutZD noprint;
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
 layout lattice / border=FALSE rows=1 columns=1;* rowgutter=8;
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

%mend ANALY_PHREG_RCS;
*** ##############################################################;
*** MACRO ends here;
*** ##############################################################;



