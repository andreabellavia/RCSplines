/* #############################################################################
** AUTHOR:  Michael Palazzolo and Jeong-Gun Park
** DATE:    August 29, 2024
** ----------------------------------------------------------------
** PURPOSE: Event Probability - Logistic Regression Model
** ----------------------------------------------------------------
** ############################################################################## */

*** ===============================;
*** Description of MACRO variables ;
*** ===============================;
*** INDAT = SAS dataset to be used for analysis;
*** TRIM = Specify (1) if predicted values should be trimmed to plot 1% to 99% percentile of predictor variable. This will be useful for data with extreme values;
*** CENSVAR = Binary outcome variable;
*** PREDVAR = Covariate to be splined;
*** RESPINTEREST = Output of interest: Probability of Failure (1) or Probability of Survival (2);
*** KNOTMETHOD = Method of knot placement: PERCENTILELIST (1) or RANGEFRACTIONS (2);
*** PTKNOTS = Locations of knots used for restricted cubic spline function for the covariate;
*** BYGRP = 1 or 2 groups to be analyzed by: assigned group variable must be (1, 2) for two groups or enter value of 0 for one group; 
*** BYGRPNM1 = The first group (BYGRP = 1) name to be used for plots: If one goup, then it can be any pseudo name;
*** BYGRPNM2 = The second group (BYGRP = 2) name to be used for plots: If one goup, then it can be any pseudo name;
*** COLOR1 = Color of spline curve for group 1 (or color for single spline if one group);
*** COLOR2 = Color of spline curve for group 2;
*** XLABEL = Label of X-axis;
*** XLABELSIZE = Size of text for X-axis;
*** XVALS = Values for the ticks of X-axis;
*** YLABELSIZE = Size of text for Y-axis;
*** YVALS = Values for the ticks of Y-axis;
*** Y2LABEL = Label of Y2-axis if overlaying histogram (hist=1);
*** Y2LABELSIZE = Size of text for Y2-axis if overlaying histogram (hist=1);
*** Y2VALS = Values for the ticks of Y2-axis if overlaying histogram (hist=1);
*** TITLE = Title for the plot;
*** KNOTLINES = Requests lines where the knots are placed (enter value of 1);
*** HIST = Request a histogram of the covariate that is splined overlayed on the plot (enter value of 1);
*** NBINS = Number of bins for the histogram (if hist=1);
*** GRID = Requests gridlines for the plot (enter a value of 1);
*** ================================;



*** ###################$$$$$$$$$$$$$$############################;
*** MACRO starts here;
*** #############################################################;
%macro LOGISTIC_PROB_RCS(indat=, trim=, censvar=, predvar=, respinterest=, knotmethod=, ptknots=, 
                             bygrp=, bygrpnm1=, bygrpnm2=, color1=, color2=, xlabel=, xlabelsize=, xvals=,  ylabelsize=, yvals=, 
                             y2label=, y2labelsize=, y2vals=, title=, knotlines=, hist=, nbins=, grid=);

*If no by group specified - specify 0;
%if &bygrp. = 0 %then %do;
data &indat.;
	set &indat.;
bygrpname = 1;
run;
%let bygrp = bygrpname;
%end;

/********************************************************************/
*** Determine the number of groups: 1 or 2 groups;
/********************************************************************/
proc freq data=&indat. noprint;
  tables &bygrp. / out=frqoutD outcum;
run;
proc means data=frqoutD noprint;
  var &bygrp;
  output out=frqoutED n=n;
run;
data _NULL_; set frqoutED;
  call symput("ngroups", n);
run;
%put &ngroups.;
*====================================================================;

*** ==============================================;
*** Logistic Regression Model;
*** ==============================================;
proc sort data=&indat.; by &bygrp.; run;
ods exclude all;
proc logistic data=&indat. descending;
	by &bygrp.;
%if &knotmethod. = 1 %then %do;
	effect splxvar = spline(&predvar. / naturalcubic knotmethod=PERCENTILELIST(&ptknots.));
%end;
%if &knotmethod. = 2 %then %do;
	effect splxvar = spline(&predvar. / naturalcubic knotmethod=RANGEFRACTIONS(&ptknots.));
%end;
	model &censvar. = splxvar / link=logit;
	output out=bloutD predicted=estimate lower=lower upper=upper;
run;
ods exclude close;
/*proc print data=bloutD(obs=200); run;*/

*** ==================================;
title2 "Estimates of Event Probability of Failure or Survival";
data bloutZD; set bloutD;
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
/*proc print data=bloutZD(obs=10); run;*/

proc sort data=bloutZD out=statsZD nodupkey; by &bygrp. &predvar.; run;

data hist_in;
	set &indat.;
predvar2 = &predvar;
keep predvar2;
run;

%if &trim. = 1 %then %do;
ods exclude all;
proc univariate data=&indat.; 
var &predvar.; 
output out=trim0 pctlpts=1,99 pctlpre=P;
run;
ods exclude close;
data _null_;
set trim0;
call symput('tlow',strip(left(P1)));
call symput('tup',strip(left(P99)));
run;
%put &tlow; 
%put &tup; 
data statsZD;
set statsZD;
if &tlow. <= &predvar. <= &tup.;
run;
%end;

data statsZD2;
	set statsZD hist_in;
run;


%if &knotmethod. = 1 %then %do;
/**** Percentiles ****/
ods exclude all;
proc univariate data=&indat.;
%if &ngroups. = 2 %then %do;
	by &bygrp.;
%end;
	var &predvar.;
	output out=stats0 median=median pctlpts=1 to 99 by 1 pctlpre=P;
run;
ods exclude close;
/*proc print data=stats0; run;*/

/**** Data for knot lines ****/
data stats1;
	set stats0;
knotloc="&ptknots.";
numknots = countw(knotloc);
k1 = scan(knotloc,1,' ');
k2 = scan(knotloc,2,' ');
k3 = scan(knotloc,3,' ');
k4 = scan(knotloc,4,' ');
k5 = scan(knotloc,5,' ');
run;
proc print data=stats1; run;
data _null_;
	set stats1;
call symputx('numknot',numknots,'g');
run;
%if &ngroups. = 1 %then %do;
%if &numknot. = 3 %then %do;
data _null_;
	set stats1;
call symputx('kn1',strip(left(k1)),'g');
call symputx('kn2',strip(left(k2)),'g');
call symputx('kn3',strip(left(k3)),'g');
run;
data _null_;
	set stats1;
call symputx('k1',P&kn1.,'g');
call symputx('k2',P&kn2.,'g');
call symputx('k3',P&kn3.,'g');
run;
%end;
%if &numknot. = 4 %then %do;
data _null_;
	set stats1;
call symputx('kn1',strip(left(k1)),'g');
call symputx('kn2',strip(left(k2)),'g');
call symputx('kn3',strip(left(k3)),'g');
call symputx('kn4',strip(left(k4)),'g');
run;
data _null_;
	set stats1;
call symputx('k1',P&kn1.,'g');
call symputx('k2',P&kn2.,'g');
call symputx('k3',P&kn3.,'g');
call symputx('k4',P&kn4.,'g');
run;
%end;
%if &numknot. = 5 %then %do;
data _null_;
	set stats1;
call symputx('kn1',strip(left(k1)),'g');
call symputx('kn2',strip(left(k2)),'g');
call symputx('kn3',strip(left(k3)),'g');
call symputx('kn4',strip(left(k4)),'g');
call symputx('kn5',strip(left(k5)),'g');
run;
data _null_;
	set stats1;
call symputx('k1',P&kn1.,'g');
call symputx('k2',P&kn2.,'g');
call symputx('k3',P&kn3.,'g');
call symputx('k4',P&kn4.,'g');
call symputx('k5',P&kn5.,'g');
run;
%end;
%end;
%if &ngroups. = 2 %then %do;
%if &numknot. = 3 %then %do;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('kn1g1',strip(left(k1)),'g');
call symputx('kn2g1',strip(left(k2)),'g');
call symputx('kn3g1',strip(left(k3)),'g');
run;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('k1g1',P&kn1g1.,'g');
call symputx('k2g1',P&kn2g1.,'g');
call symputx('k3g1',P&kn3g1.,'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('kn1g2',strip(left(k1)),'g');
call symputx('kn2g2',strip(left(k2)),'g');
call symputx('kn3g2',strip(left(k3)),'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('k1g2',P&kn1g2.,'g');
call symputx('k2g2',P&kn2g2.,'g');
call symputx('k3g2',P&kn3g2.,'g');
run;
%end;
%if &numknot. = 4 %then %do;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('kn1g1',strip(left(k1)),'g');
call symputx('kn2g1',strip(left(k2)),'g');
call symputx('kn3g1',strip(left(k3)),'g');
call symputx('kn4g1',strip(left(k4)),'g');
run;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('k1g1',P&kn1g1.,'g');
call symputx('k2g1',P&kn2g1.,'g');
call symputx('k3g1',P&kn3g1.,'g');
call symputx('k4g1',P&kn4g1.,'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('kn1g2',strip(left(k1)),'g');
call symputx('kn2g2',strip(left(k2)),'g');
call symputx('kn3g2',strip(left(k3)),'g');
call symputx('kn4g2',strip(left(k4)),'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('k1g2',P&kn1g2.,'g');
call symputx('k2g2',P&kn2g2.,'g');
call symputx('k3g2',P&kn3g2.,'g');
call symputx('k4g2',P&kn4g2.,'g');
run;
%end;
%if &numknot. = 5 %then %do;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('kn1g1',strip(left(k1)),'g');
call symputx('kn2g1',strip(left(k2)),'g');
call symputx('kn3g1',strip(left(k3)),'g');
call symputx('kn4g1',strip(left(k4)),'g');
call symputx('kn5g1',strip(left(k5)),'g');
run;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('k1g1',P&kn1g1.,'g');
call symputx('k2g1',P&kn2g1.,'g');
call symputx('k3g1',P&kn3g1.,'g');
call symputx('k4g1',P&kn4g1.,'g');
call symputx('k5g1',P&kn5g1.,'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('kn1g2',strip(left(k1)),'g');
call symputx('kn2g2',strip(left(k2)),'g');
call symputx('kn3g2',strip(left(k3)),'g');
call symputx('kn4g2',strip(left(k4)),'g');
call symputx('kn5g2',strip(left(k5)),'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('k1g2',P&kn1g2.,'g');
call symputx('k2g2',P&kn2g2.,'g');
call symputx('k3g2',P&kn3g2.,'g');
call symputx('k4g2',P&kn4g2.,'g');
call symputx('k5g2',P&kn5g2.,'g');
run;
%end;
%end;
%end;
%if &knotmethod. = 2 %then %do;
/**** Range Fraction ****/
ods exclude all;
proc univariate data=&indat.;
%if &ngroups. = 2 %then %do;
	by &bygrp.;
%end;
	var &predvar.;
	output out=stats0 median=median min=min max=max;
run;
ods exclude close;
/*proc print data=stats0; run;*/

/**** Data for knot lines ****/
data stats1;
	set stats0;
range=max-min;
knotloc="&ptknots.";
numknots = countw(knotloc);
pt1 = scan(knotloc,1,' ');
pt2 = scan(knotloc,2,' ');
pt3 = scan(knotloc,3,' ');
pt4 = scan(knotloc,4,' ');
pt5 = scan(knotloc,5,' ');
k1 = min + (pt1*range);
k2 = min + (pt2*range);
k3 = min + (pt3*range);
k4 = min + (pt4*range);
k5 = min + (pt5*range);
run;
/*proc print data=stats1; run;*/
data _null_;
	set stats1;
call symputx('numknot',numknots,'g');
run;
%if &ngroups. = 1 %then %do;
data _null_;
	set stats1;
call symputx('k1',k1,'g');
call symputx('k2',k2,'g');
call symputx('k3',k3,'g');
call symputx('k4',k4,'g');
call symputx('k5',k5,'g');
run;
%end;
%if &ngroups. = 2 %then %do;
data _null_;
	set stats1;
if &bygrp. = 1;
call symputx('k1g1',k1,'g');
call symputx('k2g1',k2,'g');
call symputx('k3g1',k3,'g');
call symputx('k4g1',k4,'g');
call symputx('k5g1',k5,'g');
run;
data _null_;
	set stats1;
if &bygrp. = 2;
call symputx('k1g2',k1,'g');
call symputx('k2g2',k2,'g');
call symputx('k3g2',k3,'g');
call symputx('k4g2',k4,'g');
call symputx('k5g2',k5,'g');
run;
%end;
%end;

/*proc print data=stats1; run;*/
/*%put &k1. &k2. &k3. &k4. &k5. &numknot.;*/


/**** Axis Options ****/
data axisopts;
xaxisvals="&xvals.";
yaxisvals="&yvals.";
y2axisvals="&y2vals.";
xmin = scan(xaxisvals,1,' ');
xmax = scan(xaxisvals,-1,' ');
ymin = scan(yaxisvals,1,' ');
ymax = scan(yaxisvals,-1,' ');
y2min = scan(y2axisvals,1,' ');
y2max = scan(y2axisvals,-1,' ');
run;
/*proc print data=axisopts; run;*/
data _null_;
	set axisopts;
call symputx('xmin',xmin,'g');
call symputx('xmax',xmax,'g');
call symputx('ymin',ymin,'g');
call symputx('ymax',ymax,'g');
call symputx('y2min',y2min,'g');
call symputx('y2max',y2max,'g');
run;

/**** Grid Displayed ****/
%if &grid. = 1 %then %do;
	%let gridval=on;
%end;
%else %do;
	%let gridval=off;
%end;

*** #####################################################;
*** PROC TEMPLATE procedure;
*** #####################################################;
proc template;
  define statgraph curveplot;
  beginGraph / border=FALSE;
    entrytitle "&title.";
  layout lattice / border=FALSE rows=1 columns=1;
      layout overlay / xaxisopts=(label=("&xlabel") 
                            labelattrs=(family="arial" size=&xlabelsize. weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin=&xmin. viewmax=&xmax. tickvaluelist=(&xvals.) minorgrid=TRUE) griddisplay=&gridval.)
                       yaxisopts=(label="&ylabeltitle"  
                            labelattrs=(family="arial" size=&ylabelsize. weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin=&ymin. viewmax=&ymax. tickvaluelist=(&yvals.) tickvalueformat=PERCENT10. minorgrid=TRUE) griddisplay=&gridval.)
%if &hist. = 1 %then %do;
					   y2axisopts=(label="&y2label"  
                            labelattrs=(family="arial" size=&y2labelsize. weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin=&y2min. viewmax=&y2max. tickvaluelist=(&y2vals.) tickvalueformat=PERCENT10.) /*griddisplay=&gridval.*/)
%end;
;

%if &ngroups eq 1 %then %do;
      bandplot x=&predvar limitlower=lowerestimA1 limitupper=upperestimA1 / fillattrs=(color=&color1. transparency=0.75);
           seriesplot x=&predvar y=estimA1 / lineattrs=(pattern=1 thickness=3 color=&color1.);
%end; 
%else %if &ngroups eq 2 %then %do;
     bandplot x=&predvar limitlower=lowerestimA1 limitupper=upperestimA1 / fillattrs=(color=&color1. transparency=0.75);
           seriesplot x=&predvar y=estimA1 / name="grp1" lineattrs=(pattern=1 thickness=3 color=&color1.);
     bandplot x=&predvar limitlower=lowerestimA2 limitupper=upperestimA2 / fillattrs=(color=&color2. transparency=0.75);
           seriesplot x=&predvar y=estimA2 / name="grp2" lineattrs=(pattern=1 thickness=3 color=&color2.);
         discretelegend "grp1" "grp2" / location=outside halign=center valign=bottom valueattrs=(family="Arial" size=9);
%end;

%if &hist. = 1 %then %do;
	histogram predvar2 / yaxis=y2 scale=PROPORTION datatransparency=0.7 NBINS=&nbins.;
%end;

%if &knotlines. = 1 and &ngroups. = 1 %then %do;
	%if &numknot. = 3 %then %do;
	referenceline x=&k1. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k2. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k3. / lineattrs=(color=black pattern=Dot);
	%end;
	%if &numknot. = 4 %then %do;
	referenceline x=&k1. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k2. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k3. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k4. / lineattrs=(color=black pattern=Dot);
	%end;
	%if &numknot. = 5 %then %do;
	referenceline x=&k1. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k2. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k3. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k4. / lineattrs=(color=black pattern=Dot);
	referenceline x=&k5. / lineattrs=(color=black pattern=Dot);
	%end;
%end;
%if &knotlines. = 1 and &ngroups. = 2 %then %do;
	%if &numknot. = 3 %then %do;
	referenceline x=&k1g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k2g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k3g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k1g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k2g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k3g2. / lineattrs=(color=&color2. pattern=Dot);
	%end;
	%if &numknot. = 4 %then %do;
	referenceline x=&k1g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k2g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k3g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k4g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k1g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k2g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k3g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k4g2. / lineattrs=(color=&color2. pattern=Dot);
	%end;
	%if &numknot. = 5 %then %do;
	referenceline x=&k1g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k2g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k3g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k4g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k5g1. / lineattrs=(color=&color1. pattern=Dot);
	referenceline x=&k1g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k2g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k3g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k4g2. / lineattrs=(color=&color2. pattern=Dot);
	referenceline x=&k5g2. / lineattrs=(color=&color2. pattern=Dot);
	%end;
%end;

     endlayout;
  endlayout;
  endGraph;
  end;
run;
*** ------------------------------;

proc sgrender data=statsZD2 template=curveplot;
  title;
  footnote;
run;

%mend LOGISTIC_PROB_RCS;
