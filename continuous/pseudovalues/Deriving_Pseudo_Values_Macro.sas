/* 
#############################################################################
** AUTHOR:  Michael Palazzolo
** DATE:    March 21, 2025
** ----------------------------------------------------------------
** PURPOSE: This macro derives pseudo-values for each observation based on Kaplan-Meier estimates for the survival curve.
Pseudo-values are based on the leave-one-out KM estimator.
This macro derives pseudo-values for a specific endpoint.
The output dataset will contain the variable 'pseudo' which contains the pseudo-values.
This macros run time is dependent on the size of the input dataset.
** ----------------------------------------------------------------
##############################################################################
*/


*** ===============================;
*** Description of MACRO variables ;
*** ===============================;
*** INDATA = SAS dataset for pseudo values;
*** ID = Unique subject ID;
*** TIME = Time variable in survival analysis;
*** CENSOR = Censoring indicator variable;
*** TIMEPOINT = Timepoint to calculate the pseudo values at (can me multiple timepoints);
*** HOWMANY = Number of observations in the dataset;
*** MODCENSOR = Censoring type ('ind' for independent, 'strat' for stratified);
*** STRATVAR = Strata for stratified censoring (must code the variable as 0,1);
*** OUTDATA = Output dataset;
*** ================================;




*** #############################################################;
*** MACRO starts here;
*** #############################################################;
%macro pseudoD (indata=,id=,time=,censor=,timepoint=,howmany=,modcensor=,stratvar=,outdata=);
%if &modcensor. = 'ind' %then %do;
/***************************************************************************************/
/************************ OVERALL KM ESTIMATE (THETA) **********************************/
/***************************************************************************************/
proc lifetest data=&indata. noprint plots=none timelist=&timepoint. reduceout outsurv=sall;
time &time.*&censor.(0);
run;
/*proc print data=sall; run;*/

data sall;
set sall;
theta = survival;
keep theta;
run;

data sout;
set &indata.;
run;
*=======================================================================================;

/***************************************************************************************/
/************ KM ESTIMATE BASED ON LEAVE-ONE-OUT ESTIMATOR (THETAMINI) *****************/
/***************************************************************************************/
%do ip=1 %to &howmany.;
proc lifetest data=&indata. noprint plots=none timelist=&timepoint. reduceout outsurv=salli;
time &time.*&censor.(0);
where &id. ne &ip.;
run;

data salli2;
set salli;
thetamini = survival;
&id. = &ip.;
keep &id. thetamini;
run;
/*proc print data=salli2; run;*/

data souti;
merge salli2 sall;
run;
/*proc print data=souti; run;*/

proc sort data=sout; by &id.; run;
proc sort data=souti; by &id.; run;
data sout;
merge sout souti;
by &id.;
run;
/*proc print data=sout; run;*/
%end;
*=======================================================================================;

/***************************************************************************************/
/*************************** PSEUDO-VALUE CALCULATION **********************************/
/***************************************************************************************/
data sfinal;
set sout;
pseudo = 1 - (&howmany. * theta - (&howmany.-1) * thetamini);
run;

data &outdata.;
	set sfinal;
run;
*=======================================================================================;
%end;


%if &modcensor. = 'strat' %then %do;
/***************************************************************************************/
/********************** DETERMINE HOWMANY FOR EACH GROUP *******************************/
/***************************************************************************************/
ods exclude all;
proc freq data=&data.;
tables &stratvar.;
ods output OneWayFreqs=howm0;
run;
ods output close;
ods exclude close;
proc print data=howm0; run;
data _null_;
set howm0;
if &stratvar. = 0 then do;
	call symput ('howmany0',trim(left(Frequency)));
end;
if &stratvar. = 1 then do;
	call symput ('howmany1',trim(left(Frequency)));
end;
run;
%put &howmany0.;
%put &howmany1.;

data grp0;
set &data.;
where &stratvar. = 0;
id0 = _n_;
run;
data grp1;
set &data.;
where &stratvar. = 1;
id1 = _n_;
run;
*=======================================================================================;


/***************************************************************************************/
/************************ OVERALL KM ESTIMATE (Group 0) ********************************/
/***************************************************************************************/
proc lifetest data=grp0 noprint plots=none timelist=&timepoint. reduceout outsurv=sall;
time &time.*&censor.(0);
run;
/*proc print data=sall; run;*/

data sall;
set sall;
theta = survival;
keep theta;
run;

data sout;
set grp0;
run;
*=======================================================================================;

/***************************************************************************************/
/************ KM ESTIMATE BASED ON LEAVE-ONE-OUT ESTIMATOR (THETAMINI) *****************/
/***************************************************************************************/
%do ip=1 %to &howmany0.;
proc lifetest data=grp0 noprint plots=none timelist=&timepoint. reduceout outsurv=salli;
time &time.*&censor.(0);
where id0 ne &ip.;
run;

data salli2;
set salli;
thetamini = survival;
id0 = &ip.;
keep id0 thetamini;
run;
/*proc print data=salli2; run;*/

data souti;
merge salli2 sall;
run;
/*proc print data=souti; run;*/

proc sort data=sout; by id0; run;
proc sort data=souti; by id0; run;
data sout;
merge sout souti;
by id0;
run;
/*proc print data=sout; run;*/
%end;
*=======================================================================================;

/***************************************************************************************/
/*************************** PSEUDO-VALUE CALCULATION **********************************/
/***************************************************************************************/
data grp0_final;
set sout;
pseudo = 1 - (&howmany0. * theta - (&howmany0.-1) * thetamini);
run;

proc datasets;
delete sout;
run;
*=======================================================================================;

/***************************************************************************************/
/************************ OVERALL KM ESTIMATE (Group 1) **********************************/
/***************************************************************************************/
proc lifetest data=grp1 noprint plots=none timelist=&timepoint. reduceout outsurv=sall;
time &time.*&censor.(0);
run;
/*proc print data=sall; run;*/

data sall;
set sall;
theta = survival;
keep theta;
run;

data sout;
set grp1;
run;
*=======================================================================================;

/***************************************************************************************/
/************ KM ESTIMATE BASED ON LEAVE-ONE-OUT ESTIMATOR (THETAMINI) *****************/
/***************************************************************************************/
%do ip=1 %to &howmany1.;
proc lifetest data=grp1 noprint plots=none timelist=&timepoint. reduceout outsurv=salli;
time &time.*&censor.(0);
where id1 ne &ip.;
run;

data salli2;
set salli;
thetamini = survival;
id1 = &ip.;
keep id1 thetamini;
run;
/*proc print data=salli2; run;*/

data souti;
merge salli2 sall;
run;
/*proc print data=souti; run;*/

proc sort data=sout; by id1; run;
proc sort data=souti; by id1; run;
data sout;
merge sout souti;
by id1;
run;
/*proc print data=sout; run;*/
%end;
*=======================================================================================;

/***************************************************************************************/
/*************************** PSEUDO-VALUE CALCULATION **********************************/
/***************************************************************************************/
data grp1_final;
set sout;
pseudo = 1 - (&howmany1. * theta - (&howmany1.-1) * thetamini);
run;

proc datasets;
delete sout;
run;
*=======================================================================================;

data sfinal;
set grp0_final grp1_final;
run;

data &outdata.;
	set sfinal;
run;
%end;

%mend pseudoD;
