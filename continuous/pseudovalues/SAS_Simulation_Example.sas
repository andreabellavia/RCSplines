/****
AUTHOR: MICHAEL PALAZZOLO
DATE: 7/24/2023
PROJECT: FLEXIBLE SPLINES MODELING OF ABSOLUTE RISK IN SURVIVAL ANALYSIS
PURPOSE: SAS CODE FOR SIMULATED DATA EXAMPLE FROM JSM PRESENTATION
****/


/**************************************************************************/
/****************** Importing Simulation SAS Dataset **********************/
/**************************************************************************/
libname in " 'insert file location of data' " access=readonly;
data dat0;
	set in.simdatafinal;
all=1; *if looking at overall (not by groups) then set all = 1;
if trt = 0 then trtx = 1; *need treatment groups to be coded at (1,2);
else if trt = 1 then trtx = 2;
run;
proc contents data=dat0; run;
*==========================================================================;

/**************************************************************************/
/*************************** MACRO CALLS **********************************/
/**************************************************************************/
%include " 'insert file path location of SAS macro' \Deriving Pseudo Values Macro FINAL.sas";
%include " 'insert file path location of SAS macro' \Macro Splines Cox FINAL.sas";
%include " 'insert file path location of SAS macro' \Macro Splines Pseudo FINAL.sas";
*==========================================================================;

/**************************************************************************/
/********************* DERIVING PSEUDO VALUES *****************************/
/**************************************************************************/
%pseudoD (indata=dat0,id=id,time=eventtime,censor=status,timepoint=3,howmany=10000,modcensor='ind',stratvar=,outdata=dat1pseudo); *50 minute run time;

data dat1; *dataset with pseudo-values already derived;
	set in.simdatawithpseudo;
run;
*==========================================================================;

/**************************************************************************/
/*************************** COX REGRESSION *******************************/
/**************************************************************************/

/******************************************/
/**************** AGE *********************/
/******************************************/
/**** Baseline Covariate Dataset - AGE ****/
proc iml;
age = do(45,65,1);
create inbaseage;
append;
quit;

%ANALY_PHREG_RCS(
indat = dat1, 
baselinecov=inbaseage,
respvar = eventtime,
timeunit = years,
censvar = status,
predvar = age, 
respinterest = 1,
ptknots = 0.1 0.5 0.9,
attime = 3, 
timeatlm = 0,
bygrp = all,
bygrpnm1 = xxx,
bygrpnm2 = xxx,
xlabel = %str(Age (years)),
ylabel = %str(Absolute Risk of Overall Mortality at 3 years),
endptnm = %str(Overall Mortality),
filenmsaved = %str(temp.rtf)
);

data ageout_cox; *Saving dataset with predictions to put into one plot;
	set statsZD;
run;
/*proc print data=ageout_cox; run;*/
*==========================================;

/******************************************/
/**************** BMI *********************/
/******************************************/
/**** Baseline Covariate Dataset - BMI ****/
proc iml;
bmi = do(28,38,0.1);
create inbasebmi;
append;
quit;

%ANALY_PHREG_RCS(
indat = dat1, 
baselinecov=inbasebmi,
respvar = eventtime,
timeunit = years,
censvar = status,
predvar = bmi, 
respinterest = 1,
ptknots = 0.1 0.5 0.9,
attime = 3, 
timeatlm = 0,
bygrp = all,
bygrpnm1 = xxx,
bygrpnm2 = xxx,
xlabel = %str(BMI (kg/m2)),
ylabel = %str(Absolute Risk of Overall Mortality at 3 years),
endptnm = %str(Overall Mortality),
filenmsaved = %str(temp.rtf)
);

data bmiout_cox; *Saving dataset with predictions to put into one plot;
	set statsZD;
run;
/*proc print data=bmiout_cox; run;*/
*==========================================;
*==========================================================================;


/**************************************************************************/
/************************* GLM PSEUDO REGRESSION **************************/
/**************************************************************************/

/******************************************/
/**************** AGE *********************/
/******************************************/
%ANALY_PSEUDO_RCS(
indat = dat1, 
respvar = pseudo,
attime = 3,
timeunit = years,
predvar = age, 
respinterest = 1,
ptknots = 0.1 0.5 0.9,
bygrp = all,
bygrpnm1 = xxx,
bygrpnm2 = xxx,
xlabel = %str(Age (years)),
ylabel = %str(Absolute Risk of Overall Mortality at 3 years),
endptnm = %str(Overall Mortality),
filenmsaved = %str(temp.rtf)
);

data ageout_pseudo; *Saving dataset with predictions to put into one plot;
	set statsZD;
run;
/*proc print data=ageout_pseudo(obs=50); run;*/
*==========================================;

/******************************************/
/**************** BMI *********************/
/******************************************/
%ANALY_PSEUDO_RCS(
indat = dat1, 
respvar = pseudo,
attime = 3,
timeunit = years,
predvar = bmi, 
respinterest = 1,
ptknots = 0.1 0.5 0.9,
bygrp = all,
bygrpnm1 = xxx,
bygrpnm2 = xxx,
xlabel = %str(BMI (kg/m2)),
ylabel = %str(Absolute Risk of Overall Mortality at 3 years),
endptnm = %str(Overall Mortality),
filenmsaved = %str(temp.rtf)
);

data bmiout_pseudo; *Saving dataset with predictions to put into one plot;
	set statsZD;
run;
/*proc print data=bmiout_pseudo(obs=50); run;*/
*==========================================;
*==========================================================================;




/*%let rundate=%sysfunc(today(),mmddyyp10.);*/
/*ods pdf file=" 'insert output pdf file path location' \Final Plots &rundate..pdf";*/
/*options nodate nonumber;*/

/**************************************************************************/
/**************************** PLOT FOR AGE ********************************/
/**************************************************************************/
/*proc print data=ageout_cox; run;*/
/*proc print data=ageout_pseudo(obs=50); run;*/

data ageout_cox2;
	set ageout_cox;
age_cox = age;
est_cox = p_Failure_m;
lower_cox = p_lower_fail_m;
upper_cox = p_upper_fail_m;
keep age_cox est_cox lower_cox upper_cox;
run;
/*proc print data=ageout_cox2; run;*/

data ageout_pseudo2;
	set ageout_pseudo;
age_pseudo = age;
est_pseudo = p_Failure_m;
lower_pseudo = p_lower_fail_m;
upper_pseudo = p_upper_fail_m;
keep age_pseudo est_pseudo lower_pseudo upper_pseudo;
run;
/*proc print data=ageout_pseudo2(obs=50); run;*/

data age_final;
	set ageout_cox2 ageout_pseudo2;
label est_cox = "Cox" est_pseudo = "GLM (pseudo-values)";
run;
/*proc print data=age_final(obs=50); run;*/

proc template;
  define statgraph curveplot;
  beginGraph / border=FALSE;
    entrytitle "Cox vs Pseudo - Age";
 layout lattice / border=FALSE rows=1 columns=1;
      layout overlay / xaxisopts=(label=("Age (years)") 
                            labelattrs=(family="arial" size=10 weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin = 45 viewmax=65) griddisplay=on)
                       yaxisopts=(label="Absolute Risk of Overall Mortality at 3 years "  
                            labelattrs=(family="arial" size=10 weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin=0.2 viewmax=0.6) griddisplay=on);
      bandplot x=age_cox limitlower=lower_cox limitupper=upper_cox / fillattrs=(color=red transparency=0.6);
           seriesplot x=age_cox y=est_cox / name="cox" lineattrs=(pattern=1 thickness=3 color=red);
      bandplot x=age_pseudo limitlower=lower_pseudo limitupper=upper_pseudo / fillattrs=(color=blue transparency=0.6);
           seriesplot x=age_pseudo y=est_pseudo / name="pseudo" lineattrs=(pattern=1 thickness=3 color=blue);
         discretelegend "cox" "pseudo" / location=inside across=1 halign=left valign=top valueattrs=(family="Arial" size=9);
	   endlayout;
 endlayout;
  endGraph;
  end;
run;

ods graphics on / reset height=650px width=850px;
proc sgrender data=age_final template=curveplot;
  title;
run;
ods graphics off;


*==========================================================================;

/**************************************************************************/
/**************************** PLOT FOR BMI ********************************/
/**************************************************************************/
/*proc print data=bmiout_cox; run;*/
/*proc print data=bmiout_pseudo(obs=50); run;*/

data bmiout_cox2;
	set bmiout_cox;
bmi_cox = bmi;
est_cox = p_Failure_m;
lower_cox = p_lower_fail_m;
upper_cox = p_upper_fail_m;
keep bmi_cox est_cox lower_cox upper_cox;
run;
/*proc print data=bmiout_cox2; run;*/

data bmiout_pseudo2;
	set bmiout_pseudo;
bmi_pseudo = bmi;
est_pseudo = p_Failure_m;
lower_pseudo = p_lower_fail_m;
upper_pseudo = p_upper_fail_m;
keep bmi_pseudo est_pseudo lower_pseudo upper_pseudo;
run;
/*proc print data=ageout_pseudo2(obs=50); run;*/

data bmi_final;
	set bmiout_cox2 bmiout_pseudo2;
label est_cox = "Cox" est_pseudo = "GLM (pseudo-values)";
run;
/*proc print data=bmi_final(obs=50); run;*/

proc template;
  define statgraph curveplot;
  beginGraph / border=FALSE;
    entrytitle "Cox vs Pseudo - BMI";
 layout lattice / border=FALSE rows=1 columns=1;
      layout overlay / xaxisopts=(label=("BMI (kg/m2)") 
                            labelattrs=(family="arial" size=10 weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin=28 viewmax=38) griddisplay=on)
                       yaxisopts=(label="Absolute Risk of Overall Mortality at 3 years "  
                            labelattrs=(family="arial" size=10 weight=bold) tickvalueattrs=(family="arial" size=9) linearopts=(viewmin=0.2 viewmax=1) griddisplay=on);
      bandplot x=bmi_cox limitlower=lower_cox limitupper=upper_cox / fillattrs=(color=red transparency=0.6);
           seriesplot x=bmi_cox y=est_cox / name="cox" lineattrs=(pattern=1 thickness=3 color=red);
      bandplot x=bmi_pseudo limitlower=lower_pseudo limitupper=upper_pseudo / fillattrs=(color=blue transparency=0.6);
           seriesplot x=bmi_pseudo y=est_pseudo / name="pseudo" lineattrs=(pattern=1 thickness=3 color=blue);
         discretelegend "cox" "pseudo" / location=inside across=1 halign=left valign=top valueattrs=(family="Arial" size=9);
	   endlayout;
 endlayout;
  endGraph;
  end;
run;

ods graphics on / reset height=650px width=850px;
proc sgrender data=bmi_final template=curveplot;
  title;
run;
ods graphics off;
*==========================================================================;

/*ods pdf close;*/
