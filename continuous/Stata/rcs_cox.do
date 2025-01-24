* Estimating and presenting nonlinear associations with Restricted Cubic Splines
*     Andrea Discacciati, Michael G. Palazzolo, Jeong-Gun Park, 
*     Giorgio E.M. Melloni, Sabina A. Murphy,  Andrea Bellavia

* Read the data. Rename variables so that the code below can be re-used. 
* t = time-to-event variable, d = observed/censored event time indicator 
* (0 = right-censored, 1 = observed event time), x = continuous variable.
clear all
cd "/Users/anddis/Library/CloudStorage/OneDrive-KarolinskaInstitutet/other/andrea_rcs/"
import excel "data_example.xlsx", clear firstrow

rename status d
rename time t
rename x x
replace x = round(x, 0.01)

* Create the restricted cubic splines and save the knots' location in a matrix
mkspline xrcs = x , nknots(4) cubic displayknots
mat knotslocation = r(knots)


*** Part 1: Model fitting and output interpretation ***
* Fit Cox proportional-hazards model with RCS transforms for variable x
stset t, fail(d) scale(365.24)
stcox xrcs1 xrcs2 xrcs3, nohr
 
* Test for overall association between x and probability of y 
* (joint Wald test on all RCS transforms).
testparm xrcs1 xrcs2 xrcs3

* Test for nonlinearity in the association between x and probability of y 
* (joint Wald test on all RCS transforms except the first one, which coincides with x itself).
testparm xrcs2 xrcs3


*** Part 2: Graphical presentation of exposure-response relationship ***
* Prepare the data to plot the Hazard Ratio for x, using 5 as the referent value (Hazard Ratio = 1).
su x, detail
levelsof x
xbrcspline xrcs, values(`r(levels)') ref(5) ///
	matknots(knotslocation) eform gen(xplot hr hrlb hrub) 
	
* Plot the Hazard Ratio for x, using 5 as the referent value (Hazard Ratio = 1). 
* Overlay an histogram for the distribution of x
twoway (histogram x, yaxis(2) yscale(alt axis(2) r(0 15)) percent color(gs12) bin(30)) ///
(line hrlb hrub hr xplot, lp(- - l) lc(black black black) yaxis(1) yscale(log axis(1) alt)), ///
	scheme(s2mono) legend(off) ///
	ylabel(0.5 1 2 4 6 10 20, angle(horiz) format(%3.2fc)) ///
	ylabel(0(5)15, angle(horiz) format(%2.0fc) axis(2)) ///
	xlabel(1/9) ///
	ytitle("Hazard Ratio (95% confidence interval)") ///
	ytitle("Percent of subjects", axis(2) placement(center)) ///
	xtitle("log(NT-proBNP)") name(graph1, replace)
graph export "stata/survival_hazardratio.pdf", replace

* Predict event probabilities at 1, 1.5, 2, 3, 3.5 years as a function of x.
predict loghr, xb
predict basesurvprob, basesurv

su basesurvprob if _d==1 & _t<1 // 1 year
gen survprob1 = 1-`r(min)'^exp(loghr) 
label var survprob1 "1 year"
su basesurvprob if _d==1 & _t<1.5 // 1.5 years
gen survprob15 = 1-`r(min)'^exp(loghr)
label var survprob15 "1.5 years"
su basesurvprob if _d==1 & _t<2 // 2 years
gen survprob2 = 1-`r(min)'^exp(loghr) 
label var survprob2 "2 years"
su basesurvprob if _d==1 & _t<3 // 3 years
gen survprob3 = 1-`r(min)'^exp(loghr)
label var survprob3 "3 years"
su basesurvprob if _d==1 & _t<3.5 // 3.5 years
gen survprob35 = 1-`r(min)'^exp(loghr) 
label var survprob35 "3.5 years"

* Plot event probabilities at 1, 1.5, 2, 3, 3.5 years as a function of x.
twoway (histogram x, yaxis(2) yscale(alt axis(2) r(0 15)) percent color(gs12) bin(30)) ///
(line survprob1 survprob15 survprob2 survprob3 survprob35 x, sort lc(black) yaxis(1) yscale(axis(1) alt)), ///
	scheme(s2mono) legend(rows(2)) ///
	ylabel(0(0.1)1, angle(horiz) format(%3.1fc)) ///
	ylabel(0(5)15, angle(horiz) format(%2.0fc) axis(2)) ///
	xlabel(1/9) ///
	ytitle("Event probability (95% confidence interval)") ///
	ytitle("Percent of subjects", axis(2) placement(center)) ///
	xtitle("log(NT-proBNP)") name(graph2, replace)
graph export "stata/survival_probability.pdf", replace

* Assess how the dose-response association changes as the number of the knots changes.
drop xrcs1 xrcs2 xrcs3 xplot hr hrlb hrub loghr basesurvprob survprob1 survprob15 survprob2 survprob3 survprob35

quietly forvalues k = 3/7 {
	noisily di "---> `k' knots" _n
	cap drop xplot
	mkspline xrcs = x , nknots(`k') cubic displayknots
	mat knotslocation = r(knots)
	stcox xrcs?
	noisily estat ic
	levelsof x
	xbrcspline xrcs, values(`r(levels)') ref(5) ///
		matknots(knotslocation) eform gen(xplot hr`k' hrlb hrub)
	label variable hr`k' "`k' knots"
	drop hrlb hrub xrcs?
}

twoway (histogram x, yaxis(2) yscale(alt axis(2) r(0 15)) percent color(gs12) bin(30)) ///
(line hr3 hr4 hr5 hr6 hr7 xplot, lc(black) yaxis(1) yscale(log axis(1) alt)), ///
	scheme(s2mono) legend(rows(2)) ///
	ylabel(1 2 4 6 10 20, angle(horiz) format(%3.2fc)) ///
	ylabel(0(5)15, angle(horiz) format(%2.0fc) axis(2)) ///
	xlabel(1/9) ///
	ytitle("Hazard Ratio (95% confidence interval)") ///
	ytitle("Percent of subjects", axis(2) placement(center)) ///
	xtitle("log(NT-proBNP)") name(graph3, replace)
graph export "stata/survival_hazardratio_knots.pdf", replace



