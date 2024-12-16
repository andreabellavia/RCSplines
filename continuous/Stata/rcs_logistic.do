* Estimating and presenting nonlinear associations with Restricted Cubic Splines
*     Andrea Discacciati, Michael G. Palazzolo, Jeong-Gun Park, 
*     Giorgio E.M. Melloni, Sabina A. Murphy,  Andrea Bellavia

* Read the data. Rename variables so that the code below can be re-used. 
* y = binary outcome, x = continuous variable 
* This example requires the xbrcspline Stata package: ssc install xbrcspline
clear all
import excel "data_example.xlsx", clear firstrow

rename status y
rename x x
replace x = round(x, 0.01)

* Create the restricted cubic splines and save the knots' location in a matrix
mkspline xrcs = x , nknots(4) cubic displayknots
mat knotslocation = r(knots)

* Fit logistic model with RCS transforms for variable x
logit y xrcs1 xrcs2 xrcs3

* Test for overall association between x and probability of y 
* (joint Wald test on all RCS transforms).
testparm xrcs1 xrcs2 xrcs3

* Test for nonlinearity in the association between x and probability of y 
* (joint Wald test on all RCS transforms except the first one, which coincides with x itself).
testparm xrcs2 xrcs3

* Prepare the data to plot the Odds Ratio for x, using 5 (median) as the referent value (Odds Ratio = 1).
su x, detail
levelsof x
xbrcspline xrcs, values(`r(levels)') ref(5) ///
	matknots(knotslocation) eform gen(xplot or orlb orub) 

* Plot the Odds Ratio for x, using 5 as the referent value (Odds Ratio = 1). 
* Overlay an histogram for the distribution of x
twoway (histogram x, yaxis(2) yscale(alt axis(2) r(0 15)) percent color(gs12) bin(30)) ///
(line orlb orub or xplot, lp(- - l) lc(black black black) yaxis(1) yscale(log axis(1) alt)), ///
	scheme(s2mono) legend(off) ///
	ylabel(0.5 1 2 4 6 10 20, angle(horiz) format(%3.2fc)) ///
	ylabel(0(5)15, angle(horiz) format(%2.0fc) axis(2)) ///
	xlabel(1/9) ///
	ytitle("Odds Ratio (95% confidence interval)") ///
	ytitle("Percent of subjects", axis(2) placement(center)) ///
	xtitle("log(NT-proBNP)") name(graph1, replace)
graph export "logistic_oddsratio.pdf", replace
	
* Predict outcome probabilities as a function of x.
predict ylogit, xb
predict yse, stdp
gen yprob = invlogit(ylogit)
gen yproblb = invlogit(ylogit - invnorm(.975)*yse)
gen yprobub = invlogit(ylogit + invnorm(.975)*yse)

* Plot outcome probabilities as a function of x.
twoway (histogram x, yaxis(2) yscale(alt axis(2) r(0 15)) percent color(gs12) bin(30)) ///
(line yproblb yprobub yprob x, sort lp(- - l) lc(black black black) yaxis(1) yscale(axis(1) alt)), ///
	scheme(s2mono) legend(off) ///
	ylabel(0(0.1)1, angle(horiz) format(%3.1fc)) ///
	ylabel(0(5)15, angle(horiz) format(%2.0fc) axis(2)) ///
	xlabel(1/9) ///
	ytitle("Outcome probability (95% confidence interval)") ///
	ytitle("Percent of subjects", axis(2) placement(center)) ///
	xtitle("log(NT-proBNP)") name(graph2, replace)
graph export "logistic_probability.pdf", replace

drop xrcs1 xrcs2 xrcs3 xplot or orlb orub ylogit yse yprob yproblb yprobub

* Let’s assess how the dose-response association changes as the number of the knots changes.
quietly forvalues k = 3/7 {
	noisily di "---> `k' knots" _n
	cap drop xplot
	mkspline xrcs = x , nknots(`k') cubic displayknots
	mat knotslocation = r(knots)
	logit y xrcs?
	noisily estat ic
	levelsof x
	xbrcspline xrcs, values(`r(levels)') ref(5) ///
		matknots(knotslocation) eform gen(xplot or`k' orlb orub)
	label variable or`k' "`k' knots"
	drop orlb orub xrcs?
}

twoway (histogram x, yaxis(2) yscale(alt axis(2) r(0 15)) percent color(gs12) bin(30)) ///
(line or3 or4 or5 or6 or7 xplot, lc(black) yaxis(1) yscale(log axis(1) alt)), ///
	scheme(s2mono) legend(row(2)) ///
	ylabel(0.25 0.5 1 2 4 6 10 20, angle(horiz) format(%3.2fc)) ///
	ylabel(0(5)15, angle(horiz) format(%2.0fc) axis(2)) ///
	xlabel(1/9) ///
	ytitle("Odds Ratio (95% confidence interval)") ///
	ytitle("Percent of subjects", axis(2) placement(center)) ///
	xtitle("log(NT-proBNP)") name(graph3, replace)
graph export "logistic_oddsratio_knots.pdf", replace
