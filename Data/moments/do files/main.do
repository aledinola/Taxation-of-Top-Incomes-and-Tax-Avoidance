** Purpose: Analyse SCF data to get distribution of LFO by wealth and income
*  Last updated by: Alessandro and Haomin July 1, 2021

clear
clear matrix
set more off

// --------------------- SET DIRECTORIES -------------------------------------//

* Change mydir
*global mydir "C:\Users\haomi\Dropbox\6_TaxAvoidance\data"
global mydir "C:\Users\aledi\Dropbox\1 - RESEARCH PROJECTS\2_ProjectTaxAvoidance\PRIVATE\Data\SCF"
global	datadir "${mydir}\data" /* raw data*/
global	savedir "${mydir}\work" /* working directory*/
global	resdir "${mydir}\results" /* results*/

global userestr = 1  //1-use data cleaning (age and female HH)



********************************************
* Table of Contents:
* 1. Data extraction, cleaning and sample restriction
*	1.1 SCF 2013
*	1.2 SUSB
* 2. Compute moments from SCF
*	2.1 Population
*		> share of ABO in the population
*		> share of ABO in the top 1%, 5%, 10% income, wealth
*		> Share of ABO and LFO by wealth quintiles
*		> share of income and wealth held by ABO
*		> Income and wealth inequality (Gini, 90-10 ratio, etc.)

*	2.2 Entrepreneurs
*		> Share of S- and C- corps
*		> Inequality of income and wealth among ABOs
*		> Employment by LFO
* 		> Emp share of businesses by size quintiles
*		> leverage by LFO
*	2.3 TAXSIM: Median dividend tax rate of C-corp owners

* 3. Compute firm size moments from SUSB
* 4. Approximate income tax function
****************************************************************************************


********************************************
* 1. Data extraction, cleaning and sample restriction 
*	1.1 SCF 2013
* codebook: https://sda.berkeley.edu/sdaweb/docs/scfcomb2019/DOC/hcbk.htm
********************************************

* Extract business-related variables from the raw data, then merge with public summary file

use	YY1 Y1 X3103 X3104 X3107 X3119 X3111 X3113 X3114 X3121 X3131 X4131 X4731 using "$datadir/p13i6.dta", clear


rename X3103 privbus /* own private business: 1 = yes, 5 = no */
rename X3104 actbus /* active business owner: 1 = yes, 5 = no, 0 = no businesses  */
rename X3107 ind /* CENSUS 2013 4-DIGIT INDUSTRY CODE */
rename X3119 lfo
rename X3111 size /* number of paid and unpaied workers, including oneself */
rename X3113 workinbus_r /* 1 = respondent works in own business */
rename X3114 workinbus_s /* 1 = spuose of respondent works in own business*/
rename X3121 busdebt /* Business debt*/
rename X3131 sales /* sales of business */
rename X4131 businc1 /*business income*/
rename X4731 businc2 /*business income*/

replace businc1 = 0 if businc1<0
replace businc2 = 0 if businc2<0
gen businc = businc1 + businc2

label define lfo 1 "Partnership" 2 "Sole Prop." 3 "S-Corp." 4 "C-Corp. and other corp." ///
	6 "Foreign" 11 "LLP" 12 "LLC" 15 "Cooperative" 40 "Not a formal bus." -7 "Other" 0 "No Business"
label value lfo lfo

label define ind 1 "Agriculture" 2 "Mining, Utilities, Construction" 3 "Manufacturing"  ///
	4 "Wholesale and Retail" 5 "Information, finance, Professional Services" ///
	6 "Educational and Health Care Services" 7 "Arts, Entertainment" ///
	8 "Other Services except Public Admin" 9 "Public Admin" 
label value ind ind

* Merge with public summary file using identifiers Y1 and YY1
ren Y1 y1 
ren YY1 yy1
merge 1:1 yy1 y1 using ${datadir}/rscfp2013.dta
drop _merge

* Data restriction
if $userestr==1 { 
	drop if hhsex==2               // keep only families where the HEAD is male 
	keep if age > 24 & age < 65 // keep individuals of working age
	*keep if occat1==1 | occat1 == 2 // keep only those who are working (for oneself or others)
}


* Create additional variables

* entre = indicator for active business owners
gen entre = (actbus == 1) 
label define entre 0 "worker" 1 "entre"
label value entre entre

* redefine lfo: partnership and sole-prop become one category
ren lfo lfo_scf
recode lfo_scf (1 2 = 1) (3 = 2) (4 = 3),gen(lfo)
replace lfo = . if lfo<1 | lfo>3
label define lfo3 1 "Sole-Prop./Partn." 2 "S-Corp." 3 "C-Corp."
label value lfo lfo3

* Create quintiles of networth
xtile nw_qt = networth [aw = wgt], nq(5)
label define nw_qt 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5"
label value nw_qt nw_qt

* Create quintiles of income
xtile income_qt = income [aw = wgt], nq(5)
label define income_qt 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5"
label value income_qt income_qt

* Create size bins
gen size_bin =.
replace size_bin = 1 if entre==1 & size>0 & size<=4
replace size_bin = 2 if entre==1 & size>4 & size<=9
replace size_bin = 3 if entre==1 & size>9 & size<=19
replace size_bin = 4 if entre==1 & size>19 & size<=99
replace size_bin = 5 if entre==1 & size>99 & size<.
label define size_bin 1 "0-4" 2 "5-9"  3 "10-19" 4 "20-99" 5 "100+"
label value size_bin size_bin
* Create indicator for wether ABO is an employer
* Definition of employer: total # working in the bus., excluding respondent and spouse.
replace workinbus_r = 0 if workinbus_r~=1
replace workinbus_s = 0 if workinbus_s~=1
gen workinbus_total = workinbus_r + workinbus_s
gen employer = .
replace employer = 1 if entre == 1 & size-workinbus_total>0 & size<.
replace employer = 0 if entre ==1 & employer~=1
* Create size excluding respondent and spouse
gen size_nf = size-workinbus_total
label var size_nf "Size excl. R and spouse"
compress

save ${savedir}/scf2013_clean.dta,replace


********************************************
*	1.2 SUSB
* Excel tables are downloaded here: 
* https://www.census.gov/programs-surveys/susb/data/tables.html
********************************************

			
* A.1 Read raw SUSB tables by year and combine them into one file
*2010
import excel using ${SUSB}/us_naicssector_lfo_2010.xls, ///
			cellrange(A8:K1380)  clear
save ${savedir}/SUSB_empdist_2010.dta,replace
* 2011
import excel using ${SUSB}/us_naicssector_lfo_2011.xls, ///
			cellrange(A8:K1377)  clear
save ${savedir}/SUSB_empdist_2011.dta,replace
* 2012
import excel using ${SUSB}/us_naicssector_lfo_2012.xls, ///
			cellrange(A9:M1378)  clear
save ${savedir}/SUSB_empdist_2012.dta,replace
* 2013
import excel using ${SUSB}/us_naicssector_lfo_2013.xlsx, ///
			cellrange(A8:K1374)  clear
save ${savedir}/SUSB_empdist_2013.dta,replace
* 2014
import excel using ${SUSB}/us_naicssector_lfo_2014.xlsx, ///
			cellrange(A9:K1383)  clear
save ${savedir}/SUSB_empdist_2014.dta,replace
* 2015
import excel using ${SUSB}/us_naicssector_lfo_2015.xlsx, ///
			cellrange(A9:K1379)  clear
save ${savedir}/SUSB_empdist_2015.dta,replace
* 2016
import excel using ${SUSB}/us_naicssector_lfo_2016.xlsx, ///
			cellrange(A10:K1380)  clear
save ${savedir}/SUSB_empdist_2016.dta,replace
* 2017
import excel using ${SUSB}/us_naicssector_lfo_2017.xlsx, ///
			cellrange(A10:K1355)  clear
save ${savedir}/SUSB_empdist_2017.dta,replace

* Combine data from all years
forvalues year = 2010/2017 {
		use ${savedir}/SUSB_empdist_`year'.dta,clear
		* Gen year
		display "year `year'"
		gen year = `year'

		* Save file
		if (`year'>2010) {
			append using "${work}/SUSB_empdist.dta"
		}
		compress
		save ${savedir}/SUSB_empdist.dta,replace
		rm ${savedir}/SUSB_empdist_`year'.dta
}

* A.2 Change variable names and add labels, create indicators for covid-impacted sectors
use ${savedir}/SUSB_empdist.dta,clear
* Gen lfo: legal form
gen lfo = substr(A,1,1)
destring lfo, replace
drop A
* NAICS code
gen naics = substr(B,1,2)
destring naics,replace force
replace naics = 0 if naics ==.
drop B C
* Employment size categories
gen emp_c = substr(D,1,1)
destring emp_c,replace
drop D


* Number of firms and establishment
ren E N_firm
ren F N_estab
* Number of total employment
ren G tot_emp
* Annual Payroll
ren J tot_payroll
* Annual Receipts
ren L tot_receipts
* Clean up
drop H I K M
order year lfo naics emp_c N_firm N_estab tot_emp tot_payroll tot_receipts
* Recode lfo
recode lfo (2=3) (3=2) (4 5 =1) (else = .),gen(lfo3)
* Labels
label define legal_form 1 "Total" 2 "Corporation" 3 "S-Corp" 4 "Partnership" ///
			5 "Sole Proprietorship" 6 "Non-Profit" 7 "Government" 8 "Other"
label define legal_form3 1 "Sole Prop/Partn" 2 "S-Corp" 3 "C-Corp" 
label define naics 0 "All" 11 "Agriculture, forestry, fishing and hunting" ///
	21 "Mining, quarrying, and oil and gas extraction" 22 "Utilities" ///
	23 "Construction" 31 "Manufacturing" 42 "Wholesale trade" 44 "Retail trade" ///
	48 "Transportation and warehousing" 51 "Information" 52 "Finance and insurance" ///
	53 "Real estate and rental and leasing" 54 "Professional, scientific, and technical services" ///
	55 "Management of companies and enterprises" 56 "Administrative and support and waste management and remediation services" ///
	61 "Educational services" 62 "Health care and social assistance" 71 "Arts, entertainment, and recreation" /// 
	72 "Accommodation and food services" 81 "Other services (except public administration)" ///
	99 "Industries not classified"

label define emp_c 1 "Total" 2 "0-4" 3 "5-9" 4 "10-19" 5 "<20" 6 "20-99" 7 "100-499" 8 "<500" 9 "500+"

label var lfo "Legal Form of Org."
label value lfo legal_form
label value lfo3 legal_form3
label var lfo3 "Legal Form of Org."
label var naics "NAICS Code"
label value naics naics
label var emp_c "Employment Size Category"
label value emp_c emp_c
label var N_firm "Number of Firms"
label var N_estab "Number of Establishment"
label var tot_emp "Number of total employment"
label var tot_payroll "Annual Payroll ($1,000)"
label var tot_receipts "Annual Receipts ($1,000)"



* Drop emp categories: total,  <20 , <500
drop if inlist(emp_c,1,5,8)



compress
save ${savedir}/SUSB_empdist.dta,replace


********************************************
* 2. Compute moments from SCF
*	2.1 Population
*		> share of ABO in the population
*		> share of ABO in the top 1%, 5%, 10% income, wealth
*		> Share of ABO and LFO by wealth quintiles
*		> Share of ABO and LFO by income quintiles
*		> share of income and wealth held by ABO
*		> Income and wealth inequality (Gini, 90-10 ratio, etc.)
*	2.2 Entrepreneurs
*		> Share of S- and C- corps
*		> Inequality of income and wealth among ABOs
*		> Employment by LFO
* 		> Emp share of businesses by size quintiles
*		> leverage by LFO
********************************************



********************************************
*	2.1 Population
********************************************

* Create log file (replacing existing file)
cap log close                                           
log using ${resdir}/data_population.txt, text replace nomsg 
log off

* Load data
use ${savedir}/scf2013_clean.dta,clear

*		> share of ABO in the population
log on
disp "Fraction of ABO in Population"
tab entre [aw = wgt]
log off

*		> Share of ABO and LFO by wealth quintiles
log on
disp "Share of entrepreneurs by wealth quintiles"
tab nw_qt entre [aw = wgt],row nofreq
log off

* Latex code
latab nw_qt entre [aw = wgt], row dec(3) tf("${resdir}/entre_networth") replace

preserve 
	collapse (sum) wgt,by(nw_qt entre) 
	bysort nw_qt: egen wgt_tot = total(wgt)
	gen share_entre = wgt/wgt_tot
	drop wgt
	reshape wide share_entre,i(nw_qt) j(entre)
	graph bar share_entre0 share_entre1, over(nw_qt) stack ///
		legend(order(1 "Worker" 2 "Entre.")) name(share_entre,replace)
	graph export ${resdir}/share_entre.png, name(share_entre) replace
restore

log on
disp "Share of LFO by wealth quintiles"
tab nw_qt lfo [aw = wgt] if lfo~=.,row nofreq
log off

* Latex code:
latab nw_qt lfo [aw = wgt] , row dec(3) tf("${resdir}/lfo_networth") replace 
preserve 
	drop if lfo==.
	collapse (sum) wgt ,by(nw_qt lfo) 
	bysort nw_qt: egen wgt_tot = total(wgt)
	gen share_lfo = wgt/wgt_tot
	drop wgt
	reshape wide share_lfo,i(nw_qt) j(lfo)
	graph bar share_lfo1 share_lfo2 share_lfo3, over(nw_qt) stack ///
		legend(order(1 "Sole Prop./Partn." 2 "S-Corp" 3 "C-Corp")) name(share_lfo,replace) ///
		bar(1, color(ebblue)) bar(2, color(sienna)) bar(3, color(gold)) 
	graph export ${resdir}/share_lfo_networth.png, name(share_lfo) replace
restore

***
*		> Share of ABO and LFO by income quintiles
log on
disp "Share of entrepreneurs by income quintiles"
tab income_qt entre [aw = wgt],row nofreq
log off

* Latex code
latab income_qt entre [aw = wgt], row dec(3) tf("${resdir}/entre_income") replace

preserve 
	collapse (sum) wgt,by(income_qt entre) 
	bysort income_qt: egen wgt_tot = total(wgt)
	gen share_entre = wgt/wgt_tot
	drop wgt
	reshape wide share_entre,i(income_qt) j(entre)
	graph bar share_entre0 share_entre1, over(income_qt) stack ///
		legend(order(1 "Worker" 2 "Entre.")) name(share_entre,replace)
	graph export ${resdir}/share_entre_income.png, name(share_entre) replace
restore

log on
disp "Share of LFO by income quintiles"
tab income_qt lfo [aw = wgt] if lfo~=.,row nofreq
log off

* Latex code:
latab income_qt lfo [aw = wgt] , row dec(3) tf("${resdir}/lfo_income") replace 
preserve 
	drop if lfo==.
	collapse (sum) wgt ,by(income_qt lfo) 
	bysort income_qt: egen wgt_tot = total(wgt)
	gen share_lfo = wgt/wgt_tot
	drop wgt
	reshape wide share_lfo,i(income_qt) j(lfo)
	graph bar share_lfo1 share_lfo2 share_lfo3, over(income_qt) stack ///
		legend(order(1 "Sole Prop./Partn." 2 "S-Corp" 3 "C-Corp")) name(share_lfo,replace) ///
		bar(1, color(ebblue)) bar(2, color(sienna)) bar(3, color(gold)) 
	graph export ${resdir}/share_lfo_income.png, name(share_lfo) replace
restore

***


*		> share of income and wealth by ABO status
foreach var in income networth {
	capture drop `var'_all
	egen `var'_all = total(`var'*wgt)
	capture drop `var'_entre
	bys entre: egen `var'_entre = total(`var'*wgt)
	capture drop `var'_share
	gen `var'_share = `var'_entre/`var'_all
}
log on
disp "share of income and wealth by ABO status"
tabstat income_share networth_share, s(mean) by(entre) 
log off


*		> share of income and wealth by LFO
preserve
drop if lfo==. 
*drop if networth <0 | networth>.
foreach var in income networth {
	capture drop `var'_all
	egen `var'_all = total(`var'*wgt)
	capture drop `var'_lfo
	bys lfo: egen `var'_lfo = total(`var'*wgt)
	capture drop `var'_share
	gen `var'_share = `var'_lfo/`var'_all
}
*log on
disp "share of income and wealth by LFO status"
tabstat income_share networth_share, s(mean) by(lfo) 
*log off
restore

*		> share of ABO in the top 1%, 5%, 10% income, wealth
capture drop networth_cdf
cumul networth [aw = wgt],gen(networth_cdf)
capture drop income_cdf
cumul income [aw = wgt],gen(income_cdf)

foreach var in networth income {
foreach x in 1 5 10 {
	
	log on
	display "Share of ABO in the top `x'% by `var'"
	tab entre [aw = wgt] if `var'_cdf>=1-`x'/100
	display "Share of lfo in the top `x'% by `var'"
	tab lfo [aw = wgt] if `var'_cdf>=1-`x'/100
	display " "
	log off
}
}

*		> Income and wealth inequality

* Gini, 90-10 ratio, etc.
log on
disp "Income Inequality"
 ineqdec0  income [aw = wgt]
disp "Wealth (networth) Inequality"
ineqdec0  networth [aw = wgt]
log off

* share of wealth and income in the top 1,10,20 and bottom 40

log on
quietly foreach var in income networth {
	noisily display "`var'"
	capture drop `var'_total
	egen `var'_total = total(`var'*wgt)
	* top 1,10,20
	foreach x in 1 10 20 {
		capture drop top`x'
		gen top`x'= (`var'_cdf>=1-`x'/100)
		capture drop temp
		gen temp = `var'*wgt*top`x'
		capture drop `var'_top`x' 
		egen `var'_top`x' = total(temp) 
		capture drop `var'_top`x'_share
		gen `var'_top`x'_share = `var'_top`x'/`var'_total if top`x'==1
		sum `var'_top`x'_share
		local `var'_top`x'_share = r(mean)
		noisily display "top `x'% share =  ``var'_top`x'_share' "
	}
	* bottom 40
	foreach x in 40 {
		capture drop bottom`x'
		gen bottom`x'= (`var'_cdf<`x'/100)
		capture drop temp
		gen temp = `var'*wgt*bottom`x'
		capture drop `var'_bottom`x' 
		egen `var'_bottom`x' = total(temp) 
		capture drop `var'_bottom`x'_share
		gen `var'_bottom`x'_share = `var'_bottom`x'/`var'_total if bottom`x'==1
		sum `var'_bottom`x'_share
		local `var'_bottom`x'_share = r(mean)
		noisily display "bottom `x'% share =  ``var'_bottom`x'_share' "
	}	
}
log off


*		> Other moments
* Ratio of mean income (entr to worker)

sum income [aw = wgt] if entre==1,d
local mean_inc_entre = r(mean)
local med_inc_entre = r(p50)

sum income [aw = wgt] if entre==0,d
local mean_inc_work = r(mean)
local med_inc_work = r(p50)
disp "median entre `med_inc_entre' worker `med_inc_work'" 

local ratiomeaninc = `mean_inc_entre'/`mean_inc_work'
local ratiomedinc = `med_inc_entre'/`med_inc_work'

log on 
disp "Ratio of mean income (entr to worker) = `ratiomeaninc' " 
disp "Ratio of median income (entr to worker) = `ratiomedinc' " 
log off

log close




****************************
*	2.2 Entrepreneurs
****************************

*		> total asset to net worth ratio by LFO
cap log close                                           
log using ${resdir}/data_entre.txt, text replace nomsg 
log off

* Load data
use ${savedir}/scf2013_clean.dta,clear
* keep only businesses: lfo = S-Corp, sole-prop/partnership, C-Corp
keep if lfo ~=.


*		> Share of S- and C- corps
log on
tab lfo [aw = wgt]
log off

*		> Inequality of income and wealth among ABOs
log on
disp "Income Inequality of ABOs"
 ineqdec0  income [aw = wgt]
disp "Wealth (networth) Inequality of ABOs"
ineqdec0  networth [aw = wgt]
log off

*		> Employment by LFO
* Share of ABOs that are employers (total and by lfo)
log on
disp "Share of ABOs that are employers"
tab lfo employer [aw = wgt],row
log off
* Emp share by LFO
preserve
	replace size = . if size<0
	replace size = size*wgt
	collapse (sum) size,by(lfo)
	egen size_total = total(size)
	gen share_lfo = size/size_total
	forvalues x = 1/3 {
		local share_`x' = share_lfo[`x']
	}
	log on
	display "Emp share of sole prop./partn. = `share_1'"
	display "Emp share of S-Corps = `share_2'"
	display "Emp share of C-Corps = `share_3'"
	log off
restore



* 		> Emp share of businesses by size quintiles
xtile emp_qt = size [aw = wgt] , nq(5)
preserve
	replace size = size*wgt
	collapse (sum) size,by(emp_qt) 
	egen size_total = total(size)
	gen share_qt = size/size_total
	log on
	quietly forvalues q = 1/5 {
		local share = share_qt[`q']
		noisily display "Emp share of size quintile `q' = `share'"
	}
	display "	"
	log off
restore

*		> leverage by LFO

*replace levratio = 0.5 if levratio>0.5 & levratio<.

log on 
display "levratio = hh debt to hh asset ratio"
tabstat levratio [aw = wgt] if levratio<=1,s(mean median) by(lfo)
log off


* aggregate business debt to aggregate net worth, by LFO
* aggregate business debt to aggregate hh asset, by LFO
* total asset to net worth ratio by LFO

replace busdebt = . if busdebt<0
preserve
	foreach var in busdebt networth asset {
		replace `var' = `var'*wgt
	}
	collapse (sum) busdebt (sum) networth (sum) asset,by(lfo)
	gen busdebt_networth_ratio = busdebt/networth
	gen busdebt_asset_ratio = busdebt/asset
	gen asset_networth_ratio = asset/networth
	log on
	display "Business debt to household net worth ratio, Business debt to household asset ratio, Household asset to net worth ratio"
	list lfo busdebt_networth_ratio busdebt_asset_ratio asset_networth_ratio,abbreviate(20)
	log off
restore



log close


*	2.3 TAXSIM: 
* 	- Median dividend tax rate of C-corp owners
*	- Average tax rate by income distribution
* 	REF: see https://taxsim.nber.org/to-taxsim/scf27-32/
* 		also see https://taxsim.nber.org/taxsim35/
* Install micombine
net sj 7-4 st0067_3
net install st0067_3
* Install TAXSIM32
net from "https://taxsim.nber.org/stata"
net describe taxsim35
net install taxsim35

* Increase maxvar
clear
set maxvar 8000

* Load SCF with tax sim data (aggregated to households) 
use ${datadir}/scftax13.dta,clear

generate weight=x42001/5
generate rep = mod(y1,10)
generate taxunit = mod(taxsimid,10)+1

taxsim35, replace 
* marginal dividend tax rate

gen drate = .
replace drate = 0 if frate<0.25
replace drate = 0.188 if frate>=0.25 & frate<=0.35
replace drate = 0.238 if frate>0.35 & frate<.

keep if taxunit ==1
keep y1 yy1 drate
duplicates drop

save ${savedir}/scf2013_drate.dta,replace

* Merge scf_drate with scf2013_clean
use ${savedir}/scf2013_clean.dta,clear
merge 1:1 y1 yy1 using ${savedir}/scf2013_drate.dta
drop if _merge == 2
drop _merge 
* Keep only c-corporation owners
keep if lfo==3
* Compute weighted average marginal dividend tax
gen intdivinc_drate = intdivinc*drate*wgt
gen intdivinc_w = intdivinc*wgt
egen tot_intdivinc_drate = total(intdivinc_drate)
egen tot_intdivinc = total(intdivinc_w)
gen tau_d = tot_intdivinc_drate/tot_intdivinc
sum tau_d
/* tau_d is 0.1806 */

*	- Average tax rate by income distribution
use ${datadir}/scftax13.dta,clear

generate weight=x42001/5
generate rep = mod(y1,10)
generate taxunit = mod(taxsimid,10)+1

taxsim35, replace 

keep if taxunit == 1
keep weight v10 v18 fiitax y1 yy1
compress
save ${savedir}/scf2013_atr.dta,replace

* Pre-tax income: v18 (federal taxable income)
* Federal income tax liability: fiitax
* Cumulative income distribution
use "${savedir}/scf2013_atr.dta",clear
merge 1:1 y1 yy1 using "${savedir}/scf2013_clean.dta"
keep if _merge==3
drop _merge
replace v18 = 0 if v18<0
sort v18
gen cdf_inc = sum(weight)
replace cdf_inc = cdf_inc/cdf_inc[_N]

gen inc_group = .
replace inc_group = 1 if cdf_inc>=0.999 & cdf_inc<. /* Top 0.1%*/
replace inc_group = 2 if cdf_inc>=0.99 & cdf_inc<0.999/* p99-p99.9*/
replace inc_group = 3 if cdf_inc>=0.9 & cdf_inc<0.99 /* P90-P99%*/
replace inc_group = 4 if cdf_inc>=0.5 & cdf_inc<0.9 /* P50-P90%*/
replace inc_group = 5 if cdf_inc>=0 & cdf_inc<0.5 /* Bottom 50% %*/

gen pretaxinc_w = v18*weight
gen tax_w = fiitax*weight

bysort inc_group: egen sum_retaxinc = sum(pretaxinc_w)
bysort inc_group: egen sum_tax = sum(tax_w)
gen atr = sum_tax/sum_retaxinc
keep atr inc_group
duplicates drop 
browse
/* Average tax rate by income group
inc_group	atr
Top 0.1		.2746097
P99-P99.9	.3062719
P90-P99		.2529508
P50-P90		.1603018
Bottom 50%	.0638568
*/

/*

//------------------Reg passthrough on networth, and log sales ------------------------------//
* Following Dyrda and Pugsley (2019), we use log-sales as proxy for entrepreneurial ability
gen passthrough = (lfo==1 | lfo==2) /*pass-through businesses*/

* we consider both total sales and sales per worker
capture drop logsales*
gen logSales = log(sales) if sales>0
gen logSales_worker = log(sales/size) if sales>0 & size>0
capture drop networth_cdf
cumul networth [aw = wgt],gen(networth_cdf)
log on
*probit passthrough  networth_cdf [pw = wgt]
probit passthrough logSales networth_cdf [pw = wgt]
*probit passthrough logSales_worker networth_cdf [pw = wgt]
log off


//---------------------- Dist of LFO by employment -------------------------//


* distribution of lfo by wealth quintile
log on 
tab size_bin lfo [aw = wgt] ,row nofreq
log off
* Latex code:
latab size_bin lfo [aw = wgt] , row dec(3) tf("${resdir}/lfo_size") replace 

preserve 
	collapse (sum) wgt ,by(size_bin lfo) 
	bysort size_bin: egen wgt_tot = total(wgt)
	gen share_lfo = wgt/wgt_tot
	drop wgt
	reshape wide share_lfo,i(size_bin) j(lfo)
	graph bar share_lfo1 share_lfo2 share_lfo3 share_lfo4, over(size_bin) stack ///
		legend(order(1 "Partnership" 2 "Sole Prop." 3 "S-Corp" 4 "C-Corp")) name(share_lfo,replace) ///
		bar(1, color(ebblue)) bar(2, color(sienna)) bar(3, color(gold)) bar(4, color(green))
	graph export ${resdir}/share_lfo_size.png, name(share_lfo) replace
restore

//---------------------- Dist of LFO by log-Sales -------------------------//
* Create quintiles of networth
capture drop sales_qt
xtile sales_qt = sales [aw = wgt] if sales>=0, nq(5)
label define sales_qt 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5",replace
label value sales_qt sales_qt

* distribution of lfo by wealth quintile
log on 
tab sales_qt lfo [aw = wgt] ,row nofreq
log off
* Latex code:
latab sales_qt lfo [aw = wgt] , row dec(3) tf("${resdir}/lfo_sales") replace 

preserve 
	collapse (sum) wgt ,by(sales_qt lfo) 
	bysort sales_qt: egen wgt_tot = total(wgt)
	gen share_lfo = wgt/wgt_tot
	drop wgt
	reshape wide share_lfo,i(sales_qt) j(lfo)
	graph bar share_lfo1 share_lfo2 share_lfo3 share_lfo4, over(sales_qt) stack ///
		legend(order(1 "Partnership" 2 "Sole Prop." 3 "S-Corp" 4 "C-Corp")) name(share_lfo,replace) ///
		bar(1, color(ebblue)) bar(2, color(sienna)) bar(3, color(gold)) bar(4, color(green))
	graph export ${resdir}/share_lfo_sales.png, name(share_lfo) replace
restore


log close

*Drop if not a worker or a self-employed household
drop if occat1>2
/*1. work for someone else 
2. self-employed/partnership
3. retired/disabled + (student/homemaker/misc. not working and age 65 or older)
4. other groups not working (mainly those under 65 and out of the labor force)*/

gen bus_own = 0
replace bus_own = 1 if actbus>0 

gen semp = 0
replace semp = 1 if occat1==2 

tab bus_own [aw = wgt]
tab semp [aw = wgt]

gen bus_all = 0

replace bus_all = 1 if bus_own== 1 | semp == 1

tab bus_all [aw = wgt]
*/

********************************************
* 3. Compute firm size moments from SUSB
********************************************

* Create log file (replacing existing file)
cap log close                                           
log using ${resdir}/SUSB_empsize.txt, text replace nomsg 
log off

use ${savedir}/SUSB_empdist.dta,clear

* keep only Sole prop, S and C corps
keep if lfo3 <.
* drop rows such that emp category is all
keep if naics ==0
* Keep 2013
keep if year == 2013


* average firm size by LFO, including large firms (500+)
bysort lfo3: egen emp = total(tot_emp)
bysort lfo3: egen firms = total(N_firm)
bysort lfo3: egen payroll = total(tot_payroll)

gen empsize = emp/firms
gen payrollsize = payroll/firms

* Adjust for non-employers using fraction of employer by LFO computed above from SCF
* Prop/Partnership: 0.4349
* S-Corps: 0.6609
* C-Corps: 0.699

replace empsize = empsize * 0.4349 if lfo3==1
replace empsize = empsize * 0.6609 if lfo3==2
replace empsize = empsize * 0.699 if lfo3==3
replace payrollsize = payrollsize * 0.4349 if lfo3==1
replace payrollsize = payrollsize * 0.6609 if lfo3==2
replace payrollsize = payrollsize * 0.699 if lfo3==3



log on 
display "Firm size by LFO"
quietly {
sum empsize if lfo3 == 1
local s_1 = r(mean)
noisily display "Average firm size of Sole Prop/Partn = `s_1'"
sum empsize if lfo3 == 2
local s_2 = r(mean)
noisily display "Average firm size of S-Corps = `s_2'"
sum empsize if lfo3 == 3
local s_3 = r(mean)
noisily display "Average firm size of C-Corps = `s_3'"
local s_2_1 = `s_2'/`s_1'
noisily display "Average firm size of S-Corps rel to Sole-Prop = `s_2_1'"
local s_3_1 = `s_3'/`s_1'
noisily display "Average firm size of C-Corps rel to Sole-Prop = `s_3_1'"
}
log off

log on 
display "Payroll size by LFO"
quietly {
sum payrollsize if lfo3 == 1
local s_1 = r(mean)
noisily display "Average payroll size of Sole Prop/Partn = `s_1'"
sum payrollsize if lfo3 == 2
local s_2 = r(mean)
noisily display "Average payroll size of S-Corps = `s_2'"
sum payrollsize if lfo3 == 3
local s_3 = r(mean)
noisily display "Average payroll size of C-Corps = `s_3'"
local s_2_1 = `s_2'/`s_1'
noisily display "Average payroll size of S-Corps rel to Sole-Prop = `s_2_1'"
local s_3_1 = `s_3'/`s_1'
noisily display "Average payroll size of C-Corps rel to Sole-Prop = `s_3_1'"
}
log off



log close
********************************************
* 4. Approximate income tax function
*	4.a share of population paying the top marginal tax rate
********************************************


* Generate net income according to the 2013 tax schedule
* Income brackets
global incbrac1 = 17850
global incbrac2 = 72500
global incbrac3 = 146400
global incbrac4 = 223050
global incbrac5 = 398350
global incbrac6 = 450000
* Marginal tax rates
global margtax1 = 0.1
global margtax2 = 0.15
global margtax3 = 0.25
global margtax4 = 0.28
global margtax5 = 0.33
global margtax6 = 0.35
global margtax7 = 0.396
* Standard tax deduction for married households filing jointly
global taxdeduct = 12200
* Mean household income 
global meaninc = 86620.32
* Income brackets: multiples of mean
global multbrac1 = 0.206
global multbrac2 = 0.837
global multbrac3 = 1.690
global multbrac4 = 2.575
global multbrac5 = 4.599
global multbrac6 = 5.195



* Load data
use ${savedir}/scf2013_clean.dta,clear

* Sample selection:
* Married, age 25-64
keep if age > 24 & age < 65 & married ==1 


* gross income = income - ssc 
gen grossinc = income - ssretinc 
replace grossinc = 0 if grossinc<0
gen grossinc_deduct = grossinc - ${taxdeduct}
replace grossinc_deduct = 0 if grossinc_deduct<0
* gen incbracX = gross income in each bracket
capture drop incbrac1
gen incbrac1 = grossinc_deduct if grossinc_deduct<=${incbrac1}
replace incbrac1 = ${incbrac1} if grossinc_deduct >${incbrac1} & grossinc_deduct<.
forvalues b = 2/6 {
	local b0 = `b'-1
	capture drop incbrac`b'
	gen incbrac`b' = 0 if grossinc_deduct<=${incbrac`b0'}
	replace incbrac`b' = grossinc_deduct - ${incbrac`b0'} if grossinc_deduct>${incbrac`b0'} & grossinc_deduct<=${incbrac`b'}
	replace incbrac`b' =   ${incbrac`b'} -  ${incbrac`b0'} if grossinc_deduct>${incbrac`b'} & grossinc_deduct <.
}
capture drop incbrac7
gen incbrac7 = 0 if grossinc_deduct<=${incbrac6}
replace incbrac7 = grossinc_deduct - ${incbrac6} if grossinc_deduct >${incbrac6} & grossinc_deduct<.

* Gen inctax: statutory income tax
capture drop inctax
gen inctax = 0 if grossinc_deduct<.
forvalues b = 1/7 {
	replace inctax = inctax + ${margtax`b'}*incbrac`b' if grossinc_deduct<.
}

* gen netinc: net income according to income tax schedule
capture drop netinc
gen netinc = grossinc - inctax

* gen avetax: average tax rate
gen avetax = inctax/grossinc
gen avenet = 1-avetax
* gen grossmult: multiples of mean gross income
gen grossmult = grossinc_deduct/${meaninc}


* Estimate tax schedule

gen lnet = log(netinc)
gen lgross = log(grossinc)
gen lavenet = log(avenet)
gen lgrossmult = log(grossmult)
* Create log file (replacing existing file)
cap log close                                           
log using ${resdir}/data_tax.txt, text replace nomsg 
display "regressing log(net income) on log(gross income)"
reg lnet lgross [aw = wgt],robust cformat(%9.5f)
*twoway (lfit lnet lgross [aw = wgt])(line lgross lgross,sort)
*(scatter lnet lgross,msymbol(x))
display "regressing log(1-ave tax rate) on log(multiples of mean income)"
reg lavenet lgrossmult [aw = wgt],robust cformat(%9.5f)
*twoway (lfit lavenet lgrossmult [aw = wgt])(scatter lavenet lgrossmult,msymbol(x))


*	4.a share of population paying the top marginal tax rate


* Load data

capture drop paytoprate
gen paytoprate = 0
replace paytoprate = 1 if grossinc_deduct >= $incbrac6
tab paytoprate [aw = wgt]

log close
