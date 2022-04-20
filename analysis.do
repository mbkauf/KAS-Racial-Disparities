*--------------------------------------------------
* Header
*--------------------------------------------------
capture log close
clear all 
set more off  
global in_path = ""  // Add file path to where clean data is stored
global out_path = ""  // Add file path to where results are stored
cd ""  // Add working directory
pause on
local date: display %tdCYND date(c(current_date), "DMY")
log using "$in_path/logs/t2tx_`date'.log", text replace 
*--------------------------------------------------

/*
Racial Disparities in Pediatric Kidney Transplants
Goal: Set up data for survival analyses, examine KM curves, and run models
Date Created: 07/17/2020
Last Updated: 08/03/2021

This following is for the time on dialysis to deceased donor transplant equation. For this equation, we have the
    following sections:

    1) Time from activation
      a) Generate any necessary variables
      b) Kaplan-Meier Curves
      c) AFT Models
    2) Time from dialysis  
      a) Kaplan-Meier Curves
      b) AFT Models
    3) Time to Graft Loss
*/
	
*--------------------------------------------------
* 1) Time since activation
*--------------------------------------------------
// Open dataset
use 4yrspeds_final_info_v3.dta, clear

// Generate transplant variable
gen     transplant = (rec_tx_dt == canhx_end_dt) & (rec_tx_dt != .)
replace transplant = 0 if don_ty=="L"  // We only want to count deceased donor transplants

// Genearte KAS variable (Activated before KAS)
gen kas = (first_active_dt >= td(04dec2014))

// Censor events that go beyond 2019 as well as pre-KAS activations into post-KAS
replace transplant = 0 if (rec_tx_dt >= td(04dec2014)) & kas==0
replace transplant = 0 if (rec_tx_dt > td(31dec2019)) & kas==1

// Drop observations that start after cut-offs
drop if kas==0 & canhx_begin_dt >= td(04dec2014)
drop if kas==1 & canhx_begin_dt > td(31dec2019)

// Generate censor date
gen     censor_time = (canhx_end_dt - first_active_dt)
replace censor_time = (td(04dec2014) - first_active_dt) ///
                        if kas == 0 & canhx_end_dt >= td(04dec2014)
replace censor_time = (td(31dec2019) - first_active_dt) ///
                        if censor_time == . & kas == 1                      
replace censor_time = (can_endwlfu - first_active_dt) ///
                        if censor_time == .
replace censor_time = (tfl_lafudate - first_active_dt) ///
                        if censor_time == .

drop if censor_time <= 0  // This is for observations from before activation
drop if censor_time == .  // This is for those who were never activated or immediately inactive
drop if don_ty=="L"  // MK: 08/03/2021 reverting back to dropping LDKT

// Generate activation year
gen active_year = year(first_active_dt)


stset censor_time, failure(transplant==1) id(pers_id) if(canhx_stat_cd != 4999) // exit(time 2000)
 
// Generate race/KAS interaction variables for KM curves
gen     race_kas = race_ethnic
replace race_kas = race_kas + 4 if kas==1

label define race_kas 1 "Pre-KAS: White Non-Hispanic" 2 "Pre-KAS: Black" ///
                      3 "Pre-KAS: Hispanic" 4 "Pre-KAS: Other" ///
                      5 "Post-KAS: White Non-Hispanic" 6 "Post-KAS: Black" ///
                      7 "Post-KAS: Hispanic" 8 "Post-KAS: Other"
label values race_kas race_kas
 
* 1b) KM Curves
// Graph curves
set scheme s2mono 

sts graph, by(race_kas) legend(on size(small) ring(0) pos(1) order(1 "White before KAS" 2 "Black before KAS" ///
    3 "Hispanic before KAS" 4 "Other before KAS" 5 "White after KAS" 6 "Black after KAS"  ///
    7 "Hispanic after KAS" 8 "Other after KAS")) ///
    title("Time Since Activation Until Transplant by Race And Era") xtitle("Months Since Activation") ///
    ytitle("Proportion Awaiting Transplant") tmax(1000)
graph export "$out_path/KM_Race_KAS_Activation.png", replace
 
sts graph if kas==0, by(race_ethnic) xtitle("Days Since Activation") ///
    legend(size(small) ring(0) pos(1) order(1 "White" 2 "Black" 3 "Hispanic" 4 "Other")) ///
    ytitle("Proportion Awaiting Transplant") title("") tmax(1000)
graph export "$out_path/KM_Race_Pre_Activation.png", replace
 
sts graph if kas==1, by(race_ethnic) xtitle("Days Since Activation") ///
    legend(size(small) ring(0) pos(1) order(1 "White" 2 "Black" 3 "Hispanic" 4 "Other")) ///
    ytitle("Proportion Awaiting Transplant") title("")  tmax(1000)
graph export "$out_path/KM_Race_Post_Activation.png", replace


stci, by(race_kas)
 
stsum if kas==0, by(race_ethnic)
stsum if kas==1, by(race_ethnic)
 

* 1c) AFT Models 
*--------------------------------------------------
* Run models by era
*--------------------------------------------------
// Unadjusted (Pre-KAS)
streg i.race_ethnic if kas == 0, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Active.xls", replace ctitle("Unadjusted - Pre-KAS") ///
        excel stats(coef ci pval) eform 

// full model (Pre-KAS) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.race_ethnic i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa active_year  if kas == 0, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Active.xls", append ctitle("Adjusted - Pre-KAS") ///
        excel stats(coef ci pval) eform 
 
// Unadjusted (Post-KAS)
streg i.race_ethnic if kas==1, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Active.xls", append ctitle("Unadjusted - Post-KAS") ///
        excel stats(coef ci pval) eform 

// full model (Post-KAS) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.race_ethnic i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa active_year  if kas == 1, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Active.xls", append ctitle("Adjusted - Post-KAS") ///
        excel stats(coef ci pval) eform 


*--------------------------------------------------
* Run models by race
*--------------------------------------------------
// Unadjusted (White)
streg i.kas if race_ethnic == 1, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Active.xls", replace ctitle("Unadjusted - White") ///
        excel stats(coef ci pval) eform 

// full model (White) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa active_year  if race_ethnic == 1, dist(logl) time tr 
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Adjusted - White") ///
        excel stats(coef ci pval) eform 

// Unadjusted (Black)
streg i.kas if race_ethnic == 2, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Unadjusted - Black") ///
        excel stats(coef ci pval) eform 


// full model (Black) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa active_year  if race_ethnic == 2, dist(logl) time tr 
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Adjusted - Black") ///
        excel stats(coef ci pval) eform 


// Unadjusted (Hispanic)
streg i.kas if race_ethnic == 3, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Unadjusted - Hispanic") ///
        excel stats(coef ci pval) eform 


// full model (Hispanic) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa active_year  if race_ethnic == 3, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Adjusted - Hispanic") ///
        excel stats(coef ci pval) eform 


// Unadjusted (Other)
streg i.kas if race_ethnic == 4, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Unadjusted - Other") ///
        excel stats(coef ci pval) eform 


// full model (Other) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa active_year if race_ethnic == 4, dist(logl) time tr 
outreg2 using "$out_path/Logl_by_Race_Active.xls", append ctitle("Adjusted - Other") ///
        excel stats(coef ci pval) eform 
 
*--------------------------------------------------
* 2) Time since Dialysis
*--------------------------------------------------
// Open dataset
use 4yrspeds_final_info_v2.dta, clear

// Drop those never on dialysis
drop if days_on_dial == 0

// Generate transplant variable
gen transplant = (rec_tx_dt != .) 
    replace transplant = 0 if don_ty=="L"  // We only want to count deceased donor transplants

gen kas = (first_dial_dt >= td(04dec2014))

replace transplant = 0 if (rec_tx_dt >= td(04dec2014)) & kas==0
replace transplant = 0 if (rec_tx_dt > td(31dec2019)) & kas==1

// Generate dialysis year
gen dialysis_year = year(first_dial_dt)

/* 
// Generate censor date
gen censor_time = transplant_dt
  replace censor_time = tfl_lafudate if censor_time == .
  replace censor_time = can_endwlfu if censor_time == .
  replace censor_time = td(04dec2014) if kas == 0 & transplant == 0 

stset censor_time, origin(time first_dial_dt) failure(transplant==1) id(pers_id) 
*/

// Generate diff time variable
gen t2censor = transplant_dt - first_dial_dt
  replace t2censor = td(04dec2014) - first_dial_dt if kas == 0 & transplant == 0
  replace t2censor = td(31dec2019) - first_dial_dt if kas == 1 & transplant == 0
  replace t2censor = tfl_lafudate - first_dial_dt if t2censor==.
  replace t2censor = can_endwlfu - first_dial_dt if t2censor==.

drop if don_ty=="L"  // MK: 08/03/2021 reverting back to dropping LDKT

stset t2censor, failure(transplant==1) id(pers_id) // exit(time 2000)

// Generate race/KAS interaction variables for KM curves
gen race_kas = race_ethnic
    replace race_kas = race_kas + 4 if kas==1

label define race_kas 1 "Pre-KAS: White Non-Hispanic" 2 "Pre-KAS: Black" ///
                      3 "Pre-KAS: Hispanic" 4 "Pre-KAS: Other" ///
                      5 "Post-KAS: White Non-Hispanic" 6 "Post-KAS: Black" ///
                      7 "Post-KAS: Hispanic" 8 "Post-KAS: Other"
label values race_kas race_kas 



* 2b) KM Curves
// Graph curves
sts graph, by(race_kas) legend(on size(small) ring(0) pos(1) order(1 "White before KAS" 2 "Black before KAS" ///
    3 "Hispanic before KAS" 4 "Other before KAS" 5 "White after KAS" 6 "Black after KAS"  ///
    7 "Hispanic after KAS" 8 "Other after KAS")) ///
    title("Time On Dialysis Until Transplant by Race And Era") xtitle("Months On Dialysis") ///
    ytitle("Proportion Awaiting Transplant")  tmax(1000)
graph export "$out_path/KM_Race_KAS_Dialysis.png", replace

sts graph if kas==0, by(race_ethnic) xtitle("Days on Dialysis") ///
    legend(size(small) ring(0) pos(1) order(1 "White" 2 "Black" 3 "Hispanic" 4 "Other")) ///
    ytitle("Proportion Awaiting Transplant") title("")  tmax(1000)
graph export "$out_path/KM_Race_Pre_Dialysis.png", replace

sts graph if kas==1, by(race_ethnic) xtitle("Days on Dialysis") ///
    legend(size(small) ring(0) pos(1) order(1 "White" 2 "Black" 3 "Hispanic" 4 "Other")) ///
    ytitle("Proportion Awaiting Transplant") title("") tmax(1000)    
graph export "$out_path/KM_Race_Post_Dialysis.png", replace


stci, by(race_kas) 


stsum if kas==0, by(race_ethnic)
stsum if kas==1, by(race_ethnic)
 

* 2c) AFT Models
*--------------------------------------------------
* Run models by era
*--------------------------------------------------
// Unadjusted (Pre-KAS)
streg i.race_ethnic if kas == 0, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Dialysis1.xls", replace ctitle("Unadjusted - Pre-KAS") ///
        excel stats(coef ci pval) eform 

// full model (Pre-KAS) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.race_ethnic i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa dialysis_year  if kas == 0, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Dialysis1.xls", append ctitle("Adjusted - Pre-KAS") ///
        excel stats(coef ci pval) eform 


// Unadjusted (Post-KAS)
streg i.race_ethnic if kas==1, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Dialysis1.xls", append ctitle("Unadjusted - Post-KAS") ///
        excel stats(coef ci pval) eform 


// full model (Post-KAS) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.race_ethnic i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa dialysis_year  if kas == 1, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Era_Dialysis1.xls", append ctitle("Adjusted - Post-KAS") ///
        excel stats(coef ci pval) eform 


*--------------------------------------------------
* Run models by race
*--------------------------------------------------
// Unadjusted (White)
streg i.kas if race_ethnic == 1, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", replace ctitle("Unadjusted - White") ///
        excel stats(coef ci pval) eform 

// full model (White) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa dialysis_year  if race_ethnic == 1, dist(logl) time tr 
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Adjusted - White") ///
        excel stats(coef ci pval) eform 


// Unadjusted (Black)
streg i.kas if race_ethnic == 2, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Unadjusted - Black") ///
        excel stats(coef ci pval) eform 

// full model (Black) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa dialysis_year if race_ethnic == 2, dist(logl) time tr 
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Adjusted - Black") ///
        excel stats(coef ci pval) eform 

// Unadjusted (Hispanic)
streg i.kas if race_ethnic == 3, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Unadjusted - Hispanic") ///
        excel stats(coef ci pval) eform 


// full model (Hispanic) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa dialysis_year  if race_ethnic == 3, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Adjusted - Hispanic") ///
        excel stats(coef ci pval) eform 

// Unadjusted (Other)
streg i.kas if race_ethnic == 4, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Unadjusted - Other") ///
        excel stats(coef ci pval) eform 

// full model (Other) that also controls peak cPRA categories, primary diagnosis, payer type, blood type, and age categories
streg i.kas i.prim_dx i.payer i.peak_cpra i.gender i.age_list_cat ///
      i.bloodtype i.dsa dialysis_year if race_ethnic == 4, dist(logl) time tr
outreg2 using "$out_path/Logl_by_Race_Dialysis1.xls", append ctitle("Adjusted - Other") ///
        excel stats(coef ci pval) eform 

*--------------------------------------------------

*--------------------------------------------------
* 3) Time to Graft Loss
*--------------------------------------------------
// Open dataset
use 4yrspeds_final_info_v2.dta, clear

// Graft Loss
gen     graft_loss = (tfl_graft_dt != .)
replace graft_loss = 0 if (tfl_graft_dt > td(31dec2019))

// KAS indicator variable (Transplanted after KAS)
gen kas = (rec_tx_dt >= td(04dec2014))

// Drop those who were not transplanted
drop if rec_tx_dt == .

// Drop if LDKT
drop if don_ty=="L"


// Time variable (days to graft loss or censor)
gen days_to_censor = tfl_graft_dt - rec_tx_dt
    replace days_to_censor = tfl_death_dt - rec_tx_dt if days_to_censor == .
    replace days_to_censor = td(31dec2019) - rec_tx_dt if days_to_censor == .
    replace days_to_censor = tfl_lafudate - rec_tx_dt if days_to_censor == . 
    replace days_to_censor = can_endwlfu - rec_tx_dt if days_to_censor == .                                                    

// Allow for up to 5 years of follow-up.
stset days_to_censor, failure(graft_loss==1) id(pers_id) // exit(time 1825)   

// White log-rank test and KM curves
sts test kas if race_ethnic==1
sts graph if race_ethnic==1, by(kas) ///
        legend(on order(1 "Pre-KAS" 2 "Post-KAS")) title("") tmax(1500) /// 
        text(.7 500 "Log-rank p-value = 0.15", size(medium))
graph export "$out_path/KM_GL_White.png", replace
 
// Black log-rank test and KM curves
sts test kas if race_ethnic==2
sts graph if race_ethnic==2, by(kas) ///
        legend(on order(1 "Pre-KAS" 2 "Post-KAS")) title("") tmax(1500) /// 
        text(.7 500 "Log-rank p-value = >0.01", size(medium)) 
graph export "$out_path/KM_GL_Black.png", replace

// Hispanic log-rank test and KM curves
sts test kas if race_ethnic==3
sts graph if race_ethnic==3, by(kas) ///
        legend(on order(1 "Pre-KAS" 2 "Post-KAS")) title("") tmax(1500) /// 
        text(.7 500 "Log-rank p-value = >0.01", size(medium))
graph export "$out_path/KM_GL_Hispanic.png", replace

// Other log-rank test and KM curves
sts test kas if race_ethnic==4
sts graph if race_ethnic==4, by(kas) ///
        legend(on order(1 "Pre-KAS" 2 "Post-KAS")) title("") tmax(1500) ///
        text(.7 500 "Log-rank p-value = 0.85", size(medium))
graph export "$out_path/KM_GL_Other.png", replace

log close


