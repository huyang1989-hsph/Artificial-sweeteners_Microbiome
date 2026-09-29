libname readfmt1 '~/fecal_questionnaire/formats';
libname readfmt  '~/formats';
options nocenter ls=125 ps=78 replace formdlim='=';
options mautosource sasautos=(channing hpstools);
options fmtsearch=(readfmt);
options fmtsearch=(readfmt1);
filename hpstools '~/sasautos';
filename channing '~/channing/sasautos';
options nocenter minoperator ls=200 ps=200 replace;

%include "~/NDS_format.sas";

libname fddr '~/raw_data';

data wfood;
	set fddr.fscrecord;
proc sort;  by dintake;
proc sort;  by id;
run;

* processing several mistaks in data;

data wfood;  set wfood;
	
	/* there are two obvious mistake in the date info 
	115553	12/10/2010
	713146	10/29/2010
	2010 should be 2012 */

	if id=115553 and dintake='10dec2010'd then dintake='10dec2012'd;
 	if id=713146 and dintake='29oct2010'd then dintake='29oct2012'd;


	/* 307838 have info of three weeks

	 307838	8/13/2012 *matched to 1st 2 facal samples
	 307838	8/14/2012
	 307838	8/15/2012
	 307838	8/16/2012
	 307838	8/17/2012
	 307838	8/18/2012
	 307838	8/19/2012

	 307838	2/11/2013 * NOT matched to facal samples
	 307838	2/12/2013
	 307838	2/13/2013
	 307838	2/14/2013
	 307838	2/15/2013
	 307838	2/16/2013
	 307838	2/17/2013

	 307838	3/29/2013 * matched to 3-4 facal samples
	 307838	3/30/2013
	 307838	3/31/2013
	 307838	4/1/2013
	 307838	4/2/2013
	 307838	4/3/2013
	 307838	4/4/2013 */

	 if id=307838 and dintake in ('11feb2013'd, '12feb2013'd, '13feb2013'd, '14feb2013'd, 
	 	'15feb2013'd, '16feb2013'd, '17feb2013'd) then delete;


	 /* 801281 has 19 food record 
	 801281	4/22/2012 * NOT matched to facal samples
	 801281	4/23/2012
	 801281	4/24/2012
	 801281	4/25/2012
	 801281	4/26/2012
	 801281	4/27/2012
	 801281	4/28/2012

	 801281	7/22/2012 * matched to 1-2 fecal samples, not 3-4 sampels available (7/26 and 7/29)
	 801281	7/23/2012
	 801281	7/24/2012
	 801281	7/25/2012
	 801281	7/26/2012

	 801281	1/20/2013 * NOT matched to facal samples
	 801281	1/21/2013
	 801281	1/22/2013
	 801281	1/23/2013
	 801281	1/24/2013
	 801281	1/25/2013
	 801281	1/26/2013 

	 in this case, in order to make it easier for later match, drop the first 7 record, 
	 repeat the 7/25/2012 & 7/26/2012 so that the Jul record can be matched to 1-2 fecals */

	 if id=801281 and dintake in ('22apr2012'd, '23apr2012'd, '24apr2012'd, '25apr2012'd, 
	 	'26apr2012'd, '27apr2012'd, '28apr2012'd) then delete;
run;

data toadd;  set wfood;
	 if id=801281 and dintake in ('25jul2012'd, '26jul2012'd) then output ;
run;

data wfood;
	set wfood toadd;
run;


* seq the data;
data wfood;
	set wfood;
proc sort;  by dintake;
proc sort;  by id;
run;

data wfood;  set wfood;
	by id;
	if first.id then SeqN=1;
	else SeqN + 1;
run;

proc freq data=wfood;
	tables SeqN;
run;

/****N     U     R     I     N     T     S****/
/* ------------------ *
*    week 1 no adj   
* ------------------ */

libname ddr '~/diet_7ddr/sas_data';

*  day 1 ;

data ntsw1d1; set ddr.dr_nts_wk1_day1; drop iamount_w1d1 cpjname_dr_w1d1;
proc sort nodupkey; by id; run;
data foodw1d1; set wfood; rename dintake=ddr_date_w1d1; if SeqN=1 then output; proc sort; by id ; run;

*  day 2 ;
data ntsw1d2; set ddr.dr_nts_wk1_day2; drop iamount_w1d2 cpjname_dr_w1d2;
proc sort nodupkey; by id; run;
data foodw1d2; set wfood; rename dintake=ddr_date_w1d2; if SeqN=2 then output; proc sort; by id ; run;

*  day 3 ;
data ntsw1d3; set ddr.dr_nts_wk1_day3; drop iamount_w1d3 cpjname_dr_w1d3;
proc sort nodupkey; by id; run;
data foodw1d3; set wfood; rename dintake=ddr_date_w1d3; if SeqN=3 then output; proc sort; by id ; run;	

*  day 4 ;
data ntsw1d4; set ddr.dr_nts_wk1_day4; drop iamount_w1d4 cpjname_dr_w1d4;
proc sort nodupkey; by id; run;
data foodw1d4; set wfood; rename dintake=ddr_date_w1d4; if SeqN=4 then output; proc sort; by id ; run;

*  day 5 ;
data ntsw1d5; set ddr.dr_nts_wk1_day5; drop iamount_w1d5 cpjname_dr_w1d5;
proc sort nodupkey; by id; run;
data foodw1d5; set wfood; rename dintake=ddr_date_w1d5; if SeqN=5 then output; proc sort; by id ; run;

*  day 6 ;
data ntsw1d6; set ddr.dr_nts_wk1_day6; drop iamount_w1d6 cpjname_dr_w1d6;
proc sort nodupkey; by id; run;
data foodw1d6; set wfood; rename dintake=ddr_date_w1d6; if SeqN=6 then output; proc sort; by id ; run;

*  day 7 ;
data ntsw1d7; set ddr.dr_nts_wk1_day7; drop iamount_w1d7 cpjname_dr_w1d7;
proc sort nodupkey; by id; run;
data foodw1d7; set wfood; rename dintake=ddr_date_w1d7; if SeqN=7 then output; proc sort; by id ; run;

*  all mean ;
data ntsw1avg; set ddr.dr_nts_wk1_mean;  
proc sort nodupkey; by id; run;

* ------------------ *
*    week 2 no adj   
* ------------------ *

*  day 1 ;
data ntsw2d1; set ddr.dr_nts_wk2_day1; drop iamount_w2d1 cpjname_dr_w2d1;
proc sort nodupkey; by id; run;
data foodw2d1; set wfood; rename dintake=ddr_date_w2d1; if SeqN=8 then output; proc sort; by id ; run;

*  day 2 ;
data ntsw2d2; set ddr.dr_nts_wk2_day2; drop iamount_w2d2 cpjname_dr_w2d2;
proc sort nodupkey; by id; run;
data foodw2d2; set wfood; rename dintake=ddr_date_w2d2; if SeqN=9 then output; proc sort; by id ; run;

*  day 3 ;
data ntsw2d3; set ddr.dr_nts_wk2_day3; drop iamount_w2d3 cpjname_dr_w2d3;
proc sort nodupkey; by id; run;
data foodw2d3; set wfood; rename dintake=ddr_date_w2d3; if SeqN=10 then output; proc sort; by id ; run;

*  day 4 ;
data ntsw2d4; set ddr.dr_nts_wk2_day4; drop iamount_w2d4 cpjname_dr_w2d4;
proc sort nodupkey; by id; run;
data foodw2d4; set wfood; rename dintake=ddr_date_w2d4; if SeqN=11 then output; proc sort; by id ; run;

*  day 5 ;
data ntsw2d5; set ddr.dr_nts_wk2_day5; drop iamount_w2d5 cpjname_dr_w2d5;
proc sort nodupkey; by id; run;
data foodw2d5; set wfood; rename dintake=ddr_date_w2d5; if SeqN=12 then output; proc sort; by id ; run;

*  day 6 ;
data ntsw2d6; set ddr.dr_nts_wk2_day6; drop iamount_w2d6 cpjname_dr_w2d6;
proc sort nodupkey; by id; run;
data foodw2d6; set wfood; rename dintake=ddr_date_w2d6; if SeqN=13 then output; proc sort; by id ; run;

*  day 7 ;
data ntsw2d7; set ddr.dr_nts_wk2_day7; drop iamount_w2d7 cpjname_dr_w2d7;
proc sort nodupkey; by id; run;
data foodw2d7; set wfood; rename dintake=ddr_date_w2d7; if SeqN=14 then output; proc sort; by id ; run;

*  all mean ;
data ntsw2avg; set ddr.dr_nts_wk2_mean;
proc sort nodupkey; by id; run;

*  2-wk mean ;
data ntswtavg; set ddr.dr_nts_wkt_mean;
proc sort nodupkey; by id; run;
/*******************************************************************************************************************************************************************/
/***************************************M      E      R      G      E***********************************************************************************************/
/*******************************************************************************************************************************************************************/
data mlvs_ddr_nts; 
merge
foodw1d1 foodw1d2 foodw1d3 foodw1d4 foodw1d5 foodw1d6 foodw1d7
foodw2d1 foodw2d2 foodw2d3 foodw2d4 foodw2d5 foodw2d6 foodw2d7

ntsw1d1 ntsw1d2 ntsw1d3 ntsw1d4 ntsw1d5 ntsw1d6 ntsw1d7 ntsw1avg 
ntsw2d1 ntsw2d2 ntsw2d3 ntsw2d4 ntsw2d5 ntsw2d6 ntsw2d7 ntsw2avg 

ntswtavg; 

by id;

run;

proc sort data=mlvs_ddr_nts; by id; run; 

libname here '~/DDR_PDI';

data mlvs_ddr_food; set here.ddr_pdi_mlvs; run; 

proc sort data=mlvs_ddr_food; by id; run; 

data mlvs_ddr; merge mlvs_ddr_nts mlvs_ddr_food; by id; run;

proc means data=mlvs_ddr;run;



data here.mlvs_ddr; set mlvs_ddr; run;	
endsas;
/*Export works bad for long colunms*/
proc export data=mlvs_ddr OUTFILE='~/food_biomarker/allddr.csv' dbms=csv replace; putnames=yes; run;

