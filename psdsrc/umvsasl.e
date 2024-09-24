/*
 * David Frey
 * University of Michigan Medicine Department of Radiology
 * Functional MRI Laboratory
 *
 * GE Medical Systems
 * Copyright (C) 1996-2003 The General Electric Company
 *
 * File Name : sparke.e
 * Language  : EPIC/ANSI C
 * Date      : 10-May-2023
 *
 */

@inline epic.h
@inline intwave.h

@global
/*********************************************************************
 *                   UMVSASL.E GLOBAL SECTION                        *
 *                                                                   *
 * Common code shared between the Host and IPG PSD processes.  This  *
 * section contains all the #define's, global variables and function *
 * declarations (prototypes).                                        *
 *********************************************************************/
#include <stdio.h>
#include <string.h>

#include "em_psd_ermes.in"
#include "grad_rf_umvsasl.globals.h"

#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "epic_error.h"
#include "epic_loadcvs.h"
#include "InitAdvisories.h"
#include "psdiopt.h"
#ifdef psdutil
#include "psdutil.h"
#endif
#include "psd_proto.h"
#include "epic_iopt_util.h"
#include "filter.h"

#include "umvsasl.h"

/* Define important values */
#define MAXWAVELEN 50000 /* Maximum wave length for gradients */
#define MAXNROTS 1000 /* Maximum number of rotations */
#define MAXNSHOTS 50 /* Maximum number of shots (waveforms) per rotation */
#define GAMMA 26754 /* Gyromagnetic ratio (rad/s/G) */
#define TIMESSI 120 /* SSP instruction time */
#define SPOIL_SEED 21001 /* rf spoiler seed */

@inline Prescan.e PSglobal
int debugstate = 1;

@ipgexport
/*********************************************************************
 *                 UMVSASL.E IPGEXPORT SECTION                       *
 *                                                                   *
 * Standard C variables of _any_ type common for both the Host and   *
 * IPG PSD processes. Declare here all the complex type, e.g.,       *
 * structures, arrays, files, etc.                                   *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG schedule_ides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
@inline Prescan.e PSipgexport
RF_PULSE_INFO rfpulseInfo[RF_FREE] = { {0,0} };

/* Define temporary error message string */
char tmpstr[200];

/* Declare sequencer hardware limit variables */
float XGRAD_max;
float YGRAD_max;
float ZGRAD_max;
float RHO_max;
float THETA_max;
int ZGRAD_risetime;
int ZGRAD_falltime;

int nrots, nshots;

/* Declare readout gradient waveform arrays */
int Gx[MAXNSHOTS][MAXWAVELEN];
int Gy[MAXNSHOTS][MAXWAVELEN];
int Gz[MAXNSHOTS][MAXWAVELEN];
int grad_len = 5000;
int acq_len = 4000;
int acq_off = 50;

/* Declare table of readout gradient transformation matrices */
long rotmattbl[MAXNROTS][9];

/* Declare receiver and Tx frequencies */
float recfreq;
float xmitfreq;

@cv
/*********************************************************************
 *                      UMVSASL.E CV SECTION                         *
 *                                                                   *
 * Standard C variables of _limited_ types common for both the Host  *
 * and IPG PSD processes. Declare here all the simple types, e.g,    *
 * int, float, and C structures containing the min and max values,   *
 * and ID description, etc.                                          *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG schedule_ides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
@inline loadrheader.e rheadercv
@inline vmx.e SysCVs

@inline Prescan.e PScvs

int numdda = 4;			/* For Prescan: # of disdaqs ps2*/

float smax = 12500.0 with {1000, 25000.0, 12500.0, VIS, "maximum allowed slew rate (G/cm/s)",};
float gmax = 4.0 with {0.5, 5.0, 4.0, VIS, "maximum allowed gradient (G/cm)",};

/* readout waveform cvs */
int ndisdaqs = 0 with {0, , 0, VIS, "number of disdaqs beginning of sequence",};
int trajid = 0 with {0, , 99999, VIS, "trajectory file id (i.e. /usr/g/bin/sparky#####.h5)",};

/* fat sat cvs */ 
int fatsup_flag = 1 with {0, 1, 1, VIS, "off (0), CHESS (1),",};
int fatsup_off = -520 with { , , -520, VIS, "fat suppression pulse frequency offset (Hz)",};
int fatsup_bw = 440 with { , , 440, VIS, "fat suppression bandwidth (Hz)",};
int rfspoil_flag = 1 with {0, 1, 1, VIS, "option to do RF phase cycling (117deg increments to rf1 phase)",};

int pgbuffertime = 248 with {100, , 248, INVIS, "gradient IPG buffer time (us)",};
float crushfac = 3.0 with {0, 10, 0, VIS, "crusher amplitude factor (a.k.a. cycles of phase/vox; dk_crush = crushfac*kmax)",};
int kill_grads = 0 with {0, 1, 0, VIS, "option to turn off readout gradients",};

/* Declare core duration variables */
int dur_rf1core = 0 with {0, , 0, INVIS, "duration of the slice selective rf1 core (us)",};
int dur_seqcore = 0 with {0, , 0, INVIS, "duration of the spiral readout core (us)",};
int deadtime1_seqcore = 0 with {0, , 0, INVIS, "pre-readout deadtime within core (us)",};
int deadtime2_seqcore = 0 with {0, , 0, INVIS, "post-readout deadtime within core (us)",};

/* inhereted from grass.e, not sure if it's okay to delete: */
float xmtaddScan;
int obl_debug = 0 with {0, 1, 0, INVIS, "On(=1) to print messages for obloptimize",};
int obl_method = 0 with {0, 1, 0, INVIS, "On(=1) to optimize the targets based on actual rotation matrices",};
int debug = 0 with {0,1,0,INVIS,"1 if debug is on ",};
float echo1bw = 16 with {,,,INVIS,"Echo1 filter bw.in KHz",};

@host
/*********************************************************************
 *                     UMVSASL.E HOST SECTION                        *
 *                                                                   *
 * Write here the code unique to the Host PSD process. The following *
 * functions must be declared here: cvinit(), cveval(), cvcheck(),   *
 * and predownload().                                                *
 *                                                                   *
 *********************************************************************/
#include <math.h>
#include <stdlib.h>
#include "grad_rf_umvsasl.h"
#include "psdopt.h"
#include "sar_pm.h"
#include "support_func.host.h"
#include "helperfuns.h"
#include "gengrads.h"

/* fec : Field strength dependency library */
#include <sysDep.h>
#include <sysDepSupport.h>      /* FEC : fieldStrength dependency libraries */

@inline loadrheader.e rheaderhost

/** Load PSD Header **/
abstract("SPGR sequence with arbitrary kspace trajectory");
psdname("sparky");

int num_conc_grad = 3;          /* always three for grass 	*/
int entry;

/* peak B1 amplitudes */
float maxB1[MAX_ENTRY_POINTS], maxB1Seq;

/* This will point to a structure defining parameters of the filter
   used for the 1st echo */
FILTER_INFO *echo1_filt; 

/* Use real time filters, so allocate space for them instead of trying
   to point to an infinite number of structures in filter.h. */
FILTER_INFO echo1_rtfilt;

/* declare function prototypes */
int gengrads(int id, float dt, float gmax, float smax,
		float ***gx, float ***gy, float ***gz,
		int *n_shots, int *n_all, int *n_wav, int *n_off);
int genrots();

/* declare function prototypes from aslprep.h */
float calc_sinc_B1(float cyc_rf, int pw_rf, float flip_rf);
float calc_hard_B1(int pw_rf, float flip_rf);
int write_scan_info();

@inline Prescan.e PShostVars            /* added with new filter calcs */

static char supfailfmt[] = "Support routine %s failed";


/************************************************************************/
/*       			CVINIT    				*/
/* Invoked once (& only once) when the PSD host process	is started up.	*/
/* Code which is independent of any OPIO button operation is put here.	*/
/************************************************************************/
STATUS cvinit( void )
{

	/* turn off bandwidth option */
	cvdef(oprbw, 500.0 / (float)GRAD_UPDATE_TIME);
	cvmin(oprbw, 500.0 / (float)GRAD_UPDATE_TIME);
	cvmax(oprbw, 500.0 / (float)GRAD_UPDATE_TIME);
	oprbw = 500.0 / (float)GRAD_UPDATE_TIME;
	pircbnub = 0;

	/* fov */
	opfov = 240;
	pifovnub = 5;
	pifovval2 = 200;
	pifovval3 = 220;
	pifovval4 = 240;
	pifovval5 = 260;
	pifovval6 = 280;

	/* tr */
	opautotr = PSD_MINIMUMTR;
	pitrnub = 2;
	pitrval2 = PSD_MINIMUMTR;
	cvmax(optr, 50s);

	/* te */
	opautote = PSD_MINTE;	
	pite1nub = 3;
	pite1val2 = PSD_MINTE;
	cvmin(opte, 0);
	cvmax(opte, 500ms);

	/* rhrecon */
	rhrecon = 6969;

	/* frequency (xres) */
	opxres = 128;
	cvmin(opxres, 16);
	cvmax(opxres, 512);
	pixresnub = 15;
	pixresval2 = 32;
	pixresval3 = 64;
	pixresval4 = 128;

	/* flip angle */
	cvmin(opflip, 0.0);
	cvmax(opflip, 360.0);
	pifanub = 2;
	pifaval2 = 90.0;

	/* hide phase (yres) option */
	piyresnub = 0;

	/* Hide inversion time */
	pitinub = 0;

	/* hide second bandwidth option */
	pircb2nub = 0;

	/* hide nex stuff */
	piechnub = 0;
	pinexnub = 0;

#ifdef ERMES_DEBUG
	use_ermes = 0;
#else /* !ERMES_DEBUG */
	use_ermes = 1;
#endif /* ERMES_DEBUG */

	configSystem();
	EpicConf();
	inittargets(&loggrd, &phygrd);

	/* Init filter slots */
	initfilter();
	
	if (_psd_rf_wait.fixedflag == 0)  { /* sets psd_grd_wait and psd_rf_wait */
		if (setsysparms() == FAILURE)  {
			epic_error(use_ermes,"Support routine setsysparams failed",
					EM_PSD_SUPPORT_FAILURE,1, STRING_ARG,"setsysparms");
			return FAILURE;
		}
	}

	if( obloptimize( &loggrd, &phygrd, scan_info, exist(opslquant),
				exist(opplane), exist(opcoax), obl_method, obl_debug,
				&opnewgeo, cfsrmode ) == FAILURE )
	{
		return FAILURE;
	}
	
	/* Get sequencer hardware limits */
	gettarget(&XGRAD_max, XGRAD, &loggrd);
	gettarget(&YGRAD_max, YGRAD, &loggrd);
	gettarget(&ZGRAD_max, ZGRAD, &loggrd);
	gettarget(&RHO_max, RHO, &loggrd);
	gettarget(&THETA_max, THETA, &loggrd);
	getramptime(&ZGRAD_risetime, &ZGRAD_falltime, ZGRAD, &loggrd);	
	ZGRAD_risetime *= 2; /* extra fluffy */
	fprintf(stderr, "ZGRAD_risetime = %d\n", ZGRAD_risetime);	

@inline Prescan.e PScvinit

#include "cvinit.in"	/* Runs the code generated by macros in preproc.*/

	return SUCCESS;
}   /* end cvinit() */

@inline InitAdvisories.e InitAdvPnlCVs

/************************************************************************/
/*       			CVEVAL    				*/
/* Called w/ every OPIO button push which has a corresponding CV. 	*/
/* CVEVAL should only contain code which impacts the advisory panel--	*/
/* put other code in cvinit or predownload				*/
/************************************************************************/
STATUS cveval( void )
{
	configSystem();
	InitAdvPnlCVs();

	pititle = 1;
	cvdesc(pititle, "Advanced pulse sequence parameters");
	piuset = 0;

	piuset += use0;
	cvdesc(opuser0, "trajectory file ID #");
	cvdef(opuser0, trajid);
	opuser0 = trajid;
	cvmin(opuser0, 0);
	cvmax(opuser0, 99999);
	trajid = opuser0;

	piuset += use1;
	cvdesc(opuser1, "number of disdaqs");
	cvdef(opuser1, ndisdaqs);
	opuser1 = ndisdaqs;
	cvmin(opuser1, 0);
	cvmax(opuser1, 500);
	ndisdaqs = opuser1;
	
	piuset += use2;
	cvdesc(opuser2, "crusher area factor (% N/fov/2)");
	cvdef(opuser2, crushfac);
	opuser2 = crushfac;
	cvmin(opuser2, 0);
	cvmax(opuser2, 10);
	crushfac = opuser2;
		
	piuset += use3;
	cvdesc(opuser3, "fat suppression (0) off or (1) on");
	cvdef(opuser3, fatsup_flag);
	opuser3 = fatsup_flag;
	cvmin(opuser3, 0);
	cvmax(opuser3, 1);	
	fatsup_flag = opuser3;

	piuset += use4;
	cvdesc(opuser4, "rf spoiling (0) off or (1) on");
	cvdef(opuser4, rfspoil_flag);
	opuser4 = rfspoil_flag;
	cvmin(opuser4, 0);
	cvmax(opuser4, 1);	
	rfspoil_flag = opuser4;

@inline Prescan.e PScveval

	return SUCCESS;
}   /* end cveval() */

void getAPxParam(optval   *min,
		optval   *max,
		optdelta *delta,
		optfix   *fix,
		float    coverage,
		int      algorithm)
{
	/* Need to be filled when APx is supported in this PSD */
}

int getAPxAlgorithm(optparam *optflag, int *algorithm)
{
	return APX_CORE_NONE;
}

/************************************************************************/
/*       			CVCHECK    				*/
/* Executed on each 'next page' to ensure prescription can proceed 	*/
/* to the next page. 							*/
/************************************************************************/
STATUS cvcheck( void )
{
	return SUCCESS;
}   /* end cvcheck() */


/************************************************************************/
/*             		    PRE-DOWNLOAD           		        */
/* Executed prior to a download--all operations not needed for the 	*/
/* advisory panel results.  Execute the	pulsegen macro expansions for	*/
/* the predownload section here.  All internal amps, slice ordering,  	*/
/* prescan slice calc., and SAT placement calculations are performed 	*/
/* in this section.  Time anchor settings for pulsegen are done in this */
/* section too.  				 			*/
/************************************************************************/
STATUS predownload( void )
{
	int echo1_freq[opslquant], rf1_freq[opslquant];
	int shot, n, slice;
	int minte, mintr;
	float rf1_b1, rffs_b1;
	int tmp_pwa, tmp_pw, tmp_pwd;
	float tmp_a, tmp_area;
	float **gx, **gy, **gz;

	/*********************************************************************/
#include "predownload.in"	/* include 'canned' predownload code */
	/*********************************************************************/	
	
	/* update sinc pulse parameters */
	pw_rf1 = 3200;
	
	/* adjust fat sup pw s.t. desired bandwidth is achieved */
	pw_rffs = 3200; /* nominal SINC1 pulse width */
	pw_rffs *= (int)round(NOM_BW_SINC1_90 / (float)fatsup_bw); /* adjust bandwidth */

	/* first, find the peak B1 for all entry points (other than L_SCAN) */
	for( entry=0; entry < MAX_ENTRY_POINTS; ++entry )
	{
		if( peakB1( &maxB1[entry], entry, RF_FREE, rfpulse ) == FAILURE )
		{
			epic_error( use_ermes, "peakB1 failed.", EM_PSD_SUPPORT_FAILURE,
					EE_ARGS(1), STRING_ARG, "peakB1" );
			return FAILURE;
		}
	}
	
	/* calculate the b1 samplitudes of each pulse */
	rf1_b1 = calc_sinc_B1(cyc_rf1, pw_rf1, opflip);
	fprintf(stderr, "predownload(): maximum B1 for rf1 pulse: %f\n", rf1_b1);
	if (rf1_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = rf1_b1;
		
	rffs_b1 = calc_sinc_B1(cyc_rffs, pw_rffs, 90.0);
	fprintf(stderr, "predownload(): maximum B1 for fatsup pulse: %f Gauss\n", rffs_b1);
	if (rffs_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = rffs_b1;
	
	/* Determine peak B1 across all entry points */
	maxB1Seq = 0.0;
	for (entry=0; entry < MAX_ENTRY_POINTS; entry++) {
		if (entry != L_SCAN) { /* since we aleady computed the peak B1 for L_SCAN entry point */
			if (peakB1(&maxB1[entry], entry, RF_FREE, rfpulse) == FAILURE) {
				epic_error(use_ermes,"peakB1 failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"peakB1");
				return FAILURE;
			}
		}
		if (maxB1[entry] > maxB1Seq)
			maxB1Seq = maxB1[entry];
	}
	fprintf(stderr, "predownload(): maxB1Seq = %f Gauss\n", maxB1Seq);
	
	/* Set xmtadd according to maximum B1 and rescale for powermon,
	   adding additional (audio) scaling if xmtadd is too big.
	   Add in coilatten, too. */
	xmtaddScan = -200 * log10( maxB1[L_SCAN] / maxB1Seq ) + getCoilAtten(); 

	if( xmtaddScan > cfdbmax )
	{
		extraScale = (float)pow( 10.0, (cfdbmax - xmtaddScan) / 200.0 );
		xmtaddScan = cfdbmax;
	} 
	else
	{
		extraScale = 1.0;
	}
	
	/* update rf amplitudes */
	a_rf1 = rf1_b1 / maxB1Seq;
	ia_rf1 = a_rf1 * MAX_PG_WAMP;
	
	a_rffs = rffs_b1 / maxB1Seq;
	ia_rffs = a_rffs * MAX_PG_WAMP;
	
	/* set ss rewinder trapezoid gradient */
	tmp_area = a_gzrf1 * (pw_gzrf1 + (pw_gzrf1a + pw_gzrf1d)/2.0);
	amppwgrad(tmp_area, gmax, 0, 0, ZGRAD_risetime, 0, &tmp_a, &tmp_pwa, &tmp_pw, &tmp_pwd); 
	tmp_a *= -0.5;
	pw_gzrf1r = tmp_pw;
	pw_gzrf1ra = tmp_pwa;
	pw_gzrf1rd = tmp_pwd;
	a_gzrf1r = tmp_a;

	/* set crusher trapezoid gradient */
	tmp_area = crushfac * 2*M_PI/GAMMA * opxres/(opfov/10.0) * 1e6; /* area under crusher s.t. dk = crushfac*kmax (G/cm*us) */
	amppwgrad(tmp_area, gmax, 0, 0, ZGRAD_risetime, 0, &tmp_a, &tmp_pwa, &tmp_pw, &tmp_pwd); 	
	
	/* apply to GRE spoiler */
	pw_gzspoil = tmp_pw;
	pw_gzspoil = tmp_pwa;
	pw_gzspoil = tmp_pwd;
	a_gzspoil = tmp_a;

	/* tmp */
	a_gzspoil = 0;

	/* generate initial gradient waveforms */
	fprintf(stderr, "predownload(): generating gradients...\n");
	if (gengrads(trajid, GRAD_UPDATE_TIME*1e-6, gmax, smax,
		&gx, &gy, &gz, &nshots, &grad_len, &acq_len, &acq_off) == 0) {
		epic_error(use_ermes,"failure to generate gradients", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	for (shot = 0; shot < nshots; shot++) {
		for (n = 0; n < grad_len; n++) {
			Gx[shot][n] = 2*round(MAX_PG_WAMP/XGRAD_max * gx[shot][n] / 2);
			Gy[shot][n] = 2*round(MAX_PG_WAMP/YGRAD_max * gy[shot][n] / 2);
			Gz[shot][n] = 2*round(MAX_PG_WAMP/ZGRAD_max * gz[shot][n] / 2);
		}
	}
	a_gxw = XGRAD_max;
	a_gyw = YGRAD_max;
	a_gzw = ZGRAD_max;
	ia_gxw = MAX_PG_WAMP;
	ia_gyw = MAX_PG_WAMP;
	ia_gzw = MAX_PG_WAMP;
	res_gxw = grad_len;
	res_gyw = grad_len;
	res_gzw = grad_len;
	pw_gxw = GRAD_UPDATE_TIME*res_gxw;
	pw_gyw = GRAD_UPDATE_TIME*res_gyw;
	pw_gzw = GRAD_UPDATE_TIME*res_gzw;	

	/* generate rotations */
	if (genrots() == 0) {
		epic_error(use_ermes,"failure to generate view transformation matrices", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	scalerotmats(rotmattbl, &loggrd, &phygrd, nrots, 0);
	
	/* calculate duration of rf1core */
	dur_rf1core = 0;
	if (fatsup_flag) {
		dur_rf1core += pgbuffertime;
		dur_rf1core += pw_rffs;
	}
	dur_rf1core += pgbuffertime;
	dur_rf1core += pw_gzspoila + pw_gzspoil + pw_gzspoild;
	dur_rf1core += pgbuffertime;
	dur_rf1core += pw_gzrf1a + pw_gzrf1 + pw_gzrf1d;
	dur_rf1core += pgbuffertime;
	dur_rf1core += pw_gzrf1ra + pw_gzrf1r + pw_gzrf1rd;
	dur_rf1core += pgbuffertime; 

	/* calculate the minimum echo time (time from center of rf1 pulse to beginning of RO grads) */
	minte = 0;
	minte += pw_gzrf1/2 + pw_gzrf1d; /* 2nd half of rf1 pulse */
	minte += pgbuffertime;
	minte += pw_gzrf1ra + pw_gzrf1r + pw_gzrf1rd;
	minte += pgbuffertime;
	minte += TIMESSI;
	minte += pgbuffertime;

	/* set deadtime1_seqcore */
	deadtime1_seqcore = opte - minte;

	/* calculate the minimum tr */
	mintr = 0;
	mintr += dur_rf1core + TIMESSI;
	mintr += pgbuffertime;
	mintr += deadtime1_seqcore;
	mintr += pw_gxw;
	mintr += pgbuffertime;

	/* set deadtime2_seqcore */
	deadtime2_seqcore = optr - mintr;	

	/* calculate duration of seqcore */
	dur_seqcore = 0;
	dur_seqcore += pgbuffertime;
	dur_seqcore += deadtime1_seqcore;
	dur_seqcore += pw_gxw;
	dur_seqcore += pgbuffertime;
	dur_seqcore += deadtime2_seqcore;

	/* set minimums */	
	cvmin(optr, mintr);
	cvmin(opte, minte);
	
	/* 
	 * Calculate RF filter and update RBW:
	 *   &echo1_rtfilt: I: all the filter parameters.
	 *   exist(oprbw): I/O: desired and final allowable bw.
	 *   exist(opxres): I: output pts generated by filter.
	 *   OVERWRITE_OPRBW: oprbw will be updated.
	 */
	if( calcfilter( &echo1_rtfilt,
				exist(oprbw),
				acq_len,
				OVERWRITE_OPRBW ) == FAILURE)
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "calcfilter:echo1" );
		return FAILURE;
	}

	echo1_filt = &echo1_rtfilt;

	/* Divide by 0 protection */
	if( (echo1_filt->tdaq == 0) || 
			floatsAlmostEqualEpsilons(echo1_filt->decimation, 0.0f, 2) ) 
	{
		epic_error( use_ermes, "echo1 tdaq or decimation = 0",
				EM_PSD_BAD_FILTER, EE_ARGS(0) );
		return FAILURE;
	}

	/* For use on the RSP schedule_ide */
	echo1bw = echo1_filt->bw;

@inline Prescan.e PSfilter

	/* For Prescan: Inform 'Auto' Prescan about prescan parameters 	*/
	pislquant = 10;	/* # of 2nd pass slices */

	/* For Prescan: Declare the entry point table 	*/
	if( entrytabinit( entry_point_table, (int)ENTRY_POINT_MAX ) == FAILURE ) 
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "entrytabinit" );
		return FAILURE;
	}

	/* For Prescan: Define the entry points in the table */
	/* Scan Entry Point */
	(void)strcpy( entry_point_table[L_SCAN].epname, "scan" );
	entry_point_table[L_SCAN].epfilter = (unsigned char)echo1_filt->fslot;
	entry_point_table[L_SCAN].epprexres = acq_len;

	(void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
	entry_point_table[L_APS2].epfilter = (unsigned char)echo1_filt->fslot;
	entry_point_table[L_APS2].epprexres = acq_len;

	(void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );
	entry_point_table[L_MPS2].epfilter = (unsigned char)echo1_filt->fslot;
	entry_point_table[L_MPS2].epprexres = acq_len;

	/* set sequence clock */
	pidmode = PSD_CLOCK_NORM;
	pitslice = optr;
	pitscan = (nshots*nrots + ndisdaqs) * optr; /* pitscan controls the clock time on the interface */	
	
	/* Set up the filter structures to be downloaded for realtime 
	   filter generation. Get the slot number of the filter in the filter rack 
	   and assign to the appropriate acquisition pulse for the right 
	   filter selection - LxMGD, RJF */
	setfilter( echo1_filt, SCAN );
	filter_echo1 = echo1_filt->fslot;
	entry_point_table[L_SCAN].epxmtadd = (short)rint( (double)xmtaddScan );

	/* APS2 & MPS2 */
	entry_point_table[L_APS2] = entry_point_table[L_MPS2] = entry_point_table[L_SCAN];	/* copy scan into APS2 & MPS2 */
	(void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
	(void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );

	/* Set up Tx/Rx frequencies */
	for (slice = 0; slice < opslquant; slice++) rsp_info[slice].rsprloc = 0;
	setupslices(rf1_freq, rsp_info, opslquant, a_gzrf1, 1.0, opfov, TYPTRANSMIT);
	setupslices(echo1_freq, rsp_info, opslquant, 0.0, 1.0, 2.0, TYPREC);

	/* Average together all slice frequencies */
	xmitfreq = 0;
	recfreq = 0;	
	for (slice = 0; slice < opslquant; slice++) {
		xmitfreq += (float)rf1_freq[slice] / (float)opslquant;
		recfreq += (float)echo1_freq[slice] / (float)opslquant;
	}

	if( orderslice( TYPNORMORDER, MAXNROTS+1, MAXNROTS+1, TRIG_INTERN ) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "orderslice" );
	}

	/* nex, exnex, acqs and acq_type are used in the rhheaderinit routine */
	/* -- to initialize recon header variables */
	if( floatsAlmostEqualEpsilons(opnex, 1.0, 2) )
	{
		baseline = 8;
		nex = 1;
		exnex = 1;
	}
	else
	{
		baseline = 0;
		nex = opnex;
		exnex = opnex;
	}

@inline loadrheader.e rheaderinit   /* Recon variables */
	
	/* Set recon header variables:
	 *   rhptsize: number of bytes per data point
	 *   rhfrsize: number of data points per acquisition
	 *   rhrawsize: total number of bytes to allocate
	 *   rhrcctrl: recon image control (bitmap)
	 *   rhexecctrl: recon executive control (bitmap)
	 */ 
	cvmax(rhfrsize, 32767);
	cvmax(rhnframes, 32767);
	cvmax(rhnslices, 32767);

	rhfrsize = acq_len;
	rhnframes = 2*ceil((float)(nshots + 1) / 2.0);
	rhnecho = 1;
	rhnslices = nrots + 1;
	rhrawsize = 2*rhptsize*rhfrsize * (rhnframes + 1) * rhnslices * rhnecho;
	
	rhrcctrl = 1; /* bit 7 (2^7 = 128) skips all recon */
	rhexecctrl = 2; /* bit 1 (2^1 = 2) sets autolock of raw files + bit 3 (2^3 = 8) transfers images to disk */

	write_scan_info();

@inline Prescan.e PSpredownload	

	return SUCCESS;
}   /* end predownload() */


@inline Prescan.e PShost


@pg
/*********************************************************************
 *                  UMVSASL.E PULSEGEN SECTION                       *
 *                                                                   *
 * Write here the functional code that loads hardware sequencer      *
 * memory with data that will allow it to play out the sequence.     *
 * These functions call pulse generation macros previously defined   *
 * with @pulsedef, and must return SUCCESS or FAILURE.               *
 *********************************************************************/
#include "support_func.h"
#include "epicfuns.h"

WF_PULSE** xGrad;     /* element (i,j) is for block j in group i */
WF_PULSE** yGrad;
WF_PULSE** zGrad;

WF_HW_WAVEFORM_PTR** x_wf;
WF_HW_WAVEFORM_PTR** y_wf;
WF_HW_WAVEFORM_PTR** z_wf;

WF_PULSE makePulse(WF_PROCESSOR wfp, char *pname, short *wave, int res, int tbeg) {
/* adapted from toppe v6 source code by Jon Fredrik Nielsen */

	WF_PULSE *ep;
	WF_PULSE proto = INITPULSE;

	ep = (WF_PULSE *) AllocNode(sizeof(WF_PULSE));
	memcpy((char*) ep, (char*) &proto, sizeof(WF_PULSE));
	pulsename(ep, pname);
	createreserve(ep, wfp, res);

	movewaveimm(wave, ep, 0, res, TOHARDWARE);
	createinstr(ep, tbeg, GRAD_UPDATE_TIME*res, MAX_PG_IAMP);
	if (wfp == TYPRHO1) {
		/* addrfbits(ep, 0, tbeg, 4*res); */
		fastAddrfbits(ep, 0, tbeg, GRAD_UPDATE_TIME*res, 50us);
	}

	return(*ep);
}
void tp_wreserve(WF_PROCESSOR wfp, WF_HW_WAVEFORM_PTR *wave_addr, int n) 
{
    SeqData seqdata;
    getWaveSeqDataWavegen(&seqdata, wfp, 0, 0, 0, PULSE_CREATE_MODE);
    *wave_addr = wreserve(seqdata, n);
}

STATUS pulsegen( void )
{
	sspinit(psd_board_type);
	int tmploc;	

	/*************************/
	/* generate readout core */
	/*************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of seqcore\n");
	tmploc = 0;
	tmploc += deadtime1_seqcore + pgbuffertime; /* add pre-readout deadtime + buffer */
	
	fprintf(stderr, "pulsegen(): generating gxw, gyw, gzw (readout gradients) and echo1 (data acquisition window)...\n");
	
	INTWAVE(XGRAD, gxw, tmploc, XGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gx[0], 1, loggrd);
	INTWAVE(YGRAD, gyw, tmploc, YGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gy[0], 1, loggrd);
	INTWAVE(ZGRAD, gzw, tmploc, ZGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gz[0], 1, loggrd);
	ACQUIREDATA(echo1, tmploc + psd_grd_wait + GRAD_UPDATE_TIME*acq_off,,,);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gxw; /* end time for readout */
	fprintf(stderr, " end: %dus\n", tmploc);

	tmploc += pgbuffertime + deadtime2_seqcore; /* add pre-readout deadtime + buffer */

	fprintf(stderr, "pulsegen(): finalizing seqcore...\n");
	fprintf(stderr, "\ttotal time: %dus (tmploc = %dus)\n", dur_seqcore, tmploc);
	SEQLENGTH(seqcore, dur_seqcore, seqcore);
	fprintf(stderr, "\tDone.\n");

	
	/*************************/
	/* generate rf1 core */
	/*************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of rf1 core\n");
	tmploc = 0;

	if (fatsup_flag) {	
		fprintf(stderr, "pulsegen(): generating rffs (fat suppresion rf pulse)...\n");
		tmploc += pgbuffertime; /* start time for rffs */
		SINC(RHO, rffs, tmploc + psd_rf_wait, 3200, 1.0, ,0.5, , , loggrd);
		fprintf(stderr, "\tstart: %dus, ", tmploc);
		tmploc += pw_rffs; /* end time for rffs */
		fprintf(stderr, " end: %dus\n", tmploc);	
	}
	
	fprintf(stderr, "pulsegen(): generating gzspoil (pre-rf1 spoiler gradient)...\n");
	tmploc += pgbuffertime; /* start time for gzspoil */
	TRAPEZOID(ZGRAD, gzspoil, tmploc + pw_gzspoila, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzspoila + pw_gzspoil + pw_gzspoild; /* end time for gzspoil pulse */
	fprintf(stderr, " end: %dus\n", tmploc);

	fprintf(stderr, "pulsegen(): generating rf1 (rf1 pulse)...\n");
	tmploc += pgbuffertime; /* start time for rf1 */
	SLICESELZ(rf1, tmploc + pw_gzrf1a, 3200, (opslthick + opslspace)*opslquant, opflip, 2, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf1a + pw_gzrf1 + pw_gzrf1d; /* end time for rf2 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);

	fprintf(stderr, "pulsegen(): generating gzrf1r (post-rf1 refocuser gradient)...\n");
	tmploc += pgbuffertime; /* start time for gzrf1r */
	TRAPEZOID(ZGRAD, gzrf1r, tmploc + pw_gzrf1ra, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf1ra + pw_gzrf1r + pw_gzrf1rd; /* end time for gzrf1r pulse */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime;

	fprintf(stderr, "pulsegen(): finalizing rf1 core...\n");
	fprintf(stderr, "\ttotal time: %dus (tmploc = %dus)\n", dur_rf1core, tmploc);
	SEQLENGTH(rf1core, dur_rf1core, rf1core);
	fprintf(stderr, "\tDone.\n");


	/**********************************/
	/* generate deadtime (empty) core */
	/**********************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of emptycore\n");

	fprintf(stderr, "pulsegen(): finalizing empty core...\n");
	SEQLENGTH(emptycore, 1000, emptycore);
	fprintf(stderr, "\tDone.\n");


@inline Prescan.e PSpulsegen

	PASSPACK(endpass, 49ms);   /* tell Signa system we're done */
	SEQLENGTH(pass, 50ms, pass);

	buildinstr();              /* load the sequencer memory       */
	fprintf(stderr, "\tDone with pulsegen().\n");

	return SUCCESS;
}   /* end pulsegen() */


/* For Prescan: Pulse Generation functions */
@inline Prescan.e PSipg


@rspvar
/*********************************************************************
 *                    UMVSASL.E RSPVAR SECTION                       *
 *                                                                   *
 * Declare here the real time variables that can be viewed and modi- *
 * fied while the IPG PSD process is running. Only limited standard  *
 * C types are provided: short, int, long, float, double, and 1D     *
 * arrays of those types.                                            *
 *                                                                   *
 * NOTE: Do not declare all real-time variables here because of the  *
 *       overhead required for viewing and modifying them.           *
 *********************************************************************/
extern PSD_EXIT_ARG psdexitarg;

/* Declare rsps */
float rfphs;
int shotn;
int rotn;
int disdaqn;
int n;
int rspfct;
int rspsct;
int view;
int slice;
int echo;

/* Inherited from grass.e: */
int dabop;
int excitation;
int rspent;
int rspdda;
int rspbas;
int rspvus;
int rspgy1;
int rspasl;
int rspesl;
int rspchp;
int rspnex;
int rspslq;

/* For Prescan: K */
int seqCount;

@inline Prescan.e PSrspvar 


@rsp
/*********************************************************************
 *                    UMVSASL.E RSP SECTION                          *
 *                                                                   *
 * Write here the functional code for the real time processing (IPG  *
 * schedule_ide). You may declare standard C variables, but of limited types *
 * short, int, long, float, double, and 1D arrays of those types.    *
 *********************************************************************/
#include <math.h>

/* For IPG Simulator: will generate the entry point list in the IPG tool */
const CHAR *entry_name_list[ENTRY_POINT_MAX] = {
	"scan", 
	"aps2",
	"mps2",
@inline Prescan.e PSeplist
};

/* Do not move the line above and do not insert any code or blank
   lines before the line above.  The code inline'd from Prescan.e
   adds more entry points and closes the list. */

long rotmat0[9]; /* Initial transformation matrix */
long zmtx[9] = {0};

STATUS psdinit( void )
{

	/* Initialize everything to a known state */
	setrfconfig( ENBL_RHO1 + ENBL_THETA );
	setssitime( TIMESSI/GRAD_UPDATE_TIME );
	rspqueueinit( 200 );	/* Initialize to 200 entries */
	scopeon( &seqcore );	/* Activate scope for core */
	syncon( &seqcore );		/* Activate sync for core */
	syncoff( &pass );		/* Deactivate sync during pass */
	seqCount = 0;		/* Set SPGR sequence counter */
	setrotatearray( 1, rsprot[0] );
	settriggerarray( 1, rsptrigger );
	setrfltrs( (int)filter_echo1, &echo1 );
			
	/* Set rf1 tx and rx frequency */
	setfrequency((int)xmitfreq, &rf1, 0);
	setfrequency((int)recfreq, &echo1, 0);

	/* Set fat sup frequency */
	if (fatsup_flag)
		setfrequency( (int)(fatsup_off / TARDIS_FREQ_RES), &rffs, 0);
	
	/* Get the original rotation matrix */
	getrotate( rotmat0, 0 );

	return SUCCESS;
}   /* end psdinit() */


@inline Prescan.e PScore


/* PLAY_DEADTIME() function for playing deadtime */
int play_deadtime(int deadtime) {
	int ttotal = 0;
	fprintf(stderr, "\tplay_deadtime(): playing deadtime (%d us)...\n", deadtime);

	/* Play empty core */
	setperiod(deadtime - TIMESSI, &emptycore, 0);
	boffset(off_emptycore);
	startseq(0, MAY_PAUSE);
	settrigger(TRIG_INTERN, 0);
	ttotal += deadtime;

	fprintf(stderr, "\tplay_deadtime(): Done.\n");	
	
	return ttotal;
}

/* function for playing GRE rf1 pulse */
int play_rf1(float phs) {
	int ttotal = 0;

	/* set rx and tx phase */
	setphase(phs, &rf1, 0);
	setphase(phs, &echo1, 0);

	/* Play the rf1 */
	fprintf(stderr, "\tplay_rf1(): playing rf1core (%d us)...\n", dur_rf1core);
	boffset(off_rf1core);
	startseq(0, MAY_PAUSE);
	settrigger(TRIG_INTERN, 0);
	ttotal += dur_rf1core + TIMESSI;

	return ttotal;	
}

/* function for playing the acquisition window */
int play_readout(int zero_grads) {
	int ttotal = 0;
	fprintf(stderr, "\tplay_readout(): playing seqcore (%d us)...\n", dur_seqcore);

	if (zero_grads) {
		setiamp(0, &gxw, 0);
		setiamp(0, &gyw, 0);
		setiamp(0, &gzw, 0);
	}

	/* play the seqcore */
	boffset(off_seqcore);
	startseq(0, MAY_PAUSE);
	settrigger(TRIG_INTERN, 0);
	ttotal += dur_seqcore + TIMESSI; 
		
	setiamp(MAX_PG_WAMP, &gxw, 0);
	setiamp(MAX_PG_WAMP, &gyw, 0);
	setiamp(MAX_PG_WAMP, &gzw, 0);

	fprintf(stderr, "\tplay_readout(): Done.\n");
	return ttotal;
}

/* function for sending endpass packet at end of sequence */
STATUS play_endscan() {
	fprintf(stderr, "\tplay_endscan(): sending endpass packet...\n");
	
	/* send SSP packet to end scan */
	boffset( off_pass );
	setwamp(SSPD + DABPASS + DABSCAN, &endpass, 2);
	settrigger(TRIG_INTERN, 0);
	startseq(0, MAY_PAUSE);  

	fprintf(stderr, "\tplay_endscan(): Done.\n");
	return SUCCESS;
}

/* function for playing prescan sequence */
STATUS prescanCore() {

	for (view = 0; view < rspvus; view++) {
		/* initialize phase */
		fprintf(stderr, "prescanCore(): playing flip pulse for prescan tr %d...\n", view);
		play_rf1(0);

		/* Load the DAB */	
		fprintf(stderr, "prescanCore(): loaddab(&echo1, 0, 0, 0, %d, DABON, PSD_LOAD_DAB_ALL)...\n", view);
		loaddab(&echo1, 0, 0, 0, view, DABON, PSD_LOAD_DAB_ALL);

		fprintf(stderr, "prescanCore(): playing readout for prescan iteration %d...\n", view);
		play_readout(1);

		fprintf(stderr, "prescanCore(): playing deadtime for prescan iteration %d...\n", view);
		play_deadtime(1s);
	}

	rspexit();

	return SUCCESS;
}

/* For Prescan: MPS2 Function */
STATUS mps2( void )
{
	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	rspent = L_MPS2;
	rspvus = 30000;
	prescanCore();
	rspexit();

	return SUCCESS;
}   /* end mps2() */


/* For Prescan: APS2 Function */
STATUS aps2( void )
{   
	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	rspent = L_APS2;
	rspvus = 1026;
	prescanCore();
	rspexit();

	return SUCCESS;
}   /* end aps2() */

STATUS scan( void )
{ 

	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	int ttotal = 0;
	fprintf(stderr, "scan(): beginning scan (t = %d / %.0f us)...\n", ttotal, pitscan);	

	/* initialize phase */
	rfphs = 0;
	
	if (ndisdaqs == 0) { /* play an empty acquisition to reset the DAB */
		/* turn the DABOFF */
		loaddab(&echo1, 0, 0, DABSTORE, 0, DABOFF, PSD_LOAD_DAB_ALL);

		play_readout(1);
		
	}	
	else
	{
		/* loop through disdaqs */
		for (disdaqn = 0; disdaqn < ndisdaqs; disdaqn++) {

			/* play the rf1 core */	
			fprintf(stderr, "scan(): playing rf1 core for disdaq %d (t = %d / %.0f us)...\n", disdaqn, ttotal, pitscan);
			ttotal += play_rf1(rfphs);
			rfphs += 117*rfspoil_flag;

			/* load the DAB */		
			fprintf(stderr, "scan(): loaddab(&echo1, %d, 0, DABSTORE, 0, DABOFF, PSD_LOAD_DAB_ALL)...\n", 0);
			loaddab(&echo1,
					0,
					0,
					DABSTORE,
					0,
					DABOFF,
					PSD_LOAD_DAB_ALL);		
			
			/* play readout */				
			fprintf(stderr, "scan(): playing seqcore for disdaq %d (%d us)...\n", disdaqn, dur_seqcore);
			play_readout(1);
		}
	}

	/* loop through rotations and shots */
	for (rotn = 0; rotn < nrots; rotn++) {
		for (shotn = 0; shotn < nshots; shotn++) {
		
			/* play the rf1 core */	
			fprintf(stderr, "scan(): playing rf1core for rot %d, shot %d (t = %d / %.0f us)...\n", rotn, shotn, ttotal, pitscan);
			ttotal += play_rf1(rfphs);
			rfphs += 117*rfspoil_flag;

			/* load the DAB */
			slice = rotn + 1;
			view = 	shotn + 1;
			echo = 0;
			fprintf(stderr, "scan(): loaddab(&echo1, %d, %d, DABSTORE, %d, DABON, PSD_LOAD_DAB_ALL)...\n", slice, echo, view);
			loaddab(&echo1,
				slice,
				echo,
				DABSTORE,
				view,
				DABON,
				PSD_LOAD_DAB_ALL);		
					
			/* play the readout */
			fprintf(stderr, "scan(): playing seqcore for rot %d, shot %d (t = %d / %.0f us)...\n", rotn, shotn, ttotal, pitscan);
			ttotal += play_readout(0);

		}
	}

	fprintf(stderr, "scan(): reached end of scan, sending endpass packet (t = %d / %.0f us)...\n", ttotal, pitscan);
	play_endscan();

	rspexit();

	return SUCCESS;
}


/********************************************
 * dummylinks
 *
 * This routine just pulls in routines from
 * the archive files by making a dummy call.
 ********************************************/
void dummylinks( void )
{
	epic_loadcvs( "thefile" );            /* for downloading CVs */
}


@host
/******************************************************
* Define the functions that will run on the host 
* during predownload operations
*****************************************************/

int genrots() {

	/* Declare values and matrices */
	int n;
	nrots = 1;

        /* Get original transformation matrix */
        for (n = 0; n < 9; n++) rotmattbl[0][n] = (float)rsprot[0][n] / MAX_PG_WAMP;

	return 1;
};

float calc_sinc_B1(float cyc_rf, int pw_rf, float flip_rf) {

	int M = 1001;
	int n;
	float w[M], x[M];
	float area = 0.0;

	/* Create an M-point symmetrical Hamming window */
	for (n = 0; n < M; n++) {
		w[n] = 0.54 - 0.46*cos( 2*M_PI*n / (M-1) );
	}	

	/* Create a sinc pulse */
	for (n = -(M-1)/2; n < (M-1)/2 + 1; n++) {
		if (n == 0)
			x[n + (M-1)/2] = 1.0;
		else
			x[n + (M-1)/2] = sin( 4 * M_PI * cyc_rf * n / (M-1) ) / ( 4 * M_PI * cyc_rf * n / (M-1) );
	}
	
	/* Calculate the area (abswidth) */
	for (n = 0; n < M; n++) {
		area += x[n] * w[n] / M;
	}

	/* Return the B1 (derived from eq. 1 on page 2-31 in EPIC manual) */
	return (SAR_ASINC1/area * 3200/pw_rf * flip_rf/90.0 * MAX_B1_SINC1_90);
}

float calc_hard_B1(int pw_rf, float flip_rf) {
	return (flip_rf / 180.0 * M_PI / GAMMA / (float)(pw_rf*1e-6));
}

int write_scan_info() {

	FILE *finfo = fopen("scaninfo.txt","w");
	fprintf(finfo, "Rx parameters:\n");
	fprintf(finfo, "\t%-50s%20f %s\n", "X/Y FOV:", (float)opfov/10.0, "cm");
	fprintf(finfo, "\t%-50s%20f %s\n", "3D slab thickness:", (float)opslquant*opslthick/10.0, "cm"); 	

	fprintf(finfo, "Hardware limits:\n");
	fprintf(finfo, "\t%-50s%20f %s\n", "Max gradient amplitude:", gmax, "G/cm");
	fprintf(finfo, "\t%-50s%20f %s\n", "Max slew rate:", smax, "G/cm/s");

	fprintf(finfo, "Readout parameters:\n");
	fprintf(finfo, "\t%-50s%20d\n", "SPARKY trajectory ID #:", trajid);
	fprintf(finfo, "\t%-50s%20f %s\n", "Flip angle:", opflip, "deg");
	fprintf(finfo, "\t%-50s%20f %s\n", "Echo time:", (float)opte*1e-3, "ms");
	fprintf(finfo, "\t%-50s%20s\n", "RF phase spoiling:", (rfspoil_flag) ? ("on") : ("off"));	
	fprintf(finfo, "\t%-50s%20f %s\n", "Shot interval (TR):", (float)optr*1e-3, "ms");
	fprintf(finfo, "\t%-50s%20d\n", "Number of disdaqs:", ndisdaqs);
	fprintf(finfo, "\t%-50s%20f %s\n", "Crusher area factor:", crushfac, "% kmax");
	fprintf(finfo, "\t%-50s%20f %s\n", "Acquisition window duration:", acq_len*GRAD_UPDATE_TIME*1e-3, "ms");
	fprintf(finfo, "\t%-50s%20s\n", "Fat suppression:", (fatsup_flag) ? ("on") : ("off"));	

	fclose(finfo);
	return 1;
}

/************************ END OF UMVSASL.E ******************************/

