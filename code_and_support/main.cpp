// Version 7.0.1 June, 2022 (built from Arvert 6.1.3 then Arvert 7.0.0 with Toffset)
// Peter Zeitler, Lehigh University

// This version includes bug fixes and improvements that mean it should be used instead of older versions (v6 and earlier).

#include "arvert.h"
#include "ZRDAAM.h"
#include  <sys/stat.h>

#define SIZE FILENAME_MAX
/*
	Unfortunately, for legacy reasons this code uses a mix of zero and one-based arrays, so be careful if you need to tinker!
	
	Following lines highlight just a few bug fixes and updates important to understanding the current version.

	15-16 November 2014:
	Got code running under C++ with changes to array allocation, few other tweaks to variable names
	to behave with RDAAM variables etc.

	29-XX November 2014:
	Yikes: went and ripped all the guts out of the motor, which is now lying on the road, Mongolian 
	style: we're going all-in: adding RDAAM/ZRDAAM, adding multiple minerals to inversion, AND moving 
	away from global variables, finally. Really a stupid thing to do - this is bad programming style: hope we can put 
	Humpty-Dumpty back together again! This will become version 5.0.0, since it involves changes to reheating, C++, RDAAM, multiple minerals, and removal 
	of global variables.

	11.24.2015:
	Now use just one diffusion_precision variable (shared by RDAAM and VD)

	~4 May 2016:
	Randomized generation of time nodes at MC stage. These are then interpolated into larger number of 
	evenly spaced nodes for CRS stage. This reduces problem where time node location acts as an implicit constraint, 
	since the location time nodes determines when a change in rate can occur.
	
	October 2016:
	Arvert 6.0.0 -- moved a number of macro parameters from <arvert.h> to <crs.in> input file, as declared variables in 
	main() and in routines. This should reduce need to recompile; it did mean lots of editing to subroutine calls and 
	function prototypes.

	October 2017:
	Version 6.1.1. Small cosmetic fixes and tweaks. No change to functionality.
	
	March 2018: Kalin reports odd convergence with poor age-spectrum fits for samples with high numbers of CRS time 
	nodes. I can't reproduce this, so for now, blame the user.

	November 2018:
	Version 6.1.2. Small bug fix: several parameters (flips, nbest, nworst, early_loss, early_weight) had been defined 
	as const and initialized; this meant values read from file were ignored. Removed const typing and adjusted function 
	calls and prototypes accordingly.

	May 9, 2020: 
	Arvert 6.1.3: cleaned up comments, tidied up location of some blocks of code, and in function MonteTt slightly 
	changed the criteria for random times being too close to start or finish, from 1% to 2%.
	
	December 28, 2021 - early January, 2022  Started Arvert 7.0.0: some bug fixes and plotting changes but main thing is
	additional of handling multiple MDD + mineral samples what can be separated by a temperature offset.
	
	January 19, 2022. Code complete: compiling, writing plot files, plot script working well. Time for testing.
	
	January 23, 2022. Code is working and plotting, for different combos of minerals and/or spectra. Built out gmtplot
	option so there is a sleeker faster option for plotting (adding all data really slows script). Will now test
	the code using the Kohistan data set.
	
	January 25-27: tracked down bugs in MontetT() Arvert 6 was passing wrong values because it filled upperlimits[] and
	lowerlimits[] with ttpoints not rttpoints number of time node info. Those arrays need to be made inside MontetT
	because each MC history has different time nodes and so different limiting arrays. Also tracked down bugs that 
	were not forcing first T data to be safe for lovera(). Jeez!

	15 June: final tweak: decoupled post processing from nbest, so user can work with more histories but still just display a smaller best set
*/
// ***************************************************************************************
// ******************* Declare variables *************************************************
// ***************************************************************************************

// NB - Many but not all lines of code that might present an issue with porting to other compilers are flagged
// with the comment PORT-ISSUE
	
int main (void) {

// Main variables (used to be global in days of yore)
	long nsamples,labsteps[MAXSAMPLES],ndomains[MAXSAMPLES],goalsteps[MAXSAMPLES],nspectra[MAXSAMPLES];
	double temperaturelab[MAXSAMPLES][MAXLABSTEPS],timelab[MAXSAMPLES][MAXLABSTEPS];
	double Eact[MAXSAMPLES][MAXDOMAINS],D0[MAXSAMPLES][MAXDOMAINS],fraction[MAXSAMPLES][MAXDOMAINS];
	double goal_loss[MAXSAMPLES][MAXLABSTEPS],goal_age[MAXSAMPLES][MAXLABSTEPS],goal_error[MAXSAMPLES][MAXLABSTEPS];
	long goal_skip[MAXSAMPLES][MAXLABSTEPS];
	char info[MAXLINE],suffix[20],isuffix[10],filename[50],finalfilename[55],recfilename[55],extent[6];	
	double timein[NODES],temperaturein[NODES], offset_temperaturein[NODES]; //T-t history, for use by Monte and CRS routines when talking to lovera()
	double MonteTime[NODES], MonteTemp[NODES];
	double Time_at_node[NODES];
	time_t RanSeedZero;
	double lowerlimit[NODES], upperlimit[NODES], lowerlimit_mc[NODES], upperlimit_mc[NODES];
	double ModelFit[MAXPOOLSIZE];
	long WorstOne,BestOne;

	double early_loss[MAXSAMPLES], early_weight[MAXSAMPLES];
	long nbest=15, nworst=15, flips=1;
	long post_set = 50;

	long restart;
	long npostprocess;
	double prange;

	double bgeom,slab, geomfac, Tc_largest;
	
	double modelduration;
	long nconstraints;
	double constraintTime[NODES], constraintTempmin[NODES], constraintTempmax[NODES];
	double heatRate, coolRate;
	double monteHeatRate,monteCoolRate;
	long rttpoints, ttpoints;
	long CRSiterations, attempts_CRS;
	double amplify, fitvalue;
	long geometry;
	long fitchoice;
	long poolsize, subsetsize, CRS_selection;
	double bulkloss39[MAXSAMPLES][MAXLABSTEPS], age[MAXSAMPLES][MAXLABSTEPS];
	long CRSselections[MAXPOOLSIZE];
	double ichange, dt;
	long fullreport;
	double maxgoodgoalage, mingoodgoalage, mingoodgoalloss, maxgoodgoalloss;
	double min_mineralage, max_mineralage;
	long goodHistory, CRSattempts;
	double avTemp[NODES], avAge[MAXSAMPLES][MAXLABSTEPS], lowEnvelope[NODES], highEnvelope[NODES];
	long gmtplot;

	long i, j, ijj, row, imonte, iCRS, imineral, n;
	long totCRS, nowCRS;
	
	long duration;
	double elapsed, runrate, elapsed2, runrate2;
	double bestCRS, worstCRS;

// for Toffset being part of CRS_selection

	double Toffset[MAXPOOLSIZE][MAXSAMPLES];
	double newToffset[MAXSAMPLES];
	double offset_min, offset_max;
	long offset_flag;

// for mineral-age machinations...
	double mineral_age[MAXSAMPLES][MAX_MINERALS], mineral_goal[MAXSAMPLES][MAX_MINERALS], mineral_goalerr[MAXSAMPLES][MAX_MINERALS];
	double mineral_eact[MAXSAMPLES][MAX_MINERALS], mineral_difzero[MAXSAMPLES][MAX_MINERALS];
	double mineral_radius[MAXSAMPLES][MAX_MINERALS];
	long mineral[MAXSAMPLES][MAX_MINERALS], use_minerals, nminerals[MAXSAMPLES];
	double mineral_weight[MAXSAMPLES][MAX_MINERALS], mineral_pool[MAXSAMPLES][MAXPOOLSIZE][MAX_MINERALS];
	double average_mineral_age[MAXSAMPLES][MAX_MINERALS];
	double mineral_U[MAXSAMPLES][MAX_MINERALS], mineral_Th[MAXSAMPLES][MAX_MINERALS], mineral_Sm[MAXSAMPLES][MAX_MINERALS];
	long diffusion_precision;
	long ndL, ndt;  // number of distance and time steps for volume-diffusion calculations
	

// Large arrays (not bothering with dynamic memory management)
	
// NOTE: 0 based, so need one extra to match old 1-based algorithms, plus one more to hold newest trial CRS history and its ages
// Thus, if actual pool size is 300 we need 302 elements in array that has members from 0 to 301

	double poolAge[MAXSAMPLES][MAXPOOLSIZE][MAXLABSTEPS];
	double poolTemp[MAXPOOLSIZE][NODES];
	double mc_Temp[MAXPOOLSIZE][NODES];	
	double poolTime[MAXPOOLSIZE][NODES];

	printf("\n");

// for recording execution speed and such
	time_t now;

// variables for talking with RDAAM / ZRDAAM routines

	TTPath path;
	TTPathPoint mypoint;
	double Heage, uncorr_age;
	double hedummy;
	bool optimize = true;

// some file-handling declarations and operations

	FILE *ifp, *ofp, *ofpscratch;
	char newFolder[SIZE], inFolder[SIZE]; // PORT-ISSUE 
	
// dummy input variable to help make input files more legible; taelk to ZRDDAM

	char dummy[50]; // 50 should be enough unless options get complex

	
// used in file I/O and examining inputs read from file
	char shorten[30];
	char line[MAXLINE];
	long proceed;

//
// ----------start doing things: find out where we are and where we are working

	getcwd( inFolder, SIZE ); // PORT-ISSUE pathname call should work but check implementation on your system
	// this call gets pathname of current working directory (e.g., where Arvert is running)
	
// change to directory from which arvert was launched
	chdir(inFolder);

// ***************************************************************************************
// ******************* Section 1. Startup message ****************************************
// ***************************************************************************************
	now = time(NULL);  	
	printf("\n");
	printf("-------------------------------------------------------------------------------\n");
	printf("                           Arvert 7.0.1 beta\n");
	printf("                            14 June, 2022\n");
	printf("\n");
	printf("                 Inversion of multiple K-feldspar age spectra\n");
	printf("                         and associated mineral ages\n");	
	printf("                   using CRS, Lovera, and RDAAM algorithms\n");
	printf("                     with randomized Mont-Carlo time nodes \n");
	printf("                        and optional post-CRS processing \n");		
	printf("\n");
	printf("                              Peter Zeitler\n");
	printf("                            Lehigh University\n");
	printf("-------------------------------------------------------------------------------\n");
	printf("                          %s",ctime(&now));
	printf("-------------------------------------------------------------------------------\n");
	
// ***************************************************************************************
// ******************* Section 2. Read data from input files *****************************
// ***************************************************************************************	
/* 		
		This section gets the following input data:
   
		(1) file domains.in (partly equivalent to Lovera's ages_me.in) 
		    which has samples' diffusion and domain parameters
		(2) file goal.in, which has the measured age spectra and flags
		    for whether or not to fit particular steps
		(3) file crs.in, which has parameters needed for the CRS algorithm and Arvert
		(4) file mineral.in, which has parameters for the mineral-age samples
*/

//  ------------ Get inversion parameters from file <crs.in> -------------

	ifp = fopen("crs.in", "r");

	if (ifp==NULL)
	{
		printf("\n FATAL ERROR: could not find (or read) file crs.in \n\n");
		exit(1);
	}

	fgets(info, MAXLINE, ifp);
	info[strlen(info) - 1] = '\0';
	
	printf("             RUN INFO:   %s\n", info);

	fgets(line, MAXLINE, ifp);
	strcpy(suffix,&line[7]);
	suffix[strlen(suffix) - 1] = '\0';
	printf("          FILE SUFFIX:   %s\n", suffix);	
	
	// ------------- Open output file to record sample and run info -------------------

	strcpy(filename,"ModelInfo.");

	now = time(NULL);
    
	fscanf(ifp, "%s%ld", &dummy,  &nsamples);
	printf("    NUMBER of SAMPLES:   %ld\n", nsamples);
	
	for (n=1; n<= nsamples; n++)
	{
		fscanf(ifp, "%s%lf", &dummy,  &Toffset[1][n]);
		if (Toffset[1][1] != 0.0)
		{
			Toffset[1][1] = 0.0;  // force first (master) sample to have zero Toffset
			printf("  !!!!!!! WARNING:   First Toffset set to zero!!");
		}
		printf("      %ld  Toffset (˚C):   %5.1f\n",n, Toffset[1][n]);	
	}
	
	fscanf(ifp, "%s%lf%lf", &dummy,  &offset_min, &offset_max);
	printf("           TOFFSET MIN, MAX:   %5.1f\t%5.1f\n", offset_min,offset_max);

	fscanf(ifp, "%s%ld", &dummy,  &offset_flag);
	assert ((offset_flag >= 0) || (offset_flag < 3 ));
	switch (offset_flag)
	{
		case 0:
			printf("     T_OFFSET TYPE?:   0 - fixed to inputs\n");			
			break;
		case 1:
			printf("     T_OFFSET TYPE?:   1 - monotonic increasing\n");		
			break;
		case 2:
			printf("     T_OFFSET TYPE?:   2 - random around reference\n");			
	}

	fscanf(ifp, "%s%lf", &dummy,  &modelduration);
	printf("\nMODEL DURATION (m.y.):   %6.1f\n", modelduration);

	assert ((modelduration > 0) && (modelduration < 4500)); //make sure times are in m.y.	
	fscanf(ifp, "%s%ld%ld", &dummy,  &rttpoints,  &ttpoints);
	assert ((ttpoints >= 5) && (ttpoints <= NODES)); // make sure arrays stay in bounds and have minimal number of points!
	assert ((rttpoints >= 3) && (rttpoints <= 20)); // make sure arrays stay in bounds and have minimal number of points!
	assert (rttpoints < ttpoints); // make sure arrays stay in bounds and have minimal number of points!
	
	printf("   MC, CRS TIME NODES:    %ld   %ld\n", rttpoints, ttpoints);
	printf("\n");	

	fscanf(ifp, "%s%ld", &dummy,  &nconstraints);
	assert ((nconstraints >= 2) && (nconstraints < 11)); // make sure arrays stay in bounds, have minimal constraints
	printf("CONSTRAINING BRACKETS tT:   %ld\n", nconstraints);
	
	printf("   TIME    TMIN   TMAX\n");  // explicit temperature constraints
	for (j = 1; j <= nconstraints; j++)
	{
		fscanf(ifp, "%s%lf%lf%lf", &dummy,  &constraintTime[j], &constraintTempmin[j], &constraintTempmax[j]);
		printf(" %6.1f   %5.1f   %5.1f\n", &dummy,  constraintTime[j], constraintTempmin[j], constraintTempmax[j]);
		assert(constraintTempmin[j] < constraintTempmax[j]);
		assert ((constraintTime[j] <= modelduration) && (constraintTime[j] >= 0.0));
		if (j == 1)
		{
			assert(constraintTime[j] == modelduration);  // first contraint has to be at start
		}
		if (j > 1)
		{
			assert(constraintTime[j] < constraintTime[j-1]); // each successive constraint must be younger
		}
		if (j == nconstraints)
		{
			assert(constraintTime[j] == 0.0); // last constraint has to be at end
		}
	}
	printf("\n");
	
	fscanf(ifp, "%s%lf%lf", &dummy,  &monteHeatRate,&monteCoolRate); // implicit rate constraints for monte-carlo histories
	printf("    MAX MONTE-CARLO HEATING RATE: %5.1f   MAX MONTE-CARLO COOLING RATE: %5.1f\n", monteHeatRate,monteCoolRate);
	assert ((monteHeatRate >= 0.0) && (monteHeatRate <= 1000)); //make sure rates are not negative or absurdly high
	assert ((monteCoolRate >= 0.0) && (monteCoolRate <= 1000));
	
	fscanf(ifp, "%s%lf%lf", &dummy,  &heatRate,&coolRate);  // implicit rate constraints for CRS histories
	printf("            MAX CRS HEATING RATE: %5.1f           MAX CRS COOLING RATE: %5.1f\n", heatRate,coolRate);
	assert ((heatRate >= 0.0) && (heatRate <= 1000)); //make sure rates are not negative or absurdly high
	assert ((coolRate >= 0.0) && (coolRate <= 1000));
	
	fscanf(ifp, "%s%lf", &dummy,  &amplify);
	printf("        CRS AMPLIFICATION FACTOR:   %4.2f\n", amplify);
	assert ((amplify >= 0.5) && (amplify <= 3)); //make sure amplify is not completely wildly absurd
	
	fscanf(ifp, "%s%ld%ld", &dummy,  &subsetsize, &poolsize);
	printf("          SUBSET SIZE, POOL SIZE:   %ld   %ld\n", subsetsize, poolsize);
	assert(poolsize <= 300);
	assert ((subsetsize < 101) && (subsetsize <= floor(poolsize/3)) && (subsetsize >= MINSUBSET)); //make sure subset size is not stupid
	printf("\n");
	
	fscanf(ifp, "%s%lf", &dummy,  &fitvalue);
	printf("               FITTING CRITERION:  %5.2f\n", fitvalue);
	assert (fitvalue >= 0); //make sure fit is not stupid
	
	fscanf(ifp, "%s%ld", &dummy,  &fitchoice);
	printf("                  FITTING OPTION:   %ld  ", fitchoice);
	assert ((fitchoice == 0) || (fitchoice == 1) || (fitchoice == 2)); // 0=step-weighted mean percent, 1=mean percent, 2=mswd-like parameter
	switch (fitchoice)
	{
		case 0:
			printf("(early-weighted mean percent deviation)\n");
			break;
		case 1:
			printf("(mean percent deviation)\n");
			break;
		case 2:
			printf("(mswd)\n");
	}
	
	fscanf(ifp, "%s%ld", &dummy,  &geometry);
	printf("          MDD DIFFUSION GEOMETRY:   %ld  ", geometry);
	assert ((geometry == 1) || (geometry == 2)); // 1 spheres, 2 slabs
	if (geometry == 1)
	{
		printf("(spherical)\n");
		geomfac = 55.;
	}
	else
	{
		printf("(infinite-slab)\n");
		geomfac = 8.7;
	}

	fscanf(ifp, "%s%ld", &dummy,  &attempts_CRS);
	printf("                MAX CRS ATTEMPTS:   %ld\n", attempts_CRS);

	fscanf(ifp, "%s%lf", &dummy,  &dt);
	printf("DISCRETIZATION DELTA-TEMPERATURE: %5.1f\n", dt);
	assert ((dt >= 1.0) && (dt <= 20)); //make sure dt value is reasonable
	
	fscanf(ifp, "%s%lf", &dummy,  &ichange);
	printf("\n           LOVERA SERIES CUT-OFF:   %5.1e\n", ichange);
	assert ((ichange <= 0.001) && (dt >= 1.0e-8)); //make sure ichange value is reasonable

	fscanf(ifp, "%s%ld", &dummy,  &fullreport);
	printf("      FLAG TO WRITE FULL REPORTS:   %ld\n", fullreport);  // 0 means only rolling-progress reports
	assert ((fullreport == 0) || (fullreport == 1)); //make sure fullreport value is reasonable	

	fscanf(ifp, "%s%ld", &dummy, &gmtplot);	
	
	switch (gmtplot)
	{
		case 1: 
			printf("           PLOT RESULTS WITH GMT:   %ld (call script, simple plot)\n", gmtplot);		
		break;
		case 2:
			printf("           PLOT RESULTS WITH GMT:   %ld (call script, full plot)\n", gmtplot);		
		break;
		case 0:
			printf("           PLOT RESULTS WITH GMT:   %ld (no script, no plot)\n", gmtplot);		
		break;
		default:
			printf("           PLOT RESULTS WITH GMT:     %ld (bad value, bad read)\n", gmtplot);
			printf("\n\n >>>> ERROR IN READING DATA; run terminated (check file crs.in)\n\n", gmtplot);			
			exit(1);
		break;
	}
	
	fscanf(ifp, "%s%ld", &dummy,  &nbest);
	printf("     # OF BEST HISTORIES TO PLOT:  %ld\n", nbest); 
	assert ((nbest > 0) || (nbest < poolsize)); // make sure nbest value is reasonable	
	
	fscanf(ifp, "%s%ld", &dummy,  &nworst);
	printf("    # OF WORST HISTORIES TO PLOT:  %ld\n", nworst);
	assert ((nworst > 0) || (nworst < poolsize)); // make sure nworst value is reasonable	
	
	fscanf(ifp, "%s%ld", &dummy,  &flips);
	printf("      MAX PERMITTED tT REVERSALS:   %ld\n", flips); 
	assert ((flips >= 0) || (flips < 10)); // make sure flips value is not insane	
	
	printf("\n");
	for(n=1; n <= nsamples; n++)
	{
		fscanf(ifp, "%s%lf%lf", &dummy,  &early_loss[n], &early_weight[n]);
		printf("     SAMPLE %ld EARLY-STEP, WEIGHT:   %5.3f\t%5.3f\n", n,early_loss[n],early_weight[n]); 
		assert ((early_loss[n] > 0.0) && (early_loss[n] < 1.0)); // make sure early_loss  is reasonable	
		assert (early_weight[n] > 0.0); // make sure early_weight is positive
	}
	printf("\n");

	fscanf(ifp, "%s%ld", &dummy,  &use_minerals);
	assert ((use_minerals >= 0) || (use_minerals < 3 ));	
	if (use_minerals == 1)
	{
		printf(" CAN MINERAL AGES BE CONSTRAINTS:   1 - yes\n");
	}
	if (use_minerals == 0)
	{
		printf(" CAN MINERAL AGES BE CONSTRAINTS:   0 - no\n");
	}
	if (use_minerals == 2)
	{
		printf(" CAN MINERAL AGES BE CONSTRAINTS:   2 - yes (and age spectra ignored!!)\n");
	}
		printf("\n");

	fscanf(ifp, "%s%ld", &dummy,  &CRSiterations);
	printf("                  CRS ITERATIONS:   %ld\n", CRSiterations);
	
	fscanf(ifp, "%s%ld", &dummy,  &npostprocess);
	printf("      POST-PROCESSING ITERATIONS:   %ld\n", npostprocess);
	
	fscanf(ifp, "%s%ld", &dummy,  &post_set);
	printf("          POST-PROCESSING SUBSET:   %ld\n", post_set);
	assert ((post_set > 1) && (post_set <= poolsize)); // make sure value makes sense	

	fscanf(ifp, "%s%lf", &dummy,  &prange);
	printf("    tT POST-PROCESSING RANGE ˚C):   %5.1f\n", prange);
	
	fscanf(ifp, "%s%ld", &dummy,  &restart);
	printf("                  RESTART OPTION:   %ld  ", restart);
	switch (restart)
	{
		case 0 :
			printf("(new start from Monte-Carlo histories)\n");
			printf("\n");
			break;
		case 1 :
			printf("(restart CRS processing) \n");
			printf("\n");
			break;
		case 2 :
			printf("(restart for post-processing\n");
			printf("\n");
			break;
	}

	fclose(ifp);
// get mineral-age information information if we are using mineral ages

	if (use_minerals > 0)
	{
		printf("\n----------------------------------- MINERAL-AGE INFO -----------------------------------\n");

		ifp=fopen("minerals.in","r");

		if (ifp==NULL)
		{
			printf("\n FATAL ERROR: could not find (or read) file mineral.in \n\n");
			exit(1);
		}

		for (n=1; n <= nsamples; n++)
		{
			fscanf(ifp, "%s%ld", &dummy, &nminerals[n]);
			printf(">>> SAMPLE %ld NUMBER OF MINERAL AGES:   %ld\n",n, nminerals[n]);
			
			for	(i = 1; i <= nminerals[n]; i++)
			{
				if (nminerals[n] > 0)
				{
					fscanf(ifp, "%s%ld", &dummy, &mineral[n][i]);
					assert ((mineral[n][i] >= 0) && (mineral[n][i] < 10 ));	
					switch (mineral[n][i])
					{
						case 0:
							printf("              ******** MINERAL %ld:   0 - apatite (volume diffusion)\n",i);
							break;
						case 1:
							printf("              ******** MINERAL %ld:   1 - zircon (volume diffusion)\n",i);
							break;
						case 2:
							printf("              ******** MINERAL %ld:   2 - titanite (volume diffusion)\n",i);
							break;
						case 3:
							printf("              ******** MINERAL %ld:   3 - other (no alpha loss)\n",i);
							break;
						case 4:
							printf("              ******** MINERAL %ld:   4 - apatite (RDAAM)\n",i);
							break;
						case 5:
							printf("              ******** MINERAL %ld:   5 - zircon (ZRDAAM)\n",i);
							break;			
					}

					fscanf(ifp, "%s%lf%lf", &dummy,  &mineral_goal[n][i],&mineral_goalerr[n][i]);
					printf("              GOAL MINERAL AGE %ld:  %6.2f Ma     ERROR IN AGE: %5.2f Ma\n",i, mineral_goal[n][i], mineral_goalerr[n][i]);
					assert ((mineral_goalerr[n][i] > 0.0) && (mineral_goal[n][i] >= 0.0)); 
	
					fscanf(ifp, "%s%lf", &dummy,  &mineral_weight[n][i]);
					printf("  MINERAL-AGE WEIGHTING FACTOR %ld:  %6.2f\n",i, mineral_weight[n][i]);
					assert ((mineral_weight[n][i] >= 0.0) && (mineral_eact[n][i] <= 500));	

					fscanf(ifp, "%s%lf", &dummy,  &mineral_radius[n][i]);
					printf("              EFFECTIVE RADIUS %ld: %7.2f microns\n",i, mineral_radius[n][i]);
					assert ((mineral_radius[n][i] >= 25.0) && (mineral_radius[n][i] <= 1000.0 ));
	
					fscanf(ifp, "%s%lf%lf", &dummy,  &mineral_eact[n][i], &mineral_difzero[n][i]);
					printf("             ACTIVATION ENERGY %ld:   %5.2f kcal/mol     DIFZERO:  %9.3e cm2/sec\n",i, mineral_eact[n][i], mineral_difzero[n][i]);
					assert ((mineral_eact[n][i] >= 20.0) && (mineral_eact[n][i] <= 60.0 ));
					assert(mineral_difzero[n][i] > 0);
		
					fscanf(ifp, "%s%lf%lf%lf", &dummy, &mineral_U[n][i], &mineral_Th[n][i], &mineral_Sm[n][i]);
					printf("         PARENT CONCENTRATIONS %ld:   %8.2f ppm U  %8.2f ppm Th  %8.2f ppm Sm\n\n",i, mineral_U[n][i], mineral_Th[n][i], mineral_Sm[n][i]);
					assert ((mineral_U[n][i] >= 0.0) && (mineral_U[n][i] <= 100000.0));		
					assert ((mineral_Th[n][i] >= 0.0) && (mineral_Th[n][i] <= 100000.0));
					assert ((mineral_Sm[n][i] >= 0.0) && (mineral_Sm[n][i] <= 100000.0));	
				} // end if-nminerals[]
			} // end of read-mineral loop
		} // end of n-loop through samples
	
		fscanf(ifp, "%s%ld", &dummy, &diffusion_precision);
		assert ((diffusion_precision >= 0) || (diffusion_precision < 3 ));	
		if (diffusion_precision == 0)
		{
			printf("  DIFFUSION and RDAAM PRECISION?:    0 - good \n");
		}
		if (diffusion_precision == 1)
		{
			printf("  DIFFUSION and RDAAM PRECISION?:    1 - better\n");
		}
		if (diffusion_precision == 2)
		{
			printf("  DIFFUSION and RDAAM PRECISION?:    2 - best\n");
		}
		
		fclose(ifp);
	} // end if_use_minerals clause
	printf("\n");
	
//  ------------ Get diffusion parameters and goal spectra from files <domains.in> and <goal.in> if we have mdd data -------------

	if (use_minerals < 2)
	{
		double Tc_max[MAXSAMPLES];

		ifp = fopen("domains.in", "r");
	
		if (ifp==NULL)
		{
			printf("\n FATAL ERROR: could not find (or read) file domains.in \n\n");
			exit(1);
		}
		for (n=1; n <= nsamples; n++)
		{	
			fscanf(ifp, "%ld", &ndomains[n]);   // number of domains
			if (ndomains[n] < 1)
			{
				nspectra[n] = 0;
			}
			else
			{
				nspectra[n] = 1;			
			}
			assert(ndomains[n] >= 0 && ndomains[n] < MAXDOMAINS);
			for (j = 1; j <= ndomains[n]; j++)
			{
				fscanf(ifp, "%lf%lf%lf", &Eact[n][j], &D0[n][j], &fraction[n][j]);
				D0[n][j] = exp(2.3025851 * D0[n][j]);  // convert log10 D0 to true values
				D0[n][j] = D0[n][j] * 3.1557600e13;    // 1/sec to 1/m.y.
			}
			if (ndomains[n] > 0)
			{
				Tc_max[n] = 473.; // get TC10 for largest domain in each sampple
				for (i = 1;i <= 20;i++)
				{
					Tc_max[n] = Eact[n][ndomains[n]] * 1000./ 1.9859 / log(geomfac * Tc_max[n] * Tc_max[n] * pow(10.,log(D0[n][ndomains[n]]/3.1557600e13)/2.3025851)/(Eact[n][ndomains[n]] * 1000./ 1.9859 * 10./(3.15576e13)));
				}				
				Tc_max[n] = Tc_max[n] - 273.15;	
			}
			else
			{
				Tc_max[n] = 0.0;			
			}
		}  // end of n-loop through samples
		fclose(ifp);


		for (n=1; n <= nsamples; n++)
		{	
			if(Tc_max[n] > Tc_largest)
			{
				Tc_largest = Tc_max[n];
			}
		}
		
		Tc_largest = Tc_largest + 200.; //hardwire temperature headroom above largest Tc
		
		//  ------------ Get goal spectra

		ifp = fopen("goal.in", "r");
	
		if (ifp==NULL)
		{
			printf("\n FATAL ERROR: could not find (or read) file goal.in \n\n");
			exit(1);
		}	

		for (n=1; n <= nsamples; n++)
		{
			if (nspectra[n] > 0)
			{
				fscanf(ifp, "%ld", &goalsteps[n]);   // number of goal-spectrum steps
				labsteps[n] = goalsteps[n];  // kludge for retro compatibility if we change heatsched()
				assert(goalsteps[n] > 0 && goalsteps[n] <= MAXLABSTEPS);
				for (j = 1; j <= goalsteps[n]; j++)
				{
					fscanf(ifp, "%lf%lf%lf%ld", &goal_loss[n][j], &goal_age[n][j], &goal_error[n][j], &goal_skip[n][j]);
	//DEBUG				printf("\t%lf\t%lf\t%lf\t%ld\n",goal_loss[n][j], goal_age[n][j], goal_error[n][j], goal_skip[n][j]);
					assert ((goal_skip[n][j] == 0) || (goal_skip[n][j] == 1 ));	
				}
				goal_loss[n][0] = 0.0; // set this to make it easier to calculate step losses with simple loop
			}
		} // end n-loop through samples
		fclose(ifp);
	}

// ******** Give user a chance to review inputs before launching model

	printf("-------------------------------------------------------------------------------\n");
	printf("Values OK?\n");
	printf("  0 -- Proceed\n");
	printf("  1 -- Abort model  --> ");
	scanf("%ld", &proceed);
	printf("-------------------------------------------------------------------------------\n");
	if (proceed == 0)
	{
		printf("\n");
		printf(">>> cntrl-C to abort <<<\n");
		printf("\n");
		printf("Model starting...\n");
	}
	else
	{
		printf("\nMODEL ABANDONED\n\n");
		exit(1);
	}

// ***************************************************************************************
// *******************  Section 3. Initialize various things *****************************
// ***************************************************************************************
	long step_flag;
	double dage;

// abandoned xls in favor of txt: this is on net faster to "open with" using Excel
	strcpy(extent,".txt");

// Set up MDD diffusion geometry: 1 for spheres, 2 for slabs, as used in code below
	if (geometry == 1) 
	{
			bgeom = 6.0;    // spheres (geometry flag is 1)
			slab = 1.;
	}
	else 
	{
			bgeom = 8.0;    // infinite slabs (geometry flag is 2)
			slab = 4.;
	}

	if (use_minerals > 0)
	{
	// set up distance and time steps for volume-diffusion routine used for mineral ages, based on user input
		switch (diffusion_precision)
		{
			case 0:
			{
				ndL = 100;  // quite fast; rather coarse for small losses and accurate integration
				ndt = 100;
				break;
			}
			case 1:
			{
				ndL = 250;  // ok but slower; still a bit coarse for small losses and accurate integration
				ndt = 250;			
				break;
			}		
			case 2:
			{
				ndL = 1000;   // more accurate but slower
				ndt = 500;			
				break;
			}
		}
	}

// seed for system rand(); we just use system rand() to initialize the much better random-generator Ranq1 from Numerical Recipes

	srand (time(NULL));	

// Shape constraints into the XXXerLimit[] array
// Array indices for upperlimit[] and lowerlimit[] map to those of timein[];

// Constraint at start:
// Lovera routine REQUIRES initial cooling before any reheating might occur. So even
// if there are no geologic constraints, at start we must have a temp well greater than Tc of largest domain.
// Experimentation shows this needs to be 150 to 200 C above largest Tc. 
// Earlier, determined max Tc for largest domain over all samples. Use this value, hardwired up 200 C.

	if (constraintTempmin[1] < (Tc_largest))  //force first minimum constraint to have headroom above Tc_largest
	{
		constraintTempmin[1] = (Tc_largest);
		if (constraintTempmin[1] >= constraintTempmax[1])  //make sure coercion doesn't foul max constraint
		{
			constraintTempmax[1] = constraintTempmin[1] + 50.;
		}
		printf("WARNING: initial temperature constraints raised to provide headroom for lovera()\n");
	}
	
// Get span constrained by observed ages. Originally this was for making timenodes. Now used mostly for plotting.
// For legacy reasons, time arrays are 1-based!

  // First, parse goal spectrum to get minimum and maximum goal-spectrum ages that are "good" (i.e., used in fitting).
  // Given that we generally expect spectra to rise, we can avoid messy search by just going up from start and down
  // from end and find first step flagged as "good". Only problem with this is if someone were to model a very odd
  // sample with differing activation energies in domains such that ages dropped at high release. But code should
  // still function in this case.
  
  // Goal arrays are 1-based
  
  // First deal with age spectra if we have them
  	if (use_minerals < 2)
  	{
		for (n=1; n <= nsamples; n++)
		{
			if (nspectra[n] > 0)
			{
				step_flag = -1;
				i = 0;
				while (step_flag < 0)
				{
					i++;
					assert(i <= labsteps[n]);  // trap infinite loop if user screwed up and no steps are flagged as good!
					if (goal_skip[n][i] > 0)
					{
						step_flag = 1;
						mingoodgoalage = goal_age[n][i];
						mingoodgoalloss = goal_loss[n][i]; 
					}
				}
				step_flag = -1;
				i = goalsteps[n]+1;
				while (step_flag < 0)
				{
					i--;
					assert(i > 0);  // trap infinite loop if user screwed up and no steps are flagged as good!
					if (goal_skip[n][i] > 0)
					{
						step_flag = 1;
						maxgoodgoalage = goal_age[n][i];
						maxgoodgoalloss = goal_loss[n][i]; 
					}
				}
			} // end if spectra[]
		} // end nsamples loop
	}
	else  // dummy out some values if no MDD data are being used
	{
		maxgoodgoalage = modelduration; 
		mingoodgoalage = 0.0;
	}
	// Next, deal with mineral ages if there are any
  	if (use_minerals > 0)
  	{
		min_mineralage = 10.e9; // set this silly high
		max_mineralage = 0.0; // set this silly low
  		for (n=1; n <= nsamples; n++)
		{
			if (nminerals[n] > 0)
			{
				for	(i = 1; i <= nminerals[n]; i++)
				{	
					if (mineral_goal[n][i] > max_mineralage) max_mineralage = mineral_goal[n][i];
					if (mineral_goal[n][i] < min_mineralage) min_mineralage = mineral_goal[n][i];	
				}
			}
		}	
	}
	else  // dummy out some values if no minerals are being used
	{
		min_mineralage = modelduration; 
		max_mineralage = 0.0;
	}
	
	// finally, handle case where both MDD and mineral data are in use; get final values

	switch (use_minerals)
	{
		case 0:  // MDD data only
		{
			min_mineralage = mingoodgoalage;
			max_mineralage = maxgoodgoalage;
			break;
		}
		case 1:  // both MDD and mineral data
		{
			if (min_mineralage > mingoodgoalage)  // youngest age is in spectrum
			{
				min_mineralage = mingoodgoalage;
			}
			if (max_mineralage < maxgoodgoalage)  // oldest age is in spectrum
			{
				max_mineralage = maxgoodgoalage;
			}
			break;
		}		
		case 2:  // mineral data only
		{
			mingoodgoalage = min_mineralage;
			maxgoodgoalage = max_mineralage;
			break;
		}
	}

// at this point, full range of MDD + mineral ages is spanned by min_mineralage and max_mineralage
// so we can use those for full range and still have mingoodgoalage and maxgoodgoalage as
// records of span of ages in spectra (if we have them)

// Check to see that nothing is or has become mung

	dage = max_mineralage - min_mineralage;

	if (dage <= 0.0)  //maybe some hoon is trying to invert a flat spectrum
	{
		max_mineralage = 1.10*min_mineralage;
	}
	assert(max_mineralage <= MAXHEADROOM * modelduration);  // insist on some timespan headroom in model
	assert(min_mineralage >= 0.0);
	assert(max_mineralage > min_mineralage);  // just to be 100% certain...

// Need to deal with model endpoints
// Use MonteTime[] array for this (these are communicated in call to MonteTt routine)
	MonteTime[1] = modelduration;
	MonteTime[ttpoints] = 0.0; 
	
// Remaining "internal" ttpoints are dealt with in MonteTt
	double node_spacing;
	long ispace;
	node_spacing = modelduration / (ttpoints - 1);
	for (ispace = 2; ispace < ttpoints; ispace++)
	{
		MonteTime[ispace] = modelduration - (ispace - 1) * node_spacing;
	}
assert(MonteTime[2] < modelduration); //one last check that I haven't shattered the fabric of time

// Get upper and lower temperature limits at each CRS time node
	upperlimit[1] = constraintTempmax[1];
	lowerlimit[1] = constraintTempmin[1];
	upperlimit[ttpoints] = constraintTempmax[nconstraints];
	lowerlimit[ttpoints] = constraintTempmin[nconstraints];	
	for (i = 2; i <= ttpoints - 1; i++)
	{
		j = 1;
		do
		{
			j++;
		}
		while (MonteTime[i] < constraintTime[j]);
		j = j - 1;
		upperlimit[i] = constraintTempmax[j] - (constraintTime[j] - MonteTime[i])/(constraintTime[j] - constraintTime[j+1])*(constraintTempmax[j] - constraintTempmax[j+1]);
		lowerlimit[i] = constraintTempmin[j] - (constraintTime[j] - MonteTime[i])/(constraintTime[j] - constraintTime[j+1])*(constraintTempmin[j] - constraintTempmin[j+1]);
	}
	
	printf("  ...finished constraints and parameter setup\n");
	
// ***************************************************************************************
// ********** Section 4. Calculate synthetic heating schedule(s) and 39Ar losses *********
// ***************************************************************************************
	// create a heating schedule based on domain parameters and 39-losses in measured spectrum
	if (use_minerals < 2)
	{
		heatsched(nsamples,nspectra,ndomains,labsteps,geometry,goalsteps,Eact,D0,fraction,goal_loss,timelab,temperaturelab);
		for (n=1; n<=nsamples; n++)
		{
			lablovera(n,ttpoints,geometry,bgeom,slab,ndomains[n], Eact[n], D0[n], fraction[n], labsteps[n], dt, ichange, timelab[n], temperaturelab[n], bulkloss39);
		}
	}
	
	printf("  ...finished heat schedules and 39Ar loss for MDD\n");
// ***************************************************************************************
// ******************* Section 5. Make starting pool of monte-carlo histories ************
// ***************************************************************************************

 if (restart == 0)
    {
    	totCRS = 0;  // we're starting a new model from scratch; initialize CRS iteration counter

    	ifp = fopen("CRScount.in", "w");
			fprintf(ifp, "%ld", totCRS);   // make sure counter file is zeroed for new run
			fclose(ifp);

// Compared to early versions, Arverts 6 and 7 create a smaller number of randomized time nodes at the Monte Carlo stage
// and then at the CRS stage interpolate regularly spaced, more frequent time nodes.
// For compatibility with old code, keep the variable ttpoints for use in CRS algorithm; introduce rttpoints as number of random MC nodes
		
   		printf("  ...starting to make %ld monte-carlo histories\n",poolsize);

// make random set of Toffsets for CRS use unless user wants to use fixed inputs

		Ranq1 ran(rand()); // reseed RNG
		
		if (offset_flag > 0)
		{
			for (imonte = 1;imonte<= poolsize;imonte++)
			{
				Toffset[imonte][1] = 0.0;
				for (n=2; n<=nsamples; n++)
				{
					if (offset_flag == 1)  // max and min values are for INCREMENT
					{
						Toffset[imonte][n] = Toffset[imonte][n-1] + ran.doub() * (offset_max - offset_min);				
					}
					else // max and min values are for ENTIRE POSSIBLE OFFSET RANGE
					{
						Toffset[imonte][n] = Toffset[imonte][1] + offset_min + ran.doub() * (offset_max - offset_min) ;	// frame as absolute offset from reference				
					}
				
				}
			}
		}
		else   // fix all Toffsets to input values (easiest way to handle later code)
		{
			for (imonte = 1;imonte<= poolsize;imonte++)
			{
				Toffset[imonte][1] = 0.0;
				for (n=2; n<=nsamples; n++)
				{
					Toffset[imonte][n] = Toffset[1][n];
				}
			}		
		}
		
		for (imonte = 1;imonte<= poolsize;imonte++)
		{
			monteTt(rttpoints, ttpoints, monteHeatRate, monteCoolRate, nconstraints, constraintTime, constraintTempmin, constraintTempmax, MonteTime, MonteTemp, temperaturein, timein, flips);
	
		// place new monte carlo history into main CRS pool
			for (ijj = 1;ijj <= ttpoints;ijj++)
			{
				poolTemp[imonte][ijj] = temperaturein[ijj];
				poolTime[imonte][ijj] = timein[ijj];
				mc_Temp[imonte][ijj]  = poolTemp[imonte][ijj]; //save a copy of monte carlo tT pool for possible plotting (will not be sorted)(so what)
			}

			if (use_minerals < 2)
			{			
	// ----- calculate an age spectrum for the new history
				for (n=1; n <= nsamples; n++)
				{
					if (nspectra[n] > 0)
					{
						// first produce offset temperature history
						for (ijj = 1; ijj <= ttpoints; ijj++)
						{
							offset_temperaturein[ijj] = temperaturein[ijj] + Toffset[imonte][n];
						}

						geolovera(n,ttpoints,geometry,bgeom,slab,ndomains[n], Eact[n], D0[n], fraction[n], labsteps[n], dt, ichange, timein, offset_temperaturein, timelab[n], temperaturelab[n], bulkloss39, age);

		// store the results of the age spectrum calculated by lovera()
						for (ijj = 1; ijj <= labsteps[n]; ijj++)
						{
							poolAge[n][imonte][ijj] = age[n][ijj];
	// DEBUG						printf("%ld\t%ld\t%ld\t%lf\t%lf\n",n,imonte,ijj,bulkloss39[n][ijj],age[n][ijj]);
						}
					}
				}
			}
			
// calculate mineral age(s) if we use them
			if (use_minerals > 0)
			{
				for (n=1; n <= nsamples; n++)
				{
					if (nminerals[n] > 0)
					{
						// first produce offset temperature history
						for (ijj = 1; ijj <= ttpoints; ijj++)
						{
							offset_temperaturein[ijj] = temperaturein[ijj] + Toffset[imonte][n];
						}
		
						for (imineral = 1; imineral <= nminerals[n]; imineral++)
						{
							switch (mineral[n][imineral])
							{
								case 0:  // volume diffusion apatite, alpha correction determined by Ft in mineral.in
								case 1:  // volume diffusion zircon, alpha correction determined by Ft in mineral.in
								case 2:  // volume diffusion titanite, alpha correction determined by Ft in mineral.in
								case 3:  // volume diffusion other mineral, no alpha correction
								{
									mineral_pool[n][imonte][imineral] = mineral_age_vd(ndL, ndt, mineral[n][imineral], mineral_radius[n][imineral], mineral_eact[n][imineral], mineral_difzero[n][imineral], modelduration, timein, offset_temperaturein, ttpoints);
									break;
								}
								case 4:   // RDAAM apatite
								{
								// set up RDAAM MODEL for current sample
									RDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius-microns, U, Th, Sm(ppm)	
			
									path.clear();
									mypoint.time = 0.0;   //   BUG??? other later calls don't do this!!!
									mypoint.temperature = offset_temperaturein[ttpoints+1];
									path.push_back(mypoint);

									for (i=ttpoints; i >= 1; i--)
									{
										mypoint.time = timein[i];
										mypoint.temperature = offset_temperaturein[i];
										path.push_back(mypoint);   // Present day is first entry in path
									}
				
									int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);		
									mineral_pool[n][imonte][imineral] = uncorr_age;
									RDAAM_FreeCalcArrays();
									break;
								}
								case 5:     // ZrDAAM zircon
								{
								// set up ZRDAAM MODEL for current sample
									ZRDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius, U, Th, Sm			
		
									path.clear();
									for (i=ttpoints; i >= 1; i--)
									{
										mypoint.time = timein[i];
										mypoint.temperature = offset_temperaturein[i];
										path.push_back(mypoint);   // Present day is first entry in path
									}
									int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);
									mineral_pool[n][imonte][imineral] = uncorr_age;
									RDAAM_FreeCalcArrays();
									break;
								}
							}		
						}
					} // end of if-nminerals[]
				}  // end of n-loop across samples	
			}
		
		// get fit this new Monte tT
			ModelFit[imonte] = fit(nsamples,nspectra,fitchoice, imonte, poolAge, goalsteps, goal_loss, goal_age, goal_error, goal_skip, use_minerals, nminerals, mineral_goal, mineral_goalerr, mineral_weight, mineral_pool,early_loss,early_weight);
		}

// leave the following double-check block alive for now until we're more confident we're making 
// consistently good single-reheating histories for the starting Monte-Carlo pool

		long flip, sign1, sign2;
		for (imonte = 1; imonte <= poolsize; imonte++)
		{
			flip = 0;	
			for (i = 1;i<= ttpoints-2;i++)
			{
				sign1 = ((poolTemp[imonte][i] - poolTemp[imonte][i+1]) <= 0 ) ? 1 : -1;
				sign2 = ((poolTemp[imonte][i+1] - poolTemp[imonte][i+2]) <= 0 ) ? 1 : -1;
				if ((sign1 * sign2) < 0) { flip = flip + 1;}
			}				
//				assert(flip < 4 && "Monte-carlo thermal history has more than two reversals in trend");	// we're just allowing one reheating pulse
			if (flip > (flips + 3))
			{
				printf("\nWARNING: at MC model %ld\tflips are %ld\n",imonte,flip);
//					for (i = 1;i<= ttpoints;i++)
//					{
//						printf("%ld\t%6.1f\t%6.1f\n",i,poolTime[imonte][i],poolTemp[imonte][i]);
//					}	
			}			
		}			
			
// write out Monte time-temperature data to MONTErestart, in tab-delimited format, for possible later plotting use
// (before changing directory!) - if user chooses restart option they might still want to overlay
// final CRS results over MC pool, thus we have to save it somewhere
	
		ofp = fopen("MONTErestart", "w");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fprintf(ofp,"%8.2f\t", timein[ijj]);
			for (imonte = 1; imonte <= poolsize - 1; imonte++)
			{
				fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
			}
			fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
		}
		fclose(ofp);			

// Change directory in preparation for writing intermediate results files 

// PORT-ISSUE - make sure these standard calls work for your system

		j = 0;	
		do
		{
			if (j == 0)
			{
				sprintf(newFolder,"RESULTS-%s", suffix);		
			}
			else     // don't overwrite older folder, so rename with digit
			{
				sprintf(newFolder,"RESULTS-%s %ld", suffix, j);
			}
			j++;
		} while (mkdir(newFolder, 0755 ) == -1); 

	/* change to new folder */

		if( chdir( newFolder) )
		{
				printf("\nCannot change to new folder\n");
				exit(EXIT_FAILURE); 
		}

	// write out monte-carlo time-temperature data, in tab-delimited format	
	// first column is time, second is average temperature, third and fourth are low and high envelope temperatures, 
	// and the n the actual histories start, best one first

		strcpy(filename,"MONTEtT_");

		strcat(filename,suffix);

		ofp = fopen(strcat(filename,extent), "w");
		
		sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
		averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);

		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fprintf(ofp,"%8.2f\t", MonteTime[ijj]);
			fprintf(ofp,"%8.2f\t", avTemp[ijj]);
			fprintf(ofp,"%8.2f\t", lowEnvelope[ijj]);
			fprintf(ofp,"%8.2f\t", highEnvelope[ijj]);
			for (imonte = 1; imonte <= poolsize - 1; imonte++)
			{

				fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
			}
			fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
		}
		fclose(ofp);

	// write out monte-carlo Toffset values, in tab-delimited format, sorted by fit
	// add initial column to make plotting quicker in Excel	
		strcpy(filename,"MONTE-offset-sort_");
		strcat(filename,suffix);
		ofp = fopen(strcat(filename,extent), "w");	
	
		for (n = 1;n<= nsamples;n++)
		{
			fprintf(ofp,"%ld\t", n);				
			for (imonte = 1; imonte <= poolsize - 1; imonte++)
			{

				fprintf(ofp,"%6.2f\t", Toffset[imonte][n]);	
			}
				fprintf(ofp,"%6.2f\n", Toffset[imonte][n]);		
		}
		fclose(ofp);	

		if (use_minerals < 2)
		{
		// write out monte-carlo spectra in age spectrum format, tab-delimited   
		// first column is fractional 39Ar loss, second is goal spectrum age, third is average of model ages, and then
		// the model predictions follow, best one first

			strcpy(filename,"MONTEage_");
			strcat(filename,suffix);
			ofp = fopen(strcat(filename,extent), "w");
		
			for(n=1; n <= nsamples; n++)
			{	
				if (nspectra[n] > 0)
				{
					fprintf(ofp,"%6.2f\t", bulkloss39[n][1] * 0.0);
					fprintf(ofp,"%6.2f\t", goal_age[n][1]);  // goal_age[]
					fprintf(ofp,"%6.2f\t", avAge[n][1]);  // average_age[]
					for (imonte = 1; imonte <= poolsize - 1; imonte++)
					{	
							fprintf(ofp,"%6.2f\t", poolAge[n][imonte][1]);
					}
					fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][1]);	
				
					for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
					{
						fprintf(ofp,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
						fprintf(ofp,"%6.2f\t", goal_age[n][ijj]);
						fprintf(ofp,"%6.2f\t", avAge[n][ijj]);	
						for (imonte = 1; imonte <= poolsize - 1; imonte++)
						{	
							fprintf(ofp,"%6.2f\t", poolAge[n][imonte][ijj]);
						}
						fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][ijj]);
			
						fprintf(ofp,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
						fprintf(ofp,"%6.2f\t", goal_age[n][ijj+1]);
						fprintf(ofp,"%6.2f\t", avAge[n][ijj+1]);
						for (imonte = 1; imonte <= poolsize - 1; imonte++)
						{	
							fprintf(ofp,"%6.2f\t", poolAge[n][imonte][ijj+1]);
						}
						fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][ijj+1]);
					}
					fprintf(ofp,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);
					fprintf(ofp,"%6.2f\t", goal_age[n][labsteps[n]]);
					fprintf(ofp,"%6.2f\t", avAge[n][labsteps[n]]);
					for (imonte = 1; imonte <= poolsize - 1; imonte++)
					{	
						fprintf(ofp,"%6.2f\t", poolAge[n][imonte][labsteps[n]]);
					}
					fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][labsteps[n]]);
				}
			}  // end of n -loop through samples
		
			fclose(ofp);
		}
	// write out monte-carlo mineral ages
		if (use_minerals > 0)
		{
			strcpy(filename,"MONTE_minerals_");
			strcat(filename,suffix);
			ofp = fopen(strcat(filename,extent), "w");
			
			for(n=1; n <= nsamples; n++)
			{				
				if (nminerals[n] > 0)
				{
					for (i = 1; i <= nminerals[n]; i++)
					{
						fprintf(ofp,"%7.2f\t", average_mineral_age[n][i]);
						for (imonte = 1; imonte <= poolsize - 1; imonte++)
						{	
							fprintf(ofp,"%7.2f\t", mineral_pool[n][imonte][i]);
						}
						fprintf(ofp,"%7.2f\n", mineral_pool[n][poolsize][i]);
					}
				}
			}	  // end of n-loop through samples
			fclose(ofp);
		}
	}
	else  // we're restarting from a pre-existing tT pool, so we don't want the monte carlo routine
	{

		ifp = fopen("CRScount.in", "r");
		fscanf(ifp, "%ld", &totCRS);   // number of previous CRS iterations
		fclose(ifp);
		printf("  ...restarting using previous pool\n");
	
	// first get saved tT data into pool array
	
		ifp = fopen("CRSrestart", "r");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fscanf(ifp,"%lf", &timein[ijj]);
			for (imonte = 1; imonte <= poolsize; imonte++)
			{
				fscanf(ifp,"%lf", &poolTemp[imonte][ijj]);
				mc_Temp[imonte][ijj]  = poolTemp[imonte][ijj]; //save a copy of monte carlo tT pool for possible plotting (will not be sorted)(so what)
			}
		}	
		fclose(ifp);
		
		ifp = fopen("MONTErestart", "r");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fscanf(ifp,"%lf", &timein[ijj]);
			for (imonte = 1; imonte <= poolsize; imonte++)
			{
				fscanf(ifp,"%lf", &mc_Temp[imonte][ijj]); // get a copy of monte carlo tT pool for possible plotting (will not be sorted)(so what)
			}
		}	
		fclose(ifp);

		ifp = fopen("TOFFSETrestart", "r");
		for (imonte = 1; imonte <= poolsize; imonte++)
		{
			for (n=1; n <= nsamples; n++)
			{
				fscanf(ifp,"%lf", &Toffset[imonte][n]); // get Toffset pool saved after earlier run
			}
		}
		fclose(ifp);
				
		printf("  ...opened restart files\n");
				
		// now determine ages for restart pool  
	
		for (imonte = 1;imonte<= poolsize;imonte++)
		{
			for (ijj = 1;ijj <= ttpoints;ijj++)
			{
				temperaturein[ijj] = poolTemp[imonte][ijj];
			}
			
			if (use_minerals < 2)
			{
				for (n=1; n <= nsamples; n++)
				{
					if(nspectra[n] > 0)
					{
						// first produce offset temperature history
						for (ijj = 1; ijj <= ttpoints; ijj++)
						{
							offset_temperaturein[ijj] = temperaturein[ijj] + Toffset[imonte][n];
						}

						geolovera(n,ttpoints,geometry,bgeom,slab,ndomains[n], Eact[n], D0[n], fraction[n], labsteps[n], dt, ichange, timein, offset_temperaturein, timelab[n], temperaturelab[n], bulkloss39, age);

		// assign lovera() age spectrum to pool	
						for (ijj = 1; ijj <= labsteps[n]; ijj++)
						{
							poolAge[n][imonte][ijj] = age[n][ijj];
						}
					}
				}
			}

			if (use_minerals > 0)
			{
				for(n=1; n <= nsamples; n++)
				{	
					if(nminerals[n] > 0)
					{
						// first produce offset temperature history
						for (ijj = 1; ijj <= ttpoints; ijj++)
						{
							offset_temperaturein[ijj] = temperaturein[ijj] + Toffset[imonte][n];
						}

						for (imineral = 1; imineral <= nminerals[n]; imineral++)
						{
							switch (mineral[n][imineral])
							{
								case 0:  // volume diffusion apatite, alpha correction determined by Ft in helium.in
								case 1:  // volume diffusion zircon, alpha correction determined by Ft in helium.in
								case 2:  // volume diffusion titanite, alpha correction determined by Ft in helium.in
								case 3:  // volume diffusion other mineral, no alpha correction
								{
									mineral_pool[n][imonte][imineral] = mineral_age_vd(ndL, ndt, mineral[n][imineral], mineral_radius[n][imineral], mineral_eact[n][imineral], mineral_difzero[n][imineral], modelduration, timein, offset_temperaturein, ttpoints);
									break;
								}
								case 4:   // RDAAM apatite
								{
								// set up RDAAM MODEL for current sample
									RDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius-microns, U, Th, Sm(ppm)	
	
									path.clear();
									for (i=ttpoints; i >= 1; i--)
									{
										mypoint.time = timein[i];
										mypoint.temperature = offset_temperaturein[i];
										path.push_back(mypoint);   // Present day is first entry in path
									}
			
									int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);				
									mineral_pool[n][imonte][imineral] = uncorr_age;
									RDAAM_FreeCalcArrays();
									break;
								}
								case 5:     // ZrDAAM zircon
								{
								// set up ZRDAAM MODEL for current sample
									ZRDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius, U, Th, Sm			
	
									path.clear();
									for (i=ttpoints; i >= 1; i--)
									{
										mypoint.time = timein[i];
										mypoint.temperature = offset_temperaturein[i];
										path.push_back(mypoint);   // Present day is first entry in path
									}
									int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);
									mineral_pool[n][imonte][imineral] = uncorr_age;
									RDAAM_FreeCalcArrays();
									break;
								}
							}		
						}
					}
				} // end of n-loop to work through samples
			} 
		} 	
			
  // Change directory in preparation for writing various intermediate-results files

	// PORT-ISSUE - make sure these standard calls work for your system
		j = 0;	
		do
		{
			if (j == 0)
			{
				sprintf(newFolder,"RESULTS-%s", suffix);		
			}
			else     // don't overwrite older folder, so rename with digit
			{
				sprintf(newFolder,"RESULTS-%s %ld", suffix, j);
			}
			j++;
		} while (mkdir(newFolder, 0755 ) == -1); 

		printf("  ...made new results folder: %s\n",newFolder);
		/* change to new folder */

		if( chdir( newFolder) )
		{
				printf("\nCannot change to new folder\n");
				exit(EXIT_FAILURE); 
		}
		
	} // done with either Monte Carlo generation or loading in histories from previous run

// Now load monte-carlo fits into ModelFit[]
	for (imonte = 1; imonte <= poolsize; imonte++)
	{
		ModelFit[imonte] = fit(nsamples,nspectra,fitchoice, imonte, poolAge, goalsteps, goal_loss, goal_age, goal_error, goal_skip, use_minerals, nminerals, mineral_goal, mineral_goalerr, mineral_weight, mineral_pool,early_loss,early_weight);
	}
	
	WorstOne = worstfit(poolsize, ModelFit); // find the worst fit }
	BestOne = bestfit(poolsize, ModelFit);  // find the best fit

// maybe superfluous -- but want to be certain
	sort(nsamples, nspectra,poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
	averages(nsamples, nspectra,ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);

	if (use_minerals < 2)
	{
	// write out goal spectra in age spectrum format, tab-delimited

		ofp = fopen("goalspec.txt", "w");
	
		for (n=1; n <= nsamples; n++)
		{
			if (nspectra[n] > 0)
			{
				fprintf(ofp,"%6.2f\t", goal_loss[n][1] * 0.0);
				fprintf(ofp,"%6.2f\n", goal_age[n][1]);	
	
				for (ijj = 1; ijj <= goalsteps[n] - 1; ijj++)
				{
					fprintf(ofp,"%6.2f\t", goal_loss[n][ijj] * 100.0);
					fprintf(ofp,"%6.2f\n", goal_age[n][ijj]);
	
					fprintf(ofp,"%6.2f\t", goal_loss[n][ijj] * 100.0);
					fprintf(ofp,"%6.2f\n", goal_age[n][ijj+1]);
				}

				fprintf(ofp,"%6.2f\t", goal_loss[n][goalsteps[n]] * 100.0);
				fprintf(ofp,"%6.2f\n", goal_age[n][goalsteps[n]]);	
			}
		}
		fclose(ofp);
	}
	
	printf("  ...ready to work with %ld thermal histories\n",poolsize);
	printf("\n");
	printf("Worst fit is %5.1f (history %ld)\n",ModelFit[WorstOne],WorstOne);
	printf(" Best fit is %5.1f (history %ld)\n",ModelFit[BestOne],BestOne);
	printf("\n");
	
	if (restart == 0)	
	{
		printf("  ...starting to make %ld new CRS histories...\n\n", CRSiterations);
	}
	else
	{
		printf("  ...starting to make %ld additional CRS histories...\n", CRSiterations);	
	}

// ***************************************************************************************
// ******************* Section 7. MAIN CRS LOOP ******************************************
// ***************************************************************************************

	duration = time(NULL); // PORT-ISSUE can be deleted as timing info is not essential

	iCRS = 0;
// do CRS algorithm until we either hit specified number of iterations, or worst fit is better than specified fit
	while ((iCRS < CRSiterations) && (ModelFit[WorstOne] > fitvalue)) 
	{
		fflush(stdout); //sets up printing of iteration monitor
		iCRS++;
		    
		CRSattempts = 0;
		goodHistory = -1;
		while ((CRSattempts < attempts_CRS) && (goodHistory < 0))  // use repeated tries to make legal CRS history
		{
			CRSattempts++;
			CRS_selection = selectsubset(subsetsize, CRSselections, poolsize);
			newhistoryCRS(ttpoints, goodHistory, CRS_selection, subsetsize, CRSselections, amplify, coolRate, heatRate, poolTemp, upperlimit, lowerlimit, timein, temperaturein, flips, Toffset, newToffset, nsamples, offset_flag, offset_max, offset_min);
		}

		if (CRSattempts >= attempts_CRS)
		{
			printf("\n\n******************************************************************\n");
			printf("* Disaster: no good CRS history after %ld tries. Model aborted. *\n",attempts_CRS);
			printf("******************************************************************\n");			
		}
		assert(CRSattempts < attempts_CRS);

		if (use_minerals < 2)
		{
			for (n=1; n <= nsamples; n++)
			{
				if(nspectra[n] > 0)
				{
					// first produce offset temperature history
					for (ijj = 1; ijj <= ttpoints; ijj++)
					{
						offset_temperaturein[ijj] = temperaturein[ijj] + newToffset[n];
					}

					geolovera(n,ttpoints,geometry,bgeom,slab,ndomains[n], Eact[n], D0[n], fraction[n], labsteps[n], dt, ichange, timein, offset_temperaturein, timelab[n], temperaturelab[n], bulkloss39, age);

	// assign lovera() age spectrum to pool	
					for (ijj = 1; ijj <= labsteps[n]; ijj++)   // put new ages into pool array's last member
					{
						poolAge[n][poolsize + 1][ijj] = age[n][ijj];
					}
				}
			}
		}
		
		if (use_minerals > 0)
		{
			for (n=1; n <= nsamples; n++)
			{
				if (nminerals[n] > 0)
				{
					// first produce offset temperature history
					for (ijj = 1; ijj <= ttpoints; ijj++)
					{
						offset_temperaturein[ijj] = temperaturein[ijj] + newToffset[n];
					}
				
					for (imineral = 1; imineral <= nminerals[n]; imineral++)
					{
						switch (mineral[n][imineral])
						{
							case 0:  // volume diffusion apatite, alpha correction determined by Ft in helium.in
							case 1:  // volume diffusion zircon, alpha correction determined by Ft in helium.in
							case 2:  // volume diffusion titanite, alpha correction determined by Ft in helium.in
							case 3:  // volume diffusion other mineral, no alpha correction
							{
		// put new ages into pool array's last member
								mineral_pool[n][poolsize + 1][imineral] = mineral_age_vd(ndL, ndt, mineral[n][imineral], mineral_radius[n][imineral], mineral_eact[n][imineral], mineral_difzero[n][imineral], modelduration, timein, offset_temperaturein, ttpoints);
								break;
							}
							case 4:   // RDAAM apatite
							{
							// set up RDAAM MODEL for current sample
								RDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius-microns, U, Th, Sm(ppm)	
							// need loop here to push t-t history to path (added nodes for before and after, 		
								path.clear();
								for (i=ttpoints; i >= 1; i--)
								{
									mypoint.time = timein[i]; 
									mypoint.temperature = offset_temperaturein[i];
									path.push_back(mypoint);   // Present day is first entry in path
								}
			
								int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);
								mineral_pool[n][poolsize + 1][imineral] = uncorr_age;
								RDAAM_FreeCalcArrays();
								break;
							}
							case 5:     // ZrDAAM zircon
							{
							// set up ZRDAAM MODEL for current sample
								ZRDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius, U, Th, Sm			
							// need loop here to push t-t history to path
								path.clear();
								for (i=ttpoints; i >= 1; i--)
								{
									mypoint.time = timein[i];
									mypoint.temperature = offset_temperaturein[i];
									path.push_back(mypoint);   // Present day is first entry in path
								}
								int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);
								mineral_pool[n][poolsize + 1][imineral] = uncorr_age;
								RDAAM_FreeCalcArrays();
								break;
							}
						}	// end switch case for mineral	
					}  // end loop through this sample's minerals
				}
			} // end n-loop through samples 
		}  // end of use-minerals block

// load new history into pool array's last element
		for  (ijj = 1; ijj <= ttpoints; ijj++) 
		{
			poolTemp[poolsize + 1][ijj] = temperaturein[ijj];
		}

// check fit of new spectrum and any ages	
		ModelFit[poolsize + 1] = fit(nsamples, nspectra, fitchoice, poolsize + 1, poolAge, goalsteps, goal_loss, goal_age, goal_error, goal_skip, use_minerals, nminerals, mineral_goal, mineral_goalerr, mineral_weight, mineral_pool,early_loss,early_weight);  // get fit for new age spectrum and He age

// if a better fit, swap element poolsize+1 with posn WorstOne	KEY STEP IN CRS ALGORITHM!!!!!
		if (ModelFit[poolsize + 1] < ModelFit[WorstOne])
		{
			for (i = 1; i <= ttpoints; i++)    // swap the temperature history into pool, and...
			{
				poolTemp[WorstOne][i] = poolTemp[poolsize + 1][i];				
			}
			for (n=1; n <= nsamples; n++)
			{
				for (i = 1; i <= labsteps[n]; i++)    // swap the age spectrum into pool
				{
					poolAge[n][WorstOne][i] = poolAge[n][poolsize + 1][i];
				}
				if (use_minerals > 0)    // swap mineral age(s) into pool
				{
					for (i = 1; i <= nminerals[n]; i++)
					{
						mineral_pool[n][WorstOne][i] = mineral_pool[n][poolsize + 1][i];
					}
				}
				
				Toffset[WorstOne][n] = newToffset[n];
				
			} // end n-loop through samples
			
			ModelFit[WorstOne] = ModelFit[poolsize + 1];
			WorstOne = worstfit(poolsize, ModelFit); // if we replaced the previous worst fit, we have to find the new worst fit }
			BestOne = bestfit(poolsize, ModelFit);
			
// ~_maybe_ not needed but I am trying to track down weird best/worst swap bug sometimes reported to console		
			sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
			averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);
		}

		if ((iCRS % 50 == 0) && (fullreport == 1))
		{
			printf("\r");  // had to put fflush(stdout) at start of main loop to make this work
			printf("  ...processed %4ld histories: best fit is %5.2f, worst fit is %11.2f",iCRS+totCRS,ModelFit[BestOne],ModelFit[WorstOne]);
		}
		else
		{
			if ((iCRS % 500 == 0) && (fullreport == 0))
			{
				printf("\r");			
				printf(" ...processed %4ld histories: best fit is %5.2f, worst fit is %11.2f",iCRS+totCRS,ModelFit[BestOne],ModelFit[WorstOne]);
			}		
		}

		// write out interim CRS time-temperature data, in tab-delimited format (but only for new runs, not restarts!)
		if (restart == 0)
		{
			if (((iCRS == 100) || (iCRS == 200) || (iCRS == 500) || (iCRS == 1000) || (iCRS == 2000)) &&  (fullreport == 1))
			{
				sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
				averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);
				nowCRS = iCRS + totCRS;
				sprintf(isuffix,"%ld",nowCRS);
				strcpy(filename,"CRStT_");
				strcat(filename,isuffix);
				strcat(filename,extent);
				ofp = fopen(filename, "w");
				for (ijj = 1;ijj<= ttpoints;ijj++)
				{
					fprintf(ofp,"%8.2f\t", timein[ijj]);
					fprintf(ofp,"%8.2f\t", avTemp[ijj]);
					fprintf(ofp,"%8.2f\t", lowEnvelope[ijj]);
					fprintf(ofp,"%8.2f\t", highEnvelope[ijj]);
					for (imonte = 1; imonte <= poolsize - 1; imonte++)
					{
						fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
					}
					fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
				}
				fclose(ofp);
			}
		}
			// write out rolling-progress CRS time-temperature data, in tab-delimited format
		if (iCRS % 500 == 0)
		{
			sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
			averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);
			strcpy(filename,"CRStT_rolling.txt");
			ofp = fopen(filename, "w");
			for (ijj = 1;ijj<= ttpoints;ijj++)
			{
				fprintf(ofp,"%8.2f\t", timein[ijj]);
				fprintf(ofp,"%8.2f\t", avTemp[ijj]);
				fprintf(ofp,"%8.2f\t", lowEnvelope[ijj]);
				fprintf(ofp,"%8.2f\t", highEnvelope[ijj]);
				for (imonte = 1; imonte <= poolsize - 1; imonte++)
				{
					fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
				}
				fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
			}
			fclose(ofp);		
			
			if (use_minerals < 2)
			{
			// write out CRS spectra in age spectrum format, tab-delimited

				strcpy(filename,"CRSage_");
				strcat(filename,isuffix);
				ofp = fopen(strcat(filename,extent), "w");
			
				for (n=1; n <= nsamples; n++)
				{	
					if (nspectra[n] > 0)
					{
						fprintf(ofp,"%6.2f\t", bulkloss39[n][1] * 0.0);
						fprintf(ofp,"%6.2f\t", goal_age[n][1]);
						fprintf(ofp,"%6.2f\t", avAge[n][1]);			
						for (imonte = 1; imonte <= poolsize - 1; imonte++)
						{	
								fprintf(ofp,"%6.2f\t", poolAge[n][imonte][1]);
						}
						fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][1]);	
				
						for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
						{
							fprintf(ofp,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							fprintf(ofp,"%6.2f\t", goal_age[n][ijj]);
							fprintf(ofp,"%6.2f\t", avAge[n][ijj]);
							for (imonte = 1; imonte <= poolsize - 1; imonte++)
							{	
								fprintf(ofp,"%6.2f\t", poolAge[n][imonte][ijj]);
							}
							fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][ijj]);
				
							fprintf(ofp,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							fprintf(ofp,"%6.2f\t", goal_age[n][ijj+1]);
							fprintf(ofp,"%6.2f\t", avAge[n][ijj+1]);
							for (imonte = 1; imonte <= poolsize - 1; imonte++)
							{	
								fprintf(ofp,"%6.2f\t", poolAge[n][imonte][ijj+1]);
							}
							fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][ijj+1]);
						}
						fprintf(ofp,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);
						fprintf(ofp,"%6.2f\t", goal_age[n][labsteps[n]]);
						fprintf(ofp,"%6.2f\t", avAge[n][labsteps[n]]);
						for (imonte = 1; imonte <= poolsize - 1; imonte++)
						{	
							fprintf(ofp,"%6.2f\t", poolAge[n][imonte][labsteps[n]]);
						}
						fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][labsteps[n]]);
					}
				}
				fclose(ofp);
			}
			
			if (use_minerals > 0)
			{			
				strcpy(filename,"CRS_minerals_");
				strcat(filename,isuffix);
				ofp = fopen(strcat(filename,extent), "w");
				for (n=1; n <= nsamples; n++)
				{	
					if(nminerals[n] > 0)
					{
						for (i = 1; i <= nminerals[n]; i++)
						{
							fprintf(ofp,"%7.2f\t", average_mineral_age[n][i]);
							for (imonte = 1; imonte <= poolsize - 1; imonte++)
							{	
								fprintf(ofp,"%7.2f\t", mineral_pool[n][imonte][i]);
							}
							fprintf(ofp,"%7.2f\n", mineral_pool[n][poolsize][i]);
						}
					}
				}
				fclose(ofp);
			}
		} //end of every-500 report block
	}  // end of main CRS loop
	
	// get final time info

		elapsed = (time(NULL) - duration)/60.0; // PORT-ISSUE can be deleted as timing not essential
		runrate = elapsed/iCRS*100; // PORT-ISSUE
		totCRS = totCRS + iCRS; // PORT-ISSUE
	
	// Begin final writing of output files and summaries

		sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
		averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);

// write out run summary

		sprintf(filename,"ModelInfo_%s.txt",suffix);
		ofp = fopen(filename, "w");
		fprintf(ofp,"Arvert 7.0.1 model run %s\n",ctime(&now));
		fprintf(ofp,"                 SAMPLE INFO:   %s\n", info);
		fprintf(ofp,"                 FILE SUFFIX:   %s\n", suffix);
		fprintf(ofp,"\n");
		fprintf(ofp,"           NUMBER OF SAMPLES:   %ld\n", nsamples);	
		for (n=1; n<= nsamples; n++)
		{
			fprintf(ofp,"%2ld              Toffset (˚C):   %5.1f\n",n, Toffset[1][n]);	
		}
		fprintf(ofp,"            TOFFSET MIN, MAX:   %5.1f\t%5.1f\n", offset_min,offset_max);
		switch (offset_flag)
		{
			case 0:
				fprintf(ofp,"              T_OFFSET TYPE?:    0 - fixed to inputs\n");			
				break;
			case 1:
				fprintf(ofp,"              T_OFFSET TYPE?:    1 - CRS monotonic increasing\n");		
				break;
			case 3:
				fprintf(ofp,"              T_OFFSET TYPE?:    2 - CRS random around reference\n");			
		}
		fprintf(ofp,"\n");

		fprintf(ofp,"           MODEL DURATION (m.y.):   %6.1f\n", modelduration);
		fprintf(ofp,"              MC, CRS TIME NODES:    %ld   %ld\n",rttpoints, ttpoints);
		fprintf(ofp,"\n        CONSTRAINING BRACKETS tT:     %ld\n", nconstraints);
		fprintf(ofp,"           TIME    TMIN   TMAX\n");
		for (j = 1; j <= nconstraints; j++)
		{
			fprintf(ofp,"         %6.1f   %5.1f   %5.1f\n", constraintTime[j], constraintTempmin[j], constraintTempmax[j]);
		}
		
		fprintf(ofp,"\n    MAX MONTE-CARLO HEATING RATE: %5.1f   MAX MONTE-CARLO COOLING RATE: %5.1f\n", monteHeatRate,monteCoolRate);
		fprintf(ofp,"            MAX CRS HEATING RATE: %5.1f           MAX CRS COOLING RATE: %5.1f\n", heatRate,coolRate);
		fprintf(ofp,"        CRS AMPLIFICATION FACTOR:   %4.2f\n", amplify);
		fprintf(ofp,"          SUBSET SIZE, POOL SIZE:   %ld   %ld\n\n", subsetsize, poolsize);
		
		fprintf(ofp,"               FITTING CRITERION:  %5.2f\n", fitvalue);
		fprintf(ofp,"                  FITTING OPTION:   %ld  ", fitchoice);
		switch (fitchoice)
		{
			case 0:
				fprintf(ofp,"(early-weighted mean percent deviation)\n");
				break;
			case 1:
				fprintf(ofp,"(mean percent deviation)\n");
				break;
			case 2:
				fprintf(ofp,"(mswd)\n");
		}

		fprintf(ofp,"          MDD DIFFUSION GEOMETRY:   %ld  ", geometry);
		if (geometry == 1)
		{
			fprintf(ofp,"(spherical)\n");
		}
		else
		{
			fprintf(ofp,"(infinite-slab)\n");
		}
		fprintf(ofp,"                MAX CRS ATTEMPTS:  %ld\n", attempts_CRS);

		fprintf(ofp,"DISCRETIZATION DELTA-TEMPERATURE: %5.1f\n", dt);
		fprintf(ofp,"\n           LOVERA SERIES CUT-OFF:   %5.1e\n", ichange);
		fprintf(ofp,"      FLAG TO WRITE FULL REPORTS:   %ld\n", fullreport);

		switch (gmtplot)
		{
			case 1: 
				fprintf(ofp,"           PLOT RESULTS WITH GMT:   %ld (script ran, simple plot)\n", gmtplot);		
			break;
			case 2:
				fprintf(ofp,"           PLOT RESULTS WITH GMT:   %ld (script ran, full plot)\n", gmtplot);		
			break;
			case 0:
				fprintf(ofp,"           PLOT RESULTS WITH GMT:   %ld (no script ran, no plot)\n", gmtplot);		
			break;
		}

		fprintf(ofp,"     # OF BEST HISTORIES TO PLOT:   %ld\n", nbest); 
		fprintf(ofp,"    # OF WORST HISTORIES TO PLOT:   %ld\n", nworst);
		fprintf(ofp,"      MAX PERMITTED tT REVERSALS:   %ld\n\n", flips); 
	
		for(n=1; n <= nsamples; n++)
		{
			fprintf(ofp,"%ld EARLY-STEP, WEIGHT:   %5.3f\t%5.3f\n", n,early_loss[n],early_weight[n]);
		}

		if (use_minerals == 1)
		{
			fprintf(ofp,"\n CAN MINERAL AGES BE CONSTRAINTS:   1 - yes\n");
		}
		if (use_minerals == 0)
		{
			fprintf(ofp,"\n CAN MINERAL AGES BE CONSTRAINTS:   0 - no\n");
		}
		if (use_minerals == 2)
		{
			fprintf(ofp,"\n CAN MINERAL AGES BE CONSTRAINTS:   2 - yes (and age spectra ignored!!)\n");
		}
		
		fprintf(ofp,"\n                  CRS ITERATIONS:  %ld\n", CRSiterations);
		fprintf(ofp,"\n      POST-PROCESSING ITERATIONS:   %ld\n", npostprocess);
		fprintf(ofp,"\n          POST-PROCESSING SUBSET:   %ld\n", post_set);		
		fprintf(ofp,"       POST-PROCESSING RANGE ˚C):   %5.1f\n", prange);
		fprintf(ofp,"                  RESTART OPTION:   %ld  ", restart);
		switch (restart)
		{
			case 0 :
				fprintf(ofp,"(new start from Monte-Carlo histories)\n");
				break;
			case 1 :
				fprintf(ofp,"(restart CRS processing)\n");
				break;
			case 2 :
				fprintf(ofp,"(restart for post-processing)\n");
				break;
		}

// send mineral information 
		fprintf(ofp,"\n----------------------------------- MINERAL-AGE INFO -----------------------------------\n");
		for (n=1; n <= nsamples; n++)
		{
			fprintf(ofp,">>> SAMPLE %ld NUMBER OF MINERAL AGES:   %ld\n",n, nminerals[n]);
			if(nminerals[n] > 0)
			{
				for	(i = 1; i <= nminerals[n]; i++)
				{
					switch (mineral[n][i])
					{
						case 0:
							fprintf(ofp,"                ***** MINERAL %ld:   0 - apatite (volume diffusion)\n",i);
							break;
						case 1:
							fprintf(ofp,"                ***** MINERAL %ld:   1 - zircon (volume diffusion)\n",i);
							break;
						case 2:
							fprintf(ofp,"                ***** MINERAL %ld:   2 - titanite (volume diffusion)\n",i);
							break;
						case 3:
							fprintf(ofp,"                ***** MINERAL %ld:   3 - other (no alpha loss\n",i);
							break;
						case 4:
							fprintf(ofp,"                ***** MINERAL %ld:   4 - apatite (RDAAM)\n",i);
							break;
						case 5:
							fprintf(ofp,"                ***** MINERAL %ld:   5 - zircon (ZRDAAM)\n",i);
							break;			
					}	
					fprintf(ofp,"             GOAL MINERAL AGE %ld:  %6.2f Ma     ERROR IN AGE: %5.2f Ma\n",i, mineral_goal[n][i], mineral_goalerr[n][i]);
					fprintf(ofp," MINERAL-AGE WEIGHTING FACTOR %ld:  %6.2f\n",i, mineral_weight[n][i]);
					fprintf(ofp,"             EFFECTIVE RADIUS %ld: %7.2f microns\n",i, mineral_radius[n][i]);
					fprintf(ofp,"            ACTIVATION ENERGY %ld:   %5.2f kcal/mol     DIFZERO:  %9.3e cm2/sec\n",i, mineral_eact[n][i], mineral_difzero[n][i]);
					fprintf(ofp,"        PARENT CONCENTRATIONS %ld:   %8.2f ppm U  %8.2f ppm Th  %8.2f ppm Sm\n\n",i, mineral_U[n][i], mineral_Th[n][i], mineral_Sm[n][i]);	
				}
			}
		} //end of n-loop through samples
	
		if (diffusion_precision == 0)
		{
			fprintf(ofp," DIFFUSION and RDAAM PRECISION?:    0 - good \n");
		}
		if (diffusion_precision == 1)
		{
			fprintf(ofp," DIFFUSION and RDAAM PRECISION?:    1 - better\n");
		}
		if (diffusion_precision == 2)
		{
			fprintf(ofp," DIFFUSION and RDAAM PRECISION?:    2 - best\n");
		}

		fprintf(ofp,"\n");
		fprintf(ofp,"DOMAIN INFO\n");
		fprintf(ofp,"-----------\n");
		for (n=1; n <= nsamples; n++)
		{
			fprintf(ofp,"Domains: %2ld\n", ndomains[n]);
			if(nspectra[n] > 0)
			{
				fprintf(ofp,"      E      D0/a2     frac.\n");  // BUG BUG BUG -- we were reporting Do/a2 not Do, and in 1/m.y.!!!  
				for (j = 1; j <= ndomains[n]; j++)
				{
					fprintf(ofp,"%2ld  %5.2f  %9.3e  %6.4f\n", j, Eact[n][j], D0[n][j]/3.1557600e13, fraction[n][j]);	
				}
			}
		} // end n-loop through samples
	
	fprintf(ofp,"\n");
	for (n=1; n <= nsamples; n++)
	{
		if(nspectra[n] > 0)
		{
			fprintf(ofp,"%ld  Goal Age Spectrum\n",n);
			fprintf(ofp,"-----------------\n");	
			fprintf(ofp,"Goal spectrum steps: %3ld\n", goalsteps[n]);
			fprintf(ofp,"      f39     age    error  skip?\n");
			for (j = 1; j <= goalsteps[n]; j++)
			{
				fprintf(ofp,"%3ld  %5.3f  %6.1f  %5.1f   %1ld\n", j, goal_loss[n][j], goal_age[n][j], goal_error[n][j], goal_skip[n][j]);
			}
			fprintf(ofp,"\n");
			fprintf(ofp,"%ld  Heating Schedule Actually Used\n",n);
			fprintf(ofp,"------------------------------\n");	
			fprintf(ofp,"Heating Steps: %3ld\n",labsteps[n]);
			fprintf(ofp,"     Temp. (C)     Time\n");
			for (i = 1; i <= labsteps[n]; i++)
			{
				fprintf(ofp,"%3ld    %5.1f       10.0\n", i, temperaturelab[n][i] - 273.15, timelab[n][i]);
			}
		}
	}	
	fclose(ofp);

// write out final CRS time-temperature data, in tab-delimited format	

	strcpy(finalfilename,"CRStT_final_");
	strcat(finalfilename,suffix);
	ofp = fopen(strcat(finalfilename,extent), "w");
	for (ijj = 1;ijj<= ttpoints;ijj++)
	{
		fprintf(ofp,"%8.2f\t", timein[ijj]);
		fprintf(ofp,"%8.2f\t", avTemp[ijj]);
		fprintf(ofp,"%8.2f\t", lowEnvelope[ijj]);
		fprintf(ofp,"%8.2f\t", highEnvelope[ijj]);
		for (imonte = 1; imonte <= poolsize - 1; imonte++)
		{
			fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
		}
		fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
	}
	fclose(ofp);	
	
// write out final CRS spectra in age spectrum format, tab-delimited

	strcpy(finalfilename,"CRSage_final_");
	strcat(finalfilename,suffix);
	ofp = fopen(strcat(finalfilename,extent), "w");

	for (n=1; n <= nsamples; n++)
	{
		if(nspectra[n] > 0)
		{
			fprintf(ofp,"%6.2f\t", bulkloss39[n][1] * 0.0);
			fprintf(ofp,"%6.2f\t", goal_age[n][1]);
			fprintf(ofp,"%6.2f\t", avAge[n][1]);
			for (imonte = 1; imonte <= poolsize - 1; imonte++)
			{	
					fprintf(ofp,"%6.2f\t", poolAge[n][imonte][1]);
			}
			fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][1]);	
	
			for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
			{
				fprintf(ofp,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
				fprintf(ofp,"%6.2f\t", goal_age[n][ijj]);
				fprintf(ofp,"%6.2f\t", avAge[n][ijj]);
				for (imonte = 1; imonte <= poolsize - 1; imonte++)
				{	
					fprintf(ofp,"%6.2f\t", poolAge[n][imonte][ijj]);
				}
				fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][ijj]);
	
				fprintf(ofp,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
				fprintf(ofp,"%6.2f\t", goal_age[n][ijj+1]);
				fprintf(ofp,"%6.2f\t", avAge[n][ijj+1]);
				for (imonte = 1; imonte <= poolsize - 1; imonte++)
				{	
					fprintf(ofp,"%6.2f\t", poolAge[n][imonte][ijj+1]);
				}
				fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][ijj+1]);
			}
			fprintf(ofp,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);
			fprintf(ofp,"%6.2f\t", goal_age[n][labsteps[n]]);
			fprintf(ofp,"%6.2f\t", avAge[n][labsteps[n]]);
			for (imonte = 1; imonte <= poolsize - 1; imonte++)
			{	
				fprintf(ofp,"%6.2f\t", poolAge[n][imonte][labsteps[n]]);
			}
			fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][labsteps[n]]);
		}
	} // end of n-loop through samples
	fclose(ofp);
	
// write out final pool of mineral ages
	if (use_minerals > 0)
	{
		strcpy(filename,"CRS_minerals_final_");
		strcat(filename,suffix);
		ofp = fopen(strcat(filename,extent), "w");
		for (n=1; n <= nsamples; n++)
		{
			if(nminerals[n] > 0)
			{
				for (i = 1; i <= nminerals[n]; i++)		
				{
					fprintf(ofp,"%7.2f\t", average_mineral_age[n][i]);
					for (imonte = 1; imonte <= poolsize - 1; imonte++)
					{	
						fprintf(ofp,"%7.2f\t", mineral_pool[n][imonte][i]);
					}
					fprintf(ofp,"%7.2f\n", mineral_pool[n][poolsize][i]);
				}
			}
		} // end of n-loop through samples
		fclose(ofp);
	}	

// write out Modelfits if user asked for full report

	if (fullreport == 1)
	{
		strcpy(finalfilename,"CRSFitsfinal_");
		strcat(finalfilename,suffix);
		ofp = fopen(strcat(finalfilename,extent), "w");
	
		for (imonte = 1; imonte <= poolsize - 1; imonte++)
		{
			fprintf(ofp,"%7.3f\n", ModelFit[imonte]);	
		}
		fprintf(ofp,"%7.3f", ModelFit[poolsize]);	
		fclose(ofp);	
	}
	
// write out age-spectrum step fits for best-fit model

	double stepweight, deviation;
	ofp = fopen("aa-stepfit.txt", "w");
	if (use_minerals < 2)  // we used the spectrum
	{
		fprintf(ofp,"Individual deviations for best-fit spectrum and mineral(s)\n");
		fprintf(ofp,"Step\tMisfit\tObserved\tModeled\tError\n");	
		for (n=1; n <= nsamples; n++)
		{	
			if(nspectra[n] > 0)
			{
				for (j = 1; j<= goalsteps[n]; j++)
				{
					if (goal_skip[n][j] > 0)
					{
						switch (fitchoice)
						{
							case 0:  // step-weighted mean percent deviation
								if (goal_loss[n][j] < early_loss[n])   // empiricially weight some early steps more given that they have more thermal relief
								{
									stepweight = early_weight[n];
								}
								else     // just use step size as weighting factor
								{
									stepweight = 1.0;
								}

								deviation = stepweight*fabs(poolAge[n][1][j] - goal_age[n][j])/goal_age[n][j];
							break;
							case 1:  // mean percent deviation
								deviation = fabs(poolAge[n][1][j] - goal_age[n][j])/goal_age[n][j];		
							break;
			
							case 2:  // mswd
								deviation = 1/(goal_error[n][j]*goal_error[n][j])*pow((poolAge[n][1][j] - goal_age[n][j]),2.0);
							break;
						}
						fprintf(ofp,"%ld\t%8.2f\t%7.2f\t%7.2f\t%7.2f\n",j,deviation,goal_age[n][j],poolAge[n][1][j],goal_error[n][j]);
					}
				}
			}
		} // end of n-loop through samples
	}
		
	// then report best-fit results for individual mineral ages, if they were used as part of inversion
	if (use_minerals > 0)
	{	
		fprintf(ofp,"Mineral\tMisfit\tObserved\tModeled\tError\n");	
		for (n=1; n <= nsamples; n++)
		{
			if (nminerals[n] > 0)
			{
				for (i = 1; i <= nminerals[n]; i++)
				{
					if (fitchoice < 2)  // 1 = mean percent deviation
					{
						deviation = mineral_weight[n][i] * fabs(mineral_pool[n][1][i] - mineral_goal[n][i]) / mineral_goal[n][i];
					}
					else                // 2 = mswd-like parameter
					{
						deviation = mineral_weight[n][i] * pow((mineral_pool[n][1][i] - mineral_goal[n][i]),2.0) / pow(mineral_goalerr[n][i],2.0);
					}
					fprintf(ofp,"%ld\t%8.2f\t%7.2f\t%7.2f\t%7.2f\n",i,deviation,mineral_goal[n][i],mineral_pool[n][1][i],mineral_goalerr[n][i]);
				}
			}
		} // end of n-loop through samples
	}
		
	fclose(ofp);
	
	// write out final CRS Toffset[] data, in tab-delimited format

	strcpy(filename,"OFFSET-CRS_FINAL_");

	strcat(filename,suffix);

	ofp = fopen(strcat(filename,extent), "w");	
	
	for (n = 1;n<= nsamples;n++)
	{
		fprintf(ofp,"%ld\t", n);
		for (imonte = 1; imonte < poolsize; imonte++)
		{
			fprintf(ofp,"%6.2f\t", Toffset[imonte][n]);	
		}
			fprintf(ofp,"%6.2f\n", Toffset[poolsize][n]);		
	}
	fclose(ofp);	

	printf("\n");
	printf("**************** CRS phase finished with no worries! *****************\n");
	printf("\n");
	switch (restart)
	{
		case 0 :
			printf("Processed %ld CRS histories \n",iCRS);
			printf("in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed,runrate);//  // PORT-ISSUE	
		break;
		case 1 :
			printf("Processed %ld additional CRS histories \n",iCRS);	
			printf("in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed,runrate);//  // PORT-ISSUE		
		break;
		case 2 :
			printf("Processed %ld CRS histories\n",iCRS);
			printf("in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed,runrate);//  // PORT-ISSUE	
			printf("  ...now post-processing an %ld additional histories using range of %5.1f ˚C\n\n",iCRS, npostprocess,prange);		
		break;
	}
	
// save copy of best and worst fits for final reporting 
	bestCRS = ModelFit[1];
	worstCRS = ModelFit[poolsize];
	
// ############################# Start post-processing here ########################
	if (restart > 1)  // post process tT (2)
	{	
		iCRS = 0;
		long pool_selection,sub_selection;
		long pnode;
		double pdTemp1, pdTemp2;
		
		duration = time(NULL); // PORT-ISSUE can be deleted as timing info is not essential		
		
	// do postprocess algorithm until we either hit specified number of iterations, or worst fit is better than specified fit
		while ( (iCRS < npostprocess) ) 
		{
			fflush(stdout); //sets up printing of iteration monitor
			iCRS++;
			
			CRSattempts = 0;
			goodHistory = -1;

	// grab a member of tT pool
			Ranq1 ran(rand());
			// our approach is to use the user-specified post_set as the subset, and do the tT post processing on members of that
			// Not that if this is too high, it will take very many iterations to decently exercise all the subset members
			// We are also depending on various arrays being properly sorted - as conceptualized postprocessing works
			// on the post_set best members of the pool. I guess if this were not true, we'd still be post-processing,
			// but rather chaotically!
			
			sub_selection = 1 + floor(ran.doub() * (post_set)); 
			if (sub_selection  > post_set) sub_selection = post_set; 
						
//***************  first, use CRS-style handling for Toffset

			double offset_centroid[MAXSAMPLES];

			if(offset_flag > 0)
			{
				for (n = 1; n <= nsamples; n++)  // zero centroid[12] array
				{
					offset_centroid[n] = 0.0;
				}

				for (n = 2; n <= nsamples; n++) // loop through the selected histories to get sum for each sample
				{
					for (j = 1; j <= post_set; j++)
					{
						offset_centroid[n] = offset_centroid[n] + Toffset[j][n];
					}
				}
				for  (n = 2; n <= nsamples; n++)   // get mean for each sample
				{
					offset_centroid[n] = offset_centroid[n] / post_set;
				}

				newToffset[1] = 0.0;

				for  (n = 2; n <= nsamples; n++)   // make new Toffset; hard-code amplify as 1.05 to explore more gently
				{
					newToffset[n] = offset_centroid[n] + 1.05 * (offset_centroid[n] - Toffset[sub_selection][n]);
				}
				// let new CRS Toffsets violate constraints, because how to repair them is complex and old code was "not fully correct"
			}
			else
			{
				for  (n = 1; n <= nsamples; n++)   // make new Toffset; added in amplify as test; might be too volatile
				{
					newToffset[n] = Toffset[1][n];
				}			
			}


// *********************** end Toffset postprocessing setup
			

		// now make deflected new trial history (maybe put this into a function, one day)

			// first choose an internal time node at random (or limit this to middle of model, or constrained part?)
	
			pnode = 2 + floor(ran.doub() * (ttpoints-3));
			if (pnode >= ttpoints) pnode = ttpoints - 1; // in the very unlikely case that random exactly = 1	
			if (pnode <= 2) pnode = 2; // in the unexpected case that we return 1

			// we will ignore implicit and flip constraints for postprocessing, but honor explicit limits
			// next, get random temperature within prange for start and end of mode
			temperaturein[1] = poolTemp[sub_selection][1] + prange/2. - ran.doub() * prange;

			temperaturein[ttpoints] = poolTemp[sub_selection][ttpoints] + prange/2. - ran.doub() * prange;

	
			temperaturein[pnode] = poolTemp[sub_selection][pnode];

			//next, reshape history using the two temperature ranges and the selected pivot node

			pdTemp1 = temperaturein[1] - poolTemp[sub_selection][1];
			for (ijj = 2; ijj < pnode; ijj++)  // working up to pivot
			{
				temperaturein[ijj] = poolTemp[sub_selection][ijj] + static_cast<double>(pnode-ijj)/static_cast<double>(pnode-1) * pdTemp1;
				if (temperaturein[ijj] > upperlimit[ijj]) temperaturein[ijj] = upperlimit[ijj];
				if (temperaturein[ijj] < lowerlimit[ijj]) temperaturein[ijj] = lowerlimit[ijj];	
			}
			pdTemp2 = temperaturein[ttpoints] - poolTemp[sub_selection][ttpoints];
			for (ijj = pnode+1; ijj < ttpoints; ijj++)  // working down to pivot
			{
				temperaturein[ijj] = poolTemp[sub_selection][ijj] + static_cast<double>(ijj-pnode)/static_cast<double>(ttpoints-pnode) * pdTemp2;	
				if (temperaturein[ijj] > upperlimit[ijj]) temperaturein[ijj] = upperlimit[ijj];
				if (temperaturein[ijj] < lowerlimit[ijj]) temperaturein[ijj] = lowerlimit[ijj];	
			}
	
			if (temperaturein[1] > upperlimit[1]) temperaturein[1] = upperlimit[1];
			if (temperaturein[1] < lowerlimit[1]) temperaturein[1] = lowerlimit[1];	
			if (temperaturein[ttpoints] > upperlimit[ttpoints]) temperaturein[ttpoints] = upperlimit[ttpoints];
			if (temperaturein[ttpoints] < lowerlimit[ttpoints]) temperaturein[ttpoints] = lowerlimit[ttpoints];	


			if (use_minerals < 2)
			{
				for (n=1; n <= nsamples; n++)
				{
					if(nspectra[n] > 0)
					{
						// first produce offset temperature history
						for (ijj = 1; ijj <= ttpoints; ijj++)
						{
							offset_temperaturein[ijj] = temperaturein[ijj] + newToffset[n]; 
						}

						geolovera(n,ttpoints,geometry,bgeom,slab,ndomains[n], Eact[n], D0[n], fraction[n], labsteps[n], dt, ichange, timein, offset_temperaturein, timelab[n], temperaturelab[n], bulkloss39, age);

		// assign lovera() age spectrum to pool	
						for (ijj = 1; ijj <= labsteps[n]; ijj++)   // put new ages into pool array's last member
						{
							poolAge[n][poolsize + 1][ijj] = age[n][ijj];
						}
					}
				}
			}
		
			if (use_minerals > 0)
			{
				for (n=1; n <= nsamples; n++)
				{
					if (nminerals[n] > 0)
					{
						// first produce offset temperature history
						for (ijj = 1; ijj <= ttpoints; ijj++)
						{
							offset_temperaturein[ijj] = temperaturein[ijj] + newToffset[n];  
						}
				
						for (imineral = 1; imineral <= nminerals[n]; imineral++)
						{
							switch (mineral[n][imineral])
							{
								case 0:  // volume diffusion apatite, alpha correction determined by Ft in helium.in
								case 1:  // volume diffusion zircon, alpha correction determined by Ft in helium.in
								case 2:  // volume diffusion titanite, alpha correction determined by Ft in helium.in
								case 3:  // volume diffusion other mineral, no alpha correction
								{
			// put new ages into pool array's last member
									mineral_pool[n][poolsize + 1][imineral] = mineral_age_vd(ndL, ndt, mineral[n][imineral], mineral_radius[n][imineral], mineral_eact[n][imineral], mineral_difzero[n][imineral], modelduration, timein, offset_temperaturein, ttpoints);
									break;
								}
								case 4:   // RDAAM apatite
								{
								// set up RDAAM MODEL for current sample
									RDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius-microns, U, Th, Sm(ppm)	
								// need loop here to push t-t history to path (added nodes for before and after, 		
									path.clear();
									for (i=ttpoints; i >= 1; i--)
									{
										mypoint.time = timein[i]; 
										mypoint.temperature = offset_temperaturein[i];
										path.push_back(mypoint);   // Present day is first entry in path
									}
			
									int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);
									mineral_pool[n][poolsize + 1][imineral] = uncorr_age;
									RDAAM_FreeCalcArrays();
									break;
								}
								case 5:     // ZrDAAM zircon
								{
								// set up ZRDAAM MODEL for current sample
									ZRDAAM_Init(diffusion_precision, mineral_radius[n][imineral], mineral_U[n][imineral], mineral_Th[n][imineral], mineral_Sm[n][imineral]); // precision, radius, U, Th, Sm			
								// need loop here to push t-t history to path
									path.clear();
									for (i=ttpoints; i >= 1; i--)
									{
										mypoint.time = timein[i];
										mypoint.temperature = offset_temperaturein[i];
										path.push_back(mypoint);   // Present day is first entry in path
									}
									int success = RDAAM_Calculate(&path, uncorr_age, hedummy, Heage, optimize);
									mineral_pool[n][poolsize + 1][imineral] = uncorr_age;
									RDAAM_FreeCalcArrays();
									break;
								}
							}	// end switch case for mineral	
						}  // end loop through this sample's minerals
					}
				} // end n-loop through samples 
			}  // end of use-minerals block

	// load new history into pool array's last element
			for  (ijj = 1; ijj <= ttpoints; ijj++) 
			{
				poolTemp[poolsize + 1][ijj] = temperaturein[ijj];
			}

	// check fit of new spectrum and any ages	
			ModelFit[poolsize + 1] = fit(nsamples, nspectra, fitchoice, poolsize + 1, poolAge, goalsteps, goal_loss, goal_age, goal_error, goal_skip, use_minerals, nminerals, mineral_goal, mineral_goalerr, mineral_weight, mineral_pool,early_loss,early_weight);  // get fit for new age spectrum and He age

	// if a better fit, swap element poolsize+1 with posn WorstOne	KEY STEP IN  ALGORITHM!!!!!
			if (ModelFit[poolsize + 1] <= ModelFit[WorstOne])  // made this <= instead of <, since we are just trying to exercise pool
			{
				for (i = 1; i <= ttpoints; i++)    // swap the temperature history into pool, and...
				{
					poolTemp[WorstOne][i] = poolTemp[poolsize + 1][i];				
				}
				for (n=1; n <= nsamples; n++)
				{
					for (i = 1; i <= labsteps[n]; i++)    // swap the age spectrum into pool
					{
						poolAge[n][WorstOne][i] = poolAge[n][poolsize + 1][i];
					}
					if (use_minerals > 0)    // swap mineral age(s) into pool
					{
						for (i = 1; i <= nminerals[n]; i++)
						{
							mineral_pool[n][WorstOne][i] = mineral_pool[n][poolsize + 1][i];
						}
					}
					Toffset[WorstOne][n] = newToffset[n];  
				} // end n-loop through samples
			
				ModelFit[WorstOne] = ModelFit[poolsize + 1];
				WorstOne = worstfit(poolsize, ModelFit); // if we replaced the previous worst fit, we have to find the new worst fit }
				BestOne = bestfit(poolsize, ModelFit);
			
	// ~_maybe_ not needed but I am still trying to track down weird best/worst swap bug sometimes reported to console		
				sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
				averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);
			}

			if ((iCRS % 50 == 0) && (fullreport == 1))
			{
				printf("\r");  // had to put fflush(stdout) at start of main loop to make this work
				printf("  ...post-processed %4ld histories: best fit is %5.2f, worst fit is %11.2f",iCRS,ModelFit[BestOne],ModelFit[WorstOne]);
			}
			else
			{
				if ((iCRS % 500 == 0) && (fullreport == 0))
				{
					printf("\r");			
					printf("  ...post-processed %4ld histories: best fit is %5.2f, worst fit is %11.2f",iCRS,ModelFit[BestOne],ModelFit[WorstOne]);
				}		
			}

				// write out rolling-progress CRS time-temperature data, in tab-delimited format
			if (iCRS % 500 == 0)
			{
				sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool, Toffset);
				averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);
				strcpy(filename,"CRStT_rolling_post.txt");
				ofp = fopen(filename, "w");
				for (ijj = 1;ijj<= ttpoints;ijj++)
				{
					fprintf(ofp,"%8.2f\t", timein[ijj]);
					fprintf(ofp,"%8.2f\t", avTemp[ijj]);
					fprintf(ofp,"%8.2f\t", lowEnvelope[ijj]);
					fprintf(ofp,"%8.2f\t", highEnvelope[ijj]);
					for (imonte = 1; imonte <= poolsize - 1; imonte++)
					{
						fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
					}
					fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
				}
				fclose(ofp);		
			} //end of every-500 report block
		}  // end of postprocessing  loop	
		
		elapsed2 = (time(NULL) - duration)/60.0; // PORT-ISSUE can be deleted as timing not essential
		runrate2 = elapsed2/iCRS*100; // PORT-ISSUE		
		printf("\n\nPost-processed %ld CRS histories\n",iCRS);
		printf("in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed2,runrate2);//  // PORT-ISSUE	

// write out Modelfits if user asked for full report

	if (fullreport == 1)
	{
		strcpy(finalfilename,"PostFitsfinal_");
		strcat(finalfilename,suffix);
		ofp = fopen(strcat(finalfilename,extent), "w");
	
		for (imonte = 1; imonte <= poolsize - 1; imonte++)
		{
			fprintf(ofp,"%7.3f\n", ModelFit[imonte]);	
		}
		fprintf(ofp,"%7.3f", ModelFit[poolsize]);	
		fclose(ofp);		
	}
	
// write out age-spectrum step fits for best-fit model

	double stepweight, deviation;
	ofp = fopen("aa-stepfit_post.txt", "w");	
	
	if (use_minerals < 2)  // we used the spectrum
	{
		fprintf(ofp,"Individual deviations for best-fit spectrum and mineral(s)\n");
		fprintf(ofp,"Step\tMisfit\tObserved\tModeled\tError\n");	
		for (n=1; n <= nsamples; n++)
		{	
			if(nspectra[n] > 0)
			{
				for (j = 1; j<= goalsteps[n]; j++)
				{
					if (goal_skip[n][j] > 0)
					{
						switch (fitchoice)
						{
							case 0:  // step-weighted mean percent deviation
								if (goal_loss[n][j] < early_loss[n])   // empiricially weight some early steps more given that they have more thermal relief
								{
									stepweight = early_weight[n];
								}
								else     // just use step size as weighting factor
								{
									stepweight = 1.0;
								}

								deviation = stepweight*fabs(poolAge[n][1][j] - goal_age[n][j])/goal_age[n][j];
							break;
							case 1:  // mean percent deviation
								deviation = fabs(poolAge[n][1][j] - goal_age[n][j])/goal_age[n][j];		
							break;
			
							case 2:  // mswd
								deviation = 1/(goal_error[n][j]*goal_error[n][j])*pow((poolAge[n][1][j] - goal_age[n][j]),2.0);
							break;
						}
						fprintf(ofp,"%ld\t%8.2f\t%7.2f\t%7.2f\t%7.2f\n",j,deviation,goal_age[n][j],poolAge[n][1][j],goal_error[n][j]);
					}
				}
			}
		} // end of n-loop through samples
	}
		
	// then report best-fit results for individual mineral ages, if they were used as part of inversion
	if (use_minerals > 0)
	{	
		fprintf(ofp,"Mineral\tMisfit\tObserved\tModeled\tError\n");	
		for (n=1; n <= nsamples; n++)
		{
			if (nminerals[n] > 0)
			{
				for (i = 1; i <= nminerals[n]; i++)
				{
					if (fitchoice < 2)  // 1 = mean percent deviation
					{
						deviation = mineral_weight[n][i] * fabs(mineral_pool[n][1][i] - mineral_goal[n][i]) / mineral_goal[n][i];
					}
					else                // 2 = mswd-like parameter
					{
						deviation = mineral_weight[n][i] * pow((mineral_pool[n][1][i] - mineral_goal[n][i]),2.0) / pow(mineral_goalerr[n][i],2.0);
					}
					fprintf(ofp,"%ld\t%8.2f\t%7.2f\t%7.2f\t%7.2f\n",i,deviation,mineral_goal[n][i],mineral_pool[n][1][i],mineral_goalerr[n][i]);
				}
			}
		} // end of n-loop through samples
	}
		
	fclose(ofp);

	printf("\n");
	printf("**************** Postprocessing phase finished with no worries! *****************\n");
	printf("\n");

	} // end of if restart= postprocessing

	printf("**************** Arvert 7.0.0 finished with no worries! *****************\n");
	printf("\n");

// ############################# post-processing ends here ########################

// wrap up reporting about the model run, now that any possible post-processing is done

	sprintf(filename,"ModelInfo_%s.txt",suffix);  // add final info to run-summary file	
	ofp = fopen(filename, "a");
	fprintf(ofp,"\n");
	fprintf(ofp,"********* Arvert 7.0.0 finished with no worries! **********\n");
	fprintf(ofp,"\n");
	switch (restart)
	{
		case 0 :
			fprintf(ofp,"Processed %ld CRS histories \n",iCRS);	
			fprintf(ofp,"in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed,runrate);//  // PORT-ISSUE	
			fprintf(ofp,"\n Best post-CRS fit is %6.2f, worst fit is %6.2f\n",bestCRS,worstCRS);	
		break;
		case 1 :
			fprintf(ofp,"Processed %ld additional CRS histories for a total of %ld\n",iCRS,totCRS);		
			fprintf(ofp,"in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed,runrate);//  // PORT-ISSUE	
			fprintf(ofp,"\n Best post-CRS fit is %6.2f, worst fit is %6.2f\n",bestCRS,worstCRS);	
		break;
		case 2 :
			fprintf(ofp,"Processed %ld CRS histories\n",iCRS);
			fprintf(ofp,"   in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed,runrate);//  // PORT-ISSUE	
			fprintf(ofp,"\nThen post-processed an %ld additional histories using range of %5.1f ˚C\n",npostprocess,prange);
			fprintf(ofp,"   in %5.2f minutes at rate of %5.3f minutes per 100 histories\n",elapsed2,runrate2);//  // PORT-ISSUE	
			fprintf(ofp,"\n        Best post-CRS fit was %6.2f, worst fit was %6.2f\n",bestCRS,worstCRS);
			fprintf(ofp,"\n Best post-processing fit was %6.2f, worst fit was %6.2f\n",ModelFit[1],ModelFit[poolsize]);
		break;
	}

	fclose(ofp);

	// write out final post-processing Toffset[] data, in tab-delimited format

	strcpy(filename,"OFFSET-POST_FINAL_");

	strcat(filename,suffix);

	ofp = fopen(strcat(filename,extent), "w");	
	
	for (n = 1;n<= nsamples;n++)
	{
		fprintf(ofp,"%ld\t", n);
		for (imonte = 1; imonte < poolsize; imonte++)
		{
			fprintf(ofp,"%6.2f\t", Toffset[imonte][n]);	
		}
			fprintf(ofp,"%6.2f\n", Toffset[poolsize][n]);		
	}
	fclose(ofp);	

// Finally, write out utility files for Arvert, into input directory

// Change directory back to input directory

	if( chdir( inFolder) )
	{
			printf("\nCannot change to new folder\n");
			exit(EXIT_FAILURE); 
	}

// set or update iteration counter

	ifp = fopen("CRScount.in", "w");
	fprintf(ifp, "%ld", totCRS);   // update counter file for new total iterations
	fclose(ifp);

// write out final CRS time-temperature data to CRSrestart, in tab-delimited format,
	
	ofp = fopen("CRSrestart", "w");
	for (ijj = 1;ijj<= ttpoints;ijj++)
	{
		fprintf(ofp,"%8.2f\t", timein[ijj]);
		for (imonte = 1; imonte <= poolsize - 1; imonte++)
		{
			fprintf(ofp,"%6.2f\t", poolTemp[imonte][ijj]);	
		}
		fprintf(ofp,"%6.2f\n", poolTemp[poolsize][ijj]);		
	}
	fclose(ofp);

// write out Monte time-temperature data to MONTErestart, in tab-delimited format, for possible later plotting use
	
	ofp = fopen("MONTErestart", "w");
	for (ijj = 1;ijj<= ttpoints;ijj++)
	{
		fprintf(ofp,"%8.2f\t", timein[ijj]);
		for (imonte = 1; imonte <= poolsize - 1; imonte++)
		{
			fprintf(ofp,"%6.2f\t", mc_Temp[imonte][ijj]);	
		}
		fprintf(ofp,"%6.2f\n", mc_Temp[poolsize][ijj]);		
	}
	fclose(ofp);
	
	ofp = fopen("TOFFSETrestart", "w");
	for (imonte = 1; imonte <= poolsize; imonte++)
	{
		for (n=1; n < nsamples; n++)
		{
			fprintf(ofp,"%5.1f\t", Toffset[imonte][n]); // get Toffset pool saved after earlier run
		}
		fprintf(ofp,"%5.1f\n", Toffset[imonte][nsamples]);
	}
	fclose(ofp);

// ***************************************************
// ****************** PLOTTING CODE ******************
// ***************************************************

// when I decided to completely generalize the possible mix of input data, this means that we wan't assume that the 
// reference sample will have any particular data type (or any at all!)
// So... this means that output to the plot script and the script itself need to change.

	if (gmtplot > 0)  // Arvert will write files that can be plotted by plot_arvert.sh, using gmt
	{
		string target, gmtcall;
		long goal_lines, goallines[MAXSAMPLES];
		double minspan;  // plot spacing for plotting mineral age

	// go back to output directory for this run
		chdir( newFolder); 

	// make the plotting directory
		sprintf(newFolder,"PLOTTING");		
		mkdir(newFolder, 0755 ); 

	// change to this plotting directory

		if( chdir( newFolder) )
		{
				printf("\nCannot change to new folder\n");
				exit(EXIT_FAILURE); 
		}

	// copy plotting script from run directory to plotting directory (easier for it to work locally)
		target = "cp ../../plot_arvert.sh .";	
		system(target.c_str());

		fflush(stdout);
		printf("Writing plot data... ");
		fflush(stdout);

// make sure arrays are sorted, one more time
		sort(nsamples, nspectra, poolsize, labsteps, ttpoints, poolAge, poolTemp, poolTime, ModelFit, use_minerals, nminerals, mineral_pool,Toffset);
		averages(nsamples, nspectra, ttpoints, poolsize, labsteps, poolAge, poolTemp, use_minerals, nminerals, mineral_pool, average_mineral_age, avTemp, avAge, lowEnvelope, highEnvelope);

// needed to automate the plot scaling
		double maxtemp = -100.0;  // needed for tT plotting
		for (ijj = 1;ijj<= nconstraints;ijj++)
		{
			if (constraintTempmax[ijj] > maxtemp)
			{
				maxtemp = constraintTempmax[ijj];
			}
		}
		double agelimit;  // needed for age plotting
		agelimit = 1.25*maxgoodgoalage;

		if (agelimit <= 50)
		{
			agelimit = ceil(agelimit/5.0) * 5.0;
		}
		else
		{
			if (agelimit <= 200)
			{
				agelimit = ceil(agelimit/10.0) * 10.0;			
			}
			else
			{
				agelimit = ceil(agelimit/25.0) * 25.0;		
			}
		}
	
	// write needed model limits and other info
		ofpscratch = fopen("limits.in", "w");			
		fprintf(ofpscratch,"%ld\t", static_cast<long>(modelduration));
		fprintf(ofpscratch,"%ld\t", static_cast<long>(maxtemp));
		fprintf(ofpscratch,"%ld\t", static_cast<long>(agelimit));
		fprintf(ofpscratch,"%ld\t", nsamples-1);
		fprintf(ofpscratch,"%ld\t", use_minerals);
		fprintf(ofpscratch,"%ld\t", nworst);
		fprintf(ofpscratch,"%ld\t", nbest);
		fprintf(ofpscratch,"%ld\t", poolsize);
		fprintf(ofpscratch,"%ld\t", offset_flag);		
		fprintf(ofpscratch,"%ld\t", gmtplot);
		fprintf(ofpscratch,"\n");	
		for (n=1; n<=nsamples; n++)
		{
			fprintf(ofpscratch,"%ld\t", nspectra[n]);			
		}	
		fclose(ofpscratch);
		
	// write backdrop for whole tT plot
		ofpscratch = fopen("backdrop.in", "w");		
		fprintf(ofpscratch,"%6.2f\t%6.2f\n", 0.0,0.0);		
		fprintf(ofpscratch,"%6.2f\t%6.2f\n", modelduration,0.0);
		fprintf(ofpscratch,"%6.2f\t%6.2f\n", modelduration,maxtemp);		
		fprintf(ofpscratch,"%6.2f\t%6.2f\n", 0.0,maxtemp);			
		fclose(ofpscratch);	
		
	// write mineral-constrained range if minerals were used
		if (use_minerals > 0) 
		{
			ofpscratch = fopen("mconstraint.in", "w");	
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", min_mineralage,0.0);		
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", max_mineralage,0.0);
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", max_mineralage,maxtemp);		
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", min_mineralage,maxtemp);			
			fclose(ofpscratch);	
		}
		
	// write spectrum-constrained range if spectra were used
		if (use_minerals < 2) 
		{		
			ofpscratch = fopen("sconstraint.in", "w");		
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", mingoodgoalage,0.0);		
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", maxgoodgoalage,0.0);
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", maxgoodgoalage,maxtemp);		
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", mingoodgoalage,maxtemp);			
			fclose(ofpscratch);
		}
		
	// write MC starting tT pool
		ofpscratch = fopen("montetT.in", "w");
		
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fprintf(ofpscratch,"%8.2f\t", MonteTime[ijj]);
			for (imonte = 1; imonte <= poolsize - 1; imonte++)
			{
				fprintf(ofpscratch,"%6.2f\t", mc_Temp[imonte][ijj]);	
			}
			fprintf(ofpscratch,"%6.2f\n", mc_Temp[poolsize][ijj]);		
		}
		fclose(ofp);
		
	// write general pool without best and worst		
		ofpscratch = fopen("pool_general.in", "w");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{		
			fprintf(ofpscratch,"%8.2f\t", timein[ijj]);
			for (imonte = nbest+1; imonte <= poolsize - nworst; imonte++)
			{	
				fprintf(ofpscratch,"%6.2f\t", poolTemp[imonte][ijj]);
			}
			fprintf(ofpscratch,"%6.2f\n", poolTemp[poolsize - nworst][ijj]);
		}
		fclose(ofpscratch);
		
	// write worst fits
		ofpscratch = fopen("pool_worst.in", "w");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{		
			fprintf(ofpscratch,"%8.2f\t", timein[ijj]);
			for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1 ; imonte++)
			{	
				fprintf(ofpscratch,"%6.2f\t", poolTemp[imonte][ijj]);	
			}
			fprintf(ofpscratch,"%6.2f\n", poolTemp[poolsize][ijj]);
		}
		fclose(ofpscratch);

	// write best fits
		ofpscratch = fopen("pool_best.in", "w");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{		
			fprintf(ofpscratch,"%8.2f\t", timein[ijj]);
			for (imonte = 1; imonte <= nbest - 1; imonte++)
			{	
				fprintf(ofpscratch,"%6.2f\t", poolTemp[imonte][ijj]);	
			}
			fprintf(ofpscratch,"%6.2f\n", poolTemp[nbest][ijj]);
		}
		fclose(ofpscratch);
		
	// write envelopes around whole pool
		ofpscratch = fopen("envelope_low.in", "w");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", timein[ijj],lowEnvelope[ijj]);	
		}
		fclose(ofpscratch);
		
		ofpscratch = fopen("envelope_high.in", "w");
		for (ijj = 1;ijj<= ttpoints;ijj++)
		{
			fprintf(ofpscratch,"%6.2f\t%6.2f\n", timein[ijj],highEnvelope[ijj]);	
		}
		fclose(ofpscratch);	
		
// Ship out the Toffset[] data
		
		double olimitl, olimitu;
		olimitu = 0.0;
		for (imonte = 1; imonte <= poolsize; imonte++)  // load 1D array with last column; should be highest
		{
			if (Toffset[imonte][nsamples] > olimitu) olimitu = Toffset[imonte][nsamples];
		}

		// figure out the plotting limits
		if (offset_flag < 2) // constraints are per sample, always positive
		{
			olimitl = 0.0;
			if (olimitu <= 50)
			{
				olimitu = ceil(olimitu/5.0) * 5.0;
			}
			else
			{
				if (olimitu <= 100)
				{
					olimitu = ceil(olimitu/10.0) * 10.0;			
				}
				else
				{
					olimitu = ceil(olimitu/25.0) * 25.0;		
				}
			}
		}
		else  // constraints are for range as specified, could be ±
		{
			if (std::abs(offset_max) <= 50)
			{
				olimitu = std::abs(offset_max)/ offset_max * ceil(std::abs(offset_max)/5.0) * 5.0;
			}
			else
			{
				if (offset_max <= 100)
				{
					olimitu = std::abs(offset_max)/ offset_max * ceil(std::abs(offset_max)/10.0) * 10.0;
				}
				else
				{
					olimitu = std::abs(offset_max)/ offset_max * ceil(std::abs(offset_max)/25.0) * 25.0;	
				}
			}		
			if (std::abs(offset_min) <= 50)
			{
				olimitl = std::abs(offset_min)/ offset_min * ceil(std::abs(offset_min)/5.0) * 5.0;
			}
			else
			{
				if (offset_max <= 100)
				{
					olimitl = std::abs(offset_min)/ offset_min * ceil(std::abs(offset_min)/10.0) * 10.0;	
				}
				else
				{
					olimitl = std::abs(offset_min)/ offset_min * ceil(std::abs(offset_min)/25.0) * 25.0;	
				}
			}	
		}
// need workaround if offsets are all zero - minor issue. And just a warning if something else went wrong.		
		if (olimitl >= olimitu)
		{
			printf(" >>>>> WARNING: olimutl %5.1f > olimutu %5.1f\n",olimitl,olimitu);
			printf(" >>>>> ...setting some blanket values to let code run to completion.\n");
			printf(" >>>>> Be warned that model results might be suspect.");
			if (offset_flag < 2)
			{
				olimitl = 0.0;
				olimitu = 100.0;
			}
			else
			{
				olimitl = -100.0;
				olimitu = 100.0;			
			}
		}
//		assert (olimitl < olimitu);
		
		minspan = 100./(nsamples + 1);
		// write nbest Toffset in pool for each sample		
		ofpscratch = fopen("Toffset_worst.in", "w");
		if (offset_flag < 2)  //first send info about plot limits
		{
			fprintf(ofpscratch,"0\t%ld\n",static_cast<long>(olimitu));		
		}
		else
		{
			fprintf(ofpscratch,"%ld\t%ld\n",static_cast<long>(olimitl),static_cast<long>(olimitu));	;		
		}		
		for (imonte = poolsize - nworst + 1; imonte <= poolsize; imonte++)
		{
			for (n = 1; n < nsamples; n++)
			{
				fprintf(ofpscratch,"%6.2f\t%6.2f\t",n * minspan,Toffset[imonte][n]);
			}
			fprintf(ofpscratch,"%6.2f\t%6.2f\n",nsamples * minspan,Toffset[imonte][nsamples]);		
		}
		fclose(ofpscratch);
		
		minspan = 100./(nsamples + 1);
		// write nbest Toffset in pool for each sample		
		ofpscratch = fopen("Toffset_best.in", "w");
		if (offset_flag < 2)  //first send info about plot limits
		{
			fprintf(ofpscratch,"0\t%ld\n",static_cast<long>(olimitu));		
		}
		else
		{
			fprintf(ofpscratch,"%ld\t%ld\n",static_cast<long>(olimitl),static_cast<long>(olimitu));	;		
		}			
		for (imonte = 1; imonte <= nbest; imonte++)
		{
			for (n = 1; n < nsamples; n++)
			{
				fprintf(ofpscratch,"%6.2f\t%6.2f\t",n * minspan,Toffset[imonte][n]);
			}
			fprintf(ofpscratch,"%6.2f\t%6.2f\n",nsamples * minspan,Toffset[imonte][nsamples]);		
		}
		fclose(ofpscratch);
		
		ofpscratch = fopen("Toffset_labels.in", "w");
		for (n = 1; n <= nsamples; n++)
		{
			if (offset_flag < 2)
			{
				fprintf(ofpscratch,"%6.2f\t%4.1f\ts%1d\n",n * minspan,0.83 * olimitu, n);
			}
			else
			{
				fprintf(ofpscratch,"%6.2f\t%4.1f\ts%1d\n",n * minspan, olimitl + 0.83 * (olimitu-olimitl), n);
			}			
		}	
		fclose(ofpscratch);	
				
		// write nworst ages in pool for each minerals		
		ofpscratch = fopen("reference_minworst.in", "w");			
		for (imonte = poolsize - nworst + 1; imonte <= poolsize;  imonte++)
		{
			for (n = 1; n < nminerals[1]; n++)
			{
				fprintf(ofpscratch,"%6.2f\t%6.2f\t",n * minspan,mineral_pool[1][imonte][n]);
			}
			fprintf(ofpscratch,"%6.2f\t%6.2f\n",nminerals[1] * minspan,mineral_pool[1][imonte][nminerals[1]]);		
		}
		fclose(ofpscratch);		

// Now work on master age plot. For legacy reasons and because it's most likely Arvert users will have MDD data
// treat first sample differently, Plus, the first sample is the reference sample.
	// age plot background  - need this for all plots, even just minerala
		ofpscratch = fopen("age_backdrop.in", "w");		
		fprintf(ofpscratch,"%6.2f\t%6.2f\n", 0.0,0.0);		
		fprintf(ofpscratch,"%6.2f\t%6.2f\n", 100.,0.0);
		fprintf(ofpscratch,"%6.2f\t%d\n", 100.,static_cast<long>(agelimit));		
		fprintf(ofpscratch,"%6.2f\t%d\n", 0.0,static_cast<long>(agelimit));			
		fclose(ofpscratch);

		if ( (use_minerals < 2) && (nspectra[1] > 0) )  // we're using spectra and we have one for the first sample
		{
		// age spectrum general
			ofpscratch = fopen("reference_spec_general.in", "w");			
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][1] * 0.0);
			for (imonte = nbest+1; imonte <= poolsize - nworst - 1; imonte++)
			{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][1]);
			}
			fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize - nworst][1]);	

			for (ijj = 1; ijj <= labsteps[1] - 1; ijj++)
			{
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				for (imonte = nbest+1; imonte <= poolsize - 1 - nworst; imonte++)
				{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][ijj]);
				}
				fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize - nworst][ijj]);
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				for (imonte = nbest+1; imonte <= poolsize - 1 - nworst; imonte++)
				{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][ijj+1]);
				}
				fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize - nworst][ijj+1]);
			}
			fprintf(ofp,"%6.2f\t", bulkloss39[1][labsteps[1]] * 100.0);
			for (imonte = nbest+1; imonte <= poolsize - 1 - nworst; imonte++)
			{	
				fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][labsteps[1]]);
			}
			fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize - nworst][labsteps[1]]);		
			fclose(ofpscratch);		

		// age spectrum worst
			ofpscratch = fopen("reference_spec_worst.in", "w");			
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][1] * 0.0);
			for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1 ; imonte++)
			{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][1]);
			}
			fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize][1]);	

			for (ijj = 1; ijj <= labsteps[1] - 1; ijj++)
			{
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1 ; imonte++)
				{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][ijj]);
				}
				fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize][ijj]);
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1 ; imonte++)
				{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][ijj+1]);
				}
				fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize][ijj+1]);
			}
			fprintf(ofp,"%6.2f\t", bulkloss39[1][labsteps[1]] * 100.0);
			for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1; imonte++)
			{	
				fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][labsteps[1]]);
			}
			fprintf(ofpscratch,"%6.2f\n", poolAge[1][poolsize][labsteps[1]]);		
			fclose(ofpscratch);	
	
		// age spectrum best
			ofpscratch = fopen("reference_spec_best.in", "w");			
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][1] * 0.0);
			for (imonte = 1; imonte <= nbest - 1; imonte++)
			{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][1]);
			}
			fprintf(ofp,"%6.2f\n", poolAge[1][nbest][1]);	

			for (ijj = 1; ijj <= labsteps[1] - 1; ijj++)
			{
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				for (imonte = 1; imonte <= nbest - 1; imonte++)
				{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][ijj]);
				}
				fprintf(ofpscratch,"%6.2f\n", poolAge[1][nbest][ijj]);
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				for (imonte = 1; imonte <= nbest - 1; imonte++)
				{	
					fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][ijj+1]);
				}
				fprintf(ofpscratch,"%6.2f\n", poolAge[1][nbest][ijj+1]);
			}
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][labsteps[1]] * 100.0);
			for (imonte = 1; imonte <= nbest - 1; imonte++)
			{	
				fprintf(ofpscratch,"%6.2f\t", poolAge[1][imonte][labsteps[1]]);
			}
			fprintf(ofpscratch,"%6.2f\n", poolAge[1][nbest][labsteps[1]]);		
			fclose(ofpscratch);	

	// write lower and upper portions of observed age spectrum
			long goal_lines = 0;
			bulkloss39[1][0] = 0.0;
			ofpscratch = fopen("goal_spectrum.in", "w");
	
			if (goal_skip[1][1] < 1)
			{
				fprintf(ofpscratch,"> -W2.0,75/75/75\n");
				goal_lines++;
			}
			else
			{
				fprintf(ofpscratch,"> -W2.0,0/0/255\n");
				goal_lines++;		
			}
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][0] * 100.0);
			fprintf(ofpscratch,"%6.2f\n", goal_age[1][1] - goal_error[1][1]);
			goal_lines++;
			for (ijj = 1; ijj <= labsteps[1] - 1; ijj++)
			{			
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj - 1] * 100.0);
				fprintf(ofpscratch,"%6.2f\n", goal_age[1][ijj] - goal_error[1][ijj]);
				goal_lines++;
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				fprintf(ofpscratch,"%6.2f\n", goal_age[1][ijj] - goal_error[1][ijj]);
				goal_lines++;
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				fprintf(ofpscratch,"%6.2f\n", goal_age[1][ijj + 1] - goal_error[1][ijj + 1]);
				goal_lines++;
				if (goal_skip[1][ijj + 1] < 1)
				{
					fprintf(ofpscratch,"> -W2.0p,75/75/75\n");
					goal_lines++;
				}
				else
				{
					fprintf(ofpscratch,"> -W2.0p,0/0/255\n");
					goal_lines++;		
				}
			}
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][labsteps[1]] * 100.0);				
			fprintf(ofpscratch,"%6.2f\n", goal_age[1][labsteps[1]] - goal_error[1][labsteps[1]]);
			goal_lines++;

	// plot upper bounds of goal spectrum - we send pen info to scratch file not direct to gmt call
			if (goal_skip[1][1] < 1)
			{
				fprintf(ofpscratch,"> -W2.0p,75/75/75\n");
				goal_lines++;
			}
			else
			{
				fprintf(ofpscratch,"> -W2.0p,0/0/255\n");
				goal_lines++;		
			}
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][0] * 100.0);
			fprintf(ofpscratch,"%6.2f\n", goal_age[1][1] + goal_error[1][1]);
			goal_lines++;
			for (ijj = 1; ijj <= labsteps[1] - 1; ijj++)
			{			
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj - 1] * 100.0);
				fprintf(ofpscratch,"%6.2f\n", goal_age[1][ijj] + goal_error[1][ijj]);
				goal_lines++;
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				fprintf(ofpscratch,"%6.2f\n", goal_age[1][ijj] + goal_error[1][ijj]);
				goal_lines++;
				fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][ijj] * 100.0);
				fprintf(ofpscratch,"%6.2f\n", goal_age[1][ijj + 1] + goal_error[1][ijj + 1]);
				goal_lines++;
				if (goal_skip[1][ijj + 1] < 1)
				{
					fprintf(ofpscratch,"> -W2.0p,75/75/75\n");
					goal_lines++;
				}
				else
				{
					fprintf(ofpscratch,"> -W2.0p,0/0/255\n");
					goal_lines++;		
				}
			}
			fprintf(ofpscratch,"%6.2f\t", bulkloss39[1][labsteps[1]] * 100.0);				
			fprintf(ofpscratch,"%6.2f", goal_age[1][labsteps[1]] + goal_error[1][labsteps[1]]);	
			goal_lines++;		
			fclose(ofpscratch);	
		}
		
	// write reference minerals data
	
		minspan = 100./(nminerals[1] + 1);
		if ((use_minerals > 0) && (nminerals[1] > 0))  // we're using minerals and have some for the first sample
		{	
			// write goal mineral age(s)
			ofpscratch = fopen("reference_mingoal.in", "w");	
			for (n = 1; n <= nminerals[1]; n++)
			{
				fprintf(ofpscratch,"%6.2f\t%6.2f\n",n * minspan,mineral_goal[1][n]);
			}
			fclose(ofpscratch);
			// write nbest ages in pool for each minerals		
			ofpscratch = fopen("reference_minbest.in", "w");			
			for (imonte = 1; imonte <= nbest; imonte++)
			{
				for (n = 1; n < nminerals[1]; n++)
				{
					fprintf(ofpscratch,"%6.2f\t%6.2f\t",n * minspan,mineral_pool[1][imonte][n]);
				}
				fprintf(ofpscratch,"%6.2f\t%6.2f\n",nminerals[1] * minspan,mineral_pool[1][imonte][nminerals[1]]);		
			}
			fclose(ofpscratch);
			// write nworst ages in pool for each minerals		
			ofpscratch = fopen("reference_minworst.in", "w");			
			for (imonte = poolsize - nworst + 1; imonte <= poolsize;  imonte++)
			{
				for (n = 1; n < nminerals[1]; n++)
				{
					fprintf(ofpscratch,"%6.2f\t%6.2f\t",n * minspan,mineral_pool[1][imonte][n]);
				}
				fprintf(ofpscratch,"%6.2f\t%6.2f\n",nminerals[1] * minspan,mineral_pool[1][imonte][nminerals[1]]);		
			}
			fclose(ofpscratch);			

		}
		
// END OF WRITING for reference sample. That was easy by comparison to...

// NOW WRITE DATA FOR MEMBER SAMPLES IF WE HAVE THEM	
		if (nsamples > 1)
		{
			if (use_minerals < 2)
			{	
			// age spectrum best
				ofpscratch = fopen("member_specbest.in", "w");
				
				for (n = 1; n <= nsamples; n++)
				{
					if (nspectra[n] > 0) // write steps per spectrum to first line
					{
						fprintf(ofpscratch,"%ld\t", 2 * labsteps[n]);
					}
					else
					{
						fprintf(ofpscratch,"%ld\t", 0);		
					}
				}

				fprintf(ofpscratch,"\n");
				for (n = 2; n <= nsamples; n++)
				{
					if (nspectra[n] > 0) // write spectrum
					{
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][1] * 0.0);				
						for (imonte = 1; imonte <= nbest - 1; imonte++)
						{	
								fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][1]);
						}
						fprintf(ofp,"%6.2f\n", poolAge[n][nbest][1]);	

						for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
						{
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							for (imonte = 1; imonte <= nbest - 1; imonte++)
							{	
								fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][ijj]);
							}
							fprintf(ofpscratch,"%6.2f\n", poolAge[n][nbest][ijj]);
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							for (imonte = 1; imonte <= nbest - 1; imonte++)
							{	
								fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][ijj+1]);
							}
							fprintf(ofpscratch,"%6.2f\n", poolAge[n][nbest][ijj+1]);
						}
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);
						for (imonte = 1; imonte <= nbest - 1; imonte++)
						{	
							fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][labsteps[n]]);
						}
						fprintf(ofpscratch,"%6.2f\n", poolAge[n][nbest][labsteps[n]]);
					}
				}	
				fclose(ofpscratch);	
		
			// age spectrum worst
				ofpscratch = fopen("member_specworst.in", "w");
				for (n = 1; n <= nsamples; n++)
				{
					if (nspectra[n] > 0) // write steps per spectrum to first line
					{
						fprintf(ofpscratch,"%ld\t", 2 * labsteps[n]);			
					}
					else
					{
						fprintf(ofpscratch,"%ld\t", 0);					
					}
				}
				fprintf(ofpscratch,"\n");
				for (n = 2; n <= nsamples; n++)
				{
					if (nspectra[n] > 0) // write spectrum
					{
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][1] * 0.0);				
						for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1; imonte++)
						{	
								fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][1]);
						}
						fprintf(ofp,"%6.2f\n", poolAge[n][poolsize][1]);	

						for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
						{
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1 ; imonte++)
							{	
								fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][ijj]);
							}
							fprintf(ofpscratch,"%6.2f\n", poolAge[n][poolsize][ijj]);
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1; imonte++)
							{	
								fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][ijj+1]);
							}
							fprintf(ofpscratch,"%6.2f\n", poolAge[n][poolsize][ijj+1]);
						}
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);
						for (imonte = poolsize - nworst + 1; imonte <= poolsize - 1; imonte++)
						{	
							fprintf(ofpscratch,"%6.2f\t", poolAge[n][imonte][labsteps[n]]);
						}
						fprintf(ofpscratch,"%6.2f\n", poolAge[n][poolsize][labsteps[n]]);
					}
				}	
				fclose(ofpscratch);	
		// write lower and upper portions of observed age spectrum
				ofpscratch = fopen("member_specgoal.in", "w");
				for (n = 2; n <= nsamples; n++)
				{
					if (nspectra[n] > 0) // write goal spectra
					{
						goal_lines = 0;
						bulkloss39[n][0] = 0.0;

	
						if (goal_skip[n][n] < 1)
						{
							fprintf(ofpscratch,"> -W2.0,75/75/75\n");
							goal_lines++;
						}
						else
						{
							fprintf(ofpscratch,"> -W2.0,0/0/255\n");
							goal_lines++;		
						}
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][0] * 100.0);
						fprintf(ofpscratch,"%6.2f\n", goal_age[n][2] - goal_error[n][1]);
						goal_lines++;
						for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
						{			
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj - 1] * 100.0);
							fprintf(ofpscratch,"%6.2f\n", goal_age[n][ijj] - goal_error[n][ijj]);
							goal_lines++;
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							fprintf(ofpscratch,"%6.2f\n", goal_age[n][ijj] - goal_error[n][ijj]);
							goal_lines++;
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							fprintf(ofpscratch,"%6.2f\n", goal_age[n][ijj + 1] - goal_error[n][ijj + 1]);
							goal_lines++;
							if (goal_skip[n][ijj + 1] < 1)
							{
								fprintf(ofpscratch,"> -W2.0p,75/75/75\n");
								goal_lines++;
							}
							else
							{
								fprintf(ofpscratch,"> -W2.0p,0/0/255\n");
								goal_lines++;		
							}
						}
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);				
						fprintf(ofpscratch,"%6.2f\n", goal_age[n][labsteps[n]] - goal_error[n][labsteps[n]]);
						goal_lines++;

				// plot upper bounds of goal spectrum - we send pen info to scratch file not direct to gmt call
						if (goal_skip[n][n] < 1)
						{
							fprintf(ofpscratch,"> -W2.0p,75/75/75\n");
							goal_lines++;
						}
						else
						{
							fprintf(ofpscratch,"> -W2.0p,0/0/255\n");
							goal_lines++;		
						}
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][0] * 100.0);
						fprintf(ofpscratch,"%6.2f\n", goal_age[n][1] + goal_error[n][1]);
						goal_lines++;
						for (ijj = 1; ijj <= labsteps[n] - 1; ijj++)
						{			
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj - 1] * 100.0);
							fprintf(ofpscratch,"%6.2f\n", goal_age[n][ijj] + goal_error[n][ijj]);
							goal_lines++;
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							fprintf(ofpscratch,"%6.2f\n", goal_age[n][ijj] + goal_error[n][ijj]);
							goal_lines++;
							fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][ijj] * 100.0);
							fprintf(ofpscratch,"%6.2f\n", goal_age[n][ijj + 1] + goal_error[n][ijj + 1]);
							goal_lines++;
							if (goal_skip[n][ijj + 1] < 1)
							{
								fprintf(ofpscratch,"> -W2.0p,75/75/75\n");
								goal_lines++;
							}
							else
							{
								fprintf(ofpscratch,"> -W2.0p,0/0/255\n");
								goal_lines++;		
							}
						}
						fprintf(ofpscratch,"%6.2f\t", bulkloss39[n][labsteps[n]] * 100.0);				
						fprintf(ofpscratch,"%6.2f\n", goal_age[n][labsteps[n]] + goal_error[n][labsteps[n]]);	
						goal_lines++;
						goallines[n] = goal_lines;
					} // end if-nspectra[]
				} // loop through nsamples for goal spectra
				fclose(ofpscratch);
		
				ofpscratch = fopen("mgoallines.in", "w");	// quick and dirty to make plotting script simpler		
				for (n=1; n <= nsamples; n++)
				{
					if (nspectra[n] > 0) // write goal spectrum lines
					{				
						fprintf(ofpscratch,"%ld\t", goallines[n]);
					}
					else
					{				
						fprintf(ofpscratch,"%ld\t", 0);
					}					
				}
				fclose(ofpscratch);	
		
			}	// end of use minerals < 2
		
			if (use_minerals > 0)
			{	
				ofpscratch = fopen("member_mingoal.in", "w");	
				for (n = 2; n <= nsamples; n++)
				{
					if (nminerals[n] > 0) // write goal ages
					{
						minspan = 100./(nminerals[n] + 1);
						// write goal ages for each minerals; 

						for (i = 1; i <= nminerals[n]; i++)
						{
							fprintf(ofpscratch,"%6.2f\t%6.2f\n",i * minspan,mineral_goal[n][i]);
						}
					}
				}
				fclose(ofpscratch);
				
				ofpscratch = fopen("member_minbest.in", "w"); // write nminerals for ALL samples
				for (n = 1; n <= nsamples; n++)
				{
					fprintf(ofpscratch,"%ld\t",nminerals[n]);				
				}				
				fprintf(ofpscratch,"\n");

				for (n = 2; n <= nsamples; n++)
				{				
					if (nminerals[n] > 0) // write best ages for members
					{
						minspan = 100./(nminerals[n] + 1);
						for (imonte = 1; imonte <= nbest; imonte++)
						{
							for (i = 1; i < nminerals[n]; i++)
							{
								fprintf(ofpscratch,"%6.2f\t%6.2f\t",i * minspan,mineral_pool[n][imonte][i]);
							}
							fprintf(ofpscratch,"%6.2f\t%6.2f\n",nminerals[n] * minspan,mineral_pool[n][imonte][nminerals[n]]);		
						}
					}
				}
				fclose(ofpscratch);
				
				ofpscratch = fopen("member_minworst.in", "w"); // write nminerals for ALL samples
				for (n = 1; n <= nsamples; n++)
				{
					fprintf(ofpscratch,"%ld\t",nminerals[n]);				
				}				
				fprintf(ofpscratch,"\n");

				for (n = 2; n <= nsamples; n++)
				{				
					if (nminerals[n] > 0) // write worst ages for members
					{
						minspan = 100./(nminerals[n] + 1);
						for (imonte = poolsize - nworst + 1; imonte <= poolsize;  imonte++)
						{
							for (i = 1; i < nminerals[n]; i++)
							{
								fprintf(ofpscratch,"%6.2f\t%6.2f\t",i * minspan,mineral_pool[n][imonte][i]);
							}
							fprintf(ofpscratch,"%6.2f\t%6.2f\n",nminerals[n] * minspan,mineral_pool[n][imonte][nminerals[n]]);		
						}
					}
				}
				fclose(ofpscratch);
			}
		} // END OF MEMBER IF STATEMENT
		
// announce status of plot commands
		printf("\n  ...plot data written.\n");		

 // create a pdf plot for viewing
 		printf("\nLaunching plot script - plotting might take a minute...\n");
		gmtcall = "./plot_arvert.sh";
		system(gmtcall.c_str());
 		printf("  ...plot created. Model done.\n\n");		
	}  // end of gmt-plotting block

// cleanup

    return 0;
}
