// this is for arvert 7.0.x
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
//#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <algorithm>
#include <vector>

#include "nr3.h"
#include "ran_pz.h"

// following block is for lovera() routine
	#define     nxi 50001
	#define     nch  8005
	#define     ntp   100
	#define      pi     3.141592653589793
	#define  xlambd     0.000543
	#define       R     0.001987
	#define 	  ecam  50
	#define 	 elab  60
	#define 	euplus  50
	#define    NODES 101    //number of tT nodes plus spare 

// stopping distances for U-He volume-diffusion ages in mineral_age_vd() 
	#define stopping_apatite 18.81e-4   // for 238U, ketcham et al., 2011
	#define stopping_zircon 15.55e-4
	#define stopping_titanite 17.46e-4

// miscellaneous macros for program operation
	#define MAXHEADROOM	0.85  // maximum step or mineral age as fraction of modelduration
	#define MAX_MINERALS 51  // protect the foolish from themselves
	#define MAXLINE 101  // for file reading
	#define MINSUBSET  5  // this is very low but just protecting against disaster
	
	#define MAXSAMPLES  11  // max number of samples to invert (10; +1 due to 1-based arrays)
	#define MAXDOMAINS  11  // max number of domains per MDD sample (10; +1 due to 1-based arrays)
	#define MAXLABSTEPS  101  // max number of age-spectrum steps(100) plus one more
	#define MAXPOOLSIZE 302 // max size of MC and CRS pool (300; plus 2, for swapping and for 1-based arrays)
	
// function prototypes

void averages(long nsamples, long nspectra[], long ttpoints, long poolsize, long labsteps[MAXSAMPLES], double poolAge[MAXSAMPLES][MAXPOOLSIZE][MAXLABSTEPS], double poolTemp[MAXPOOLSIZE][NODES], long use_minerals, long nminerals[MAXSAMPLES], double mineral_pool[MAXSAMPLES][MAXPOOLSIZE][MAX_MINERALS], double average_mineral_age[MAXSAMPLES][MAX_MINERALS], double avTemp[NODES], double avAge[][MAXLABSTEPS], double lowEnvelope[NODES], double highEnvelope[NODES]);

long bestfit(long poolsize, double ModelFit[MAXPOOLSIZE]);

double fit(long nsamples, long nspectra[], long fitchoice, long index, const double poolAge[][MAXPOOLSIZE][MAXLABSTEPS], long goalsteps[], const double goal_loss[][MAXLABSTEPS], const double goal_age[][MAXLABSTEPS], const double goal_error[][MAXLABSTEPS], const long goal_skip[][MAXLABSTEPS], long use_minerals, long nminerals[], const double mineral_goal[MAXSAMPLES][MAX_MINERALS], const double mineral_goalerr[MAXSAMPLES][MAX_MINERALS], const double mineral_weight[MAXSAMPLES][MAX_MINERALS], const double mineral_pool[][MAXPOOLSIZE][MAX_MINERALS], const double early_loss[], const double early_weight[]);

void heatsched(long nsamples, long nspectra[], long ndomains[], long labsteps[], long geometry, long goalsteps[], double Eact[][MAXDOMAINS], double D0[][MAXDOMAINS], double fraction[][MAXDOMAINS], double goal_loss[][MAXLABSTEPS], double timelab[][MAXLABSTEPS], double temperaturelab[][MAXLABSTEPS]);

void geolovera(long sample, long ttpoints, long geometry, double bgeom, double slab, long ndomains, double Eact[MAXDOMAINS], double D0[MAXDOMAINS], double fraction[MAXDOMAINS], long labsteps, double dt, double ichange, double timein[], double temperaturein[], const double timelab[MAXLABSTEPS],const double temperaturelab[MAXLABSTEPS], double bulkloss39[][MAXLABSTEPS], double age[][MAXLABSTEPS]);

void lablovera(long sample, long ttpoints, long geometry, double bgeom, double slab, long ndomains, double Eact[MAXDOMAINS], double D0[MAXDOMAINS], double fraction[MAXDOMAINS], long labsteps, double dt, double ichange, const double timelab[MAXLABSTEPS],const double temperaturelab[MAXLABSTEPS], double bulkloss39[][MAXLABSTEPS]);

double mineral_age_vd(long ndL, long ndt, long mineral, double radius, double eact, double difzero, double modelduration, const double timein[], const double temperaturein[], long ttpoints); 

void monteTt(long rttpoints, long ttpoints, double monteHeatRate, double monteCoolRate, long nconstraints, double constraintTime[NODES], double constraintTempmin[NODES], double constraintTempmax[NODES], double MonteTime[NODES], double MonteTemp[NODES], double temperaturein[NODES], double timein[NODES], long flips);

long randomselect(long poolsize);

void newhistoryCRS(long ttpoints, long &goodHistory, long CRS_selection, long subsetsize, const long CRSselections[], double amplify, double coolRate, double heatRate, const double poolTemp[][NODES], double upperlimit[], double lowerlimit[], double timein[], double temperaturein[], long flips, double Toffset[][MAXSAMPLES], double newToffset[], long nsamples, long offset_flag, double offset_max, double offset_min);

long selectsubset(long subsetsize, long CRSselections[MAXPOOLSIZE], long poolsize);

void sort(long nsamples, long nspectra[], long poolsize, long labsteps[], long ttpoints, double poolAge[][MAXPOOLSIZE][MAXLABSTEPS], double poolTemp[][NODES], double poolTime[][NODES], double ModelFit[], long use_minerals, long nminerals[], double mineral_pool[][MAXPOOLSIZE][MAX_MINERALS], double Toffset[][MAXSAMPLES]);

long worstfit(long poolsize, const double ModelFit[]);

