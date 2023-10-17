/* 	

		Functions selectsubset() and randomselect() of Arvert 7.0.1
		P. Zeitler, Lehigh University
		8/6/2000
		
		This is THE key routine for the constrained-random-search algorithm. It creates a new
		thermal history using the subset of thermal histories plus one additional thermal
		history. The new thermal history is created by projecting the temperature values of
		the one selected history through the centroid of the subset of histories; this projection
		is amplified by an amount specified by the input variable <amplify>, which typically
		has a value of about 1.3 (values too close to 1.00 don't sufficiently explore parameter
		space, and values much larger slow convergence because they lead to wild histories).
		Example: if amplify is 1.3, the temperature at a time node is 100, and the centroid (mean) value for the pool is
		130 at the same time node, then the temperature for the new CRS history at this time
		mode will be (130 - 100)* 1.3 + 130, or 169.
		
		Note that (as pointed out by Sean Willet in his CRS code for fission-track annealing), 
		it might make sense to use a weighted mean in creating the centroid for the subset, the
		weighting being based on the exponential temperature dependence of diffusion. At present,
		however, the code below uses just an arithmetic mean to create the centroid values.
		
		This routine works in four parts:
		
		(1) determines centroid values for subset
		
		(2) makes new CRS history

	  	(3) checks to make sure this history is reasonable
	  		(a) repairs any temperatures that have sneaked across explicit bounds
	  		(b) makes no attempt to repair violations of implicit constraints -- history is rejected
	  	
		(4) loads this new history into array temperaturein[nn] for use in lovera()
		
		3/15/03 -- narrowed explicit constraints
		
		3/30/03 -- restored implicit rate constraints. No repair attempt if history illegal--just retry whole process.

		April 2 2012: Fixed "bug" (meaning weird, scary, non-understood, probably cancerous issue) in which
		     slow heating rates near -1C/m.y. were not being caught out as bad if no heating allowed.
		     Code was: if ((rate < 0.0) && abs(rate) > heatRate)); this was NOT trapping the bad rate if
		     it was between -1 and 0. Eh??? Fixed by just saying if (rate < -heatRate). But why?????

		24 October 2014 to November 15, 2014. Note that above "bug" was a real bug, because I used abs() not fabs()... 
		abs() returs int not real! You need to use labs() for long and fabs() for double. So if rate was > -1 and < 0, 
		we were not getting the expected result.
		     
		Made changes to parse thermal history and reject it if there are > 3 reversals of direction. 
		This might lead to problems with higher numbers of nodes because too many attempts will be taken:
		code will be slower and might break limit on # of CRS attempts (may have to loosen this).
		
		11/2014 removed global variables
		
		January 2022 vetted for Arvert 7.. Note that a more recent 2021 paper might suggest that the new history
		should be built not from a random member of pool, but the best fit. Good for convergence, but iffy for 
		exploration? Quick test suggested this is not a useful innovation, so things left as they are. 
		
		February 2022. added code for Toffset being part of CRS.
		
		3/11/2022 experimented with exponential average of subset temperatures (used kspar E, actual temps in K).
		NO big diff - dropped in June 2022.
		
		15 June, 2022: Realized that code written to honor Toffset constraints was mung, Deleted it but introduced a local
		off_amplify that's hard-coded to 1.05, to make the offset CRS routine be less pushy so less likely to
		produce strong excursions in value.
*/

#include "arvert.h"

// ------- Function newhistoryCRS() -------------------------

void newhistoryCRS(long ttpoints, long &goodHistory, long CRS_selection, long subsetsize, const long CRSselections[MAXPOOLSIZE], double amplify, double coolRate, double heatRate, const double poolTemp[MAXPOOLSIZE][NODES], double upperlimit[NODES], double lowerlimit[NODES], double timein[NODES], double temperaturein[NODES], const long flips, double Toffset[MAXPOOLSIZE][MAXSAMPLES], double newToffset[MAXSAMPLES], long nsamples, long offset_flag, double offset_max, double offset_min)
{
	long i, j, n, index;
	double centroid[NODES],rate;
	double offset_centroid[MAXSAMPLES];
	double off_amplify = 1.05;
		
	long flip, sign1, sign2;

// first, determine the centroid values for the subset

	for (i = 1; i <= ttpoints; i++)  // zero centroid[12] array
	{
		centroid[i] = 0.0;
	}
	
	for (n = 2; n <= nsamples; n++)  // zero centroid[12] array
	{
		offset_centroid[n] = 0.0;
	}

	for (i = 1; i <= ttpoints; i++) // loop through the selected histories to get sum at each time node
	{
		for (j = 1; j <= subsetsize; j++)
		{
			index = CRSselections[j];
			centroid[i] = centroid[i] + poolTemp[index][i]; 
//			centroid[i] = centroid[i] + exp(-7500./1.987/(273.15 + poolTemp[index][i]));		// TEST 031122	
		}
	}

	for (n = 2; n <= nsamples; n++) // loop through the selected histories to get sum for each sample
	{
		for (j = 1; j <= subsetsize; j++)
		{
			index = CRSselections[j];
			offset_centroid[n] = offset_centroid[n] + Toffset[index][n];
		}
	}
	
	for  (i = 1; i <= ttpoints; i++)   // get mean at each time node
	{
//		centroid[i] = centroid[i] / subsetsize;  
		centroid[i] = centroid[i] / subsetsize;  		// TEST 031122	
//		centroid[i] = -7500./1.987/log(centroid[i]) - 273.15;  		// TEST 031122	
	}

	for  (n = 2; n <= nsamples; n++)   // get mean for each sample
	{
		offset_centroid[n] = offset_centroid[n] / subsetsize;
	}

// now, actually build the new CRS history
	for  (i = 1; i <= ttpoints; i++) // loop through the time nodes
	{
		temperaturein[i] = centroid[i] + amplify * (centroid[i] - poolTemp[CRS_selection][i]);
	}
	
// build the new CRS offsets if Toffsets are part of inversion
// hard code amplify (off_amplify) to a more modest 1.05 to avoid stomping too badly on Toffset constraints
	newToffset[1] = 0.0;
	if (offset_flag > 0)
	{
		for  (n = 2; n <= nsamples; n++)   // new CRS Toffset (added amplify back in (might be too large so that new histories get smacked against limits)
		{
			newToffset[n] = offset_centroid[n] + off_amplify * (offset_centroid[n] - Toffset[CRS_selection][n]);
		}
// decided to drop following code, and dishonor user constraints
/*
	// honor Toffset constraints according to how offsets used (mono progressive or random). Cap at limit if limit exceeded.
		if (offset_flag == 1) // mono progressive
		{
			for (n=2; n <= nsamples; n++)		
			{
				if ((newToffset[n] - newToffset[n-1]) > offset_max) newToffset[n] = newToffset[n-1] + offset_max;
				if ((newToffset[n] - newToffset[n-1]) < offset_min) newToffset[n] = newToffset[n-1] + offset_min;					
			}			
		}
		else  // offset_flag must be 2 (random)
		{
			for (n=2; n <= nsamples; n++)		
			{
				if (newToffset[n] > offset_max) newToffset[n] = offset_max;
				if (newToffset[n] < offset_min) newToffset[n] = offset_min;			
			}
		}
*/
	}
	else
	{
		for (n = 2; n <= nsamples; n++)   // just keep the input values
		{
			newToffset[n] = Toffset[1][n];
		}
	}

// next, check the new history for points that exceed explicit bounds
// Deal with these extreme points by setting them just inside the explicit limits

	for  (i = 1; i <= ttpoints; i++) // loop through the time nodes
	{
		if (temperaturein[i] > upperlimit[i])
		{
			temperaturein[i]= upperlimit[i] - 1.0;
		}
		if (temperaturein[i] < lowerlimit[i])
		{
			temperaturein[i]= lowerlimit[i] + 1.0;
		}
	}
	
// now, check for violations of implicit rate constraints, and either try to fix them or exit routine in failure and retry from scratch

	goodHistory = 1;

	for  (i = 1; i <= ttpoints - 1; i++) // loop through the time nodes
	{
		rate = (temperaturein[i] - temperaturein[i+1])/(timein[i] - timein[i+1]);   // positive is cooling
		
		if (rate > coolRate)
		{
			goodHistory = -1;
		}	
		if (rate < -heatRate)
		{
			goodHistory = -1;
		}
	}
	
// now check for too many flips in the sign of the temperature rate
	flip = 0;
	for  (i = 1; i <= ttpoints-2; i++) // loop through the time nodes
	{	
		sign1 = ((temperaturein[i] - temperaturein[i+1]) <= 0 ) ? 1 : -1;
		sign2 = ((temperaturein[i+1] - temperaturein[i+2]) <= 0 ) ? 1 : -1;
		if ((sign1 * sign2) < 0) { flip = flip + 1;}
	}

	if (flip > (flips + 2) )  // key line here: allow more oscillation with larger values if you want them
	{
		goodHistory = -flip;
	}

	return;
}