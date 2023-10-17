/* 	

		Function monteTt() of Arvert 7.0.x
		P. Zeitler, Lehigh University
		11/2014 Major revision to older code, to damp tT oscillation during any reheating
		
		5/2016  Major revision to generate temps at randomized time nodes, and then pass regularly sampled version
		of these to the CRS algorithm. I.e., we make a thermal history with a limited number of random time nodes
		and then at a much finer time grid, we let the CRS process recombine these. This better honors the match between
		the number of constraints, the number of histories in the subset, and the number of time nodes, but produces
		smoother and more realistic final CRS tT histories.
		
		Generates a random time-temperature history that obeys heating/cooling contraints (explicit temperature limits,
		implicit heating/cooling rate). For each trial, time nodes are chosen at random (other than start and end).
		
		23 October 2014: add changes that allow history with at most one reheating episode (e.g. cool-heat-cool)
		
		20 October 2014:  mostly working but buggy: not working with nodes above 12 or so. WTF?
		
		Early November 2014: This is proving to be a tremendous pain in the ass to get working properly. I grow old,
		I grow old, I shall wear the bottom of my trousers rolled... 
		
		15 November 2014: Gotta move on to other work. I am sure there is some cleverer and far more compact way of
		constraining tT reversals, but I'm going the route of more code so what's happening is more clear. One 
		problem is dealing with the start and end nodes, if we randomly choose to start with one of those. So we'll
		just handle each set of cases separately (start nodes, end nodes, and then the general case for nodes that
		are in between). I am also seeing REALLY weird behavior where logical tests based on longs are not working
		which leads me to think that (??) some sort of implicit type conversion is happening?? Whatever that it is, 
		it is NOT helping with debugging this code.
		
		Will need to do much testing as we allow this to enter the main Arvert fork (at least, make sure it does no harm to 
		models that only use monotonic cooling!). Seems to work ok.
		
		Switched to newer Ranq1 random-number generator from Numerical Recipes -- now need to include nr.h and ran_pz.h
		and declare an instance of the Ranq1 struct, then call doub() function from this struct to get a random between 
		0 and 1.
		
		May 2016 - Major change to introduce random location of time nodes, and then output these histories to be  
		projected into a finer equally spaced time mesh for the CRS work. Thus, user now needs to tell Arvert how many 
		MC and how many CRS nodes they want (there should be more CRS than MC modes; I guess they could be equal but 
		have not tested this).
		
		May 9, 2020: 
		Arvert 6.1.3: cleaned up comments, tidied up location of some blocks of code, and in this routine slightly 
		changed the criteria for random times being too close to start or finish, from 1% to 2%.
		
		January 2022 checked as part of Arvert 7.x multi-sample. No changes needed other than NODES macro.
		January 26, 2022. Discovered through experimentation with lovera() that start temp needs to be _well_ above
		max Tc to avoid failure of lovera() routine. Decided to remove protection from this routine and handle it up
		front in main().
		
		Also discovered that this routine is not honoring constraints at start (and elsewhere). Short of debugging
		the whole thing in its complexity, tweaked so that first point always honors external limits for first bounding 
		point.
		
		Not enough, this fix... discovered mistaken use of || instead of && in a compound if-clause. Jeez - somehow
		this swine of a routine worked in Arvert 6. 
		
		Also added coded to main() that determines Tc10 of largest domain in all samples and then adds 200 C to that
		to ensure safe start for lovera() (based on some testing with the lovera() forward model). This overrides user model-
		start constraint if necessary and prints a warning.
		
		Let as few as 3 MC points be used (on realization that since time nodes are randomized, using fewer
		MC nodes might be a better seed for CRS algorithm, with less hassles related to preventing extra flips).
		
		Seems to be working now (January 28th).
		
		10 March 2022. Taking hint from Dale Issler, experiment with using tangent transform in 
		translating rate constraints to constraining creating of new MC temperatures. For slow rates, not a
		big issue, but if higher rates allowed, there is large difference between mean rate allowed versus
		mean rate created; adding this unskews this bias towards faster cooling. Didn't get this working in a way that improved things.
*/

#include "arvert.h"

void monteTt(long rttpoints, long ttpoints, double monteHeatRate, double monteCoolRate, long nconstraints, double constraintTime[NODES], double constraintTempmin[NODES], double constraintTempmax[NODES], double MonteTime[NODES], double MonteTemp[NODES], double temperaturein[NODES], double timein[NODES], const long flips)
{
	long i, j, start, index;
	double highlimit,lowlimit, ratelimit;
	double lowerlimit[NODES], upperlimit[NODES];  // NODES is mild overkill but simpler then yet another MACRO
	
	bool started_heating, now_heating, previously_heating;
	long flip;
	
	Ranq1 ran(rand());

// NOTE: MonteTime[] stores time coords of Monte Carlo nodes in problem frame of reference.

// For this Arvert monteTt.cpp variant, we want at least 5 tT nodes. This should be checked in main(), but double-check
// here:

	assert(rttpoints >= 3 && "insufficient tT points (< 3) in MonteTt for this Arvert variant");	

// For Arvert 6+, we now implement random time nodes at the MC stage. This goes a long way to eliminating the location 
// of time nodes as being an obscure but significant implicit constraint on when the tT history can change slope.
// E.g., // imagine needing a short fast cooling pulse at 15 Ma, but only having fixed nodes at 20 and 10 Ma.

// That means we have to make the random nodes. The simplest way seems to be to first generate the desired number, 
// then sort what we've created.

// We need to make rttpoints - 2 time nodes, since the ends are fixed at modelduration and 0.0.
// Looks like we're using 1-based arrays (thought I had purged those...dream on)

	double localtime[rttpoints+1], localtemp[rttpoints+1];
	
	double swap;
	long duplicate;
	
	localtime[1] = MonteTime[1];
	localtime[rttpoints] = MonteTime[ttpoints];

//make rttpoints-2 time nodes, making sure they are different than start or end

	duplicate = -1;
	do
	{
		for (i = 2; i < rttpoints; i++)
		{
			localtime[i] = ran.doub() * MonteTime[1];
			if (localtime[i] < 0.02 * MonteTime[1]) localtime[i] = 0.02 * MonteTime[1];
			if (localtime[i] > 0.98 * MonteTime[1]) localtime[i] = 0.98 * MonteTime[1];
		}

	// sort the time nodes	
		for (i = 1; i <= rttpoints-1; i++)
		{
			for (j = i+1; j <= rttpoints; j++)
			{
				if (localtime[j] > localtime[i])
				{
					swap = localtime[i];
					localtime[i] = localtime[j];
					localtime[j] = swap;
				}
			}		
		}

// check for nodes being too close (might need to finesse what "too close" means!)
// Right now, we keep nodes apart by at least 2% of total duration

		for (i = 1; i < rttpoints; i++)
		{
			if ((localtime[i] - localtime[i+1]) < 0.02 * localtime[1])
			{
				duplicate = 1;
			}
		}
	} while (duplicate < 0);

// Now need to set up upper and lower limits for this local set of MC time nodes

	upperlimit[1] = constraintTempmax[1];
	lowerlimit[1] = constraintTempmin[1];
	upperlimit[rttpoints] = constraintTempmax[nconstraints];
	lowerlimit[rttpoints] = constraintTempmin[nconstraints];	
	for (i = 2; i < rttpoints; i++)
	{
		j = 1;
		do
		{
			j++;
		}
		while (localtime[i] < constraintTime[j]);
		
		j = j - 1;
		
		upperlimit[i] = constraintTempmax[j] - (constraintTime[j] - localtime[i])/(constraintTime[j] - constraintTime[j+1])*(constraintTempmax[j] - constraintTempmax[j+1]);
		lowerlimit[i] = constraintTempmin[j] - (constraintTime[j] - localtime[i])/(constraintTime[j] - constraintTime[j+1])*(constraintTempmin[j] - constraintTempmin[j+1]);
	}

/*	
Making a thermal history: to avoid bias, this version of the original code retains Sean Willet's scheme where we 
start at a randomly selected node and work in both directions.	However, the MAJOR DIFFERENCE in this code is that
if reheating is allowed, we limit reversals in the tT path to basically one reheating event: we disallow
oscillations. The user needs to be EXTREMELY CAREFUL in watching how this added constraint on the original
Monte-Carlo pool interacts with explicit and implicit constraints, potentially producing bias to the solution.
Of course, the knuckleheaded user should ALWAYS be EXTREMELY CAREFUL when specifying constraints!

As noted above, we eschew elegance for clarity, so we use separate blocks of similar code to work through cases where
the start node is at or within one step of start or end of model - these positions complicate initializing initial 
trends. 
*/

// lock in first temperature point no matter what since we need to protect lovera()
	localtemp[1] = lowerlimit[1] + ran.doub() * (upperlimit[1] - lowerlimit[1]);

// first get a randomly selected starting node

	start = 1 + floor(ran.doub() * (rttpoints));
	if (start > rttpoints) start = rttpoints; // in the very unlikely case that random exactly = 1

// get a  random temperature for this starting point
	localtemp[start] = lowerlimit[start] + ran.doub() * (upperlimit[start] - lowerlimit[start]);


// In the interests of pool diversity, let's declare that it's more important to start randomly from
// all nodes than it is to slavishly lock in a single reheating pulse. This simplifies the code, but
// might allow a bit of a tT jiggle at the beginning or end.

// node-1 case
	if (start == 1)  // already made node 1, now just work to model end
	{
		// generate first younger point out of loop, in order to get initial trend efficiently

		ratelimit = localtemp[1] + monteHeatRate*(localtime[1] - localtime[2]);  //first, use implicit limit
		highlimit = (ratelimit < upperlimit[2]) ? ratelimit : upperlimit[2];

		ratelimit = localtemp[1] - monteCoolRate*(localtime[1] - localtime[2]);  //next, use implicit limit
		lowlimit = (ratelimit > lowerlimit[2]) ? ratelimit : lowerlimit[2];
		localtemp[2] = lowlimit + ran.doub() * (highlimit - lowlimit);

		// now get initial trend
		started_heating = (localtemp[2] >= localtemp[1]) ? true : false; //if segment dT is positive, we are heating; set flag to true
		previously_heating = started_heating;

		flip = 0;
		for (index = 3; index <= rttpoints; index++)
		{
			ratelimit = localtemp[index-1] + monteHeatRate*(localtime[index-1] - localtime[index]);  //first, determine limit set by heating rate
			highlimit = (ratelimit < upperlimit[index]) ? ratelimit : upperlimit[index]; // choose smaller of ratelimit and explicit temperature boundary

			ratelimit = localtemp[index-1] - monteCoolRate*(localtime[index-1] - localtime[index]);  //as above but do lower limit
			lowlimit = (ratelimit > lowerlimit[index]) ? ratelimit : lowerlimit[index];

			if (flip <= flips)  // generate temperature anywhere in full range
			{
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);
			}
			else  // generate temperature in a range that is restricted by the requirement to stop flipping; details complex
			{
				(started_heating) ? (lowlimit = localtemp[index-1]) : (highlimit = localtemp[index-1]);
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);			
			}
			previously_heating = (localtemp[index-1] >= localtemp[index-2]) ? true : false;
			now_heating = (localtemp[index] >= localtemp[index-1]) ? true : false;

			flip = ((previously_heating)==(now_heating)) ? (flip+0) : (flip+1); // check for flip to heating/cooling from cooling/heating; set flag
		}
//		if (localtemp[1] < lowerlimit[1]) {printf(" >>>>>>>>> lowerlimit violation created in block 1: %ld\t%lf\t%lf\n",start,lowerlimit[1],localtemp[1]);}
	}
	
// node-2 case
	if (start == 2)  // make up any old node[1] then work to model end, any extra flip be damned
	{
		// first make node[ttpoints] and get trend data

		ratelimit = localtemp[2] + monteCoolRate*(localtime[1] - localtime[2]);  //first, use implicit limit
		highlimit = (ratelimit < upperlimit[1]) ? ratelimit : upperlimit[1];

		ratelimit = localtemp[2] - monteHeatRate*(localtime[1] - localtime[2]);  //next, use implicit limit
		lowlimit = (ratelimit > lowerlimit[1]) ? ratelimit : lowerlimit[1];
//		localtemp[1] = lowlimit + ran.doub() * (highlimit - lowlimit);
		localtemp[1] = lowerlimit[1] + ran.doub() * (upperlimit[1] - lowerlimit[1]);
	
		// now get initial trend
		started_heating = (localtemp[2] >= localtemp[1]) ? true : false; //if segment dT is positive, we are heating; set flag to true
		previously_heating = started_heating;

		flip = 0;
		for (index = 3; index <= rttpoints; index++)
		{
			ratelimit = localtemp[index-1] + monteHeatRate*(localtime[index-1] - localtime[index]);  //first, determine limit set by heating rate
			highlimit = (ratelimit < upperlimit[index]) ? ratelimit : upperlimit[index]; // choose smaller of ratelimit and explicit temperature boundary

			ratelimit = localtemp[index-1] - monteCoolRate*(localtime[index-1] - localtime[index]);  //as above but do lower limit
			lowlimit = (ratelimit > lowerlimit[index]) ? ratelimit : lowerlimit[index];

			if (flip <= flips)  // generate temperature anywhere in full range
			{
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);
			}
			else  // generate temperature in a range that is restricted by the requirement to stop flipping; details complex
			{
				(started_heating) ? (lowlimit = localtemp[index-1]) : (highlimit = localtemp[index-1]);
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);			
			}
			previously_heating = (localtemp[index-1] >= localtemp[index-2]) ? true : false;
			now_heating = (localtemp[index] >= localtemp[index-1]) ? true : false;

			flip = ((previously_heating)==(now_heating)) ? (flip+0) : (flip+1); // check for flip to heating/cooling from cooling/heating; set flag
		}
//		if (localtemp[1] < lowerlimit[1]) {printf(" >>>>>>>>> lowerlimit violation created in block 2: %ld\t%lf\t%lf\n",start,lowerlimit[1],localtemp[1]);}
	}	

// next-to-last node case
	if (start == (rttpoints - 1) )  // make up any old node[ttpoints] then work to model start, any extra flip be damned
	{
		// generate first older point out of loop, in order to get initial trend efficiently

		ratelimit = localtemp[rttpoints - 1] + monteHeatRate*(localtime[rttpoints-1] - localtime[rttpoints]);  //first, use implicit  limit
		highlimit = (ratelimit < upperlimit[rttpoints]) ? ratelimit : upperlimit[rttpoints];

		ratelimit = localtemp[rttpoints - 1] - monteCoolRate*(localtime[rttpoints-1] - localtime[rttpoints]);  //next, use implicit  limit
		lowlimit = (ratelimit > lowerlimit[rttpoints]) ? ratelimit : lowerlimit[rttpoints];
		localtemp[rttpoints-1] = lowlimit + ran.doub() * (highlimit - lowlimit);

		// now get initial trend
		started_heating = (localtemp[rttpoints] >= localtemp[rttpoints-1]) ? true : false; //if segment dT is positive, we are heating; set flag to 1
		previously_heating = started_heating;

		flip = 0;
		for (index = rttpoints - 2; index > 1; index--)  //working backward, implicit cooling limit determines implicit max limit, and vice versa
		{
			ratelimit = localtemp[index+1] + monteCoolRate*(localtime[index] - localtime[index+1]);  //first, determine limit set by heating rate
			highlimit = (ratelimit < upperlimit[index]) ? ratelimit : upperlimit[index]; // choose smaller of ratelimit and explicit temperature boundary

			ratelimit = localtemp[index+1] - monteHeatRate*(localtime[index] - localtime[index+1]);  //as above but do lower limit
			lowlimit = (ratelimit > lowerlimit[index]) ? ratelimit : lowerlimit[index];

			if (flip <= flips)  // generate temperature anywhere in full range
			{
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);
			}
			else  // generate temperature in a range that is restricted by the requirement to stop flipping; also force monotonic cooling
			{
				(started_heating) ? (highlimit = localtemp[index+1]) : (lowlimit = localtemp[index+1]); // this forces monotonic cooling
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);			
			}
			previously_heating = (localtemp[index+1] <= localtemp[index+2]) ? true : false;
			now_heating = (localtemp[index] <= localtemp[index+1]) ? true : false;

			flip = ((previously_heating)==(now_heating)) ? (flip+0) : (flip+1); // check for flip to heating/cooling from cooling/heating; set flag
		}
//		if (localtemp[1] < lowerlimit[1]) {printf(" >>>>>>>>> lowerlimit violation created in block 3: %ld\t%lf\t%lf\n",start,lowerlimit[1],localtemp[1]);}
	}

// last-node case
	if ((start == rttpoints) && (start != 3) )  // just work back to model start
	{
	// generate temperatures at times older than starting point (LOWER node numbers)
		flip = 0;
		// generate first older point out of loop, in order to get initial trend efficiently

		ratelimit = localtemp[rttpoints] + monteCoolRate*(localtime[rttpoints-1] - localtime[rttpoints]);  //first, use implicit  limit
		highlimit = (ratelimit < upperlimit[rttpoints-1]) ? ratelimit : upperlimit[rttpoints-1];

		ratelimit = localtemp[rttpoints] - monteHeatRate*(localtime[rttpoints-1] - localtime[rttpoints]);  //next, use implicit  limit
		lowlimit = (ratelimit > lowerlimit[rttpoints-1]) ? ratelimit : lowerlimit[rttpoints-1];
		localtemp[rttpoints-1] = lowlimit + ran.doub() * (highlimit - lowlimit);

		// now get initial trend
		started_heating = (localtemp[rttpoints] >= localtemp[rttpoints-1]) ? true : false; //if segment dT is positive, we are heating; set flag to 1
		previously_heating = started_heating;

		for (index = rttpoints - 2; index > 1; index--)  //working backward, implicit cooling limit determines implicit max limit, and vice versa
		{
			ratelimit = localtemp[index+1] + monteCoolRate*(localtime[index] - localtime[index+1]);  //first, determine limit set by heating rate
			highlimit = (ratelimit < upperlimit[index]) ? ratelimit : upperlimit[index]; // choose smaller of ratelimit and explicit temperature boundary

			ratelimit = localtemp[index+1] - monteHeatRate*(localtime[index] - localtime[index+1]);  //as above but do lower limit
			lowlimit = (ratelimit > lowerlimit[index]) ? ratelimit : lowerlimit[index];

			if (flip <= flips)  // generate temperature anywhere in full range
			{
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);
			}
			else  // generate temperature in a range that is restricted by the requirement to stop flipping; also force monotonic cooling
			{
				(started_heating) ? (highlimit = localtemp[index+1]) : (lowlimit = localtemp[index+1]); // this forces monotonic cooling
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);			
			}
			previously_heating = (localtemp[index+1] <= localtemp[index+2]) ? true : false;
			now_heating = (localtemp[index] <= localtemp[index+1]) ? true : false;

			flip = ((previously_heating)==(now_heating)) ? (flip+0) : (flip+1); // check for flip to heating/cooling from cooling/heating; set flag
		}
//		if (localtemp[1] < lowerlimit[1]) {printf(" >>>>>>>>> lowerlimit violation created in block 4: %ld\t%lf\t%lf\n",start,lowerlimit[1],localtemp[1]);}
	}	

// general case	where there are more than five nodes and we are in the middle
	if ((start > 2) && (start < (rttpoints - 1) ))   // changed from ||
	{
	// generate temperatures at times older than starting point (LOWER node numbers)
		flip = 0;
		// generate first older point out of loop, in order to get initial trend efficiently

		ratelimit = localtemp[start] + monteCoolRate*(localtime[start-1] - localtime[start]);  //first, use implicit cooling limit
		highlimit = (ratelimit < upperlimit[start-1]) ? ratelimit : upperlimit[start-1];

		ratelimit = localtemp[start] - monteHeatRate*(localtime[start-1] - localtime[start]);  //next, use implicit heating limit
		lowlimit = (ratelimit > lowerlimit[start-1]) ? ratelimit : lowerlimit[start-1];
		localtemp[start-1] = lowlimit + ran.doub() * (highlimit - lowlimit);

		// now get initial trend
		started_heating = (localtemp[start] >= localtemp[start-1]) ? true : false; //if segment dT is positive, we are heating; set flag to 1
		previously_heating = started_heating;

		for (index = start - 2; index > 1; index--)  //working backward, implicit cooling limit determines implicit max limit, and vice versa
		{
			ratelimit = localtemp[index+1] + monteCoolRate*(localtime[index] - localtime[index+1]);  //first, determine limit set by heating rate
			highlimit = (ratelimit < upperlimit[index]) ? ratelimit : upperlimit[index]; // choose smaller of ratelimit and explicit temperature boundary

			ratelimit = localtemp[index+1] - monteHeatRate*(localtime[index] - localtime[index+1]);  //as above but do lower limit
			lowlimit = (ratelimit > lowerlimit[index]) ? ratelimit : lowerlimit[index];

			if (flip <= flips)  // generate temperature anywhere in full range
			{
			  localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);
			}
			else  // generate temperature in a range that is restricted by the requirement to stop flipping; also force monotonic cooling
			{
			  (started_heating) ? (highlimit = localtemp[index+1]) : (lowlimit = localtemp[index+1]); // this forces monotonic cooling
			  localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);			
			}
			previously_heating = (localtemp[index+1] <= localtemp[index+2]) ? true : false;
			now_heating = (localtemp[index] <= localtemp[index+1]) ? true : false;

			flip = ((previously_heating)==(now_heating)) ? (flip+0) : (flip+1); // check for flip to heating/cooling from cooling/heating; set flag
		}
		
	// now return to generate temperatures at times younger than starting point (LARGER node numbers)
		for (index = start + 1; index <= rttpoints; index++)
		{
			ratelimit = localtemp[index-1] + monteHeatRate*(localtime[index-1] - localtime[index]);  //first, determine limit set by heating rate
			highlimit = (ratelimit < upperlimit[index]) ? ratelimit : upperlimit[index]; // choose smaller of ratelimit and explicit temperature boundary

			ratelimit = localtemp[index-1] - monteCoolRate*(localtime[index-1] - localtime[index]);  //as above but do lower limit
			lowlimit = (ratelimit > lowerlimit[index]) ? ratelimit : lowerlimit[index];

			if (flip <= flips)  // generate temperature anywhere in full range
			{
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);
			}
			else  // generate temperature in a range that is restricted by the requirement to stop flipping; details complex
			{
				(started_heating) ? (lowlimit = localtemp[index-1]) : (highlimit = localtemp[index-1]);
				localtemp[index] = lowlimit + ran.doub() * (highlimit - lowlimit);			
			}
			previously_heating = (localtemp[index-1] >= localtemp[index-2]) ? true : false;
			now_heating = (localtemp[index] >= localtemp[index-1]) ? true : false;

			flip = ((previously_heating)==(now_heating)) ? (flip+0) : (flip+1); // check for flip to heating/cooling from cooling/heating; set flag
		}
//		if (localtemp[1] < lowerlimit[1]) {printf(" >>>>>>>>> lowerlimit violation created in block 5: %ld\t%lf\t%lf\n",start,lowerlimit[1],localtemp[1]);}
	}
	
// add some protection for Lovera routine even if it makes for a few fugly histories
// protection for lovera() now handle in main(), when setting up model limits.
	
// Need to interpolate our thermal history defined by rttpoints MC nodes to ttpoints CRS nodes. Ugh - this is never fun 
// - how many times have I struggled to write a routine like this?

	double node_spacing, temp_slope;

	node_spacing = localtime[1] / (ttpoints - 1);

	MonteTemp[1] = localtemp[1];
//	if (MonteTemp[1] < lowerlimit[1]) {printf(" >>>>>>>>> lowerlimit violation snuck through %lf\t%lf\n",MonteTemp[1],lowerlimit[1]);}
	MonteTemp[ttpoints] = localtemp[rttpoints];

	index = 2;
	for (i=2; i<=rttpoints; i++)  // step through segments
	{
		temp_slope = (localtemp[i] - localtemp[i-1])/(localtime[i-1] - localtime[i]);
		while (MonteTime[index] > localtime[i])
		{
			MonteTemp[index] = localtemp[i-1] + (localtime[i-1] - MonteTime[index]) * temp_slope;
			index++;
		}
	}

// need to interface history with input arrays for Lovera()
	for (index = 1; index <= ttpoints; index++)
	{
		timein[index] = MonteTime[index];
		temperaturein[index] = MonteTemp[index];
	}
	
	return;
}