/* 	
		Functions selectsubset() and randomselect() of Arvert 7.0 and later
		P. Zeitler, Lehigh University
		8/6/2000
		
		The former is a key routine for the constrained-random-search algorithm. It selects
		a random subset of thermal histories from the main pool of histories. The
		selected histories are unique--there is no duplication of the selections.
		
		Function randomselect() just picks a random number in a range equivalent to the
		number of histories in the pool.
		
		11/2014  - removed global variables
		
		January 2022   vetted for Arvert 7
		
		A recent 2021 paper suggests ( I think ) that rather than a random selection like this, that
		the best-fit history be used. This should speed convergence... but would it hamper exploration?
		Quick test suggests this does not work well, so this routine left unchanged.
*/

#include "arvert.h"
//#include "nr3.h"
//#include "ran_pz.h"

// ------- Function selectsubset() -------------------------

long selectsubset(long subsetsize, long CRSselections[MAXPOOLSIZE], long poolsize)
{
	long i, j, duplicate, CRS_selection;

	CRSselections[1] = randomselect(poolsize);
	for (i = 2; i <= subsetsize + 1; i++)
	{
		do
		{
			CRSselections[i] = randomselect(poolsize);
			duplicate = 0;
			for (j = 1; j <= i - 1; j++)
			{
				if (CRSselections[i] == CRSselections[j])
				{
					duplicate = 1;
					break;
				}
			}
		}
		while (duplicate == 1);
	}
	
	CRS_selection = CRSselections[subsetsize + 1];  // index of history to be modified using CRS algorithm

	return CRS_selection;
}


// ------- Function randomselect() -------------------------

long randomselect(long poolsize)
{
	long selection;
	Ranq1 ran(rand());
		
	selection = poolsize + 1;
	while (selection > poolsize)  // guard against unlikely case that ran3() = 1, making selection i_carlos + 1
	{
		selection = 1 + floor(ran.doub() * (poolsize));
	}
	
	return selection;
}