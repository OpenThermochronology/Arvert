/* 	

		Sorting function sort() of Arvert 7.0
		P. Zeitler, Lehigh University
		based on 6.1.3 code
		
		sort() arranges PoolAge, PoolTemp, Hepool, and ModelFit in order: best to worst fit
		
		4/2004 -- added sorting of mineral-age data
		
		11/2014 - switched away from global variables; added option for multiple mineral ages
		
		December 2021  altered to work with multi-sample Arvert 7
  	
*/

#include "arvert.h"

void sort(long nsamples, long nspectra[],long poolsize, long labsteps[MAXLABSTEPS], long ttpoints, double poolAge[][MAXPOOLSIZE][MAXLABSTEPS], double poolTemp[MAXPOOLSIZE][NODES], double poolTime[MAXPOOLSIZE][NODES], double ModelFit[MAXPOOLSIZE], long use_minerals, long nminerals[MAXSAMPLES], double mineral_pool[MAXSAMPLES][MAXPOOLSIZE][MAX_MINERALS], double Toffset[MAXPOOLSIZE][MAXSAMPLES])
{
	long i, j, k, imineral, n;
	double swap;
	
	for (i = 1; i <= poolsize - 1; i++)
	{
		for (j = i + 1; j <= poolsize; j ++)
		{
			if (ModelFit[j] < ModelFit[i])
			{
				// swap fit
				swap = ModelFit[i];
				ModelFit[i] = ModelFit[j];
				ModelFit[j] = swap;

				// swap all the mineral ages in the mineral-age pools		
				if (use_minerals > 0)
				{					
					for (n=1; n <= nsamples; n++)
					{
						if(nminerals[n] > 0)
						{
							for (imineral = 1; imineral <= nminerals[n]; imineral++)
							{
								swap = mineral_pool[n][i][imineral];
								mineral_pool[n][i][imineral] = mineral_pool[n][j][imineral];
								mineral_pool[n][j][imineral] = swap;
							}
						}
					}
				}
				if (use_minerals < 2)
				{					
					for (n=1; n <= nsamples; n++)
					{				
						if(nspectra[n] > 0)
						{
							// swap all the spectra in the age-spectrum pool
							for (k = 1; k <= labsteps[n]; k++)
							{
								swap = poolAge[n][i][k];
								poolAge[n][i][k] = poolAge[n][j][k];
								poolAge[n][j][k] = swap;
							}
						}
					}
				}
				
				// swap all the temperature nodes in the history
				for (k = 1; k <= ttpoints; k++)
				{
					swap = poolTemp[i][k];
					poolTemp[i][k] = poolTemp[j][k];
					poolTemp[j][k] = swap;
				}
				// swap all the time nodes in the history
				for (k = 1; k <= ttpoints; k++)
				{
					swap = poolTime[i][k];
					poolTime[i][k] = poolTime[j][k];
					poolTime[j][k] = swap;
				}
				
				// swap all the Toffsets

				for (n=2; n <= nsamples; n++)
				{
					swap = Toffset[i][n];
					Toffset[i][n] = Toffset[j][n];
					Toffset[j][n] = swap;
				}				
			}
		}
	}
	
	return;
}