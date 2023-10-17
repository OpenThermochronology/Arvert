/* 	

		Averaging routine averages() of Arvert 7.0,x
		P. Zeitler, Lehigh University
		December 2021 modified from 6.1.3. code
		
		averages() determines average and bracketing thermal history and average age spectrum
		
		11/2014 now pass parameters by reference and value, not globals; modified for multiple mineral ages
  	
  		December 2021 altered to handle multiple samples
*/

#include "arvert.h"

void averages(long nsamples, long nspectra[], long ttpoints, long poolsize, long labsteps[MAXSAMPLES], double poolAge[MAXSAMPLES][MAXPOOLSIZE][MAXLABSTEPS], double poolTemp[MAXPOOLSIZE][NODES], long use_minerals, long nminerals[MAXSAMPLES], double mineral_pool[MAXSAMPLES][MAXPOOLSIZE][MAX_MINERALS], double average_mineral_age[MAXSAMPLES][MAX_MINERALS], double avTemp[NODES], double avAge[][MAXLABSTEPS], double lowEnvelope[NODES], double highEnvelope[NODES])
{
	long i, j, n;
	
	for (i = 1; i <= ttpoints; i++)
	{
		lowEnvelope[i] = 1000.;
		highEnvelope[i] = -100.;
		avTemp[i] = 0.0; 
		for (j = 1; j <= poolsize; j++)
		{
			avTemp[i] = avTemp[i] + poolTemp[j][i];
			if (lowEnvelope[i] > poolTemp[j][i])
			{
				lowEnvelope[i] = poolTemp[j][i];
			}
			if (highEnvelope[i] < poolTemp[j][i])
			{
				highEnvelope[i] = poolTemp[j][i];
			}
		}
		avTemp[i] = avTemp[i]/(poolsize);
	}

	if (use_minerals < 2)
	{
		for (n = 1; n <= nsamples; n++)
		{
			if (nspectra[n] > 0)
			{
				for (i = 1; i <= labsteps[n]; i++)
				{
					avAge[n][i] = 0.0;
					for (j = 1; j <= poolsize; j++)
					{
						avAge[n][i] = avAge[n][i] + poolAge[n][j][i];
					}
					avAge[n][i] = avAge[n][i]/(poolsize);
				}
			}
		}
	}

	if (use_minerals > 0)
	{	
		for (n = 1; n <= nsamples; n++)
		{	
			if (nminerals[n] > 0)
			{
				for (i = 1; i <= nminerals[n]; i++)
				{
					average_mineral_age[n][i] = mineral_pool[n][1][i];
					for (j = 2; j <= poolsize; j++)
					{
						average_mineral_age[n][i] = average_mineral_age[n][i] + mineral_pool[n][j][i];
					}
					average_mineral_age[n][i] = average_mineral_age[n][i] / poolsize;
				}
			}
		}
	}

	return;
}