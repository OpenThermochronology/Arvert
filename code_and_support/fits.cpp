/* 	

		Fitting functions fit(), worstfit(), bestfit() of Arvert 7.0.x
		P. Zeitler, Lehigh University
		2/16/2004  Update January 2022
		
		fit() can use either mean percent deviation or mswd
		
		mean percent deviation is just the arithmetic average of the absolute value
		of the percent deviation between measured and goals steps, for each fitted step;
		for minerals it just measured and predicted age.
		
		2/16/04 -- routine changed to calculate mswd rather than chi-square as the second fitting option,
		and added bestfit() routine
  		4/04 -- added code for optional fitting of He ages
  		4/12 -- fitmax variable is an unused legacy of attempt to initially focus on early heating steps and
  		then open up to all heating steps. Didn't work well. Maybe better to do this using different
  		version of goal-spectrum fitting and restart option
  		11/2014 -- removed use of global variables Accommodated multiple mineral ages.
  		11.24.2015 -- removed fitmax variable
  		11.28.2015  need to modify this to better deal with weighting. Let's try this. Treat each age-spectrum step and mineral age
  		equally. Then calculate sum of deviations with mineral ages weighted such that 1 is baseline; greater than 1
  		amplifies weight, less than 1 de-amplifies. Then get a modified MSWD where mean is weighted mean that includes weights.
  		Also tried a repair to option 0 - just weight all early steps below a cut-off higher (rather than using step size as well).
  		
  		January 2022: Changes to accomodate new options to have data combos of different data types, from multiple 
  		samples. There are two main options: fitting each sample separately and combining them, or fitting all steps and
  		ages as one pool of observables (related by Toffset). 
  		
  		Decided to start with lumping all observables into one single set, rather than sample by sample. Mineral
  		weighting intended to factor relative to one heating step.
  		
  		Note that this means that samples with more MDD steps and more mineral ages get more weight. That sort
  		of makes sense, since we know more about it.
  		
  		As of version 7.0.x this routine is functioning but has NOT been deeply tested for correctness w.r.t different
  		combos of weighting, data types.
*/

#include "arvert.h"

// ********** function fit()
// key routine for determining objective function for inversion

double fit(long nsamples, long nspectra[MAXSAMPLES], long fitchoice, long index, const double poolAge[][MAXPOOLSIZE][MAXLABSTEPS], long goalsteps[], const double goal_loss[][MAXLABSTEPS], const double goal_age[][MAXLABSTEPS], const double goal_error[][MAXLABSTEPS], const long goal_skip[][MAXLABSTEPS], long use_minerals, long nminerals[], const double mineral_goal[][MAX_MINERALS], const double mineral_goalerr[][MAX_MINERALS], const double mineral_weight[][MAX_MINERALS], const double mineral_pool[][MAXPOOLSIZE][MAX_MINERALS], const double early_loss[], const double early_weight[])
{
	long i, j, fittedSteps, n;
	double spec_fit, mineral_fit, total_fit, deviation, mineral_deviation, stepweight, sumweight, sum_mineralweight;

	deviation = 0.0;
	fittedSteps = 0;
	stepweight = 0.0;
	sumweight = 0.0;
	sum_mineralweight = 0.0;

	for (n=1; n <= nsamples; n++)
	{

	// first do fits if we are using MDD age spectra
		if (use_minerals < 2)  // we want to fit an age spectrum but using only steps flagged as ok
		{
			if (nspectra[n] > 0)
			{
				for (j = 1; j<= goalsteps[n]; j++)
				{
					if (goal_skip[n][j] > 0)
					{
						fittedSteps = fittedSteps + 1;
						switch (fitchoice)
						{
							case 0:  // step-weighted mean percent deviation
								if (goal_loss[n][j] < early_loss[n])   // empiricially weight some early steps more given that they have more thermal relief
								{
									stepweight = early_weight[n];
								}
								else     // just weight all steps the same
								{
									stepweight = 1.0;
								}
								sumweight = sumweight + stepweight;
								deviation = deviation + stepweight*fabs(poolAge[n][index][j] - goal_age[n][j])/goal_age[n][j];
							break;
							case 1:  // mean percent deviation
								deviation = deviation + fabs(poolAge[n][index][j] - goal_age[n][j])/goal_age[n][j];		
							break;
				
							case 2:  // mswd
								deviation = deviation + 1/(goal_error[n][j]*goal_error[n][j])*pow((poolAge[n][index][j] - goal_age[n][j]),2.0);
							break;
						}
					}
				}
				assert(fittedSteps >= 1);
			} // end if-nspectra[]
		}

	// now do fits for any minerals ages if they're being used
		if (use_minerals > 0)
		{	
			if (nminerals[n] > 0)
				{
				if (fitchoice < 2)  // 1 = mean percent deviation
				{
					for (i = 1; i <= nminerals[n]; i++)
					{
						mineral_deviation = mineral_deviation + mineral_weight[n][i] * fabs(mineral_pool[n][index][i] - mineral_goal[n][i]) / mineral_goal[n][i];
						sum_mineralweight = sum_mineralweight + mineral_weight[n][i];
					}
				}
				else                // 2 = mswd-like parameter
				{
					for (i = 1; i <= nminerals[n]; i++)
					{
						mineral_deviation = mineral_deviation + mineral_weight[n][i] * pow((mineral_pool[n][index][i] - mineral_goal[n][i]),2.0) / pow(mineral_goalerr[n][i],2.0);
						sum_mineralweight = sum_mineralweight + mineral_weight[n][i];
					}
				}
			}	
		}	
	}
	
	switch (use_minerals)
	{
		case 0:   // only age spectrum
			if (fitchoice < 2) // mpd
			{
				total_fit = 100.0 * deviation / sumweight;		
			}
			else
			{
				total_fit = (double) (deviation / fittedSteps);		
			}
			break;
		case 1:   // spectra and mineral ages
			if (fitchoice < 2) // mpd
			{
				total_fit = 100.0 * (deviation + mineral_deviation) / (sumweight + sum_mineralweight);		
			}
			else
			{
				total_fit = (double) ((deviation + mineral_deviation) / (fittedSteps + sum_mineralweight));		
			}				
			break;
		case 2:  // mineral ages
			if (fitchoice < 2) // mpd
			{
				total_fit = 100.0 * mineral_deviation / sum_mineralweight;		
			}
			else
			{
				total_fit = (double) (mineral_deviation / sum_mineralweight);		
			}			
			break;	
	}
	
	return total_fit;
}

// ********** function worstfit()
// important routine for CRS algorithm since we have to continually find out which thermal history is worst 

long worstfit(long poolsize, const double ModelFit[MAXPOOLSIZE])
{	
	long ii, WorstOne;
	
	WorstOne = 1;
	for (ii = 2; ii <= poolsize; ii++)
	{
		if (ModelFit[ii] > ModelFit[WorstOne])
		{
			WorstOne = ii;
		}
	}
	return WorstOne;
}


// ********** function bestfit()
// important routine for CRS algorithm since we have to continually find out which thermal history is best 

long bestfit(long poolsize, double ModelFit[MAXPOOLSIZE])
{	
	long ii, BestOne;
	
	BestOne = 1;
	for (ii = 2; ii <= poolsize; ii++)
	{
		if (ModelFit[ii] < ModelFit[BestOne])
		{
			BestOne = ii;
		}
	}
	return BestOne;
}
