/* 	

		Function heatsched() of Arvert 7.0.x
		P. Zeitler, Lehigh University
		modified from version 6.1.3 code
		
   Generates a heating schedule based on the losses observed for the goal spectrum.
   
		It is easiest to create model age spectra which have step sizes that match the goal spectrum, rather than
		use the actual lab heating schedule but end up with steps which don't match. In an ideal world, of course,
		the lab schedule would give you the same steps. But it's not an ideal world because the domain distribution and
		kinetics are not perfect (this opens a philosophical question,  left for another day, about whether this
		inversion should start with the domains determination as part of it. Ye gods.).
		
		It's visually easier to compare results, and fitting is far more straightforward, if the goal and
		model spectra have the same step end points.
	
  		~2014: Removed global variables and now pass parameters using mix of by value and by reference, as needed.
  		
  		CRITICAL SUMMARY THAT IT IS VITAL TO UNDERSTAND: Arvert 7.0.x calculates a heating schedule that exactly predicts
  		the observed losses. It DOES NOT USE the sample's actual heating schedule.
  		
  		June 2020: patch, to address a malfunction in odd cases (e.g. old Newark Basin data that jump from 0,55 to 0.99 for
  		final two steps. Raised second guess[2] value from 3000 to 10000 for all steps but first.
		
		December 2021: as part of Avert 7, alter routine to process multiple samples _within_ this routine.
		DIfferent MDD samples need different synthetic heating schedules.
		
		Note that this only needs to be called once, at start, since all tT models will share the same set of lab
		heating schedules.
		
		In a future release (Arvert 7.1.x?) will add a routine that gives the user the option to use the actual losses
		predicted from the domain structure and actual lab heating schedule. This will be complicated because of the
		need to change how fitting of spectra is done, plus will need to keep an eye on what the loss predictions
		actually look like.
*/

#include "arvert.h"

void heatsched(long nsamples, long nspectra[MAXSAMPLES], long ndomains[MAXDOMAINS], long labsteps[MAXLABSTEPS], long geometry, long goalsteps[MAXLABSTEPS], double Eact[][MAXDOMAINS], double D0[][MAXDOMAINS], double fraction[][MAXDOMAINS], double goal_loss[][MAXLABSTEPS], double timelab[][MAXLABSTEPS], double temperaturelab[][MAXLABSTEPS])
{
	long step, i, j, n;
	double target[MAXLABSTEPS],cumDta2[21],guess[4],Arloss[4];
	double Dta2, a_term, volume, delta;

	for (n=1; n <= nsamples; n++)
	{
		if (nspectra[n] > 0)
		{
			labsteps[n] = goalsteps[n]; // legacy from olden times when people were asked to input the heating schedule
	
			for (i = 1; i <= 10; i++)
			{
				cumDta2[i] = 0.0;
			}

			for (i = 1; i <= labsteps[n]; i++)
			{
				target[i] = goal_loss[n][i];
				assert((target[i] <= 1.0) && (target[i] != target[i-1]));
				if (target[i] == 1.0) //   protect against case of total loss being entered by a dodo.
				{
					target[i] = 0.9999;
				}
			}

			delta = 1.0;
			guess[1] = 10; // these guesses MUST bracket real temperature value for first step!! Let's hope these do!
			guess[2] = 3000;
			do                    //  deal with first step, for which there is no previous 39Ar loss.
			{
				guess[3] = (guess[1] + guess[2])/2;
				for (i = 1; i <= 3; i++) 
				{
					Arloss[i] = 0.0;
					volume = 0;
					for (j = 1; j <= ndomains[n]; j++)
					{
						volume = volume + fraction[n][j];
		//arbitrarily make simulated steps 10 minutes long				
						a_term = D0[n][j]/31557600e6 * 600;  // Do/a2 back to 1/sec, then times 10 minutes to get Dot/a2

						cumDta2[j] = a_term * exp(-Eact[n][j] / R / (guess[i] + 273.15));
						if (geometry == 1) // spherical
						{
							if (cumDta2[j] <= 0.000001)
							{
								Arloss[i] = Arloss[i] + fraction[n][j]*(6 * sqrt(cumDta2[j]/pi));
							}
							else
							{
								long index;
								double sumseries, previous, change;
								sumseries = 0.0;
								previous = 1.0;
								index = 0;
								do
								{
									index = index + 1;
									change = fabs(previous - sumseries)/previous;
									previous = sumseries;
									sumseries = sumseries + 6/pi/pi/index/index*exp(-index*index*pi*pi*cumDta2[j]);
								}
								while (change > 0.000000001);
								Arloss[i] = Arloss[i] + fraction[n][j]*(1 - sumseries);
							}
						}
						else  // slab
						{
							if (cumDta2[j] < 0.159)
							{
								Arloss[i] = Arloss[i] + fraction[n][j]*(1.12838 * sqrt(cumDta2[j]));
							}
							else
							{
								Arloss[i] = Arloss[i] + fraction[n][j]*(1 - 0.810569 * exp(-2.467401 * cumDta2[j]));
							}
						}  // end if-geometry
					}   // end for-j to domains					
			
					Arloss[i] = Arloss[i]/volume;

				} // for loop calculating all guess temperatures 
		
				if (Arloss[3] <= target[1])
				{
					guess[1] = guess[3];
				}
				else
				{
					guess[2] = guess[3];
				}
				delta = fabs((Arloss[3] - target[1]))/target[1];
			}
			while (delta > 0.001);  // end conditional loop iterating guess temperatures 
	
			temperaturelab[n][1] = guess[3] + 273.15;
	
			for (j = 1; j <= ndomains[n]; j++)
			{
				a_term = D0[n][j]/31557600e6 * 600;   // Dt/a2 for ten minutes
				cumDta2[j] = a_term * exp(-Eact[n][j] / R / (guess[3] + 273.15));
			}
		
			for (step = 2; step <= labsteps[n]; step++)   /// now deal with all other steps
			{
				delta = 1.0;
				guess[1] = 10.;  // for the record, guess temps work in C
				guess[2] = 10000.;
				do
				{
					guess[3] = (guess[1] + guess[2])/2;
					for (i = 1; i <= 3; i++)
					{
						Arloss[i] = 0.0;
						for (j = 1; j <= ndomains[n]; j++)
						{
							a_term = D0[n][j]/31557600e6 * 600;
							Dta2 = cumDta2[j] + a_term * exp(-Eact[n][j] / R / (guess[i] + 273.15));
		//	printf("step %d   domain %d   Dta2 %e\n",step,j,Dta2);
						if (geometry == 1)
						{
							if (Dta2 <= 0.000001)
							{
								Arloss[i] = Arloss[i] + fraction[n][j]*(6 * sqrt(Dta2/pi));
							}
							else
							{
								long index;
								double sumseries, previous, change;
								sumseries = 0.0;
								previous = 1.0;
								index = 0;
								do
								{
									index = index + 1;
									change = fabs(previous - sumseries)/previous;
									previous = sumseries;
									sumseries = sumseries + 6/pi/pi/index/index*exp(-index*index*pi*pi*Dta2);
								}
								while (change > 0.000000001);
								Arloss[i] = Arloss[i] + fraction[n][j]*(1 - sumseries);
							}
						}
						else
							{
								if (Dta2 < 0.159)
								{
									Arloss[i] = Arloss[i] + fraction[n][j]*(1.12838 * sqrt(Dta2));
								}
								else
								{
									Arloss[i] = Arloss[i] + fraction[n][j]*(1 - 0.810569 * exp(-2.46740 * Dta2));
								}
							}
						}
						Arloss[i] = Arloss[i]/volume;
					} //  loop through guess temperatures
			
		// printf("step %d Arloss[1 - 3] is %f   %f   %f\n",step, Arloss[3],target[step],guess[3]);
					if (Arloss[3] <= target[step])
					{
						guess[1] = guess[3];
					}
						else
					{
						guess[2] = guess[3];
					}
					delta = fabs(Arloss[3] - target[step])/target[step];
				}
				while (delta > 0.001);
	
				temperaturelab[n][step] = guess[3] + 273.15;  // convert degrees C to K
				for (j = 1; j <= ndomains[n]; j++)
				{
					a_term = D0[n][j]/31557600e6 * 600;   // Dt/a2 for ten minutes
					cumDta2[j] = cumDta2[j] + a_term * exp(-Eact[n][j] / R / (guess[3] + 273.15));
				}
			}  

			for (i = 1; i <= labsteps[n]; i++)
			{
				timelab[n][i] = 10/5.2596e11;  // 10 minutes converted to units of m.y. 
			}
		}
	} // end of n-loop through samples
return;
}