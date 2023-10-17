/* 	

		Function geolovera() of Arvert 7.0.x
		P. Zeitler, Lehigh University
		modified from 6.1.3 code
		
   C then C++ translation of Oscar Lovera's fortan program ages.f (translated through intermediate
   Pascal translation, also by Zeitler!). Somehow it (still) works.
   
   Handles heating as well as cooling as long as, after the model start
   there is complete cooling from well above closure to below closure.
   
   2013 to Nov 2014: small changes to work under C++. Note use of std::vector. 
   
   11/2014 - Arvert now passes parameters by mix of value and reference, instead of globals
   
   N.B. a more recent standalone function, loverafunc(), is written in plain C. This code has included updates 
   added to that code. MOST importantly, unlike Oscar's original code, this tests both the 40Ar and 39Ar infinite
   series for fractional convergence (Oscar tested only 40Ar). Without this, low losses of perhaps 0.1% or less will
   not be accurately predicted. With this fix, all should be fine with the variable ichange set to 1e-7. With this value,
   code will be slower so if user is not modeling first few percent anyway, higher values can be tolerated. For
   intensive work it would pay to test values of ichange before commencing on lots of work. Thing to look for is that 
   predicted Ar losses match input losses in file goal.in over the range that is being modeled.
   
   December 2021 - in setting up for Arvert 7 and multiple samples, realized that at present (e.g., Arvert 6.1.3), where
   we use a fake heating schedule instead of the observed one, bulkloss values from Lovera() are identical to the input 			
   observed losses. So, there is no need to pass bulkloss39[] back to main(). This will simplify code. BUT if we 
   eventually move towards using both observed and predicted losses, we will need to communicate with lovera() about
   39Ar or extract this routine as a separate function. 
   We need to allow lovera() to calculate bulkloss39[] because its components are used for 40/39 ratio (but can we 
   excise this to save compute cycles?).
   
   Geolovera() is based on lovera() from Arvert 6.1.3.
   This routine just calculates the geological portion of 40Ar evolution, and needs to be paired with lablovera(), 
   which determines the laboratory release, which only needs to be calculated once per MDD samples within this inverse
   code.
*/

#include "arvert.h"
#include <vector>
using namespace std;

void geolovera(long sample, long ttpoints, long geometry, double bgeom, double slab, long ndomains, double Eact[MAXDOMAINS], double D0[MAXDOMAINS], double fraction[MAXDOMAINS], long labsteps, double dt, double ichange, double timein[], double temperaturein[], const double timelab[MAXLABSTEPS], const double temperaturelab[MAXLABSTEPS], double bulkloss39[][MAXLABSTEPS], double age[][MAXLABSTEPS])
{
	double r40[MAXLABSTEPS], tra40[MAXLABSTEPS];
	double sra40[MAXLABSTEPS], sra39[MAXLABSTEPS], scura[MAXLABSTEPS];
	double c40, c39wd;
	
	long j, ns, i;
	
//	long dummy;

// following are from chist procedure of Lovera

	double Tdiff, deltaTemp, deltaTime;
	long newpoints;
	long ij, nt;
	
// following are from geosec procedure of Lovera
	double zit, avtemperature, a1, uplus, cam, sum;
	long	n, jj, nj;
	long	madj, mfake, m;

// following are from lab procedure of Lovera
	double  kpie, prev_r40, change40, zii, ab;
	double zita[MAXLABSTEPS];
	long jjj, k;
	
//  clock_t started, ended;
//  double elapsed;

	std::vector<double> xi(50003);
	std::vector<double> temperature(nch);
	std::vector<double> time(nch);
	std::vector<double> dzita(nch);
	std::vector<double> d(nch);	

	//  ------------ Set some needed parameters -------------

	c40 = 1.0 - exp(-xlambd * timein[1]);
	c39wd = exp(-xlambd * timein[1]);

	ns = 1;
	temperature[1] = temperaturein[1];
	time[1] = timein[1];

// discretize the input thermal history--new PZ algorithm that uses Oscar's approach of discretizing
// over constant temperature intervals, not time intervals. This solves two problems: two many time nodes if
// there is little T change over large time spans, or very large T change over a fixed small time interval, that is,
// this approach is far better for highly variable cooling history with very high and low rates.

	for (k = 1; k <= ttpoints - 1; k++)
	{
		Tdiff = temperaturein[k] - temperaturein[k + 1];
		if (fabs(Tdiff) <= dt)  //if delta-T small have at least number of time and temp points given in input history
		{
			ns = ns + 1;
			temperature[ns] = temperaturein[k + 1];
			time[ns] = timein[k + 1];
		}
		else
		{
			newpoints = (long) floor(fabs(Tdiff)/dt + 0.5) + 1;
			deltaTemp = Tdiff/newpoints;
			deltaTime = (timein[k] - timein[k+1])/newpoints;
			for (ij = 1; ij <= newpoints; ij++)
			{
			temperature[ij + ns] = temperaturein[k] - ij*deltaTemp;
			time[ij + ns] = timein[k] - ij*deltaTime;
			}
			ns = ns + newpoints;
		}
		assert(ns <= nch - 2);
	}

// fix last point to that of input data (this avoids problems with interpolated non-rational times ending up just slightly negative)

	temperature[ns] = temperaturein[ttpoints];
	time[ns] = timein[ttpoints];


	for (j = 0; j <= labsteps; j++) // added this-- Lovera's code did not initialize arrays to zero for summation! }
	{
		tra40[j] = 0.0;
	}

// -------------------------- MAIN LOOP ------------------------------------

	for (j = 1; j <= ndomains; j++)   // calculate the age spectrum for each domain
	{

	// following code is Lovera's GEOSEC procedure
		zit = 0.0;
		for (jj = 1; jj <= ns; jj++)
		{
			avtemperature = 273.15 + (temperature[jj + 1] + temperature[jj]) / 2; //include conversion to K
			d[jj] = D0[j]/slab * exp(-Eact[j] / R / avtemperature);
		}

		dzita[ns] = 0.0;

		for (nj = 2; nj <= ns; nj++)
		{
			jj = nj - 1;
			dzita[ns - jj] = dzita[ns + 1 - jj] + fabs(time[ns + 1 - jj] - time[ns - jj]) * d[ns - jj];
		}

/* Computation of HXI(M)
   This series is VERY slow to converge for some values. To handle this, follow Oscar's approach
   of breaking series into three blocks and only evaluating at increasingly more widely separated indices. Every n,
   then every 5 n, and then every 200 n. 
*/
		for (m = 0; m <= 49; m++)
		{
			sum = 0.0;
			for (n = 1; n <= ns-1; n++)
			{
				a1 = pow(((m*geometry + 1) * pi),2.0) * d[n] - xlambd; 
				uplus = pow(((m*geometry + 1) * pi),2.0) * dzita[n + 1] + xlambd * (time[1] - time[n + 1]);
				if (uplus <= euplus) // was 50 in original code; now set as macro in arvert.h
				{
					if ((a1 * fabs(time[n + 1] - time[n])) > ecam)
					{
						cam = 1.0;
					}
					else
					{
						cam = 1.0 - exp(-a1 * fabs(time[n + 1] - time[n]));
					}     
					sum = sum + d[n] * pow(((m*geometry + 1) * pi),2) / a1 * cam * exp(-uplus);
				}     
			}  
			xi[m+1] = sum;
		}        // end of for m loop

		for (m = 1; m <= 310; m++)
		{
			mfake = m*5 + 49;
			sum = 0.0;
			for (n = 1; n <= ns-1; n++)
			{
				a1 = pow(((mfake*geometry + 1) * pi),2.0) * d[n] - xlambd;
				uplus = pow(((mfake*geometry + 1) * pi),2.0) * dzita[n + 1] + xlambd * (time[1] - time[n + 1]);
				if (uplus <= euplus)  // was 50 in original code; euplus now set as macro in arvert.h
				{
					if ((a1 * fabs(time[n + 1] - time[n])) > ecam)
					{
							cam = 1.0;
					}
					else
					{
							cam = 1.0 - exp(-a1 * fabs(time[n + 1] - time[n]));
					} 
					sum = sum + d[n] * pow(((mfake*geometry + 1) * pi),2.0) / a1 * cam * exp(-uplus);
				} 
			}
			xi[mfake+1] = sum;
		
			for (madj = 1;madj <= 4;madj++)
			{
				xi[mfake+1-madj] = xi[mfake+1] - (xi[mfake+1] - xi[mfake+1-5])/5*madj;
			}
		}         // end for m loop

		for (m = 1; m <= 242; m++)
		{
			mfake = m*200 + 1599;
			sum = 0.0;
			for (n = 1; n <= ns-1; n++)
			{
				a1 = pow(((mfake*geometry + 1) * pi),2.0) * d[n] - xlambd;
				uplus = pow(((mfake*geometry + 1) * pi),2.0) * dzita[n + 1] + xlambd * (time[1] - time[n + 1]);
				if (uplus <= euplus) // was 50 in original code; euplus now set as macro in arvert.h
				{
					if ((a1 * fabs(time[n + 1] - time[n])) > ecam)
					{
						cam = 1.0;
					}
					else
					{
						cam = 1.0 - exp(-a1 * fabs(time[n + 1] - time[n]));
					} 
					sum = sum + d[n] * pow(((mfake*geometry + 1) * pi),2.0) / a1 * cam * exp(-uplus);
				}
			}        //  end for n loop
			xi[mfake+1] = sum;
			for (madj = 1;madj <= 199;madj++)
			{
				xi[mfake+1-madj] = xi[mfake+1] - (xi[mfake+1] - xi[mfake+1-200])/200*madj;
			}

		}   // end for m loop

//  following code is Lovera's LAB procedure. I've snipped out the 39-related pieces into lablovera() and just
// use the array it passes.

		zii= 0.0;			

		for (nt = 1;nt <= labsteps;nt++)
		{
			zita[nt] = zii + D0[j]/slab * timelab[nt] * exp(-Eact[j] / R / temperaturelab[nt]);
			zii = zita[nt];
		}

		for (jjj = 1;jjj <= labsteps;jjj++)
		{
			r40[jjj] = 0.0;
			k = -1;
			change40 = 1.0;
			prev_r40 = 0.0;

			while ( (change40 > ichange) && (k < 49999)) // loop only until values stop changing
			{
				k = k + 1;
				kpie = pow(((k*geometry + 1) * pi),2.0);
				if (zita[jjj] * kpie > elab)  // was 100 in original code; elab now set as macro in arvert.h
				{
					ab = 0.0;
				}
				else
				{
					ab = exp(-zita[jjj] * kpie) / kpie;
				}
				r40[jjj] = r40[jjj] + (c40 - 1 + xi[k+1]) * (1.0 / kpie - ab);
				change40 = (r40[jjj] - prev_r40) / r40[jjj];
				prev_r40 = r40[jjj];
			}
				r40[jjj] = bgeom * r40[jjj];
		}
		for (k = 1; k <= labsteps;k++)
		{
			tra40[k] = fraction[j] * r40[k] + tra40[k];
		}
	}   // of main loop in j

	for (j = 1; j <= labsteps;j++)      // final data reduction 
	{
		sra40[j] = tra40[j] - tra40[j - 1];
		sra39[j] = bulkloss39[sample][j] - bulkloss39[sample][j - 1];
		if ((sra40[j] <= 0.0) || (sra39[j] <= 0.0))
		{
			age[sample][j] = 0.0;
		}
		else
		{
			scura[j] = sra40[j] / sra39[j];
			age[sample][j] = log(1.0 + scura[j] / c39wd) / xlambd;
		}
	}

	return;
}