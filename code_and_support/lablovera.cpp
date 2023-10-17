/* 	
		Function lablovera() of Arvert 7.0.x
		P. Zeitler, Lehigh University
		modified from 6.1.3 code
		
   Based on C then C++ translation of Oscar Lovera's fortan program ages.f (translated through intermediate
   Pascal translation, also by Zeitler!). Somehow it (still) works.
   
   Based on lovera() from Arvert 6.1.3
   This routine just calculates the 39Ar release in the lab, which only needs to be done once per MDD sample.

*/

#include "arvert.h"
#include <vector>
using namespace std;

void lablovera(long sample, long ttpoints, long geometry, double bgeom, double slab, long ndomains, double Eact[MAXDOMAINS], double D0[MAXDOMAINS], double fraction[MAXDOMAINS], long labsteps, double dt, double ichange, const double timelab[MAXLABSTEPS],const double temperaturelab[MAXLABSTEPS], double bulkloss39[][MAXLABSTEPS])
{
	double r39[MAXLABSTEPS];
	
	long j, ns, i;
	long nt;

// following are from lab procedure of Lovera
	double  kpie, prev_r40, prev_r39, change39, change40, zii, ab;
	double zita[MAXLABSTEPS];
	long jjj, k;


	std::vector<double> temperature(nch);
	std::vector<double> time(nch);
	std::vector<double> dzita(nch);
	std::vector<double> d(nch);	

	//  ------------ Set some needed parameters -------------

	for (j = 0; j <= labsteps; j++) // added this-- Lovera's code did not initialize arrays to zero for summation! }
	{
		bulkloss39[sample][j] = 0.0;
	}

// -------------------------- MAIN LOOP ------------------------------------

	for (j = 1; j <= ndomains; j++)   // calculate the age spectrum for each domain
	{

//  following code is Lovera's LAB procedure

		zii= 0.0;			

		for (nt = 1;nt <= labsteps;nt++)
		{
			zita[nt] = zii + D0[j]/slab * timelab[nt] * exp(-Eact[j] / R / temperaturelab[nt]);
			zii = zita[nt];
		}

		for (jjj = 1;jjj <= labsteps;jjj++)
		{
			r39[jjj] = 0.0;
			k = -1;
			change39 = 1.0;
			prev_r39 = 0.0;

			while ((change39 > ichange)  && (k < 49999)) // loop only until values stop changing or terms is high
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
					r39[jjj] = r39[jjj]  + ab;
				}

				change39 = (r39[jjj] - prev_r39) / r39[jjj];
				prev_r39 = r39[jjj];
			}
				r39[jjj] = 1.0 - bgeom * r39[jjj];
		}
		for (k = 1; k <= labsteps;k++)
		{
			bulkloss39[sample][k] = fraction[j] * r39[k] + bulkloss39[sample][k];
		}
	}   // of main loop in j

	return;
}