/* 

		Function mineral of Arvert 7.0.x
		P. Zeitler, Lehigh University
		3/31/2004 - variable changes
		11/2014 - pass parameters, not globals
		
		calculates a mineral age given a thermal history, using simple volume diffusion and radiogenic production

	IMPORTANT NOTES:
	
	1. This routine uses SPHERICAL diffusion geometry. This has two ramifications. First, the kinetic data you
	supply for the sample MUST have been derived for this geometry. Second, the grain dimension you specify
	for your sample MUST be the "effective" radius of the sphere having the same surface/volume ratio of your
	unknown (see discussion in paper by Meesters and Dunai), or the Ft-equivalent sphere (see discussion in paper by
	Ketcham et al.)
	
	2. This inversion works with ages UNCORRECTED for alpha loss. That is because the primary data you measure
	include the effects of alpha recoil, which can have some impacts on diffusion, in that alpha loss modifies
	the expected diffusion profile at grain edges (making this profile less steep). No matter what one does, it's
	an imperfect world: one either alpha-corrects a cooling age, which is not an ideal thing to do, especially 
	for a young but slowly cooled samples (these are rare unless you mess with boreholes). On the other hand, if you
	leave the samples uncorrected, you have to make the model incorporate alpha-loss. I have chosen to leave the
	data as primary as possible and put calculations into the model.
	
	One outcome of this is that constraining and calculated He ages, being "raw", will not directly equate to Tc in 
	cooling systems. This makes no difference to the math or the constraining power, but can be jarring to people
	when they see what appears to be a constraint appearing at a time that doesn't seem to be in the data.
	
	3. This model only uses He production from 238-U. The contributions from 235-U and 232-Th are ignored. To first
	and probably second order this makes no difference since it is the diffusion of He that is most important.
	There could be situations where either the different alpha range for Th and/or the different decay timescales for
	235-U and 232-Th relative to the diffusional timescale could conceivably cause small differences in results.
	But this is not likely to be important; if it were, it would show up in samples within the PRZ, where the
	system is most sensitive to the balance between diffusion rates and production rates. Keep in mind that
	the modern 235-U contribution of He is only like 5% of that from 238-U.
	
	4. April 2012. Allowed modeling of titanite, and also recast this routine as a general mineral-age constraint (can 
	use a stopping distance of zero and therefore model argon or U-Pb systems (albeit with U lamda!! this should have 
	only a tiny effect for PRZ-stagnant samples).
	
	11/24/2015 now pass along a precision parameter that can scale to execution speed, similar to the classic
	HeFTy good / better / best precisions.
	
	January 2022 Vetted for use in Arvert 7.x multi-mineral. No real changes since this is being called once per mineral
*/

#include "arvert.h"

#define	Rgas         1.987        // gas constant 
#define	lamda238  1.55125e-10
#define Umass     238e6
#define Ualphas   8

double mineral_age_vd(long ndL, long ndt, long mineral, double radius, double eact, double difzero, double modelduration, const double timein[NODES], const double temperaturein[NODES], long ttpoints)
{	
	double age;
	double stopping_distance;

//	long nodes;  // nodes is now a macro fixed in arvert.h
	long i, dtLast, dLLast;
	long nx, timestep, segment;
	long steptime, jkl;
	long iL, LL, LK, K2;
	long array_size;
	
	double U;

	double dt, dL, surfconc; 
	double dtsec;
	double gastotal;
	double H;
	double arg, proportion, temp, time;
	double Xo, YearTime, U_back_then;
	double rrlam, lam;

	double *prevconc, *Hestar, *dif;

// arrays for use in tridiagonal algorithm 
	double *A, *B, *C, *dd, *dv, *beta, *gamma;
// for alpha-loss correction
	double *alpha_loss;

// ***** ALLOCATE MEMORY

	array_size = ndL+2;

// arrays involving distance
	prevconc = (double *) calloc(array_size, sizeof(double));
	dd = (double *) calloc(array_size, sizeof(double));
	dv = (double *) calloc(array_size, sizeof(double));
	beta = (double *) calloc(array_size, sizeof(double));
	gamma = (double *) calloc(array_size, sizeof(double));
	A = (double *) calloc(array_size, sizeof(double));
	B = (double *) calloc(array_size, sizeof(double));
	C = (double *) calloc(array_size, sizeof(double));
	alpha_loss = (double *) calloc(array_size, sizeof(double));

// arrays involving time
	array_size = ndt+2;

	Hestar = (double *) calloc(array_size, sizeof(double));
	dif = (double *) calloc(array_size, sizeof(double));
	
// ***** Initialize some things

	U = 50;  // value doesn't matter for this purely volume-diffusion code

// do some things needed for U-He ages

	switch (mineral)
	{
		case 0:
			stopping_distance = stopping_apatite; // mean stopping distance in cm for 238-U alpha-recoil, apatite
			break;
		case 1:
			stopping_distance = stopping_zircon; // mean stopping distance in cm for 238-U alpha-recoil, zircon
			break;
		case 2:
			stopping_distance = stopping_titanite; // mean stopping distance in cm for 238-U alpha-recoil, titanite
			break;
		case 3:
			stopping_distance =   0.0; // otherwise assume that no alpha correction will be modeled (e.g. for Ar-Ar biotite...)
			break;				
	}
		
		radius = radius / 1.0e4;  // convert microns to cm for use in this volume-diffusion routine

	// Set up dLLast and get fractional volume for each domain
	// For legacy reasons much of this code uses 1-based arrays, sorry
		dLLast = ndL + 1;  // need one more grid points than there are distance steps, timesteps
		dtLast = ndt + 1;

	// make certain that 1-based arrays stay in bounds
		assert(dLLast <= array_size-1);
		assert(dtLast <= array_size-1);

	// Make the node spacings, etc.
		dL = radius/ndL;
		dt = modelduration/ndt;

	// Set surface boundary condition
		surfconc = 0.0;

	// ***** Calc diffusivities
	// Need to calculate diffusivities at each timestep, which will vary since temperature is a function of time.  
	// Temperatures are obtained by linear interpolation through specified T-t pairs.

		dif[1] = difzero * exp(-eact * 1000 / Rgas / (temperaturein[1] + 273.15)); // at start
		timestep = 2;
		time = modelduration - (timestep-1)*dt;
		for (segment = 1; segment <= ttpoints-1; segment++)
		{
			while (time >= timein[segment+1])
			{
				proportion = (time - timein[segment+1]) / (timein[segment] - timein[segment+1]);
				temp = proportion * (temperaturein[segment] - temperaturein[segment + 1]) + temperaturein[segment + 1] + 273.15;
			 	arg = -eact*1000 / Rgas / temp;
				dif[timestep] = difzero * exp(arg);
				timestep = timestep + 1;
				time = modelduration - (timestep-1)*dt;
			}
		}

	// ***** code from ancient procedure GeolParams;

	/* Set radiogenic helium production array.  Dt in years!!
	Note: the amount of 238-U was greater a long time ago. The non-linearity in the age equation reflects this. Because 	
	we linearly break elapsed time into little chunks, we would miss this if we didn't make an adjustment. So, we need 
	to nudge up the U content to reflect higher 238-U. This effect should be negligible for Tertiary models, bit is 
	increasingly important for longer runs. */

	// ignoring Th and certainly Sm for now

			U_back_then = U / exp(-lamda238 * modelduration*1.0e6);
			Hestar[1] = U_back_then/Umass * (exp(Ualphas * lamda238 * dt * 1.0e6) - 1);
			dtsec = dt*1.0e6*365.35*24*3600;
			timestep = 2;
			for (i = timestep; i <= dtLast; i++)
			{
				YearTime = modelduration - dt*(i-1);
				U_back_then = U / exp(-lamda238 * YearTime*1.0e6);
				Hestar[i] = U_back_then/Umass * (exp(Ualphas*lamda238 * dt*1.0e6) - 1);
			}
					
		for (i = 1; i <= dLLast; i++)
		{
			Xo = (i-1)*dL;
			if (Xo < (radius - stopping_distance))
			{
				alpha_loss[i] = 1.0;
			}
			else
			{
				alpha_loss[i] = 0.5 + ((0.5/Xo*(Xo*Xo+radius*radius-stopping_distance*stopping_distance))-Xo)/(2*stopping_distance);
			}
		}

	// ***** code from ancient procedure NUMDFSetup;

		prevconc[dLLast] = surfconc;
		for (i = 2;i <= ndL; i++) // assign starting concentration profile
		{
			prevconc[i] = 0.0;
		}
		prevconc[1] = prevconc[2]; // VITAL: this is the no-flow boundary condition at CENTER of crystal!

	// ***** code from ancient procedure NUMDF - Crank-Nicolson finite-difference for spherical geometry

		  // --------- Main time loop -------------------------------------- 
		for (steptime = 2; steptime <= dtLast; steptime++)
		{
			rrlam = dif[steptime] * dtsec/(dL * dL);				
	      // ----------- Main depth loop --------------------------
	      		// Deal with boundary conditions:
	         // note that CENTER is at nx = 0 !!
	        
			A[2] = 0.0;
			B[2] = -(2.0*rrlam + 2.0);
			C[2] = 2.0*rrlam;  
	        
			dd[2] = (2.0 * rrlam - 2.0) * prevconc[1] - 2.0 * rrlam * prevconc[2];
			dd[2] = dd[2] - 2.0 * alpha_loss[2] * Hestar[steptime];

			for (nx = 3; nx <= ndL; nx++)
			{
				lam = rrlam/(nx - 1);
				A[nx] = (nx - 2.0)*lam;
				B[nx] = -(2.0 + 2.0*(nx - 1) * lam);
				C[nx] = nx * lam;
				
				dd[nx] = -lam * (nx - 2.0) * prevconc[nx - 1] + (2.0 * (nx - 1) * lam - 2.0) * prevconc[nx] - 2 * alpha_loss[nx] * Hestar[steptime];
				dd[nx] = dd[nx] - lam * nx * prevconc[nx + 1];
			}
	      // ----------- end depth loop ------------------------------
	// ***** code from extremely ancient procedure TRIDAG, From Numerical Recipes text by Carnahan, Luther and Wilkes

			beta[2] = B[2];
			gamma[2] = dd[2] / B[2];
			LL = 3;
			for (iL = LL;iL <= ndL;iL++)
			{
				lam = rrlam/(iL - 1);
				beta[iL] = B[iL] - (A[iL] * C[iL - 1])/beta[iL - 1];
				gamma[iL] = (dd[iL] - A[iL] * gamma[iL - 1]) / beta[iL];
			}
			dv[ndL] = gamma[ndL];
			LK = ndL - 2;
			for (K2 = 1;K2 <= LK;K2++)
			{
				jkl = ndL - K2;
				dv[jkl] = gamma[jkl] - (C[jkl] * dv[jkl + 1]) / beta[jkl];
			}
	//	 end of tridiagonal code

	// Assign concentrations calculated by TRIDAG code
			for (nx = 2;nx <= ndL;nx++)
			{
				prevconc[nx] = dv[nx];
			}
			prevconc[1] = prevconc[2];  // no-flow boundary condition at center
// Note that more sophisticated implementations might have a fake node at -1 to nuance this no-flow condition at center
			prevconc[dLLast] = surfconc; // surface boundary condition (O for this purpose!)
		}

	//	end of crank-nicolson code


	// FOLLOWING IS SUMMARY ROOUTINE - very simple trapezoidal integration of concentration profile, of course
	// adjusted for spherical geometry
		gastotal = 0.0;
		for (nx = 1; nx <= ndL;nx++)  // get amount of gas in crystal at end of last step (note adjustment for spherical geometry)
		{
			H = (prevconc[nx] + prevconc[nx + 1]) / 2;
			gastotal = gastotal + H * 4.1889*((nx*dL)*(nx*dL)*(nx*dL) - ((nx-1)*dL)*((nx-1)*dL)*((nx-1)*dL));
		}
		age = 1/(lamda238*1.0e6) * log((gastotal/(4.1889*radius*radius*radius)/Ualphas) / (U/Umass) + 1); 

// ************************** Wrap things up... Doom if you omit the following.  **************************

	free(prevconc);
	free(Hestar);
	free(dif);
	free(dd);
	free(dv);
	free(beta);
	free(gamma);
	free(A);
	free(B);
	free(C);
	free(alpha_loss);
	
	return age;
}