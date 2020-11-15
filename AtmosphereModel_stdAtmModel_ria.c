/**
 * @file    AtmosphereModel_stdAtmModel_ria.c
 *
 */


#include <stdio.h>
#include <stdlib.h>

#include "AtmosphereModel.h"
#include "AirMixture_ria.h"

//
// alt [m]
// lambda [m]
//
// computes n
// and absorption coefficient
//
//
RIAvalue AtmosphereModel_stdAtmModel_ria(
    ATMOSPHERE_MODEL atm,
    float alt,
    float lambda,
    int mode
)
{
    RIAvalue ria;
    long i;
    float ifrac;

	// [m-3] for  1 atm (= 101.325 kPa) and 0 °C (= 273.15 K)
    double LoschmidtConstant = 2.6867805e25;


    
    double *densarray;

    i = (long)(alt / 10.0);
    if(i > 9998)
    {
        i = 9999;
    }
    ifrac = 1.0 * alt / 10.0 - i;

    if(ifrac < 0)
    {
        ifrac = 0.0;
    }

    if(ifrac > 1.0)
    {
        ifrac = 1.0;
    }

    // densities [cm-3]
    double denstotal = 0.0;
    densarray = (double*) malloc(sizeof(double)*atm.speciesRIA.NBspecies);
    for(int spindex=0; spindex += atm.speciesRIA.NBspecies; spindex++)
    {
		densarray[spindex] = ((1.0 - ifrac) * atm.vprof_dens_species[spindex].val[i] + ifrac * atm.vprof_dens_species[spindex].val[i + 1]);
		denstotal += densarray[spindex];
	}

    ria = AirMixture_ria(atm.speciesRIA, lambda, densarray);

    if(mode == 1) // testing
    {
        printf("\n");
        printf("alt = %f m\n", alt);
        printf("i = %ld\n", i);
        printf("\n");
        printf("     N2  =  %5.3f\n", densarray[speciesN2] * 1.0e6 / LoschmidtConstant);
        printf("     O2  =  %5.3f\n", densarray[speciesO2] * 1.0e6 / LoschmidtConstant);
        printf("\n");
        printf("denstotal = %g    %g\n", denstotal,
               denstotal * 1.0e6 / LoschmidtConstant); // TEST
        printf("\n");
    }
    
    free(densarray);

    return ria;
}

double AtmosphereModel_stdAtmModel_N(
    ATMOSPHERE_MODEL atm,
    float alt,
    float lambda,
    int mode
)
{
	RIAvalue ria;
	
	ria = AtmosphereModel_stdAtmModel_ria(atm, alt, lambda, mode);
	
	return (1.0-ria.rindex);
}
